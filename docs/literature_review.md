# Literature review: co-occurrence and mutual exclusivity of cancer alterations

What published work does, where it agrees (the standard to follow), where it
differs (the choices we have to make), and what this project takes from each.
Every design choice at the end was **tested on GENIE** before being adopted; the
scripts are in `scripts/analyses/`.

---

## 1. Papers reviewed

| # | Paper | Data | What it does | Key finding for us |
|---|---|---|---|---|
| 1 | Canisius, Martens & Wessels 2016, *Genome Biol* — **DISCOVER** | TCGA pan-cancer | Pair test with a Poisson-Binomial null that gives every tumour its own alteration rate | Assuming equal rates across tumours (Fisher) creates spurious co-occurrence and misses exclusivity; "biology drives mutual exclusivity but chance explains most co-occurrence" |
| 2 | Ferrer-Bonsoms et al. 2022, *Bioinformatics* — **Rediscover** | — | Fast R implementation of DISCOVER | The tool most GENIE-scale studies could use; we benchmarked it (it inherits DISCOVER's bias towards exclusivity, 69% false positives on our null) |
| 3 | Kim et al. 2017, *Bioinformatics* — **WeSME** | TCGA | P-values that control each patient's and each gene's mutation rate, approximating a margin-preserving permutation test | Same idea as DISCOVER, reached by permutation; exclusivity can also reveal mutational processes (APOBEC) |
| 4 | Leiserson et al. 2016, *Bioinformatics* — **WExT**; Leiserson et al. 2015 — **CoMEt** | TCGA colorectal / endometrial | Exact exclusivity test conditioning on per-sample, per-gene probabilities; extends to sets of >2 genes | Highly variable mutation rates (hypermutated colorectal / endometrial) confound exclusivity unless per-sample rates are used |
| 5 | **van de Haar** et al. 2019, *Cell* — "Identifying epistasis in cancer genomes: a delicate affair" | TCGA | Re-examines published exclusivity maps | Most exclusivity is explained by **disease subtype and tumour mutation load**, not pathway structure; driver genes are mutated preferentially in low-load tumours, so they look exclusive with many genes |
| 6 | **Park & Lehner** 2015, *Mol Syst Biol* | >3,000 TCGA tumours, 22 cancer types | Tests within each cancer type; odds-ratio heterogeneity test (with permutations) for whether an interaction differs between types | Interactions **change between cancer types** ("plasticity of epistasis"); 30 of 52 re-testable interactions were cancer-type-specific |
| 7 | **Mina** et al. 2017, *Cancer Cell* — **SELECT**; Mina et al. 2020, *Nat Genet* | 6,456 / >9,000 TCGA tumours | Evolutionary dependencies between *functional* alterations; validation against CRISPR/drug screens in 2,000 cell lines | Exclusivity ≈ functional redundancy, co-occurrence ≈ synergy (e.g. PIK3CA + NFE2L2 in squamous cancers); dependencies predict therapeutic response |
| 8 | Iranzo et al. 2022, *Cell Rep* — **Coselens** | TCGA, 16 cancer types | Conditional selection measured from excess non-synonymous mutations (dN/dS-style) | 296 conditionally selected pairs; epistasis crosses pathways and shapes subtypes |
| 9 | **Scharpf** et al. 2022, *Cancer Res* — mutant RAS in **GENIE** | 66,372 GENIE tumours, 51 cancer types | Bayesian hierarchical models of co-mutation, modelling **TMB and mutational signatures** | Co-mutation is lineage- and **allele-specific**: KRAS G12C lung cancer enriched for NTRK3 and chromatin regulators |
| 10 | Cook et al. 2021, *Nat Commun* | 13,492 tumours | KRAS allele-level co-mutation networks per tissue | Each KRAS allele has its own tissue-specific co-mutation network |
| 11 | Vaeyens et al. 2023, *medRxiv* | 64,807 cBioPortal samples + cell lines | Variant-class-level exclusivity of BRAF / EGFR / KRAS | Exclusivity depends on the variant class (class I BRAF, "classical" EGFR), not just the gene |
| 12 | Sanchez-Vega et al. 2018, *Cell* | 9,125 TCGA tumours, 64 subtypes | Curated pathway alterations; co-occurrence/exclusivity within and between 10 pathways | Standard "functional alteration" definitions; exclusivity within pathways, co-occurrence between them |
| 13 | El Tekle et al. 2021, *Trends Cancer* | review | Categories of co-occurrence and incompatibility | Two kinds of exclusivity: same-pathway redundancy vs divergent-pathway incompatibility (e.g. oncogene-induced senescence) |
| 14 | Skoulidis & Heymach 2019, *Nat Rev Cancer* | review | Co-mutations in lung cancer biology and therapy | Co-mutations (TP53, STK11, KEAP1) change prognosis and drug response — the clinical reason for Aim 1 |
| 15 | Kuipers et al. 2021, *PLoS Comput Biol* — **GeneAccord** | single-cell AML | Clonal co-occurrence / exclusivity inside a tumour | Sample-level co-occurrence can hide subclonal structure; needs clonal reconstruction |
| 16 | Campbell & Reyna 2026, *bioRxiv* — **MOSAIC** | — | Subgroup-aware exclusivity | Existing methods cannot tell exclusivity from interaction apart from exclusivity from subgroup structure |
| 17 | Recent **GENIE disease papers** (bladder, gastric, GIST, angiosarcoma, chordoma, ... 2024-26) | GENIE subsets | Descriptive landscape + pairwise co-occurrence p-values | Mostly raw pairwise tests without burden, subtype or panel handling (see §3) |

## 2. What (almost) everyone agrees on — the standard to follow

1. **Never assume every tumour has the same chance of every alteration.**
   DISCOVER, WeSME and WExT all build the null so that each tumour's own
   alteration count (or rate) is respected. Fisher's test (and cBioPortal's
   co-occurrence tab, which uses it) does not.
2. **Test within a cancer type.** Interactions change between cancer types
   (Park & Lehner); pan-cancer results are summaries of within-type results.
3. **Use functional (driver-level) alterations**, not every variant call
   (SELECT, Sanchez-Vega).
4. **Correct for multiple testing** (FDR) and expect exclusivity to be the more
   reliable signal; co-occurrence is easily produced by chance (DISCOVER).
5. **Treat mutation load and subtype as confounders**, not background
   (van de Haar, MOSAIC, Scharpf).

## 3. The common "old" approach, and why it falls short

Most single-disease GENIE papers report pairwise co-occurrence p-values on raw
counts, typically through cBioPortal, whose co-occurrence tab uses a one-sided
Fisher's test. Problems, each shown on our data:

| Problem | Evidence on GENIE |
|---|---|
| Ignores mutation load | On null data built from real **non-hypermutated** patients, Fisher calls 6-15% of pairs significant at p < 0.05 (should be 5%); the conditional test 1.3-2.9% |
| Hypermutated tumours pooled in | Colorectal: 85% of pairs "significant" with Fisher on all patients vs 8% in non-hypermutated. Example from the literature: POLD1–KMT2D–ARID1A "co-occurrence" in gastric cancer — three favourite targets of hypermutation |
| Subtypes pooled | "Pancreatic cancer" mixes adenocarcinoma (SMAD4) and neuroendocrine tumours (MEN1): MEN1 vs SMAD4 looks strongly exclusive (q = 6e-22) and disappears once subtype is accounted for |
| Panels ignored | A gene absent from a panel is counted as "not mutated"; 20 panels claim copy-number data they never released |

## 4. Where methods differ — and our choice

| Question | Options in the literature | Our choice and why |
|---|---|---|
| How to model mutation load | Plug-in per-tumour rates (DISCOVER, Rediscover); margin-preserving permutation (WeSME); per-sample weights (WExT); Bayesian hierarchical with TMB/signatures (Scharpf) | **Conditional (Rasch) test** — the exact version of the margin-preserving permutation null. The plug-in estimate biases towards exclusivity on panel-sized gene sets (46% false positives on our benchmark; Rediscover 69%) |
| Hypermutated tumours | Ignore; exclude; model TMB as covariate (Scharpf) | **Analyse non-hypermutated patients as the main analysis** (Jason's suggestion), hypermutated (TMB 10-100) separately as a description, ultramutated (>= 100) excluded. Splitting alone is not enough (Fisher still inflated), so the conditional test is kept |
| Subtype | Usually ignored; van de Haar and MOSAIC show it matters | **Gene rates fitted within detailed subtypes** (>= 50 patients) — see §5 |
| Gene vs variant | Gene level (most); allele level (Scharpf, Cook, Vaeyens) | Gene level in the screen; **variant level in the lookup tool**, and a variant-level screen for key hotspot genes in Phase 4 |
| What to rank by | p-value (DISCOVER); effect-based statistics (SELECT); selection strength (Coselens) | Keep q-value, report fold enrichment and **flag weak effects**: in the previous analysis 25% of candidates differed from chance by < 1.5-fold (KRAS–TP53 in pancreas: fold 1.08, q = 2e-24) |
| Validation | None (most GENIE papers); cell-line screens (Mina 2020); cross-cohort (Park & Lehner) | **Out-of-centre replication**: discover in MSK, re-test in all other centres (§5); plus simulation-based null calibration |
| Multi-gene sets | CoMEt, WExT, MEMo | Not in the screen; the lookup tool counts exact combinations directly |
| Cancer-type specificity | Park & Lehner heterogeneity test | Adopt for Aim 2 (Phase 4) |
| Clonality | GeneAccord (needs subclonal data) | **Not adopted** — panel data has no clonal reconstruction; noted as a limitation |
| Selection-based (dN/dS) | Coselens | **Not adopted** — needs exome/genome-wide passenger mutations that panels do not have |

What is unusual about this project compared with all of the above: it handles
**panel testability** explicitly (which genes were actually tested, for
mutations and for copy number, per patient). None of the reviewed methods
does this; most were built on whole-exome TCGA data where every gene is tested.

## 5. Tests we ran on GENIE to decide

| Test | Result | Script |
|---|---|---|
| Fisher vs conditional on null data from non-hypermutated patients | Fisher 6-15% false positives at p < 0.05 (up to 1.8% of pairs after FDR); conditional 1.3-2.9%, 0% after FDR | `04_null_calibration_nonhyper.py` |
| TMB split: Fisher / conditional, all vs non-hypermutated | Splitting cuts Fisher's calls sharply (colorectal 85% -> 8%) but Fisher stays above the conditional test (lung 31% vs 24%) | `02_tmb_split_fisher_vs_conditional.py` |
| Ultramutated tumours inside the hypermutated group | Endometrial: 21.6% significant pooled; 0.3% (TMB 10-100) and 0% (>= 100) when separated | `03_ultramutated_split.py` |
| Subtype-specific gene rates | Removes 390 of 1,571 hits in 10 cancer types, concentrated in melanoma, endometrial, ovarian and glioma; the removed ones are textbook subtype markers (glioma CIC vs PTEN = oligodendroglioma vs glioblastoma; breast CDH1 vs ERBB2 = lobular vs ductal; ovarian PIK3CA vs TP53 = clear-cell vs serous; pancreatic MEN1 vs SMAD4) | `05_subtype_and_replication.py` |
| Replication MSK -> other centres | 96% of MSK hits keep the same direction elsewhere, none flips significantly. Hits that survive subtype adjustment replicate **73%** of the time, hits that do not **35%** | `05_subtype_and_replication.py` |
| Whole-design comparison by replication | Discover in MSK patients, re-test the same pairs in all other centres; "replicated" = p < 0.05, same direction (10 cancer types). A: current; B: non-hypermutated only; C: B + subtype-specific gene rates; D: C + gene filter counted among tested patients; E: C on mutations only. Replication rates: **A 44%, B 66%, C 66%, D 64%, E 84%**. Significant pairs: A 3,762, B 1,571, C 1,229, D 1,525, E 419. Caveat: artefacts shared by all centres (hypermutation, subtype) also "replicate", so replication is a floor on reliability, not proof | `06_design_comparison.py` |
| Filters (gene frequency, joint coverage, no gene filter) | Dropping the gene-frequency filter multiplies pairs ~10x and run time to 1-20 h per cancer type. The joint-coverage filter removes few pairs (287 in 5 cancer types, 222 in breast), and those rest on one centre's panel: only 1 of them could even be re-tested outside MSK | `p2_variants.py` grid, `07_joint_coverage_filter.py` |
| Mutation reporting vs panel gene lists | 1,144 (panel, gene) pairs report zero mutations where >= 10 are expected; 2 panels release no mutations at all | `01_mutation_panel_reporting.py` |

## 6. What this changes (implemented 2026-10-03)

Result of the redesigned Phase 2: 6,636 pairs in non-hypermutated patients
(51 cancer types), 1,783 significant, **1,543 interaction candidates**, of which
**464 replicate outside MSK** (65% of testable MSK hits replicate, none flips).


**Phase 2 (co-occurrence screen)**

| Change | Why (evidence) |
|---|---|
| Main analysis on **non-hypermutated patients** (TMB < 10 on panels >= 1 Mb) | Jason's suggestion; replication 44% -> 66%; the pooled analysis spent most tests on genes frequent only because of hypermutation (melanoma 6,873 pooled pairs vs 488) |
| Keep the **conditional test** | Fisher remains miscalibrated in non-hypermutated patients (6-15% false positives) |
| Gene rates fitted **within detailed subtypes** (>= 50 patients) | van de Haar 2019, Park & Lehner 2015; removes textbook subtype artefacts (MEN1 vs SMAD4, CIC vs PTEN, CDH1 vs ERBB2); subtype-robust hits replicate 73% vs 35% |
| Keep the four filters; each has a different job | F1 enough patients; F2 gene altered often enough (and, as defined over all patients, widely tested); F3 pair not resting on one centre's panel; F4 enough power. Counting F2 among tested patients adds pairs that replicate worse; F3 removes only single-centre pairs |
| **Hypermutated patients (TMB 10-100) as their own analysis**, reported as a description | Their patterns are largely MSI-vs-MSS subtype structure; ultramutated (>= 100) excluded |
| Add a **mutation-only analysis** on all patients with mutation data | Uses more patients (no copy-number restriction) and replicates best (84%); complements the mutation-or-copy-number analysis rather than replacing it (copy-number biology such as CDKN2A deletion, ERBB2 amplification is only in the latter) |
| Add **out-of-centre replication** columns per pair | Standard validation (Park & Lehner re-testing); gives the lookup tool a confidence level per pair |
| **Flag weak effects** (< 1.5-fold either way, `weak_effect`) | Large cohorts make tiny effects significant: 25% of the previous candidates and 22% of the new ones (346 of 1,543) are within 1.5-fold of chance |

**Phase 1:** mutation testability read partly from the data, as for copy number —
genes count as tested only where the panel has exon probes, and samples on panels
that released no mutations (VICC-02-XFV2, CHOP-FUSIP) count as not tested for
mutations. Gene-level gaps that can only be found statistically are left as they
are (documented in the README).

**Phase 3 (lookup tool):** already variant-level (Scharpf 2022, Cook 2021, Vaeyens
2023 show co-mutation is allele-specific); add the subtype filter, the replication
flag and the effect size from Phase 2.

**Phase 4 (Aim 2):** cancer-type-specific interactions with an odds-ratio
heterogeneity test (Park & Lehner 2015); variant-level screen for hotspot genes
(KRAS alleles, BRAF V600E vs non-V600, EGFR L858R / exon 19); interpret exclusive
pairs as same-pathway redundancy vs divergent-pathway incompatibility (El Tekle 2021).

**Not adopted:** clonality tests (GeneAccord — panel data has no subclonal
reconstruction); dN/dS-based selection (Coselens — needs genome-wide passengers);
cell-line validation (Mina 2020 — out of scope, a possible extension).

## References

1. Canisius S, Martens JWM, Wessels LFA. *Genome Biol* 2016;17:261. doi:10.1186/s13059-016-1114-x
2. Ferrer-Bonsoms JA et al. *Bioinformatics* 2022;38:844-5. doi:10.1093/bioinformatics/btab709
3. Kim YA et al. *Bioinformatics* 2017;33:814-21. doi:10.1093/bioinformatics/btw242
4. Leiserson MDM et al. *Bioinformatics* 2016;32:i736-45. doi:10.1093/bioinformatics/btw462 — and *Genome Biol* 2015;16:160 (CoMEt)
5. van de Haar J et al. *Cell* 2019;177:1375-83. doi:10.1016/j.cell.2019.05.005
6. Park S, Lehner B. *Mol Syst Biol* 2015;11:824. doi:10.15252/msb.20156102
7. Mina M et al. *Cancer Cell* 2017;32:155-68. doi:10.1016/j.ccell.2017.06.010 — and *Nat Genet* 2020;52:1198-207. doi:10.1038/s41588-020-0703-5
8. Iranzo J et al. *Cell Rep* 2022;40:111272. doi:10.1016/j.celrep.2022.111272
9. Scharpf RB et al. *Cancer Res* 2022;82:4058-78. doi:10.1158/0008-5472.CAN-22-1731
10. Cook JH et al. *Nat Commun* 2021;12:1808. doi:10.1038/s41467-021-22125-z
11. Vaeyens F et al. *medRxiv* 2023. doi:10.1101/2023.10.21.23297089
12. Sanchez-Vega F et al. *Cell* 2018;173:321-37. doi:10.1016/j.cell.2018.03.035
13. El Tekle G et al. *Trends Cancer* 2021;7:823-36. doi:10.1016/j.trecan.2021.04.009
14. Skoulidis F, Heymach JV. *Nat Rev Cancer* 2019;19:495-509. doi:10.1038/s41568-019-0179-8
15. Kuipers J et al. *PLoS Comput Biol* 2021;17:e1009036. doi:10.1371/journal.pcbi.1009036
16. Campbell K, Reyna MA. *bioRxiv* 2026. doi:10.64898/2026.04.29.721672
17. AACR Project GENIE Consortium. *Cancer Discov* 2017;7:818-31. doi:10.1158/2159-8290.CD-17-0151
