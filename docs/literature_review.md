# Literature review: co-occurrence and mutual exclusivity of cancer alterations

This review maps what has been done on the question this project asks, places the project in
that work, and records how each published method family performed **when we ran it (or its core
idea) on our own GENIE data**. It is organised to serve the six purposes of a literature review:
a map of the field (§1-5), positioning and gaps (§9), avoiding past errors (§3, §6), choosing
methods (§6-8), synthesised insight (§8, §10), and documented scholarly method (§0).

Scripts for every test are in `scripts/analyses/`; numbers are from GENIE v19.0 unless stated.

---

## 0. How this review was done

**Question.** How should co-occurrence and mutual exclusivity of somatic alterations be tested
in large, multi-centre, targeted-panel cohorts such as GENIE -- and what has gone wrong before?

**Sources searched** (2026-09-30 to 2026-10-22):
- Europe PMC REST API (covers PubMed/MEDLINE, PubMed Central and preprints), sorted by citations;
- Firecrawl research index (PubMed, bioRxiv, medRxiv, arXiv; semantic search);
- tool documentation and code (cBioPortal manual, maftools, PathwayMapper, CIViC) for methods
  actually used in practice.

**Search terms** (combined in several framings, because one query misses papers another finds):
"mutual exclusivity" / "co-occurrence" + "somatic" / "cancer" + "mutation rate" / "heterogeneous";
"epistasis cancer genomes confounders subtype mutation load"; "GENIE co-occurrence panel";
"allele-specific co-mutation KRAS"; "conditional selection driver"; "chromosomal instability
copy number co-occurrence"; "fusion mutual exclusivity driver"; "patient similarity panel
sequencing"; "clinical co-mutation NSCLC"; "driver pathway gene sets mutual exclusivity".

**Included:** methods for testing pairs or sets of alterations in tumour sequencing; studies of
confounders (mutation load, subtype, tissue, copy-number complexity); GENIE and large-panel
co-mutation studies; clinical reviews of co-mutations; data resources and tools used in practice.
**Excluded:** germline co-occurrence; non-cancer studies; single case reports, except where they
document a rare combination.

**Appraisal.** Each method paper was read for (1) what its null hypothesis assumes, (2) which
confounders it handles, (3) what data it was built on, (4) how it was validated -- and then, where
possible, run on our data (§7). Every reference was checked against Europe PMC metadata
(authors, journal, volume, pages, DOI); three errors were caught this way (a wrong first author,
a follow-up paper returned instead of the original Dendrix paper, a correction notice instead of
the PCAWG paper).

**Limitations.** Narrative, single-reviewer review, not a PRISMA systematic review; English only;
full texts were read for the papers the design rests on, including van de Haar 2019 and Scharpf
2022 (PMC6816465, PMC9627127) and SELECT (Mina 2017, Methods S1), checked statement by statement on
26 Oct 2026; a few papers cited for a single point (e.g. El Tekle 2021, Vaeyens 2023) were read as
abstracts only.

---

## 1. The problem

A tumour usually carries a few **driver alterations** -- mutations, copy-number changes or gene
fusions that help it grow. Drivers do not combine at random:

| Pattern | Meaning | Typical explanation | Example |
|---|---|---|---|
| **Mutual exclusivity** | two genes are rarely altered in the same tumour | both hit the same pathway, so one is enough (*redundancy*); or both together harm the cell (*incompatibility*) | KRAS vs EGFR in lung cancer -- both activate RAS signalling |
| **Co-occurrence** | two genes are altered together more often than chance | they cooperate, or one makes the other tolerable | KRAS + STK11 in lung cancer |

**Why clinicians care:** a co-alteration can change prognosis or drug response -- *STK11*/*KEAP1*
alongside *KRAS* predicts poor immunotherapy response (Skoulidis & Heymach 2019); *BRAF* V600E +
*NRAS* predicts resistance to BRAF inhibitors (CIViC); *TP53* or *CDKN2A/B* loss in ALK/RET/ROS1
fusion-positive lung cancer predicts poorer targeted-therapy outcomes (Jiang 2023).

**Why it is hard:** most apparent patterns are not biology. They come from
- **mutation load** -- a few tumours carry most alterations, so everything "co-occurs" in them;
- **subtype** -- a "cancer type" can mix diseases with different drivers, which then look exclusive;
- **genomic position / copy-number complexity** -- one amplification hits neighbouring genes, and
  chromosomally unstable tumours gain many genes at once;
- **testing gaps** -- a gene not on the panel is read as "not altered".

Most methods below are attempts to remove one or more of these.

### Glossary

| term | meaning |
|---|---|
| null model | what "chance" means in the test -- the most important choice in this field |
| Fisher's exact test | classic 2x2 test; assumes every tumour has the same chance of every alteration |
| Poisson-Binomial | distribution of a count when each tumour has its own probability (DISCOVER, ours) |
| conditional (Rasch) test | our test: each tumour's chance is computed *given its own number of alterations* |
| permutation null | reshuffle the data keeping each tumour's and each gene's totals fixed (WeSME); our conditional test computes this null exactly |
| TMB | tumour mutational burden, mutations per million DNA letters; >= 10 = hypermutated |
| copy-number burden / CIN | how many genes a tumour has gained or lost / ongoing chromosomal instability |
| amplicon | a stretch of DNA amplified as one unit, carrying several genes together |
| subtype confounding | two alterations look exclusive only because each marks a different disease subtype |
| dN/dS | ratio of protein-changing to silent mutations; > 1 means a gene is under positive selection |
| FDR / q-value | false discovery rate; q < 0.05 = at most ~5% of called pairs expected to be false |
| replication | the same result in independent patients (other centres, or another cohort such as TCGA) |
| shrinkage / partial pooling | pulling a noisy estimate towards an average across related groups (hierarchical models) |
| over-adjustment | adjusting for a variable that is partly caused by the outcome, which hides real effects |

---

## 2. Data sources used in previous work

| Dataset | Type | Size | Strength | Weakness for this question |
|---|---|---|---|---|
| **TCGA PanCanAtlas** (Hoadley 2018) | whole exome + copy number, research cohort | ~10,000 tumours, 33 types | every gene tested in every tumour; matched normal; uniform processing | few patients per type, so rare pairs cannot be tested; mostly primary, untreated tumours |
| **ICGC/TCGA PCAWG** (2020) | whole genome | ~2,600 tumours | non-coding regions, structural variants | small |
| **MSK-IMPACT** (Zehir 2017) | clinical panel, one centre | ~10,000 at publication | real clinical samples incl. metastases | one centre, one panel family |
| **AACR GENIE** (2017-, this project) | clinical panels pooled from many centres | ~270,000 samples (v19) | by far the largest; real-world; enough patients for rare combinations | 167 panels, centre-specific pipelines, coarse clinical data |
| **GENIE BPC** | GENIE subset with curated treatment and outcome | thousands, 6 cancer types | treatment and survival | separate access, small |
| **Cell lines** (CRISPR / drug screens; used by Mina 2020) | functional screens | ~1-2,000 lines | can test whether a dependency is functional | not patients |

**Key point:** almost all method papers were built and validated on **TCGA**, where every gene is
tested in every tumour. GENIE breaks that assumption -- hence this project's per-patient
testability tracking for mutations, copy number and fusions (Phase 1). TCGA is used here the
other way round: as an **independent validation cohort** (§10).

---

## 3. What (almost) everyone agrees on -- and the common shortcut that ignores it

The standard to follow:
1. **Never assume every tumour has the same chance of every alteration** (DISCOVER, WeSME, WExT).
2. **Test within a cancer type** -- interactions change between tissues (Park & Lehner 2015).
3. **Use functional (driver-level) alterations**, not every variant call (SELECT; Sanchez-Vega 2018).
4. **Correct for multiple testing**, and expect exclusivity to be the more reliable signal (DISCOVER).
5. **Treat mutation load, subtype and copy-number complexity as confounders** (van de Haar 2019;
   Kumar 2011; Zack 2013).

The common shortcut: most single-disease GENIE papers (bladder, gastric, GIST, angiosarcoma,
chordoma, ... 2024-26) report pairwise p-values on raw counts, typically via cBioPortal, whose
co-occurrence tab uses Fisher's test. What that costs, shown on our data:

| Problem | Evidence on GENIE |
|---|---|
| Ignores mutation load | On null data built from real non-hypermutated patients, Fisher calls 6-15% of pairs significant at p < 0.05 (should be 5%); the conditional test 1.3-2.9% |
| Hypermutated tumours pooled in | Colorectal: 85% of pairs "significant" with Fisher on all patients vs 8% in non-hypermutated. Literature example: POLD1-KMT2D-ARID1A "co-occurrence" in gastric cancer -- three favourite targets of hypermutation |
| Subtypes pooled | "Pancreatic cancer" mixes adenocarcinoma (SMAD4) and neuroendocrine tumours (MEN1): MEN1 vs SMAD4 looks strongly exclusive (q = 6e-22) and disappears once subtype is accounted for |
| Panels ignored | A gene absent from a panel is counted as "not mutated"; 20 panels claim copy-number data they never released, 9 claim fusion detection they never released |

---

## 4. Method families

Each family: what it does, the data it was built on, and its **status here** --
*adopted*, *tested on our data*, *idea borrowed* (when the method cannot run on panel data as
published), or *not used*. Results of running them are in §7.

### A. Simple pairwise tests -- the common default
| Method | Paper | Idea | Status here |
|---|---|---|---|
| Fisher's exact / hypergeometric test | -- | 2x2 table per pair | tested: miscalibrated on our data (§7) |
| cBioPortal "Mutual Exclusivity" tab | Cerami 2012; Gao 2013 | one-sided Fisher + log odds ratio | = Fisher |
| maftools `somaticInteractions()` | Mayakonda 2018 | pairwise Fisher on a mutation table | = Fisher |

### B. Tests that respect each tumour's own alteration count -- the current standard
| Method | Paper | Data | Idea | Status here |
|---|---|---|---|---|
| **DISCOVER** | Canisius 2016 | TCGA pan-cancer | each tumour gets its own alteration probability (estimated), Poisson-Binomial null | tested: biased towards exclusivity on panel-sized gene sets (§7) |
| **Rediscover** | Ferrer-Bonsoms 2022 | -- | fast R implementation of DISCOVER | tested: same bias, 69% false positives on our null benchmark |
| **WeSME** | Kim 2017 | TCGA | approximates a permutation null keeping each patient's and gene's totals | idea adopted: our conditional test is this null computed exactly |
| **WExT** | Leiserson 2016 | TCGA colorectal/endometrial | exact exclusivity test with per-sample, per-gene weights | same family; superseded by the exact conditional test |
| **Conditional (Rasch) test** | this project | GENIE | each tumour's pair probability *given its own alteration count*; exact Poisson-Binomial tail | **adopted** (Phase 2) |

### C. Sets of mutually exclusive genes (more than two)
| Method | Paper | Idea | Status here |
|---|---|---|---|
| **MEMo** | Ciriello 2012 | exclusive modules within an interaction network | idea tested (§7) |
| **Dendrix / Multi-Dendrix** | Vandin 2012; Leiserson 2013 | gene sets maximising coverage minus overlap | tested on our data (§7) |
| **CoMEt** | Leiserson 2015 | exact test for exclusive sets | same aim |
| **Mutex** | Babur 2015 | exclusivity + shared signalling target | same aim |

Built on TCGA; used for pathway discovery rather than pair-level clinical questions. The lookup
tool counts exact multi-gene combinations instead (Phase 3).

### D. Selection- and evolution-based methods
| Method | Paper | Data | Idea | Status here |
|---|---|---|---|---|
| **SELECT** | Mina 2017; Mina 2020 | TCGA (>9,000 tumours in 2020) + ~2,000 cell lines (2020) | weighted mutual information with a per-class sample-count correction; empirical null from margin-preserving randomisations within tumour (sub)types; effect size and empirical FDR; 2020: 852 dependencies (476 exclusive, 376 co-occurring), robust to TMB, signatures and intratumour heterogeneity; in cell-line knockout screens co-occurrence meant greater dependence on the altered genes (synergy) and exclusivity meant one alteration rescues loss of the other (redundancy) | **re-implemented** from Methods S1 and benchmarked (§7); redundancy interpretation used in Phase 4 |
| **SelectSim** | Iyer 2026 [42] | GENIE v15 (MSK + DFCI only, per panel or shared genes) + TCGA; >70,000 tumours, 119 subtypes; OncoKB-oncogenic mutations only | simulate mutations under independence keeping gene frequency and TMB (within subtype); TMB-weighted co-mutation effect size; empirical FDR; validated across cohorts | compared pair by pair (notebook 04 §9): every robust pair it reports for the same cancer type agrees in direction; it has no copy number or fusions, so most of our copy-number pairs are outside its scope |
| **Coselens** | Iranzo 2022 | TCGA, 16 types | is a gene's dN/dS higher when another gene is mutated? (conditional selection) | idea borrowed: conditional nonsynonymous/silent ratio on GENIE panels (§7) |

### E. Confounders and context
| Paper | Data | Finding | Status here |
|---|---|---|---|
| **van de Haar 2019** | 13 TCGA + MSKCC cohorts, 8 tumour types, 341 cancer genes | five methods (DISCOVER, MEMo, WExT, Fisher, MEGSA) re-run: most mutual exclusivity they find tracks subtype and mutation load, not pathway structure -- even for load-aware DISCOVER; 78-100% of hits involve the 15 genes most mutated in low-load tumours | adopted: subtype strata + load-conditional test |
| **Park & Lehner 2015** | >3,000 TCGA tumours, 22 types | interactions change between cancer types; odds-ratio heterogeneity test | adopted: Phase 4 §2 |
| **Sanchez-Vega 2018** | 9,125 TCGA tumours | 10 curated oncogenic pathways; exclusivity common within pathways | adopted: pathway test, Phase 4 §1 |
| **MOSAIC** (Campbell & Reyna 2026) | -- | methods cannot separate true exclusivity from subgroup structure | addressed by subtype strata |
| **CDCOCA** (Kumar 2011); **Zack 2013**; **Watkins 2020** | copy-number cohorts | copy-number co-occurrence must be judged against each tumour's copy-number complexity | tested: burden sensitivity check (§8) |
| **GeneAccord** (Kuipers 2021) | single-cell AML | clonal co-occurrence inside a tumour | not used: panels have no clonal reconstruction |

### F. Allele-level, fusion and GENIE-scale studies (closest to this project)
| Paper | Data | Method | Finding | Status here |
|---|---|---|---|---|
| **Scharpf 2022** | 66,372 GENIE v6.1 tumours (27,257 tumour-normal for sequence co-mutations) | Bayesian hierarchical log-linear model pooled across cancer types; TMB-quintile and dominant-signature stratified sensitivity models | RAS co-mutation is lineage- and allele-specific (KRAS G12C lung enriched for NTRK3, chromatin regulators) | ideas borrowed: smoking-signature adjustment (Phase 4 §3); partial pooling across cancer types (§7) |
| **Cook 2021** | 13,492 tumours | per-KRAS-allele co-mutation networks | each allele has its own partners in each tissue | adopted: allele screen (Phase 4 §3) |
| **Vaeyens 2023** | 64,807 cBioPortal samples | exclusivity by variant class | BRAF/EGFR exclusivity depends on variant class | adopted: variant-level lookup |
| **Jiang 2023; Tang 2018; Zhang 2020** | fusion-positive cohorts / review | co-alterations with ALK/RET/ROS1 fusions | TP53 / CDKN2A loss common and prognostic; double-driver cases rare but real | adopted: fusion screen (§8) |

**Scharpf 2022 is the closest prior work** -- same database, same concern with TMB, and panel-aware
denominators like ours (counted only on platforms covering both genes). Differences: it covers RAS
genes only; it handles mutation load by stratifying into TMB quintiles rather than conditioning on
each tumour's own count; it does not adjust for copy-number burden; and it has no lookup tool. It also
benchmarked itself against SELECT and DISCOVER (3,513 lung cancers) and found its TMB-stratified and
pooled models concordant -- the G12C/STK11 and G12C/KEAP1 links are strongest in low-TMB tumours.

---

## 5. Public tools a clinician would use today

| Tool | Paper | What it gives | What it lacks for this question |
|---|---|---|---|
| **cBioPortal** (incl. GENIE's instance) | Cerami 2012; Gao 2013 | interactive counts, Fisher-based co-occurrence | no adjustment for load, subtype or panel testability; one query at a time |
| **OncoKB** | Chakravarty 2017 | evidence level per variant (FDA, guidelines) | single variants, not combinations; licence needed |
| **CIViC** | Griffith 2017 | open, curated clinical evidence, some combination entries | no frequencies |
| **Cancer Hotspots**, **COSMIC** | Chang 2018; Tate 2019 | where recurrent mutations sit and how often | no co-occurrence statistics |

Today, "how common is KRAS G12C + STK11 + KEAP1 in lung cancer, and is that more than chance?"
means querying cBioPortal by hand and reading literature. Phase 3 automates it.

---

## 6. Where methods differ -- and the choice made here

| Question | Options in the literature | Choice here and why |
|---|---|---|
| How to model mutation load | plug-in per-tumour rates (DISCOVER); permutation (WeSME, SELECT); per-sample weights (WExT); TMB-quintile strata + Bayesian pooling (Scharpf) | **conditional (Rasch) test** -- the exact margin-preserving null; the plug-in estimate biases towards exclusivity on panel-sized gene sets |
| Hypermutated tumours | ignore; exclude; TMB strata (Scharpf) | **non-hypermutated patients as the main analysis** (Jason's suggestion); hypermutated (TMB 10-100) separately; ultramutated (>= 100) excluded. Splitting alone is not enough, so the conditional test is kept |
| Subtype | usually ignored | **gene rates fitted within detailed subtypes** (>= 50 patients) |
| Copy-number complexity | ignored by most; CDCOCA / Zack adjust | **sensitivity check, not default**: building it in over-adjusts (§8) |
| Gene vs variant | gene level (most); allele level (Scharpf, Cook, Vaeyens) | gene level in the screen; variant level in the lookup tool and Phase 4 allele screen |
| Fusions | rarely in co-occurrence screens | **added** as a third alteration type with its own testability (§8) |
| What to rank by | p-value (DISCOVER); effect scores (SELECT); selection strength (Coselens) | q-value + fold enrichment; **weak effects flagged** (< 1.5-fold) |
| Validation | none (most GENIE papers); cell lines (Mina 2020); cross-cohort (Park & Lehner) | **out-of-centre replication built in**; **TCGA as external cohort** (§10) |
| Multi-gene sets | MEMo, Dendrix, CoMEt, Mutex | not in the screen; exact combination counts in the tool (§7 compares) |
| Clonality | GeneAccord | not possible on panel data (limitation) |

---

## 7. Head-to-head: each method family run on the same GENIE data

Every family was run on the same patients and pairs (non-hypermutated, panels covering every
gene used, 9 cancer types) and scored the same way: **(1)** false positives on null data that keep
every patient's alteration count and every gene's frequency but destroy gene-gene links, and
**(2)** replication of each method's findings in TCGA (`15_method_benchmark.py`).
All §7 tables are reproduced in `notebooks/05_validation_and_methods.ipynb`.

9 cancer types, 2,553 pairs, identical patients for every method:

| Family / method | null pairs p < 0.05 (target 5%) | null pairs after FDR (target ~0) | significant on real data (q < 0.05) | of those, replicate in TCGA | significant in TCGA, opposite direction |
|---|---|---|---|---|---|
| A Fisher (cBioPortal, maftools) | 11.4% | 0.6% | 990 | 28.6% | 0 |
| B1 DISCOVER-style plug-in (Rediscover) | **32.7%** | **18.1%** | 1,101 | 30.5% | 0 |
| B2 conditional test, no strata | 2.0% | 0% | 853 | 33.8% | 0 |
| **B3 conditional + subtype strata (Phase 2)** | **1.9%** | **0%** | 677 | **37.0%** | 0 |
| D SELECT (Mina 2017; re-implemented, empirical FDR) | 3.8% | 0.2% | 770 | 33.7% | 0 |

Ranking only -- of each method's K strongest pairs testable in TCGA, % replicating:

| K | A Fisher | B1 DISCOVER | B2 conditional | B3 + subtype | D SELECT (effect size) |
|---|---|---|---|---|---|
| 50 | 84% | 92% | 92% | 86% | **96%** |
| 100 | 77% | 79% | 83% | 78% | **84%** |
| 200 | 57% | 64% | 63% | 56% | 61% |

**How SELECT was re-implemented** (from Mina 2017, Methods S1): weighted mutual information with the
per-class sample-count correction (sigmoid, l = 15, r = 0.1; classes = detailed subtypes); an empirical
null from 1,000 randomisations that keep every sample's and every gene's counts within each subtype
(curveball trades; naive column shuffles for blocks with < 30 alterations, as in the paper); effect size
= observed wMI minus the null's upper quartile; Tusher-style empirical FDR per interaction type. Two
steps were not reproduced: randomising each alteration type separately (all methods here see one
gene-level matrix) and the correction for indirect links (a network-building step). A first,
simplified "SELECT-style" score without a null had ranked worst; that score was not SELECT and is
withdrawn.

What it shows:
- **Calibration needs a margin-preserving null.** The conditional test (exact) and SELECT (by
  randomisation) are both calibrated. Fisher over-calls on null data; the DISCOVER-style plug-in is
  worst -- 18% of null pairs survive FDR, all as false exclusivity. (van de Haar 2019 describe Fisher
  as *conservative* in raw counts; both views hold: it makes few calls on their data, but here, where
  the null keeps each tumour's load, its calls are inflated, and in both studies they concentrate on
  low-load genes.)
- **Our calls replicate best** (37% with subtype strata vs 34% for SELECT and 29% for Fisher), with
  fewer calls -- the extra calls of the other methods are mostly the ones that fail.
- **The strongest pairs are found by every method** (84-92% of each top-50 replicate); methods differ
  on the marginal calls, which is where false positives live.
- **Caveat:** subtype strata rank slightly lower at K = 200 (56% vs 63%) because the pairs they remove
  are subtype markers, which TCGA reproduces too (it mixes subtypes as well). External replication
  cannot penalise an artefact the two cohorts share -- the case for subtype strata rests on the
  removed pairs being textbook subtype markers (§8), not on replication.
- **SELECT ranks the very strongest pairs best** (96% of its top 50 replicate): its effect size is a
  good ranking score. Our test's advantage is the precision of the significant list, not the top of
  the ranking.
- WeSME's permutation null is the null the conditional test computes exactly, so it is represented
  by B2 (the null data above are generated by exactly that permutation).

**Mutation-load confounding (van de Haar 2019 diagnostic; `21_mutation_load_confounding.py`).**
Each gene's mutation load association (MLA: z-score of log TMB in a logistic regression of "gene
altered"; our reading of their standardised score, as their STAR Methods are not in the open
manuscript) was correlated with the number of significant exclusive pairs each method calls for it.
van de Haar found r = -0.60 to -0.87 for every method on TCGA/MSKCC exomes. On GENIE panels:

| Patients | Fisher | DISCOVER-style | SELECT | Conditional | Conditional + subtype |
|---|---|---|---|---|---|
| non-hypermutated only (our design) | +0.16 | +0.25 | +0.13 | +0.16 | +0.20 |
| hypermutated (TMB 10-100) included | -0.17 | -0.07 | -- | -0.16 | -0.17 |

- The strong confounding van de Haar report does not appear for any method on GENIE panels.
- A weak version appears once hypermutated tumours are included -- and it affects the conditional test
  as much as Fisher, so the test alone does not remove it.
- Analysing non-hypermutated tumours separately (Prof. Wong's suggestion) removes it for every method:
  direct evidence for that design decision, complementing the calibration result above (the split and
  the conditional test fix different problems).

**Exclusive gene sets (family C).** Dendrix's coverage-minus-overlap score versus a load-aware
version using the conditional model's expected overlap, scored by whether the top sets form one
known pathway (`16_exclusive_gene_sets.py`):

| | Dendrix (coverage - overlap) | load-aware (observed vs expected overlap) |
|---|---|---|
| top-20 sets with all 3 genes in one pathway | 5% | **15%** |
| top-20 sets with >= 2 genes sharing a pathway | 46% | **73%** |
| mean frequency of the genes chosen | 0.24-0.30 | 0.19-0.22 |

Top set per cancer type -- load-aware: BRAF/KRAS/NRAS (colorectal), AKT1/PIK3CA/PIK3R1
(endometrial), CDK4/CDKN2A/RB1 (glioma), EGFR/ERBB2/KRAS (lung) -- textbook pathways. Dendrix:
TP53 bundled with other common genes (ATM/SMAD2/TP53, AKT1/PTEN/TP53, KRAS/MDM2/TP53), because a
coverage score that ignores mutation load rewards whatever is frequent. Same lesson as for pairs:
**the load-aware null matters for sets too.** (6 cancer types, 35 most altered genes each, all
triples scored.)

**Conditional selection (family D, Coselens idea).** Coselens needs genome-wide passenger
mutations, which panels lack, but panels do report silent mutations in the genes they cover
(85,959 in non-hypermutated patients). Borrowed idea: is gene B's ratio of protein-changing to
silent mutations higher (co-selection) or lower (selection against) when gene A is mutated?
Silent mutations measure each group's background mutation rate, so this controls mutation load
by a completely different route than the conditional test (`17_conditional_selection_dnds.py`):
- 219 testable directions among mutation-driven replicated candidates;
- direction agrees with Phase 2 in **84%**; among the 53 that are themselves significant, **100%**
  agree -- e.g. IDH1-mutant gliomas carry 8x more protein-changing (relative to silent) CIC
  mutations than IDH1-wild-type gliomas (oligodendroglioma biology).
- Limitation: selection-based checks share subtype confounding (KRAS vs MEN1 in pancreatic cancer
  appears in both, partly the adenocarcinoma/neuroendocrine split).

**Partial pooling across cancer types (family F, Scharpf idea).** Each cancer type's effect is
pulled towards the pair's pooled effect across types, more when it is noisy (empirical-Bayes,
DerSimonian-Laird; `18_hierarchical_shrinkage.py`) -- an approximation of Scharpf's fully Bayesian
hierarchical model (fitted in Stan). Against TCGA effects:
- all 2,152 testable pairs: correlation 0.691 -> 0.700, squared error 0.206 -> 0.195;
- noisy estimates (< 30 co-altered patients): correlation 0.558 -> 0.577, error -7%;
- direction agreement unchanged (~97%).
A modest, consistent gain -> pooled estimates kept as an extra column for small cancer types.

---

## 8. Tests that shaped the design (GENIE)

| Test | Result | Script |
|---|---|---|
| Fisher vs conditional on null data (non-hypermutated) | Fisher 6-15% false positives at p < 0.05 (up to 1.8% of pairs after FDR); conditional 1.3-2.9%, 0% after FDR | `04_null_calibration_nonhyper.py` |
| TMB split | splitting cuts Fisher's calls sharply (colorectal 85% -> 8%) but Fisher stays above the conditional test (lung 31% vs 24%) | `02_tmb_split_fisher_vs_conditional.py` |
| Ultramutated tumours | endometrial: 21.6% significant pooled; 0.3% (TMB 10-100) and 0% (>= 100) when separated | `03_ultramutated_split.py` |
| Subtype-specific gene rates | removes 390 of 1,571 hits in 10 types; removed ones are textbook subtype markers (CIC vs PTEN in glioma, CDH1 vs ERBB2 in breast, PIK3CA vs TP53 in ovarian, MEN1 vs SMAD4 in pancreas) | `05_subtype_and_replication.py` |
| Replication MSK -> other centres | 96% keep direction, none flips; subtype-robust hits replicate 73%, others 35% | `05_subtype_and_replication.py` |
| Whole-design comparison (replication) | all patients 44%; non-hypermutated 66%; + subtype 66% (artefacts removed); mutation-only 84% (fewer pairs) | `06_design_comparison.py` |
| Filters | dropping the gene-frequency filter: ~10x pairs, 1-20 h per cancer type; the joint-coverage filter removes only single-centre pairs | `p2_variants.py`, `07_joint_coverage_filter.py` |
| Copy-number burden (CDCOCA / Zack) | per-pair logistic re-test: 21 of 25 copy-number-driven candidates survive, 4 are artefacts (BCL2L1/MYC colorectal, CCND1/CDKN2A bladder, CDK4/TERT and FOXA1/TERT lung). In the final run 1,506 of 1,824 significant pairs (1,276 of 1,584 candidates) pass the burden sensitivity check. As default strata: counting all genes over-adjusts (CCND1/FGFR1 lost); counting genes outside the analysis does not improve replication (64% vs 66%) -> **sensitivity column `burden_robust`, not default** | `12_burden_adjusted_retest.py`, `13_burden_strata.py` |
| Fusions | 9,619 driver-fusion calls (Phase 1; 35 panels actually release fusions, while 9 panels claim fusion detection but release none). Fusion screen: 255 pairs in 13 cancer types, 109 significant (92 exclusive). Recovers known biology: ALK/ROS1/RET fusions avoid EGFR, KRAS, STK11, KEAP1 in lung cancer (ALK + EGFR: 4 vs 140 expected); ERG fusions avoid SPOP, FOXA1, CDK12 and co-occur with PTEN loss in prostate cancer; FGFR2 fusions avoid IDH1/KRAS and co-occur with BAP1 in cholangiocarcinoma; ALK fusions co-occur with CDKN2A/B loss (Jiang 2023). Note frequency vs enrichment: TP53 is present in 24% of ALK-fusion lung cancers (Jiang 2023 reports ~31% in fusion-positive NSCLC) yet is *under*-represented relative to expectation (124 vs 191) -- ALK fusions arise mostly in never-smokers. Artefact flagged: breast 'FGFR1 fusions' co-occur with FGFR1 amplification -- rearrangements inside the 8p11 amplicon, not driver fusions | Phase 1 §4c, Phase 2 §6e |
| Mutation reporting vs panel lists | 1,144 (panel, gene) pairs report zero mutations where >= 10 are expected; 8 panels release no mutations | `01_mutation_panel_reporting.py` |

**Phases 3-4 methods** (patient matching, therapy evidence, pathway, heterogeneity, allele and
smoking tests) are summarised with their results in §10.

---

## 9. Opposing views, and where this project sits

| Debate | One side | Other side | What our data say |
|---|---|---|---|
| Is co-occurrence biology or chance? | mostly chance once heterogeneity is modelled (Canisius 2016) | co-selection is widespread (Mina 2020: 852 dependencies in >9,000 tumours, 44% co-occurring, robust to TMB; Iranzo 2022) | both: exclusive candidates outnumber co-occurring ~3:1 and are pathway-structured; some co-occurrences are strong, mutation-driven, replicated in TCGA and supported by conditional selection (KEAP1/STK11, IDH1/CIC) |
| What does exclusivity mean? | pathway redundancy (Sanchez-Vega 2018) | artefact of subtype and load (van de Haar 2019); incompatibility (El Tekle 2021); functional redundancy shown in cell-line knockouts (Mina 2020) | all three occur: 3-4x enrichment within pathways; subtype strata remove a quarter of hits (van de Haar attribute 58% of colorectal hits to CMS subtypes); cross-pathway exclusive pairs (EGFR vs STK11/KEAP1) fit incompatibility |
| Universal or tissue-specific? | core interactions recur pan-cancer | interactions rewire between tissues (Park & Lehner 2015) | 30% of comparable pairs differ between cancer types, 47 reverse direction; 21 universal |
| Gene or variant level? | gene-level screens | allele-specific partners (Cook 2021; Scharpf 2022) | allele differences real but few (10 of 201) and partly explained by smoking |
| How to handle hypermutation? | TMB strata (Scharpf) or condition on load (DISCOVER family) | split the cohort (Jason) | both needed |
| Copy-number co-occurrence | count overlapping events | judge against complexity (Kumar 2011; Zack 2013) | complexity matters (4 of 25 artefacts) but naive adjustment over-adjusts |
| Is anything left to discover at the driver-pair level? | large screens keep reporting new pairs | most significant exclusivity calls involve a few low-load driver genes and are likely confounded (van de Haar 2019) | of our 386 most robust pairs, 41 are in neither the literature nor the SelectSim atlas, and those are amplicon passengers or subtype markers; the strongest candidate (CTCF vs TP53, endometrial) is in SelectSim's atlas; all pairs SelectSim reports agree in direction (notebook 04 §9) |
| Significance or effect size? | p/q-values | effect scores (SELECT) | with >10,000 patients tiny effects are significant; 22% of candidates (350 of 1,584) < 1.5-fold, flagged |

**The gap this project fills:**
1. **Panel-aware statistics at GENIE scale.** The rigorous methods (family B) were built on TCGA;
   GENIE disease papers mostly use family A. This project brings the rigorous method to GENIE and
   tracks testability per patient for mutations, copy number and fusions -- read from the data,
   because the panel metadata is wrong for all three.
2. **Confounders handled together.** Load, hypermutation, subtype, one-DNA-event artefacts,
   copy-number complexity and centre effects -- each tested before adoption, rejected options
   documented.
3. **Validation built in.** Out-of-centre replication for every pair, plus an independent cohort
   (TCGA) and an independent route (conditional selection from silent mutations).
4. **A clinical lookup tool** combining variant-level frequencies, exact combination counts, an
   on-demand test, partner prediction and CIViC evidence -- existing tools offer one of these each.
5. **A novelty audit, not just a hit list.** Every robust pair is checked against the literature
   and the SelectSim atlas (`19_novelty_screen.py`) and the unreported ones are stress-tested for proxies (`20_candidate_followup.py`),
   so claims of "new" biology are separated from passengers and subtype markers.

---

## 10. Validation and the later phases

**External validation in TCGA** (`14_tcga_external_validation.py`). The same pairs re-tested with
the same test in TCGA PanCanAtlas (independent patients, whole exome, non-hypermutated, rules of
Phase 1 applied as closely as possible):

| GENIE result | testable in TCGA | same direction | TCGA p < 0.05, same direction | significant opposite |
|---|---|---|---|---|
| replicated candidates | 316 | 94% | 46% | 0 |
| all candidates | 765 | 88% | 32% | 0 |
| not significant in GENIE | 1,317 | 69% | 7% (≈ chance) | 13 |

GENIE's replicated candidates validate 6.5x more often than GENIE-null pairs, none reverses, and
KEAP1/STK11, EGFR/KRAS, CCND1/FGFR1 and CDKN2A/EGFR are all confirmed. TCGA is ~20x smaller, so a
non-significant TCGA result mostly reflects low power. Concordance is lowest in lung cancer
(79% same direction), where TCGA is half squamous and heavily smoking-related while GENIE is
mostly adenocarcinoma.

**Phases 3 and 4: methods tried, and what won**

| Question | Options tried (from the literature) | Result | Adopted |
|---|---|---|---|
| Match a new patient to similar cases (registered aim) | nearest neighbours (Jaccard; rarity-weighted Jaccard as in Foote 2025 / TumorComparer; + subtype); summed pairwise ratios; per-gene logistic models | held-out AUROC: base rate 0.82, similarity 0.81-0.82, summed ratios 0.77, **logistic 0.855** (better in every cancer type, calibrated) | logistic (`predict_partners`); exact matches via `combo_prevalence` |
| Therapy evidence | OncoKB (licence); CIViC (open) | 3,433 accepted predictive/prognostic items incl. combinations | CIViC |
| Rare pairs below Phase 2's filters | counts only; on-demand test | on-demand test reproduces Phase 2 (KRAS/STK11 expected 460.5 vs 460.4) | `test_pair` |
| Are frequencies reproducible? | out-of-centre comparison | correlation 0.89-0.996 | reported |
| Pathway structure (Sanchez-Vega; Mina) | within / across the 10 pathways | exclusivity 3.3x (replicated 4.5x) within pathways; co-occurrence not enriched across (OR ~1) | Phase 4 §1 |
| Cancer-type specificity (Park & Lehner) | Cochran's Q across types | 30% heterogeneous; 47 reverse; 21 universal | Phase 4 §2 |
| Allele-specific partners (Scharpf; Cook; Vaeyens) | allele groups + Cochran's Q | 10 of 201 differ (EGFR L858R vs ex19 for RBM10, q = 4e-22) | Phase 4 §3 |
| Smoking confounding of KRAS G12C | signature-aware (Scharpf) | G12C/STK11 2.00 -> 1.74; G12C/RBM10 explained entirely; G12D unchanged | Phase 4 §3 |

---

## References

**Methods -- pairwise and load-aware tests**
1. Canisius S, Martens JWM, Wessels LFA. A novel independence test for somatic alterations in cancer shows that biology drives mutual exclusivity but chance explains most co-occurrence. *Genome Biol* 2016;17:261. doi:10.1186/s13059-016-1114-x
2. Ferrer-Bonsoms JA, Jareno L, Rubio A. Rediscover: an R package to identify mutually exclusive mutations. *Bioinformatics* 2022;38:844-5. doi:10.1093/bioinformatics/btab709
3. Kim YA, Madan S, Przytycka TM. WeSME: uncovering mutual exclusivity of cancer drivers and beyond. *Bioinformatics* 2017;33:814-21. doi:10.1093/bioinformatics/btw242
4. Leiserson MDM, Reyna MA, Raphael BJ. A weighted exact test for mutually exclusive mutations in cancer. *Bioinformatics* 2016;32:i736-45. doi:10.1093/bioinformatics/btw462
5. Mayakonda A et al. Maftools: efficient and comprehensive analysis of somatic variants in cancer. *Genome Res* 2018;28:1747-56. doi:10.1101/gr.239244.118

**Methods -- gene sets**
6. Ciriello G, Cerami E, Sander C, Schultz N. Mutual exclusivity analysis identifies oncogenic network modules. *Genome Res* 2012;22:398-406. doi:10.1101/gr.125567.111
7. Vandin F, Upfal E, Raphael BJ. De novo discovery of mutated driver pathways in cancer. *Genome Res* 2012;22:375-85. doi:10.1101/gr.120477.111
8. Leiserson MDM, Blokh D, Sharan R, Raphael BJ. Simultaneous identification of multiple driver pathways in cancer. *PLoS Comput Biol* 2013;9:e1003054. doi:10.1371/journal.pcbi.1003054
9. Leiserson MDM, Wu HT, Vandin F, Raphael BJ. CoMEt: a statistical approach to identify combinations of mutually exclusive alterations in cancer. *Genome Biol* 2015;16:160. doi:10.1186/s13059-015-0700-7
10. Babur Ö et al. Systematic identification of cancer driving signaling pathways based on mutual exclusivity of genomic alterations. *Genome Biol* 2015;16:45. doi:10.1186/s13059-015-0612-6

**Methods -- selection and evolution**
11. Mina M et al. Conditional selection of genomic alterations dictates cancer evolution and oncogenic dependencies. *Cancer Cell* 2017;32:155-68. doi:10.1016/j.ccell.2017.06.010
12. Mina M et al. Discovering functional evolutionary dependencies in human cancers. *Nat Genet* 2020;52:1198-207. doi:10.1038/s41588-020-0703-5
13. Iranzo J, Gruenhagen G, Calle-Espinosa J, Koonin EV. Pervasive conditional selection of driver mutations and modular epistasis networks in cancer. *Cell Rep* 2022;40:111272. doi:10.1016/j.celrep.2022.111272
14. Kuipers J et al. Statistical tests for intra-tumour clonal co-occurrence and exclusivity. *PLoS Comput Biol* 2021;17:e1009036. doi:10.1371/journal.pcbi.1009036

**Confounders and context**
15. van de Haar J et al. Identifying epistasis in cancer genomes: a delicate affair. *Cell* 2019;177:1375-83. doi:10.1016/j.cell.2019.05.005
16. Park S, Lehner B. Cancer type-dependent genetic interactions between cancer driver alterations indicate plasticity of epistasis across cell types. *Mol Syst Biol* 2015;11:824. doi:10.15252/msb.20156102
17. Sanchez-Vega F et al. Oncogenic signaling pathways in The Cancer Genome Atlas. *Cell* 2018;173:321-37. doi:10.1016/j.cell.2018.03.035
18. Campbell K, Reyna MA. MOSAIC: model-based, subgroup-aware identification of driver mutations in cancer. *bioRxiv* 2026. doi:10.64898/2026.04.29.721672
19. Kumar N, Rehrauer H, Cai H, Baudis M. CDCOCA: a statistical method to define complexity dependence of co-occuring chromosomal aberrations. *BMC Med Genomics* 2011;4:21. doi:10.1186/1755-8794-4-21
20. Zack TI et al. Pan-cancer patterns of somatic copy number alteration. *Nat Genet* 2013;45:1134-40. doi:10.1038/ng.2760
21. Watkins TBK et al. Pervasive chromosomal instability and karyotype order in tumour evolution. *Nature* 2020;587:126-32. doi:10.1038/s41586-020-2698-6
22. El Tekle G et al. Co-occurrence and mutual exclusivity: what cross-cancer mutation patterns can tell us. *Trends Cancer* 2021;7:823-36. doi:10.1016/j.trecan.2021.04.009

**Allele-level, fusion and GENIE studies**
23. Scharpf RB et al. Genomic landscapes and hallmarks of mutant RAS in human cancers. *Cancer Res* 2022;82:4058-78. doi:10.1158/0008-5472.CAN-22-1731
24. Cook JH et al. The origins and genetic interactions of KRAS mutations are allele- and tissue-specific. *Nat Commun* 2021;12:1808. doi:10.1038/s41467-021-22125-z
25. Vaeyens F et al. Variant-specific landscape of mutual exclusivity among BRAF, EGFR, and KRAS oncogenes in human cancer. *medRxiv* 2023. doi:10.1101/2023.10.21.23297089
26. Jiang B et al. TP53 or CDKN2A/B covariation in ALK/RET/ROS1-rearranged NSCLC is associated with a high TMB, tumor immunosuppressive microenvironment and poor prognosis. *J Cancer Res Clin Oncol* 2023;149:10041-52. doi:10.1007/s00432-023-04924-7
27. Tang Z et al. Coexistent genetic alterations involving ALK, RET, ROS1 or MET in 15 cases of lung adenocarcinoma. *Mod Pathol* 2018;31:307-12. doi:10.1038/modpathol.2017.109
28. Zhang R, Dong L, Yu J. Concomitant pathogenic mutations and fusions of driver oncogenes in tumors. *Front Oncol* 2020;10:544579. doi:10.3389/fonc.2020.544579
29. Skoulidis F, Heymach JV. Co-occurring genomic alterations in non-small-cell lung cancer biology and therapy. *Nat Rev Cancer* 2019;19:495-509. doi:10.1038/s41568-019-0179-8

**Data sources**
30. AACR Project GENIE Consortium. AACR Project GENIE: powering precision medicine through an international consortium. *Cancer Discov* 2017;7:818-31. doi:10.1158/2159-8290.CD-17-0151
31. Hoadley KA et al. Cell-of-origin patterns dominate the molecular classification of 10,000 tumors from 33 types of cancer. *Cell* 2018;173:291-304. doi:10.1016/j.cell.2018.03.022
32. ICGC/TCGA Pan-Cancer Analysis of Whole Genomes Consortium. Pan-cancer analysis of whole genomes. *Nature* 2020;578:82-93. doi:10.1038/s41586-020-1969-6
33. Zehir A et al. Mutational landscape of metastatic cancer revealed from prospective clinical sequencing of 10,000 patients. *Nat Med* 2017;23:703-13. doi:10.1038/nm.4333

**Tools and knowledge bases**
34. Cerami E et al. The cBio cancer genomics portal: an open platform for exploring multidimensional cancer genomics data. *Cancer Discov* 2012;2:401-4. doi:10.1158/2159-8290.CD-12-0095
35. Gao J et al. Integrative analysis of complex cancer genomics and clinical profiles using the cBioPortal. *Sci Signal* 2013;6:pl1. doi:10.1126/scisignal.2004088
36. Chakravarty D et al. OncoKB: a precision oncology knowledge base. *JCO Precis Oncol* 2017;2017:PO.17.00011. doi:10.1200/PO.17.00011
37. Griffith M et al. CIViC is a community knowledgebase for expert crowdsourcing the clinical interpretation of variants in cancer. *Nat Genet* 2017;49:170-4. doi:10.1038/ng.3774
38. Chang MT et al. Accelerating discovery of functional mutant alleles in cancer. *Cancer Discov* 2018;8:174-83. doi:10.1158/2159-8290.CD-17-0321
39. Tate JG et al. COSMIC: the Catalogue Of Somatic Mutations In Cancer. *Nucleic Acids Res* 2019;47:D941-7. doi:10.1093/nar/gky1015

**Patient similarity**
40. Foote MB et al. Analysis of shared variants between cancer biospecimens. *Clin Cancer Res* 2025;31:376-86. doi:10.1158/1078-0432.CCR-24-1583
41. Sinha R, Schultz N, Sander C. Comparing cancer cell lines and tumor samples by genomic profiles (TumorComparer). *bioRxiv* 2015. doi:10.1101/028159

**Closest recent atlas (added 25 Oct 2026)**
42. Iyer A, Petrovic M, Sesia D, Nanni L, Mina M, Ciriello G. Evolving patterns of co-mutations from tumor initiation to metastatic progression. *Nat Genet* 2026. doi:10.1038/s41588-026-02661-4
