# Learning Cancer Mutation Landscapes from Large-Scale Somatic Mutation Testing Data

Final Year Project (HKU, B.Sc. Bioinformatics) analysing the [AACR Project
GENIE](https://www.aacr.org/professionals/research/aacr-project-genie/)
public dataset (v19.0, ~270,000 clinical tumour sequencing samples) to find
which gene alterations occur together, or avoid each other, within and across
cancer types -- and to turn that into a lookup tool for the clinic.

**Supervisor:** Prof. Jason Wong

## Project aims

Refocused with Jason on 2026-09-24 (supersedes the original cluster-and-survival brief):

1. **A co-mutation lookup tool.** Clinicians meet rare combinations of
   alterations and currently need someone to search GENIE by hand. The tool
   takes a cancer type and a patient's known alterations and returns what else
   is typically altered alongside them, how common the exact combination is,
   and (planned) the therapy evidence for it. Co-occurrence is the focus;
   exclusivity is a bonus insight. This is also the registered aim of matching
   a new patient to similar cases and reporting how common the combination is.
2. **Scientific insight.** Which combinations are unusually common or rare in a
   cancer type, which are specific to one cancer type or to hypermutated
   tumours, and why.

Why the refocus: Phase 2's gates (gene altered in >= 3% of patients, >= 5 expected
co-altered) remove exactly the rare combinations clinicians ask about (for *EGFR*
in lung cancer, 651 partner genes occur in >= 3 patients but only 28 are testable
pairs), so the tool reads the Phase 1 tables directly and uses Phase 2 as
annotation. GENIE survival data is coarse (year-level, a subset of centres), so
multi-gene clusters and survival analysis were dropped to optional.

### Registered project description

*Official title:* **Learning Cancer Mutation Landscapes from Large-Scale Somatic
Mutation Testing Data.** Analyse a database of over 200,000 clinical tumour
sequencing samples to investigate patterns of somatic mutations across cancer
types. *Registered Aim 1:* identify genes whose mutations co-occur more
frequently than expected within and across cancers. *Registered Aim 2:* develop
a computational method that matches a new patient's mutation profile to
genetically similar cases and reports the prevalence of the observed
combination (delivered by the Phase 3 lookup tool; similarity matching itself
was tested and rejected, see Phase 3). Learning outcomes: apply statistical and
bioinformatics methods to large genomic datasets; develop and evaluate
computational approaches; interpret findings in precision oncology; critically
evaluate methods and communicate strengths and limitations.

## Repository structure

Organised along the usual lines for computational biology projects (Noble 2009; Wilson et al.
2017): raw data are read-only, every derived file is produced by code in the repository, and each
folder has its own README.

```
README.md, LICENSE        This file (aims, pipeline, results, data problems); MIT licence for the code
requirements.txt, setup.sh  Python environment (R packages: scripts/figures/install_r_packages.R)

notebooks/                Main pipeline, run in numeric order
  01_phase1_preprocessing.ipynb      Phase 1: clean the data, record what was tested
  02_phase2_comutation_matrix.ipynb  Phase 2: pairwise co-occurrence / exclusivity tests
  03_phase3_lookup_tool.ipynb        Phase 3: the lookup tool (Aim 1)
  04_phase4_insight.ipynb            Phase 4: what the patterns tell us (Aim 2), incl. novelty audit
  05_validation_and_methods.ipynb    TCGA validation, method benchmark, borrowed methods
  06_phase2_deep_dive.ipynb          Trust checks on Phase 2 hits (subtype, burden, proximity, passengers)

scripts/                  See scripts/README.md (index of every script and the run order)
  analyses/    One question per script (01-20), behind each design decision
  validation/  Checks of the Phase 2 test itself: simulation, exact tail, Rediscover (R)
  figures/     R publication figures; run_all.sh rebuilds all of them
  fetch_external.py, build_notebook05.py

data/                     See data/README.md (sources and a guide to every processed file)
  raw/         GENIE v19.0, unmodified (not in git -- see Data access)
  external/    Annotation resources and external cohorts (not in git)
  processed/   Derived tables produced by the notebooks and scripts (in git)

results/
  figures/     Publication figures, PDF + 300-dpi PNG (in git)
  analyses/, validation/, logs/   Script outputs and logs (not in git; regenerate by running the scripts)

docs/
  FYP_onboarding_guide.docx  End-to-end walkthrough for someone new to the project
  supervisor_summary.md      Progress, results and open decisions for supervisor review
  literature_review.md       Published methods, how each performed on our data, what we adopt
```

Earlier exploratory notebooks (old clusters, survival, SNV/CNV exploration)
were removed; they are in the git history (last present in commit `32a6c97`).
The original project plan (`PLAN.md`: clusters, survival, Fisher's test, a
data-file inventory) was deleted once it no longer matched the work; it is in
the git history (last present in commit `92e8234`).

Every notebook is self-contained, starts with a "data flow at a glance"
diagram, and ends with a Discussion section.

## Data access

GENIE data cannot be redistributed. `data/raw/` and `data/external/` are
gitignored; only the derived, aggregate outputs in `data/processed/` are
tracked. To reproduce this pipeline, request GENIE v19.0 access via
[AACR Project GENIE](https://www.aacr.org/professionals/research/aacr-project-genie/)
and place the files under `data/raw/`.

## Setup

```
cd fyp
./setup.sh
```

Creates a Python virtual environment, installs dependencies from
`requirements.txt`, and registers a Jupyter kernel named
**"Python (fyp-genie)"** — select it before running any notebook.

The figures and the Rediscover comparison need R (tested with 4.5.1):
`Rscript scripts/figures/install_r_packages.R`. The full run order is in
`scripts/README.md`.

## Status

- **Phase 1 and Phase 2 are current** (re-run end to end, no errors). Phase 2 was
  redesigned after a literature review and head-to-head tests
  ([`docs/literature_review.md`](docs/literature_review.md)): non-hypermutated
  patients as the main analysis, gene rates per detailed subtype, and
  out-of-centre replication built in. Gene fusions added as a third alteration
  type; copy-number burden checked as a sensitivity analysis.
- **Validated externally** in TCGA (independent patients, whole exome) and
  benchmarked against every published method family on the same data
  (literature review §7, §10).
- **Phase 3 (lookup tool)** and **Phase 4 (Aim 2 insight)** are built and run;
  their methods were chosen by testing alternatives (literature review §7).
- **Known data limit:** some panels report mutations in far fewer genes than
  they list; this can only be detected statistically, so it is documented and
  left as is (the directly checkable parts are fixed -- see *Problems found in
  the GENIE data itself*).
- **Not yet supervisor-reviewed:** the conditional pair test, the germline
  cutoff change (gnomAD 1e-4 -> 1e-3), one sample per patient, TMB-based
  hypermutation status, the expected-count gate, the shared-DNA filter, and
  the v3 design (non-hypermutated main analysis, subtype strata, replication),
  fusions, and the burden sensitivity check. One-page summary for review:
  [`docs/supervisor_summary.md`](docs/supervisor_summary.md).

### Headline numbers

| | |
|---|---|
| Samples / patients | 271,837 samples, 167 panels -> **227,696 patients** (one sample each) |
| Samples with copy-number data | 172,874 (63.6%) |
| Mutations kept (likely harmful) | 3,458,550 -> **1,644,226** |
| Deep copy-number calls | 377,216 -> 291,962 after dropping direction-inconsistent bystanders |
| Driver gene fusions | 56,356 structural variants -> **9,619** driver-fusion calls (ALK, ROS1, RET, ERG, FGFR2 ...) |
| Combined alteration table | 1,945,807 rows (232,386 samples): mutations, copy number, fusions |
| Phase 2 main analysis | non-hypermutated patients: 6,697 gene pairs (227 genes), 51 cancer types, 112,779 patients |
| Significant (q < 0.05) | 1,824 (27.2%) |
| ...and not one copy-number event hitting neighbouring genes | **1,584 interaction candidates** -- 426 co-occurring, 1,158 exclusive; 350 weak (< 1.5-fold) |
| ...also replicated outside MSK (significant in MSK, confirmed in other centres) | **486**; of all 839 MSK hits testable elsewhere, 65% replicate, none in the opposite direction |
| ...survive the copy-number burden sensitivity check | 1,276 of 1,584 candidates |
| External validation in TCGA (replicated candidates) | 94% same direction, 46% significant (vs 7% for GENIE-null pairs), **0 reversed** |
| Fusion screen | 255 pairs, 13 cancer types, 109 significant (e.g. ALK/ROS1/RET fusions vs EGFR/KRAS/STK11) |
| Hypermutated patients on their own (TMB 10-100) | 881 of 98,907 pairs significant (0.9%), 18 cancer types |
| Mutation-only analysis (all patients with mutation data) | 2,403 pairs, 57 cancer types, 155,020 patients, 573 significant |
| Mechanism-specific tests (mutation / copy number per side) | 8,445 tests, 2,239 significant, 687 cross-mechanism |
| Pan-cancer pairs (>= 5 cancer types, same direction) | 13 |

Jason's question -- *"half being significant sounds quite high"* -- was right
about the method: Fisher's test called 72.3% of pairs significant. The share is
now 27.2% of a much smaller, driver-only pair list (hypermutation no longer
pushes passenger genes into the analysis), and only 1,584 pass as candidates.

## How the pipeline works

### Phase 1 -- clean the data (no analysis)

Turns the raw GENIE files into one reliable table of *which gene is altered,
how, in which patient* -- and records which genes each patient was actually
tested for.

1. **Patients.** Join the sample and patient files (271,837 samples). Keep
   **one sample per patient** (227,696): a primary tumour and its metastasis
   share their early mutations, so counting both would count the same pair
   twice. The sample kept is the one with copy-number data, then the primary
   tumour, then the larger panel.
   Each patient also gets a **TMB group** from GENIE's official TMB, trusted
   only on panels that sequence >= 1 Mb: on smaller panels ~61% of samples come
   out as "TMB >= 10", which is noise, not biology.
2. **What was tested.** 167 different gene panels are used across hospitals,
   so "no mutation reported" can mean *tested and clean* or *never tested*.
   Two tables record this separately: genes each panel tests for mutations
   (exon probes in `genomic_information.txt`) and genes with copy-number values
   (read from `data_CNA.txt` itself, because the panel metadata turned out to
   be wrong -- see below). Per patient, `CNA_TESTED` and `MUT_TESTED` record
   whether any copy-number / mutation data was released at all. Every frequency downstream divides only by patients
   actually tested for the gene.
3. **Mutations: keep the likely harmful ones.** 3.46M calls -> 1,644,226.
   - somatic only, and protein-changing only (silent / intron / UTR dropped);
   - rare in healthy people: gnomAD max frequency <= 1e-3 (a 1e-4 cutoff was
     tried first, but it removed ~104k real cancer mutations -- see below);
   - missense changes need evidence of harm: `AlphaMissense pathogenic OR
     (PolyPhen damaging AND SIFT deleterious) OR a cancerhotspots.org hotspot`
     (rule approved by Jason). AlphaMissense is newer and better calibrated
     than PolyPhen/SIFT; hotspots rescue well-known recurrent drivers.
4. **Copy number: deep calls only (+2 amplification, -2 deep deletion).** Why
   not the shallow +1 / -1 calls:
   - a one-copy gain or loss is usually a passenger -- large chunks of a
     chromosome arm gained or lost together, sweeping up many genes. There are
     2.58M shallow calls against 0.38M deep ones;
   - one lost copy rarely switches a tumour suppressor off by itself (the other
     copy still works), and one extra copy rarely activates an oncogene;
   - shallow calls depend heavily on tumour purity and on each centre's
     pipeline, and are reported **inconsistently between centres**: MSK panels
     report +/-1 in only ~5% of samples, DFCI panels in 98-99% (median 25-69
     per sample). Including them would make "altered" mean different things
     depending on the hospital;
   - Jason's point: a shallow deletion is not a signal on its own, but a
     shallow deletion **plus** a damaging mutation in the same tumour
     suppressor suggests both copies are gone (two-hit). That is kept as a
     separate feature (`putative_biallelic_tsg_events.parquet`), not in the
     main table.

   Deep calls whose direction contradicts the gene's role (an oncogene deleted,
   a tumour suppressor amplified -- usually a bystander next to the real
   target, e.g. *RAD21* on the *MYC* amplicon) are dropped using OncoKB gene
   roles: 377,216 -> 291,962. Genes OncoKB does not review are kept.
5. **Gene fusions.** `data_sv.txt` holds 56,356 structural variants; only events
   joining two *different* genes where one is a known fusion driver (OncoKB /
   Mitelman-recurrent list, antisense events excluded) are kept: **9,619**
   driver-fusion calls, named after the driver (e.g. `ALK_fusion`). A sample
   counts as fusion-tested only if its panel actually released structural
   variants (35 panels do; read from the data, `SV_TESTED`).
6. **Output.** `alterations_long.parquet` (1,945,807 rows, 232,386 samples:
   1,644,226 mutations + 291,962 deep copy-number calls + 9,619 fusions),
   `clinical_tidy`, `panel_gene_coverage`, `cna_gene_panel_coverage` and
   `sv_gene_panel_coverage`. Phase 1 does **not** choose genes -- every
   analysis filter lives in Phase 2.

### Phase 2 -- find gene pairs that go together or avoid each other

Run separately for each cancer type. All four filters are defined in one
settings cell and nowhere else.

0. **Who.** Non-hypermutated patients (TMB < 10 on a panel >= 1 Mb): 159,390
   patients. Jason's suggestion, and the design whose findings replicate best
   across centres (44% -> 66%). Hypermutated patients are analysed separately
   (step 8).
1. **Filter 1** -- cancer types with >= 100 patients with both mutation and
   copy-number data: **51 types with testable pairs, 112,779 patients**.
   "Altered" means mutation *or* copy-number change, so "not altered" must be
   checkable for both.
2. **Filter 2** -- genes altered in >= 3% (and >= 5) of that type's patients.
   Counted over all patients, so it also requires the gene to be widely tested;
   counting only tested patients was tried and adds pairs that replicate worse.
3. **Filter 3** -- a pair is tested only if >= 50% of patients were tested for
   both genes. It rarely binds; the pairs it removes rest on one centre's panel
   (only 1 of them could be re-tested outside MSK).
4. **Filter 4** -- and only if >= 5 patients are *expected* to carry both
   (decided from the gene frequencies, before looking at the result).
5. **The test.** Each patient is compared with their **own number of
   alterations**: a tumour with 30 alterations will carry almost any pair by
   chance, one with 2 will not. (Conditional / Rasch model, exact
   Poisson-Binomial tail, then Benjamini-Hochberg FDR within each cancer type.)
   **6,697 pairs -> 1,824 significant (27.2%).**
6. **Subtypes.** Gene rates are fitted within each detailed subtype (>= 50
   patients): broad cancer types mix subtypes with different drivers, which
   otherwise looks like exclusivity -- pancreatic *MEN1* vs *SMAD4*
   (neuroendocrine vs adenocarcinoma), glioma *CIC* vs *PTEN*, breast *CDH1*
   vs *ERBB2* all disappear. Hits that survive this replicate 73% of the time,
   hits that do not 35%.
7. **One-DNA-event check.** Removed: pairs where one copy-number event hit both
   genes (a gained chromosome arm, or neighbours on one amplicon). Neighbouring
   genes (<= 10 Mb) are removed only when copy number drives the pair -- two
   point mutations are two separate events (Jason), so *KEAP1* / *STK11* /
   *SMARCA4* (19p13, lung) and *PBRM1* / *SETD2* (3p21, kidney) stay.
   **1,584 interaction candidates**; `weak_effect` flags the 350 within
   1.5-fold of chance (with 10,000+ patients, tiny effects are significant).
8. **Validation and side analyses.**
   - *Out-of-centre replication:* every pair re-tested in MSK patients and in
     all other centres; 486 candidates are `replicated`.
   - *Copy-number burden sensitivity check* (CDCOCA, Kumar 2011; Zack 2013):
     each pair re-tested with gene rates fitted within groups of copy-number
     burden; `burden_robust` -- 1,276 of 1,584 candidates pass. Built into the
     main test it over-adjusts (loses the known CCND1/FGFR1 co-amplification),
     so it is a flag, not the default.
   - *Fusion screen*: each fusion with >= 20 carriers tested against its cancer
     type's genes -- 109 significant pairs, e.g. ALK/ROS1/RET fusions avoid
     EGFR, KRAS, STK11 in lung cancer; ERG fusions avoid SPOP/FOXA1 and go with
     PTEN loss in prostate cancer.
   - *Hypermutated patients on their own* (TMB 10-100): ultramutated tumours
     (TMB >= 100, mostly *POLE*) are excluded -- pooled with the rest they made
     21.6% of endometrial pairs "significant", 0.3% and 0% when separated.
     881 of 98,907 pairs significant; much of it is subtype structure (MSI vs
     MSS), so it is reported as a description, not as candidates.
   - *Mutation-only analysis* on every patient with mutation data (155,020):
     mutation calls are the most consistent between centres (replication 84%).

**Why both the TMB split and the conditional test.** Splitting alone is not
enough: on null data built from real non-hypermutated patients (no true
interactions), Fisher's test still calls 6-15% of pairs significant at p < 0.05
(should be 5%) against 1.3-2.9% for the conditional test, because copy-number
burden still varies and TMB does not measure it.
The conditional test alone is not enough either -- on all patients it leaves
25% of colorectal pairs significant vs 5.8% in non-hypermutated patients.

**Why not Fisher's test or the earlier rate-adjusted test.** Fisher's test
assumes every patient has the same chance of every alteration; with a few
patients carrying most alterations it calls co-occurrence everywhere (72.3% of
pairs "significant" originally). A DISCOVER-style rate-adjusted test fixes that
but estimates each patient's rate from their own few alterations, which pushes
pairs towards "exclusive". On simulated data with **no** interactions
(MSK-IMPACT468 benchmark, `scripts/validation/`):

| test | null pairs with p < 0.05 (target 5%) | significant after FDR |
|---|---|---|
| plug-in rate model (DISCOVER-style) | 46.1% | 33.36% |
| Rediscover R package | 69.4% | 598 of 990 pairs |
| **conditional test (used)** | **2.9%** | **0.00%** |

The conditional test recovers 5/5 planted interactions; its Poisson-Binomial
tail matches Rediscover's on 990 of 990 benchmark pairs and 60-digit
arithmetic to ~1e-12.

**Known biology recovered:** *KRAS*/*EGFR*, *KRAS*/*BRAF*, *EGFR*/*IDH1*,
*PIK3CA*/*PTEN* (breast) exclusive; *MDM2* amplification vs *TP53* mutation
exclusive across cancer types; the breast-specific 11q13/8p12
co-amplification; *KEAP1*/*STK11* co-mutation in lung cancer.

### Phase 3 -- the lookup tool (Aim 1)

`notebooks/03_phase3_lookup_tool.ipynb` answers a clinician's five questions:

| Question | Function | How |
|---|---|---|
| What else is usually altered with this? | `lookup()` | panel-aware counts in carriers vs non-carriers (95% interval), by TMB group and optionally subtype, with Phase 2's verdict (`candidate, replicated` / ...) |
| How common is this exact combination? | `combo_prevalence()` | exact counts among patients tested for every gene |
| Is this (rare) pair really over- or under-represented? | `test_pair()` | Phase 2's conditional test on demand, for variants and for pairs below Phase 2's filters (exploratory when < 5 expected) |
| Given the whole profile, what else is likely? | `predict_partners()` | one regularised logistic model per gene |
| What does it mean for treatment? | `therapy_evidence()` | CIViC evidence, including combinations |

Queries can be genes, variants (`("KRAS", "p.G12C")`), copy-number changes
(`("MET", "AMP")`) or fusions (`"ALK_fusion"`).

How the methods were chosen (all on held-out or independent data):
- **Matching a patient to similar cases** (the registered aim) was tested four
  ways. Whole-profile similarity -- plain or rarity-weighted Jaccard, with or
  without the same subtype -- does *not* beat the cancer-type base rate
  (AUROC 0.81-0.82 vs 0.82), because panel profiles carry only 2-4 alterations.
  **Per-gene logistic models** do, in every cancer type (0.855), and are well
  calibrated. Exact matches are still counted by `combo_prevalence()`.
- **On-demand test** reproduces Phase 2 on pairs Phase 2 tested (KRAS/STK11:
  expected 460.5 vs 460.4).
- **Frequencies are reproducible:** MSK vs other centres correlate at
  0.89-0.996 (median difference < 1 percentage point).
- **Therapy evidence:** OncoKB needs a licence token, so the open CIViC
  knowledgebase is used (3,433 accepted predictive / prognostic items).

Worked example in the notebook: *BRAF* V600E + *NRAS* in melanoma -- 7 of 4,144
patients, 27x rarer than expected, and CIViC evidence (level B) that the
combination confers resistance to BRAF inhibitors.

Limit: co-occurrence is not treatment prediction -- GENIE's main release has no
treatment or response data; therapy information comes from CIViC.

### Phase 4 -- what the patterns tell us (Aim 2)

`notebooks/04_phase4_insight.ipynb` tests claims from the literature on the
Phase 2 results:

1. **Exclusivity is pathway redundancy; co-occurrence is not pathway
   cooperation.** Exclusive pairs are 3.3x (replicated: 4.5x) enriched within
   one of the 10 TCGA oncogenic pathways (EGFR/KRAS, BRAF/NRAS, MDM2/TP53,
   PIK3CA/PIK3R1); co-occurring pairs are *not* enriched across pathways
   (OR ~1). Exclusive pairs across pathways (EGFR vs STK11 / KEAP1 in lung)
   are candidate incompatibilities.
2. **A third of interactions depend on the cancer type** (Cochran's Q, 538
   pairs tested in >= 3 types): 47 reverse direction -- EGFR/TP53 (co-occurring
   in lung, exclusive in glioma), KRAS/TP53 (co-occurring in pancreas,
   exclusive in lung and colorectal), CCND1/TP53 (head and neck vs breast);
   21 are universal (ATM/TP53, CDKN2A/RB1, MDM2/RB1).
3. **Partners depend on the variant, and smoking confounds some of it.**
   EGFR L858R vs exon 19 deletion differ in RBM10 (q = 4e-22); PIK3CA helical
   vs kinase in GATA3 and PTEN. Adjusting for each patient's smoking signature
   (C>A share) shrinks KRAS G12C/STK11 from 2.00x to 1.74x and removes
   G12C/RBM10 entirely; the non-smoking allele G12D is unchanged.
4. **Rare combinations exist and are counted:** 351 replicated exclusive
   pairs, 130 seen in 1-20 patients (ALK fusion + STK11 in lung cancer 0 vs 60
   expected; GNAQ+GNA11 in melanoma 0 vs 44; H3F3A+IDH1 in glioma 1 vs 43).
5. **Hypermutated tumours mainly add subtype structure** (MSI vs MSS in
   colorectal cancer), not new interactions.
6. **Mutation-only and mutation-or-copy-number analyses agree** in direction on
   98.9% of pairs significant in both.
7. **What is new?** Every robust pair (386) was checked against the literature
   (`19_novelty_screen.py`, Europe PMC titles/abstracts) and against the SelectSim
   atlas (Iyer ... Ciriello, Nat Genet 2026; 892 validated mutation-only pairs).
   41 are in neither, and they are amplicon passengers (GLI1~CDK4, KDR~PDGFRA,
   RARA~ERBB2, MTAP~CDKN2A, CD79B) or subtype markers. The strongest candidate
   after subtype, third-gene, burden and TMB adjustment (`20_candidate_followup.py`),
   CTCF vs TP53 in endometrial cancer, is in SelectSim's atlas -- we confirm it under
   stricter adjustment. Every robust pair SelectSim reports for the same cancer
   type agrees in direction. The distinct contribution is copy-number and fusion
   partners (outside SelectSim's mutation-only scope), panel-aware statistics and
   the artefact audits.

Also in the notebook: the copy-number burden confound (CD79B, 4 of 25
copy-number-driven pairs are artefacts) and NF1+PTPN11 as a cooperating
exception to the pathway rule.

Outputs: `cancer_type_heterogeneity.parquet`, `rare_combinations.parquet`,
`allele_*.parquet` (`scripts/analyses/08_allele_screen.py`), `novelty_screen.parquet`,
`candidate_followup.parquet`.

## How far the results can be trusted

All results below are collected, with tables, in `notebooks/05_validation_and_methods.ipynb`.

- **External cohort (TCGA PanCanAtlas, whole exome, independent patients;**
  `scripts/analyses/14_tcga_external_validation.py`): replicated candidates keep
  their direction in 94% of testable pairs, are significant in TCGA 46% of the
  time (GENIE-null pairs: 7%, about chance) and **none reverses**.
- **Method benchmark** (`15_method_benchmark.py`): on the same data, Fisher
  (cBioPortal, maftools) calls 11% of null pairs significant and a
  DISCOVER-style test 33% (18% even after FDR); the conditional test 2% (0%
  after FDR), and its calls replicate in TCGA best (37% vs 29% for Fisher).
- **Independent route** (`17_conditional_selection_dnds.py`, Coselens idea):
  gene B's ratio of protein-changing to silent mutations shifts with gene A in
  the direction Phase 2 predicts (84%; 100% of the significant cases).
- **Exclusive gene sets** (`16_exclusive_gene_sets.py`): a load-aware set score
  finds whole pathways 3x as often as Dendrix's (15% vs 5% of top sets).
- **Partial pooling across cancer types** (`18_hierarchical_shrinkage.py`,
  Scharpf idea): modestly improves agreement with TCGA for noisy estimates.

## Figures

Each notebook has a "Publication figures" section, just before its Discussion, showing the figures
that come from it (last column below).
Publication figures are drawn in R (ggplot2, patchwork, ComplexHeatmap), following the figure
types used in the co-mutation literature. Rebuild them all with `bash scripts/figures/run_all.sh`
(a Python step exports small CSVs, because R cannot read parquet here). Output:
`results/figures/*.pdf` (vector, for the thesis) and `*.png` (300 dpi, for slides).
Colours: green = co-occurring, purple = exclusive (ColorBrewer PRGn); methods in Okabe-Ito;
oncoprints in cBioPortal colours.

| Figure | Panels | Modelled on | Shown in notebook |
|---|---|---|---|
| 1 Data | patient flow; cohort by centre; panel metadata vs released data; alterations per patient | GENIE 2017 Figs 1-3 | 01 |
| 2 Can the test be trusted? | null p-value histograms; head-to-head benchmark; GENIE vs TCGA effects | DISCOVER 2016 Fig 2; SelectSim 2026 Fig 1; Park & Lehner 2015 | 05 (pointer in 02) |
| 3 Landscape | triangle heatmaps (lung, breast, colorectal); volcano; lung network | maftools; Sanchez-Vega 2018 Fig 6; Park & Lehner 2015 Fig 3 | 02 |
| 4 Biology | pathway x pathway matrix; within vs between pathways; pairs across cancer types; burden confound | Sanchez-Vega 2018 Fig 6; van de Haar 2019 Fig 3 | 04 |
| 5 Fusions and variants | lung oncoprint with ALK/ROS1/RET; fusion partners; smoking adjustment | Skoulidis 2015 Fig 2 | 04 |
| 6 The tool | KRAS G12C partners; partner prediction; MSK vs other-centre frequencies | GENIE 2017 Fig 3 | 03 |
| S1-S6 | filter funnel; UpSet of lung drivers; KRAS/TP53 forest plot; rare combinations; dN/dS route; gene sets | Lex 2014; maftools Fig 4A; CoMEt 2015 Fig 3 | S1: 02; S2: 03; S4: 04; S3, S5, S6: 05 |

## Problems found in the GENIE data itself

Checked directly against the data rather than trusting documentation. Each
one would have silently biased results.

| Problem | Evidence | What we do |
|---|---|---|
| **Panel metadata overstates copy-number testing** | `assay_information.txt` says 48 panels report copy number; 20 of them have no samples in `data_CNA.txt`. 98,963 samples were being counted "no copy-number change" when never tested. *CDKN2A* deletion in glioma: 17.4% -> 34.4% once fixed (published ~30-40%). | Copy-number testing read from `data_CNA.txt` per sample (`CNA_TESTED`) |
| **Copy-number gene lists differ from mutation gene lists** | 16.2% of (panel, gene) pairs sequenced for mutations have no copy-number values; one panel (`YALE-HSM-V1`) reports copy number for none of its 50 genes. ~199,000 sample-gene cells affected. | Separate copy-number coverage table, read from `data_CNA.txt` |
| **Panel gene lists overstate mutation reporting too** | Comparing reported mutations with each gene's rate on other panels: 1,144 (panel, gene) pairs in 47 panels expect >= 10 mutated patients but report **zero**. Worst panels: `JHU-500STP` lists 760 genes, reports mutations in 46; `VICC-02-XFV2` lists 105, reports none; `CHOP-FUSIP` is a fusion-only panel whose genes are listed anyway; the exome `UHN-WGS-V1` lists 18,849 genes, reports 25. Also 699 pairs are "covered" only by intron (fusion) probes -- 96% report nothing. 6.9% of patient-gene cells across all patients; 3.8% in Phase 2's patients (mostly spared by the copy-number restriction; worst Phase 2 genes: *TEK*, *PIK3C2G*, *PPARG* in melanoma, 6-9% of tested patients). `clinicalReported = False` is **not** a problem (those genes report normally). | **Partly fixed, from the data:** samples on the 8 panels that released no mutations (2,974 samples, e.g. `VICC-02-XFV2`, `CHOP-FUSIP`, `PROV-MPNC`) count as not tested for mutations (`MUT_TESTED`), and genes count as tested only where the panel has exon probes (52,562 instead of 53,291 panel-gene pairs). Gene-level gaps that only statistics can reveal (e.g. `JHU-500STP`) are left as they are |
| **Panel metadata misstates fusion detection** | 9 panels claim structural-variant detection but released none; 3 released fusions without claiming to | Fusion testability read from `data_sv.txt` (`SV_TESTED`, `sv_gene_panel_coverage`) |
| **Germline filter removed real cancer mutations** | At gnomAD 1e-4, ~104k somatic calls were removed -- mostly acquired drivers that reach "healthy" gnomAD donors through clonal haematopoiesis (*JAK2* V617F lost 90% of its calls; *DNMT3A* R882, *SF3B1* K700E, *MYD88* L265P, *EGFR* T790M). | Cutoff 1e-3 |
| **Copy-number bystanders** | 85,254 deep calls go the wrong way for the gene's role (oncogene deleted / tumour suppressor amplified) -- genes riding next to the real target. | Dropped (OncoKB roles) |
| **Shallow copy-number calls not comparable across centres** | MSK panels report +/-1 in ~5% of samples, DFCI panels in 98-99%. | Deep calls only (Phase 1, step 4) |
| **TMB meaningless on small panels** | On panels < 1 Mb, ~61% of samples come out as TMB >= 10. | TMB trusted only on panels >= 1 Mb; others "unknown" |
| **Corrupted gene coordinates** | 33 of 19,605 genes in `genomic_information.txt` have impossible spans (e.g. *ZNF703* ~97 Mb). | Excluded from the gene-distance check |

## Validation scripts

`scripts/README.md` lists every script with the question it answers, its output and where the
result is discussed, plus the full run order. Validation scripts load the test functions directly
from the Phase 2 notebook, so what is validated is exactly what the notebook runs.

## License

Code is released under the MIT License (`LICENSE`). It does not cover GENIE data or
other third-party data, which remain under their own terms (see `data/README.md`).
