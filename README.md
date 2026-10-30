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
  02b_phase2_validation.ipynb        Phase 2b: validation (precision and recall, known biology, TCGA)
  02c_phase2_extensions.ipynb        Phase 2c: side analyses (hypermutated, mutations only, burden, fusions, method comparison)
  03_phase3_lookup_tool.ipynb        Phase 3: the lookup tool (Aim 1)
  04_phase4_insight.ipynb            Phase 4: what the patterns tell us (Aim 2), incl. novelty audit
  05_phase2_deep_dive.ipynb          Trust checks on Phase 2 hits (subtype, burden, proximity, passengers)

scripts/                  See scripts/README.md (index of every script and the run order)
  analyses/    One question per script (01-20), behind each design decision
  validation/  Checks of the Phase 2 test itself: simulation, exact tail, Rediscover (R)
  figures/     R publication figures; run_all.sh rebuilds all of them
  fetch_external.py, build_notebook_02b_validation.py

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
were removed; they are in the git history (last present in commit `0ce282a`).
The original project plan (`PLAN.md`: clusters, survival, Fisher's test, a
data-file inventory) was deleted once it no longer matched the work; it is in
the git history (last present in commit `7009281`).

Every notebook is self-contained, starts with a "data flow at a glance"
diagram, and ends with a Discussion section.

## Data access

GENIE data cannot be redistributed. `data/raw/` and `data/external/` are
gitignored, and so are the individual-level tables in `data/processed/`
(`clinical_tidy`, `alterations_long`, `cohort`, `events`, `putative_biallelic_tsg_events`, rebuilt by
notebook 01); only aggregate outputs are tracked, and notebooks print table
schemas rather than patient rows. To reproduce this pipeline, request GENIE v19.0 access via
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
| Mutations kept (likely harmful) | 3,458,550 -> **1,659,562** (incl. 15,103 TERT promoter hotspot and 233 MET exon 14 splice calls) |
| Deep copy-number calls | 377,216 -> 291,962 after dropping direction-inconsistent bystanders |
| Driver gene fusions | 56,356 structural variants -> **15,847** driver-fusion calls in 432 genes (published pair in Mitelman, or a recurrent fusion oncogene; ERG, BRAF, ALK, EWSR1, RET ...) |
| Combined alteration table | 1,967,371 rows (232,820 samples): mutations, copy number, fusions |
| Phase 2 main analysis (SNV+CNV) | non-hypermutated patients: 6,809 gene pairs (224 genes, incl. the TERT promoter), 52 cancer types, 113,042 patients |
| Significant (q < 0.05) | 1,811 (26.6%) |
| ...and not one copy-number event hitting neighbouring genes | **1,569 interaction candidates** -- 419 co-occurring, 1,150 exclusive |
| ...also replicated outside MSK (significant in MSK, confirmed in other centres) | **472**; of all 814 MSK hits testable elsewhere, 65.5% replicate, none in the opposite direction |
| ...survive the copy-number burden sensitivity check | 1,264 of 1,569 candidates |
| External validation in TCGA (replicated candidates) | 94% same direction, 46% significant (vs 7% for GENIE-null pairs), **0 reversed** |
| Fusion screen | 375 pairs, 18 cancer types, 132 significant (e.g. ALK/ROS1/RET fusions vs EGFR/KRAS/STK11) |
| Hypermutated patients on their own (TMB 10-100) | 881 of 99,429 pairs significant (0.9%), 18 cancer types |
| Other alteration sets, same patients and test | SNV only: 2,174 pairs, 466 candidates (76% of MSK hits replicate); SNV+CNV+SV: 6,871 pairs, 1,609 candidates (65%) |
| Mechanism-specific tests (mutation / copy number per side) | 8,661 tests, 2,263 significant, 701 cross-mechanism |
| Pan-cancer pairs (>= 5 cancer types, same direction) | 13 |

Jason's question -- *"half being significant sounds quite high"* -- was right
about the method: Fisher's test called 72.3% of pairs significant. The share is
now 26.6% of a much smaller, driver-only pair list (hypermutation no longer
pushes passenger genes into the analysis), and only 1,569 pass as candidates.

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
3. **Mutations: keep the likely harmful ones.** 3.46M calls -> 1,659,562 (including two exceptions to "protein-altering only": TERT promoter hotspots C228T / C250T as their own feature, and MET exon 14 splice changes).
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
5. **Gene fusions.** `data_sv.txt` holds 56,356 structural variants. An event joining two *different* genes
   (not antisense) is a driver fusion if its gene pair is **published in the Mitelman database** (CC BY 4.0) **or**
   one partner is a **recurrent fusion oncogene** (>= 5 published partners in Mitelman and an OncoKB oncogene):
   **15,847** driver-fusion calls in 432 genes, named after the driver (the OncoKB oncogene; e.g. `ALK_fusion`).
   Frame is not used: GENIE's frame labels are missing for most events and call famous drivers such as
   TMPRSS2-ERG "out of frame". A sample
   counts as fusion-tested only if its panel actually released structural
   variants (35 panels do; read from the data, `SV_TESTED`).
6. **Output.** `alterations_long.parquet` (1,967,371 rows, 232,820 samples:
   1,659,562 mutations + 291,962 deep copy-number calls + 15,847 fusions),
   `clinical_tidy`, `panel_gene_coverage`, `cna_gene_panel_coverage` and
   `sv_gene_panel_coverage`. Phase 1 does **not** choose genes -- every
   analysis filter lives in Phase 2.
7. **Hand-off to Phase 2** (Section 9). Three tables, checked before they are written
   (`scripts/validation/check_phase1_handoff.py`), and the only Phase 1 files Phase 2 reads:
   `cohort` (one row per patient: 227,696, with a TMB group -- 159,390 not hypermutated, 24,452
   hypermutated, 658 ultramutated, 43,196 TMB unknown), `events` (their 1,704,134 driver events:
   mutations, copy-number changes, fusions) and `testability` (which panel tests which gene feature
   for mutation, copy number and fusion). 99.8% of events were testable on their patient's panel.

### Phase 2 -- find gene pairs that go together or avoid each other

`notebooks/02_phase2_comutation_matrix.ipynb`, run separately for each cancer type; all filters are
set in one settings cell. The list is built three times, with three definitions of "altered" -- **SNV**
(mutations, incl. TERT promoter and MET exon 14), **SNV+CNV** (the main table) and **SNV+CNV+SV** (+ driver
fusions, kept separate until fusions are agreed) -- on the same patients with the same test; notebook 02b
compares them (SNV replicates best, SNV+CNV sees copy-number biology at the same precision).

1. **Collect the data; keep non-hypermutated patients** (TMB < 10 on a panel >= 1 Mb; 159,390
   patients). Prof. Wong's suggestion, and the design whose findings replicate best across centres
   (44% -> 66%).
2. **Filter.** Cancer types with >= 100 patients tested for mutations and copy number (**51 types,
   113,042 patients**); genes altered in >= 3% (and >= 5) of patients; pairs with >= 50% of patients
   tested for both genes and >= 5 patients *expected* to carry both (decided before seeing the result).
3. **The conditional test.** Each patient is compared with their **own number of alterations**: a
   tumour with 30 alterations will carry almost any pair by chance, one with 2 will not (conditional /
   Rasch model, exact Poisson-Binomial tail). Gene rates are fitted within detailed subtypes, which
   removes subtype artefacts (pancreatic *MEN1* vs *SMAD4*, glioma *CIC* vs *PTEN*, breast *CDH1* vs
   *ERBB2*). *Side comparisons (not used):* Fisher's test called 72.3% of pairs significant; the
   DISCOVER-style test is biased towards exclusivity on panel data (table below).
4. **FDR** (Benjamini-Hochberg within each cancer type): **6,809 pairs -> 1,811 significant (26.6%)** (main set).
5. **Remove one-DNA-event pairs.** Removed: pairs where one copy-number event hit
   both genes (a gained chromosome arm, or neighbours on one amplicon). Neighbours (<= 10 Mb) are removed
   only when copy number drives the pair -- two point mutations are two separate events, so *KEAP1* /
   *STK11* / *SMARCA4* (lung) and *PBRM1* / *SETD2* (kidney) stay. This comes after FDR because it
   re-labels significant pairs without changing any p-value. **1,569 interaction candidates**.
6. **Replication and the numbers.** Every pair re-tested in MSK patients and in all other centres: 65.5%
   of MSK hits confirmed, none reversed, 472 candidates `replicated`; the funnel and the three sets side by
   side (Phase 2 figure). All evaluation -- known biology, strongest pairs, precision and recall, TCGA, which
   set is better -- is in notebook 02b.
7. **Save** `comutation_pairs.parquet` for Phases 3 and 4.

**Why not Fisher's test or the earlier rate-adjusted test.** Fisher's test assumes every patient has
the same chance of every alteration; with a few patients carrying most alterations it calls
co-occurrence everywhere. A DISCOVER-style rate-adjusted test fixes that but estimates each patient's
rate from their own few alterations, which pushes pairs towards "exclusive". On simulated data with
**no** interactions (MSK-IMPACT468 benchmark, `scripts/validation/`):

| test | null pairs with p < 0.05 (target 5%) | significant after FDR |
|---|---|---|
| plug-in rate model (DISCOVER-style) | 46.1% | 33.36% |
| Rediscover R package | 69.4% | 598 of 990 pairs |
| **conditional test (used)** | **2.9%** | **0.00%** |

The conditional test recovers 5/5 planted interactions; its Poisson-Binomial tail matches Rediscover's
on 990 of 990 benchmark pairs and 60-digit arithmetic to ~1e-12. **Why both the TMB split and the
conditional test:** on null data from real non-hypermutated patients Fisher still calls 6-15% of pairs
at p < 0.05 (conditional test 1.3-2.9%), and the conditional test alone on all patients leaves 25% of
colorectal pairs significant vs 5.8% in non-hypermutated patients.

**Known biology recovered:** *KRAS*/*EGFR*, *KRAS*/*BRAF*, *EGFR*/*IDH1*, *PIK3CA*/*PTEN* (breast)
exclusive; *MDM2* amplification vs *TP53* mutation exclusive across cancer types; *KEAP1*/*STK11*
co-mutation in lung cancer.

### Phase 2c -- side analyses of the pair test

`notebooks/02c_phase2_extensions.ipynb` re-uses the Phase 2 test; none of it changes the 1,569
candidates. (Hypermutated tumours are now described in Phase 4, Section 5; mutations only is the SNV set of notebook 02.)
- *Hypermutated patients on their own* (TMB 10-100; ultramutated TMB >= 100 excluded -- pooled they made
  21.6% of endometrial pairs "significant", 0.3% and 0% when separated): 881 of 99,429 pairs
  significant, mostly subtype structure (MSI vs MSS); a description, not candidates.
- *Mutations only* (155,020 patients with mutation data): the most consistent between centres
  (replication 84%).
- *Copy-number burden flag* (Kumar 2011; Zack 2013): pairs re-tested within groups of genome
  instability; `burden_robust` -- 1,264 of 1,569 candidates pass. Built into the main test it
  over-adjusts (loses the known CCND1/FGFR1 co-amplification), so it is a flag.
- *Fusion screen* (fusions with >= 20 carriers): 132 significant pairs, e.g. ALK/ROS1/RET fusions avoid
  EGFR, KRAS, STK11 in lung cancer; ERG fusions avoid SPOP/FOXA1 and go with PTEN loss in prostate cancer.
- *Mutation-load check* (van de Haar 2019) and *method comparison on real data* (significant lists in
  TCGA, ranking by effect size).

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
7. **What is new?** Every robust pair (387) was checked against the literature
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

`notebooks/02b_phase2_validation.ipynb` measures the test's **precision** (of the pairs it calls, how many are
real) and **recall** (of the real pairs, how many it finds), equally weighted, comparing it with Fisher's test (cBioPortal, maftools), a DISCOVER-style test and SELECT (Mina 2017,
re-implemented) on identical data.

- **Data with no real links** (`15_method_benchmark.py`): every call is wrong here. Fisher calls
  11% of pairs at p < 0.05 and the DISCOVER-style test 33% (18% after FDR); ours 2% (0% after FDR), SELECT 4%.
  No candidate is significant in the opposite direction in TCGA.
- **Precision and recall** (`22_planted_signal_simulation.py`): links of known strength planted into
  data without links, in three scenarios (our model's assumption; real data with links removed; a
  misspecified stress test). On realistic data ours has precision 99% and recall 49%; SELECT 50% and 42%;
  Fisher and DISCOVER-style reach recall 60-63% only at a precision of 4-8%. **At the same precision ours
  finds the most real links** (recall 58% at 95% precision vs 13% for SELECT and 2% for the others). Weak co-occurrence and part of the exclusive links are missed. A hidden patient group the
  strata do not capture breaks every method, ours included -- the subtype strata and the separate
  hypermutated analysis are part of the method.
- **Known biology:** 26 of 29 textbook pairs chosen from the literature beforehand
  (`data/curated/textbook_pairs.csv`) are found (Fisher 28, DISCOVER-style and SELECT 27); the misses are weak
  pairs. The 10 strongest exclusive pairs are all established biology; of the 10 strongest co-occurring pairs
  5 are established, 1 is a suspected artefact and 4 are not yet explained (`data/curated/top_pairs_biology.csv`).
  Of 387 robust pairs, 235 are well known; every pair also in the SelectSim atlas agrees in direction.
- **Independent patients (TCGA PanCanAtlas,** `14_tcga_external_validation.py`): candidates keep their
  direction in 88% of testable pairs (control pairs: 69%), effect sizes correlate (Spearman 0.74; 0.82 for
  replicated candidates), and agreement rises with the strength of the GENIE signal.
- **Mutation-load confounding** (notebook 02c, Section E; van de Haar 2019 diagnostic): absent once
  hypermutated tumours are analysed separately, for every method.

## Figures

One combined figure per phase (up to four panels), titled with the phase's main message. Each
notebook draws its own figure by running the phase's R script from a code cell, so the figure is
always rebuilt from that phase's saved outputs; `bash scripts/figures/run_all.sh` rebuilds all of them.
R (ggplot2, patchwork, ComplexHeatmap); a Python step first exports small CSVs, because R cannot read
parquet here. Each figure is saved as PDF (thesis) and 300-dpi PNG (slides). Green = co-occurring,
purple = exclusive.

| Figure (`results/figures/`) | Main message | Panels | Notebook |
|---|---|---|---|
| `phase1_data/figure_phase1` | 113,042 patients and their driver events, ready for Phase 2 | patient flow and every filter · TMB groups · filters per alteration type · events Phase 2 receives | 01 |
| `phase2_pair_test/figure_phase2` | The gene-pair list: 1,569 candidates, 472 replicated | three alteration sets · funnel to candidates (main set) · candidates per cancer type | 02 |
| `phase2b_validation/figure_phase2b` | The highest precision, and the most real links found at any given precision | precision in simulation · recall by link strength · precision-recall curves · recall of textbook pairs · GENIE vs TCGA | 02b |
| `phase2c_extensions/figure_phase2c` | Fusions avoid the other drivers | lung oncoprint · fusion partners | 02c |
| `phase3_lookup_tool/figure_phase3` | The tool's answers hold up | KRAS G12C partners · partner prediction · MSK vs other centres · exact combinations | 03 |
| `phase4_biology/figure_phase4` | Exclusivity follows pathways; pairs change with the cancer type | pathways · across cancer types · variants and smoking · rare combinations | 04 |

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
