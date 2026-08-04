# FYP Plan — Co-mutation Cluster Discovery in AACR GENIE

**Data:** AACR GENIE public release v19.0 (271,837 samples, 3.46M mutation records, CNA matrix for 172,875 samples, 117 cancer types, 166 gene panels). Raw files live in `data/raw/`; cleaned/derived tables produced by this project live in `data/processed/`.

## Setup

Notebooks need a Python kernel with the right packages installed. One-time setup:

```
cd fyp
./setup.sh
```

This creates `.venv/`, installs everything in `requirements.txt`, and registers
a Jupyter kernel named **"Python (fyp-genie)"**. Open a notebook and select
that kernel (top-right in JupyterLab / VS Code's kernel picker) before
running cells. Re-running `./setup.sh` later is safe — it reuses the existing
`.venv` and just re-syncs packages.

## Official FYP project title & registered description

**Title:** *Learning Cancer Mutation Landscapes from Large-Scale Somatic Mutation Testing Data*

**Project description** (as registered):
- Analyse a database of over 200,000 clinical tumour sequencing samples to investigate patterns of somatic mutations across different cancer types.
- **Aim 1:** identify genes whose mutations co-occur more frequently than expected within and across cancers, providing insights into tumour biology.
- **Aim 2:** develop a computational method that matches a new patient's mutation profile to genetically similar cases in the database and reports the prevalence of the observed mutation combination.
- Combines statistical analysis, bioinformatics, and large-scale genomic data mining.

**Learning activities / methodology:**
- Explore and preprocess a large clinical cancer genomics dataset.
- Perform statistical analyses of mutation co-occurrence within and across cancer types.
- Develop algorithms to measure similarity between tumour mutation profiles.
- Evaluate different similarity metrics and computational approaches.
- Validate findings using independent subsets of data and appropriate statistical methods.
- Interpret results in the context of cancer biology and precision oncology.

**Intended student learning outcomes:**
- Apply bioinformatics and statistical methods to analyse large-scale genomic datasets.
- Develop and evaluate computational approaches for investigating complex biological data.
- Interpret genomic findings in the context of cancer biology and precision medicine.
- Critically evaluate analytical methods and communicate their strengths and limitations.
- Present scientific findings effectively through written reports and oral presentations.
- Demonstrate independent research, problem-solving and data analysis skills.

**How this maps onto the phased plan below:** Aim 1 = Phase 2 (pairwise co-occurrence) + Phase 3 (multi-gene clusters). Aim 2 is phrased as *similarity search* — match a new sample to its most genetically similar existing cases and report how common that combination is — which is slightly more general than "assign to one discrete cluster." Jason's cluster-based approach (Phase 5) is one concrete way to implement Aim 2: clusters act as pre-computed groups of similar cases, and assignment + reporting cluster prevalence satisfies the aim. If clusters turn out to be poorly separated in practice, a direct k-nearest-neighbours-style similarity search over individual cases (without discrete clustering) is the fallback that still satisfies Aim 2 as officially written — worth flagging as an explicit checkpoint in Phase 5 rather than committing to clustering only.

## Background & clinical motivation

The originating problem is the molecular tumor board (MTB) workflow at HKSH: a
patient's somatic NGS panel report comes back with a list of mutations, and
the board has to decide which ones are clinically actionable and how. The
original project framing (before it was narrowed to the concrete deliverable
below) had three parts:

1. **ESCAT-based actionability annotation** — classify each mutation by the
   [ESMO Scale for Clinical Actionability of molecular Targets (ESCAT)](https://www.esmo.org/scales-and-tools/esmo-scale-for-clinical-actionability-of-molecular-targets-escat),
   a six-tier framework (I–V, X) that ranks genomic alterations by strength of
   evidence for a matched targeted therapy, so cancer centres can grade
   mutations consistently instead of ad hoc literature review at every board.
2. **Co-mutation detection and highlighting** — flag clinically significant
   co-occurring mutation pairs, because single-gene actionability tiering
   misses interactions that change prognosis or drug response. Known examples
   from GENIE-scale analyses: *KRAS*/*EGFR* and *KRAS*/*BRAF* are largely
   mutually exclusive in NSCLC/colorectal cancer, and *PIK3CA* is mutually
   exclusive with *AKT1*/*PTEN* in breast cancer — but some co-mutation pairs
   are not just mutually exclusive, they carry independent prognostic or
   predictive weight (e.g. *STK11* co-mutation with *KRAS* in NSCLC is
   associated with reduced immunotherapy response). A board reading each gene
   in isolation can miss this.
3. **Similarity-based recommendations for non-actionable mutations** — most
   mutations on a panel report have no approved matched drug (ESCAT Tier X).
   For these, infer a "next-best" option (off-label use, pathway/domain
   similarity to a druggable target, or a relevant clinical trial) rather than
   leaving the board with nothing.

**Why AACR GENIE specifically, and why clusters over single genes:** a single
hospital's annual caseload for any specific cancer type is too small to
establish statistically robust co-mutation associations — you might see a
handful of patients per year with a given rare combination. [AACR Project
GENIE](https://pubmed.ncbi.nlm.nih.gov/28572459/) is a multi-institution
consortium pooling clinical-grade panel sequencing with outcome data
specifically to give researchers the numbers needed to validate genomic
biomarkers and find these patterns (~200k samples here — enough to look at
combinations of *more than two* mutated genes per cancer type, not just
pairs).

The supervisor's brief operationalizes point 2 (and implicitly extends point
3) into a concrete, buildable pipeline: **define co-mutation clusters from
GENIE, test them for survival association, then build an algorithm that
assigns a new patient sample to one of these clusters.** A cluster assignment
is itself a form of "similarity-based recommendation" — instead of reasoning
from abstract pathway/domain similarity, the board would see "this sample's
mutation profile matches empirically-derived cluster X, which in GENIE
patients was associated with [prognosis / enrichment in cancer type Y]" as a
decision-support signal, particularly useful for the majority of mutations
that have no direct ESCAT Tier I/II match. The ESCAT annotation piece (point
1) is not part of the current phased plan — it would sit alongside cluster
assignment as a separate, complementary annotation layer on top of the same
panel report.

## Aims (from supervisor's brief)

1. Define co-mutation clusters (>2 genes) based on AACR GENIE data, leveraging the ~200k sample size.
2. Analyze co-mutation clusters for clinical (mainly survival) associations, and cross-check interesting clusters in other datasets.
3. Develop an algorithm/software to assign a new sample to a co-mutation cluster defined in (1).

Explicit constraints from Jason:
- Must account for CNVs, not just point mutations.
- Must restrict to **likely pathogenic** mutations — silent/benign variants should not count as "mutated."

---

## Data inventory (confirmed 2026-07-06)

| File | Rows / size | Role |
|---|---|---|
| `data_mutations_extended.txt` | 3,458,551 rows / 1.0 GB | MAF — somatic mutation calls, one row per variant per sample |
| `data_CNA.txt` | 1,004 genes × 172,875 sample columns / 421 MB | GISTIC-style CNA calls (-2/-1/0/1/2), `NA` where gene not on that sample's panel |
| `data_clinical_sample.txt` | 271,837 rows | Sample-level metadata: `CANCER_TYPE`, `SEQ_ASSAY_ID`, `ONCOTREE_CODE` |
| `data_clinical_patient.txt` | 1 row/patient | `SEX`, `DEAD`, `YEAR_DEATH`, `INT_DOD`, `YEAR_CONTACT` (survival fields) |
| `genomic_information.txt` | 71 MB | Authoritative panel coverage: `SEQ_ASSAY_ID` × `Hugo_Symbol` × `includeInPanel` (True/False) — tells us which genes were actually assayed for each panel |
| `data_gene_panel_*.txt` (166 files) | small | Redundant with `genomic_information.txt`, simpler gene-list format per panel |
| `assay_information.txt` | 36 KB | Panel technical metadata (platform, gene count, alteration types covered — some panels are CNA-only or SV-only) |
| `tmb_19.0-public.tsv` | 271,838 rows | Precomputed tumor mutational burden per sample |
| `data_sv.txt` | 16 MB | Structural variants (fusions) — secondary priority |
| `meta_gene_matrix.txt` | tiny | Points to the gene-panel matrix file structure |
| `release-notes.pdf` | — | GENIE v19.0 release documentation |

**Key data-quality fact driving the whole design:** not every sample was sequenced on the same panel, so "gene X not reported for sample Y" can mean either *tested, wild-type* or *never tested*. In `data_CNA.txt` this is explicit (`NA` = not tested). In the mutation MAF it is implicit — must be reconstructed from `genomic_information.txt`. Every downstream frequency/co-occurrence calculation must use panel coverage as the denominator, not the full cohort, or results will be biased toward whatever genes are on the most common panels (e.g., MSK-IMPACT).

---

## Phase 1 — Data tidying & alteration matrix construction (current phase)

**Goal:** this phase has three jobs, not one:
1. **Organize** — collapse the five raw GENIE tables into one tidy, panel-aware **alteration table** that replaces the raw files for everything downstream.
2. **Filter** — reduce the raw mutation calls down to likely-pathogenic ones only: somatic (not inherited), protein-changing (not silent), and rare in the healthy population (not a common germline polymorphism). Phase 2's co-occurrence statistics are only meaningful if "altered" means "likely a real cancer-driving change," not "any raw variant call."
3. **Track testing coverage** — build `panel_gene_coverage.parquet` so every gene's "no mutation reported" can be disambiguated into either "tested, and clean" or "never tested." Example: if gene X is only on a small panel used for 800 of 5,000 Lung Cancer samples, and 40 of those 800 show it altered, the true rate is 40/800 = 5% — but treating all 5,000 as tested would wrongly give 40/5,000 = 0.8%, a number that just reflects which panel was used, not real biology. Without this step, Phase 2's frequency/co-occurrence math would silently use the wrong denominator, and results would just reflect which panels happen to be common rather than real biology.

1. **Clinical join**: merge `data_clinical_sample.txt` + `data_clinical_patient.txt` on patient ID → one row per sample with cancer type, panel ID, and survival fields attached.
2. **Panel coverage map**: parse `genomic_information.txt` (`includeInPanel == True`) into a `SEQ_ASSAY_ID → set of covered genes` lookup. This is what makes later denominators correct.
3. **Pathogenicity filter on mutations**:
   - Keep `Mutation_Status` ∈ {Somatic, SOMATIC}.
   - Keep protein-altering `Variant_Classification` only (Missense, Nonsense, Frame_Shift_Ins/Del, In_Frame_Ins/Del, Splice_Site, Translation_Start_Site, Nonstop_Mutation). Drop Silent, Intron, 5'/3'UTR, Flank, RNA, IGR.
   - Drop likely-germline calls using population frequency: exclude variants with `gnomAD_AF > 1e-4` (any subpopulation).
   - For missense variants specifically, require supporting evidence of deleteriousness (Polyphen "probably/possibly damaging" OR SIFT "deleterious") OR presence in a cancer hotspot list (cancerhotspots.org — to be downloaded separately) — missense without either signal is dropped or flagged low-confidence rather than silently kept.
   - Document the exact filter thresholds and counts-before/after in the notebook so the choice is auditable and revisable after discussion with Jason.
4. **CNA filter**: threshold `data_CNA.txt` to deep/high-confidence calls only (`-2` = deep deletion, `2` = amplification); treat `-1`/`1` as a separate, lower-confidence category, kept out of the primary matrix pending discussion. `NA` entries recorded as "not tested," not "wild-type."
5. **Unify into one long alteration table**: `(Sample_ID, Hugo_Symbol, Alteration_Type ∈ {MUT, CNV}, detail cols)`, tagged with `Cancer_Type` and `SEQ_ASSAY_ID` from the clinical join.
6. **Outputs** (Parquet, in `data/processed/`):
   - `clinical_tidy.parquet`
   - `panel_gene_coverage.parquet`
   - `alterations_long.parquet`
   - Summary QC report (counts before/after each filter, per-cancer-type sample counts, per-panel gene counts) so filter choices can be sanity-checked with Jason before Phase 2.

Deliverable this session: `notebooks/01_phase1_preprocessing.ipynb`, self-contained (all loading/filtering logic defined in the notebook itself, no separate script dependency).

---

## Phase 2 — Pairwise co-mutation matrix

1. Restrict to cancer types with adequate sample size (define a minimum N, e.g. ≥100, after seeing the Phase 1 distribution).
2. Per cancer type, build a binary gene × sample alteration matrix from `alterations_long`, masked by `panel_gene_coverage` so a gene-pair's contingency table only counts samples where **both** genes were actually assayed.
3. Fisher's exact test per gene pair (co-occurrence vs. mutual exclusivity), analogous to `maftools::somaticInteractions()`, but with the panel-aware denominator correction layered on top (maftools does not handle partial panel coverage natively).
4. Benjamini-Hochberg FDR correction across all pairs tested, per cancer type.
5. Output: long table of (Cancer_Type, Gene_A, Gene_B, odds_ratio/log-OR, p, q, direction).

## Phase 3 — Multi-gene (>2) cluster discovery

1. Frequent itemset mining (Apriori/FP-growth via `mlxtend`) on the panel-masked binary matrices to find recurrent 3+ gene combinations above a support threshold, OR graph community detection (e.g. Louvain) on the Phase 2 pairwise co-occurrence graph, thresholded by q-value.
2. Compare both approaches; pick whichever gives more clinically interpretable, stable clusters (stability checked via bootstrap resampling).
3. Name/characterize each cluster by its cancer-type distribution and gene composition.

## Phase 4 — Clinical association

1. Merge cluster assignments with survival fields (`DEAD`, `YEAR_DEATH`/`INT_DOD` relative to `YEAR_CONTACT` — note GENIE's survival data is coarse, year-level, and only available for a subset of centers).
2. Kaplan-Meier + log-rank test per cluster vs. rest; Cox proportional hazards adjusting for cancer type and age.
3. Flag clusters with significant associations for validation in an external cohort (e.g. TCGA, cBioPortal public studies) — Jason noted GENIE alone lacks rich clinical annotation.

## Phase 5 — New-sample cluster assignment algorithm

1. Given the clusters from Phase 3, define a similarity/distance metric between a new sample's alteration profile and each cluster's defining gene set (e.g. Jaccard similarity, or a trained classifier — logistic regression / gradient boosting — using cluster membership as labels).
2. Handle the case where the new sample's panel doesn't cover all cluster-defining genes (very likely in practice — HKSH panels differ from GENIE's) — restrict comparison to the intersection of tested genes, or flag reduced confidence.
3. Package as a small, documented function/CLI: input = sample's alteration list + panel ID, output = best-matching cluster(s) + confidence.

---

## Open questions to confirm with Jason before Phase 2

- Minimum sample-size cutoff per cancer type to include in matrix building.
- Whether `-1`/`1` (shallow) CNA calls should be included at all, or only `-2`/`2`.
- Acceptable pathogenicity filter (this plan uses gnomAD AF + Polyphen/SIFT + hotspot; Jason may prefer OncoKB annotations if API/offline data becomes available).
- Whether hematologic/heme panels (different alteration types, e.g. `CHOP-HEMEP`) should be modeled separately from solid tumor panels.
