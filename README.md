# Learning Cancer Mutation Landscapes from Large-Scale Somatic Mutation Testing Data

Final Year Project (HKU, B.Sc. Bioinformatics) analysing the [AACR Project
GENIE](https://www.aacr.org/professionals/research/aacr-project-genie/)
public dataset (v19.0, ~270,000 clinical tumour sequencing samples) to find
patterns of somatic mutation co-occurrence across and within cancer types,
and lay the groundwork for matching a new patient's mutation profile to
genetically similar cases.

**Supervisor:** Prof. Jason Wong

## Project aims

1. Identify genes whose mutations co-occur more frequently than expected,
   both **within** a cancer type and **across** cancer types.
2. Develop a computational method that matches a new patient's mutation
   profile to genetically similar cases in the database, and reports how
   common that combination is.

See [`PLAN.md`](PLAN.md) for the full phased methodology, supervisor
correspondence context, and open questions.

## Repository structure

```
notebooks/
  01_phase1_preprocessing.ipynb      Data cleaning, pathogenicity filtering,
                                      panel-coverage tracking, gene lists
  02_phase2_comutation_matrix.ipynb  Pairwise co-mutation matrix (within and
                                      across cancer types) + mechanism-specific
                                      tests (SNV/CNV on each side of a pair)
  03_phase3_gene_clusters.ipynb      Multi-gene (3+) cluster discovery
  04_phase4_survival_association.ipynb  Cluster-survival association
  05_phase2b_snv_cnv_exploration.ipynb  Exploratory: SNV/CNV split, CNV driver
                                      filter, and the CNA-testability discovery
  06_phase2c_alteration_states.ipynb  Exploratory: per-gene alteration states,
                                      mechanism-specific tests, subtype checks

data/
  raw/         Original GENIE files (not tracked in git -- see Data access)
  external/    Downloaded reference data: AlphaMissense, cancerhotspots.org,
               OncoKB gene list (not tracked in git)
  processed/   Cleaned outputs produced by the notebooks (tracked in git)

validation/    Checks of the Phase 2 test against simulation, scipy/mpmath
               and the Rediscover R package (outputs in data/validation/)
PLAN.md        Full methodology, phase-by-phase plan, open questions
requirements.txt / setup.sh   Environment setup
```

Every notebook is self-contained (no shared script dependency), starts with
a "data flow at a glance" diagram, and ends with a Discussion section
(insight / result / limitations / open questions).

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

## Status

**Phase 1 (preprocessing) and Phase 2 (co-mutation matrix) are current.**
Notebooks 05 and 06 are the exploratory work behind earlier Phase 2 changes,
kept for the record. Phases 3 and 4 have a first pass but have **not** been
re-run on the current Phase 1/2 outputs -- Phase 3 should be rebuilt from
`interaction_candidate` pairs only.

**Not yet supervisor-reviewed** (flag at the next check-in): the conditional
pair test, the germline cutoff change (gnomAD 1e-4 -> 1e-3), one sample per
patient, TMB-based hypermutation status, the expected-count gate, and the
shared-DNA filter. Each is documented with its evidence in the notebooks.

### Headline numbers

| | |
|---|---|
| Samples / patients | 271,837 samples, 167 panels -> **227,696 patients** (one sample each) |
| Samples with copy-number data | 172,874 (63.6%) -- 20 of the 48 panels claiming CNA support released none |
| Mutations kept (pathogenicity filter) | 3,458,550 -> **1,644,226** |
| Deep CNV calls | 377,216 -> 291,962 after dropping direction-inconsistent bystanders |
| Combined alteration events | 1,936,188 (230,886 samples) |
| Gene pairs tested (Phase 2) | 33,092 across 53 cancer types (134,087 patients) |
| Significant (q < 0.05) | 4,657 (14.1%) |
| ...also significant in non-hypermutated patients | 2,096 (6.3%) |
| ...and not neighbouring genes / one shared DNA event | **1,847 interaction candidates (5.6%)** -- 556 co-occurring, 1,291 exclusive |
| Mechanism-specific tests (SNV/CNV per side) | 38,406 tests, 5,378 significant, 1,092 cross-mechanism |
| Pan-cancer candidates (>= 5 cancer types, same direction) | 13 |

Jason's question -- *"half being significant sounds quite high"* -- was right:
72.3% (Fisher) -> 49.6% (rate-adjusted) -> **14.1% significant, 5.6% candidates** now.

## Methods in brief

**Phase 1 -- one clean, panel-aware table.**
- *Mutations*: somatic, protein-changing, gnomAD max AF <= 1e-3, and for
  missense `AlphaMissense pathogenic OR (Polyphen damaging AND SIFT
  deleterious) OR cancerhotspots.org residue` (supervisor-approved). The
  germline cutoff was 1e-4 until checked: it removed ~104k somatic calls,
  dominated by acquired drivers seen in gnomAD through clonal haematopoiesis
  (JAK2 V617F lost 90% of its calls; DNMT3A R882, SF3B1 K700E, MYD88 L265P,
  APC E1309Dfs, EGFR T790M).
- *Copy number*: deep calls only; a TSG deleted or oncogene amplified is kept,
  the reverse (a bystander of a large event) dropped (OncoKB roles).
- *Testability*: which genes each panel sequences (mutations) and where
  copy-number values actually exist (read from `data_CNA.txt`, which disagrees
  with panel metadata for 20 panels and 16.2% of panel-gene pairs).
- *Patients*: one representative sample per patient (has CNA data > primary >
  larger panel) -- a primary and its metastasis share trunk mutations.
- *Hypermutation*: official GENIE TMB, trusted only on panels >= 1 Mb
  (on smaller panels it calls ~61% of samples "TMB >= 10", which is noise).

**Phase 2 -- a pair test that compares like with like.** For each pair, the
chance a patient carries both genes is computed *given that patient's own
number of alterations* (conditional / Rasch likelihood; gene parameters by
conditional ML within TMB strata), summed into a Poisson-Binomial null with
an exact tail, then BH-FDR within each cancer type. A pair is tested only if
>= 50% of patients were tested for both genes and >= 5 are *expected* to
carry both (decided before seeing the result). Significant pairs are
re-tested in non-hypermutated patients, and pairs that are neighbours
(<= 10 Mb) or one copy-number event on one chromosome are excluded from the
interaction candidates.

**Why not the earlier tests.** Fisher's test assumes every patient has the
same chance of every alteration; with a few patients carrying most
alterations it calls co-occurrence everywhere. The DISCOVER-style
rate-adjusted test fixes that but estimates each patient's rate from their
own few alterations, which biases pairs towards "exclusive". On simulated
data with **no** interactions (MSK-IMPACT468 benchmark, `validation/`):

| test | null pairs with p < 0.05 (target 5%) | significant after FDR |
|---|---|---|
| plug-in rate model (DISCOVER-style) | 46.1% | 33.36% |
| Rediscover R package | 69.4% | 598 of 990 pairs |
| **conditional test (used)** | **2.9%** | **0.00%** |

The conditional test recovers 5/5 planted interactions; its Poisson-Binomial
tail matches Rediscover's on 990 of 990 benchmark pairs and 60-digit
arithmetic to ~1e-12.

**Known biology recovered**: KRAS/EGFR, KRAS/BRAF, EGFR/IDH1, PIK3CA/PTEN
(breast) exclusive; MDM2 amplification vs TP53 mutation exclusive across
cancer types; the breast-specific 11q13/8p12 co-amplification; across
cancers CDKN2A/RB1, ATM/TP53 and MDM2/TP53 exclusivity.

## Validation scripts

`validation/` reproduces every check above from the processed data
(`build_msk468_matrix.py` first; `rediscover_run.R` needs R + Rediscover).
They load the test functions directly from the Phase 2 notebook, so what is
validated is exactly what the notebook runs.
