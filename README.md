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
                                      across cancer types)
  03_phase3_gene_clusters.ipynb      Multi-gene (3+) cluster discovery
  04_phase4_survival_association.ipynb  Cluster-survival association

data/
  raw/         Original GENIE files (not tracked in git -- see Data access)
  external/    Downloaded reference data: AlphaMissense, cancerhotspots.org,
               OncoKB gene list (not tracked in git)
  processed/   Cleaned outputs produced by the notebooks (tracked in git)

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

**Phase 1 (preprocessing) is current.** Phase 2 (co-mutation matrix) is
complete but its committed outputs predate Phase 1's latest revision
(copy-number testability fix + bystander filter) and are being re-run.
Phases 3 and 4 have a working first pass but are not yet re-run either.

Current headline numbers (from `data/processed/qc_summary.json`):

| | |
|---|---|
| Samples | 271,837 (167 sequencing panels) |
| Samples with copy-number data | 172,874 (63.6%) -- 20 of the 48 panels claiming CNA support have none |
| Mutations kept after pathogenicity filter | 1,597,106 |
| Deep CNV calls (+2/-2) | 377,216 -> 291,962 after dropping direction-inconsistent bystanders |
| Combined alteration events | 1,889,068 (229,185 samples) |
| Gene pairs tested (Phase 2, per-cancer-type) | 19,454 |
| Significant pairs (q < 0.05) | 9,161 (47.1%) |
| Pan-cancer recurrent pairs (≥5 cancer types, consistent direction) | 215 |

**Pathogenicity filter** (missense mutations): `AlphaMissense = pathogenic
OR (Polyphen = damaging AND SIFT = deleterious) OR cancerhotspots.org
residue match` — approved by the supervisor after comparing 4 evidence
sources and multiple candidate rules (see Phase 1, Section 3).

**CNV filter**: deep amplification / deep deletion only, then each call is
checked against the gene's OncoKB role -- a tumour suppressor *deleted* or an
oncogene *amplified* is kept; the reverse (almost always a bystander on a
large arm-level event) is dropped. Genes OncoKB doesn't curate are kept.
Copy-number testability is read from `data_CNA.txt` directly rather than
trusted from panel metadata, which is wrong for 20 panels (see Phase 1,
Sections 2 and 4).

**Gene list for comparison across samples**: built per cancer type (not
globally pooled) — a gene qualifies for a cancer type if ≥80% of that
cancer type's own samples were tested for it, with ≥100 tested samples as
a floor (see Phase 1, Section 8).
