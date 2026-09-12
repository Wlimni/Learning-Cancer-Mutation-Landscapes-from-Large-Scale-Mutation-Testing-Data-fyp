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

**Phase 1 (preprocessing) and Phase 2 (co-mutation matrix) are current**,
including supervisor-approved methodology revisions. Notebooks 05 and 06 are
the exploratory work that led to Phase 2's latest changes; kept for the
record, not part of the pipeline. Phases 3 and 4 have a working first pass
but are not yet re-run against the latest Phase 1/2 outputs.

Current headline numbers (from `data/processed/qc_summary.json`):

| | |
|---|---|
| Samples | 271,837 (167 sequencing panels) |
| Samples with copy-number data | 172,874 (63.6%) -- 20 of the 48 panels claiming CNA support have none |
| Mutations kept after pathogenicity filter | 1,597,106 |
| Deep CNV calls (+2/-2) | 377,216 -> 291,962 after dropping direction-inconsistent bystanders |
| Combined alteration events | 1,889,068 (229,185 samples) |
| Gene pairs tested (Phase 2, per-cancer-type, CNA-tested patients) | 20,415 across 58 cancer types |
| Significant pairs (q < 0.05) | 8,725 (42.7%) |
| Mechanism-specific tests (SNV/CNV on each side, ≥10 co-occurring) | 12,116 tests, 8,472 significant |
| …of which cross-mechanism (invisible to any single-mechanism table) | 918 |
| Pan-cancer recurrent pairs (≥5 cancer types, consistent direction) | 174 |

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

**Phase 2 denominator**: the main co-mutation table is restricted to
patients with copy-number data, so "altered" (mutation OR copy-number
change) means the same thing for every patient in a contingency table.
Each pair is additionally tested four ways -- mutation vs mutation,
mutation vs copy-number, and the reverse, copy-number vs copy-number --
so a pair's signal can be attributed to a specific alteration-type
combination (see Phase 2, Section 6b). Structural variants are excluded
by supervisor decision (partial overlap with SNV/CNV events; see Phase 1
Discussion).

**Gene list for comparison across samples**: built per cancer type (not
globally pooled) — a gene qualifies for a cancer type if ≥80% of that
cancer type's own samples were tested for it, with ≥100 tested samples as
a floor (see Phase 1, Section 8).
