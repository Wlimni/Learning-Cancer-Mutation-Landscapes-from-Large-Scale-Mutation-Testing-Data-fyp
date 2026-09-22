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

**Not yet supervisor-reviewed**: the switch from Fisher's exact test to a
DISCOVER-style rate-adjusted test (see "Statistical test" below) was built
in direct response to Jason's question about the significant-pair fraction,
with the evidence and controlled comparison documented in Phase 2's Section
4 -- but the change itself hasn't been signed off yet. Flag this
specifically at the next check-in.

Current headline numbers (from `data/processed/qc_summary.json`):

| | |
|---|---|
| Samples | 271,837 (167 sequencing panels) |
| Samples with copy-number data | 172,874 (63.6%) -- 20 of the 48 panels claiming CNA support have none |
| Mutations kept after pathogenicity filter | 1,597,106 |
| Deep CNV calls (+2/-2) | 377,216 -> 291,962 after dropping direction-inconsistent bystanders |
| Combined alteration events | 1,889,068 (229,185 samples) |
| Gene pairs tested (Phase 2, per-cancer-type, CNA-tested patients) | 25,609 across 52 cancer types |
| Significant pairs (q < 0.05) | 12,711 (49.6%) |
| Mechanism-specific tests (SNV/CNV on each side, ≥10 co-occurring) | 27,504 tests, 13,915 significant |
| …of which cross-mechanism (invisible to any single-mechanism table) | 1,637 |
| Pan-cancer recurrent pairs (≥5 cancer types, consistent direction) | 187 |

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

**Statistical test (Phase 2)**: gene pairs are tested with a DISCOVER-style
rate-adjusted test (Canisius et al. 2016), not Fisher's exact test. Fisher's
test assumes every patient has the same baseline chance of any gene being
altered, which is badly wrong when mutation burden is skewed within a
cancer type (e.g. Melanoma: median 6 altered genes/patient, max 494) --
checked directly, this inflated apparent co-occurrence specifically (90-94%
of Fisher-significant pairs were "co-occurring" vs. only 6-9% "exclusive").
The fix fits a per-gene, per-patient rate model (iterative proportional
fitting) and tests each pair against the Poisson-Binomial null those rates
imply. Controlled comparison (same gene lists, same floors, test statistic
only): significant pairs 72.3% → 49.6%, co-occurring:exclusive split
94.4%:5.6% → 82.9%:17.1%. Known biology (KRAS/EGFR, KRAS/BRAF, EGFR/IDH1,
11q13/8p12) survives the switch and is if anything clearer (see Phase 2,
Section 4).

**Gene list for comparison across samples**: built per cancer type (not
globally pooled). Two versions exist: `consensus_genes_per_cancer_type_strict80.parquet`
(the original rule -- ≥80% of that cancer type's own samples tested, ≥100
floor -- kept for comparison) and `consensus_genes_per_cancer_type.parquet`
(a deliberately loose candidate pool -- ≥50 tested samples, no coverage-
fraction requirement -- which Phase 2 actually reads). The 80% single-gene
rule was replaced after a concrete failure case: Melanoma's own list had
only 51 genes under it despite 10,203 patients, because its samples are
split across 87 different panels, so almost no gene reaches 80% of
Melanoma's own population even though famous drivers (BRAF/NRAS/KRAS) are
each tested in 97-100% of it. Reliability is now enforced per-*pair* in
Phase 2 (a joint-coverage gate: enough patients tested for both genes in a
specific pair, ≥50% of that cancer type plus an absolute floor) instead of
per-gene here (see Phase 1, Section 8b, and Phase 2, Section 4).
