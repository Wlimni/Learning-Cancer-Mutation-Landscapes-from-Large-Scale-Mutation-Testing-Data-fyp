# Scripts

Notebooks (`notebooks/01`-`05`, with `02b`, `02c`) are the main pipeline; these scripts hold the checks and analyses
behind each design decision, the validation, and the figures. Each file opens with a header giving
its purpose, inputs and outputs. **Run everything from the repository root.** Scripts that test the
Phase 2 method load its functions directly from `notebooks/02_phase2_comutation_matrix.ipynb`
(`validation/common.py`), so what is checked is exactly what the notebook runs.

## Run order (full reproduction)

| Step | Run | Produces |
|---|---|---|
| 1 | `notebooks/01_phase1_preprocessing.ipynb` | clean tables (`data/processed/`) |
| 2 | `notebooks/02_phase2_comutation_matrix.ipynb`, then `02c_phase2_extensions.ipynb` | Phase 2 pair table (core), then the side analyses (the longest step) |
| 3 | `scripts/fetch_external.py` | CIViC + pathways (`data/external/`) |
| 4 | `scripts/analyses/08`, `10`, `12` | allele screen, partner models, burden re-test |
| 5 | `notebooks/03_phase3_lookup_tool.ipynb`, `04_phase4_insight.ipynb` | tool and insight results |
| 6 | `scripts/analyses/14`, `15`, `17`, `18`, `21`, `22`, `23` | validation tables (14 needs the internet; 21 and 22 after 15; 23 after 14 and 22; 22 takes 1-2 hours) |
| 7 | `scripts/analyses/19`, `20` | novelty screen (19 needs the internet), follow-up |
| 8 | `notebooks/02b_phase2_validation.ipynb`, `05_phase2_deep_dive.ipynb` | read-only summaries |
| 9 | `bash scripts/figures/run_all.sh` | all figures (`results/figures/`) |

Long jobs are best run detached, e.g. `nohup caffeinate -i .venv/bin/python scripts/analyses/15_method_benchmark.py > results/logs/15.log 2>&1 &`.

## `analyses/` -- one question per script

Two roles. **Pipeline** scripts produce files that the notebooks or figures read, so they must be
re-run when the data change. **Design evidence** scripts were run once to choose a design (for
example, splitting off hypermutated tumours, subtype strata, the filters); the notebooks and the
literature review quote their results, nothing reads their outputs, and their numbers describe the
design at the time they were run (30 Sep - 4 Oct 2026, 21 Oct for 13).

| # | Role | Question | Main output | Discussed in |
|---|---|---|---|---|
| 01 | design evidence | Do panels report mutations in every gene they claim to cover? | `results/analyses/` | notebook 01 §10b; review §8 |
| 02 | design evidence | Does splitting off hypermutated tumours fix Fisher's test? | printed | review §8 (TMB split) |
| 03 | design evidence | Is the hypermutated group one population (MSI vs POLE)? | printed | review §8 (ultramutated) |
| 04 | design evidence | Is Fisher calibrated in non-hypermutated patients? | printed | review §8 (null calibration) |
| 05 | design evidence | Do subtype strata and out-of-centre replication help? | `results/analyses/` | notebook 02, Step 3a |
| 06 | design evidence | Which whole Phase 2 design replicates best across centres? | `results/analyses/` | notebook 02, Discussion |
| 07 | design evidence | Is the joint-coverage filter needed? | printed | review §8 (filters) |
| 08 | pipeline | Do variants of one gene have different partners? Is it smoking? | `allele_*.parquet` | notebook 04 §3 |
| 09 | design evidence | Which similarity metric matches a new patient best? | `results/analyses/` | notebook 03 §6 |
| 10 | pipeline (Fig 6b) | Do per-gene models predict partners better than base rates? | `results/analyses/predict_partners.csv` | notebook 03 §6 |
| 12 | pipeline | Systematic burden re-test of copy-number-driven pairs | `burden_adjusted_retest.parquet` | notebook 04 §7 |
| 13 | design evidence | Should burden strata be built into the main test? | `results/analyses/` | notebook 02c, Section C |
| 14 | pipeline | Do the pairs replicate in TCGA? | `tcga_validation.parquet` | notebook 02c Section G (figure, panel e) |
| 15 | pipeline | How does each published method family perform on the same data? | `method_benchmark*.parquet` | notebook 02c Sections F, G |
| 17 | pipeline | Does an independent dN/dS route agree? | `conditional_selection_dnds.parquet` | notebook 04 (novelty audit, via script 19) |
| 18 | pipeline | Does shrinkage across cancer types help? | `shrunk_effects.parquet` | Phase 4 figure, panel b |
| 19 | pipeline | Which robust pairs are unreported (literature + SelectSim)? | `novelty_screen.parquet` | notebook 04 §9 |
| 20 | pipeline | Do the unreported pairs survive third-gene, burden and TMB adjustment? | `candidate_followup.parquet` | notebook 04 §9 |
| 21 | pipeline | Is any method confounded by mutation load (van de Haar diagnostic)? | `mutation_load_confounding*.parquet` | notebook 02c, Section E |
| 22 | pipeline | Sensitivity and precision of each method: simulation with planted links (3 scenarios) | `simulation_planted.parquet` | notebook 02c Section G |
| 23 | pipeline | Precision, recall and TCGA agreement of each alteration set (SNV, SNV+CNV, SNV+CNV+SV) | `alteration_set_validation.parquet` | notebook 02b |
| `p2_variants.py` | helper | Helper: runs Phase 2 design variants with the notebook's own functions (used by 05-08 and 13) | `results/analyses/variants/` | -- |

## `validation/` -- is the Phase 2 test itself correct?

| Script | Checks |
|---|---|
| `build_msk468_matrix.py` | Builds the benchmark (one panel, every gene tested in every sample). **Run first.** |
| `null_calibration.py` | False-positive rate on simulated data with no interactions |
| `tail_ground_truth.py` | Exact Poisson-Binomial tail vs scipy and 60-digit arithmetic |
| `rediscover_run.R`, `rediscover_compare.py` | Same input through Rediscover (R) and through our test |
| `common.py` | Loads the Phase 2 functions from the notebook |
| `check_phase1_handoff.py` | Checks the three tables Phase 1 hands to Phase 2 (one row per patient, events belong to the cohort, testable on the panel); run by notebook 01, Section 9 |

Outputs go to `results/validation/`.

## `figures/` -- figures by phase (R)

`run_all.sh` runs `00_export_figure_data.py` (parquet -> small CSVs, since R here cannot read parquet)
and then one script per phase: `phase1_data.R`, `phase2_pair_test.R`, `phase3_lookup_tool.R`,
`phase4_biology.R`, plus `phase2b_validation.R`, `phase2c_extensions.R` and `phase2c_methods.R`. Each writes one combined figure for its phase to `results/figures/<phase>/`; the notebooks call these scripts themselves.
Shared style is in `_theme.R`. Install the R packages once with
`Rscript scripts/figures/install_r_packages.R`.

## Other

| Script | Purpose |
|---|---|
| `fetch_external.py` | Downloads CIViC evidence and the TCGA pathway gene lists |
| `build_notebook_02b_validation.py` | Generates `notebooks/02b_phase2_validation.ipynb` |

Retired: `11_novel_insights.py` (first CD79B burden check, replaced by script 12; last version in commit `7811201`).
