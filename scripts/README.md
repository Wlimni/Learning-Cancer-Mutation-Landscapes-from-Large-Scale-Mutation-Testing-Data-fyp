# Scripts

Notebooks (`notebooks/01`-`06`) are the main pipeline; these scripts hold the checks and analyses
behind each design decision, the validation, and the figures. Each file opens with a header giving
its purpose, inputs and outputs. **Run everything from the repository root.** Scripts that test the
Phase 2 method load its functions directly from `notebooks/02_phase2_comutation_matrix.ipynb`
(`validation/common.py`), so what is checked is exactly what the notebook runs.

## Run order (full reproduction)

| Step | Run | Produces |
|---|---|---|
| 1 | `notebooks/01_phase1_preprocessing.ipynb` | clean tables (`data/processed/`) |
| 2 | `notebooks/02_phase2_comutation_matrix.ipynb` | Phase 2 pair tables (the longest step) |
| 3 | `scripts/fetch_external.py` | CIViC + pathways (`data/external/`) |
| 4 | `scripts/analyses/08`, `10`, `12` | allele screen, partner models, burden re-test |
| 5 | `notebooks/03_phase3_lookup_tool.ipynb`, `04_phase4_insight.ipynb` | tool and insight results |
| 6 | `scripts/analyses/14`-`18` | validation tables (14 needs the internet) |
| 7 | `scripts/analyses/19`, `20` | novelty screen (19 needs the internet), follow-up |
| 8 | `notebooks/05_validation_and_methods.ipynb`, `06_phase2_deep_dive.ipynb` | read-only summaries |
| 9 | `bash scripts/figures/run_all.sh` | all figures (`results/figures/`) |

Long jobs are best run detached, e.g. `nohup caffeinate -i .venv/bin/python scripts/analyses/15_method_benchmark.py > results/logs/15.log 2>&1 &`.

## `analyses/` -- one question per script

| # | Question | Main output | Discussed in |
|---|---|---|---|
| 01 | Do panels report mutations in every gene they claim to cover? | `results/analyses/` | notebook 01 §9b; review §8 |
| 02 | Does splitting off hypermutated tumours fix Fisher's test? | printed | review §8 (TMB split) |
| 03 | Is the hypermutated group one population (MSI vs POLE)? | printed | review §8 (ultramutated) |
| 04 | Is Fisher calibrated in non-hypermutated patients? | printed | review §8 (null calibration) |
| 05 | Do subtype strata and out-of-centre replication help? | `results/analyses/` | notebook 02 §4b |
| 06 | Which whole Phase 2 design replicates best across centres? | `results/analyses/` | notebook 02 §12 |
| 07 | Is the joint-coverage filter needed? | printed | review §8 (filters) |
| 08 | Do variants of one gene have different partners? Is it smoking? | `allele_*.parquet` | notebook 04 §3 |
| 09 | Which similarity metric matches a new patient best? | `results/analyses/` | notebook 03 §6 |
| 10 | Do per-gene models predict partners better than base rates? | `results/analyses/predict_partners.csv` | notebook 03 §6 |
| 11 | First search for new insight (found the burden confound) | `cnv_burden_confound_check.parquet` | notebook 04 §7 |
| 12 | Systematic burden re-test of copy-number-driven pairs | `burden_adjusted_retest.parquet` | notebook 04 §7 |
| 13 | Should burden strata be built into the main test? | `results/analyses/` | notebook 02 §6d |
| 14 | Do the pairs replicate in TCGA? | `tcga_validation.parquet` | notebook 05 §1 |
| 15 | How does each published method family perform on the same data? | `method_benchmark*.parquet` | notebook 05 §2 |
| 16 | Exclusive gene sets: Dendrix vs load-aware | `exclusive_gene_sets.parquet` | notebook 05 §3 |
| 17 | Does an independent dN/dS route agree? | `conditional_selection_dnds.parquet` | notebook 05 §4 |
| 18 | Does shrinkage across cancer types help? | `shrunk_effects.parquet` | notebook 05 §5 |
| 19 | Which robust pairs are unreported (literature + SelectSim)? | `novelty_screen.parquet` | notebook 04 §9 |
| 20 | Do the unreported pairs survive third-gene, burden and TMB adjustment? | `candidate_followup.parquet` | notebook 04 §9 |
| `p2_variants.py` | Helper: runs Phase 2 design variants with the notebook's own functions (used by 05-08 and 13) | `results/analyses/variants/` | -- |

## `validation/` -- is the Phase 2 test itself correct?

| Script | Checks |
|---|---|
| `build_msk468_matrix.py` | Builds the benchmark (one panel, every gene tested in every sample). **Run first.** |
| `null_calibration.py` | False-positive rate on simulated data with no interactions |
| `tail_ground_truth.py` | Exact Poisson-Binomial tail vs scipy and 60-digit arithmetic |
| `rediscover_run.R`, `rediscover_compare.py` | Same input through Rediscover (R) and through our test |
| `common.py` | Loads the Phase 2 functions from the notebook |

Outputs go to `results/validation/`.

## `figures/` -- publication figures (R)

`run_all.sh` runs `00_export_figure_data.py` (parquet -> small CSVs, since R here cannot read parquet)
and then every `fig*.R`. Shared style is in `_theme.R`. Install the R packages once with
`Rscript scripts/figures/install_r_packages.R`.

## Other

| Script | Purpose |
|---|---|
| `fetch_external.py` | Downloads CIViC evidence and the TCGA pathway gene lists |
| `build_notebook05.py` | Generates `notebooks/05_validation_and_methods.ipynb` |
