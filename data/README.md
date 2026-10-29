# Data

Three layers, following the usual raw -> external -> processed convention (plus a small `curated/` folder): raw files are never
edited, every processed file is produced by a notebook or script, and anything can be rebuilt by
re-running them in order (see `scripts/README.md` for the run order).

| Folder | What | In git? |
|---|---|---|
| `raw/` | AACR Project GENIE v19.0 public release, unmodified | No (GENIE terms forbid redistribution) |
| `external/` | Public annotation resources and external cohorts | No (re-downloadable) |
| `curated/` | Hand-curated reference lists with literature sources (`textbook_pairs.csv`: 29 pairs with a published direction, chosen before looking at our results; `top_pairs_biology.csv`: biology notes for the strongest pairs) | Yes |
| `processed/` | Tables produced by this project | Aggregate tables only; individual-level tables (clinical_tidy, alterations_long, cohort, events, putative_biallelic_tsg_events) and the AlphaMissense extract are not tracked |

## raw/ -- AACR Project GENIE v19.0

Request access at [AACR Project GENIE](https://www.aacr.org/professionals/research/aacr-project-genie/)
(Synapse) and place the release files here unchanged.

| File | Content |
|---|---|
| `data_mutations_extended.txt` | Somatic mutation calls (MAF), 3.46M rows |
| `data_CNA.txt` | Gene-level copy number (-2 ... +2; NA = not tested), 1,004 genes x 172,875 samples |
| `data_sv.txt` | Structural variants / fusions, 56,356 rows |
| `data_clinical_sample.txt`, `data_clinical_patient.txt` | Cancer type, panel, sample type; sex, vital status |
| `genomic_information.txt` | Genomic regions covered by each panel |
| `assay_information.txt` | Panel metadata (claimed alteration types) -- not trusted, see README "Problems found" |
| `tmb_19.0-public.tsv` | Tumour mutational burden per sample |
| `data_gene_matrix.txt`, `data_cna_hg19.seg`, `gene_panels/` | Release files not used directly |

## external/ -- annotation resources and external cohorts

| File / folder | Source | Used for | Fetched |
|---|---|---|---|
| `AlphaMissense_aa_substitutions.tsv.gz` | [AlphaMissense](https://storage.googleapis.com/dm_alphamissense/AlphaMissense_aa_substitutions.tsv.gz) (Cheng 2023) | Missense pathogenicity filter (Phase 1) | Jul 2026 |
| `hotspots_v2.xls` | [cancerhotspots.org](https://www.cancerhotspots.org/files/hotspots_v2.xls) v2 | Hotspot rescue of missense calls (Phase 1) | Jul 2026 |
| `hgnc_complete_set.txt` | [HGNC](https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt) | Gene symbol harmonisation | Jul 2026 |
| `oncokb_curated_genes.json` | [OncoKB](https://www.oncokb.org/api/v1/utils/allCuratedGenes) public gene list | Oncogene / tumour-suppressor roles for copy-number direction | Jul 2026 |
| `civic/` | [CIViC](https://civicdb.org) nightly summaries | Treatment evidence in the lookup tool (Phase 3) | Oct 2026 |
| `pathways/sanchez_vega_2018_pathways.tsv` | Sanchez-Vega 2018 pathways via [PathwayMapper](https://github.com/iVis-at-Bilkent/pathway-mapper) | Pathway structure (Phase 4) | Oct 2026 |
| `tcga/` | TCGA PanCancer Atlas via the public cBioPortal REST API (cache) | External validation (script 14) | Oct 2026 |
| `selectsim/supp_tables.xlsx` | Supplementary tables of Iyer et al., *Nat Genet* 2026 (SelectSim) | Novelty screen (script 19) | Oct 2026 |
| `europepmc_counts.json` | [Europe PMC](https://europepmc.org) REST API (cache of paper counts) | Novelty screen (script 19) | Oct 2026 |

`scripts/fetch_external.py` downloads CIViC and the pathway list; scripts 14 and 19 fill their own
caches; the Phase 1 notebook documents the other downloads.

## processed/ -- outputs of this project

Patients are one sample per patient; "pairs" are gene pairs within one cancer type unless stated.
Phase 2 pair tables share a core set of columns: `n_tested_both`, `n_both_altered`, `expected_both`,
`fold_enrichment` (observed / expected), `direction`, `p_value`, `q_value` (BH-FDR per cancer type).

### Phase 1 -- clean data (`notebooks/01_phase1_preprocessing.ipynb`)

| File | One row per | Notes |
|---|---|---|
| `clinical_tidy.parquet` | sample (271,837) | Clinical fields + `IS_REPRESENTATIVE` (one per patient), `MUT_TESTED` / `CNA_TESTED` / `SV_TESTED` (read from the data), TMB, `HYPERMUTATION_STATUS` |
| `alterations_long.parquet` | alteration (1,967,371) | Likely-harmful mutations (incl. `TERT_promoter` hotspots and MET exon 14 splice changes), deep copy-number calls, driver fusions (`{gene}_fusion`) |
| `panel_gene_coverage.parquet` | panel x gene, mutations (52,632) | Which genes each panel sequences (exon probes), plus `TERT_promoter` where the promoter is tested |
| `cna_gene_panel_coverage.parquet` | panel x gene, copy number (11,680) | Read from `data_CNA.txt`, not metadata |
| `sv_gene_panel_coverage.parquet` | panel x fusion gene (6,244) | Fusion genes on panels that released fusions (any probe) |
| **`cohort.parquet`** | patient (227,696) | **Hand-off to Phase 2** (Section 5b): one sample per patient, clinical and testing fields, `TMB_GROUP` (non_hypermutated / hypermutated / ultramutated / unknown) |
| **`events.parquet`** | driver event (1,723,487) | **Hand-off to Phase 2:** the driver events of those patients (same columns as `alterations_long`) |
| **`testability.parquet`** | panel x gene feature (62,250) | **Hand-off to Phase 2:** `MUT_TESTED` / `CNA_TESTED` / `SV_TESTED` in one table; `FEATURE_KIND` = gene / fusion / promoter |
| `putative_biallelic_tsg_events.parquet` | event (24,193) | Shallow deletion + damaging mutation in a tumour suppressor (kept separate) |
| `qc_summary.json` | -- | Filter funnel and data-quality counts quoted in the README |
| `alphamissense_filtered.parquet` (+ `.ids.json`) | protein variant | AlphaMissense extract for GENIE proteins; large, not in git |

### Phase 2 -- pair tests (`notebooks/02_phase2_comutation_matrix.ipynb`; side tables from `02c_phase2_extensions.ipynb`)

| File | One row per | Notes |
|---|---|---|
| `comutation_pairs.parquet` | pair (6,871) | **Main result.** Adds replication (`rep_*`, `replicated`), burden flag (`burden_*`, added by notebook 02c), mechanism and proximity columns, `interaction_candidate`, `weak_effect` |
| `comutation_funnel.csv` | filter step | Tested -> significant -> candidates -> replicated |
| `comutation_pairs_hypermutated.parquet` | pair (99,429) | Hypermutated patients (TMB 10-100), descriptive  (notebook 02c) |
| `comutation_pairs_mutation_only.parquet` | pair (2,554) | Mutations only, all mutation-tested patients  (notebook 02c) |
| `comutation_pairs_fusions.parquet` | fusion x gene pair (375) | Fusion screen  (notebook 02c) |
| `mechanism_pairs.parquet` | pair x mechanism (8,661) | Mutation vs copy number on each side |
| `hypermutation_status_by_cancer_type.parquet` | cancer type (56) | Counts per TMB group  (notebook 02c) |

### Phase 2 deep dive (`notebooks/05_phase2_deep_dive.ipynb`)

| File | One row per | Notes |
|---|---|---|
| `mechanism_pairs_checked.parquet` | cross-mechanism hit (687) | Subtype, burden, proximity and passenger checks |

### Phase 4 -- insight (`notebooks/04_phase4_insight.ipynb` and scripts)

| File | Written by | Notes |
|---|---|---|
| `cancer_type_heterogeneity.parquet` | notebook 04 | Pair heterogeneity across cancer types (Cochran's Q); pattern = universal / specific / opposite |
| `rare_combinations.parquet` | notebook 04 | Replicated exclusive pairs, how rare, CIViC combination evidence |
| `allele_pairs.parquet`, `allele_heterogeneity.parquet`, `allele_smoking_adjustment.parquet` | script 08 | Variant-level partners; smoking-signature adjustment |
| `burden_adjusted_retest.parquet` | script 12 | 25 copy-number-driven pairs re-tested with burden |
| `novelty_screen.parquet` | script 19 | 387 robust pairs vs Europe PMC and the SelectSim atlas |
| `candidate_followup.parquet` | script 20 | Unreported candidates adjusted for subtype, third genes, burden, TMB |

### Validation (`scripts/analyses/14-22`; summarised in `notebooks/02b_phase2_validation.ipynb` and `02c_phase2_extensions.ipynb`)

| File | Written by | Notes |
|---|---|---|
| `tcga_validation.parquet` | script 14 | Phase 2 pairs + TCGA re-test (`tcga_*` columns) |
| `method_benchmark.parquet`, `method_benchmark_null.parquet` | script 15 | p-values of each method family on real and curveball-null data |
| `conditional_selection_dnds.parquet` | script 17 | Protein-changing : silent ratios by partner status |
| `shrunk_effects.parquet` | script 18 | Effects shrunk across cancer types (DerSimonian-Laird) |
| `simulation_planted.parquet` | script 22 | Per pair, scenario and replicate: planted truth and every method's p, q and direction |
| `mutation_load_confounding.parquet`, `..._with_hyper.parquet` | script 21 | Per gene: mutation load association and significant exclusive calls per method |
