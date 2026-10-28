"""
Script Name: 00_export_figure_data.py
Purpose : Export the small tables the R figure scripts need (R cannot read parquet here, so
          every figure reads a CSV). Only reshapes saved Phase 1-4 outputs; no new statistics.
Inputs  : data/processed/*.parquet, results/analyses/predict_partners.csv,
          data/external/pathways/sanchez_vega_2018_pathways.tsv
Outputs : results/figures/data/*.csv (one or more per figure panel)
Note    : run from the repo root, then the R scripts in scripts/figures/.
"""
from pathlib import Path

import numpy as np
import pandas as pd

P = Path("data/processed")
OUT = Path("results/figures/data"); OUT.mkdir(parents=True, exist_ok=True)
cl = pd.read_parquet(P / "clinical_tidy.parquet")
al = pd.read_parquet(P / "alterations_long.parquet")
cp = pd.read_parquet(P / "comutation_pairs.parquet")
CTS = sorted(cp.Cancer_Type.unique())

# Main cohort exactly as Phase 2 uses it: one sample per patient, non-hypermutated,
# copy-number tested, in a qualifying cancer type
rep = cl[cl.IS_REPRESENTATIVE]
main = rep[(rep.HYPERMUTATION_STATUS == "not_hypermutated") & rep.CNA_TESTED & rep.CANCER_TYPE.isin(CTS)]
print(f"main cohort: {len(main):,} patients in {main.CANCER_TYPE.nunique()} cancer types")

# ---- Figure 1 -------------------------------------------------------------
flow = [("Samples in GENIE v19.0", len(cl)),
        ("One sample per patient", len(rep)),
        ("Copy-number tested", int(rep.CNA_TESTED.sum())),
        ("Not hypermutated (TMB < 10)", int((rep.CNA_TESTED & (rep.HYPERMUTATION_STATUS == "not_hypermutated")).sum())),
        (f"In {len(CTS)} cancer types with >= 100 patients", len(main))]
pd.DataFrame(flow, columns=["step", "n"]).to_csv(OUT / "fig1a_flow.csv", index=False)

top_centres = main.CENTER.value_counts().index[:7]
(main.assign(centre=main.CENTER.where(main.CENTER.isin(top_centres), "Other"))
     .groupby(["CANCER_TYPE", "centre"]).size().rename("patients").reset_index()
     .to_csv(OUT / "fig1b_cohort.csv", index=False))

# Panel metadata (what the panel says it tests) vs data (what it actually released)
pg = pd.read_parquet(P / "panel_gene_coverage.parquet")
claims = pg.groupby("SEQ_ASSAY_ID").alteration_types.agg(lambda s: ";".join(s))
released = rep.groupby("SEQ_ASSAY_ID")[["MUT_TESTED", "CNA_TESTED", "SV_TESTED"]].any()
panels = rep.SEQ_ASSAY_ID.value_counts().rename("patients").to_frame().join(released)
panels["claims_MUT"] = claims.reindex(panels.index).fillna("").str.contains("snv")
panels["claims_CNA"] = claims.reindex(panels.index).fillna("").str.contains("cna")
panels["claims_SV"] = claims.reindex(panels.index).fillna("").str.contains("structural")
panels["genes"] = pg.groupby("SEQ_ASSAY_ID").size().reindex(panels.index)
panels.rename_axis("panel").reset_index().to_csv(OUT / "fig1c_panels.csv", index=False)

genes = set(cp.Gene_A) | set(cp.Gene_B)
n_alt = (al[al.Sample_ID.isin(main.SAMPLE_ID) & al.Hugo_Symbol.isin(genes)]
         .groupby("Sample_ID").Hugo_Symbol.nunique())
(main[["SAMPLE_ID", "CANCER_TYPE"]].assign(n_altered=main.SAMPLE_ID.map(n_alt).fillna(0).astype(int))
     .to_csv(OUT / "fig1d_alterations_per_patient.csv", index=False))

# ---- Figure 2 -------------------------------------------------------------
pd.read_parquet(P / "method_benchmark_null.parquet").to_csv(OUT / "fig2a_null_pvalues.csv", index=False)
keep = ["Cancer_Type", "Gene_A", "Gene_B", "n_both_altered", "expected_both", "fold_enrichment",
        "log2_fold_enrichment", "direction", "q_value", "interaction_candidate", "replicated",
        "weak_effect", "disc_q_value", "disc_direction", "rep_p_value",
        "rep_direction"]
cp[keep].to_csv(OUT / "fig2_pairs.csv", index=False)
pd.read_csv(P / "comutation_funnel.csv").to_csv(OUT / "fig2f_funnel.csv", index=False)

# Lung oncoprint: drivers + fusions, patients tested for all three alteration types
LUNG = "Non-Small Cell Lung Cancer"
ONC = ["EGFR", "KRAS", "ALK_fusion", "ROS1_fusion", "RET_fusion", "BRAF", "MET", "ERBB2", "STK11", "KEAP1", "TP53"]
lung = main[(main.CANCER_TYPE == LUNG) & main.MUT_TESTED & main.SV_TESTED]
la = al[al.Sample_ID.isin(lung.SAMPLE_ID) & al.Hugo_Symbol.isin(ONC)].copy()
la["type"] = np.select([la.Alteration_Type == "SV", la.CNA_Direction.eq("Amplification"), la.CNA_Direction.eq("Deep_Deletion")],
                       ["FUSION", "AMP", "DEL"], "MUT")
la = la.groupby(["Sample_ID", "Hugo_Symbol"]).type.agg(lambda s: ";".join(sorted(set(s)))).reset_index()
la.to_csv(OUT / "fig2e_lung_oncoprint.csv", index=False)
pd.DataFrame({"n_patients": [len(lung)]}).to_csv(OUT / "fig2e_lung_n.csv", index=False)
pd.read_parquet(P / "comutation_pairs_fusions.parquet").to_csv(OUT / "fig2e_fusion_pairs.csv", index=False)

# ---- Figure 3 -------------------------------------------------------------
ap = pd.read_parquet(P / "allele_pairs.parquet")
ap[(ap.allele == "KRAS_G12C") & (ap.Cancer_Type == LUNG)].to_csv(OUT / "fig3a_kras_g12c.csv", index=False)
# exact combinations of lung drivers, for the UpSet plot (patients with >= 1 of them)
UPS = ["TP53", "KRAS", "EGFR", "STK11", "KEAP1", "CDKN2A", "SMARCA4", "ALK_fusion"]
ua = al[al.Sample_ID.isin(lung.SAMPLE_ID) & al.Hugo_Symbol.isin(UPS)][["Sample_ID", "Hugo_Symbol"]].drop_duplicates()
ua.assign(v=1).pivot_table(index="Sample_ID", columns="Hugo_Symbol", values="v", fill_value=0) \
  .reindex(columns=UPS, fill_value=0).to_csv(OUT / "fig3b_lung_combinations.csv")
pd.read_csv("results/analyses/predict_partners.csv").to_csv(OUT / "fig3c_predict_partners.csv", index=False)

# Gene frequency in MSK vs other centres (patients whose panel covers the gene)
cov = pg[["SEQ_ASSAY_ID", "Hugo_Symbol"]]
rows = []
for ct in ["Non-Small Cell Lung Cancer", "Breast Cancer", "Colorectal Cancer", "Prostate Cancer",
           "Pancreatic Cancer", "Glioma"]:
    m = main[main.CANCER_TYPE == ct][["SAMPLE_ID", "SEQ_ASSAY_ID", "CENTER"]]
    g_ct = sorted(set(cp.loc[cp.Cancer_Type == ct, "Gene_A"]) | set(cp.loc[cp.Cancer_Type == ct, "Gene_B"]))
    g_ct = [g for g in g_ct if not g.endswith("_fusion")]
    t = m.merge(cov[cov.Hugo_Symbol.isin(g_ct)], on="SEQ_ASSAY_ID")
    hit = al[al.Sample_ID.isin(m.SAMPLE_ID) & al.Hugo_Symbol.isin(g_ct)][["Sample_ID", "Hugo_Symbol"]].drop_duplicates()
    t = t.merge(hit.assign(alt=1), left_on=["SAMPLE_ID", "Hugo_Symbol"], right_on=["Sample_ID", "Hugo_Symbol"], how="left")
    t["group"] = np.where(t.CENTER == "MSK", "MSK", "other")
    f = t.groupby(["Hugo_Symbol", "group"]).alt.agg(lambda s: s.notna().mean()).unstack()
    n = t.groupby(["Hugo_Symbol", "group"]).size().unstack()
    rows.append(f.join(n, rsuffix="_n").assign(Cancer_Type=ct).reset_index())
pd.concat(rows).to_csv(OUT / "fig3d_centre_frequency.csv", index=False)

# ---- Figure 4 -------------------------------------------------------------
pw = pd.read_csv("data/external/pathways/sanchez_vega_2018_pathways.tsv", sep="\t")
pw[["gene", "pathway"]].to_csv(OUT / "fig4a_pathways.csv", index=False)
se = pd.read_parquet(P / "shrunk_effects.parquet")
se[["Cancer_Type", "Gene_A", "Gene_B", "key", "n_both_altered", "expected_both", "q_value", "direction",
    "theta", "var", "theta_shrunk", "tau2"]].to_csv(OUT / "fig4c_shrunk.csv", index=False)
pd.read_parquet(P / "allele_smoking_adjustment.parquet").to_csv(OUT / "fig4d_smoking.csv", index=False)
rc = pd.read_parquet(P / "rare_combinations.parquet")
rc.rename(columns={"civic_combination_evidence": "civic"}).to_csv(OUT / "fig4f_rare.csv", index=False)

# ---- Phase 2b (validation) ----------------------------------------------
tv = pd.read_parquet(P / "tcga_validation.parquet")
tv[tv.tcga_p.notna()][keep + ["tcga_n_both", "tcga_expected", "tcga_direction", "tcga_p"]] \
  .to_csv(OUT / "fig5a_tcga.csv", index=False)
pd.read_parquet(P / "method_benchmark.parquet").to_csv(OUT / "fig5b_benchmark_real.csv", index=False)
pd.read_parquet(P / "novelty_screen.parquet")[["Cancer_Type", "Gene_A", "Gene_B", "novelty", "selectsim"]] \
  .to_csv(OUT / "fig5f_known.csv", index=False)
pd.read_parquet(P / "mutation_load_confounding.parquet").to_csv(OUT / "fig5e_load.csv", index=False)
pd.read_parquet(P / "mutation_load_confounding_with_hyper.parquet").to_csv(OUT / "fig5e_load_hyper.csv", index=False)
print("exported:", sorted(p.name for p in OUT.glob("*.csv")))
