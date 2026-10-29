"""
Script Name: 00_export_figure_data.py
Purpose : Export the small tables the R figure scripts need (R cannot read parquet here, so
          every figure reads a CSV). Only reshapes saved Phase 1-4 outputs; no new statistics.
Inputs  : data/processed/*.parquet, results/analyses/predict_partners.csv,
          data/external/pathways/sanchez_vega_2018_pathways.tsv
Outputs : results/figures/data/*.csv (one or more per figure panel)
Note    : run from the repo root, then the R scripts in scripts/figures/.
"""
import json
from pathlib import Path

import numpy as np
from scipy.stats import false_discovery_control
import pandas as pd

P = Path("data/processed")
OUT = Path("results/figures/data"); OUT.mkdir(parents=True, exist_ok=True)
cl = pd.read_parquet(P / "clinical_tidy.parquet")
al = pd.read_parquet(P / "alterations_long.parquet")
cp = pd.read_parquet(P / "comutation_pairs.parquet")
CTS = sorted(cp.Cancer_Type.unique())

# ---- Phase 1: what is handed to Phase 2 (cohort / events, notebook 01 Section 9) ----------
co = pd.read_parquet(P / "cohort.parquet")
ev = pd.read_parquet(P / "events.parquet")
both = co[co.MUT_TESTED & co.CNA_TESTED]
main = both[(both.TMB_GROUP == "non_hypermutated") & both.CANCER_TYPE.isin(CTS)]   # Phase 2's main analysis
print(f"main cohort: {len(main):,} patients in {main.CANCER_TYPE.nunique()} cancer types")
g = both.TMB_GROUP.value_counts()
flow = [("Samples in GENIE v19.0", len(cl), ""),
        ("One sample per patient", len(co), f"{len(cl) - len(co):,} repeat samples of the same patient"),
        ("Tested for mutation and copy number", len(both), f"{len(co) - len(both):,} on panels without copy-number data"),
        ("Not hypermutated", int(g["non_hypermutated"]),
         f"{g['hypermutated']:,} hypermutated (Phase 2c), {g['ultramutated']:,} ultramutated,\n  {g['unknown']:,} TMB unknown"),
        (f"In the {len(CTS)} cancer types with >= 100 patients", len(main),
         f"{int(g['non_hypermutated']) - len(main):,} in smaller cancer types")]
pd.DataFrame(flow, columns=["step", "n", "removed"]).to_csv(OUT / "fig1a_flow.csv", index=False)

top = both[both.CANCER_TYPE.isin(CTS)].CANCER_TYPE.value_counts().index[:15]
(both[both.CANCER_TYPE.isin(top)].groupby(["CANCER_TYPE", "TMB_GROUP"]).size().rename("patients").reset_index()
     .to_csv(OUT / "fig1b_tmb_groups.csv", index=False))

q = json.load(open(P / "qc_summary.json"))
mf, cf, ff = q["mutation_filter_funnel"], q["cna_filter_funnel"], q["fusion_filter_funnel"]
pd.DataFrame([("Mutations", "Raw calls", mf["raw_rows"]), ("Mutations", "Somatic", mf["after_somatic_filter"]),
              ("Mutations", "Protein-altering", mf["after_protein_altering_filter"]),
              ("Mutations", "Not common germline", mf["after_germline_af_filter"]),
              ("Mutations", "Likely pathogenic", mf["after_missense_evidence_filter_FINAL_AM_hotspot_rule"]),
              ("Copy number", "Deep calls (+2 / -2)", cf["deep_calls_(+2/-2)"]),
              ("Copy number", "Direction fits the gene's role", cf["after_dropping_non_consistent"]),
              ("Fusions", "Structural variants", ff["raw_structural_variants"]),
              ("Fusions", "Known driver fusions", ff["driver_fusion_calls"])],
             columns=["type", "step", "n"]).to_csv(OUT / "fig1c_filters.csv", index=False)

em = ev[ev.Sample_ID.isin(main.SAMPLE_ID)].merge(main[["SAMPLE_ID", "CANCER_TYPE"]].rename(columns={"CANCER_TYPE": "ct"}),
                                                  left_on="Sample_ID", right_on="SAMPLE_ID")
em["event"] = np.select([em.Alteration_Type == "MUT", em.CNA_Direction == "Amplification", em.Alteration_Type == "CNV"],
                        ["Mutation", "Amplification", "Deep deletion"], "Fusion")
top_m = main.CANCER_TYPE.value_counts().index[:15]
comp = em.groupby(["ct", "event"]).size().rename("events").reset_index()
comp["patients"] = comp.ct.map(main.CANCER_TYPE.value_counts())
comp.loc[comp.ct.isin(top_m)].to_csv(OUT / "fig1d_composition.csv", index=False)
em.event.value_counts().rename_axis("event").rename("events").reset_index().assign(patients=len(main)) \
  .to_csv(OUT / "fig1d_composition_total.csv", index=False)

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
tb_cov = pd.read_parquet(P / "testability.parquet")
cov = tb_cov.loc[tb_cov.MUT_TESTED, ["SEQ_ASSAY_ID", "Hugo_Symbol"]]   # mutation-tested genes per panel
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
# Simulation with planted links (script 22): precision, recall and precision-recall curves per method
SIM = {"Fisher": ("q_fisher", "p_fisher", "dir_fisher"), "DISCOVER-style": ("q_discover", "p_discover", "dir_discover"),
       "SELECT": ("q_select", "es_select", "dir_select"), "Ours": ("q_cond_subtype", "p_cond_subtype", "dir_cond_subtype")}
sim = pd.read_parquet(P / "simulation_planted.parquet")
pos = sim.truth != "none"
prec_rec, by_strength, pr = [], [], []
for m, (qc, sc, dc) in SIM.items():
    hit = pos & (sim[dc] == sim.truth)                       # a planted pair found in the right direction
    for sc_name, s in sim.groupby("scenario"):
        call = s[qc] < .05; h = hit.loc[s.index]
        prec_rec.append(dict(scenario=sc_name, method=m, calls=int(call.sum()),
                             precision=(call & h).sum() / max(call.sum(), 1), recall=(call & h).sum() / pos.loc[s.index].sum()))
        for (d, f), g in s[pos.loc[s.index]].groupby(["truth", "fold"]):
            by_strength.append(dict(scenario=sc_name, method=m, direction=d, fold=f, n=len(g),
                              recall=((g[qc] < .05) & (g[dc] == g.truth)).mean()))
        score = -s[sc] if sc == "es_select" else s[sc]        # SELECT ranks by effect size, the others by p
        o = np.argsort(score.to_numpy(), kind="stable"); hh = h.to_numpy()[o]
        tp = np.cumsum(hh); k = np.arange(1, len(hh) + 1)
        keep = np.unique(np.linspace(0, len(hh) - 1, 300).astype(int))
        pr += [dict(scenario=sc_name, method=m, recall=tp[i] / hh.sum(), precision=tp[i] / k[i]) for i in keep]
pd.DataFrame(prec_rec).to_csv(OUT / "fig5g_sim_precision_recall.csv", index=False)
pd.DataFrame(by_strength).to_csv(OUT / "fig5h_sim_recall_by_strength.csv", index=False)
pd.DataFrame(pr).to_csv(OUT / "fig5i_sim_pr.csv", index=False)

# Textbook pairs (data/curated): does each method find them in the real benchmark data?
mb = pd.read_parquet(P / "method_benchmark.parquet")
tb = pd.read_csv("data/curated/textbook_pairs.csv")
rows = []
for r in tb.itertuples():
    m = mb[(mb.Cancer_Type == r.Cancer_Type) & (((mb.Gene_A == r.Gene_A) & (mb.Gene_B == r.Gene_B))
                                                | ((mb.Gene_A == r.Gene_B) & (mb.Gene_B == r.Gene_A)))]
    for name, (qc, _, dc) in SIM.items():
        pc = qc.replace("q_", "p_")
        q = mb["q_select"] if pc == "p_select" else pd.Series(false_discovery_control(mb[pc].clip(1e-300, 1)), mb.index)
        found = len(m) > 0 and q[m.index[0]] < .05 and m[dc].iloc[0] == r.expected_direction
        rows.append(dict(Cancer_Type=r.Cancer_Type, pair=f"{r.Gene_A} + {r.Gene_B}", expected=r.expected_direction,
                         method=name, found=bool(found)))
pd.DataFrame(rows).to_csv(OUT / "fig5j_textbook.csv", index=False)
pd.read_parquet(P / "novelty_screen.parquet")[["Cancer_Type", "Gene_A", "Gene_B", "novelty", "selectsim"]] \
  .to_csv(OUT / "fig5f_known.csv", index=False)
pd.read_parquet(P / "mutation_load_confounding.parquet").to_csv(OUT / "fig5e_load.csv", index=False)
pd.read_parquet(P / "mutation_load_confounding_with_hyper.parquet").to_csv(OUT / "fig5e_load_hyper.csv", index=False)
print("exported:", sorted(p.name for p in OUT.glob("*.csv")))
