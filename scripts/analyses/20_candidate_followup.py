"""
Script Name: 20_candidate_followup.py
Purpose : Stress-test the driver-level pairs that the novelty screen (19) found rarely or never
          reported. For each pair, ask whether a third gene, the subtype or copy-number burden
          explains it: logistic regression of gene B on gene A, crude and adjusted for
          (a) detailed subtype, (b) the obvious third-gene drivers, (c) copy-number burden,
          (d) TMB within the non-hypermutated range, and (e) restricted to TMB < 5 -- residual
          microsatellite instability just below the TMB 10 cut-off can create pairs of MSI-target genes.
Inputs  : data/processed/{clinical_tidy,alterations_long,panel_gene_coverage,cna_gene_panel_coverage}.parquet
Outputs : data/processed/candidate_followup.parquet; printed table
Notes   : patients = Phase 2 main cohort (one per patient, non-hypermutated, copy-number tested)
          whose panels cover every gene in the model. Run from the repo root.
"""
import numpy as np
import pandas as pd
import statsmodels.api as sm

P = "data/processed/"
CANDIDATES = [  # cancer type, A, B, third-gene covariates (chosen from known biology of that tumour)
    ("Endometrial Cancer", "CTCF", "PIK3R1", ["PTEN", "PIK3CA", "TP53", "ARID1A", "CTNNB1", "KRAS"]),
    ("Endometrial Cancer", "CTCF", "TP53", ["PTEN", "PIK3CA", "PIK3R1", "ARID1A", "CTNNB1"]),
    ("Endometrial Cancer", "CCNE1", "CTNNB1", ["TP53", "PTEN", "PIK3CA", "ARID1A"]),
    ("Non-Small Cell Lung Cancer", "MDM2", "SMARCA4", ["TP53", "KRAS", "STK11", "KEAP1", "EGFR", "CDK4"]),
    ("Glioma", "CDKN2B", "H3F3A", ["IDH1", "TP53", "ATRX", "EGFR", "PDGFRA"]),
    ("Esophagogastric Cancer", "KRAS", "RARA", ["ERBB2", "TP53", "CDH1"]),          # RARA ~ ERBB2 amplicon?
    ("Colorectal Cancer", "BCL2L1", "SMAD4", ["TP53", "APC", "KRAS", "BRAF", "PIK3CA"]),
]

cl = pd.read_parquet(P + "clinical_tidy.parquet")
al = pd.read_parquet(P + "alterations_long.parquet", columns=["Sample_ID", "Hugo_Symbol", "Alteration_Type"])
mcov = pd.read_parquet(P + "panel_gene_coverage.parquet")[["SEQ_ASSAY_ID", "Hugo_Symbol"]]
ccov = pd.read_parquet(P + "cna_gene_panel_coverage.parquet")
main = cl[cl.IS_REPRESENTATIVE & cl.CNA_TESTED & (cl.HYPERMUTATION_STATUS == "not_hypermutated")]
burden = al[al.Alteration_Type == "CNV"].groupby("Sample_ID").size()


def fit(y, X):
    try:
        r = sm.Logit(y, sm.add_constant(X, has_constant="add")).fit(disp=0, maxiter=200)
        return np.exp(r.params["A"]), r.pvalues["A"]
    except Exception:
        return np.nan, np.nan


rows = []
for ct, A, B, cov in CANDIDATES:
    genes = [A, B] + cov
    m = main[main.CANCER_TYPE == ct]
    # panel must cover every gene for mutations (and for copy number, since altered = MUT or CNV)
    ok_m = mcov[mcov.Hugo_Symbol.isin(genes)].groupby("SEQ_ASSAY_ID").Hugo_Symbol.nunique() == len(genes)
    ok_c = ccov[ccov.Hugo_Symbol.isin(genes)].groupby("SEQ_ASSAY_ID").Hugo_Symbol.nunique() == len(genes)
    panels = set(ok_m[ok_m].index) & set(ok_c[ok_c].index)
    m = m[m.SEQ_ASSAY_ID.isin(panels)]
    X = (al[al.Sample_ID.isin(m.SAMPLE_ID) & al.Hugo_Symbol.isin(genes)].assign(v=1)
         .pivot_table(index="Sample_ID", columns="Hugo_Symbol", values="v", aggfunc="max")
         .reindex(index=m.SAMPLE_ID, columns=genes, fill_value=0).fillna(0))
    det = m.set_index("SAMPLE_ID").CANCER_TYPE_DETAILED.fillna("NA")
    det = det.where(det.map(det.value_counts()) >= 50, "other")
    sub = pd.get_dummies(det, drop_first=True, dtype=float)
    own = X.drop(columns=[A, B]).sum(1)                                   # exclude the pair's own CNVs
    bur = np.log1p(burden.reindex(X.index).fillna(0) - X[[A, B]].sum(1).clip(upper=burden.reindex(X.index).fillna(0)))
    y, a = X[B], X[A].rename("A")
    res = dict(Cancer_Type=ct, A=A, B=B, patients=len(X), n_A=int(a.sum()), n_B=int(y.sum()), both=int((a * y).sum()),
               expected=a.sum() * y.sum() / len(X))
    for name, Z in [("crude", a.to_frame()), ("subtype", pd.concat([a, sub], axis=1)),
                    ("subtype+genes", pd.concat([a, sub, X[cov]], axis=1)),
                    ("subtype+genes+burden", pd.concat([a, sub, X[cov], bur.rename("burden")], axis=1))]:
        res[f"OR_{name}"], res[f"p_{name}"] = fit(y, Z)
    tmb = np.log1p(m.set_index("SAMPLE_ID").TMB_MUT_PER_MB.reindex(X.index)).rename("logTMB")
    full = pd.concat([a, sub, X[cov], bur.rename("burden"), tmb], axis=1)
    res["OR_all+TMB"], res["p_all+TMB"] = fit(y, full)
    low = tmb.index[tmb < np.log1p(5)]
    res["both_TMB<5"] = int((a[low] * y[low]).sum())
    res["OR_TMB<5"], res["p_TMB<5"] = fit(y[low], full.loc[low].drop(columns="logTMB"))
    rows.append(res)

out = pd.DataFrame(rows)
out.to_parquet(P + "candidate_followup.parquet", index=False)
pd.set_option("display.width", 250)
print(out.round(4).to_string(index=False))
