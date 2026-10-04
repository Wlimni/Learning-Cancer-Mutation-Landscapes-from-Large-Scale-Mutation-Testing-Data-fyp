"""
Script Name: 10_predict_partners.py
Purpose : Can a patient's known alterations predict their other alterations better than
          the cancer-type base rate? (the useful core of "matching a new patient to
          similar cases"). Same held-out task as 09_similarity_metrics.py, but on an
          80/20 patient split, comparing:
            baseline   -- each gene's alteration rate in the training patients
            subtype    -- the rate within the patient's detailed subtype
            logistic   -- one L2-regularised logistic model per gene: other genes'
                          alteration status + alteration load + subtype
Inputs  : data/processed/*.parquet
Outputs : printed table; results/analyses/predict_partners.csv
Notes   : untested genes are 0 in the inputs and excluded from training targets.
"""
import numpy as np, pandas as pd
from sklearn.linear_model import LogisticRegression
rng = np.random.default_rng(0)
P = "data/processed/"
cl = pd.read_parquet(P + "clinical_tidy.parquet")
cl = cl[cl.IS_REPRESENTATIVE & cl.CNA_TESTED & cl.MUT_TESTED & (cl.HYPERMUTATION_STATUS == "not_hypermutated")].set_index("SAMPLE_ID")
al = pd.read_parquet(P + "alterations_long.parquet", columns=["Sample_ID", "Hugo_Symbol"]).drop_duplicates()
cov = pd.read_parquet(P + "panel_gene_coverage.parquet")[["SEQ_ASSAY_ID", "Hugo_Symbol"]].merge(
      pd.read_parquet(P + "cna_gene_panel_coverage.parquet"))
rows = []
for ct in ["Non-Small Cell Lung Cancer", "Breast Cancer", "Colorectal Cancer", "Glioma", "Melanoma", "Pancreatic Cancer"]:
    c = cl[cl.CANCER_TYPE == ct]
    a = al[al.Sample_ID.isin(c.index)]
    freq = a.Hugo_Symbol.value_counts(); genes = sorted(freq[freq >= 20].index)
    X = (pd.crosstab(a.Sample_ID, a.Hugo_Symbol).reindex(index=c.index, columns=genes, fill_value=0) > 0).to_numpy()
    T = (cov.assign(v=True).pivot_table(index="SEQ_ASSAY_ID", columns="Hugo_Symbol", values="v", aggfunc="any")
         .reindex(index=c.SEQ_ASSAY_ID, columns=genes).fillna(False).astype(bool).to_numpy())
    X &= T
    det = c.CANCER_TYPE_DETAILED.fillna("NA"); top = det.value_counts(); det = det.where(det.map(top) >= 50, "other").to_numpy()
    S = pd.get_dummies(det).to_numpy(float)
    train = rng.random(len(X)) < 0.8; test = np.flatnonzero(~train & (X.sum(1) >= 3))
    rate = X[train].sum(0) / np.maximum(T[train].sum(0), 1)
    sub_rate = {s: X[train & (det == s)].sum(0) / np.maximum(T[train & (det == s)].sum(0), 1) for s in np.unique(det)}
    # Held-out query: hide one real alteration per test patient (fixed draw for all methods).
    hidden = {i: rng.choice(np.flatnonzero(X[i])) for i in test}
    Q = X.astype(float).copy()
    for i, h in hidden.items(): Q[i, h] = 0
    load = np.log1p(Q.sum(1, keepdims=True))
    pred = np.zeros((len(X), len(genes)))
    for g in range(len(genes)):
        tr = train & T[:, g]
        if X[tr, g].sum() < 5: pred[:, g] = rate[g]; continue
        feats = np.hstack([np.delete(Q, g, axis=1), load, S])
        m = LogisticRegression(C=0.5, max_iter=300).fit(feats[tr], X[tr, g])
        pred[:, g] = m.predict_proba(feats)[:, 1]
    res = {"baseline": [], "subtype": [], "logistic": []}
    for i, h in hidden.items():
        neg = np.flatnonzero(T[i] & ~X[i])
        if len(neg) < 5: continue
        for name, score in [("baseline", rate), ("subtype", sub_rate[det[i]]), ("logistic", pred[i])]:
            res[name].append((score[h] > score[neg]).mean() + 0.5 * (score[h] == score[neg]).mean())
    for k, v in res.items(): rows.append(dict(cancer_type=ct, method=k, n_test=len(v), auroc=np.mean(v)))
    print(ct, {k: round(np.mean(v), 3) for k, v in res.items()}, flush=True)
t = pd.DataFrame(rows); t.to_csv("results/analyses/predict_partners.csv", index=False)
print(t.groupby("method").auroc.mean().round(3).to_string())
