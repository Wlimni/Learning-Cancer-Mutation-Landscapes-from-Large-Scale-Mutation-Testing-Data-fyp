"""
Script Name: 09_similarity_metrics.py
Purpose : Choose the similarity metric for matching a new patient to GENIE patients
          (registered Aim: "match a new patient's profile to similar cases"). Task: hide
          one of a patient's real alterations, rank every gene the patient was tested for
          by how often the patient's k nearest neighbours carry it, and score how high the
          hidden gene ranks (AUROC vs the patient's truly unaltered genes).
Metrics : baseline (cancer-type frequency, no similarity); Jaccard; IDF-weighted Jaccard
          (rare alterations count more -- cf. shared-variant similarity scores, Clin Cancer
          Res 2025; TumorComparer); IDF-weighted restricted to the same detailed subtype.
          All are panel-aware: only genes tested in BOTH patients enter the comparison.
Inputs  : data/processed/{clinical_tidy,alterations_long,panel_gene_coverage,cna_gene_panel_coverage}.parquet
Outputs : printed table; results/analyses/similarity_metrics.csv
"""
import numpy as np, pandas as pd
rng = np.random.default_rng(0)
P = "data/processed/"
cl = pd.read_parquet(P + "clinical_tidy.parquet")
cl = cl[cl.IS_REPRESENTATIVE & cl.CNA_TESTED & cl.MUT_TESTED].set_index("SAMPLE_ID")
al = pd.read_parquet(P + "alterations_long.parquet", columns=["Sample_ID", "Hugo_Symbol"]).drop_duplicates()
cov = pd.read_parquet(P + "panel_gene_coverage.parquet")[["SEQ_ASSAY_ID", "Hugo_Symbol"]].merge(
      pd.read_parquet(P + "cna_gene_panel_coverage.parquet"))
K, N_TEST = 50, 800
rows = []
for ct in ["Non-Small Cell Lung Cancer", "Breast Cancer", "Colorectal Cancer", "Glioma", "Melanoma", "Pancreatic Cancer"]:
    c = cl[cl.CANCER_TYPE == ct]
    a = al[al.Sample_ID.isin(c.index)]
    freq = a.Hugo_Symbol.value_counts()
    genes = sorted(freq[freq >= 10].index)                       # genes seen in >= 10 patients
    X = (pd.crosstab(a.Sample_ID, a.Hugo_Symbol).reindex(index=c.index, columns=genes, fill_value=0) > 0).to_numpy()
    T = (cov.assign(v=True).pivot_table(index="SEQ_ASSAY_ID", columns="Hugo_Symbol", values="v", aggfunc="any")
         .reindex(index=c.SEQ_ASSAY_ID, columns=genes).fillna(False).astype(bool).to_numpy())
    X &= T
    rate = X.sum(0) / np.maximum(T.sum(0), 1)                    # per-gene alteration rate among tested
    idf = -np.log(np.clip(rate, 1e-4, 1))                        # rarer alteration -> larger weight
    sub = c.CANCER_TYPE_DETAILED.fillna("NA").to_numpy()
    cand = np.flatnonzero(X.sum(1) >= 3)
    test = rng.choice(cand, size=min(N_TEST, len(cand)), replace=False)
    aucs = {m: [] for m in ["baseline", "jaccard", "idf_jaccard", "idf_jaccard_subtype", "pairwise_logfold"]}
    # Pairwise co-alteration counts among patients tested for both genes (for "pairwise_logfold").
    XT = X.astype(float); TT = T.astype(float)
    n_both = XT.T @ XT                                           # patients with both altered
    n_a_tb = XT.T @ TT                                           # altered in a, tested for b
    for i in test:
        hidden = rng.choice(np.flatnonzero(X[i]))
        q = X[i].copy(); q[hidden] = False                       # the query profile without the hidden gene
        negatives = np.flatnonzero(T[i] & ~X[i])                 # tested, truly unaltered
        if len(negatives) < 5: continue
        others = np.ones(len(X), bool); others[i] = False
        for m in aucs:
            if m == "baseline":
                score = rate.copy()
            elif m == "pairwise_logfold":
                # Each known alteration a shifts gene g by its observed/expected co-alteration ratio
                # (expected = patients with a, tested for g, times g's rate); counts exclude patient i.
                known = np.flatnonzero(q)
                nb_ = n_both[known] - np.outer(X[i, known], X[i])
                ex_ = (n_a_tb[known] - np.outer(X[i, known], T[i])) * rate
                score = np.log(rate + 1e-6) + np.log((nb_ + 1) / (ex_ + 1)).sum(0)
            else:
                w = np.ones(len(genes)) if m == "jaccard" else idf
                inter = (X & q) @ w                               # shared alterations (query side has q only)
                union = (X @ (w * T[i])) + (T @ (w * q)) - inter  # altered in either, tested in both
                sim = np.divide(inter, union, out=np.zeros(len(X)), where=union > 0)
                ok = others.copy()
                if m.endswith("subtype"): ok &= sub == sub[i]
                nb = np.flatnonzero(ok)[np.argsort(-sim[ok])[:K]]
                # neighbour frequency per gene, among neighbours tested for it; shrink to the cancer-type rate
                score = (X[nb].sum(0) + 5 * rate) / (T[nb].sum(0) + 5)
            s_pos, s_neg = score[hidden], score[negatives]
            aucs[m].append((s_pos > s_neg).mean() + 0.5 * (s_pos == s_neg).mean())
    for m, v in aucs.items():
        rows.append(dict(cancer_type=ct, metric=m, n_test=len(v), auroc=np.mean(v)))
    print(ct, {m: round(np.mean(v), 3) for m, v in aucs.items()}, flush=True)
t = pd.DataFrame(rows); t.to_csv("results/analyses/similarity_metrics.csv", index=False)
print(t.pivot(index="cancer_type", columns="metric", values="auroc").round(3).to_string())
print(t.groupby("metric").auroc.mean().round(3).to_string())
