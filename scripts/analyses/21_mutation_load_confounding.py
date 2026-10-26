"""
Script Name: 21_mutation_load_confounding.py
Purpose : van de Haar et al. (Cell 2019) diagnostic for "mutation load confounding": genes mutated mostly
          in low-load tumours (negative mutation load association, MLA) collect spurious mutual
          exclusivity. For each method in the benchmark (15), correlate each gene's MLA with the number
          of significant exclusive pairs the method calls for it; a strongly negative correlation means
          the method is confounded by mutation load.
Inputs  : data/processed/{clinical_tidy,alterations_long,panel_gene_coverage,cna_gene_panel_coverage,
          comutation_pairs,method_benchmark}.parquet
Outputs : data/processed/mutation_load_confounding.parquet (one row per gene x cancer type: MLA and
          significant exclusive calls per method); printed summary
Notes   : as in van de Haar's STAR Methods (PMC6816465, supplement S3): each gene's profile is regressed on
          tumour load (+ intercept) with statsmodels Logit, and the score is the fitted load coefficient
          divided by its standard error. MLA uses mutation load (GENIE's TMB, mut/Mb, untransformed); ALA,
          their variant for mutation + copy-number profiles like ours, uses alteration load (each tumour's
          total number of alterations). Patients and genes are exactly those of the benchmark.
"""
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import false_discovery_control, pearsonr, spearmanr

P = "data/processed/"
cl = pd.read_parquet(P + "clinical_tidy.parquet")
cl = cl[cl.IS_REPRESENTATIVE & cl.CNA_TESTED & cl.MUT_TESTED & (cl.HYPERMUTATION_STATUS == "not_hypermutated")]
al = pd.read_parquet(P + "alterations_long.parquet", columns=["Sample_ID", "Hugo_Symbol", "Alteration_Type"])
al = al[al.Alteration_Type.isin(["MUT", "CNV"])][["Sample_ID", "Hugo_Symbol"]].drop_duplicates()
ALT_LOAD = al.groupby("Sample_ID").size()                        # alteration load: all altered genes per tumour
cov = pd.read_parquet(P + "panel_gene_coverage.parquet")[["SEQ_ASSAY_ID", "Hugo_Symbol"]].merge(
      pd.read_parquet(P + "cna_gene_panel_coverage.parquet"))
cp = pd.read_parquet(P + "comutation_pairs.parquet")
mb = pd.read_parquet(P + "method_benchmark.parquet")
M = {"Fisher": "p_fisher", "DISCOVER-style": "p_discover", "SELECT": "p_select",
     "Conditional": "p_conditional", "Conditional + subtype": "p_cond_subtype"}
DIR = {"p_fisher": "dir_fisher", "p_discover": "dir_discover", "p_select": "dir_select",
       "p_conditional": "dir_conditional", "p_cond_subtype": "dir_cond_subtype"}
for name, pc in M.items():                    # q-values exactly as in the benchmark (SELECT: its own empirical FDR)
    mb["q_" + pc] = mb["q_select"] if pc == "p_select" else false_discovery_control(mb[pc].clip(1e-300, 1))

rows = []
for ct in mb.Cancer_Type.unique():
    # same patients and genes as 15_method_benchmark.py
    c = cl[cl.CANCER_TYPE == ct]
    genes = sorted(set(cp.loc[cp.Cancer_Type == ct, "Gene_A"]) | set(cp.loc[cp.Cancer_Type == ct, "Gene_B"]))
    per_panel = c.SEQ_ASSAY_ID.value_counts()
    big = per_panel.index[per_panel.cumsum() / per_panel.sum() <= 0.85].tolist() or per_panel.index[:1].tolist()
    common = set.intersection(*[set(cov.loc[cov.SEQ_ASSAY_ID == p, "Hugo_Symbol"]) for p in big])
    genes = [g for g in genes if g in common]
    c = c[c.SEQ_ASSAY_ID.isin(big)]
    a_c = al[al.Sample_ID.isin(c.SAMPLE_ID) & al.Hugo_Symbol.isin(genes)]
    X = pd.crosstab(a_c.Sample_ID, a_c.Hugo_Symbol).reindex(index=c.SAMPLE_ID, columns=genes, fill_value=0) > 0
    X = X.loc[:, X.sum() >= 5]
    load = c.set_index("SAMPLE_ID").TMB_MUT_PER_MB.reindex(X.index).astype(float)
    aload = ALT_LOAD.reindex(X.index).fillna(0).astype(float)
    ok = load.notna()
    m = mb[mb.Cancer_Type == ct]
    for g in X.columns:
        r = sm.Logit(X.loc[ok, g].astype(float), sm.add_constant(load[ok])).fit(disp=0)
        ra = sm.Logit(X[g].astype(float), sm.add_constant(aload)).fit(disp=0)
        row = dict(Cancer_Type=ct, gene=g, freq=X[g].mean(), MLA=r.tvalues.iloc[1], ALA=ra.tvalues.iloc[1])
        for name, pc in M.items():
            hit = m[((m.Gene_A == g) | (m.Gene_B == g)) & (m["q_" + pc] < 0.05)]
            row["ex_" + name] = int((hit[DIR[pc]] == "exclusive").sum())
            row["co_" + name] = int((hit[DIR[pc]] == "co-occurring").sum())
        rows.append(row)

t = pd.DataFrame(rows)
t["low_MLA"] = t.groupby("Cancer_Type").MLA.transform(lambda s: s <= s.quantile(0.25))
t.to_parquet(P + "mutation_load_confounding.parquet", index=False)
print(f"{len(t)} gene x cancer-type rows, {t.Cancer_Type.nunique()} cancer types; "
      f"MLA median {t.MLA.median():.2f} (range {t.MLA.min():.2f} to {t.MLA.max():.2f})")
print("Lowest-MLA genes:", t.nsmallest(8, "MLA")[["Cancer_Type", "gene", "MLA"]].round(2).to_string(index=False))

# van de Haar's diagnostic: correlation between MLA and number of significant exclusive calls per gene
print("\nMLA vs significant exclusive calls per gene (negative = confounded by mutation load)")
pairs_low = {}
for name in M:
    r, p = pearsonr(t.MLA, t["ex_" + name]); rho, _ = spearmanr(t.MLA, t["ex_" + name])
    r_ala, _ = pearsonr(t.ALA, t["ex_" + name])
    per = t.groupby("Cancer_Type").apply(lambda g: pearsonr(g.MLA, g["ex_" + name])[0] if g["ex_" + name].std() > 0 else np.nan)
    # share of exclusive pairs (counted once) touching a bottom-quartile MLA gene
    m = mb[mb["q_" + M[name]] < 0.05]; m = m[m[DIR[M[name]]] == "exclusive"]
    low = set(zip(t.loc[t.low_MLA, "Cancer_Type"], t.loc[t.low_MLA, "gene"]))
    touch = [(ct, a) in low or (ct, b) in low for ct, a, b in zip(m.Cancer_Type, m.Gene_A, m.Gene_B)]
    print(f"  {name:<22} MLA r = {r:+.2f} (p = {p:.1e}), Spearman {rho:+.2f} | ALA r = {r_ala:+.2f} | per cancer type median r = "
          f"{per.median():+.2f} | exclusive calls {len(m):>4}, touching a low-MLA gene {np.mean(touch) if touch else np.nan:.0%}")
base = np.mean([((ct, a) in low) or ((ct, b) in low) for ct, a, b in zip(mb.Cancer_Type, mb.Gene_A, mb.Gene_B)])
print(f"  (all tested pairs touching a low-MLA gene: {base:.0%})")

# ---------------------------------------------------------------------------------------------
# Part 2. The benchmark above contains only non-hypermutated tumours. Does load confounding appear
# once hypermutated tumours (TMB 10-100) are kept in, and does the conditional test still avoid it?
# Same cancer types, genes and pairs; Fisher, DISCOVER-style and the conditional tests are re-run on
# all tumours with a measured TMB below 100 (SELECT is skipped: 1,000 randomisations per type).
import sys
sys.path.insert(0, "scripts/validation")
from common import load_phase2_functions
from scipy.stats import fisher_exact
F = load_phase2_functions()
cl_all = pd.read_parquet(P + "clinical_tidy.parquet")
cl_all = cl_all[cl_all.IS_REPRESENTATIVE & cl_all.CNA_TESTED & cl_all.MUT_TESTED
                & cl_all.HYPERMUTATION_STATUS.isin(["not_hypermutated", "hypermutated"]) & (cl_all.TMB_MUT_PER_MB < 100)]
M2 = {"Fisher": "p_fisher", "DISCOVER-style": "p_discover", "Conditional": "p_conditional",
      "Conditional + subtype": "p_cond_subtype"}
rows2, genes2 = [], []
for ct in mb.Cancer_Type.unique():
    c = cl_all[cl_all.CANCER_TYPE == ct]
    genes = sorted(set(mb.loc[mb.Cancer_Type == ct, "Gene_A"]) | set(mb.loc[mb.Cancer_Type == ct, "Gene_B"]))
    per_panel = c.SEQ_ASSAY_ID.value_counts()
    big = per_panel.index[per_panel.cumsum() / per_panel.sum() <= 0.85].tolist() or per_panel.index[:1].tolist()
    common = set.intersection(*[set(cov.loc[cov.SEQ_ASSAY_ID == p, "Hugo_Symbol"]) for p in big])
    genes = [g for g in genes if g in common]; c = c[c.SEQ_ASSAY_ID.isin(big)]
    a_c = al[al.Sample_ID.isin(c.SAMPLE_ID) & al.Hugo_Symbol.isin(genes)]
    X = (pd.crosstab(a_c.Sample_ID, a_c.Hugo_Symbol).reindex(index=c.SAMPLE_ID, columns=genes, fill_value=0) > 0).values
    det = c.CANCER_TYPE_DETAILED.fillna("NA"); strata = det.where(det.map(det.value_counts()) >= 50, "other").to_numpy()
    load = c.TMB_MUT_PER_MB.astype(float).to_numpy()
    aload = ALT_LOAD.reindex(c.SAMPLE_ID).fillna(0).astype(float).to_numpy()
    T = np.ones_like(X); cond = F["ConditionalModel"](X, T); cond_s = F["ConditionalModel"](X, T, strata); pi = F["fit_pi"](X, T, None)
    gi = {g: i for i, g in enumerate(genes)}
    pr = mb[(mb.Cancer_Type == ct) & mb.Gene_A.isin(gi) & mb.Gene_B.isin(gi)]
    for A, B in zip(pr.Gene_A, pr.Gene_B):
        a, b = gi[A], gi[B]; xa, xb = X[:, a], X[:, b]
        n11 = int((xa & xb).sum()); n10 = int(xa.sum()) - n11; n01 = int(xb.sum()) - n11; n00 = len(X) - n11 - n10 - n01
        res = dict(Cancer_Type=ct, Gene_A=A, Gene_B=B, p_fisher=fisher_exact([[n11, n10], [n01, n00]])[1],
                   dir_fisher="co-occurring" if n11 * n00 > n10 * n01 else "exclusive")
        res["p_discover"], res["dir_discover"], _ = F["discover_pair_test"](pi[:, a], pi[:, b], n11)
        res["p_conditional"], res["dir_conditional"], _ = F["pair_test"](cond.joint(a, b), n11)
        res["p_cond_subtype"], res["dir_cond_subtype"], _ = F["pair_test"](cond_s.joint(a, b), n11)
        rows2.append(res)
    for g, i in gi.items():
        if X[:, i].sum() >= 5:
            r = sm.Logit(X[:, i].astype(float), sm.add_constant(load)).fit(disp=0)
            ra = sm.Logit(X[:, i].astype(float), sm.add_constant(aload)).fit(disp=0)
            genes2.append(dict(Cancer_Type=ct, gene=g, MLA=r.tvalues[1], ALA=ra.tvalues[1]))
    print(f"{ct:<28} all tumours {len(X):>6,} (hypermutated {(c.HYPERMUTATION_STATUS == 'hypermutated').sum():,})", flush=True)

r2, g2 = pd.DataFrame(rows2), pd.DataFrame(genes2)
g2["low_MLA"] = g2.groupby("Cancer_Type").MLA.transform(lambda s: s <= s.quantile(0.25))
low2 = set(zip(g2.loc[g2.low_MLA, "Cancer_Type"], g2.loc[g2.low_MLA, "gene"]))
print("\nWith hypermutated tumours included: MLA vs significant exclusive calls per gene")
for name, pc in M2.items():
    r2["q"] = false_discovery_control(r2[pc].clip(1e-300, 1))
    ex = r2[(r2.q < .05) & (r2[DIR[pc]] == "exclusive")]
    g2["ex_" + name] = [int(((ex.Cancer_Type == ct) & ((ex.Gene_A == g) | (ex.Gene_B == g))).sum()) for ct, g in zip(g2.Cancer_Type, g2.gene)]
    r, p = pearsonr(g2.MLA, g2["ex_" + name]); r_ala, _ = pearsonr(g2.ALA, g2["ex_" + name])
    touch = np.mean([(ct, a) in low2 or (ct, b) in low2 for ct, a, b in zip(ex.Cancer_Type, ex.Gene_A, ex.Gene_B)]) if len(ex) else np.nan
    print(f"  {name:<22} MLA r = {r:+.2f} (p = {p:.1e}) | ALA r = {r_ala:+.2f} | exclusive calls {len(ex):>4}, touching a low-MLA gene {touch:.0%}")
print(f"  (all tested pairs touching a low-MLA gene: "
      f"{np.mean([(ct, a) in low2 or (ct, b) in low2 for ct, a, b in zip(r2.Cancer_Type, r2.Gene_A, r2.Gene_B)]):.0%})")
g2.to_parquet(P + "mutation_load_confounding_with_hyper.parquet", index=False)
