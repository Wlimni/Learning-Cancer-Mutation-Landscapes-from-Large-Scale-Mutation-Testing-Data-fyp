"""
Script Name: 02_tmb_split_fisher_vs_conditional.py
Purpose : Test Jason's "split by TMB" idea: does Fisher's test stop over-calling
          once hypermutated patients are removed? Compares % significant
          (BH q<0.05) for Fisher on all vs non-hypermutated patients, against
          the conditional test, on the same pairs.
Inputs  : data/processed/{alterations_long,clinical_tidy,panel_gene_coverage,
          cna_gene_panel_coverage,comutation_pairs}.parquet
Outputs : printed table; alteration-burden spread inside the non-hyper group
"""
import numpy as np, pandas as pd
from scipy.stats import fisher_exact
from scipy.stats import false_discovery_control
import sys; sys.path.insert(0, "scripts/validation"); from common import load_phase2_functions
F = load_phase2_functions()
P = "data/processed/"
cl = pd.read_parquet(P + "clinical_tidy.parquet")
cl = cl[cl.IS_REPRESENTATIVE & cl.CNA_TESTED]
al = pd.read_parquet(P + "alterations_long.parquet", columns=["Sample_ID", "Hugo_Symbol", "Alteration_Type"])
mcov = pd.read_parquet(P + "panel_gene_coverage.parquet")[["SEQ_ASSAY_ID", "Hugo_Symbol"]]
ccov = pd.read_parquet(P + "cna_gene_panel_coverage.parquet")
cp = pd.read_parquet(P + "comutation_pairs.parquet")

rows = []
for ct in ["Non-Small Cell Lung Cancer", "Breast Cancer", "Colorectal Cancer", "Endometrial Cancer", "Melanoma", "Bladder Cancer"]:
    pairs = cp[cp.Cancer_Type == ct]
    genes = sorted(set(pairs.Gene_A) | set(pairs.Gene_B))
    c = cl[cl.CANCER_TYPE == ct].set_index("SAMPLE_ID")
    X = pd.crosstab(al.Sample_ID, al.Hugo_Symbol).reindex(index=c.index, columns=genes, fill_value=0) > 0
    # tested = gene on the sample's panel for mutation AND for copy number (as Phase 2)
    both = mcov.merge(ccov).assign(v=True).pivot_table(index="SEQ_ASSAY_ID", columns="Hugo_Symbol", values="v", aggfunc="any")
    T = both.reindex(index=c.SEQ_ASSAY_ID, columns=genes).fillna(False).astype(bool).set_axis(c.index)
    nonhyper = (c.HYPERMUTATION_STATUS == "not_hypermutated").values
    for grp, m in [("all", np.ones(len(c), bool)), ("non-hyper", nonhyper)]:
        ps = []
        for a, b in zip(pairs.Gene_A, pairs.Gene_B):
            t = m & T[a].values & T[b].values
            xa, xb = X[a].values[t], X[b].values[t]
            n11 = (xa & xb).sum(); n10 = (xa & ~xb).sum(); n01 = (~xa & xb).sum(); n00 = t.sum() - n11 - n10 - n01
            ps.append(fisher_exact([[n11, n10], [n01, n00]])[1])
        q = false_discovery_control(ps)
        rows.append((ct, "Fisher " + grp, len(pairs), f"{(q < .05).mean():.1%}"))
    rows.append((ct, "conditional (all, Phase 2)", len(pairs), f"{(pairs.q_value < .05).mean():.1%}"))
    # conditional test on the non-hyper group only, same pairs, same matrix as the Fisher rows
    alt, tst = X.values[nonhyper], T.values[nonhyper]
    M = F["ConditionalModel"](alt, tst); gi = {g: i for i, g in enumerate(genes)}; ps = []
    for a, b in sorted(zip(pairs.Gene_A.map(gi), pairs.Gene_B.map(gi))):
        jt = tst[:, a] & tst[:, b]
        ps.append(F["pair_test"](M.joint(a, b)[jt], int((alt[:, a] & alt[:, b] & jt).sum()))[0])
    rows.append((ct, "conditional (non-hyper)", len(pairs), f"{(false_discovery_control(ps) < .05).mean():.1%}"))
    # burden spread among non-hypermutated patients: TMB does not see copy-number burden
    k = X[nonhyper].sum(axis=1)
    print(f"{ct}: non-hyper alterations/patient median {k.median():.0f}, 90th pct {k.quantile(.9):.0f}, max {k.max()}; "
          f"top 10% carry {k.sort_values(ascending=False).head(max(1, len(k)//10)).sum()/k.sum():.0%} of alterations")
print(pd.DataFrame(rows, columns=["cancer_type", "test", "pairs", "% significant"]).to_string(index=False))
