"""
Script Name: 04_null_calibration_nonhyper.py
Purpose : Jason's option "non-hypermutated patients only, then plain Fisher is fine":
          is Fisher's test calibrated there? Builds null data from the REAL non-hyper
          matrix by random swaps that keep every patient's alteration count and every
          gene's frequency (curveball algorithm) but destroy any gene-gene link, then
          counts how often Fisher and the conditional test call a pair significant.
Inputs  : data/processed/*.parquet, notebook test functions (validation/common.py)
Outputs : printed false-positive rates per cancer type
Notes   : patients restricted to those tested for every gene used (no mask needed).
"""
import sys, numpy as np, pandas as pd
from scipy.stats import fisher_exact, false_discovery_control
sys.path.insert(0, "scripts/validation"); from common import load_phase2_functions
F = load_phase2_functions(); rng = np.random.default_rng(1)
P = "data/processed/"
cl = pd.read_parquet(P + "clinical_tidy.parquet"); cl = cl[cl.IS_REPRESENTATIVE & cl.CNA_TESTED]
al = pd.read_parquet(P + "alterations_long.parquet", columns=["Sample_ID", "Hugo_Symbol"]).drop_duplicates()
cov = pd.read_parquet(P + "panel_gene_coverage.parquet")[["SEQ_ASSAY_ID", "Hugo_Symbol"]].merge(
      pd.read_parquet(P + "cna_gene_panel_coverage.parquet"))
cp = pd.read_parquet(P + "comutation_pairs.parquet")


def curveball(X, n_trades):
    rows = [set(np.flatnonzero(r)) for r in X]
    n = len(rows)
    for _ in range(n_trades):
        i, j = rng.integers(n, size=2)
        if i == j: continue
        a, b = rows[i], rows[j]
        only_a, only_b = list(a - b), list(b - a)
        if not only_a or not only_b: continue
        pool = only_a + only_b; rng.shuffle(pool)
        k = len(only_a); both = a & b
        rows[i], rows[j] = both | set(pool[:k]), both | set(pool[k:])
    Y = np.zeros_like(X)
    for r, s in enumerate(rows): Y[r, list(s)] = True
    return Y


def run_tests(X):
    M = F["ConditionalModel"](X, np.ones_like(X)); G = X.shape[1]; pf, pc = [], []
    for a in range(G):
        for b in range(a + 1, G):
            q = M.joint(a, b)
            if q.sum() < 5: continue                  # same outcome-free gate as Phase 2
            n11 = int((X[:, a] & X[:, b]).sum()); n10 = int(X[:, a].sum()) - n11
            n01 = int(X[:, b].sum()) - n11; n00 = len(X) - n11 - n10 - n01
            pf.append(fisher_exact([[n11, n10], [n01, n00]])[1]); pc.append(F["pair_test"](q, n11)[0])
    return np.array(pf), np.array(pc)


print(f"{'cancer type':<28}{'patients':>9}{'genes':>6}{'pairs':>7} | {'Fisher p<.05':>12}{'q<.05':>7} | {'cond p<.05':>10}{'q<.05':>7}")
for ct in ["Non-Small Cell Lung Cancer", "Breast Cancer", "Colorectal Cancer", "Endometrial Cancer",
           "Melanoma", "Bladder Cancer", "Glioma", "Prostate Cancer"]:
    c = cl[(cl.CANCER_TYPE == ct) & (cl.HYPERMUTATION_STATUS == "not_hypermutated")]
    genes = sorted(set(cp.loc[cp.Cancer_Type == ct, "Gene_A"]) | set(cp.loc[cp.Cancer_Type == ct, "Gene_B"]))
    # keep genes and patients so that every kept patient was tested for every kept gene
    per_panel = c.SEQ_ASSAY_ID.value_counts()
    big = per_panel.index[per_panel.cumsum() / per_panel.sum() <= 0.8].tolist() or per_panel.index[:1].tolist()
    common = set.intersection(*[set(cov.loc[cov.SEQ_ASSAY_ID == p, "Hugo_Symbol"]) for p in big])
    genes = [g for g in genes if g in common]
    c = c[c.SEQ_ASSAY_ID.isin(big)]
    X = (pd.crosstab(al.Sample_ID, al.Hugo_Symbol).reindex(index=c.SAMPLE_ID, columns=genes, fill_value=0) > 0).values
    X = X[:, X.sum(0) >= 5]
    fp_f, fp_c, q_f, q_c, npairs = [], [], [], [], 0
    for rep in range(5):
        Y = curveball(X, 20 * len(X))
        pf, pc = run_tests(Y); npairs = len(pf)
        fp_f.append((pf < .05).mean()); fp_c.append((pc < .05).mean())
        q_f.append((false_discovery_control(pf) < .05).mean()); q_c.append((false_discovery_control(pc) < .05).mean())
    print(f"{ct:<28}{len(X):>9,}{X.shape[1]:>6}{npairs:>7} | {np.mean(fp_f):>12.1%}{np.mean(q_f):>7.1%} | "
          f"{np.mean(fp_c):>10.1%}{np.mean(q_c):>7.1%}")
