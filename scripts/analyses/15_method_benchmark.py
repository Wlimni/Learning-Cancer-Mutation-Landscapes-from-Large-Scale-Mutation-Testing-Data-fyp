"""
Script Name: 15_method_benchmark.py
Purpose : Run each method family from the literature on the SAME GENIE data and score them
          the same way (docs/literature_review.md, Section 3):
            A  Fisher's exact test (cBioPortal, maftools somaticInteractions)
            B1 DISCOVER-style plug-in Poisson-Binomial (Canisius 2016; Rediscover 2022)
            B2 conditional (Rasch) test, no strata (this project's core test)
            B3 conditional test + subtype strata (this project's Phase 2)
            D  SELECT (Mina 2017), core steps as in its Methods S1: weighted mutual information
               with the per-class sample-count correction (sigmoid w, l = 15, r = 0.1; classes =
               detailed subtypes), an empirical null from 1,000 block-wise randomisations that keep
               every sample's and gene's counts within each subtype (curveball / edge flipping; naive
               column shuffles for blocks with < 30 alterations), effect size = wMI minus the null's
               upper quartile, and Tusher-style empirical FDR per interaction type. Not reproduced:
               separate randomisation per alteration type (every method here sees one gene-level
               'altered' matrix) and the ASC correction for indirect links (a network step)
          (WeSME's margin-preserving permutation is the exact null that B2 computes analytically;
           it is checked separately below on a subset, because permutations are slow.)
Scores  : (1) null calibration -- % of pairs called at p < 0.05 on null data that keep every
              patient's alteration count and every gene's frequency (curveball swaps);
          (2) external replication -- of each method's top-K pairs (same K for all methods),
              the % that replicate in TCGA (p < 0.05, same direction; data/processed/tcga_validation)
Inputs  : data/processed/*.parquet, Phase 2 test code (scripts/validation/common.py)
Outputs : data/processed/method_benchmark.parquet (real data), method_benchmark_null.parquet
          (curveball null p-values, for the calibration figure); printed tables
Notes   : per cancer type, patients restricted to panels covering every gene used, so that all
          methods see an identical, fully tested matrix.
"""
import sys
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, false_discovery_control
sys.path.insert(0, "scripts/validation")
from common import load_phase2_functions

F = load_phase2_functions()
rng = np.random.default_rng(7)
rng_sel = np.random.default_rng(2017)   # SELECT's own randomisations, so the shared curveball null stays as before
P = "data/processed/"
cl = pd.read_parquet(P + "clinical_tidy.parquet")
cl = cl[cl.IS_REPRESENTATIVE & cl.CNA_TESTED & cl.MUT_TESTED & (cl.HYPERMUTATION_STATUS == "not_hypermutated")]
al = pd.read_parquet(P + "alterations_long.parquet", columns=["Sample_ID", "Hugo_Symbol", "Alteration_Type"])
al = al[al.Alteration_Type.isin(["MUT", "CNV"])][["Sample_ID", "Hugo_Symbol"]].drop_duplicates()
cov = pd.read_parquet(P + "panel_gene_coverage.parquet")[["SEQ_ASSAY_ID", "Hugo_Symbol"]].merge(
      pd.read_parquet(P + "cna_gene_panel_coverage.parquet"))
cp = pd.read_parquet(P + "comutation_pairs.parquet")
tv = pd.read_parquet(P + "tcga_validation.parquet")[["Cancer_Type", "Gene_A", "Gene_B", "tcga_p", "tcga_direction"]]
CTS = ["Non-Small Cell Lung Cancer", "Breast Cancer", "Colorectal Cancer", "Glioma", "Endometrial Cancer",
       "Bladder Cancer", "Esophagogastric Cancer", "Ovarian Cancer", "Head and Neck Cancer"]


def curveball(X, n_trades):
    rows = [set(np.flatnonzero(r)) for r in X]; n = len(rows)
    for _ in range(n_trades):
        i, j = rng.integers(n, size=2)
        if i == j: continue
        a, b = rows[i], rows[j]; oa, ob = list(a - b), list(b - a)
        if not oa or not ob: continue
        pool = oa + ob; rng.shuffle(pool); both = a & b
        rows[i], rows[j] = both | set(pool[:len(oa)]), both | set(pool[len(oa):])
    Y = np.zeros_like(X)
    for r, s in enumerate(rows): Y[r, list(s)] = True
    return Y


def curveball_blocks(X, blocks, rounds=10):
    """Margin-preserving randomisation within each block: vectorised curveball trades (each round pairs
    every row with another row of its block). Blocks with < 30 alterations get a naive per-gene shuffle."""
    Y = X.copy()
    for b in np.unique(blocks):
        rows = np.flatnonzero(blocks == b)
        if Y[rows].sum() < 30:
            for g in range(Y.shape[1]):
                Y[rows, g] = rng_sel.permutation(Y[rows, g])
            continue
        for _ in range(rounds):
            r = rng_sel.permutation(rows); r = r[: len(r) // 2 * 2]
            ia, ib = r[0::2], r[1::2]
            A, B = Y[ia], Y[ib]
            D = A ^ B
            ka = (A & ~B).sum(1)                                    # A keeps its number of private genes
            keys = np.where(D, rng_sel.random(D.shape), 2.0)
            rank = np.argsort(np.argsort(keys, axis=1), axis=1)
            toA = D & (rank < ka[:, None])
            both = A & B
            Y[ia], Y[ib] = both | toA, both | (D & ~toA)
    return Y


def wmi(X, blocks, a, b, l=15.0, r=0.1):
    """SELECT's weighted mutual information for gene pairs (a[k], b[k]) (Mina 2017, Methods S1, eq. 3-6)."""
    Xf = X.astype(float); N = len(X)
    n11 = (Xf[:, a] * Xf[:, b]).sum(0); nx, ny = Xf[:, a].sum(0), Xf[:, b].sum(0)
    n10, n01 = nx - n11, ny - n11; n00 = N - n11 - n10 - n01
    Nstar = np.zeros(len(a))
    for blk in np.unique(blocks):
        m = blocks == blk; f = Xf[m].mean(0)
        Nstar += m.sum() / (1 + np.exp(-l * (np.maximum(f[a], f[b]) - r)))
    Nstar = np.maximum(Nstar, n11 + n10 + n01)                    # side condition: N* >= events of x or y
    J = np.stack([n00 - (N - Nstar), n01, n10, n11]) / Nstar    # f*(0,0), f*(0,1), f*(1,0), f*(1,1)
    fx = np.stack([J[0] + J[1], J[2] + J[3]]); fy = np.stack([J[0] + J[2], J[1] + J[3]])
    E = np.stack([fx[0] * fy[0], fx[0] * fy[1], fx[1] * fy[0], fx[1] * fy[1]])
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.nansum(np.where(J > 0, J * np.log(J / E), 0.0), axis=0), n11


def select_test(X, blocks, pairs_idx, n_rand=1000):
    a = np.array([i for i, _ in pairs_idx]); b = np.array([j for _, j in pairs_idx])
    obs, n11 = wmi(X, blocks, a, b)
    null, null11 = np.empty((n_rand, len(a))), np.empty((n_rand, len(a)))
    for k in range(n_rand):
        null[k], null11[k] = wmi(curveball_blocks(X, blocks), blocks, a, b)
    p = (null >= obs).mean(0)                                       # empirical p (SELECT: % random >= observed)
    es = obs - np.quantile(null, 0.75, axis=0)                      # max of the null within its IQR = upper quartile
    direction = np.where(n11 > null11.mean(0), "co-occurring", "exclusive")
    # Tusher empirical FDR, per interaction type: expected false calls at threshold t = mean number of
    # randomised pairs whose own empirical p is <= t
    srt = np.sort(null, axis=0)                                      # p of each randomised score vs the null
    p_null = np.column_stack([1 - np.searchsorted(srt[:, j], null[:, j], side="left") / n_rand for j in range(len(a))])
    q = np.ones(len(a))
    for d in ("co-occurring", "exclusive"):
        m = direction == d
        if not m.any(): continue
        ps = np.sort(p[m]); fdr = np.array([(p_null[:, m] <= t).sum(1).mean() / max((ps <= t).sum(), 1) for t in ps])
        fdr = np.minimum.accumulate(fdr[::-1])[::-1]               # monotone q-values
        q[m] = np.minimum(fdr[np.searchsorted(ps, p[m], side="right") - 1], 1)
    return p, q, es, direction


def methods(X, strata, pairs_idx):
    """p-values (or scores) of every method for the given gene-index pairs."""
    T = np.ones_like(X)
    cond = F["ConditionalModel"](X, T)
    cond_s = F["ConditionalModel"](X, T, strata)
    pi = F["fit_pi"](X, T, None)
    blocks = pd.factorize(strata)[0]
    p_sel, q_sel, es_sel, dir_sel = select_test(X, blocks, pairs_idx)
    out = []
    for a, b in pairs_idx:
        xa, xb = X[:, a], X[:, b]
        n11 = int((xa & xb).sum()); n10 = int(xa.sum()) - n11; n01 = int(xb.sum()) - n11; n00 = len(X) - n11 - n10 - n01
        p_f = fisher_exact([[n11, n10], [n01, n00]])[1]
        dir_f = "co-occurring" if n11 * n00 > n10 * n01 else "exclusive"
        p_d, dir_d, _ = F["discover_pair_test"](pi[:, a], pi[:, b], n11)
        p_c, dir_c, e_c = F["pair_test"](cond.joint(a, b), n11)
        p_s, dir_s, _ = F["pair_test"](cond_s.joint(a, b), n11)
        k = len(out)
        out.append(dict(p_fisher=p_f, dir_fisher=dir_f, p_discover=p_d, dir_discover=dir_d,
                        p_conditional=p_c, dir_conditional=dir_c, p_cond_subtype=p_s, dir_cond_subtype=dir_s,
                        p_select=p_sel[k], q_select=q_sel[k], es_select=es_sel[k], dir_select=dir_sel[k], expected=e_c))
    return out


null_rows, real_rows = [], []
for ct in CTS:
    c = cl[cl.CANCER_TYPE == ct]
    genes = sorted(set(cp.loc[cp.Cancer_Type == ct, "Gene_A"]) | set(cp.loc[cp.Cancer_Type == ct, "Gene_B"]))
    per_panel = c.SEQ_ASSAY_ID.value_counts()
    big = per_panel.index[per_panel.cumsum() / per_panel.sum() <= 0.85].tolist() or per_panel.index[:1].tolist()
    common = set.intersection(*[set(cov.loc[cov.SEQ_ASSAY_ID == p, "Hugo_Symbol"]) for p in big])
    genes = [g for g in genes if g in common]
    c = c[c.SEQ_ASSAY_ID.isin(big)]
    a_c = al[al.Sample_ID.isin(c.SAMPLE_ID) & al.Hugo_Symbol.isin(genes)]          # this cancer type only (fast)
    X = (pd.crosstab(a_c.Sample_ID, a_c.Hugo_Symbol).reindex(index=c.SAMPLE_ID, columns=genes, fill_value=0) > 0).values
    keep = X.sum(0) >= 5; X = X[:, keep]; genes = [g for g, k in zip(genes, keep) if k]
    det = c.CANCER_TYPE_DETAILED.fillna("NA"); vc = det.value_counts()
    strata = det.where(det.map(vc) >= 50, "other").to_numpy()
    gi = {g: i for i, g in enumerate(genes)}
    pr = cp[(cp.Cancer_Type == ct) & cp.Gene_A.isin(gi) & cp.Gene_B.isin(gi)]
    idx = [(gi[a], gi[b]) for a, b in zip(pr.Gene_A, pr.Gene_B)]
    # real data
    res = pd.DataFrame(methods(X, strata, idx))
    res[["Cancer_Type", "Gene_A", "Gene_B"]] = pr[["Cancer_Type", "Gene_A", "Gene_B"]].to_numpy()
    real_rows.append(res)
    # null data: same patient totals and gene frequencies, no gene-gene links
    Y = curveball(X, 5 * len(X))         # 5 swaps per patient mixes the matrix well
    nres = pd.DataFrame(methods(Y, strata, idx))
    nres["Cancer_Type"] = ct
    null_rows.append(nres)
    print(f"{ct:<28} patients {len(X):>6,} genes {len(genes):>3} pairs {len(idx):>4}", flush=True)

real = pd.concat(real_rows, ignore_index=True)
null = pd.concat(null_rows, ignore_index=True)
real = real.merge(tv, on=["Cancer_Type", "Gene_A", "Gene_B"], how="left")
real.to_parquet(P + "method_benchmark.parquet", index=False)
null.to_parquet(P + "method_benchmark_null.parquet", index=False)   # for the calibration figure

M = {"A Fisher (cBioPortal, maftools)": ("p_fisher", "dir_fisher"),
     "B1 DISCOVER-style plug-in": ("p_discover", "dir_discover"),
     "B2 conditional (no strata)": ("p_conditional", "dir_conditional"),
     "B3 conditional + subtype (Phase 2)": ("p_cond_subtype", "dir_cond_subtype"),
     "D SELECT (Mina 2017)": ("p_select", "dir_select")}
QCOL = {"p_select": "q_select"}                         # SELECT uses its own empirical FDR, the others BH


def qvals(df, pc_):
    return df[QCOL[pc_]].to_numpy() if pc_ in QCOL else false_discovery_control(df[pc_].clip(1e-300, 1))


print("\n(1) Null calibration: % of pairs with p < 0.05 when no gene pair interacts (target 5%)")
for name, (pc_, _) in M.items():
    q = qvals(null, pc_)
    print(f"  {name:<36} {(null[pc_] < .05).mean():6.1%}   after FDR: {(q < .05).mean():.1%}")

print("\n(2) Significant calls on real data (q < 0.05) and their replication in TCGA")
t = real[real.tcga_p.notna()]
for name, (pc_, dc) in M.items():
    q = qvals(real, pc_)
    real[f"q_{pc_}"] = q
    s = real[(q < .05) & real.tcga_p.notna()]
    rep = ((s.tcga_p < .05) & (s.tcga_direction == s[dc])).mean()
    opp = ((s.tcga_p < .05) & (s.tcga_direction != s[dc])).sum()
    print(f"  {name:<36} significant {(q < .05).sum():>5}  | testable in TCGA {len(s):>4}  replicate {rep:5.1%}  opposite {opp}")

print("\n(3) Top-K ranking: of each method's K strongest pairs (testable in TCGA), % replicating in TCGA")
ranks = {name: (t[pc_], t[dc]) for name, (pc_, dc) in M.items()}
ranks["D SELECT (effect size)"] = (-t["es_select"], t["dir_select"])        # SELECT ranks by effect size
for name in ("D SELECT (Mina 2017)",):
    ranks.pop(name)
for K in [50, 100, 200]:
    line = f"  K={K:<4}"
    for name, (score, d) in ranks.items():
        top = t.assign(_s=score.to_numpy(), _d=d.to_numpy()).nsmallest(K, "_s")
        line += f" | {name.split(' ')[0]} {((top.tcga_p < .05) & (top.tcga_direction == top._d)).mean():5.1%}"
    print(line)
