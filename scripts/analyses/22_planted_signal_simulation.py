"""
Script Name: 22_planted_signal_simulation.py
Purpose : Sensitivity and precision of each method, on data where the truth is known. Background
          data with no gene-gene links are generated three ways, a few pairs get a planted link of
          known strength, and every method of the benchmark (script 15) is run on the same data:
            S1 model-based  -- each cell drawn from the fitted gene + patient rate model within
                               subtypes (the null our test assumes; a sanity check that favours us)
            S2 real data    -- the real matrix, shuffled within subtypes so every patient keeps its
                               alteration count and every gene its frequency (curveball); realistic
            S3 misspecified -- S1 plus what our test does not model: a hidden patient group that
                               shifts gene rates, and 5% hypermutated-like patients with extra,
                               gene-specific alterations (a stress test)
          Planted links: per cancer type, 4 pairs with no gene in common (true links are rare in real data;
          planting many makes most genes part of a link and creates knock-on associations); co-occurring at
          1.5x / 2x / 3x the background count of patients with both, or exclusive at 1/1.5 / 1/2 / 1/3.
          Gene B's alterations are moved between patients with and without gene A, so B's frequency
          is unchanged. All other tested pairs are true negatives.
Inputs  : as scripts/analyses/15_method_benchmark.py (same patients, genes, pairs and method code)
Outputs : data/processed/simulation_planted.parquet -- one row per pair, scenario and replicate:
          truth, planted effect and every method's p, q and direction (pair-level only, no patients)
Notes   : 3 scenarios x 6 replicates x 9 cancer types; SELECT uses 200 randomisations (1,000 in
          script 15) to keep the run to a few hours. q-values: Benjamini-Hochberg per cancer type
          (as in Phase 2); SELECT its own empirical FDR. Seeded.
"""
import time
import numpy as np
import pandas as pd
from scipy.special import expit, logit
from scipy.stats import false_discovery_control

# Benchmark set-up and method code, executed from script 15 itself so both scripts test the same thing
src = open("scripts/analyses/15_method_benchmark.py").read()
ns = {"__name__": "benchmark"}
exec(src[:src.index("null_rows, real_rows = [], []")], ns)
ns["select_test"].__defaults__ = (200,)          # fewer SELECT randomisations (see Notes)
cl, al, cov, cp, CTS, F = (ns[k] for k in ("cl", "al", "cov", "cp", "CTS", "F"))
methods, curveball_blocks = ns["methods"], ns["curveball_blocks"]

rng = np.random.default_rng(22)
n_planted = 0                                      # running count, so all effect sizes are used equally
EFFECTS = [("co-occurring", 1.5), ("co-occurring", 2.0), ("co-occurring", 3.0),
           ("exclusive", 1.5), ("exclusive", 2.0), ("exclusive", 3.0)]
N_REP, MAX_PLANTED, OUT = 6, 4, "data/processed/simulation_planted.parquet"
PCOLS = ["p_fisher", "p_discover", "p_conditional", "p_cond_subtype"]


def matrix(ct):
    """Same patients, genes and pairs as script 15 for this cancer type."""
    c = cl[cl.CANCER_TYPE == ct]
    genes = sorted(set(cp.loc[cp.Cancer_Type == ct, "Gene_A"]) | set(cp.loc[cp.Cancer_Type == ct, "Gene_B"]))
    per_panel = c.SEQ_ASSAY_ID.value_counts()
    big = per_panel.index[per_panel.cumsum() / per_panel.sum() <= 0.85].tolist() or per_panel.index[:1].tolist()
    common = set.intersection(*[set(cov.loc[cov.SEQ_ASSAY_ID == p, "Hugo_Symbol"]) for p in big])
    genes = [g for g in genes if g in common]
    c = c[c.SEQ_ASSAY_ID.isin(big)]
    a_c = al[al.Sample_ID.isin(c.SAMPLE_ID) & al.Hugo_Symbol.isin(genes)]
    X = (pd.crosstab(a_c.Sample_ID, a_c.Hugo_Symbol).reindex(index=c.SAMPLE_ID, columns=genes, fill_value=0) > 0).values
    keep = X.sum(0) >= 5; X = X[:, keep]; genes = [g for g, k in zip(genes, keep) if k]
    det = c.CANCER_TYPE_DETAILED.fillna("NA"); vc = det.value_counts()
    strata = det.where(det.map(vc) >= 50, "other").to_numpy()
    gi = {g: i for i, g in enumerate(genes)}
    pr = cp[(cp.Cancer_Type == ct) & cp.Gene_A.isin(gi) & cp.Gene_B.isin(gi)]
    return X, strata, genes, [(gi[a], gi[b]) for a, b in zip(pr.Gene_A, pr.Gene_B)]


def background(X, strata, scenario):
    """Data with no planted links (see Purpose)."""
    if scenario == "S2 real data, links removed":
        return curveball_blocks(X, pd.factorize(strata)[0])
    pi = np.clip(F["fit_pi"](X, np.ones_like(X), strata), 1e-6, 1 - 1e-6)
    if scenario == "S1 model-based":
        return rng.random(X.shape) < pi
    n, g = X.shape                                                  # S3: what our test does not model
    z = rng.choice([-1.0, 1.0], n)                                  # hidden patient group
    hyper = rng.random(n) < 0.05                                    # hypermutated-like patients
    lp = logit(pi) + np.outer(z, rng.normal(0, 0.7, g)) + np.outer(hyper, 1.5 + rng.normal(0, 0.7, g))
    return rng.random(X.shape) < expit(lp)


def plant(Y, pairs_idx):
    """Plant links in pairs with no gene in common; returns the new matrix and the truth per pair."""
    global n_planted
    Y = Y.copy(); used = set(); truth = {}
    freq = Y.mean(0)
    for k in rng.permutation(len(pairs_idx)):
        a, b = pairs_idx[k]
        if a in used or b in used or min(freq[a], freq[b]) < 0.03 or len(truth) >= MAX_PLANTED:
            continue
        direction, fold = EFFECTS[n_planted % len(EFFECTS)]
        xa, xb = Y[:, a], Y[:, b]
        n11 = int((xa & xb).sum())
        target = round(n11 * fold) if direction == "co-occurring" else round(n11 / fold)
        gain = np.flatnonzero(xa & ~xb) if direction == "co-occurring" else np.flatnonzero(xa & xb)
        lose = np.flatnonzero(~xa & xb) if direction == "co-occurring" else np.flatnonzero(~xa & ~xb)
        m = min(abs(target - n11), len(gain), len(lose))
        if m == 0:
            continue
        on, off = rng.choice(gain, m, replace=False), rng.choice(lose, m, replace=False)
        if direction == "co-occurring":
            Y[on, b], Y[off, b] = True, False                       # B moves into A-carriers ...
        else:
            Y[on, b], Y[off, b] = False, True                       # ... or out of them; B's total is unchanged
        used |= {a, b}; truth[k] = (direction, fold); n_planted += 1
    return Y, truth


rows, t0 = [], time.time()
for ct in CTS:
    X, strata, genes, idx = matrix(ct)
    for scenario in ["S1 model-based", "S2 real data, links removed", "S3 misspecified"]:
        for rep in range(N_REP):
            Y, truth = plant(background(X, strata, scenario), idx)
            res = pd.DataFrame(methods(Y, strata, idx))
            for pc in PCOLS:                                        # BH per cancer type, as in Phase 2
                res["q" + pc[1:]] = false_discovery_control(res[pc].clip(1e-300, 1))
            res["Gene_A"] = [genes[a] for a, _ in idx]; res["Gene_B"] = [genes[b] for _, b in idx]
            res["truth"] = [truth.get(k, ("none", 1.0))[0] for k in range(len(idx))]
            res["fold"] = [truth.get(k, ("none", 1.0))[1] for k in range(len(idx))]
            res["Cancer_Type"], res["scenario"], res["rep"] = ct, scenario, rep
            rows.append(res)
            print(f"{ct:<28} {scenario:<28} rep {rep}  planted {len(truth):>2}  {time.time() - t0:7.0f}s", flush=True)
    pd.concat(rows, ignore_index=True).to_parquet(OUT, index=False)     # save as we go

sim = pd.concat(rows, ignore_index=True)
sim.to_parquet(OUT, index=False)

# Quick summary (the notebook draws the full picture): power and false discovery rate at q < 0.05
M = {"Fisher": ("q_fisher", "dir_fisher"), "DISCOVER-style": ("q_discover", "dir_discover"),
     "SELECT": ("q_select", "dir_select"), "Ours": ("q_cond_subtype", "dir_cond_subtype")}
for scenario, s in sim.groupby("scenario"):
    print(f"\n{scenario}: {(s.truth != 'none').sum()} planted, {(s.truth == 'none').sum()} negatives")
    for name, (qc, dc) in M.items():
        call = s[qc] < .05
        tp = call & (s.truth != "none") & (s[dc] == s.truth)
        print(f"  {name:<15} power {tp.sum() / (s.truth != 'none').sum():5.1%}   "
              f"false discovery rate {(call & ~tp).sum() / max(call.sum(), 1):5.1%}   calls {call.sum()}")
