"""
Script Name: 18_hierarchical_shrinkage.py
Purpose : Family F (Scharpf et al. 2022, GENIE): a Bayesian hierarchical model lets each cancer
          type's estimate of a pair's effect borrow strength from the same pair in other cancer
          types. Borrowed as empirical-Bayes shrinkage (DerSimonian-Laird random effects): each
          cancer type's log(observed/expected) is pulled towards the pair's pooled effect, more
          strongly when that cancer type's estimate is noisy.
Test    : do shrunk effects agree better with TCGA (independent patients) than raw effects?
          Measured by correlation of effect sizes and direction agreement, overall and for the
          noisiest estimates (few co-altered patients), where borrowing should help most.
Inputs  : data/processed/tcga_validation.parquet (Phase 2 pairs + TCGA re-test)
Outputs : data/processed/shrunk_effects.parquet; printed comparison
"""
import numpy as np
import pandas as pd

m = pd.read_parquet("data/processed/tcga_validation.parquet")
m["theta"] = np.log((m.n_both_altered + 0.5) / (m.expected_both + 0.5))
m["var"] = 1 / (m.n_both_altered + 0.5)
m["key"] = m.Gene_A + "/" + m.Gene_B
out = []
for k, g in m.groupby("key"):
    if len(g) < 2:
        g = g.assign(theta_shrunk=g.theta, tau2=np.nan); out.append(g); continue
    w = 1 / g["var"]; mu_fe = (w * g.theta).sum() / w.sum()
    Q = (w * (g.theta - mu_fe) ** 2).sum(); df = len(g) - 1
    tau2 = max(0.0, (Q - df) / (w.sum() - (w ** 2).sum() / w.sum()))        # DerSimonian-Laird
    w_re = 1 / (g["var"] + tau2); mu = (w_re * g.theta).sum() / w_re.sum()
    shrink = g["var"] / (g["var"] + tau2) if tau2 > 0 else pd.Series(1.0, index=g.index)
    out.append(g.assign(theta_shrunk=shrink * mu + (1 - shrink) * g.theta, tau2=tau2))
s = pd.concat(out)
s.to_parquet("data/processed/shrunk_effects.parquet", index=False)

t = s[s.tcga_p.notna()].copy()
t["theta_tcga"] = np.log((t.tcga_n_both + 0.5) / (t.tcga_expected + 0.5))
t["tcga_sig"] = t.tcga_p < 0.05
for name, sel in [("all pairs testable in TCGA", t.index == t.index),
                  ("noisy GENIE estimates (< 30 co-altered patients)", t.n_both_altered < 30),
                  ("pairs tested in >= 3 GENIE cancer types", t.key.map(m.key.value_counts()) >= 3)]:
    d = t[sel]
    r_raw, r_sh = np.corrcoef(d.theta, d.theta_tcga)[0, 1], np.corrcoef(d.theta_shrunk, d.theta_tcga)[0, 1]
    sig = d[d.tcga_sig]
    agree_raw = (np.sign(sig.theta) == np.sign(sig.theta_tcga)).mean()
    agree_sh = (np.sign(sig.theta_shrunk) == np.sign(sig.theta_tcga)).mean()
    mse_raw, mse_sh = ((d.theta - d.theta_tcga) ** 2).mean(), ((d.theta_shrunk - d.theta_tcga) ** 2).mean()
    print(f"{name:<50} n={len(d):>5} | correlation raw {r_raw:.3f} -> shrunk {r_sh:.3f} | "
          f"squared error raw {mse_raw:.3f} -> shrunk {mse_sh:.3f} | direction (TCGA-significant) {agree_raw:.1%} -> {agree_sh:.1%}")
