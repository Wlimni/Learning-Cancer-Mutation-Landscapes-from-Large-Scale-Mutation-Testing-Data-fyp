"""
Script Name: rediscover_compare.py
Purpose : Compare the Phase 2 test with Rediscover on identical input, one
          half at a time:
            B vs A  -- my tail on THEIR probabilities   (isolates the tail)
            C vs A  -- my rate fit + my tail, end to end
          plus how well each rate model reproduces the observed totals.
Inputs  : data/validation/msk468_matrix.npy, rediscover_pm.csv, rediscover_p.csv
Outputs : printed agreement table
Notes   : run rediscover_run.R first.
"""
import numpy as np
import pandas as pd

from common import OUT, load_msk468, load_phase2_functions

f = load_phase2_functions()
alt, genes = load_msk468()
tested = np.ones_like(alt)
pm_theirs = pd.read_csv(OUT / "rediscover_pm.csv").to_numpy().T        # -> samples x genes
pi_mine = f["fit_pi"](alt, tested)

print("Rate models -- max error reproducing the observed totals")
for name, pm in [("Rediscover getPM", pm_theirs), ("Phase 2 fit_pi", pi_mine)]:
    print(f"  {name:<18} genes {np.abs(pm.sum(0) - alt.sum(0)).max():.2e} | "
          f"patients {np.abs(pm.sum(1) - alt.sum(1)).max():.2e}")
print(f"  correlation between the two: {np.corrcoef(pm_theirs.ravel(), pi_mine.ravel())[0, 1]:.4f}")

ref = pd.read_csv(OUT / "rediscover_p.csv")
rows = []
for i, j, p_ref in ref.itertuples(index=False):
    i, j = int(i), int(j)
    both = int((alt[:, i] & alt[:, j]).sum())
    tail = f["poisson_binomial_tail"]
    rows.append((p_ref, tail(pm_theirs[:, i] * pm_theirs[:, j], both, upper=True),
                 tail(pi_mine[:, i] * pi_mine[:, j], both, upper=True)))
A, B, C = np.array(rows).T

for label, x in [("B vs A: my tail, their probabilities", B), ("C vs A: full pipeline vs Rediscover", C)]:
    ok = (A > 1e-300) & (x > 1e-300)
    rel = np.abs(x[ok] - A[ok]) / A[ok]
    same = ((x < 0.05) == (A < 0.05))
    print(f"\n{label}")
    print(f"  median relative difference {np.median(rel):.1e}")
    print(f"  same call at p < 0.05: {same.sum()} of {len(A)} ({same.mean():.1%})")
print(f"\nsignificant at p < 0.05 -- Rediscover {(A < 0.05).sum()}, Phase 2 {(C < 0.05).sum()}")
print(f"Spearman correlation of the ranking: {pd.Series(A).rank().corr(pd.Series(C).rank()):.3f}")
print(f"negative p-values returned by Rediscover: {(A < 0).sum()}")
