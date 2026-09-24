"""
Script Name: null_calibration.py
Purpose : Are the p-values honest? Simulate data with NO interactions (every
          cell independent given the logistic patient x gene model fitted to
          the benchmark), run each test, and count false alarms. A correct
          test gives ~5% of p < 0.05 and ~0 pairs significant after FDR.
            plug-in      -- per-patient rate estimated from the same data
                            (DISCOVER-style; what Phase 2 used before)
            conditional  -- ConditionalModel, the test Phase 2 uses now
            true pi      -- the generating probabilities (ideal reference)
          Plus a positive control: 5 planted interactions must be recovered.
Inputs  : data/validation/msk468_matrix.npy (build_msk468_matrix.py)
Outputs : printed calibration table
Notes   : ~20 min; seeded. The Rediscover arm is rediscover_run.R on a matrix
          written here (data/validation/null_genes_by_samples.csv).
"""
import itertools

import numpy as np
import pandas as pd

from common import OUT, bh_fdr, load_msk468, load_phase2_functions

f = load_phase2_functions()
pair_test = f["pair_test"]
rng = np.random.default_rng(31337)
alt_real, _ = load_msk468()
tested = np.ones_like(alt_real)
PAIRS = list(itertools.combinations(range(alt_real.shape[1]), 2))
N_REPS = 4

pi_true = f["fit_logistic_pi"](alt_real, tested)          # the null the simulations come from


def pvals(alt, arm):
    if arm == "conditional":
        m = f["ConditionalModel"](alt, tested)
        q_of = lambda i, j: m.joint(i, j)
    else:
        pi = pi_true if arm == "true pi" else f["fit_pi"](alt, tested)
        q_of = lambda i, j: pi[:, i] * pi[:, j]
    return np.array([pair_test(q_of(i, j), int((alt[:, i] & alt[:, j]).sum()))[0] for i, j in PAIRS])


arms = ["plug-in", "conditional", "true pi"]
acc = {a: [] for a in arms}
for rep in range(N_REPS):
    sim = rng.random(pi_true.shape) < pi_true
    if rep == 0:
        pd.DataFrame(sim.T.astype(int)).to_csv(OUT / "null_genes_by_samples.csv", index=False, header=False)
    for a in arms:
        p = pvals(sim, a)
        acc[a].append(((p < 0.05).mean(), (bh_fdr(p) < 0.05).mean()))
    print(f"  replicate {rep + 1}/{N_REPS} done", flush=True)

print(f"\n{'NULL (no interactions exist)':<30}{'frac p<0.05':>13}{'pairs sig. after FDR':>22}")
for a, v in acc.items():
    v = np.array(v)
    print(f"{a:<30}{v[:, 0].mean():>13.3f}{v[:, 1].mean():>21.2%}")
print(f"{'target':<30}{0.05:>13.3f}{0.0:>21.2%}")

sim = rng.random(pi_true.shape) < pi_true
planted = {(0, 1), (2, 3), (4, 5), (6, 7), (8, 9)}
for i, j in planted:
    idx = np.flatnonzero(sim[:, i])
    sim[rng.choice(idx, size=len(idx) // 2, replace=False), j] = True
q = bh_fdr(pvals(sim, "conditional"))
hits = sum(1 for k, pr in enumerate(PAIRS) if pr in planted and q[k] < 0.05)
print(f"\nPOSITIVE CONTROL (conditional): {hits}/5 planted interactions recovered; "
      f"{(q < 0.05).sum()} of {len(q)} pairs significant overall")
