"""
Script Name: tail_ground_truth.py
Purpose : Check the exact Poisson-Binomial tail against two sources of truth
          that are independent of it:
            A. identical probabilities -> it IS a Binomial; scipy's tail
               (regularised incomplete beta) is accurate to ~1e-308
            B. heterogeneous probabilities -> 60-digit mpmath DP, no truncation
Inputs  : the notebook's poisson_binomial_tail (via common.py)
Outputs : printed relative errors
Notes   : needs `pip install mpmath`; test B takes a few minutes.
"""
import numpy as np
from mpmath import mp, mpf
from scipy.stats import binom

from common import load_phase2_functions

tail = load_phase2_functions()["poisson_binomial_tail"]

print("A. identical probabilities vs exact Binomial (scipy)")
n, p = 36841, 100 / 36841
probs = np.full(n, p)
for k in [100, 120, 150, 200, 250, 300, 400, 500]:
    truth, mine = binom.sf(k - 1, n, p), tail(probs, k, upper=True, force="exact")
    print(f"  k={k:4d}  truth {truth:.6e}  mine {mine:.6e}  rel err {abs(mine - truth) / truth:.1e}")

print("\nB. heterogeneous probabilities vs 60-digit mpmath DP")
mp.dps = 60
probs2 = np.clip(np.random.default_rng(7).beta(0.3, 40, 2000), 1e-9, 0.99)
pmf = [mpf(0)] * (len(probs2) + 1)
pmf[0] = mpf(1)
for pr in probs2:                                   # full-support exact convolution
    q = mpf(float(pr))
    for m in range(len(probs2), 0, -1):
        pmf[m] = pmf[m] * (1 - q) + pmf[m - 1] * q
    pmf[0] *= (1 - q)
for k in [40, 60, 80, 100, 120, 150, 200]:
    truth, mine = sum(pmf[k:]), tail(probs2, k, upper=True, force="exact")
    print(f"  k={k:4d}  truth {float(truth):.6e}  mine {mine:.6e}  rel err {float(abs(mpf(mine) - truth) / truth):.1e}")
