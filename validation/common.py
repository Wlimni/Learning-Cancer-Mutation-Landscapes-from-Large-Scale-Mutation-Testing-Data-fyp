"""
Script Name: common.py
Purpose : Shared helpers for the validation scripts -- loads the Phase 2 test
          functions straight out of the notebook, so what is validated is
          exactly what the notebook runs (no second copy that can drift).
Inputs  : notebooks/02_phase2_comutation_matrix.ipynb (cell id "disc0v3rcd")
Outputs : none (imported by the other scripts)
"""
import json
from pathlib import Path

import numpy as np
from scipy.special import expit, logit
from scipy.stats import norm

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "data" / "validation"
OUT.mkdir(parents=True, exist_ok=True)
PHASE2_NB = ROOT / "notebooks" / "02_phase2_comutation_matrix.ipynb"
TEST_CELL_ID = "disc0v3rcd"


def load_phase2_functions() -> dict:
    """Execute the notebook's rate-model + Poisson-Binomial cell; return its namespace."""
    nb = json.loads(PHASE2_NB.read_text())
    cell = next(c for c in nb["cells"] if c.get("id") == TEST_CELL_ID)
    ns = {"np": np, "norm": norm, "expit": expit, "logit": logit}
    exec("".join(cell["source"]), ns)
    return ns


def load_msk468():
    """The validation matrix built by build_msk468_matrix.py: (samples x genes bool, gene list)."""
    return np.load(OUT / "msk468_matrix.npy"), json.loads((OUT / "msk468_genes.json").read_text())


def bh_fdr(p: np.ndarray) -> np.ndarray:
    n = len(p)
    order = np.argsort(p)
    ranks = np.empty(n, int)
    ranks[order] = np.arange(1, n + 1)
    q = np.minimum.accumulate((p * n / ranks)[order][::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(q, 0, 1)
    return out
