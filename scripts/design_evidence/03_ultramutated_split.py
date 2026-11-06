"""
Script Name: 03_ultramutated_split.py
Purpose : Is the hypermutated group one population? Re-test its pairs with the
          ultramutated (TMB >= 100, typically POLE) patients removed, and within
          each sub-group, to see how much of the signal is MSI-vs-POLE mixing.
Inputs  : data/processed/*.parquet, notebook test functions (validation/common.py)
Outputs : printed table
"""
import sys, numpy as np, pandas as pd
from scipy.stats import false_discovery_control
sys.path.insert(0, "scripts/validation"); from common import load_phase2_functions
F = load_phase2_functions()
P = "data/processed/"
cl = pd.read_parquet(P + "clinical_tidy.parquet"); cl = cl[cl.IS_REPRESENTATIVE & cl.CNA_TESTED]
al = pd.read_parquet(P + "alterations_long.parquet", columns=["Sample_ID", "Hugo_Symbol"])
cov = pd.read_parquet(P + "panel_gene_coverage.parquet")[["SEQ_ASSAY_ID", "Hugo_Symbol"]].merge(
      pd.read_parquet(P + "cna_gene_panel_coverage.parquet"))
cp = pd.read_parquet(P + "comutation_pairs.parquet")

for ct in ["Endometrial Cancer", "Colorectal Cancer"]:
    pairs = cp[(cp.Cancer_Type == ct) & cp.hyper_p_value.notna()]
    genes = sorted(set(pairs.Gene_A) | set(pairs.Gene_B)); gi = {g: i for i, g in enumerate(genes)}
    c = cl[(cl.CANCER_TYPE == ct) & (cl.HYPERMUTATION_STATUS == "hypermutated")].set_index("SAMPLE_ID")
    X = (pd.crosstab(al.Sample_ID, al.Hugo_Symbol).reindex(index=c.index, columns=genes, fill_value=0) > 0).values
    T = (cov.assign(v=True).pivot_table(index="SEQ_ASSAY_ID", columns="Hugo_Symbol", values="v", aggfunc="any")
         .reindex(index=c.SEQ_ASSAY_ID, columns=genes).fillna(False).astype(bool).values)
    ultra = (c.TMB_MUT_PER_MB >= 100).values
    print(f"\n{ct}: hypermutated {len(c):,}, of which TMB >= 100: {ultra.sum():,} ({ultra.mean():.0%})")
    for name, m in [("all hypermutated (Phase 2)", np.ones(len(c), bool)), ("TMB 10-100 only", ~ultra), ("TMB >= 100 only", ultra)]:
        alt, tst = X[m], T[m]
        M = F["ConditionalModel"](alt, tst); ps = []
        for a, b in sorted(zip(pairs.Gene_A.map(gi), pairs.Gene_B.map(gi))):
            jt = tst[:, a] & tst[:, b]
            ps.append(F["pair_test"](M.joint(a, b)[jt], int((alt[:, a] & alt[:, b] & jt).sum()))[0])
        print(f"  {name:<28} n={m.sum():>5,}  significant {(false_discovery_control(ps) < .05).mean():.1%} of {len(ps):,} pairs")
