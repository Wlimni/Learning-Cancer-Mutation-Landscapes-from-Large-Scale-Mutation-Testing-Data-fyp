"""
Script Name: 06_design_comparison.py
Purpose : Score whole Phase 2 designs by out-of-centre replication: discover in MSK
          patients, re-test the same pairs in all other centres (independent
          patients, different panels and pipelines). A design is better if more of
          its discoveries replicate (p < 0.05, same direction).
Designs : A current   -- all patients, gene rates per TMB group, gene filter vs all patients
          B split     -- non-hypermutated only (Jason), same gene filter
          C split+sub -- B + gene rates per detailed subtype
          D final     -- C + gene filter counted among patients TESTED for the gene
          E mut-only  -- D on mutations only, all mutation-tested patients
Inputs  : p2_variants.py harness
Outputs : printed table; results/analyses/design_compare.parquet
"""
import sys, pandas as pd
sys.path.insert(0, "scripts/analyses")
import p2_variants as V

CTS = ["Non-Small Cell Lung Cancer", "Breast Cancer", "Colorectal Cancer", "Glioma", "Melanoma", "Prostate Cancer",
       "Pancreatic Cancer", "Ovarian Cancer", "Endometrial Cancer", "Bladder Cancer"]
DESIGNS = {"A current": dict(group="all"), "B split": dict(group="nonhyper"),
           "C split+sub": dict(group="nonhyper", subtype="yes"),
           "D final": dict(group="nonhyper", subtype="yes", gene_rule="of_tested"),
           "E mut-only": dict(mode="snv", group="nonhyper", subtype="yes", gene_rule="of_tested")}
rows = []
for name, kw in DESIGNS.items():
    for ct in CTS:
        full = V.run(ct, **kw)
        if full is None or not len(full):
            continue
        genes = sorted(set(full.Gene_A) | set(full.Gene_B))
        m = V.run(ct, center="MSK", genes_fixed=genes, **kw)
        o = V.run(ct, center="nonMSK", genes_fixed=genes, **kw)
        k = ["Gene_A", "Gene_B"]
        r = m[m.q < .05].merge(o, on=k, suffixes=("", "_o")) if (m is not None and o is not None and len(m) and len(o)) else pd.DataFrame()
        rows.append(dict(design=name, cancer_type=ct, patients=full.attrs["n_patients"], pairs=len(full),
                         significant=int((full.q < .05).sum()),
                         msk_sig_testable=len(r),
                         replicated=int(((r.p_o < .05) & (r.direction_o == r.direction)).sum()) if len(r) else 0,
                         opposite=int(((r.p_o < .05) & (r.direction_o != r.direction)).sum()) if len(r) else 0))
        print(rows[-1], flush=True)
t = pd.DataFrame(rows)
t.to_parquet("results/analyses/design_comparison.parquet")
s = t.groupby("design")[["pairs", "significant", "msk_sig_testable", "replicated", "opposite"]].sum()
s["replication_rate"] = (s.replicated / s.msk_sig_testable).round(3)
print(s.to_string())
