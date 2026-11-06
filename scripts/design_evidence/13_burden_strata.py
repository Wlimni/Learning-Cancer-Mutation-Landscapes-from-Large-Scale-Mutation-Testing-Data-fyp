"""
Script Name: 13_burden_strata.py
Purpose : Should Phase 2 fit gene rates within copy-number-burden groups (0 / 1-3 / >= 4
          deep CNV calls) as well as subtypes? Precedent: CDCOCA (Kumar 2011) and Zack et al.
          2013 show copy-number co-occurrence must account for each tumour's overall
          copy-number complexity. Compared on the same 10 cancer types as
          06_design_comparison.py, by (1) out-of-centre replication and (2) what happens to
          the 24 copy-number-driven candidates re-tested in 12_burden_adjusted_retest.py.
Inputs  : p2_variants.py harness; data/processed/burden_adjusted_retest.parquet
Outputs : results/analyses/burden_strata.parquet; printed comparison
"""
import sys
import pandas as pd
sys.path.insert(0, "scripts/analyses")
import p2_variants as V

CTS = ["Non-Small Cell Lung Cancer", "Breast Cancer", "Colorectal Cancer", "Glioma", "Melanoma", "Prostate Cancer",
       "Pancreatic Cancer", "Ovarian Cancer", "Endometrial Cancer", "Bladder Cancer"]
DESIGNS = {"C subtype (current)": dict(group="nonhyper", subtype="yes"),
           "F subtype + CNV burden": dict(group="nonhyper", subtype="yes", burden="yes"),
           "G subtype + CNV burden (genes outside the analysis)": dict(group="nonhyper", subtype="yes", burden="outside")}
k = ["Gene_A", "Gene_B"]
rows, full_res = [], []
for name, kw in DESIGNS.items():
    for ct in CTS:
        full = V.run(ct, **kw)
        genes = sorted(set(full.Gene_A) | set(full.Gene_B))
        m = V.run(ct, center="MSK", genes_fixed=genes, **kw)
        o = V.run(ct, center="nonMSK", genes_fixed=genes, **kw)
        r = m[m.q < .05].merge(o, on=k, suffixes=("", "_o"))
        rows.append(dict(design=name, cancer_type=ct, pairs=len(full), significant=int((full.q < .05).sum()),
                         msk_sig_testable=len(r), replicated=int(((r.p_o < .05) & (r.direction_o == r.direction)).sum())))
        full_res.append(full.assign(design=name, Cancer_Type=ct))
        print(rows[-1], flush=True)
t = pd.DataFrame(rows)
s = t.groupby("design")[["pairs", "significant", "msk_sig_testable", "replicated"]].sum()
s["replication_rate"] = (s.replicated / s.msk_sig_testable).round(3)
print(s.to_string())

fr = pd.concat(full_res)
fr.to_parquet("results/analyses/burden_strata.parquet", index=False)
chk = pd.read_parquet("data/processed/burden_adjusted_retest.parquet")
chk = chk[chk.Cancer_Type.isin(CTS)]
both = fr.pivot_table(index=["Cancer_Type", "Gene_A", "Gene_B"], columns="design", values="q").reset_index()
chk = chk.merge(both, on=["Cancer_Type", "Gene_A", "Gene_B"], how="left")
print("\nCopy-number-driven candidates: verdict of the logistic re-test vs significance with burden strata")
print(chk[["Cancer_Type", "Gene_A", "Gene_B", "verdict"] + list(DESIGNS)].round(4).to_string(index=False))
