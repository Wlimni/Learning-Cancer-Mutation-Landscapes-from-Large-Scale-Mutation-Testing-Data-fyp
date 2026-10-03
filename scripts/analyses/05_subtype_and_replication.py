"""
Script Name: 05_subtype_and_replication.py
Purpose : Test literature-suggested changes on Phase 2 (non-hypermutated patients,
          conditional test): (1) subtype confounding -- gene rates fitted within
          detailed subtypes (van de Haar 2019; Park & Lehner 2015); (2) replication
          -- discover in MSK patients, re-test the same pairs in other centres.
Inputs  : p2_variants.py (notebook functions), data/processed
Outputs : results/analyses/lit_tests.parquet (pair-level results for every variant)
"""
import sys, pandas as pd
sys.path.insert(0, "scripts/analyses")
import p2_variants as V          # executes the notebook cells once (chdir into notebooks/)

CTS = ["Non-Small Cell Lung Cancer", "Breast Cancer", "Colorectal Cancer", "Glioma", "Melanoma", "Prostate Cancer",
       "Pancreatic Cancer", "Ovarian Cancer", "Endometrial Cancer", "Bladder Cancer"]
out = []
for ct in CTS:
    base = V.run(ct, group="nonhyper")
    sub = V.run(ct, group="nonhyper", subtype="yes")
    genes = sorted(set(base.Gene_A) | set(base.Gene_B))
    msk = V.run(ct, group="nonhyper", center="MSK", genes_fixed=genes)
    rep = V.run(ct, group="nonhyper", center="nonMSK", genes_fixed=genes)
    for tag, d in [("base", base), ("subtype", sub), ("MSK", msk), ("nonMSK", rep)]:
        if d is not None and len(d):
            out.append(d.assign(Cancer_Type=ct, variant=tag, n_patients=d.attrs["n_patients"]))
    s = lambda d: 0 if d is None or not len(d) else int((d.q < .05).sum())
    print(f"{ct:<28} base {s(base):>4}/{len(base):<5} subtype {s(sub):>4}/{len(sub):<5} "
          f"MSK {s(msk):>4}/{0 if msk is None else len(msk):<5} (n={0 if msk is None else msk.attrs['n_patients']}) "
          f"other {s(rep):>4}/{0 if rep is None else len(rep):<5} (n={0 if rep is None else rep.attrs['n_patients']})", flush=True)
pd.concat(out).to_parquet("results/analyses/subtype_and_replication.parquet")
