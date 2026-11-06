"""
Script Name: 07_joint_coverage_filter.py
Purpose : Is Filter 3 (>= 50% of patients tested for BOTH genes) needed on top of
          Filter 2 (gene altered in >= 3% of ALL patients)? Compare the pairs it
          removes against the pairs it keeps, by out-of-centre replication.
Inputs  : p2_variants.py harness (non-hypermutated patients, subtype strata)
Outputs : printed table
"""
import sys, pandas as pd
sys.path.insert(0, "scripts/analyses")
import p2_variants as V

k = ["Gene_A", "Gene_B"]
tot = {"kept": [0, 0], "removed_by_F3": [0, 0]}
for ct in ["Breast Cancer", "Colorectal Cancer", "Non-Small Cell Lung Cancer", "Melanoma", "Bladder Cancer"]:
    kw = dict(group="nonhyper", subtype="yes")
    with_f3, no_f3 = V.run(ct, **kw), V.run(ct, joint_cov=0, **kw)
    removed = no_f3.merge(with_f3[k], on=k, how="left", indicator=True)
    removed = removed[removed._merge == "left_only"][k]
    genes = sorted(set(no_f3.Gene_A) | set(no_f3.Gene_B))
    m = V.run(ct, center="MSK", genes_fixed=genes, joint_cov=0, **kw)
    o = V.run(ct, center="nonMSK", genes_fixed=genes, joint_cov=0, **kw)
    r = m[m.q < .05].merge(o, on=k, suffixes=("", "_o"))
    r["rep"] = (r.p_o < .05) & (r.direction_o == r.direction)
    r = r.merge(removed.assign(removed=True), on=k, how="left").fillna({"removed": False})
    for grp, g in r.groupby("removed"):
        key = "removed_by_F3" if grp else "kept"
        tot[key][0] += len(g); tot[key][1] += int(g.rep.sum())
    print(f"{ct:<28} pairs with F3 {len(with_f3):>4}, without {len(no_f3):>4} (+{len(removed)}) | "
          f"MSK hits among removed pairs: {int(r.removed.sum())}, replicated {int(r[r.removed].rep.sum())}", flush=True)
for key, (n, rep) in tot.items():
    print(f"{key:<14} MSK hits testable elsewhere {n:>4}, replicated {rep:>4} ({rep / max(n, 1):.0%})")
