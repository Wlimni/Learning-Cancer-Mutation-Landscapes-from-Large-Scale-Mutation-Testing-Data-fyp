"""
Script Name: 14_tcga_external_validation.py
Purpose : External validation of Phase 2 in TCGA PanCancer Atlas (independent patients,
          whole-exome sequencing, different centres and calling pipelines), via the public
          cBioPortal API. The same pairs are re-tested with the same conditional test, after
          applying Phase 1/2 rules as closely as the data allow:
            - mutations: protein-changing; a missense change counts only if the same protein
              change is kept by Phase 1's pathogenicity filter in GENIE (truncating always count)
            - copy number: deep GISTIC calls (+2/-2), direction-consistent with OncoKB roles
            - patients: TMB < 10 mut/Mb (non-hypermutated, as Phase 2's main analysis)
            - strata: TCGA study (e.g. LUAD vs LUSC, GBM vs LGG) ~ detailed subtype
Inputs  : cBioPortal REST API (cached to data/external/tcga/), data/processed/{comutation_pairs,
          alterations_long}.parquet, data/external/oncokb_curated_genes.json, Phase 2 test code
Outputs : data/processed/tcga_validation.parquet; printed concordance summary
"""
import json, sys, time, urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, "scripts/validation")
from common import load_phase2_functions, bh_fdr

F = load_phase2_functions()
API = "https://www.cbioportal.org/api"
CACHE = Path("data/external/tcga"); CACHE.mkdir(parents=True, exist_ok=True)
STUDIES = {  # TCGA PanCancer Atlas study -> GENIE cancer type
    "Non-Small Cell Lung Cancer": ["luad", "lusc"], "Breast Cancer": ["brca"], "Colorectal Cancer": ["coadread"],
    "Glioma": ["gbm", "lgg"], "Melanoma": ["skcm"], "Prostate Cancer": ["prad"], "Pancreatic Cancer": ["paad"],
    "Ovarian Cancer": ["ov"], "Endometrial Cancer": ["ucec"], "Bladder Cancer": ["blca"],
    "Esophagogastric Cancer": ["stad", "esca"], "Head and Neck Cancer": ["hnsc"],
    "Renal Cell Carcinoma": ["kirc", "kirp", "kich"], "Hepatobiliary Cancer": ["lihc", "chol"],
    "Thyroid Cancer": ["thca"], "Soft Tissue Sarcoma": ["sarc"], "Mesothelioma": ["meso"], "Germ Cell Tumor": ["tgct"],
}
TRUNC = {"Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site", "Nonstop_Mutation"}
INFRAME = {"In_Frame_Del", "In_Frame_Ins"}


def call(path, body=None, tries=4):
    for k in range(tries):
        try:
            req = urllib.request.Request(API + path, data=None if body is None else json.dumps(body).encode(),
                                         headers={"Content-Type": "application/json", "Accept": "application/json"},
                                         method="GET" if body is None else "POST")
            return json.load(urllib.request.urlopen(req, timeout=300))
        except Exception as e:
            if k == tries - 1: raise
            time.sleep(5 * (k + 1))


def cached(name, fn):
    f = CACHE / f"{name}.json"
    if f.exists(): return json.loads(f.read_text())
    d = fn(); f.write_text(json.dumps(d)); return d


cp = pd.read_parquet("data/processed/comutation_pairs.parquet")
genie_mut = pd.read_parquet("data/processed/alterations_long.parquet", columns=["Hugo_Symbol", "Alteration_Type", "HGVSp_Short"])
kept_missense = set(zip(genie_mut.Hugo_Symbol, genie_mut.HGVSp_Short.astype(str).str.replace("p.", "", regex=False)))
roles = json.loads(Path("data/external/oncokb_curated_genes.json").read_text())
role = {}
for r in roles:
    og, tsg = r.get("oncogene", False), r.get("tsg", False)
    role[r["hugoSymbol"]] = "BOTH" if og and tsg else "ONC" if og else "TSG" if tsg else None

rows = []
for ct, studies in STUDIES.items():
    pairs = cp[cp.Cancer_Type == ct]
    if pairs.empty: continue
    genes = sorted(set(pairs.Gene_A) | set(pairs.Gene_B))
    gmap = cached(f"genes_{ct}", lambda: call("/genes/fetch?geneIdType=HUGO_GENE_SYMBOL", genes))
    entrez = {g["hugoGeneSymbol"]: g["entrezGeneId"] for g in gmap}
    ids = [entrez[g] for g in genes if g in entrez]
    alt_rows, samples, strata, tmb_all = [], [], [], {}
    for st in studies:
        sid = f"{st}_tcga_pan_can_atlas_2018"
        smp = cached(f"samples_{st}", lambda: call(f"/sample-lists/{sid}_cnaseq/sample-ids"))
        tmb = cached(f"tmb_{st}", lambda: call(f"/studies/{sid}/clinical-data?clinicalDataType=SAMPLE&attributeId=TMB_NONSYNONYMOUS&projection=SUMMARY"))
        tmb_all.update({d["sampleId"]: float(d["value"]) for d in tmb})
        muts = cached(f"mut_{st}_{ct}", lambda: call(f"/molecular-profiles/{sid}_mutations/mutations/fetch?projection=DETAILED",
                                                     {"sampleListId": f"{sid}_cnaseq", "entrezGeneIds": ids}))
        cnas = cached(f"cna_{st}_{ct}", lambda: call(f"/molecular-profiles/{sid}_gistic/discrete-copy-number/fetch?discreteCopyNumberEventType=HOMDEL_AND_AMP&projection=DETAILED",
                                                     {"sampleListId": f"{sid}_cnaseq", "entrezGeneIds": ids}))
        for m in muts:
            g, mt, pc = m["gene"]["hugoGeneSymbol"], m.get("mutationType", ""), str(m.get("proteinChange", ""))
            if mt in TRUNC or mt in INFRAME or (mt == "Missense_Mutation" and (g, pc) in kept_missense):
                alt_rows.append((m["sampleId"], g))
        for c in cnas:
            g, v = c["gene"]["hugoGeneSymbol"], c["alteration"]
            ok = {"TSG": v == -2, "ONC": v == 2}.get(role.get(g), True)
            if ok: alt_rows.append((c["sampleId"], g))
        samples += smp; strata += [st] * len(smp)
    s = pd.Series(strata, index=samples)
    s = s[~s.index.duplicated()]
    s = s[[tmb_all.get(x, np.nan) < 10 for x in s.index]]           # non-hypermutated, as Phase 2
    alt = pd.DataFrame(alt_rows, columns=["s", "g"]).drop_duplicates()
    X = (pd.crosstab(alt.s, alt.g).reindex(index=s.index, columns=genes, fill_value=0) > 0).to_numpy()
    T = np.ones_like(X, dtype=bool)                                   # exome + genome-wide CNA: all tested
    M = F["ConditionalModel"](X, T, s.to_numpy())
    gi = {g: i for i, g in enumerate(genes)}
    res = []
    for a, b in sorted(zip(pairs.Gene_A.map(gi), pairs.Gene_B.map(gi))):
        q = M.joint(a, b)
        if q.sum() < 5: continue                                      # same power gate
        both = int((X[:, a] & X[:, b]).sum())
        p, d, e = F["pair_test"](q, both)
        res.append((genes[a], genes[b], both, e, d, p))
    r = pd.DataFrame(res, columns=["Gene_A", "Gene_B", "tcga_n_both", "tcga_expected", "tcga_direction", "tcga_p"])
    if len(r): r["tcga_q"] = bh_fdr(r.tcga_p.to_numpy())
    r["Cancer_Type"] = ct; r["tcga_patients"] = len(s)
    rows.append(r)
    print(f"{ct:<28} TCGA patients {len(s):>5}  pairs testable {len(r):>4} of {len(pairs)}", flush=True)

tv = pd.concat(rows, ignore_index=True)
m = cp.merge(tv, on=["Cancer_Type", "Gene_A", "Gene_B"], how="left")
m.to_parquet("data/processed/tcga_validation.parquet", index=False)
t = m[m.tcga_p.notna()]
for name, sel in [("replicated candidates (GENIE)", t.interaction_candidate & t.replicated),
                  ("all candidates", t.interaction_candidate),
                  ("not significant in GENIE", t.q_value >= 0.05)]:
    g = t[sel]
    same = g.tcga_direction == g.direction
    print(f"{name:<32} testable in TCGA {len(g):>4} | same direction {same.mean():.0%} | "
          f"TCGA p<0.05 same dir {((g.tcga_p < .05) & same).mean():.0%} | TCGA q<0.05 same dir {((g.tcga_q < .05) & same).mean():.0%} | "
          f"TCGA significant opposite {((g.tcga_p < .05) & ~same).sum()}")
