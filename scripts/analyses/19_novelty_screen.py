"""
Script Name: 19_novelty_screen.py
Purpose : Which robust Phase 2 interactions has the literature NOT reported? For every strong
          candidate (replicated outside MSK, burden-robust, not weak), count Europe PMC papers
          whose title or abstract names both genes (any context; and with co-occurrence /
          exclusivity wording), with and without the cancer type. Full-text search was tried first
          and rejected: pan-cancer papers listing hundreds of genes made every pair look "known".
          A second source is the SelectSim atlas (Iyer ... Ciriello, Nat Genet 2026; 892 validated
          mutation-only co-mutations, supplementary table 7): atlas papers list hundreds of pairs only
          in supplements, which a title/abstract search cannot see.
Inputs  : data/processed/tcga_validation.parquet (Phase 2 pairs + TCGA re-test),
          data/external/selectsim/supp_tables.xlsx (SelectSim supplementary tables),
          data/processed/conditional_selection_dnds.parquet, Europe PMC REST API (cached)
Outputs : data/processed/novelty_screen.parquet; printed shortlist
Notes   : a literature count is a triage signal, not proof of novelty -- short gene symbols that
          are also English words (MET, ATM, KIT, REL) inflate counts, so they err towards "known".
          The shortlist must be checked by reading. Run from the repo root.
"""
import json, time, urllib.parse, urllib.request
from pathlib import Path

import pandas as pd

CACHE = Path("data/external/europepmc_counts.json")
cache = json.loads(CACHE.read_text()) if CACHE.exists() else {}
API = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?format=json&resultType=idlist&pageSize=1&query="
WORDS = '("co-occur*" OR "cooccur*" OR "mutually exclusive" OR "mutual exclusivity" OR "co-mutat*" OR "concurrent")'
CT_TERM = {  # cancer type -> search term
    "Non-Small Cell Lung Cancer": '(lung OR NSCLC)', "Breast Cancer": "breast", "Colorectal Cancer": '(colorectal OR colon)',
    "Glioma": '(glioma OR glioblastoma OR astrocytoma)', "Bladder Cancer": '(bladder OR urothelial)',
    "Endometrial Cancer": '(endometrial OR uterine)', "Esophagogastric Cancer": '(gastric OR esophageal OR oesophageal)',
    "Melanoma": "melanoma", "Prostate Cancer": "prostate", "Hepatobiliary Cancer": '(cholangiocarcinoma OR hepatocellular OR biliary)',
    "Pancreatic Cancer": "pancreatic", "Ovarian Cancer": "ovarian", "Soft Tissue Sarcoma": "sarcoma",
    "Head and Neck Cancer": '("head and neck" OR HNSCC)', "Renal Cell Carcinoma": '(renal OR kidney)', "Thyroid Cancer": "thyroid",
}


def ta(term):  # title-or-abstract field search
    return f'(TITLE:{term} OR ABSTRACT:{term})'


def hits(q):
    if q not in cache:
        for attempt in range(4):
            try:
                with urllib.request.urlopen(API + urllib.parse.quote(q), timeout=30) as r:
                    cache[q] = json.load(r)["hitCount"]; break
            except Exception:
                time.sleep(2 ** attempt)
        time.sleep(0.1)
        if len(cache) % 50 == 0:
            CACHE.write_text(json.dumps(cache))
    return cache.get(q)


t = pd.read_parquet("data/processed/tcga_validation.parquet")
c = t[t.interaction_candidate & t.replicated & t.burden_robust & ~t.weak_effect].copy()
g = lambda s: s.str.replace("_fusion", "", regex=False)
c["gA"], c["gB"] = g(c.Gene_A), g(c.Gene_B)
ct_term = lambda ct: CT_TERM.get(ct, '"' + ct.split()[0] + '"')
c["papers_pair"] = [hits(f'{ta(chr(34) + a + chr(34))} AND {ta(chr(34) + b + chr(34))}') for a, b in zip(c.gA, c.gB)]
c["papers_pair_in_type"] = [hits(f'{ta(chr(34) + a + chr(34))} AND {ta(chr(34) + b + chr(34))} AND {ta(ct_term(ct))}')
                            for a, b, ct in zip(c.gA, c.gB, c.Cancer_Type)]
c["papers_interaction_in_type"] = [hits(f'{ta(chr(34) + a + chr(34))} AND {ta(chr(34) + b + chr(34))} AND {ta(ct_term(ct))} AND {ta(WORDS)}')
                                   for a, b, ct in zip(c.gA, c.gB, c.Cancer_Type)]
CACHE.write_text(json.dumps(cache))

# Independent support: TCGA (same direction, p < 0.05) and the dN/dS route (same direction)
c["tcga_support"] = (c.tcga_p < 0.05) & (c.tcga_direction == c.direction)
dn = pd.read_parquet("data/processed/conditional_selection_dnds.parquet")
dn_ok = {(r.Cancer_Type, frozenset((r.A, r.B))) for r in dn.itertuples() if r.dnds_direction == r.phase2_direction and r.p < 0.05}
c["dnds_support"] = [(ct, frozenset((a, b))) in dn_ok for ct, a, b in zip(c.Cancer_Type, c.Gene_A, c.Gene_B)]
c["mechanism"] = c[["n_MUT_MUT", "n_CNV_CNV", "n_mixed"]].idxmax(axis=1).str.replace("n_", "")
# SelectSim validated pairs: "A - B_CO" / "A - B_ME" per tumour class (or pan_can)
ss = pd.read_excel("data/external/selectsim/supp_tables.xlsx", sheet_name="supplementary table 7 (part 2)")
ss[["pair", "kind"]] = ss.pair_type.str.rsplit("_", n=1, expand=True)
ss["key"] = [frozenset(p.split(" - ")) for p in ss.pair]
ss["dir"] = ss.kind.map({"CO": "co-occurring", "ME": "exclusive"})
in_type = {(k, cl): d for k, cl, d in zip(ss.key, ss.Tumor_class, ss.dir)}
pan = {k: d for k, cl, d in zip(ss.key, ss.Tumor_class, ss.dir) if cl == "pan_can"}
cls = c.Cancer_Type.str.lower().str.replace(" ", "_").str.replace("hepatobiliary_cancer", "biliary_tract_cancer")
keys = [frozenset((a, b)) for a, b in zip(c.gA, c.gB)]
c["selectsim"] = [("same type, same direction" if in_type.get((k, t)) == d else
                   "same type, opposite direction" if (k, t) in in_type else
                   "pan-cancer only" if k in pan else "not reported")
                  for k, t, d in zip(keys, cls, c.direction)]
c["novelty"] = pd.cut(c.papers_pair_in_type, [-1, 2, 19, 1e9], labels=["rarely or never reported", "some reports", "well known"])

c["new_candidate"] = (c.novelty == "rarely or never reported") & (c.selectsim == "not reported")

cols = ["Cancer_Type", "Gene_A", "Gene_B", "direction", "n_both_altered", "expected_both", "fold_enrichment", "q_value",
        "mechanism", "tcga_support", "dnds_support", "papers_pair", "papers_pair_in_type", "papers_interaction_in_type", "novelty",
        "selectsim", "new_candidate"]
c[cols].to_parquet("data/processed/novelty_screen.parquet", index=False)
print(c.groupby(["direction", "novelty"], observed=True).size().unstack(fill_value=0))
print("\nSelectSim atlas status, by what drives the pair (SelectSim analysed mutations only):")
print(pd.crosstab(c.mechanism, c.selectsim))
print(f"\nNot in titles/abstracts AND not in SelectSim: {c.new_candidate.sum()} of {len(c)}")
print("\nShortlist: rarely/never reported, supported by TCGA or dN/dS")
sl = c[c.new_candidate & (c.tcga_support | c.dnds_support)].sort_values("q_value")
print(sl[cols].round(4).to_string(index=False))
print("\nRarely/never reported, strongest without external support (TCGA untestable or not significant)")
print(c[c.new_candidate & ~(c.tcga_support | c.dnds_support)]
      .sort_values("q_value").head(20)[cols].round(4).to_string(index=False))
