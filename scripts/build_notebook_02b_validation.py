"""
Script Name: build_notebook_02b_validation.py
Purpose : Generate notebooks/02b_phase2_validation.ipynb -- which alteration set (SNV, SNV+CNV, SNV+CNV (drivers), SNV+CNV+SV) gives
          the better gene-pair list? Three checks: known biology, precision and recall where the truth is known,
          and independent patients (TCGA). The comparison of methods (Fisher, DISCOVER-style, SELECT) is in 02c.
Inputs  : none (the notebook reads data/processed/*.parquet and data/curated/textbook_pairs.csv)
Outputs : notebooks/02b_phase2_validation.ipynb (unexecuted; run with nbconvert --execute)
Note    : the notebook only reads saved outputs; re-run scripts/analyses/23_alteration_set_validation.py to
          regenerate them.
"""
import nbformat as nbf

md, code = nbf.v4.new_markdown_cell, nbf.v4.new_code_cell
cells = [
md("""# Phase 2b · Which alteration set gives the better gene-pair list?

Notebook 02 builds the gene-pair list four times, with the same test, on the same patients -- only the definition
of "altered" differs:

| Set | A gene counts as altered if it has ... |
|---|---|
| **SNV** | a likely-harmful mutation (incl. TERT promoter and MET exon 14) |
| **SNV+CNV** | a mutation **or** a deep copy-number change -- **the main table** (for now) |
| **SNV+CNV (drivers)** | a mutation **or** a deep copy-number change **in an established driver gene** (COSMIC CGC / Vogelstein / CIViC); other copy-number calls count as not altered |
| **SNV+CNV+SV** | a mutation, a deep copy-number change **or** a driver fusion |

Each set is judged three ways:

| Part | Question | Data |
|---|---|---|
| 1 | **Known biology:** how many published gene pairs does it find? | 58 pairs from the literature (`data/curated/textbook_pairs.csv`) |
| 2 | **Precision and recall:** are its calls real, and how many real links does it find? | real data with links of known strength planted |
| 3 | **Independent patients:** do its candidates hold up in TCGA? | TCGA PanCancer Atlas |

How our test compares with other methods (Fisher, DISCOVER-style, SELECT) is in notebook 02c, Section G."""),
code("""import numpy as np
import pandas as pd

P = "../data/processed/"
C = "../data/curated/"
pd.set_option("display.width", 160)
SETS = {"SNV": "comutation_pairs_snv.parquet", "SNV+CNV": "comutation_pairs.parquet",
        "SNV+CNV (drivers)": "comutation_pairs_snv_cnv_drivers.parquet",
        "SNV+CNV+SV": "comutation_pairs_snv_cnv_sv.parquet"}
pairs = {n: pd.read_parquet(P + f) for n, f in SETS.items()}
v = pd.read_parquet(P + "alteration_set_validation.parquet")      # scripts/analyses/23_alteration_set_validation.py
pd.DataFrame({n: {"pairs tested": len(p), "interaction candidates": int(p.interaction_candidate.sum())}
              for n, p in pairs.items()})"""),

md("""---
## 1 · Known biology: how many published pairs does each set find?

58 gene pairs with a published direction (15 co-occurring, 43 exclusive), each with a reference and DOI, chosen
**before** looking at the results. They cover mutations, copy-number changes and fusions in 17 cancer types -- e.g.
*EGFR* vs *KRAS* and *ALK* fusion vs *EGFR* in lung cancer, *MDM2* amplification vs *TP53* in glioma, *SF3B1* vs
*SRSF2* in myelodysplasia. A pair is **found** if it is an interaction candidate in the published direction;
otherwise the reason is given."""),
code("""tb = pd.read_csv(C + "textbook_pairs.csv")


def verdict(p, r):
    m = p[(p.Cancer_Type == r.Cancer_Type) & (((p.Gene_A == r.Gene_A) & (p.Gene_B == r.Gene_B))
                                             | ((p.Gene_A == r.Gene_B) & (p.Gene_B == r.Gene_A)))]
    if m.empty: return "not tested"
    m = m.iloc[0]
    if m.q_value >= .05: return "not significant"
    if m.direction != r.expected_direction: return "opposite direction"
    if not m.interaction_candidate: return "removed: one DNA event"
    return "found"


known = tb[["Cancer_Type", "Gene_A", "Gene_B", "expected_direction", "reference"]].copy()
for n, p in pairs.items():
    known[n] = [verdict(p, r) for r in tb.itertuples()]
pd.DataFrame({n: known[n].value_counts() for n in pairs}).fillna(0).astype(int)"""),
md("""Only the pairs that some set does not find:"""),
code("""known[known[list(pairs)].ne("found").any(axis=1)].drop(columns="reference")"""),
md("""> **Reading.** Adding alteration types finds more known biology: **35 of 58** pairs with mutations only, **44** with
> copy number, **43** with copy number restricted to established drivers, **50** with fusions. Most of the difference
> is pairs a set cannot see at all ("not tested"): copy-number pairs (*MDM2* vs *TP53*, *CCNE1* vs *BRCA1/2*) need
> copy number, fusion pairs (*ALK* vs *EGFR*, *ERG* vs *SPOP*, *BRAF* vs *RET*) need fusions. No set finds a pair in
> the opposite direction.
>
> Restricting copy number to established drivers costs **one** known pair, *ERBB2* vs *FGFR2* in oesophagogastric
> cancer (the pair is still tested, but no longer significant). Two other misses are worth knowing: *KIT* vs *PDGFRA*
> in GIST is found by SNV but removed once copy number is added -- the two genes are neighbours on 4q12, so the
> one-DNA-event rule treats them as one amplicon; and a few pairs are tested but not significant in every set (*BRAF*
> vs *NRAS* in colorectal cancer, *SF3B1* vs *SRSF2*, *RHOA* + *TET2*) -- the recall cost of a strict test (Part 2)."""),

md("""---
## 2 · Precision and recall where the truth is known

On real data nobody knows every true link, so precision and recall are measured on data where the answer is known:
each set's own real data (9 large cancer types), **shuffled within subtypes** so that every patient keeps its number of
alterations and every gene its frequency but no gene-gene link survives, then **~215 links of known strength
planted** (co-occurring or exclusive, 1.5x to 3x; 6 replicates). A planted link is **found** if called (q < 0.05) in
the right direction; every other call is wrong. The same shuffled data without planted links shows how often a set
calls something that is not there."""),
code("""rows = []
for n in SETS:
    nl = v[(v.alteration_set == n) & (v.check == "no links")]
    s = v[(v.alteration_set == n) & (v.check == "planted links")]
    call = s.q < .05; hit = call & (s.truth != "none") & (s.direction == s.truth)
    rows.append({"alteration set": n, "no links: called (q < 0.05)": (nl.q < .05).mean(),
                 "planted links": int((s.truth != "none").sum()), "calls": int(call.sum()),
                 "precision": hit.sum() / max(call.sum(), 1), "recall": hit.sum() / (s.truth != "none").sum()})
pd.DataFrame(rows).set_index("alteration set").round(3)"""),
md("""Recall by the strength of the planted link:"""),
code("""s = v[(v.check == "planted links") & (v.truth != "none")]
s.assign(found=(s.q < .05) & (s.direction == s.truth)).groupby(["truth", "fold", "alteration_set"]).found.mean() \\
 .unstack().round(2)"""),
md("""> **Reading.** The test is **precise on every set**: nothing is called on data without links, and 95-99% of calls on
> planted data are real. Recall is 46-52% -- strong links are found, weak ones (1.5x) mostly not. **SNV+CNV (drivers)
> has the highest precision (99.1% vs 96.2%) and a higher recall than SNV+CNV (49.1% vs 46.3%)**, on a similar number
> of tested pairs (2,339 vs 2,552), so the gain is not just from testing fewer pairs. Copy-number calls in genes that
> are not established drivers add noise, not signal."""),

md("""---
## 3 · Independent patients: TCGA

Each set's candidates re-tested with the same test in TCGA (whole exome, other centres, non-hypermutated patients;
mutations only for the SNV set). TCGA has about 1/20 of GENIE's patients, so **direction** is the main measure.
**Control** = pairs *not* significant in GENIE: how much agreement to expect without a real link. TCGA has no fusion
calls, so fusion pairs cannot be tested there."""),
code("""rows = []
for n in SETS:
    t = v[(v.alteration_set == n) & (v.check == "TCGA")]
    for grp, sel in [("candidates", t.interaction_candidate == True), ("control: not significant", t.q_value >= .05)]:
        d = t[sel]; same = d.tcga_direction == d.direction
        rows.append({"alteration set": n, "group": grp, "testable in TCGA": len(d), "same direction": round(same.mean(), 3),
                     "significant, same direction": round(((d.tcga_p < .05) & same).mean(), 3),
                     "significant, opposite": int(((d.tcga_p < .05) & ~same).sum())})
pd.DataFrame(rows).set_index(["alteration set", "group"])"""),
md("""> **Reading.** Candidates of every set keep their direction in TCGA far more often than control pairs, and none is
> significant the other way. SNV candidates hold up best (93% vs 88% once copy number is added) -- mutations are
> called more consistently than copy number between cohorts -- but the SNV set has a fifth as many candidates. The
> driver-list set keeps the same agreement as SNV+CNV (87% vs 88%) while its control pairs agree less often (64% vs
> 69%), so its margin over chance is the larger one (24 vs 19 points)."""),

md("""---
## 4 · What does the copy-number driver list change?

In an amplified or deleted region usually only one gene is the real driver; its neighbours are passengers or
biomarkers (supervisor). **SNV+CNV (drivers)** keeps copy-number calls only in established driver genes (Phase 1,
Section 3b). Comparing it with SNV+CNV: which candidates disappear, which remain, and do the ones that disappear look
like real biology?"""),
code("""key = ["Cancer_Type", "Gene_A", "Gene_B"]
a, b = pairs["SNV+CNV"], pairs["SNV+CNV (drivers)"]
ca, cb = a[a.interaction_candidate], b[b.interaction_candidate]
m = a[key + ["interaction_candidate", "replicated", "fold_enrichment", "direction", "q_value", "n_CNV_CNV", "n_mixed", "n_MUT_MUT"]].merge(
    b[key + ["interaction_candidate", "replicated"]], on=key, how="outer", suffixes=("", "_drivers"), indicator=True)
m["status"] = "neither"
m.loc[m.interaction_candidate.eq(True) & m.interaction_candidate_drivers.eq(True), "status"] = "candidate in both"
m.loc[m.interaction_candidate.eq(True) & ~m.interaction_candidate_drivers.eq(True), "status"] = "only with all copy-number calls"
m.loc[~m.interaction_candidate.eq(True) & m.interaction_candidate_drivers.eq(True), "status"] = "only with the driver list"
m = m[m.status != "neither"]
t = v[(v.alteration_set == "SNV+CNV") & (v.check == "TCGA")].set_index(key)
m["tcga_same_direction"] = [t.tcga_direction.get(tuple(r), None) == r2 for r, r2 in zip(m[key].values.tolist(), m.direction)]
m["copy_number_share"] = (m.n_CNV_CNV + m.n_mixed.fillna(0)) / (m.n_CNV_CNV + m.n_mixed.fillna(0) + m.n_MUT_MUT).clip(lower=1)
display(m.groupby("status").agg(pairs=("status", "size"), replicated=("replicated", "mean"),
                                 copy_number_driven=("copy_number_share", lambda s: (s >= .5).mean())).round(2))
lost = m[m.status == "only with all copy-number calls"].assign(effect=lambda d: d.fold_enrichment.apply(lambda x: abs(np.log2(max(x, 1e-9)))))
lost.sort_values("effect", ascending=False)[key + ["direction", "fold_enrichment", "q_value", "replicated"]].head(20).round(3)"""),
md("""> **Reading.** 1,151 candidates are the same under both rules; **418 appear only with all copy-number calls** and 81
> only with the driver list.
>
> - The 418 are mostly **amplicon neighbours**: 87% are exclusive, 80% are copy-number driven, and only **5.5%
>   replicate** outside MSK (37% for the 1,151 shared candidates). Typical examples are *FGF19*, *FGF4*, *WHSC1L1* and
>   *PAK1* against *TP53* or *PIK3CA* in breast cancer (neighbours of *CCND1* and *FGFR1* on 11q13 / 8p11), and
>   *BCL2L1* against *KRAS* or *BRAF* in colorectal cancer (a 20q gain that marks chromosomally unstable tumours).
>   For 241 of them the pair is no longer tested at all: the gene fell below the 3% frequency filter once its
>   uncertain copy-number calls stopped counting.
> - The 81 that appear only with the list are weak (1.2-1.5-fold) and replicate no better (5%), so nothing real is
>   gained there.
> - The list does not remove everything: *FGF3* is an established driver, so it still pairs with *MDM2*, *MYC* or
>   *CD79B* in breast cancer. A "lead gene per amplicon" step could trim that; it is not needed for the conclusion."""),

md("""---
### Figure -- Phase 2b at a glance
**a** published pairs found by each set, or why not · **b** precision and recall on data with planted links ·
**c** candidates keeping their direction in TCGA, against control pairs.  \\n*Drawn by `scripts/figures/phase2b_validation.R`
-> `results/figures/phase2b_validation/figure_phase2b.pdf`*"""),
code("""import subprocess, sys
from IPython.display import Image, display

# Rebuild this figure from the saved outputs: export small CSVs, then run the R script
subprocess.run([sys.executable, "scripts/figures/00_export_figure_data.py"], cwd="..", check=True, capture_output=True)
subprocess.run(["Rscript", "scripts/figures/phase2b_validation.R"], cwd="..", check=True, capture_output=True)
display(Image(filename="../results/figures/phase2b_validation/figure_phase2b.png", width=1000))"""),

md("""---
## Summary

| | SNV | SNV+CNV (current) | **SNV+CNV (drivers)** | SNV+CNV+SV |
|---|---|---|---|---|
| Interaction candidates | 466 | 1,569 | **1,232** | 1,609 |
| MSK hits replicated in other centres | 76% | 65% | **68%** | 65% |
| Known pairs found (of 58) | 35 | 44 | **43** | 50 |
| Precision / recall (planted links) | 95.5% / 49.3% | 96.2% / 46.3% | **99.1% / 49.1%** | 95.7% / 51.9% |
| TCGA: same direction (control) | 93% (71%) | 88% (69%) | **87% (64%)** | 88% (69%) |
| Significant in the opposite direction (TCGA) | 0 | 0 | **0** | 0 |

1. **The test is precise on every set.** The choice of set does not trade away trustworthiness.
2. **More alteration types find more known biology:** copy number adds 9 known pairs, fusions 6 more.
3. **Copy-number calls in non-driver genes are mostly passengers.** The 418 candidates they create replicate at 5.5%.
   Counting only established drivers keeps 43 of the 44 known pairs, raises precision and recall in the simulation,
   and keeps the TCGA agreement.
4. **SNV calls hold up slightly better in independent data**, because mutations are called more consistently than
   copy number.

**Proposal: make SNV+CNV (drivers) the main table.** Fusions (SNV+CNV+SV) add the most known biology at no cost to
precision; they stay a separate set until it is decided whether fusions belong in the main analysis."""),
]
nb = nbf.v4.new_notebook(cells=cells, metadata={"kernelspec": {"name": "fyp-genie", "display_name": "Python (fyp-genie)", "language": "python"}})
nbf.write(nb, "notebooks/02b_phase2_validation.ipynb")
