"""
Script Name: build_notebook05.py
Purpose : Generate notebooks/05_validation_and_methods.ipynb, which collects the validation and
          literature-method results (scripts/analyses/12-18) in one readable place.
Inputs  : none (the notebook itself reads data/processed/*.parquet)
Outputs : notebooks/05_validation_and_methods.ipynb (unexecuted; run with nbconvert --execute)
Note    : the notebook only reads saved outputs; re-run scripts 12-18 to regenerate them.
"""
import nbformat as nbf

md, code = nbf.v4.new_markdown_cell, nbf.v4.new_code_cell
cells = [
md("""# Phase 2 validation -- does the method hold up, and what do the published methods add?

Phases 2-4 rest on one choice: a **conditional pair test** that respects each tumour's own
alteration count, fitted within subtypes. This notebook collects the checks of that choice and
the results of borrowing each method family from the literature review
(`docs/literature_review.md` §4-§8):

| § | Check | Literature family | Script |
|---|---|---|---|
| 1 | Independent patients (TCGA) | external validation | `14_tcga_external_validation.py` |
| 2 | Head-to-head benchmark on the same data | A Fisher, B DISCOVER, D SELECT | `15_method_benchmark.py` |
| 3 | Exclusive gene *sets*, not pairs | C Dendrix / MEMo | `16_exclusive_gene_sets.py` |
| 4 | A route that ignores co-occurrence counts | E selection (dN/dS, Coselens) | `17_conditional_selection_dnds.py` |
| 5 | Borrowing strength across cancer types | F hierarchical models | `18_hierarchical_shrinkage.py` |
| 6 | Copy-number burden confound | -- (found here) | `12_burden_adjusted_retest.py`, Phase 2 §6d |
| 7 | Gene fusions | -- | Phase 2 §6e |

The scripts take minutes to hours (TCGA download, curveball nulls); this notebook only reads
their saved outputs."""),
code("""import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control

P = "../data/processed/"
pd.set_option("display.width", 160)
tv = pd.read_parquet(P + "tcga_validation.parquet")"""),

md("""## 1. External validation in TCGA

TCGA PanCancer Atlas: different patients, whole-exome sequencing, different centres and callers.
The same pairs are re-tested with the same conditional test after applying the Phase 1 filters
as closely as TCGA allows (non-hypermutated patients; TCGA study as the subtype stratum).
TCGA has about 1/20 of GENIE's patients, so many pairs cannot be tested there; the question is
whether the ones that can **keep their direction**."""),
code("""t = tv[tv.tcga_p.notna()]
rows = []
for name, sel in [("replicated candidates", t.interaction_candidate & t.replicated),
                  ("all candidates", t.interaction_candidate),
                  ("not significant in GENIE (control)", t.q_value >= 0.05)]:
    g = t[sel]; same = g.tcga_direction == g.direction
    rows.append(dict(group=name, testable_in_TCGA=len(g), same_direction=same.mean(),
                     TCGA_p05_same_dir=((g.tcga_p < .05) & same).mean(),
                     TCGA_p05_opposite=int(((g.tcga_p < .05) & ~same).sum())))
pd.DataFrame(rows).round(3)"""),
code("""# Textbook pairs: GENIE vs TCGA effect side by side
known = [("Non-Small Cell Lung Cancer", "KEAP1", "STK11"), ("Non-Small Cell Lung Cancer", "EGFR", "KRAS"),
         ("Breast Cancer", "CCND1", "FGFR1"), ("Glioma", "CDKN2A", "EGFR")]
k = pd.concat([t[(t.Cancer_Type == c) & (t.Gene_A.isin([a, b])) & (t.Gene_B.isin([a, b]))] for c, a, b in known])
k.assign(tcga_fold=(k.tcga_n_both / k.tcga_expected).round(2))[
    ["Cancer_Type", "Gene_A", "Gene_B", "fold_enrichment", "q_value", "tcga_n_both", "tcga_expected", "tcga_fold", "tcga_p"]].round(4)"""),
md("""**Reading.** Replicated candidates keep their direction in TCGA almost always and are
significant there several times more often than the control pairs; none is significant the other
way. The out-of-centre replication flag (Phase 2 §6c) is therefore a useful confidence label, and
the control row shows the agreement is not an artefact of shared gene frequencies."""),

md("""## 2. Head-to-head benchmark of the method families

All methods run on the **same** patients, genes and pairs (9 large cancer types):

- **A Fisher** -- what cBioPortal and maftools report;
- **B1 DISCOVER-style** -- per-tumour rates plugged in as known (Canisius 2016);
- **B2 conditional** -- our test without subtype strata; **B3** -- with strata (Phase 2);
- **D SELECT-style** -- weighted mutual information (Mina 2017), a ranking score only.

**(1) Null calibration** needs a curveball-permuted null (matrices that keep every gene's and every
tumour's alteration count but break any real pairing). It is computed in the script, not
saved; from `results/logs/15_method_benchmark.log` -- share of null pairs with p < 0.05 (target 5%):

| Method | p < 0.05 | after FDR |
|---|---|---|
| A Fisher | 11.4% | 0.6% |
| B1 DISCOVER-style | 32.7% | 18.1% |
| B2 conditional | 2.0% | 0.0% |
| B3 conditional + subtype | 1.9% | 0.0% |

**(2)-(3)** below are recomputed from the saved real-data results."""),
code("""mb = pd.read_parquet(P + "method_benchmark.parquet")
M = {"A Fisher": ("p_fisher", "dir_fisher"), "B1 DISCOVER-style": ("p_discover", "dir_discover"),
     "B2 conditional": ("p_conditional", "dir_conditional"), "B3 conditional + subtype": ("p_cond_subtype", "dir_cond_subtype")}
rows = []
for name, (pc, dc) in M.items():
    q = false_discovery_control(mb[pc].clip(1e-300, 1))
    s = mb[(q < .05) & mb.tcga_p.notna()]
    rows.append(dict(method=name, significant=int((q < .05).sum()), testable_in_TCGA=len(s),
                     replicate_in_TCGA=((s.tcga_p < .05) & (s.tcga_direction == s[dc])).mean(),
                     significant_opposite=int(((s.tcga_p < .05) & (s.tcga_direction != s[dc])).sum())))
pd.DataFrame(rows).round(3)"""),
code("""# (3) Ranking: of each method's K strongest pairs, how many replicate in TCGA?
t = mb[mb.tcga_p.notna()]
scores = {name: (t[pc], t[dc]) for name, (pc, dc) in M.items()}
scores["D SELECT-style"] = (-t["mi_select"], t["dir_select"])
pd.DataFrame({name: [((top.tcga_p < .05) & (top.tcga_direction == top._d)).mean()
                     for K in (50, 100, 200)
                     for top in [t.assign(_s=s.to_numpy(), _d=d.to_numpy()).nsmallest(K, "_s")]]
              for name, (s, d) in scores.items()}, index=["top 50", "top 100", "top 200"]).round(3)"""),
md("""**Reading.**
- **Calibration is where the methods differ.** Fisher roughly doubles the false-positive rate
  and DISCOVER-style plug-in is far worse on panel data (it treats estimated rates as known);
  the conditional test is conservative, and multiple-testing correction then removes the noise.
- **Precision follows calibration.** Our significant list is the shortest but the most
  reproducible in TCGA (37% vs 29% for Fisher).
- **Ranking does not.** The strongest 50-100 pairs replicate at a similar rate whatever the method;
  the very strong signals are found by everyone. The method choice matters for the *long tail* of
  moderate effects, which is exactly where a clinical lookup tool is queried."""),

md("""## 3. Exclusive gene sets (Dendrix / MEMo idea)

Pairs miss groups of three or more genes that substitute for each other (one pathway, one hit).
Dendrix scores a set by coverage minus overlap, which favours frequent genes. The load-aware
version asks whether the set overlaps **less than expected given each tumour's alteration
count** (same logic as the pair test). Benchmark: how often the top-20 sets per cancer type fall
in one Sanchez-Vega (2018) pathway."""),
code("""gs = pd.read_parquet(P + "exclusive_gene_sets.parquet")
summ = []
for ct, g in gs.groupby("Cancer_Type"):
    d = g.nlargest(20, "dendrix_W"); l = g[g.coverage >= 0.10].nsmallest(20, "z_exclusivity")
    summ.append(dict(Cancer_Type=ct, dendrix_all3=d.same_pathway_all3.mean(), loadaware_all3=l.same_pathway_all3.mean(),
                     dendrix_any2=d.same_pathway_any2.mean(), loadaware_any2=l.same_pathway_any2.mean(),
                     dendrix_top=d.iloc[0].genes, loadaware_top=l.iloc[0].genes))
summ = pd.DataFrame(summ)
display(summ.round(2))
summ[["dendrix_all3", "loadaware_all3", "dendrix_any2", "loadaware_any2"]].mean().round(2)"""),
md("""**Reading.** Load-aware sets are about three times as often whole-pathway modules; Dendrix's
top sets are dominated by the most frequent genes (e.g. TP53 with anything rare). The pair-level
finding -- exclusivity is pathway-structured -- holds for sets as well."""),

md("""## 4. An independent route: does gene A change selection on gene B? (dN/dS idea)

Co-occurrence counts can be distorted by panels and copy-number artefacts. The ratio of
protein-changing to silent mutations in gene B does not depend on those: if B is **less useful**
once A is mutated, B's protein-changing share drops in A-mutant tumours (Coselens, Iranzo 2022).
Only mutation-driven replicated candidates with enough silent mutations can be tested."""),
code("""dn = pd.read_parquet(P + "conditional_selection_dnds.parquet")
agree = dn.dnds_direction == dn.phase2_direction
print(f"testable directions: {len(dn)}; agree with Phase 2: {agree.mean():.0%}; "
      f"p < 0.05: {(dn.p < .05).sum()}, of which agree: {((dn.p < .05) & agree).sum()}")
dn.sort_values("p").head(10).round(4)"""),
md("""**Reading.** Direction agrees in most testable cases and in every significant one, from data
the pair test never used. The interactions are therefore selection, not counting artefacts.
Limitation: GENIE releases few silent mutations (panels filter them), so only a minority of
pairs can be checked."""),

md("""## 5. Borrowing strength across cancer types (hierarchical idea)

Each pair is estimated separately in every cancer type; small cancer types give noisy
estimates. A DerSimonian-Laird random-effects model per gene pair pulls each estimate towards
the pair's average across cancer types, by an amount set by its own noise and by how much the
pair truly varies. Judged against the TCGA estimate of the same pair."""),
code("""s = pd.read_parquet(P + "shrunk_effects.parquet")
t = s[s.tcga_p.notna()].copy()
t["theta_tcga"] = np.log((t.tcga_n_both + 0.5) / (t.tcga_expected + 0.5))
n_types = s.key.value_counts()
rows = []
for name, sel in [("all testable", t.index == t.index), ("noisy (< 30 co-altered)", t.n_both_altered < 30),
                  ("pair in >= 3 cancer types", t.key.map(n_types) >= 3)]:
    d = t[sel]
    rows.append(dict(group=name, n=len(d), corr_raw=np.corrcoef(d.theta, d.theta_tcga)[0, 1],
                     corr_shrunk=np.corrcoef(d.theta_shrunk, d.theta_tcga)[0, 1],
                     sq_err_raw=((d.theta - d.theta_tcga) ** 2).mean(), sq_err_shrunk=((d.theta_shrunk - d.theta_tcga) ** 2).mean()))
pd.DataFrame(rows).round(3)"""),
md("""**Reading.** A small, consistent gain, largest for the noisy estimates. Not adopted for the
headline calls (it would blur the cancer-type-specific and reversing pairs that Phase 4 studies),
but `theta_shrunk` is a better *effect size to display* for rare pairs in the lookup tool."""),

md("""## 6. Copy-number burden: the confound not covered by the literature methods

Tumours with many amplifications and deletions co-amplify unrelated genes. The conditional test
counts alterations per tumour, but within the tested genes only. Phase 2 §6d adds a burden
stratum as a sensitivity column; script 12 re-tests every copy-number-driven co-occurring
candidate with burden as a covariate."""),
code("""cp = pd.read_parquet(P + "comutation_pairs.parquet")
c = cp[cp.interaction_candidate]
print(f"candidates {len(c):,}; burden-robust {c.burden_robust.sum():,} ({c.burden_robust.mean():.0%})")
br = pd.read_parquet(P + "burden_adjusted_retest.parquet")
display(br.verdict.str.split(":").str[0].value_counts())
br[~br.verdict.str.startswith("survives")].round(4)"""),

md("""## 7. Gene fusions

Fusions are tested only on panels that both detect them and release them (read from the data).
Fusion screen = pairs between a driver fusion and the Phase 2 genes."""),
code("""fu = pd.read_parquet(P + "comutation_pairs_fusions.parquet")
sig = fu[fu.q_value < .05]
print(f"fusion pairs tested {len(fu)}; significant {len(sig)} "
      f"({(sig.direction == 'exclusive').sum()} exclusive, {(sig.direction == 'co-occurring').sum()} co-occurring)")
sig.sort_values("q_value").head(15)[["Cancer_Type", "Gene_A", "Gene_B", "n_both_altered", "expected_both",
                                     "fold_enrichment", "direction", "q_value"]].round(4)"""),

md("""## Discussion -- what trying the published methods taught us

1. **The main design choice is confirmed from three independent directions**: null calibration
   (only the conditional test stays below 5%), independent patients (TCGA: 94% same direction,
   none opposite) and an independent signal (dN/dS agrees in every significant case).
2. **The advantage is in the long tail, not the headline pairs.** All methods find KEAP1/STK11 or
   EGFR/KRAS. They disagree on moderate effects, which is where Fisher and DISCOVER-style tests
   produce false calls -- and where a clinician's lookup query usually lands.
3. **Exclusivity is a pathway property; co-occurrence is not.** Seen for pairs (Phase 4) and now for
   gene sets: load-aware sets are whole-pathway three times as often as Dendrix sets.
4. **Copy number needs its own artefact checks.** Burden, neighbouring genes and shared events
   explain a share of co-occurrence that no published method addresses for panel data.
5. **Frequency and enrichment answer different questions.** ALK-fusion lung cancers often carry
   TP53 mutations, yet less often than expected: the tool must show both.
6. **Borrowed but not adopted:** shrinkage (small gain, blurs specificity) and SELECT-style
   ranking (no better than p-value ranking). Both are reported, not used for headline calls."""),
]
nb = nbf.v4.new_notebook(cells=cells, metadata={"kernelspec": {"name": "python3", "display_name": "Python 3", "language": "python"}})
nbf.write(nb, "notebooks/05_validation_and_methods.ipynb")
