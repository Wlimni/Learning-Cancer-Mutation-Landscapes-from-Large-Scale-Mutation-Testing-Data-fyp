# Progress summary for supervisor review

**Project:** co-mutation patterns in AACR GENIE v19.0 -- Aim 1: a lookup tool for clinicians;
Aim 2: scientific insight into rare and cancer-type-specific combinations.
**Status (25 Oct 2026):** all four phases built, validated and committed; the v3 design has not yet
been reviewed by you. Details: `README.md`; method justification: `docs/literature_review.md`;
end-to-end walkthrough: `docs/FYP_onboarding_guide.docx`.

## 1. Bottom line

- **The method is sound.** Of all published method families, only our conditional test is calibrated
  on GENIE panel data, and its findings replicate in TCGA (94% same direction, none reversed).
- **The tool works** for genes, variants, copy-number changes and fusions, including rare combinations
  that the statistical screen cannot test.
- **The biology is mostly confirmation, and we checked that rather than assumed it.** Mutation-driven
  driver pairs in common cancers are almost all already reported; our distinct contribution is
  copy-number and fusion partners, panel-aware statistics across all centres, the artefact audits and
  the tool. This positioning needs your view (decision 1).

## 2. Since our last meeting (Phase 2 results)

### Your suggestions, implemented

| Suggestion | What we did | Effect on GENIE |
|---|---|---|
| Exclude hypermutated tumours from the main analysis | Main analysis on non-hypermutated patients (TMB < 10 on panels >= 1 Mb); hypermutated tumours (TMB 10-100) analysed separately and reported descriptively; ultramutated (TMB >= 100) excluded | Replication across centres rose from 44% to 66%; pooling ultramutated tumours had made 21.6% of endometrial pairs "significant" (0.3% once separated) |
| Remove neighbouring genes only when copy number drives the pair (two point mutations are two events) | Neighbour rule restricted to copy-number pairs | Restores KEAP1/STK11/SMARCA4 (lung) and PBRM1/SETD2 (kidney) |
| Fisher's test seemed to call too much ("half significant sounds high") | Replaced by the conditional test; FDR per cancer type | Significant pairs 72% -> 27% on a driver-only list; on null data Fisher 6-15% at p < 0.05, ours 1.3-2.9% |

Splitting off hypermutated tumours alone was not enough (Fisher stays inflated afterwards), so the
conditional test is kept as well.

### Other changes to Phase 2

| Change | Why (evidence on GENIE) |
|---|---|
| Gene rates fitted within detailed subtypes | Removes textbook subtype artefacts (pancreatic MEN1 vs SMAD4, glioma CIC vs PTEN, breast CDH1 vs ERBB2) |
| Testability read from the data, not panel metadata -- copy number, mutations, fusions | Metadata is wrong for all three: 20 panels claim copy number they never released, 8 released no mutations, 9 claim fusion detection with no fusions released |
| Gene fusions added as a third alteration type | 985 lung-cancer patients with ALK/ROS1/RET fusions (15,847 driver fusions overall, defined from Mitelman + OncoKB) were invisible before |
| Copy-number burden checked as a sensitivity flag | 4 of 25 copy-number-driven pairs are artefacts of unstable genomes; building burden into the main test over-adjusts (loses known CCND1/FGFR1) |
| DISCOVER-style test and the Rediscover R package tested and rejected | Biased towards exclusivity on panel-sized gene sets (46-69% false positives on null data) |

### New work since the meeting

- **Phase 3 (Aim 1):** the lookup tool, five functions (Section 3).
- **Phase 4 (Aim 2):** pathway structure, cancer-type differences, variant-level partners and smoking,
  rare combinations, the copy-number burden confound.
- **Validation (notebook 02b):** precision and recall -- on data without links, in a
  simulation with planted links (three scenarios), recovery of 29 textbook pairs, and TCGA replication --
  each compared with Fisher, DISCOVER-style and SELECT on identical data.
- **Novelty audit:** every robust pair checked against the literature and the SelectSim atlas
  (notebook 04, section 9).
- **Figures:** one combined R figure per phase (up to four panels, titled with its main message), drawn from inside each notebook.

## 3. Results by aim

### Aim 1 -- lookup tool (Phase 3 figure)
- Five functions: counts with intervals, exact combinations, an on-demand test for rare pairs, partner
  prediction (beats base rates in every cancer type; similarity-based matching does not) and CIViC
  treatment evidence.
- Example: BRAF V600E + NRAS in melanoma -- 7 patients, 27x rarer than expected, CIViC level B
  resistance to BRAF inhibitors.
- Gene frequencies agree between MSK and other centres (r = 0.90-1.00).

### Aim 2 -- biology (Phase 2 and Phase 4 figures)
- **6,871 gene pairs** tested in 52 cancer types (113,042 non-hypermutated patients) ->
  **1,609 interaction candidates** (436 co-occurring, 1,173 exclusive); **494 replicated** in centres
  other than MSK; 1,300 also pass the copy-number burden check. New positive control: the TERT promoter avoids
  ATRX / IDH1 and co-occurs with EGFR / PTEN / CDKN2A in glioma, and co-occurs with BRAF in thyroid cancer.
- **Pathways and context:** exclusivity is 3-4x concentrated within one pathway, co-occurrence is not
  pathway-structured; a third of interactions change between cancer types (47 reverse, e.g. KRAS/TP53
  together in pancreas, apart in lung).
- **Variants and smoking:** the smoking signature explains part of KRAS G12C/STK11 and all of
  G12C/RBM10.
- **Fusions:** 109 significant fusion pairs -- ALK/ROS1/RET avoid EGFR/KRAS/STK11 (ALK + STK11:
  0 patients vs 60 expected); ERG avoids SPOP/FOXA1 and goes with PTEN loss in prostate cancer.
- **Is anything new?** All 387 robust pairs were checked against the literature and the SelectSim
  atlas (Nat Genet 2026; GENIE v15 MSK + DFCI, mutations only). The 41 found in neither are amplicon
  passengers or subtype markers; our best candidate (CTCF vs TP53, endometrial) is in SelectSim's
  atlas and we confirm it under stricter adjustment. All robust pairs SelectSim reports agree in
  direction.

### Can the results be trusted? (Phase 2b figure; notebook 02b)
- **Data with no real links** (every call is wrong): Fisher (cBioPortal) calls 11% of pairs
  significant and a DISCOVER-style test 33% (18% after FDR); ours 2% (0%), SELECT (Mina 2017,
  re-implemented) 4%. No candidate is reversed in TCGA.
- **Precision and recall (new simulation):** links of known strength planted into data without links,
  three scenarios. On realistic data ours has precision 99% and recall 49%; Fisher and DISCOVER-style reach
  recall 60-63% only at a precision of 4-8%. At the same precision ours finds the most real links (recall
  58% at 95% precision vs 13% for SELECT). It misses weak
  co-occurrence and part of the exclusive links. A hidden patient group the strata do not capture breaks
  every method, ours included.
- **Known biology:** 26 of 29 textbook pairs found (Fisher 28, DISCOVER-style and SELECT 27). The 10
  strongest exclusive pairs are all established; of the 10 strongest co-occurring pairs, 4 are not yet
  explained (BRAF+SETD2 and CTNNB1+SMAD4 in lung, NF1+PTPN11 in glioma, APC+SPOP in prostate).
- **TCGA** (independent patients, whole exome): candidates keep their direction in 88% (control pairs
  69%), effect sizes correlate (Spearman 0.74), and agreement rises with GENIE signal strength.
- **Mutation-load confounding** (van de Haar 2019's diagnostic): absent once hypermutated tumours are
  analysed separately -- direct evidence for your suggestion to exclude them (notebook 02c, Section E).

## 4. Decisions that need your view

1. **Positioning against SelectSim (Nat Genet 2026).** It covers mutation-driven pairs in MSK and
   DFCI. Should the thesis lead with what it does not cover -- copy-number and fusion partners,
   all-centre panel-aware statistics, artefact audits and the clinical tool?
2. **Thesis emphasis:** the tool, or the biology? This decides what to polish first.
3. **Main analysis = non-hypermutated patients.** Agree that hypermutated tumours are reported
   descriptively only? (Both replication and the mutation-load check support it.)
4. **Subtype labels.** GENIE's detailed labels differ by centre ("Glioblastoma" vs "Glioblastoma
   Multiforme"). Harmonise via OncoTree codes?
5. **Weak-effect threshold** for reporting (currently 1.5-fold).
6. **Data access:** apply for GENIE BPC (treatment and outcome) and an OncoKB academic token?
7. **Driver events beyond protein changes (for information).** Besides mutations, amplifications/deletions and
   fusions, Phase 1 now keeps two well-established drivers that the "protein-altering" label misses: **TERT promoter
   hotspots** (C228T/C250T; switch TERT on, like an amplification; a separate feature, testable on the 80 panels
   whose design covers the promoter) and **MET exon 14 splice changes** beyond the canonical sites. Fusions are now
   defined from external sources (Mitelman published pairs, or a recurrent fusion oncogene in Mitelman + OncoKB)
   instead of a hand-typed gene list. OK with this scope?

## 5. Limitations and risks

- **Overlap with published work:** SelectSim already reports the mutation-driven driver pairs; novelty
  has to come from copy number, fusions, rare cancer types and the tool (decision 1).
- Panels that list more genes than they report (detectable only statistically); no clonality on panel
  data.
- GENIE has no treatment data, so the tool reports evidence, not predicted response.
- TCGA is ~20x smaller than GENIE, so many pairs cannot be tested there.

## 6. Next steps (proposed order)

1. Your review of the v3 design and the decisions above.
2. Run the novelty audit on copy-number x mutation pairs, fusion partners and rare cancer types, where
   new biology is most likely.
3. A simple web interface for the lookup tool once we agree what it shows.
4. Optional: OncoKB as a second evidence source; GENIE BPC for outcomes.

## 7. Where to look

| What | Where |
|---|---|
| Figures by phase (PDF and PNG) | `results/figures/phase1_data/` ... `validation/` |
| Validation and method comparison | `notebooks/02b_phase2_validation.ipynb` (side analyses: `02c_phase2_extensions.ipynb`) |
| Biology, including the novelty audit | `notebooks/04_phase4_insight.ipynb` (section 9) |
| The lookup tool | `notebooks/03_phase3_lookup_tool.ipynb` |
| Why each method was chosen | `docs/literature_review.md` |
