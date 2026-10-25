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

## 2. Results by aim

### Aim 1 -- lookup tool (Phase 3; Fig 6)
- Five functions: counts with intervals, exact combinations, an on-demand test for rare pairs, partner
  prediction (beats base rates in every cancer type; similarity-based matching does not) and CIViC
  treatment evidence.
- Example: BRAF V600E + NRAS in melanoma -- 7 patients, 27x rarer than expected, CIViC level B
  resistance to BRAF inhibitors.
- Gene frequencies agree between MSK and other centres (r = 0.90-1.00).

### Aim 2 -- biology (Phases 2 and 4; Figs 3-5)
- **6,697 gene pairs** tested in 51 cancer types (112,779 non-hypermutated patients) ->
  **1,584 interaction candidates** (426 co-occurring, 1,158 exclusive); **486 replicated** in centres
  other than MSK; 1,276 also pass the copy-number burden check.
- **Pathways and context:** exclusivity is 3-4x concentrated within one pathway, co-occurrence is not
  pathway-structured; a third of interactions change between cancer types (47 reverse, e.g. KRAS/TP53
  together in pancreas, apart in lung).
- **Variants and smoking:** the smoking signature explains part of KRAS G12C/STK11 and all of
  G12C/RBM10.
- **Fusions:** 109 significant fusion pairs -- ALK/ROS1/RET avoid EGFR/KRAS/STK11 (ALK + STK11:
  0 patients vs 60 expected); ERG avoids SPOP/FOXA1 and goes with PTEN loss in prostate cancer.
- **Is anything new?** All 386 robust pairs were checked against the literature and the SelectSim
  atlas (Nat Genet 2026; GENIE v15 MSK + DFCI, mutations only). The 41 found in neither are amplicon
  passengers or subtype markers; our best candidate (CTCF vs TP53, endometrial) is in SelectSim's
  atlas and we confirm it under stricter adjustment. All robust pairs SelectSim reports agree in
  direction.

### Can the results be trusted? (Fig 2; notebook 05)
- **External validation in TCGA** (independent patients, whole exome): replicated candidates keep the
  same direction in **94%**, are significant in TCGA **46%** of the time vs 7% for GENIE-null pairs,
  and **none reverses** (KEAP1/STK11, EGFR/KRAS, CCND1/FGFR1, CDKN2A/EGFR confirmed).
- **Method benchmark** on identical data: Fisher (cBioPortal) calls 11% of null pairs significant, a
  DISCOVER-style test 33% (18% after FDR), ours 2% (0%); our calls also replicate best in TCGA (37% vs
  29% for Fisher).
- **Independent route:** gene B's ratio of protein-changing to silent mutations depends on gene A in
  the direction Phase 2 predicts (84%; 100% of the significant cases).

## 3. What changed since the last meeting, and why

| Change | Why (evidence on GENIE) |
|---|---|
| Pair test respects each tumour's own alteration count (conditional test) instead of Fisher | On null data from real patients, Fisher calls 6-15% of pairs significant at p < 0.05 (should be 5%); ours 1.3-2.9%. The DISCOVER-style alternative (and the Rediscover R package) is biased towards exclusivity on panel-sized gene sets (46-69% false positives) |
| Main analysis on non-hypermutated patients (your suggestion); hypermutated analysed separately | Replication of findings across centres rose from 44% to 66%; splitting alone is not enough, so the conditional test is kept too |
| Ultramutated tumours (TMB >= 100) excluded from the hypermutated group | Mixing them in made 21.6% of endometrial pairs "significant"; separated, 0.3% |
| Gene rates fitted within detailed subtypes | Removes textbook subtype artefacts (pancreatic MEN1 vs SMAD4, glioma CIC vs PTEN, breast CDH1 vs ERBB2) |
| Neighbouring genes removed only when copy number drives the pair (your suggestion) | Restores KEAP1/STK11/SMARCA4 (lung) and PBRM1/SETD2 (kidney) |
| Testability read from the data, not panel metadata -- copy number, mutations, now fusions | Metadata is wrong for all three: 20 panels claim copy number they never released, 8 released no mutations, 9 claim fusion detection with no fusions released |
| Gene fusions added as a third alteration type | 985 lung-cancer patients with ALK/ROS1/RET fusions (9,619 driver fusions overall) were invisible before; clinically one of the most common questions |
| Copy-number burden checked as a sensitivity analysis | 4 of 25 copy-number-driven pairs are artefacts of unstable genomes (cf. CDCOCA 2011, Zack 2013); building burden into the main test over-adjusts (loses known CCND1/FGFR1), so it is a flag, not the default |
| Novelty audit against the literature and the SelectSim atlas | Separates genuinely unreported pairs from amplicon passengers and subtype markers before anything is called "new" |

## 4. Decisions that need your view

1. **Positioning against SelectSim (Nat Genet 2026).** It covers mutation-driven pairs in MSK and
   DFCI. Should the thesis lead with what it does not cover -- copy-number and fusion partners,
   all-centre panel-aware statistics, artefact audits and the clinical tool?
2. **Thesis emphasis:** the tool, or the biology? This decides what to polish first.
3. **Main analysis = non-hypermutated patients.** Agree that hypermutated tumours are reported
   descriptively only?
4. **Subtype labels.** GENIE's detailed labels differ by centre ("Glioblastoma" vs "Glioblastoma
   Multiforme"). Harmonise via OncoTree codes?
5. **Weak-effect threshold** for reporting (currently 1.5-fold).
6. **Data access:** apply for GENIE BPC (treatment and outcome) and an OncoKB academic token?

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
| Figures (6 main, 6 supplementary; PDF and PNG) | `results/figures/` |
| Validation and method comparison | `notebooks/05_validation_and_methods.ipynb` |
| Biology, including the novelty audit | `notebooks/04_phase4_insight.ipynb` (section 9) |
| The lookup tool | `notebooks/03_phase3_lookup_tool.ipynb` |
| Why each method was chosen | `docs/literature_review.md` |
