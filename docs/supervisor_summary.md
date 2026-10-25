# Progress summary for supervisor review

**Project:** co-mutation patterns in AACR GENIE v19.0 -- (1) a lookup tool for clinicians,
(2) scientific insight into rare and cancer-type-specific combinations.
Details: `README.md`; method justification: `docs/literature_review.md`.

## 1. What changed since the last meeting, and why

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

## 2. Headline results

- **6,697 gene pairs** tested in 51 cancer types (112,779 non-hypermutated patients) ->
  **1,584 interaction candidates** (426 co-occurring, 1,158 exclusive); **486 replicated** in centres other than MSK;
  1,276 also pass the copy-number burden check.
- **External validation in TCGA** (independent patients, whole exome): replicated candidates keep
  the same direction in **94%**, are significant in TCGA **46%** of the time vs 7% for GENIE-null
  pairs, and **none reverses** (KEAP1/STK11, EGFR/KRAS, CCND1/FGFR1, CDKN2A/EGFR confirmed).
- **Method benchmark:** every published method family run on the same data -- Fisher (cBioPortal)
  calls 11% of null pairs significant, a DISCOVER-style test 33% (18% after FDR), ours 2% (0%);
  our calls also replicate best in TCGA (37% vs 29% for Fisher).
- **Fusions:** 109 significant fusion pairs -- ALK/ROS1/RET avoid EGFR/KRAS/STK11 (ALK + STK11:
  0 patients vs 60 expected); ERG avoids SPOP/FOXA1 and goes with PTEN loss in prostate cancer.
- **Independent route:** gene B's ratio of protein-changing to silent mutations depends on gene A
  in the direction Phase 2 predicts (84%; 100% of the significant cases).
- **Biology (Phase 4):** exclusivity is 3-4x concentrated within one pathway, co-occurrence is not
  pathway-structured; a third of interactions change between cancer types (47 reverse, e.g.
  KRAS/TP53 together in pancreas, apart in lung); smoking signature explains part of KRAS G12C/STK11
  and all of G12C/RBM10.
- **Is anything new?** A literature screen of the 386 most robust pairs: 343 already reported; the
  43 unreported are mostly amplicon passengers or subtype markers. One candidate survives every
  adjustment and is backed by TCGA and the dN/dS route: CTCF mutation is avoided in TP53-mutant
  endometrial cancer (22 vs 121 expected). Needs molecular-class labels to confirm.
- **Lookup tool (Phase 3):** counts, exact combinations, on-demand test for rare pairs, partner
  prediction (beats base rates in every cancer type; similarity-based matching does not) and CIViC
  treatment evidence -- e.g. BRAF V600E + NRAS in melanoma: 7 patients, 27x rarer than expected,
  CIViC level B resistance to BRAF inhibitors.

## 3. Decisions that need your view

1. **Main analysis = non-hypermutated patients.** Agree that hypermutated tumours are reported
   descriptively only?
2. **Subtype labels.** GENIE's detailed labels differ by centre ("Glioblastoma" vs "Glioblastoma
   Multiforme"). Harmonise via OncoTree codes?
3. **Weak-effect threshold** for reporting (currently 1.5-fold).
4. **Thesis emphasis:** the tool, or the biology? This decides what to polish first.
5. **Data access:** apply for GENIE BPC (treatment and outcome) and an OncoKB academic token?

## 4. Limitations to discuss

Panels that list more genes than they report (detectable only statistically); no clonality on
panel data; GENIE has no treatment data, so the tool reports evidence, not predicted response;
TCGA is ~20x smaller than GENIE, so many pairs cannot be tested there.
