"""
Script Name: 17_conditional_selection_dnds.py
Purpose : Family D, Coselens (Iranzo et al. 2022) -- borrowed idea. Selection on a gene shows up
          as an excess of protein-changing (nonsynonymous) over silent (synonymous) mutations;
          silent mutations are not selected, so they measure each tumour group's background
          mutation rate. "Conditional selection" of gene B on gene A: B's nonsyn/syn ratio is
          higher (co-selection) or lower (selection against) in A-mutant than in A-wild-type
          tumours. Because the comparison is within B and against B's own silent mutations, it
          is protected from mutation-load confounding by a completely different route than the
          conditional test.
Test    : for Phase 2's replicated mutation-mutation candidates, compare B's nonsyn:syn ratio
          by A status (Fisher's exact test on the 2x2 table of mutation counts), and check the
          direction agrees with Phase 2. Panels report few silent mutations per gene, so many
          pairs will be too sparse -- reported as the main limitation.
Inputs  : data/raw/data_mutations_extended.txt (incl. Silent), data/processed/{clinical_tidy,comutation_pairs}
Outputs : data/processed/conditional_selection_dnds.parquet; printed summary
"""
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

NONSYN = {"Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del",
          "In_Frame_Ins", "Splice_Site", "Nonstop_Mutation", "Translation_Start_Site"}
maf = pd.read_csv("data/raw/data_mutations_extended.txt", sep="\t", low_memory=False,
                  usecols=["Tumor_Sample_Barcode", "Hugo_Symbol", "Variant_Classification", "Mutation_Status"])
maf = maf[maf.Mutation_Status.astype(str).str.upper() == "SOMATIC"]
maf["cls"] = np.where(maf.Variant_Classification.isin(NONSYN), "ns",
                      np.where(maf.Variant_Classification == "Silent", "syn", None))
maf = maf[maf.cls.notna()]
cl = pd.read_parquet("data/processed/clinical_tidy.parquet")
cl = cl[cl.IS_REPRESENTATIVE & (cl.HYPERMUTATION_STATUS == "not_hypermutated")]
maf = maf[maf.Tumor_Sample_Barcode.isin(cl.SAMPLE_ID)].merge(cl[["SAMPLE_ID", "CANCER_TYPE"]],
                                                             left_on="Tumor_Sample_Barcode", right_on="SAMPLE_ID")
print(f"silent mutations available: {(maf.cls == 'syn').sum():,} vs {(maf.cls == 'ns').sum():,} nonsynonymous")

cp = pd.read_parquet("data/processed/comutation_pairs.parquet")
cand = cp[cp.interaction_candidate & cp.replicated & (cp.n_MUT_MUT >= 0.8 * cp.n_both_altered)]
rows = []
for r in cand.itertuples():
    m = maf[maf.CANCER_TYPE == r.Cancer_Type]
    for A, B in [(r.Gene_A, r.Gene_B), (r.Gene_B, r.Gene_A)]:
        a_pos = set(m.loc[(m.Hugo_Symbol == A) & (m.cls == "ns"), "SAMPLE_ID"])
        b = m[m.Hugo_Symbol == B]
        ns_pos, syn_pos = ((b.cls == "ns") & b.SAMPLE_ID.isin(a_pos)).sum(), ((b.cls == "syn") & b.SAMPLE_ID.isin(a_pos)).sum()
        ns_neg, syn_neg = ((b.cls == "ns") & ~b.SAMPLE_ID.isin(a_pos)).sum(), ((b.cls == "syn") & ~b.SAMPLE_ID.isin(a_pos)).sum()
        if syn_pos + syn_neg < 5:          # too few silent mutations to say anything
            continue
        odds, p = fisher_exact([[ns_pos, syn_pos], [ns_neg, syn_neg]])
        ratio = ((ns_pos + 0.5) / (syn_pos + 0.5)) / ((ns_neg + 0.5) / (syn_neg + 0.5))
        rows.append(dict(Cancer_Type=r.Cancer_Type, A=A, B=B, phase2_direction=r.direction,
                         ns_Amut=ns_pos, syn_Amut=syn_pos, ns_Awt=ns_neg, syn_Awt=syn_neg,
                         selection_ratio=ratio, p=p,
                         dnds_direction="co-occurring" if ratio > 1 else "exclusive"))
t = pd.DataFrame(rows)
t.to_parquet("data/processed/conditional_selection_dnds.parquet", index=False)
print(f"mutation-driven replicated candidates: {len(cand)}; testable directions (>= 5 silent mutations in B): {len(t)}")
if len(t):
    agree = (t.dnds_direction == t.phase2_direction)
    print(f"direction agrees with Phase 2: {agree.mean():.0%}; nominally significant (p < 0.05): {(t.p < .05).sum()}, "
          f"of which same direction: {((t.p < .05) & agree).sum()}")
    print(t.sort_values("p").head(12)[["Cancer_Type", "A", "B", "phase2_direction", "ns_Amut", "syn_Amut", "ns_Awt",
                                       "syn_Awt", "selection_ratio", "p"]].round(4).to_string(index=False))
