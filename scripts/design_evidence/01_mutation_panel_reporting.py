"""
Script Name: 01_mutation_panel_reporting.py
Purpose : Does the mutation side have the same problem as copy number -- genes the
          panel metadata says were tested, but whose mutations were never reported?
          Checks (1) genes 'covered' only by intron/intergenic rows, (2) genes marked
          clinicalReported=False, (3) observed vs expected mutated samples per
          (panel, gene), (4) mutations reported in genes the panel does not list.
Inputs  : data/raw/{genomic_information,data_mutations_extended,data_clinical_sample}.txt
Outputs : printed summary; results/analyses/mutation_panel_reporting.parquet (per panel-gene table)
"""
import numpy as np, pandas as pd
R = "data/raw/"
g = pd.read_csv(R + "genomic_information.txt", sep="\t", low_memory=False,
                usecols=["Hugo_Symbol", "SEQ_ASSAY_ID", "Feature_Type", "includeInPanel", "clinicalReported"])
g = g[g.includeInPanel == True]
pg = g.groupby(["SEQ_ASSAY_ID", "Hugo_Symbol"]).agg(
    has_exon=("Feature_Type", lambda s: (s == "exon").any()),
    clin_true=("clinicalReported", lambda s: (s == True).any()),
    clin_false=("clinicalReported", lambda s: (s == False).any())).reset_index()
print(f"covered (panel, gene) pairs: {len(pg):,}")
print(f"(1) covered only by intron/intergenic rows (no exon sequenced): {(~pg.has_exon).sum():,} pairs in "
      f"{pg.loc[~pg.has_exon, 'SEQ_ASSAY_ID'].nunique()} panels")
print(f"(2) clinicalReported: True {pg.clin_true.sum():,}, only False {(pg.clin_false & ~pg.clin_true).sum():,}, "
      f"unset {(~pg.clin_true & ~pg.clin_false).sum():,}")

cs = pd.read_csv(R + "data_clinical_sample.txt", sep="\t", comment="#", usecols=["SAMPLE_ID", "SEQ_ASSAY_ID", "CANCER_TYPE"])
maf = pd.read_csv(R + "data_mutations_extended.txt", sep="\t", low_memory=False,
                  usecols=["Tumor_Sample_Barcode", "Hugo_Symbol", "Variant_Classification"])
PROT = {"Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del",
        "In_Frame_Ins", "Splice_Site", "Translation_Start_Site", "Nonstop_Mutation"}
m = maf[maf.Variant_Classification.isin(PROT)].merge(cs, left_on="Tumor_Sample_Barcode", right_on="SAMPLE_ID")
m = m[["SAMPLE_ID", "Hugo_Symbol", "SEQ_ASSAY_ID", "CANCER_TYPE"]].drop_duplicates(["SAMPLE_ID", "Hugo_Symbol"])

# (4) calls in genes the sample's panel does not list at all
inpanel = m.merge(pg[["SEQ_ASSAY_ID", "Hugo_Symbol"]], how="left", indicator=True)
out = inpanel[inpanel._merge == "left_only"]
print(f"(4) protein-changing calls in genes NOT on the sample's panel: {len(out):,} of {len(m):,} "
      f"({len(out)/len(m):.2%}); panels: {out.SEQ_ASSAY_ID.nunique()}")
print(out.groupby("SEQ_ASSAY_ID").size().sort_values(ascending=False).head(6).to_string())

# (3) observed vs expected mutated samples per (panel, gene), expected from the gene's
# cancer-type-specific rate on all OTHER panels that cover it
n_ct = cs.groupby(["SEQ_ASSAY_ID", "CANCER_TYPE"]).size().rename("n").reset_index()
cov_ct = pg[["SEQ_ASSAY_ID", "Hugo_Symbol"]].merge(n_ct)                              # tested samples
mut_ct = m.merge(pg[["SEQ_ASSAY_ID", "Hugo_Symbol"]]).groupby(["SEQ_ASSAY_ID", "Hugo_Symbol", "CANCER_TYPE"]).size().rename("k").reset_index()
t = cov_ct.merge(mut_ct, how="left").fillna({"k": 0})
tot = t.groupby(["Hugo_Symbol", "CANCER_TYPE"])[["n", "k"]].sum().rename(columns={"n": "N", "k": "K"})
t = t.join(tot, on=["Hugo_Symbol", "CANCER_TYPE"])
t["rate_other"] = (t.K - t.k) / (t.N - t.n).replace(0, np.nan)
t["exp"] = t.n * t.rate_other.fillna(0)
pgk = t.groupby(["SEQ_ASSAY_ID", "Hugo_Symbol"]).agg(n=("n", "sum"), obs=("k", "sum"), exp=("exp", "sum")).reset_index()
pgk = pgk.merge(pg)
silent = pgk[(pgk.exp >= 10) & (pgk.obs == 0)]
under = pgk[(pgk.exp >= 10) & (pgk.obs < 0.1 * pgk.exp)]
print(f"\n(3) (panel, gene) pairs with >= 10 mutated samples expected: {(pgk.exp >= 10).sum():,}")
print(f"    ZERO observed: {len(silent):,} pairs in {silent.SEQ_ASSAY_ID.nunique()} panels, "
      f"{silent.n.sum():,} sample-gene cells")
print(f"    < 10% of expected: {len(under):,} pairs in {under.SEQ_ASSAY_ID.nunique()} panels")
print("    of the zero-observed: no exon rows", (~silent.has_exon).sum(), "| clinicalReported only False",
      (silent.clin_false & ~silent.clin_true).sum())
print(silent.sort_values("exp", ascending=False).head(25)[["SEQ_ASSAY_ID", "Hugo_Symbol", "n", "exp", "obs", "has_exon", "clin_true", "clin_false"]].to_string(index=False))
print("\nzero-observed pairs per panel:")
print(silent.groupby("SEQ_ASSAY_ID").agg(genes=("Hugo_Symbol", "size"), exp=("exp", "sum")).sort_values("exp", ascending=False).head(12).to_string())
# the other direction: does the intron-only / clinicalReported flag predict silence?
e10 = pgk[pgk.exp >= 10]
for name, mask in [("no exon rows", ~e10.has_exon), ("clinicalReported only False", e10.clin_false & ~e10.clin_true),
                   ("exon + reported/unset", e10.has_exon & ~(e10.clin_false & ~e10.clin_true))]:
    s = e10[mask]
    print(f"  {name:<28} pairs {len(s):>6,}  obs/exp {s.obs.sum() / max(s.exp.sum(), 1):.2f}  zero-observed {(s.obs == 0).mean():.1%}")
pgk.to_parquet("results/analyses/mutation_panel_reporting.parquet")
