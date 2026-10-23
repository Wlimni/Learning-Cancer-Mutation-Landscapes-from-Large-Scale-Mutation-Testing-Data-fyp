"""
Script Name: 16_exclusive_gene_sets.py
Purpose : Family C (MEMo, Ciriello 2012; Dendrix, Vandin 2012; Multi-Dendrix, Leiserson 2013;
          Mutex, Babur 2015): find SETS of genes whose alterations rarely overlap -- usually one
          pathway. Compared on GENIE, per cancer type, for sets of 3 genes:
            Dendrix weight  W = 2*|patients covered| - sum of each gene's patients
                            (high coverage, little overlap; ignores mutation load)
            load-aware      observed vs EXPECTED pairwise overlap within the set, expected from
                            the Phase 2 conditional model (each patient's own alteration count)
Score   : pathway coherence of the top 20 sets (all 3 genes in one TCGA pathway; >= 2 sharing
          one), and how frequent the genes are (Dendrix is known to favour common genes).
Inputs  : data/processed/*.parquet; data/external/pathways/sanchez_vega_2018_pathways.tsv
Outputs : data/processed/exclusive_gene_sets.parquet; printed comparison
"""
import itertools, sys
import numpy as np
import pandas as pd
sys.path.insert(0, "scripts/validation")
from common import load_phase2_functions

F = load_phase2_functions()
P = "data/processed/"
cl = pd.read_parquet(P + "clinical_tidy.parquet")
cl = cl[cl.IS_REPRESENTATIVE & cl.CNA_TESTED & cl.MUT_TESTED & (cl.HYPERMUTATION_STATUS == "not_hypermutated")]
al = pd.read_parquet(P + "alterations_long.parquet", columns=["Sample_ID", "Hugo_Symbol", "Alteration_Type"])
al = al[al.Alteration_Type.isin(["MUT", "CNV"])][["Sample_ID", "Hugo_Symbol"]].drop_duplicates()
cov = pd.read_parquet(P + "panel_gene_coverage.parquet")[["SEQ_ASSAY_ID", "Hugo_Symbol"]].merge(
      pd.read_parquet(P + "cna_gene_panel_coverage.parquet"))
cp = pd.read_parquet(P + "comutation_pairs.parquet")
pw = pd.read_csv("data/external/pathways/sanchez_vega_2018_pathways.tsv", sep="\t").groupby("gene").pathway.apply(set)


def coherence(genes):
    sets = [pw.get(g, set()) for g in genes]
    all3 = bool(set.intersection(*sets)) if all(sets) else False
    any2 = any(sets[i] & sets[j] for i, j in itertools.combinations(range(3), 2))
    return all3, any2


rows = []
for ct in ["Non-Small Cell Lung Cancer", "Breast Cancer", "Colorectal Cancer", "Glioma", "Endometrial Cancer", "Bladder Cancer"]:
    c = cl[cl.CANCER_TYPE == ct]
    genes = sorted(set(cp.loc[cp.Cancer_Type == ct, "Gene_A"]) | set(cp.loc[cp.Cancer_Type == ct, "Gene_B"]))
    per_panel = c.SEQ_ASSAY_ID.value_counts()
    big = per_panel.index[per_panel.cumsum() / per_panel.sum() <= 0.85].tolist() or per_panel.index[:1].tolist()
    common = set.intersection(*[set(cov.loc[cov.SEQ_ASSAY_ID == p, "Hugo_Symbol"]) for p in big])
    genes = [g for g in genes if g in common]
    c = c[c.SEQ_ASSAY_ID.isin(big)]
    a_c = al[al.Sample_ID.isin(c.SAMPLE_ID) & al.Hugo_Symbol.isin(genes)]          # this cancer type only (fast)
    X = (pd.crosstab(a_c.Sample_ID, a_c.Hugo_Symbol).reindex(index=c.SAMPLE_ID, columns=genes, fill_value=0) > 0).values
    freq = X.mean(0)
    top = np.argsort(-freq)[:35]                                  # 35 most altered genes -> 6,545 triples
    det = c.CANCER_TYPE_DETAILED.fillna("NA"); vc = det.value_counts()
    M = F["ConditionalModel"](X, np.ones_like(X), det.where(det.map(vc) >= 50, "other").to_numpy())
    E = np.zeros((X.shape[1], X.shape[1]))
    for a, b in itertools.combinations(sorted(top), 2):
        E[a, b] = E[b, a] = M.joint(a, b).sum()                   # expected co-altered patients
    O = X[:, top].T.astype(int) @ X[:, top].astype(int)
    pos = {g: k for k, g in enumerate(top)}
    for trip in itertools.combinations(sorted(top), 3):
        cols = X[:, list(trip)]
        covered = cols.any(1).sum(); w_dendrix = 2 * covered - cols.sum()
        obs = sum(O[pos[a], pos[b]] for a, b in itertools.combinations(trip, 2))
        exp = sum(E[a, b] for a, b in itertools.combinations(trip, 2))
        z = (obs - exp) / np.sqrt(max(exp, 1e-9))                 # negative = more exclusive than expected
        all3, any2 = coherence([genes[i] for i in trip])
        rows.append(dict(Cancer_Type=ct, genes="/".join(genes[i] for i in trip), coverage=covered / len(X),
                         dendrix_W=w_dendrix, obs_overlap=obs, exp_overlap=exp, z_exclusivity=z,
                         mean_gene_freq=freq[list(trip)].mean(), same_pathway_all3=all3, same_pathway_any2=any2))
    print(f"{ct:<28} patients {len(X):,}  triples scored {len(top) * (len(top) - 1) * (len(top) - 2) // 6:,}", flush=True)

t = pd.DataFrame(rows)
t.to_parquet(P + "exclusive_gene_sets.parquet", index=False)
summary = []
for ct, g in t.groupby("Cancer_Type"):
    d = g.nlargest(20, "dendrix_W"); l = g[g.coverage >= 0.10].nsmallest(20, "z_exclusivity")
    summary.append(dict(Cancer_Type=ct,
                        dendrix_all3=d.same_pathway_all3.mean(), loadaware_all3=l.same_pathway_all3.mean(),
                        dendrix_any2=d.same_pathway_any2.mean(), loadaware_any2=l.same_pathway_any2.mean(),
                        dendrix_gene_freq=d.mean_gene_freq.mean(), loadaware_gene_freq=l.mean_gene_freq.mean(),
                        dendrix_top=d.iloc[0].genes, loadaware_top=l.iloc[0].genes))
s = pd.DataFrame(summary)
print("\nTop-20 sets per cancer type: pathway coherence and gene frequency, Dendrix vs load-aware")
print(s.round(2).to_string(index=False))
print("\nmean over cancer types:", s[["dendrix_all3", "loadaware_all3", "dendrix_any2", "loadaware_any2"]].mean().round(2).to_dict())
