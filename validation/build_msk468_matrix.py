"""
Script Name: build_msk468_matrix.py
Purpose : Build the benchmark matrix every validation script uses: all
          MSK-IMPACT468 samples x the 45 genes most often altered on them.
          One panel means every gene was tested in every sample, so no
          testability masking is involved and external tools (Rediscover)
          see exactly the same input as the Phase 2 test.
Inputs  : data/processed/alterations_long.parquet, clinical_tidy.parquet
Outputs : data/validation/msk468_matrix.npy   (samples x genes, bool)
          data/validation/msk468_genes.json   (gene order, alphabetical)
          data/validation/msk468_genes_by_samples.csv  (0/1, for R)
Notes   : a technical benchmark of the METHOD, so all samples on the panel are
          used (not one per patient) -- it must match the numbers reported in
          Phase 2 Section 4c: 36,841 samples, 105,014 altered cells.
"""
import json

import numpy as np
import pandas as pd

from common import OUT, ROOT

PANEL, N_GENES = "MSK-IMPACT468", 45

clin = pd.read_parquet(ROOT / "data/processed/clinical_tidy.parquet", columns=["SAMPLE_ID", "SEQ_ASSAY_ID"])
samples = clin.loc[clin["SEQ_ASSAY_ID"] == PANEL, "SAMPLE_ID"]
alt = pd.read_parquet(ROOT / "data/processed/alterations_long.parquet", columns=["Sample_ID", "Hugo_Symbol"])
alt = alt[alt["Sample_ID"].isin(set(samples))].drop_duplicates()

genes = sorted(alt["Hugo_Symbol"].value_counts().index[:N_GENES])
M = (alt[alt["Hugo_Symbol"].isin(genes)].assign(v=True)
     .pivot(index="Sample_ID", columns="Hugo_Symbol", values="v")
     .reindex(index=samples, columns=genes).fillna(False).astype(bool))

np.save(OUT / "msk468_matrix.npy", M.to_numpy())
(OUT / "msk468_genes.json").write_text(json.dumps(genes))
pd.DataFrame(M.to_numpy().T.astype(int)).to_csv(OUT / "msk468_genes_by_samples.csv", index=False, header=False)
print(f"{M.shape[0]:,} samples x {M.shape[1]} genes, {int(M.to_numpy().sum()):,} altered cells")
