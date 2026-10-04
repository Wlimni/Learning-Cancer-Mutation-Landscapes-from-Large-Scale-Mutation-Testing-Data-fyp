"""
Script Name: fetch_external.py
Purpose : Download the external reference data used by Phases 3-4 (not tracked in git):
          CIViC clinical evidence (nightly release) and the 10 TCGA oncogenic pathways
          (Sanchez-Vega et al. 2018, as distributed with cBioPortal's PathwayMapper).
Inputs  : internet access
Outputs : data/external/civic/nightly-ClinicalEvidenceSummaries.tsv
          data/external/civic/nightly-MolecularProfileSummaries.tsv
          data/external/pathways/sanchez_vega_2018_pathways.tsv  (pathway, gene)
Notes   : CIViC changes nightly -- record the download date when reporting results.
"""
import json
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
civic = ROOT / "data" / "external" / "civic"; civic.mkdir(parents=True, exist_ok=True)
for f in ["nightly-ClinicalEvidenceSummaries.tsv", "nightly-MolecularProfileSummaries.tsv"]:
    urllib.request.urlretrieve(f"https://civicdb.org/downloads/nightly/{f}", civic / f)
    print("downloaded", f)

url = "https://raw.githubusercontent.com/iVis-at-Bilkent/pathway-mapper/master/packages/pathway-mapper/src/data/pathways.json"
data = json.load(urllib.request.urlopen(url))
rows = [(pw, line.split("\t")[0])
        for pw in ["Cell Cycle", "HIPPO", "MYC", "NOTCH", "NRF2", "PI3K", "RTK-RAS", "TGF-Beta", "TP53", "WNT"]
        for line in data[pw] if len(line.split("\t")) > 2 and line.split("\t")[2] == "GENE"]
t = pd.DataFrame(rows, columns=["pathway", "gene"]).drop_duplicates()
t = t[~t["gene"].str.contains(" ")]                       # drop group nodes such as "Activin ligands"
out = ROOT / "data" / "external" / "pathways"; out.mkdir(parents=True, exist_ok=True)
t.to_csv(out / "sanchez_vega_2018_pathways.tsv", sep="\t", index=False)
print(f"pathways: {t['pathway'].nunique()} pathways, {t['gene'].nunique()} genes")
