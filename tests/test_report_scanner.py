"""Script Name: test_report_scanner.py.
Purpose: verify scanner extraction and paired upload behavior.
Inputs: synthetic report text and isolated Flask client. Outputs: test results.
Notes: no private report is bundled in the test suite.
"""

import io
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import pandas as pd

from scripts import report_lookup
from scripts.report_scanner.app import create_app
from scripts.report_scanner.extract import parse_findings


class ScannerTests(unittest.TestCase):
    def test_two_columns_and_directional_copy_number(self):
        report = "\n".join([
            "SUMMARY OF GENOMIC FINDINGS",
            " " * 70 + "GENOMIC FINDINGS WITH UNCERTAIN CLINICAL ACTIONABILITY",
            "NF1 F1247Ifs*18   10.10%      None" + " " * 39 + "MSH2 R406* 77.60%",
            " " * 70 + "CYP2D6 Copy Number 6.21 Fold",
            " " * 70 + "Loss",
            "Tumor Mutation Burden: TMB Intermediate (7.73 mut/Mb)",
            "No FDA-approved therapies are available in Glioma associated with Genomic Findings.",
        ])
        self.assertEqual(parse_findings(report), {
            "drivers": ["NF1 F1247Ifs*18"],
            "uncertain": ["MSH2 R406*", "CYP2D6 LOSS"],
            "cancer_type": "Glioma",
        })

    def test_search_stops_when_observed_results_exceed_the_budget(self):
        with TemporaryDirectory() as temp:
            root = Path(temp)
            raw = root / "data/raw"
            raw.mkdir(parents=True)
            pd.DataFrame({
                "IS_REPRESENTATIVE": [True], "CANCER_TYPE": ["Glioma"],
                "SAMPLE_ID": ["sample"], "SEQ_ASSAY_ID": ["panel"], "MUT_TESTED": [True],
            }).to_parquet(root / "clinical_tidy.parquet")
            genes = ["ABC", "BCD", "CDE"]
            pd.DataFrame({"SEQ_ASSAY_ID": ["panel"] * 3, "Hugo_Symbol": genes}).to_parquet(
                root / "panel_gene_coverage.parquet")
            pd.DataFrame({
                "Hugo_Symbol": genes, "HGVSp_Short": ["p.A1V"] * 3,
                "Tumor_Sample_Barcode": ["sample"] * 3, "Mutation_Status": ["Somatic"] * 3,
            }).to_csv(raw / "data_mutations_extended.txt", sep="\t", index=False)
            findings = [(g, "A1V") for g in genes]
            with patch.multiple(report_lookup, ROOT=root, PROCESSED_DIR=root):
                with self.assertRaisesRegex(ValueError, "Increase the minimum"):
                    report_lookup.scan_report("Glioma", findings, [], max_combinations=2)

    def test_pairing_and_single_use_upload(self):
        app = create_app("192.0.2.10", 5050)
        client = app.test_client()
        self.assertEqual(client.get("/").status_code, 200)
        status = client.get("/api/status").get_json()
        path = "/" + status["phone_url"].split("/", 3)[3]
        phone = {"REMOTE_ADDR": "192.0.2.12"}
        self.assertEqual(client.get("/", environ_overrides=phone).status_code, 403)
        self.assertEqual(client.get("/api/status", environ_overrides=phone).status_code, 403)
        self.assertEqual(client.get(path, environ_overrides=phone).status_code, 200)
        self.assertEqual(client.get("/capture/invalid", environ_overrides=phone).status_code, 404)
        url = path + "/upload"
        # The server checks the file format before OCR and permits another try after a bad upload.
        def upload(contents):
            return client.post(url, data={"report": (io.BytesIO(contents), "case.pdf")},
                               environ_overrides=phone)
        self.assertEqual(upload(b"not a PDF").status_code, 400)
        def extract_during_reset(_contents):
            self.assertEqual(client.post("/api/reset", json={}).status_code, 409)
            self.assertIn(b"Reading report", client.get(path, environ_overrides=phone).data)
            return {
                "drivers": ["TP53 R273H"], "uncertain": [], "cancer_type": "Glioma", "method": "PDF text"
            }
        with patch("scripts.report_scanner.app.extract_report", side_effect=extract_during_reset):
            self.assertEqual(upload(b"%PDF-1.7").status_code, 200)
            self.assertEqual(upload(b"%PDF-1.7").status_code, 409)
        self.assertEqual(client.get("/api/status").get_json()["draft"]["drivers"], ["TP53 R273H"])
        self.assertEqual(client.post("/api/reset", json={}).status_code, 200)
        self.assertEqual(client.get(path, environ_overrides=phone).status_code, 404)


if __name__ == "__main__":
    unittest.main()
