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

    def test_manual_lookup_can_start_without_an_upload(self):
        client = create_app("192.0.2.10", 5050).test_client()
        pairing = client.get("/api/status").get_json()["pairing"]
        body = {"pairing": pairing, "cancer_type": "Glioma", "drivers": "TP53 R273H",
                "uncertain": "NF1", "min_observed": 1}
        with patch("scripts.report_scanner.app.Thread") as worker:
            self.assertEqual(client.post("/api/lookup", json=body).status_code, 202)
            self.assertEqual(worker.call_args.kwargs["args"][1:],
                             ("Glioma", [("TP53", "R273H")], ["NF1"], 1))
            self.assertEqual(client.post("/api/lookup", json=body).status_code, 409)
        self.assertEqual(client.get("/api/status").get_json()["state"], "running")

    def test_report_download_requires_a_completed_local_lookup_and_contains_all_rows(self):
        import pandas as pd

        client = create_app("192.0.2.10", 5050).test_client()
        self.assertEqual(client.get("/api/report").status_code, 409)
        pairing = client.get("/api/status").get_json()["pairing"]
        body = {"pairing": pairing, "cancer_type": "Glioma", "drivers": "TP53 R273H",
                "uncertain": "NF1", "min_observed": 1}
        tables = {"findings": pd.DataFrame(), "phase2_pairs": pd.DataFrame(),
                  "combinations": pd.DataFrame({"n_observed": range(1, 202)})}
        with patch("scripts.report_scanner.app.scan_report", return_value=tables), \
             patch("scripts.report_scanner.app.Thread", side_effect=lambda target, args, daemon:
                   type("Worker", (), {"start": lambda self: target(*args)})()):
            self.assertEqual(client.post("/api/lookup", json=body).status_code, 202)
        report = client.get("/api/report")
        self.assertEqual(report.status_code, 200)
        self.assertIn("attachment", report.headers["Content-Disposition"])
        self.assertEqual(report.get_json()["cancer_type"], "Glioma")
        self.assertEqual(len(report.get_json()["combinations"]), 201)
        self.assertEqual(len(client.get("/api/status").get_json()["results"]["combinations"]), 200)
        self.assertEqual(client.get("/api/report", environ_overrides={"REMOTE_ADDR": "192.0.2.12"}).status_code, 403)
        client.post("/api/reset", json={})
        self.assertEqual(client.get("/api/report").status_code, 409)

    def test_laptop_file_upload_checks_pairing_and_shares_review_state(self):
        client = create_app("192.0.2.10", 5050).test_client()
        status = client.get("/api/status").get_json()
        phone_url = "/" + status["phone_url"].split("/", 3)[3] + "/upload"

        def upload(token, data, remote="127.0.0.1"):
            return client.post("/api/upload", data={"report": (io.BytesIO(data), "report.pdf")},
                               headers={"X-Report-Pairing": token}, environ_overrides={"REMOTE_ADDR": remote})

        self.assertEqual(upload(status["pairing"], b"%PDF-1.7", remote="192.0.2.12").status_code, 403)
        self.assertEqual(upload("wrong", b"%PDF-1.7").status_code, 409)
        self.assertEqual(upload(status["pairing"], b"invalid").status_code, 400)
        with patch("scripts.report_scanner.app.extract_report", return_value={
            "drivers": ["TP53 R273H"], "uncertain": [], "cancer_type": "Glioma", "method": "PDF text"
        }):
            self.assertEqual(upload(status["pairing"], b"%PDF-1.7").status_code, 200)
            self.assertEqual(client.post(phone_url).status_code, 409)
        self.assertEqual(client.get("/api/status").get_json()["state"], "review")

    def test_phone_upload_passes_files_over_previous_limit_as_a_stream(self):
        client = create_app("192.0.2.10", 5050).test_client()
        phone_url = client.get("/api/status").get_json()["phone_url"]
        path = "/" + phone_url.split("/", 3)[3] + "/upload"
        large_file = io.BytesIO(b"%PDF-1.7" + b"\0" * (13 * 1024 * 1024))

        def extract(stream):
            self.assertEqual(stream.read(8), b"%PDF-1.7")
            self.assertGreater(stream.seek(0, io.SEEK_END), 12 * 1024 * 1024)
            return {"drivers": [], "uncertain": [], "cancer_type": "", "method": "PDF text"}

        with patch("scripts.report_scanner.app.extract_report", side_effect=extract):
            response = client.post(path, data={"report": (large_file, "large.pdf")})
        response.request.environ["wsgi.input"].close()
        self.assertEqual(response.status_code, 200)

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
