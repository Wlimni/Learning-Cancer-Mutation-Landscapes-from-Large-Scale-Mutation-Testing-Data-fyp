"""Script Name: app.py.
Purpose: pair a phone with a laptop, review extracted findings, and run GENIE lookup.
Inputs: local report uploads and reviewed findings. Outputs: aggregate results in browser.
Notes: local research demo for a trusted network. Files and session data are transient.
"""

import argparse
import hmac
import ipaddress
import re
import secrets
import socket
import subprocess
import sys
import time
from dataclasses import dataclass, field
from datetime import datetime, timezone
from functools import cache
from io import BytesIO
from threading import Lock, Thread

import pandas as pd
import qrcode
from flask import Flask, Response, abort, jsonify, render_template, request
from qrcode.image.svg import SvgPathImage

from scripts.report_lookup import PROCESSED_DIR, scan_report
from scripts.report_scanner.extract import GENE, PROTEIN, extract_report

PAIR_TIMEOUT = 15 * 60
COMBINATION_LIMIT = 10_000
FINDING_LINE = re.compile(
    rf"({GENE}(?:_promoter|_fusion)?)(?:\s+({PROTEIN}|AMP|DEL|GAIN|LOSS|C228T|C250T))?"
)


@dataclass
class Session:
    token: str = field(default_factory=lambda: secrets.token_urlsafe(24))
    expires: float = field(default_factory=lambda: time.monotonic() + PAIR_TIMEOUT)
    status: str = "waiting"
    draft: dict | None = None
    results: dict | None = None
    report: dict | None = None
    error: str = ""


def reviewed_findings(text):
    if not isinstance(text, str) or len(text) > 4000:
        raise ValueError("Enter one finding per line.")
    findings = []
    for line in text.splitlines():
        if not line.strip():
            continue
        match = FINDING_LINE.fullmatch(line.strip())
        if not match:
            raise ValueError(f"Check the finding format: {line[:50]}")
        gene, detail = match.groups()
        if gene.endswith("_fusion") and detail:
            raise ValueError("Enter a fusion as GENE_fusion.")
        if gene.endswith("_promoter") and (gene != "TERT_promoter" or detail not in ("C228T", "C250T")):
            raise ValueError("Enter a promoter hotspot as TERT_promoter C228T or C250T.")
        if detail in ("C228T", "C250T") and gene != "TERT_promoter":
            raise ValueError("Use TERT_promoter for a TERT promoter hotspot.")
        findings.append((gene, detail) if detail else gene)
    return findings


def create_app(phone_host, port):
    app = Flask(__name__)
    app.config.update(
        TRUSTED_HOSTS=["localhost", "127.0.0.1", phone_host],
        TEMPLATES_AUTO_RELOAD=True,
    )
    lock = Lock()
    session = Session()

    def desktop_only():
        if request.remote_addr != "127.0.0.1" or request.host.split(":")[0] not in {"localhost", "127.0.0.1"}:
            abort(403)

    def phone_session(token):
        if not token.isascii() or not hmac.compare_digest(token, session.token) or time.monotonic() >= session.expires:
            abort(404)
        return session

    @app.after_request
    def response_headers(response):
        response.headers["Cache-Control"] = "no-store"
        response.headers["Referrer-Policy"] = "no-referrer"
        response.headers["X-Content-Type-Options"] = "nosniff"
        response.headers["Content-Security-Policy"] = (
            "default-src 'self'; frame-ancestors 'none'; form-action 'self'; base-uri 'none'"
        )
        return response

    @app.get("/")
    def desktop():
        desktop_only()
        return render_template("desktop.html")

    def phone_url():
        return f"http://{phone_host}:{port}/capture/{session.token}"

    @app.get("/api/status")
    def status():
        desktop_only()
        with lock:
            return jsonify(
                state=session.status,
                pairing=session.token,
                phone_url=phone_url(),
                expired=time.monotonic() >= session.expires,
                draft=session.draft,
                results=session.results,
                error=session.error,
            )

    @cache
    def report_options():
        genes = set()
        for table in ("panel_gene_coverage", "cna_gene_panel_coverage", "sv_gene_panel_coverage"):
            try:
                genes.update(pd.read_parquet(PROCESSED_DIR / f"{table}.parquet",
                                             columns=["Hugo_Symbol"])["Hugo_Symbol"].dropna())
            except FileNotFoundError:
                continue
        try:
            cancer_types = pd.read_parquet(PROCESSED_DIR / "clinical_tidy.parquet",
                                           columns=["CANCER_TYPE"])["CANCER_TYPE"].dropna().unique()
        except FileNotFoundError:
            cancer_types = []
        return {"genes": sorted(genes), "cancer_types": sorted(cancer_types)}

    @app.get("/api/options")
    def options():
        desktop_only()
        return jsonify(report_options())

    @app.get("/qr.svg")
    def qr():
        desktop_only()
        with lock:
            url = phone_url()
        qr_image = qrcode.make(url, image_factory=SvgPathImage, border=4)
        stream = BytesIO()
        qr_image.save(stream)
        return Response(stream.getvalue(), mimetype="image/svg+xml")

    @app.get("/capture/<token>")
    def capture(token):
        with lock:
            current = phone_session(token)
            if current.status == "extracting":
                return render_template("phone.html", extracting=True, token=token)
            sent = current.status != "waiting"
        return render_template("phone.html", sent=sent, token=token)

    def process_upload(token, *, from_phone):
        with lock:
            if from_phone:
                current = phone_session(token)
            elif token != session.token:
                return "Start a new session and try again.", 409
            else:
                current = session
            if current.status != "waiting":
                return "A report is already in this session. Start a new session to send another.", 409
            current.status = "extracting"
        try:
            report = request.files.get("report")
            if not report:
                raise ValueError("Choose a PDF or image.")
            draft = extract_report(report.stream)
        except ValueError as exc:
            return str(exc), 400
        else:
            with lock:
                current.draft = draft
                current.status = "review"
            return "", 200
        finally:
            with lock:
                if current.status == "extracting":
                    current.status = "waiting"

    @app.post("/capture/<token>/upload")
    def upload(token):
        error, code = process_upload(token, from_phone=True)
        if code == 409:
            abort(409, description=error)
        return render_template("phone.html", sent=not error, token=token, error=error), code

    @app.post("/api/upload")
    def upload_from_laptop():
        desktop_only()
        error, code = process_upload(request.headers.get("X-Report-Pairing"), from_phone=False)
        return (jsonify(error=error) if error else jsonify(ok=True)), code

    @app.get("/api/report")
    def download_report():
        desktop_only()
        with lock:
            report = session.report if session.status == "complete" else None
        if report is None:
            return jsonify(error="Run GENIE lookup before saving a report."), 409
        response = jsonify(report)
        response.headers["Content-Disposition"] = 'attachment; filename="genie-lookup-report.json"'
        return response

    @app.post("/api/reset")
    def reset():
        nonlocal session
        desktop_only()
        request.get_json()
        with lock:
            if session.status in {"extracting", "running"}:
                return jsonify(error="Wait for the current operation to finish."), 409
            session = Session()
        return jsonify(ok=True)

    def lookup(current, cancer, drivers, uncertain, minimum):
        try:
            tables = scan_report(cancer, drivers, uncertain, min_observed=minimum,
                                 max_combinations=COMBINATION_LIMIT)
            combinations = tables["combinations"]
            report = {
                "dataset": "AACR Project GENIE",
                "cancer_type": cancer,
                "min_observed": minimum,
                "generated_at": datetime.now(timezone.utc).isoformat(),
                "findings": tables["findings"].to_dict("records"),
                "combinations": combinations.to_dict("records"),
                "total_combinations": len(combinations),
                "phase2_pairs": tables["phase2_pairs"].to_dict("records"),
            }
            preview = {**report, "combinations": report["combinations"][:200]}
            with lock:
                current.report, current.results, current.status = report, preview, "complete"
        except (FileNotFoundError, ValueError) as exc:
            error = ("Local GENIE data is missing. Follow the project data setup before running lookup."
                     if isinstance(exc, FileNotFoundError) else str(exc))
            with lock:
                current.error, current.status = error, "review" if current.draft else "waiting"
        except Exception:
            app.logger.exception("GENIE lookup failed")
            with lock:
                current.error = "Lookup failed. Check the laptop terminal for details."
                current.status = "review" if current.draft else "waiting"

    @app.post("/api/lookup")
    def run_lookup():
        desktop_only()
        body = request.get_json()
        if not isinstance(body, dict):
            abort(400)
        try:
            cancer = body.get("cancer_type")
            if not isinstance(cancer, str) or not cancer.strip() or len(cancer) > 100:
                raise ValueError("Enter the GENIE cancer type, for example Glioma.")
            drivers = reviewed_findings(body.get("drivers"))
            uncertain = reviewed_findings(body.get("uncertain"))
            if not 2 <= len(set(drivers + uncertain)) <= 50:
                raise ValueError("Enter between 2 and 50 distinct findings.")
            minimum = body.get("min_observed")
            if type(minimum) is not int or not 1 <= minimum <= 10000:
                raise ValueError("Minimum observed must be an integer from 1 to 10,000.")
        except ValueError as exc:
            return jsonify(error=str(exc)), 400
        with lock:
            current = session
            if body.get("pairing") != current.token or current.status not in {"waiting", "review", "complete"}:
                return jsonify(error="Wait for the upload or lookup to finish, or start a new session."), 409
            current.status, current.error, current.results, current.report = "running", "", None, None
        Thread(target=lookup, args=(current, cancer.strip(), drivers, uncertain, minimum), daemon=True).start()
        return jsonify(ok=True), 202

    return app


def laptop_ipv4():
    # Prefer Wi-Fi on macOS. The default route can point through a VPN.
    if sys.platform == "darwin":
        ports = subprocess.run(["networksetup", "-listallhardwareports"], capture_output=True,
                               text=True, timeout=5).stdout
        wifi = re.search(r"Hardware Port: Wi-Fi\s+Device: (\S+)", ports)
        if wifi:
            address = subprocess.run(["ipconfig", "getifaddr", wifi[1]], capture_output=True,
                                     text=True, timeout=5).stdout.strip()
            if address:
                return address
    with socket.socket(socket.AF_INET, socket.SOCK_DGRAM) as route:
        route.connect(("192.0.2.1", 9))
        return route.getsockname()[0]


def main():
    parser = argparse.ArgumentParser(description="Run the local genomic report scanner.")
    parser.add_argument("--phone-host", help="Laptop Wi-Fi IPv4 address for the QR link")
    parser.add_argument("--port", type=int, default=5050)
    args = parser.parse_args()
    try:
        host = str(ipaddress.IPv4Address(args.phone_host or laptop_ipv4()))
    except (ipaddress.AddressValueError, OSError, subprocess.TimeoutExpired) as exc:
        parser.error(f"Specify the laptop Wi-Fi IPv4 address with --phone-host: {exc}")
    print(f"Laptop review: http://127.0.0.1:{args.port}", flush=True)
    create_app(host, args.port).run(host="0.0.0.0", port=args.port, debug=False, threaded=True)


if __name__ == "__main__":
    main()
