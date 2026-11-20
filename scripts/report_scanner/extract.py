"""Script Name: extract.py.
Purpose: extract draft report findings using local PDF text and OCR.
Inputs: PDF, JPEG or PNG bytes. Outputs: editable draft findings.
Notes: supports the demo layout and common gene/protein notation.
"""

import re
from shutil import copyfileobj, which
import subprocess
from pathlib import Path
from tempfile import TemporaryDirectory

GENE = r"[A-Z][A-Za-z0-9-]{1,19}"
PROTEIN = r"(?:p\.)?[A-Z]\d+(?:[A-Z*]|Ter)(?:fs(?:\*?\d+|Ter\d+)?)?"
TERT = r"(?:TERT(?:_promoter|\s+promoter)\s+)(?P<hotspot>C228T|C250T)"
FINDING = re.compile(
    rf"(?<![\w-])(?:{TERT}|(?P<gene>{GENE})[ \t]+(?P<detail>{PROTEIN}))(?![\w*])"
)
NUMBER = re.compile(rf"\b(?P<gene>{GENE})[ \t]+(?:copy[- ]number|CNV)\b", re.I)
DIRECTION = re.compile(r"\b(loss|deletion|gain|amplification)\b", re.I)
MAP = {"loss": "LOSS", "deletion": "LOSS", "gain": "GAIN", "amplification": "GAIN"}
UNCERTAIN = re.compile(r"GENOMIC FINDINGS WITH UNCERTAIN CLINICAL ACTIONABILITY", re.I)
CLINICAL_DECLARATION = re.compile(
    r"available in ([A-Za-z0-9 -]{2,80}?) associated with genomic findings", re.I
)


def parse_findings(text):
    """Read protein findings and copy-number directions without using VAFs as inputs."""
    drivers, uncertain = [], []
    summary = False
    boundary = None
    lines = text.splitlines()
    for index, line in enumerate(lines):
        if "SUMMARY OF GENOMIC FINDINGS" in line.upper():
            summary = True
        if section := UNCERTAIN.search(line):
            boundary = section.start()
        for match in FINDING.finditer(line):
            if match["hotspot"]:
                finding = f"TERT_promoter {match['hotspot']}"
            else:
                finding = f"{match['gene']} {match['detail'].removeprefix('p.')}"
            driver_row = summary and (boundary is None or (
                boundary >= 25 and match.start() < 24
                and re.search(r"\b(?:None|Yes)\b", line[match.end():boundary], re.I)
            ))
            (drivers if driver_row else uncertain).append(finding)
        for match in NUMBER.finditer(line):
            direction = DIRECTION.match(line[match.end():].lstrip(" :-"))
            if direction is None and index + 1 < len(lines):
                # In the demo PDF the direction can wrap beneath the gene.
                direction = next((m for m in DIRECTION.finditer(lines[index + 1])
                                  if abs(m.start() - match.start()) <= 12), None)
            if direction:
                uncertain.append(f"{match['gene']} {MAP[direction[1].lower()]}")
    cancer = CLINICAL_DECLARATION.search(" ".join(text.split()))
    return {
        "drivers": list(dict.fromkeys(drivers)),
        "uncertain": list(dict.fromkeys(uncertain)),
        "cancer_type": cancer[1] if cancer else "",
    }


def _run(command):
    try:
        return subprocess.run(command, capture_output=True, check=True, timeout=30).stdout
    except FileNotFoundError as exc:
        raise ValueError(f"Install {command[0]} on the laptop and try again.") from exc
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as exc:
        raise ValueError("The report could not be read. Try a clear image or an unlocked PDF.") from exc


def extract_report(stream):
    """Extract from a temporary local file and discard it when finished."""
    header = stream.read(64)
    stream.seek(0)
    is_pdf = header.startswith(b"%PDF-")
    is_heic = header[4:8] == b"ftyp" and any(
        brand in header[8:] for brand in (b"heic", b"heix", b"hevc", b"mif1", b"msf1")
    )
    if not (is_pdf or is_heic or header.startswith((b"\xff\xd8\xff", b"\x89PNG\r\n\x1a\n"))):
        raise ValueError("Choose a PDF, JPEG, PNG or HEIC image.")
    with TemporaryDirectory(prefix="report-scanner-") as temp:
        directory = Path(temp)
        extension = ".pdf" if is_pdf else ".heic" if is_heic else ""
        report = directory / ("report" + extension)
        with report.open("wb") as output:
            copyfileobj(stream, output)
        if is_pdf:
            text = _run(["pdftotext", "-f", "1", "-l", "4", "-layout", str(report), "-"]).decode("utf-8", "replace")
            draft = parse_findings(text)
            if draft["drivers"] or draft["uncertain"]:
                return {**draft, "method": "PDF text"}
            _run(["pdftoppm", "-f", "1", "-l", "4", "-scale-to", "2400",
                  "-png", str(report), str(directory / "page")])
            images = sorted(directory.glob("page-*.png"))
        else:
            if is_heic:
                if not which("sips"):
                    raise ValueError("Convert the HEIC image to JPEG or PNG and try again.")
                converted = directory / "photo.png"
                _run(["sips", "-s", "format", "png", str(report), "--out", str(converted)])
                report = converted
            images = [report]
        text = "\n".join(_run(["tesseract", str(image), "stdout", "--psm", "6",
                                "-c", "preserve_interword_spaces=1"]).decode("utf-8", "replace")
                         for image in images)
        return {**parse_findings(text), "method": "Photo OCR"}
