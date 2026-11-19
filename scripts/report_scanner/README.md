# Phone report scanner

The laptop shows a QR code for a one-time phone upload. It extracts a draft
locally and lets you correct the calls before using the same lookup as the
[Phase 3 notebook](../../notebooks/03_lookup_tool/03_phase3_lookup_tool.ipynb).

1. Set up the Python environment and local [GENIE data](../../data/README.md).
   Install the scanner dependencies from the repository root.

   ```sh
   .venv/bin/python -m pip install -r scripts/report_scanner/requirements.txt
   brew install poppler tesseract
   ```

   On systems without Homebrew, install Poppler (`pdftotext` and `pdftoppm`)
   and Tesseract using the system package manager.

2. Connect the laptop and phone to the same trusted Wi-Fi, then start the
   app from the repository root.

   ```sh
   .venv/bin/python -m scripts.report_scanner.app
   ```

   Open **http://127.0.0.1:5050** on the laptop. Scan the QR code with the
   phone camera. If the QR opens the wrong address, rerun with
   `--phone-host <laptop-Wi-Fi-IPv4>`. Use `--port` if 5050 is in use.

3. In the phone page, take a photo or choose a PDF, JPEG or PNG. On the
   laptop, check the extracted gene names, variants and categories. Correct
   anything the extraction missed. Enter the exact GENIE cancer type,
   then choose **Run GENIE lookup**.

The demo parser is tuned to the photographed layout. PDF text is used
first, with OCR for images and for PDFs whose text produces no findings.
The scanner reads up to four PDF pages and accepts uploads under 12 MB.
Keep MSI, TMB and VAF as report context. Enter an exact protein change
as `TP53 R273H`, a directional copy-number call as `CYP2D6 LOSS`, or
a bare gene to query the curated gene-level feature. Use `ATP1A1_fusion`
for a curated fusion feature if needed.

The QR link expires after 15 minutes and accepts one successful upload.
The review and results are available only at the laptop loopback address.
Uploads are temporary and session state is cleared when the laptop starts
a new session or exits. The phone connection uses local HTTP. Use a
synthetic demo report for presentations and keep real reports under the
ignored `data/reports/` folder. This app is intended for a local demo.

Results show observed combinations separately from gene-level Phase 2
evidence. The web view displays up to the first 200 combinations and gives
the total. If more than 10,000 combinations meet the count, the app asks
you to raise the minimum before running the lookup again.

Run scanner tests with `.venv/bin/python -m unittest discover -s tests -v`.
