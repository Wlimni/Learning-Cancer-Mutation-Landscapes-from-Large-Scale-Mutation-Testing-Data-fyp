# Phone report scanner

Scan a report directly on the laptop or send it from a phone using the QR code.
The right-hand editor accepts manually entered findings or a scanned draft.
Both paths use the same lookup as the
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

   Open **http://127.0.0.1:5050** on the laptop. If you started the app
   with `--port 5057`, open **http://127.0.0.1:5057** instead. Open the
   Flask URL, since the HTML file and VS Code Live Server do not provide
   the lookup API.

   Enter the cancer type and findings in the **Findings** tab. Gene
   suggestions come from the local GENIE panel coverage tables. To import
   instead, open **Import**, choose a PDF, JPEG, PNG or HEIC (including a saved
   screenshot), and click **Scan report**. To use a phone, scan the QR code
   under **From your phone** in the same tab.
   If the QR opens the wrong address, rerun with
   `--phone-host <laptop-Wi-Fi-IPv4>`. Use `--port` if 5050 is in use.

3. On the phone page, take a photo or choose a PDF or image. The phone uploads it after you confirm.
   With either upload method, check the extracted gene names, variants and
   categories on the laptop. Correct
   anything the extraction missed. Enter the exact GENIE cancer type,
   then choose **Run GENIE lookup**.

The demo parser is tuned to the photographed layout. PDF text is used
first, with OCR for images and for PDFs whose text produces no findings.
The scanner reads up to four PDF pages. HEIC photos can be read on macOS.
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


After a lookup, the report preview is laid out on A4 pages. Tables that
continue onto another page repeat their column headings. **Export PDF**
prints those same pages, including the summary and the displayed results.
Select **Save as PDF** in the print dialog. **Start new session** clears the current findings and results.
