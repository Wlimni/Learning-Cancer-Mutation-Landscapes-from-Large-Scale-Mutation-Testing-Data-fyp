"use strict";

const byId = (id) => document.getElementById(id);
const form = byId("review-form");
const findingRows = byId("finding-rows");
const geneSymbols = new Map();
let edited = false;
let draft = null;
const uploadForm = byId("upload-form");
const filePicker = byId("report-file");
let selectedReport = null;
let pairing = "";
let draftLoaded = false;
let resultsLoaded = false;
let reportSnapshot = null;
let actionPending = false;
let requestRevision = 0;
let errorTimer = null;
const exportPdf = byId("export-pdf");

function message(id, text, timeout = 0) {
  if (id === "error") window.clearTimeout(errorTimer);
  const element = byId(id);
  if (element.textContent !== text) element.textContent = text;
  element.hidden = !text;
  if (id === "error" && text && timeout) {
    errorTimer = window.setTimeout(() => message(id, ""), timeout);
  }
}

async function post(path, data) {
  const response = await fetch(path, {
    method: "POST",
    headers: {"Content-Type": "application/json"},
    body: JSON.stringify(data),
  });
  if (!response.ok) {
    const body = await response.json().catch(() => ({}));
    throw new Error(body.error || "The request failed. Please try again.");
  }
}

function addTable(container, title, headings, rows, empty) {
  const section = document.createElement("section");
  container.append(section);
  const heading = document.createElement("h3");
  heading.textContent = title;
  section.append(heading);
  if (!rows.length) {
    const text = document.createElement("p");
    text.textContent = empty;
    section.append(text);
    return;
  }
  const wrap = document.createElement("div");
  wrap.className = "table-wrap";
  const table = document.createElement("table");
  table.createTHead().insertRow().append(...headings.map((value) => {
    const cell = document.createElement("th");
    cell.scope = "col";
    cell.textContent = value;
    return cell;
  }));
  const body = table.createTBody();
  for (const values of rows) {
    const row = body.insertRow();
    values.forEach((value) => {
      const cell = row.insertCell();
      cell.textContent = value;
    });
  }
  wrap.append(table);
  section.append(wrap);
}

function setExportAvailable(available) {
  if (exportPdf) exportPdf.hidden = exportPdf.disabled = !available;
}

function setRunCompleted(completed) {
  byId("run").dataset.settled = String(completed);
  byId("run-label").textContent = completed ? "Re-Run GENIE lookup" : "Run GENIE lookup";
}

exportPdf?.addEventListener("click", () => {
  const previousTitle = document.title;
  document.title = "GENIE Lookup Report";
  try {
    window.print();
  } finally {
    document.title = previousTitle;
  }
});

const reportStage = byId("ep-report");
const firstPage = reportStage.querySelector(".page");

// Keep the paper at A4 dimensions so the preview and print use the same line breaks.
new ResizeObserver(([entry]) => {
  reportStage.style.setProperty("--preview-scale", Math.min(1, entry.contentRect.width / firstPage.offsetWidth));
}).observe(reportStage);

function clearReportPages() {
  for (const page of reportStage.querySelectorAll(".continuation-page")) page.remove();
  firstPage.querySelector(".page-number").textContent = "";
}

function continuationPage() {
  const page = document.createElement("article");
  page.className = "page continuation-page";
  const body = document.createElement("div");
  body.className = "report-body result-tables";
  const footer = firstPage.querySelector(".report-foot").cloneNode(true);
  footer.querySelector("[id]").removeAttribute("id");
  page.append(body, footer);
  reportStage.append(page);
  return body;
}

function reportOverflows(container) {
  const body = container.closest(".report-body");
  return body.scrollHeight > body.clientHeight + 1;
}

function paginateResults(panel) {
  clearReportPages();
  const sections = [...panel.children];
  panel.replaceChildren();
  let current = panel;
  for (const section of sections) {
    current.append(section);
    if (!reportOverflows(current)) continue;
    const table = section.querySelector("table");
    if (!table) {
      current = continuationPage();
      current.append(section);
      continue;
    }

    let tbody = table.tBodies[0];
    const rows = [...tbody.rows];
    tbody.replaceChildren();
    let fragment = section;
    for (const [index, row] of rows.entries()) {
      tbody.append(row);
      if (!reportOverflows(current)) continue;

      // Leave complete rows on this page and repeat the heading and columns on the next.
      row.remove();
      if (!tbody.rows.length) fragment.remove();
      current = continuationPage();
      fragment = section.cloneNode(false);
      const heading = section.querySelector("h3").cloneNode(true);
      if (index) heading.textContent += " (continued)";
      const wrap = table.parentElement.cloneNode(false);
      const nextTable = table.cloneNode(false);
      nextTable.append(table.tHead.cloneNode(true), document.createElement("tbody"));
      wrap.append(nextTable);
      fragment.append(heading, wrap);
      current.append(fragment);
      tbody = nextTable.tBodies[0];
      tbody.append(row);
    }
  }
  const pages = [...reportStage.querySelectorAll(".page")];
  pages.forEach((page, i) => {
    page.querySelector(".page-number").textContent = `Page ${i + 1} of ${pages.length}`;
  });
}

function showResults(results) {
  reportSnapshot = results;
  updateSummary();
  const panel = byId("result-tables");
  panel.replaceChildren();
  const date = new Date(results.generated_at);
  byId("foot-date").textContent = Number.isNaN(date.getTime()) ? "" : `Generated ${date.toLocaleDateString()}`;
  addTable(panel, "Observed combinations", ["Size", "Findings", "Patients observed", "Jointly tested", "Frequency"],
    results.combinations.map((r) => [r.size, r.combination.join(" + "), r.n_observed, r.n_tested, r.pct.toFixed(2) + "%"]),
    "No combinations met the minimum observed count.");
  if (results.total_combinations > results.combinations.length) {
    const text = document.createElement("p");
    text.textContent = `Showing ${results.combinations.length} of ${results.total_combinations} combinations. Increase the minimum count to narrow the results.`;
    panel.append(text);
  }
  addTable(panel, "Phase 2 gene-pair evidence", ["Gene pair", "Direction", "q-value", "Replicated"],
    results.phase2_pairs.map((r) => [r.Gene_A + " + " + r.Gene_B, r.direction, r.q_value.toPrecision(3), r.replicated ? "Yes" : "No"]),
    "No Phase 2 candidates matched the submitted genes.");
  addTable(panel, "Individual findings", ["Finding", "Report group", "Patients observed", "Tested"],
    results.findings.map((r) => [r.finding, r.report_group, r.n_observed, r.n_tested]),
    "No findings were available.");
  byId("empty-state").hidden = true;
  byId("results").hidden = false;
  const previewScale = reportStage.style.getPropertyValue("--preview-scale");
  reportStage.style.setProperty("--preview-scale", "1");
  try {
    paginateResults(panel);
  } finally {
    if (previewScale) reportStage.style.setProperty("--preview-scale", previewScale);
    else reportStage.style.removeProperty("--preview-scale");
  }
}

let rowCounter = 0;
function addFinding(group = "driver", gene = "", alteration = "") {
  const row = byId("finding-row").content.firstElementChild.cloneNode(true);
  row.querySelector(".gene").value = gene;
  row.querySelector(".alteration").value = alteration;
  const name = `group-${rowCounter++}`; // one radio set per row
  for (const radio of row.querySelectorAll(".grp")) {
    radio.name = name;
    radio.checked = radio.value === group;
  }
  findingRows.append(row);
  updateSummary();
  return row;
}

const rowGroup = (row) => row.querySelector(".grp:checked")?.value ?? "driver";

// Live preview of what will be searched, shown on the report page.
function updateSummary() {
  const cancer = byId("cancer").value.trim();
  const minimum = byId("minimum").value || "1";
  const inputFindings = [...findingRows.children].flatMap((row) => {
    const gene = row.querySelector(".gene").value.trim();
    if (!gene) return [];
    const detail = row.querySelector(".alteration").value.trim();
    return [{finding: [gene, detail].filter(Boolean).join(" "), report_group: rowGroup(row)}];
  });
  const distinct = new Set(inputFindings.map(({finding}) => finding.replace(/ p\./, " ").toUpperCase()));
  byId("run").dataset.ready = String(Boolean(cancer) && distinct.size >= 2 && byId("minimum").validity.valid);

  const reportCancer = reportSnapshot?.cancer_type ?? cancer;
  const findings = reportSnapshot?.findings ?? inputFindings;
  const sumCancer = byId("sum-cancer");
  sumCancer.textContent = reportCancer || "Not set";
  sumCancer.classList.toggle("empty", !reportCancer);
  byId("sum-min").textContent = reportSnapshot?.min_observed ?? minimum;
  const chips = byId("sum-chips");
  chips.replaceChildren();
  for (const {finding, report_group} of findings) {
    const chip = document.createElement("span");
    chip.className = `chip chip-${report_group.includes("driver") ? "driver" : "uncertain"}`;
    const findingLabel = document.createElement("span");
    findingLabel.textContent = finding;
    const dot = document.createElement("i");
    dot.setAttribute("aria-hidden", "true");
    chip.append(dot, findingLabel);
    chips.append(chip);
  }
  byId("sum-count").textContent = findings.length;
}

function clearEditor() {
  form.reset();
  findingRows.replaceChildren();
  addFinding();
  addFinding("uncertain");
  edited = false;
  draft = null;
  byId("import-draft").hidden = true;
  message("error", "");
  message("review-note", "");
  updateSummary();
}

function reviewMessage(text, warning = false) {
  byId("review-note").classList.toggle("note-warning", warning);
  message("review-note", text);
}

function importDraft() {
  if (!draft) return;
  findingRows.replaceChildren();
  byId("cancer").value = draft.cancer_type || byId("cancer").value;
  for (const [group, findings] of [["driver", draft.drivers], ["uncertain", draft.uncertain]]) {
    for (const finding of findings) {
      const [gene, ...detail] = finding.trim().split(/\s+/);
      addFinding(group, gene, detail.join(" "));
    }
  }
  if (!findingRows.children.length) addFinding();
  updateSummary();
  edited = false;
  byId("import-draft").hidden = true;
  const hasFindings = draft.drivers.length || draft.uncertain.length;
  reviewMessage(hasFindings
    ? `${draft.method} extracted a draft. Check each finding and report group before running the lookup.`
    : "No findings were recognized. Enter them here or start a new session to scan a clearer report.", !hasFindings);
}

async function loadOptions() {
  try {
    const response = await fetch("/api/options");
    if (!response.ok) return;
    const {genes, cancer_types} = await response.json();
    for (const gene of genes) {
      geneSymbols.set(gene.toUpperCase(), gene);
    }
    for (const cancer of cancer_types) byId("cancer-types").append(new Option(cancer, cancer));
    if (document.activeElement.matches(".gene")) suggestGenes(document.activeElement);
  } catch {
    // The editor still accepts typed gene symbols when suggestions are unavailable.
  }
}

function suggestGenes(input) {
  const prefix = input.value.trim().toUpperCase();
  const options = byId("gene-options");
  options.replaceChildren();
  if (prefix.length < 2) return;
  for (const [key, gene] of geneSymbols) {
    if (key.startsWith(prefix)) options.append(new Option(gene, gene));
    if (options.children.length >= 12) break;
  }
}

findingRows.addEventListener("focusin", (event) => {
  if (event.target.matches(".gene")) suggestGenes(event.target);
});

function collectFindings() {
  const findings = {driver: [], uncertain: []};
  for (const row of findingRows.children) {
    const geneInput = row.querySelector(".gene");
    const entered = geneInput.value.trim();
    const detail = row.querySelector(".alteration").value.trim();
    if (!entered && !detail) continue;
    if (!entered) {
      message("error", "Enter the gene for this alteration.");
      if (byId("pane-findings").classList.contains("ep-on")) {
        geneInput.setCustomValidity("Enter the gene for this alteration.");
        geneInput.reportValidity();
      }
      return null;
    }
    const gene = geneSymbols.get(entered.toUpperCase()) ?? entered;
    findings[rowGroup(row)].push([gene, detail].filter(Boolean).join(" "));
  }
  return findings;
}

form.addEventListener("input", (event) => {
  edited = true;
  updateSummary();
  if (event.target.matches(".gene")) suggestGenes(event.target);
  event.target.closest(".finding-row")?.querySelector(".gene").setCustomValidity("");
});

byId("add-finding").addEventListener("click", () => {
  edited = true;
  addFinding().querySelector(".gene").focus();
});

findingRows.addEventListener("click", (event) => {
  const button = event.target.closest(".remove-row");
  if (!button) return;
  edited = true;
  const row = button.closest(".finding-row");
  const next = row.nextElementSibling || row.previousElementSibling;
  row.remove();
  updateSummary();
  (next?.querySelector(".gene") || byId("add-finding")).focus();
});

byId("import-draft").addEventListener("click", () => {
  importDraft();
  byId("cancer").focus();
});

// Panel tabs: 1 · Import, 2 · Findings
const tabs = [...document.querySelectorAll(".ep-tab")];
function showTab(name) {
  for (const tab of tabs) {
    const on = tab.dataset.t === name;
    tab.classList.toggle("ep-on", on);
    tab.setAttribute("aria-selected", String(on));
    byId(tab.getAttribute("aria-controls")).classList.toggle("ep-on", on);
  }
}
for (const tab of tabs) tab.addEventListener("click", () => showTab(tab.dataset.t));

// Drag the panel edge to resize it, like the report editor.
const panel = byId("ep-panel");
byId("ep-drag").addEventListener("pointerdown", (event) => {
  event.preventDefault();
  const handle = event.currentTarget;
  handle.classList.add("ep-dragging");
  handle.setPointerCapture(event.pointerId);
  const move = (e) => panel.style.setProperty("--pw", `${Math.round(e.clientX)}px`);
  const stop = () => {
    handle.classList.remove("ep-dragging");
    handle.removeEventListener("pointermove", move);
    handle.removeEventListener("pointerup", stop);
  };
  handle.addEventListener("pointermove", move);
  handle.addEventListener("pointerup", stop);
});

const loadingPanel = document.createElement("div");
loadingPanel.id = "lookup-loading";
loadingPanel.className = "lookup-loading";
loadingPanel.hidden = true;
loadingPanel.setAttribute("aria-hidden", "true");
const loadingCard = document.createElement("div");
loadingCard.className = "lookup-loading-card";
const loadingLabel = document.createElement("strong");
loadingCard.append(loadingLabel);
loadingPanel.append(loadingCard);
document.body.append(loadingPanel);

fetch("/static/dna-trail.svg")
  .then((response) => response.ok ? response.text() : "")
  .then((source) => {
    const svg = new DOMParser().parseFromString(source, "image/svg+xml").documentElement;
    if (svg.localName !== "svg") return;
    svg.setAttribute("aria-hidden", "true");
    for (const [selector, className] of [[".ep-hd", "header-dna"], [".page-head", "report-banner-dna"]]) {
      const backdrop = svg.cloneNode(true);
      backdrop.classList.add(className);
      // Each decorative copy has its own SVG IDs and no animated highlights.
      for (const element of backdrop.querySelectorAll("[id]")) element.id = className + "-" + element.id;
      for (const element of backdrop.querySelectorAll("use")) {
        if (element.classList.contains("dna-trail")) element.remove();
        else element.setAttribute("href", "#" + className + "-" + element.getAttribute("href").slice(1));
      }
      document.querySelector(selector).prepend(backdrop);
    }
    loadingCard.prepend(svg);
  })
  .catch(() => { /* The loading text remains available. */ });

function setLoading(label = "") {
  if (label) message("status", label);
  loadingLabel.textContent = label;
  loadingPanel.hidden = !label;
  byId("results").setAttribute("aria-busy", String(Boolean(label)));
}

async function refresh() {
  const revision = requestRevision;
  try {
    const response = await fetch("/api/status", {cache: "no-store"});
    if (!response.ok) throw new Error("The laptop server is unavailable.");
    const data = await response.json();
    if (revision !== requestRevision) return;
    if (pairing !== data.pairing) {
      pairing = data.pairing;
      byId("qr").src = "/qr.svg?" + encodeURIComponent(pairing);
      byId("phone-link").href = data.phone_url;
      reportSnapshot = null;
      clearReportPages();
      byId("foot-date").textContent = "";
      clearEditor();
      showTab("import");
      byId("empty-state").hidden = false;
      uploadForm.reset();
      selectReport(null);
      byId("results").hidden = true;
      draftLoaded = resultsLoaded = false;
      setRunCompleted(false);
      setExportAvailable(false);
    }
    const status = {
      waiting: "Import a report or enter findings manually.",
      extracting: "Extracting report findings…",
      review: "Ready for review",
      running: "Searching GENIE…",
      complete: "Lookup complete",
    };
    if (!actionPending) {
      message("status", status[data.state]);
      setLoading(["running", "extracting"].includes(data.state) ? status[data.state] : "");
    }
    const acceptingReport = ["waiting", "extracting"].includes(data.state);
    byId("phone-intake").hidden = data.state !== "waiting";
    byId("phone-link").hidden = data.expired;
    message("phone-help", data.expired ? "Phone link expired. Select New for a new QR code." : "Scan on the same Wi-Fi.");
    uploadForm.hidden = !acceptingReport;
    byId("import-done").hidden = acceptingReport ||
      !(data.draft?.drivers.length || data.draft?.uncertain.length);
    byId("upload-fields").disabled = data.state !== "waiting" || actionPending;
    byId("reset").disabled = ["extracting", "running"].includes(data.state) || actionPending;
    const busy = ["extracting", "running"].includes(data.state) || actionPending;
    byId("review-fields").disabled = busy;
    byId("run").disabled = busy;
    byId("import-draft").disabled = busy;
    setExportAvailable(data.state === "complete" && Boolean(data.results) && !busy);
    if (data.draft && !draftLoaded) {
      draft = data.draft;
      if (edited) {
        if (draft.drivers.length || draft.uncertain.length) {
          byId("import-draft").hidden = false;
          reviewMessage("A scan is ready. Your entries are still here. Replace them with the scan if you want to review it.");
        } else {
          reviewMessage("No findings were recognized in the scan. You can keep entering them here.", true);
        }
      } else {
        importDraft();
      }
      draftLoaded = true;
      showTab("findings");
    }
    if (data.results && !resultsLoaded) {
      showResults(data.results);
      setRunCompleted(true);
      reviewMessage("Lookup complete. You can export the report as a PDF.");
      resultsLoaded = true;
      byId("empty-state").hidden = true;
    }
    if (data.error) message("error", data.error);
  } catch (err) {
    if (revision === requestRevision && !actionPending) {
      setLoading();
      message("status", err.message);
    }
  } finally {
    window.setTimeout(refresh, 900);
  }
}

function selectReport(file) {
  selectedReport = file;
  message("selected-report", file ? file.name : "PDF, JPEG, PNG or HEIC.");
  byId("scan-report").disabled = !file;
  message("upload-error", "");
}

filePicker.addEventListener("change", () => selectReport(filePicker.files[0] ?? null));

uploadForm.addEventListener("submit", async (event) => {
  event.preventDefault();
  if (!selectedReport || actionPending) return;
  message("upload-error", "");
  const data = new FormData();
  data.append("report", selectedReport);
  actionPending = true;
  requestRevision++;
  setLoading("Extracting report findings…");
  byId("upload-fields").disabled = true;
  try {
    const response = await fetch("/api/upload", {
      method: "POST",
      headers: {"X-Report-Pairing": pairing},
      body: data,
    });
    const result = await response.json().catch(() => ({}));
    if (!response.ok) throw new Error(result.error || "Upload failed. Try choosing the file again.");
    uploadForm.hidden = true;
    message("status", "Ready for review");
  } catch (err) {
    message("upload-error", err.message);
    byId("upload-fields").disabled = false;
  } finally {
    setLoading();
    requestRevision++;
    actionPending = false;
  }
});

byId("reset").addEventListener("click", async () => {
  actionPending = true;
  byId("reset").disabled = true;
  setExportAvailable(false);
  message("error", "");
  try {
    await post("/api/reset", {});
  } catch (err) {
    message("error", err.message);
  } finally {
    actionPending = false;
  }
});

form.addEventListener("submit", async (event) => {
  event.preventDefault();
  if (!pairing || actionPending) return;
  const findings = collectFindings();
  if (!findings) return;
  if (!findings.driver.length && !findings.uncertain.length) {
    message("error", "No information yet. Import a report or enter findings first.", 4000);
    return;
  }
  if (!form.checkValidity()) {
    if (byId("pane-findings").classList.contains("ep-on")) form.reportValidity();
    message("error", "Check the cancer type and minimum patient count.");
    return;
  }
  actionPending = true;
  requestRevision++;
  setLoading("Searching GENIE…");
  setExportAvailable(false);
  message("status", "Searching GENIE…");
  byId("run").disabled = true;
  message("error", "");
  try {
    await post("/api/lookup", {
      pairing,
      cancer_type: byId("cancer").value,
      drivers: findings.driver.join("\n"),
      uncertain: findings.uncertain.join("\n"),
      min_observed: Number(byId("minimum").value),
    });
    reportSnapshot = null;
    clearReportPages();
    byId("foot-date").textContent = "";
    byId("results").hidden = true;
    byId("empty-state").hidden = false;
    resultsLoaded = false;
    updateSummary();
  } catch (err) {
    setLoading();
    message("error", err.message);
  } finally {
    requestRevision++;
    actionPending = false;
  }
});

clearEditor();
loadOptions();
refresh();
