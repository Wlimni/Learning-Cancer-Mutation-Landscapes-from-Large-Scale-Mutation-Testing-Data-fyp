"use strict";

const byId = (id) => document.getElementById(id);
const form = byId("review-form");
let pairing = "";
let draftLoaded = false;
let resultsLoaded = false;
let actionPending = false;

function message(id, text) {
  const element = byId(id);
  element.textContent = text;
  element.hidden = !text;
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
  const heading = document.createElement("h3");
  heading.textContent = title;
  container.append(heading);
  if (!rows.length) {
    const text = document.createElement("p");
    text.textContent = empty;
    container.append(text);
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
    for (const value of values) row.insertCell().textContent = value;
  }
  wrap.append(table);
  container.append(wrap);
}

function showResults(results) {
  const panel = byId("result-tables");
  panel.replaceChildren();
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
  byId("results").hidden = false;
}

async function refresh() {
  try {
    const response = await fetch("/api/status", {cache: "no-store"});
    if (!response.ok) throw new Error("The laptop server is unavailable.");
    const data = await response.json();
    if (pairing !== data.pairing) {
      pairing = data.pairing;
      byId("qr").src = "/qr.svg?" + encodeURIComponent(pairing);
      byId("phone-link").href = data.phone_url;
      form.reset();
      form.hidden = true;
      byId("results").hidden = true;
      draftLoaded = resultsLoaded = false;
      message("review-note", "Your draft will appear here after upload.");
    }
    const status = {
      waiting: data.expired ? "Pairing expired. Start a new session." : "Waiting for the phone",
      extracting: "Extracting report findings…",
      review: "Ready for review",
      running: "Searching GENIE…",
      complete: "Lookup complete",
    };
    message("status", status[data.state]);
    byId("phone-link").hidden = data.state !== "waiting" || data.expired;
    byId("reset").disabled = ["extracting", "running"].includes(data.state) || actionPending;
    byId("run").disabled = data.state === "running" || actionPending;
    if (data.draft && !draftLoaded) {
      byId("cancer").value = data.draft.cancer_type;
      byId("drivers").value = data.draft.drivers.join("\n");
      byId("uncertain").value = data.draft.uncertain.join("\n");
      const noFindings = !data.draft.drivers.length && !data.draft.uncertain.length;
      message("review-note", noFindings
        ? "No findings were recognized. Add them by hand or start a new session with a clearer report."
        : `${data.draft.method} extracted a draft. Check every call and category against the report.`);
      form.hidden = false;
      draftLoaded = true;
    }
    if (data.results && !resultsLoaded) {
      showResults(data.results);
      resultsLoaded = true;
    }
    if (data.error) message("error", data.error);
  } catch (err) {
    message("status", err.message);
  }
  window.setTimeout(refresh, 900);
}

byId("reset").addEventListener("click", async () => {
  actionPending = true;
  byId("reset").disabled = true;
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
  actionPending = true;
  byId("run").disabled = true;
  message("error", "");
  try {
    await post("/api/lookup", {
      pairing,
      cancer_type: byId("cancer").value,
      drivers: byId("drivers").value,
      uncertain: byId("uncertain").value,
      min_observed: Number(byId("minimum").value),
    });
    message("status", "Searching GENIE…");
    byId("results").hidden = true;
    resultsLoaded = false;
  } catch (err) {
    message("error", err.message);
  } finally {
    actionPending = false;
  }
});

refresh();
