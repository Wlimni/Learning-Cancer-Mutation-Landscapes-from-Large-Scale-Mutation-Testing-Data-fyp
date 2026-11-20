"use strict";

// A native camera confirmation and a file selection both fire change.
const actions = document.querySelector(".phone-actions");
actions?.addEventListener("change", (event) => {
  const picker = event.target;
  if (picker.type !== "file" || !picker.files.length || actions.classList.contains("is-uploading")) return;
  actions.classList.add("is-uploading");
  document.getElementById("send-status").hidden = false;
  picker.form.requestSubmit();
});
window.addEventListener("pageshow", () => {
  if (!actions) return;
  actions.classList.remove("is-uploading");
  document.getElementById("send-status").hidden = true;
  for (const form of actions.querySelectorAll("form")) form.reset();
});
