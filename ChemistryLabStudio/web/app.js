/* SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright © 2026 Jean-François Barthélémy (Cerema, UMR MCD)
 *
 * The page: forms that edit the model, the script the server writes from it,
 * the status of the compute process, and the results of a run.
 *
 * The model lives in the server. Every edit is sent there (debounced), and the
 * answer brings back the script and the problems; the forms are redrawn only
 * when the model changes from elsewhere (a library entry, an opened file), so
 * that typing is never interrupted. Server data reaches the page through
 * `textContent` only.
 */

(() => {
  const $ = (s) => document.querySelector(s);
  const S = {
    model: null, problems: [], path: null, stem: "study", script: "",
    catalog: null, library: [], results: null, status: { state: "starting" },
  };

  /* ── helpers ──────────────────────────────────────────────────── */

  function el(tag, attrs, ...children) {
    const e = document.createElement(tag);
    for (const [k, v] of Object.entries(attrs || {})) {
      if (v === undefined || v === null || v === false) continue;
      if (k === "class") e.className = v;
      else if (k === "text") e.textContent = v;
      else if (k.startsWith("on")) e.addEventListener(k.slice(2), v);
      else if (k in e && typeof v !== "string") e[k] = v;
      else e.setAttribute(k, v === true ? "" : v);
    }
    for (const c of children.flat()) if (c !== null && c !== undefined) e.append(c);
    return e;
  }

  function fmt(x, digits) {
    if (x === null || x === undefined || Number.isNaN(x)) return "—";
    const d = digits === undefined ? 3 : digits;
    if (x !== 0 && (Math.abs(x) < 1e-3 || Math.abs(x) >= 1e5)) return x.toExponential(d - 1);
    return x.toFixed(d);
  }

  const css = (name) => getComputedStyle(document.documentElement).getPropertyValue(name).trim();

  /* ── model and script ─────────────────────────────────────────── */

  let pushTimer = null;

  function edited() {
    clearTimeout(pushTimer);
    pushTimer = setTimeout(push, 250);
  }

  async function push() {
    const { status, data } = await Studio.post("/api/model", { model: S.model });
    if (status === 200) applyState(data, false);
    else toast(data.error || t("error"), true);
  }

  function applyState(d, redraw) {
    S.model = d.model;
    S.problems = d.problems;
    S.path = d.path;
    S.stem = d.stem;
    renderScript(d.script);
    renderPath();
    renderProblems();
    if (redraw) renderForms();
  }

  function renderScript(source) {
    S.script = source;
    const code = $("#code-text");
    code.removeAttribute("data-highlighted");
    code.className = "language-julia";
    code.textContent = source;
    if (window.hljs) hljs.highlightElement(code);
  }

  function renderPath() {
    const label = $("#path-label");
    label.textContent = S.path || t("untitled");
    label.title = S.path || t("untitled_hint");
  }

  function problemText(p) {
    const key = "problem_" + p.code;
    const s = t(key, p.values);
    return s === key ? p.text : s;
  }

  function renderProblems() {
    const list = $("#problems");
    list.replaceChildren(
      ...(S.problems.length ? S.problems.map((p) => el("li", { class: "problem", text: problemText(p) }))
        : [el("li", { class: "ok", text: t("no_problem") })])
    );
    updateRunButton();
  }

  /* ── forms ────────────────────────────────────────────────────── */

  function renderForms() {
    if (!S.model) return;
    $("#title").value = S.model.title;
    $("#w-b").value = S.model.w_b;
    $("#temperature").value = S.model.temperature_c;
    renderSystem();
    renderMaterials();
  }

  function renderSystem() {
    const m = S.model;
    const dbs = S.catalog ? S.catalog.databases.map((d) => d.name).filter((n) => n.endsWith(".json")) : [m.database];
    $("#database").replaceChildren(...dbs.map((n) => el("option", { value: n, text: n, selected: n === m.database })));

    const lists = S.catalog ? S.catalog.phase_lists : [{ name: m.phase_list, source: "" }];
    $("#phase-list").replaceChildren(...lists.map((l) => el("option", { value: l.name, text: l.name, selected: l.name === m.phase_list })));
    const current = lists.find((l) => l.name === m.phase_list);
    $("#phase-list-source").textContent = current && current.source ? t("source", { source: current.source }) : "";

    const models = S.catalog ? S.catalog.activity_models : ["cemdata18_KOH", "cemdata18_NaOH"];
    $("#activity-model").replaceChildren(
      ...models.map((key) =>
        el("label", { class: "choice" },
          el("input", {
            type: "radio", name: "activity-model", value: key, checked: key === m.activity_model,
            onchange: () => { m.activity_model = key; edited(); },
          }),
          el("span", {},
            el("b", { text: t("am_" + key) }),
            el("span", { class: "hint", text: t("am_" + key + "_text") })
          )
        )
      )
    );
  }

  function templateEntry(name) {
    return S.catalog ? S.catalog.templates.find((e) => e.name === name) : null;
  }

  function renderMaterials() {
    const m = S.model;
    const cards = m.materials.map((mat, i) => materialCard(mat, i));
    $("#materials").replaceChildren(...cards);
    renderFractionSum();
  }

  function renderFractionSum() {
    const total = S.model.materials.reduce((s, x) => s + (Number(x.fraction) || 0), 0);
    const note = $("#fraction-sum");
    note.textContent = t("fraction_sum", { total: fmt(total, 3) });
    note.classList.toggle("problem", Math.abs(total - 1) > 1e-9);
  }

  function materialCard(mat, i) {
    const m = S.model;
    const templates = S.catalog ? S.catalog.templates : [{ name: mat.template }];
    const select = el("select", {
      onchange: (e) => {
        mat.template = e.target.value;
        mat.extents = [];
        mat.material_extent = null;
        renderMaterials();
        edited();
      },
    }, templates.map((e) => el("option", {
      value: e.name, text: e.name, selected: e.name === mat.template,
      disabled: !!e.unavailable, title: e.unavailable || "",
    })));

    const remove = m.materials.length > 1
      ? el("button", { class: "small", title: t("remove_material"), text: "−",
          onclick: () => { m.materials.splice(i, 1); renderMaterials(); edited(); } })
      : null;

    const fraction = el("input", {
      type: "number", step: "0.01", min: "0", max: "1", value: mat.fraction,
      oninput: (e) => { mat.fraction = Number(e.target.value); renderFractionSum(); edited(); },
    });

    const entry = templateEntry(mat.template);
    const rows = [];
    const names = entry && entry.constituents ? entry.constituents.map((c) => c.name)
      : mat.extents.map((x) => x.constituent);
    for (const name of names) {
      const c = entry && entry.constituents ? entry.constituents.find((x) => x.name === name) : null;
      const set = mat.extents.find((x) => x.constituent === name);
      rows.push(el("tr", {},
        el("td", { text: name }),
        el("td", { class: "num", text: c ? fmt(100 * c.mass_fraction, 1) + " %" : "" }),
        el("td", {}, el("input", {
          type: "number", step: "0.05", min: "0", max: "1", placeholder: "1",
          value: set ? set.extent : "",
          oninput: (e) => {
            const v = e.target.value;
            mat.extents = mat.extents.filter((x) => x.constituent !== name);
            if (v !== "") mat.extents.push({ constituent: name, extent: Number(v) });
            edited();
          },
        }))
      ));
    }
    const table = el("table", { class: "data compact" },
      el("thead", {}, el("tr", {},
        el("th", { text: t("constituent") }), el("th", { text: t("mass_fraction") }), el("th", { text: t("degree") }))),
      el("tbody", {}, rows)
    );

    const source = entry && entry.source ? el("p", { class: "hint", text: t("source", { source: entry.source }) }) : null;
    return el("div", { class: "card" },
      el("header", {}, select, remove),
      source,
      el("div", { class: "field" }, el("label", { text: t("binder_fraction") }), fraction),
      rows.length ? table : el("p", { class: "hint", text: t("constituents_pending") }),
      el("p", { class: "hint", text: t("degree_hint") })
    );
  }

  function wireStaticForms() {
    $("#title").addEventListener("input", (e) => { S.model.title = e.target.value; edited(); });
    $("#w-b").addEventListener("input", (e) => { S.model.w_b = Number(e.target.value); edited(); });
    $("#temperature").addEventListener("input", (e) => { S.model.temperature_c = Number(e.target.value); edited(); });
    $("#phase-list").addEventListener("change", (e) => { S.model.phase_list = e.target.value; renderSystem(); edited(); });
    $("#material-add").addEventListener("click", () => {
      const first = S.catalog ? S.catalog.templates.find((e) => !e.unavailable) : null;
      S.model.materials.push({ template: first ? first.name : S.model.materials[0].template,
        fraction: 0, extents: [], material_extent: null });
      renderMaterials();
      edited();
    });
  }

  /* ── library ──────────────────────────────────────────────────── */

  async function loadLibrary() {
    const { status, data } = await Studio.get("/api/library");
    if (status !== 200) return;
    S.library = data;
    $("#library").replaceChildren(...data.map((e) =>
      el("div", { class: "card" },
        el("header", {}, el("b", { text: e.title }), e.designation ? el("span", { class: "tag", text: e.designation }) : null),
        el("p", { class: "hint", text: e.summary || "" }),
        el("p", { class: "hint", text: t("sources", { sources: (e.sources || []).join(", ") }) }),
        el("button", { class: "small primary", text: t("load"), onclick: () => loadEntry(e.id) })
      )
    ));
  }

  async function loadEntry(id) {
    const { status, data } = await Studio.post("/api/library", { id });
    if (status !== 200) return toast(data.error || t("error"), true);
    applyState(data, true);
    showTab("mix");
    toast(t("loaded"));
  }

  /* ── status of the compute process ───────────────────────────── */

  const BADGE = {
    preparing: "status_preparing", starting: "status_starting", ready: "status_ready",
    busy: "status_busy", failed: "status_failed", stopped: "status_stopped",
  };

  async function poll() {
    const { status, data } = await Studio.get("/api/status").catch(() => ({ status: 0, data: null }));
    if (status !== 200) {
      setBadge("bad", t("status_offline"));
      return;
    }
    S.status = data;
    const cls = data.state === "ready" ? "ok" : data.state === "failed" ? "bad" : "busy";
    setBadge(cls, t(BADGE[data.state] || data.state));
    $("#status").title = data.message || "";
    $("#stop").hidden = data.state !== "busy";
    $("#ready-card").hidden = !!S.results || data.state === "ready";
    $("#ready-message").textContent = data.message || t(BADGE[data.state] || data.state);
    $("#log").textContent = data.log.slice(-60).join("\n");
    $("#log").scrollTop = $("#log").scrollHeight;
    renderBanner(data);
    updateRunButton();
    if (!S.catalog && (data.state === "ready" || data.state === "busy")) loadCatalog();
  }

  function setBadge(cls, text) {
    const b = $("#status");
    b.className = "badge " + cls;
    b.textContent = text;
  }

  function renderBanner(data) {
    const banner = $("#banner");
    if (data.state !== "failed") { banner.hidden = true; return; }
    banner.hidden = false;
    banner.replaceChildren(
      el("b", { text: t("failed_title") }), " ",
      el("button", { class: "small", text: t("restart"), onclick: restart }),
      el("pre", { text: data.message })
    );
  }

  async function restart() {
    await Studio.post("/api/compute/restart");
    toast(t("restarting"));
  }

  let catalogLoading = false;
  async function loadCatalog() {
    if (catalogLoading) return;
    catalogLoading = true;
    const { status, data } = await Studio.get("/api/catalog");
    catalogLoading = false;
    if (status !== 200 || data.pending) return;
    if (data.error) { toast(t("catalog_failed") + " " + data.error, true); S.catalog = null; return; }
    S.catalog = data;
    renderForms();
  }

  function updateRunButton() {
    const ready = S.status.state === "ready";
    $("#run").disabled = !ready || S.problems.length > 0;
    $("#run").title = !ready ? t("run_wait") : S.problems.length ? t("run_problems") : "";
  }

  /* ── running ──────────────────────────────────────────────────── */

  async function run() {
    clearTimeout(pushTimer);
    await push();
    $("#run").disabled = true;
    setBadge("busy", t("status_busy"));
    const { status, data } = await Studio.post("/api/run");
    poll();
    if (status === 422) {
      toast(t("run_problems"), true);
      showTab("calculation");
      return;
    }
    if (status !== 200 || !data.ok) {
      $("#output").textContent = (data && (data.output || "")) + "\n" + ((data && data.error) || t("error"));
      showResult("output");
      toast(t("run_failed"), true);
      return;
    }
    S.results = data.results;
    $("#output").textContent = data.output || t("no_output");
    $("#ready-card").hidden = true;
    renderResults();
    toast(t("run_done", { seconds: fmt(data.elapsed, 2) }));
  }

  async function stop() {
    await Studio.post("/api/compute/stop");
    toast(t("stopping"));
  }

  /* ── results ──────────────────────────────────────────────────── */

  function plotLayout(extra) {
    const ink = css("--ink"), muted = css("--muted"), line = css("--line");
    const axis = { color: muted, gridcolor: line, zerolinecolor: line, automargin: true };
    return Object.assign({
      paper_bgcolor: "rgba(0,0,0,0)", plot_bgcolor: "rgba(0,0,0,0)",
      font: { color: ink, size: 12, family: getComputedStyle(document.body).fontFamily },
      margin: { l: 10, r: 20, t: 36, b: 40 },
      xaxis: Object.assign({}, axis), yaxis: Object.assign({}, axis, { autorange: "reversed" }),
      showlegend: false,
    }, extra || {});
  }

  /* A horizontal bar per label, every label shown: the plot grows with the
     number of bars rather than letting Plotly thin the axis. What has not
     reacted is drawn in the muted color, apart from what formed. */
  function barh(id, labels, values, title, xtitle, xtype) {
    const height = Math.max(220, 20 * labels.length + 100);
    const layout = plotLayout({ title: { text: title, font: { size: 13 } }, height });
    layout.xaxis.title = { text: xtitle };
    layout.yaxis.dtick = 1;
    if (xtype) layout.xaxis.type = xtype;
    const colors = labels.map((l) => (l.startsWith("unreacted ") ? css("--muted") : css("--accent")));
    const trace = { type: "bar", orientation: "h", x: values, y: labels, marker: { color: colors },
      hovertemplate: "%{y}: %{x:.4g}<extra></extra>" };
    document.getElementById(id).style.height = height + "px";
    Plotly.react(id, [trace], layout, { displaylogo: false, responsive: true });
  }

  function renderResults() {
    const r = S.results;
    if (!r) return;
    const binder = r.binder_mass_g ? fmt(r.binder_mass_g, 0) : "100";
    const pm = r.phase_masses;
    barh("plot-masses", pm.labels, pm.values, t("plot_masses"), t("plot_masses_axis", { binder }));
    const vf = r.volume_fractions;
    if (vf.labels) {
      $("#plot-volumes").hidden = false;
      $("#volumes-note").textContent = "";
      barh("plot-volumes", vf.labels, vf.values.map((v) => (v === null ? null : 100 * v)), t("plot_volumes"), t("plot_volumes_axis"));
    } else {
      $("#plot-volumes").hidden = true;
      $("#volumes-note").textContent = vf.reason;
    }

    const ps = r.pore_solution;
    $("#pore-summary").replaceChildren(
      summaryItem(t("ph"), fmt(ps.pH, 2)),
      summaryItem(t("ionic_strength"), withUnit(fmt(r.ionic_strength, 3), "mol/kg")),
      summaryItem(t("porosity"), withUnit(fmt(percent(r.porosity.total), 1), "%"), t("porosity_unknown")),
      summaryItem(t("bound_water"), withUnit(fmt(r.bound_water, 3), "g/g"))
    );
    // What is dissolved below a trace would squeeze the log axis to nothing.
    const mm = ps.molalities.values.map((v) => (v === null ? null : 1000 * v));
    const shown = ps.molalities.labels.map((k, i) => [k, mm[i]]).filter(([, v]) => v !== null && v >= TRACE);
    barh("plot-pore", shown.map((x) => x[0]), shown.map((x) => x[1]), t("plot_pore"), t("plot_pore_axis"), "log");
    fillTable("#pore-table", [t("element"), t("plot_pore_axis")],
      ps.molalities.labels.map((k, i) => [k, { text: mm[i] !== null && mm[i] < TRACE ? "< 1e-9" : fmt(mm[i], 3), cls: "" }]));

    fillTable("#si-table", [t("phase"), t("si")],
      r.saturation_indices.labels.map((k, i) => {
        const v = r.saturation_indices.values[i];
        return [k, { text: fmt(v, 3), cls: v === null ? "" : Math.abs(v) < 1e-6 ? "ok" : v > 0 ? "warn" : "" }];
      }));

    fillTable("#certificate-table", [t("cert_field"), t("cert_value")],
      Object.entries(r.certificate).map(([k, v]) => [t("cert_" + k), typeof v === "number" ? fmt(v, 3) : String(v)]));
  }

  // mmol per kg of water below which an element is reported as a trace
  const TRACE = 1e-9;

  // A missing value stays missing: in JavaScript `100 * null` is 0.
  const percent = (x) => (x === null || x === undefined ? null : 100 * x);
  const withUnit = (text, unit) => (text === "—" ? text : text + " " + unit);

  function summaryItem(label, value, missingHint) {
    return el("div", { class: "item", title: value === "—" && missingHint ? missingHint : null },
      el("span", { class: "hint", text: label }), el("b", { text: value }));
  }

  function fillTable(sel, head, rows) {
    $(sel).replaceChildren(
      el("thead", {}, el("tr", {}, head.map((h) => el("th", { text: h })))),
      el("tbody", {}, rows.map((row) => el("tr", {}, row.map((c) =>
        typeof c === "object" ? el("td", { class: "num " + c.cls, text: c.text }) : el("td", { text: c })))))
    );
  }

  /* ── files ────────────────────────────────────────────────────── */

  const dirOf = (p) => (p ? p.replace(/[\\/][^\\/]*$/, "") : "");

  async function save() {
    if (!S.path) return saveAs();
    await push();
    const { status, data } = await Studio.post("/api/save", {});
    if (status === 200) { applyState(data, false); toast(t("saved", { path: data.path })); }
    else toast(data.error || t("error"), true);
  }

  async function saveAs() {
    const path = await Picker.pick("save", dirOf(S.path), S.stem + ".jl");
    if (!path) return;
    await push();
    let { status, data } = await Studio.post("/api/save", { path });
    if (status === 409) {
      if (!confirm(t("overwrite", { path: data.path }))) return;
      ({ status, data } = await Studio.post("/api/save", { path, overwrite: true }));
    }
    if (status === 200) { applyState(data, false); toast(t("saved", { path: data.path })); }
    else toast(data.error || t("error"), true);
  }

  async function open() {
    const path = await Picker.pick("open", dirOf(S.path));
    if (!path) return;
    const { status, data } = await Studio.post("/api/open", { path });
    if (status !== 200) return toast(data.error || t("error"), true);
    applyState(data, true);
    toast(data.message || t("opened"), !!data.message);
  }

  /* ── tabs, theme, language ────────────────────────────────────── */

  function showTab(name) {
    document.querySelectorAll("#steps button").forEach((b) => b.classList.toggle("on", b.dataset.tab === name));
    document.querySelectorAll("#left .tab").forEach((d) => d.classList.toggle("on", d.id === "tab-" + name));
  }

  function showResult(name) {
    document.querySelectorAll("#result-tabs button").forEach((b) => b.classList.toggle("on", b.dataset.result === name));
    document.querySelectorAll("#results .result").forEach((d) => d.classList.toggle("on", d.id === "result-" + name));
    window.dispatchEvent(new Event("resize"));
  }

  const THEMES = ["auto", "light", "dark"];
  function applyTheme(theme) {
    if (theme === "auto") document.documentElement.removeAttribute("data-theme");
    else document.documentElement.setAttribute("data-theme", theme);
    $("#theme").textContent = theme === "dark" ? "☾" : theme === "light" ? "☀" : "◐";
    try { localStorage.setItem("clstudio.theme", theme); } catch (e) { /* not remembered */ }
    renderResults();
  }
  function savedTheme() {
    try { return localStorage.getItem("clstudio.theme") || "auto"; } catch (e) { return "auto"; }
  }

  async function toggleLanguage() {
    await Studio.setLanguage(Studio.language() === "en" ? "fr" : "en");
    $("#lang").textContent = Studio.language() === "en" ? "FR" : "EN";
    renderPath();
    renderProblems();
    renderForms();
    renderResults();
    if (S.library.length) loadLibrary();
  }

  /* ── start ────────────────────────────────────────────────────── */

  document.addEventListener("DOMContentLoaded", async () => {
    await Studio.setLanguage(Studio.savedLanguage());
    $("#lang").textContent = Studio.language() === "en" ? "FR" : "EN";
    applyTheme(savedTheme());
    document.querySelectorAll("#steps button").forEach((b) => b.addEventListener("click", () => showTab(b.dataset.tab)));
    document.querySelectorAll("#result-tabs button").forEach((b) => b.addEventListener("click", () => showResult(b.dataset.result)));
    $("#run").addEventListener("click", run);
    $("#stop").addEventListener("click", stop);
    $("#save").addEventListener("click", save);
    $("#saveas").addEventListener("click", saveAs);
    $("#open").addEventListener("click", open);
    $("#open2").addEventListener("click", open);
    $("#copy").addEventListener("click", async () => {
      try { await navigator.clipboard.writeText(S.script); toast(t("copied")); }
      catch (e) { toast(t("copy_failed"), true); }
    });
    $("#theme").addEventListener("click", () => applyTheme(THEMES[(THEMES.indexOf(savedTheme()) + 1) % THEMES.length]));
    $("#lang").addEventListener("click", toggleLanguage);
    wireStaticForms();

    const { status, data } = await Studio.get("/api/state");
    if (status !== 200) {
      setBadge("bad", t("status_denied"));
      toast(data.error || t("error"), true);
      return;
    }
    applyState(data, true);
    loadLibrary();
    poll();
    setInterval(poll, 1500);
  });
})();
