/* SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright © 2026 Jean-François Barthélémy (Cerema, UMR MCD)
 *
 * The page: forms that edit the model, the script the server writes from it,
 * the findings of the diagnosis, the status of the compute process, and the
 * results of a run.
 *
 * The model lives in the server. Every edit changes the page's copy and is
 * sent there (debounced); the answer brings back the script and the problems.
 * The forms are redrawn only when the model changes from elsewhere (a library
 * entry, an opened file) or when an edit changes what a form offers (another
 * template, another mode), so that typing is never interrupted. Server data
 * reach the page through `textContent` only.
 */

(() => {
  const $ = (s) => document.querySelector(s);
  const el = Forms.make;
  const S = {
    model: null, problems: [], path: null, stem: "study", script: "",
    catalog: null, schema: null, library: [], results: null, status: { state: "starting" },
    species: null, speciesFor: null, diagnostics: null, filter: { text: "", state: "crystal", elements: [] },
  };

  /* ── helpers ──────────────────────────────────────────────────── */

  function fmt(x, digits) {
    if (x === null || x === undefined || Number.isNaN(x)) return "—";
    const d = digits === undefined ? 3 : digits;
    if (x !== 0 && (Math.abs(x) < 1e-3 || Math.abs(x) >= 1e5)) return x.toExponential(d - 1);
    return x.toFixed(d);
  }

  const css = (name) => getComputedStyle(document.documentElement).getPropertyValue(name).trim();
  // `replaceChildren` would write a null as the text "null": the absent parts are dropped.
  const fill = (node, ...kids) => node.replaceChildren(...kids.flat().filter((c) => c !== null && c !== undefined && c !== false));
  const unique = (xs) => [...new Set(xs)];

  /* ── model and script ─────────────────────────────────────────── */

  let pushTimer = null;
  let pushing = Promise.resolve();

  function edited() {
    clearTimeout(pushTimer);
    pushTimer = setTimeout(() => { pushing = push(); }, 250);
    scheduleDiagnosis();
  }

  async function push() {
    const { status, data } = await Studio.post("/api/model", { model: S.model });
    if (status === 200) applyState(data, false);
    else toast(data.error || t("error"), true);
  }

  async function flush() {
    clearTimeout(pushTimer);
    await pushing;
    await push();
  }

  /* The page keeps its own copy while the user edits: the server's answer to
     an edit is the same model, and replacing the copy would detach the forms
     from what they edit. */
  function applyState(d, redraw) {
    if (redraw) S.model = d.model;
    S.problems = d.problems;
    S.path = d.path;
    S.stem = d.stem;
    renderScript(d.script);
    renderPath();
    renderProblems();
    if (redraw) { renderForms(); scheduleDiagnosis(); }
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

  // The parts of the model a problem names, in the page's language.
  const WHERE = { "the system": "where_system", "the activity model": "where_activity", "the mix": "where_mix",
    "the calculation": "where_calculation", "the constraint": "where_constraint" };

  function problemText(p) {
    const key = "problem_" + p.code;
    const v = Object.assign({}, p.values);
    if (v.field) v.field = Forms.fieldLabel(v.field);
    if (v.where) v.where = v.where.split(", ").map((w) => (WHERE[w] ? t(WHERE[w]) : w)).join(", ");
    for (const k of ["min", "max"]) if (v[k] === null || v[k] === undefined) v[k] = "—";
    const s = t(key, v);
    return s === key ? p.text : s;
  }

  function renderProblems() {
    fill($("#problems"),
      ...(S.problems.length ? S.problems.map((p) => el("li", { class: "problem", text: problemText(p) }))
        : [el("li", { class: "ok", text: t("no_problem") })]));
    updateRunButton();
  }

  /* ── diagnosis ────────────────────────────────────────────────── */

  let diagnosisTimer = null;
  let diagnosing = Promise.resolve();

  function scheduleDiagnosis() {
    clearTimeout(diagnosisTimer);
    diagnosisTimer = setTimeout(() => { diagnosing = diagnose(); }, 1500);
  }

  async function diagnose() {
    if (S.status.state !== "ready" || S.problems.length) {
      S.diagnostics = null;
      renderDiagnostics();
      return;
    }
    await pushing;
    const { status, data } = await Studio.post("/api/diagnose");
    if (status !== 200) return;
    if (data.ok) S.diagnostics = data.diagnostics;
    else if (data.error && !/already running/.test(data.error)) S.diagnostics = [{ level: "blocking", code: "failed", message: data.error, fix: "", subjects: [] }];
    renderDiagnostics();
  }

  const LEVEL_ICON = { blocking: "⛔", warning: "⚠", info: "ℹ" };

  function renderDiagnostics() {
    const items = S.diagnostics === null
      ? [el("li", { class: "hint", text: t(S.problems.length ? "diagnosis_after_problems" : "diagnosis_pending") })]
      : S.diagnostics.length
        ? S.diagnostics.map((d) => el("li", { class: "diag " + d.level },
            el("span", { class: "icon", text: LEVEL_ICON[d.level] || "•" }),
            el("span", {}, el("span", { text: diagnosticText(d) }), d.fix ? el("span", { class: "hint fix", text: " " + d.fix } ) : null)))
        : [el("li", { class: "ok", text: t("diagnosis_clean") })];
    for (const id of ["#diagnostics", "#diagnostics-system"]) fill($(id), ...items.map((x) => x.cloneNode(true)));
  }

  // The messages are ChemistryLab's own, in English; a translation keyed by
  // the code replaces them when it exists.
  function diagnosticText(d) {
    const key = "diag_" + d.code;
    const s = t(key, { subjects: d.subjects.join(", ") });
    return s === key ? d.message : s;
  }

  /* ── forms ────────────────────────────────────────────────────── */

  function renderForms() {
    if (!S.model || !S.schema) return;
    $("#title").value = S.model.title;
    renderSystem();
    renderMaterials();
    renderMix();
    renderCalculation();
  }

  /* System: the database, the system and its species, the activity model. */
  function renderSystem() {
    const m = S.model;
    fill($("#database-form"), Forms.kinded("database", m.database, (v) => {
      const changed = v.kind !== m.database.kind || v.name !== m.database.name;
      m.database = v;
      edited();
      if (changed) { S.species = null; renderSpecies(); }
    }));
    renderDatabaseStatus();
    fill($("#system-form"), Forms.kinded("system", m.system, (v) => {
      const redraw = v.kind !== m.system.kind || v.phase_list !== m.system.phase_list || v.gel !== m.system.gel;
      m.system = v;
      edited();
      if (redraw) renderPhases();
    }, { ctx: { skip: ["add", "remove", "phases", "solid_solutions"] } }));
    fill($("#activity-form"), Forms.kinded("activity", m.activity, (v) => { m.activity = v; edited(); }));
    renderPhases();
    renderSpecies();
  }

  function renderDatabaseStatus() {
    const dbs = S.catalog ? S.catalog.databases : [];
    fill($("#database-status"), ...dbs.map((d) =>
      el("li", {}, el("span", { class: "tag " + (d.status === "local" || d.status === "cached" ? "ok" : ""), text: t("db_" + d.status) }),
        " ", el("span", { text: d.name }))));
  }

  function currentList() {
    const s = S.model.system;
    return s.kind === "phase_list" && S.catalog ? S.catalog.phase_lists.find((l) => l.name === s.phase_list) : null;
  }

  /* The phases of the system as chips: where each comes from, and a button to
     leave it out or take it back. */
  function renderPhases() {
    const s = S.model.system;
    const chips = [];
    const chip = (name, cls, why, action, actionText) => el("span", { class: "chip " + cls, title: why },
      el("span", { text: name }), action ? el("button", { class: "x", text: actionText || "×", title: why, onclick: action }) : null);
    if (s.kind === "phase_list") {
      const list = currentList();
      const base = list ? unique(list.reactants.concat(list.products)) : [];
      const sss = list ? list.solid_solutions.map((x) => (x === "CSHQ" ? s.gel : x)) : [];
      for (const p of base) {
        if (s.remove.includes(p)) chips.push(chip(p, "removed", t("chip_removed"), () => { s.remove = s.remove.filter((x) => x !== p); changedPhases(); }, "↺"));
        else chips.push(chip(p, list.reactants.includes(p) ? "reactant" : "", t(list.reactants.includes(p) ? "chip_reactant" : "chip_product"), () => { s.remove.push(p); changedPhases(); }));
      }
      for (const p of sss) chips.push(chip(p, "ss", t("chip_solid_solution")));
      for (const p of s.add) chips.push(chip(p, "added", t("chip_added"), () => { s.add = s.add.filter((x) => x !== p); changedPhases(); }));
      if (!list) chips.push(el("span", { class: "hint", text: t("constituents_pending") }));
    } else if (s.kind === "custom") {
      for (const p of s.phases) chips.push(chip(p, "", t("chip_custom"), () => { s.phases = s.phases.filter((x) => x !== p); changedPhases(); }));
      for (const p of s.solid_solutions) chips.push(chip(p, "ss", t("chip_solid_solution"), () => { s.solid_solutions = s.solid_solutions.filter((x) => x !== p); changedPhases(); }));
    }
    fill($("#phases"), ...chips);
    renderSolidSolutionPicker();
  }

  function changedPhases() { renderPhases(); renderSpeciesList(); edited(); }

  function systemPhases() {
    const s = S.model.system;
    if (s.kind === "custom") return new Set(s.phases);
    const list = currentList();
    const base = list ? list.reactants.concat(list.products).filter((p) => !s.remove.includes(p)) : [];
    return new Set(base.concat(s.add));
  }

  function addPhase(p) {
    const s = S.model.system;
    if (s.kind === "custom") { if (!s.phases.includes(p)) s.phases.push(p); }
    else if (s.remove.includes(p)) s.remove = s.remove.filter((x) => x !== p);
    else if (!s.add.includes(p)) s.add.push(p);
    changedPhases();
  }

  function renderSolidSolutionPicker() {
    const s = S.model.system;
    const box = $("#ss-picker");
    if (s.kind !== "custom" || !S.catalog) { fill(box); return; }
    const names = S.catalog.solid_solutions.map((x) => x.name).filter((n) => !s.solid_solutions.includes(n));
    const select = el("select", {}, names.map((n) => el("option", { value: n, text: n })));
    fill(box, el("div", { class: "row" }, select,
      el("button", { class: "small", text: t("add_solid_solution"), onclick: () => { s.solid_solutions.push(select.value); changedPhases(); } })));
  }

  /* The species of the database: a search, a filter by state and by the
     elements a species may be made of, and a button to add one. */
  async function renderSpecies() {
    const db = S.model.database;
    const box = $("#species-browser");
    if (db.kind !== "builtin") { fill(box, el("p", { class: "hint", text: t("species_builtin_only") })); return; }
    if (!S.species || S.speciesFor !== db.name) {
      fill(box, el("p", { class: "hint", text: t("species_loading") }));
      if (S.status.state !== "ready") return;
      const name = db.name;
      const { status, data } = await Studio.get("/api/species?database=" + encodeURIComponent(name));
      if (status !== 200 || S.model.database.name !== name) return;
      S.species = data;
      S.speciesFor = name;
    }
    const elements = unique(S.species.flatMap((x) => x.elements)).sort(byPeriod);
    fill(box,
      el("div", { class: "row" },
        el("input", { type: "search", placeholder: t("species_search"), value: S.filter.text, spellcheck: "false",
          oninput: (e) => { S.filter.text = e.target.value; renderSpeciesList(); } }),
        el("select", { onchange: (e) => { S.filter.state = e.target.value; renderSpeciesList(); } },
          ["crystal", "aqueous", "gas", "all"].map((x) => el("option", { value: x, text: t("state_" + x), selected: x === S.filter.state })))),
      el("div", { class: "elements", title: t("elements_hint") }, elements.map((e) =>
        el("button", { class: "element" + (S.filter.elements.includes(e) ? " on" : ""), text: e, onclick: (ev) => {
          S.filter.elements = S.filter.elements.includes(e) ? S.filter.elements.filter((x) => x !== e) : S.filter.elements.concat([e]);
          ev.target.classList.toggle("on");
          renderSpeciesList();
        } }))),
      el("div", { id: "species-list", class: "species-list" }));
    renderSpeciesList();
  }

  // Elements by period, then group, as a periodic table reads.
  const ORDER = "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am".split(" ");
  const byPeriod = (a, b) => (ORDER.indexOf(a) + 1 || 999) - (ORDER.indexOf(b) + 1 || 999);

  function renderSpeciesList() {
    const box = document.getElementById("species-list");
    if (!box || !S.species) return;
    const f = S.filter;
    const text = f.text.trim().toLowerCase();
    const inSystem = systemPhases();
    const hits = S.species.slice().sort((a, b) => a.symbol.localeCompare(b.symbol)).filter((x) =>
      (f.state === "all" || x.state === f.state) &&
      (!f.elements.length || x.elements.every((e) => f.elements.includes(e))) &&
      (!text || x.symbol.toLowerCase().includes(text) || x.formula.toLowerCase().includes(text) || (x.cement || "").toLowerCase().includes(text)));
    const T = S.model.mix.temperature + 273.15;
    fill(box,
      el("p", { class: "hint", text: t("species_count", { n: hits.length }) }),
      el("table", { class: "data compact" }, el("tbody", {}, hits.slice(0, 150).map((x) => {
        const outside = x.T_min !== null && x.T_max !== null && (T < x.T_min || T > x.T_max);
        const there = inSystem.has(x.symbol);
        return el("tr", { title: t("species_card", { formula: x.formula, M: fmt(x.molar_mass, 2), lo: fmt(x.T_min - 273.15, 0), hi: fmt(x.T_max - 273.15, 0) }) },
          el("td", {}, el("b", { text: x.symbol })),
          el("td", { text: x.formula }),
          // a cement notation with an improper fraction (H₃//₂) reads worse than the formula
          el("td", { class: "hint", text: x.cement && !x.cement.includes("//") ? x.cement : "" }),
          el("td", { class: "warn", text: outside ? "⚠" : "", title: outside ? t("outside_range") : null }),
          el("td", {}, x.state === "crystal"
            ? el("button", { class: "small", text: there ? "✓" : "+", disabled: there, title: t(there ? "in_system" : "add_to_system"), onclick: () => addPhase(x.symbol) })
            : null));
      }))),
      hits.length > 150 ? el("p", { class: "hint", text: t("species_more", { n: hits.length - 150 }) }) : null);
  }

  /* Materials: a card per material, its kind, its fields, a law per
     constituent. */
  const CLINKER = ["C3S", "C2S", "C3A", "C4AF"];

  function templateEntry(name) {
    return S.catalog ? S.catalog.templates.find((e) => e.name === name) : null;
  }

  function constituentsOf(mat) {
    if (mat.kind === "template") {
      const e = templateEntry(mat.template);
      return e && e.constituents ? e.constituents.map((c) => c.name) : [];
    }
    if (mat.kind === "oxides") return mat.route === "bogue" ? CLINKER.concat([mat.sulfate]) : [mat.name];
    return mat.minerals.map((r) => r.formula).filter((x) => x);
  }

  function materialContext(mat) {
    return {
      keys: () => constituentsOf(mat),
      note: (f, key) => {
        if (mat.kind !== "template") return "";
        const e = templateEntry(mat.template);
        const c = e && e.constituents ? e.constituents.find((x) => x.name === key) : null;
        return c ? fmt(100 * c.mass_fraction, 1) + " %" : "";
      },
      noneLabel: () => t("reacts_fully"),
      units: () => ["%"],
    };
  }

  function renderMaterials() {
    const m = S.model;
    fill($("#materials"), ...m.materials.map((mat, i) => materialCard(mat, i)));
    renderFractionSum();
  }

  function materialCard(mat, i) {
    const m = S.model;
    const card = el("div", { class: "card material" });
    function draw() {
      const before = { kind: mat.kind, template: mat.template, route: mat.route, sulfate: mat.sulfate, name: mat.name, rows: JSON.stringify(mat.minerals) };
      const entry = mat.kind === "template" ? templateEntry(mat.template) : null;
      fill(card,
        el("header", {}, el("b", { text: t("material_n", { n: i + 1 }) }),
          m.materials.length > 1 ? el("button", { class: "small", title: t("remove_material"), text: "−",
            onclick: () => { m.materials.splice(i, 1); renderMaterials(); edited(); } }) : null),
        Forms.kinded("material", mat, (v) => {
          m.materials[i] = v;
          mat = v;
          renderFractionSum();
          edited();
          // Another template, route or list of minerals offers other constituents.
          const now = { kind: v.kind, template: v.template, route: v.route, sulfate: v.sulfate, name: v.name, rows: JSON.stringify(v.minerals) };
          if (now.kind !== before.kind || now.template !== before.template) { v.extents = []; v.material_extent = null; draw(); }
          else if (now.route !== before.route || now.sulfate !== before.sulfate || (now.name !== before.name && v.route === "oxides") || now.rows !== before.rows) draw();
        }, { ctx: materialContext(mat) }),
        entry && entry.source ? el("p", { class: "hint refs" }, t("source_label"), " ", reference(entry.source)) : null);
    }
    draw();
    return card;
  }

  function renderFractionSum() {
    const total = S.model.materials.reduce((s, x) => s + (Number(x.fraction) || 0), 0);
    const note = $("#fraction-sum");
    note.textContent = t("fraction_sum", { total: fmt(total, 3) });
    note.classList.toggle("problem", Math.abs(total - 1) > 1e-9);
  }

  function addMaterial() {
    const first = S.catalog ? S.catalog.templates.find((e) => !e.unavailable && e.kind === "scm") : null;
    const mat = Forms.defaults("material", "template");
    if (first) mat.template = first.name;
    mat.fraction = 0;
    S.model.materials.push(mat);
    renderMaterials();
    edited();
  }

  function renderMix() {
    fill($("#mix-form"), Forms.fieldsOf("mix", S.model.mix, (v) => { S.model.mix = v; edited(); renderSpeciesList(); }));
  }

  /* Calculation: the three modes, the calculations of the chosen one, and the
     constraint on the equilibrium. */
  const modeOf = (kind) => { const m = Forms.memberOf("calculation", kind); return m ? m.mode : "equilibrium"; };

  function renderCalculation() {
    const m = S.model;
    const mode = modeOf(m.calculation.kind);
    document.querySelectorAll("#modes button").forEach((b) => {
      b.classList.toggle("on", b.dataset.mode === mode);
      b.disabled = !Forms.membersOf("calculation").some((x) => x.mode === b.dataset.mode);
    });
    fill($("#calculation-form"), Forms.kinded("calculation", m.calculation, (v) => {
      m.calculation = v;
      edited();
    }, { only: (x) => x.mode === mode, hideKind: mode === "equilibrium", ctx: { units: () => ["mmol/L", "mol/L"] } }));
    fill($("#constraint-form"), Forms.kinded("constraint", m.constraint, (v) => { m.constraint = v; edited(); }));
  }

  function setMode(mode) {
    if (modeOf(S.model.calculation.kind) === mode) return;
    const first = Forms.membersOf("calculation").find((x) => x.mode === mode);
    if (!first) return;
    S.model.calculation = Forms.defaults("calculation", first.kind);
    renderCalculation();
    edited();
  }

  /* ── library ──────────────────────────────────────────────────── */

  async function loadLibrary() {
    const { status, data } = await Studio.get("/api/library");
    if (status !== 200) return;
    S.library = data;
    renderLibrary();
  }

  /* A published work, by its key: its citation linked to its DOI once the
     catalog of ChemistryLab is read, its key until then. */
  function reference(key) {
    const r = S.catalog && S.catalog.references ? S.catalog.references[key] : null;
    if (!r) return el("span", { text: key });
    return r.doi ? el("a", { href: "https://doi.org/" + r.doi, target: "_blank", rel: "noopener", text: r.citation })
      : el("span", { text: r.citation });
  }

  const references = (keys) => keys.flatMap((k, i) => [i ? el("br") : null, reference(k)]);

  function renderLibrary() {
    fill($("#library"), ...S.library.map((e) =>
      el("div", { class: "card" },
        el("header", {}, el("b", { text: e.title }), e.designation ? el("span", { class: "tag", text: e.designation }) : null,
          e.process ? el("span", { class: "tag", text: Forms.kindLabel(e.process) }) : null),
        el("p", { class: "hint", text: e.summary || "" }),
        el("p", { class: "hint refs" }, ...references(e.sources || [])),
        el("button", { class: "small primary", text: t("load"), onclick: () => loadEntry(e.id) }))));
  }

  async function loadEntry(id) {
    const { status, data } = await Studio.post("/api/library", { id });
    if (status !== 200) return toast(data.error || t("error"), true);
    applyState(data, true);
    showTab("materials");
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
    const was = S.status.state;
    S.status = data;
    const cls = data.state === "ready" ? "ok" : data.state === "failed" ? "bad" : "busy";
    setBadge(cls, t(BADGE[data.state] || data.state));
    $("#status").title = data.message || "";
    $("#stop").hidden = data.state !== "busy";
    $("#ready-card").hidden = !!S.results || data.state === "ready" || S.running;
    $("#ready-message").textContent = data.message || t(BADGE[data.state] || data.state);
    $("#log").textContent = data.log.slice(-60).join("\n");
    $("#log").scrollTop = $("#log").scrollHeight;
    renderBanner(data);
    updateRunButton();
    if (!S.catalog && (data.state === "ready" || data.state === "busy")) loadCatalog();
    if (data.state === "ready" && was !== "ready") {
      if (!S.species) renderSpecies();
      if (S.diagnostics === null) scheduleDiagnosis();
    }
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
    fill(banner,
      el("b", { text: t("failed_title") }), " ",
      el("button", { class: "small", text: t("restart"), onclick: restart }),
      el("pre", { text: data.message }));
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
    Forms.setCatalog(data);
    renderForms();
    renderLibrary();
  }

  function updateRunButton() {
    const ready = S.status.state === "ready";
    $("#run").disabled = !ready || S.problems.length > 0;
    $("#run").title = !ready ? t("run_wait") : S.problems.length ? t("run_problems") : "";
  }

  /* ── running ──────────────────────────────────────────────────── */

  async function run() {
    clearTimeout(diagnosisTimer);
    await flush();
    await diagnosing;
    $("#run").disabled = true;
    setBadge("busy", t("status_busy"));
    $("#stop").hidden = false;
    S.running = true;
    $("#ready-card").hidden = true;
    const { status, data } = await Studio.post("/api/run");
    S.running = false;
    poll();
    if (status === 422) {
      toast(t("run_problems"), true);
      showTab("calculation");
      return;
    }
    if (status !== 200 || !data.ok) {
      $("#output").textContent = ((data && data.output) || "") + "\n" + ((data && data.error) || t("error"));
      showResult("output");
      toast(data && data.stopped ? t("run_stopped") : t("run_failed"), true);
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
      xaxis: Object.assign({}, axis), yaxis: Object.assign({}, axis),
      showlegend: false,
    }, extra || {});
  }

  const PLOT_CONFIG = { displaylogo: false, responsive: true, toImageButtonOptions: { format: "svg" } };

  /* A horizontal bar per label, every label shown: the plot grows with the
     number of bars rather than letting Plotly thin the axis. What has not
     reacted is drawn in the muted color, apart from what formed. */
  function barh(id, labels, values, title, xtitle, xtype) {
    const height = Math.max(220, 20 * labels.length + 100);
    const layout = plotLayout({ title: { text: title, font: { size: 13 } }, height });
    layout.xaxis.title = { text: xtitle };
    layout.yaxis.dtick = 1;
    layout.yaxis.autorange = "reversed";
    if (xtype) layout.xaxis.type = xtype;
    const colors = labels.map((l) => (l.startsWith("unreacted ") ? css("--muted") : css("--accent")));
    const trace = { type: "bar", orientation: "h", x: values, y: labels, marker: { color: colors },
      hovertemplate: "%{y}: %{x:.4g}<extra></extra>" };
    show(id, height);
    Plotly.react(id, [trace], layout, PLOT_CONFIG);
  }

  /* Curves along a path: stacked areas for what adds up (masses, volumes), or
     plain lines; a legend beside them. */
  function curves(id, x, series, title, xtitle, ytitle, opts) {
    const o = opts || {};
    const layout = plotLayout({ title: { text: title, font: { size: 13 } }, height: 360, showlegend: series.length > 1,
      legend: { font: { size: 10 } } });
    layout.xaxis.title = { text: xtitle };
    layout.yaxis.title = { text: ytitle };
    if (o.log) layout.yaxis.type = "log";
    const traces = series.map((s) => Object.assign({
      type: "scatter", mode: o.stack ? "lines" : "lines+markers", x, y: s.values, name: s.label,
      line: { width: o.stack ? 0.5 : 2 }, hovertemplate: s.label + ": %{y:.4g}<extra></extra>",
    }, o.stack ? { stackgroup: "one" } : {}, s.label.startsWith("unreacted ") ? { line: { color: css("--muted"), width: 0.5 } } : {}));
    show(id, 360);
    Plotly.react(id, traces, layout, PLOT_CONFIG);
  }

  function show(id, height) {
    const d = document.getElementById(id);
    d.hidden = false;
    d.style.height = height + "px";
  }

  function hide(id) { const d = document.getElementById(id); d.hidden = true; Plotly.purge(d); }

  function renderResults() {
    const r = S.results;
    if (!r) return;
    const path = r.kind === "path";
    document.querySelectorAll("#result-tabs button").forEach((b) => {
      const only = b.dataset.only;
      b.hidden = !!only && only !== r.kind;
    });
    if (path) renderPath_(r); else renderEquilibrium(r);
    const current = document.querySelector("#result-tabs button.on");
    if (current && current.hidden) showResult("phases");
  }

  function renderEquilibrium(r) {
    const binder = r.binder_mass_g ? fmt(r.binder_mass_g, 0) : "100";
    const pm = r.phase_masses;
    barh("plot-masses", pm.labels, pm.values, t("plot_masses"), t("plot_masses_axis", { binder }));
    const vf = r.volume_fractions;
    if (vf.labels) {
      $("#volumes-note").textContent = "";
      barh("plot-volumes", vf.labels, vf.values.map((v) => (v === null ? null : 100 * v)), t("plot_volumes"), t("plot_volumes_axis"));
    } else {
      hide("plot-volumes");
      $("#volumes-note").textContent = vf.reason;
    }

    const ps = r.pore_solution;
    fill($("#pore-summary"),
      summaryItem(t("ph"), fmt(ps.pH, 2)),
      summaryItem(t("ionic_strength"), withUnit(fmt(r.ionic_strength, 3), "mol/kg")),
      summaryItem(t("porosity"), withUnit(fmt(percent(r.porosity.total), 1), "%"), t("porosity_unknown")),
      summaryItem(t("bound_water"), withUnit(fmt(r.bound_water, 3), "g/g")));
    // What is dissolved below a trace would squeeze the log axis to nothing.
    const mm = ps.molalities.values.map((v) => (v === null ? null : 1000 * v));
    const shown = ps.molalities.labels.map((k, i) => [k, mm[i]]).filter(([, v]) => v !== null && v >= TRACE);
    barh("plot-pore", shown.map((x) => x[0]), shown.map((x) => x[1]), t("plot_pore"), t("plot_pore_axis"), "log");
    hide("plot-ph");
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

  /* The axis of a path, from its label: an age, a fraction, an amount of a
     species, a step, a temperature, a ratio. */
  function pathAxis(label) {
    const key = "axis_" + label;
    const s = t(key);
    return s === key ? t("axis_added", { species: label }) : s;
  }

  function renderPath_(r) {
    const x = r.parameter;
    const xt = pathAxis(r.label);
    const binder = r.binder_mass_g ? fmt(r.binder_mass_g, 0) : "100";
    const pm = r.phase_masses;
    curves("plot-masses", x, pm.labels.map((l, k) => ({ label: l, values: pm.series[k] })), t("plot_masses"), xt,
      t("plot_masses_axis", { binder }), { stack: true });
    const vf = r.volume_fractions;
    if (vf.labels) {
      $("#volumes-note").textContent = "";
      curves("plot-volumes", x, vf.labels.map((l, k) => ({ label: l, values: vf.series[k].map((v) => (v === null ? null : 100 * v)) })),
        t("plot_volumes"), xt, t("plot_volumes_axis"), { stack: true });
    } else {
      hide("plot-volumes");
      $("#volumes-note").textContent = vf.reason;
    }
    fill($("#pore-summary"));
    curves("plot-ph", x, [{ label: "pH", values: r.pH }], t("ph"), xt, "pH");
    const mol = r.molalities;
    const lines = mol.labels.map((l, k) => ({ label: l, values: mol.series[k].map((v) => (v === null ? null : 1000 * v)) }))
      .filter((s) => s.values.some((v) => v !== null && v >= TRACE));
    curves("plot-pore", x, lines, t("plot_pore"), xt, t("plot_pore_axis"), { log: true });
    fill($("#pore-table"));
    curves("plot-porosity", x, [
      { label: t("porosity") + " (%)", values: r.porosity.map(percent) },
      { label: t("bound_water") + " (g/g × 100)", values: r.bound_water.map(percent) },
    ], t("porosity_water"), xt, "%");
    const head = [xt, t("certified"), "pH", t("porosity") + " (%)", t("bound_water") + " (g/g)"].concat(pm.labels);
    fillTable("#path-table", head, x.map((v, i) => [
      { text: fmt(v, 3), cls: "" }, { text: r.certified[i] ? "✓" : "✗", cls: r.certified[i] ? "ok" : "warn" },
      { text: fmt(r.pH[i], 2), cls: "" }, { text: fmt(percent(r.porosity[i]), 1), cls: "" }, { text: fmt(r.bound_water[i], 3), cls: "" },
    ].concat(pm.series.map((s) => ({ text: fmt(s[i], 2), cls: "" })))));
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
    fill($(sel),
      el("thead", {}, el("tr", {}, head.map((h) => el("th", { text: h })))),
      el("tbody", {}, rows.map((row) => el("tr", {}, row.map((c) =>
        typeof c === "object" ? el("td", { class: "num " + c.cls, text: c.text }) : el("td", { text: c }))))));
  }

  /* Copy a table as tab-separated text, for a spreadsheet. */
  async function copyTable(sel) {
    const rows = [...document.querySelectorAll(sel + " tr")].map((tr) => [...tr.children].map((c) => c.textContent).join("\t"));
    try { await navigator.clipboard.writeText(rows.join("\n")); toast(t("copied_rows")); }
    catch (e) { toast(t("copy_failed"), true); }
  }

  /* ── files ────────────────────────────────────────────────────── */

  const dirOf = (p) => (p ? p.replace(/[\\/][^\\/]*$/, "") : "");

  async function save() {
    if (!S.path) return saveAs();
    await flush();
    const { status, data } = await Studio.post("/api/save", {});
    if (status === 200) { applyState(data, false); toast(t("saved", { path: data.path })); }
    else toast(data.error || t("error"), true);
  }

  async function saveAs() {
    const path = await Picker.pick("save", dirOf(S.path), S.stem + ".jl");
    if (!path) return;
    await flush();
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

  /* ── tabs, theme, language, expert fields ─────────────────────── */

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

  function applyExpert(on) {
    document.body.classList.toggle("expert-on", on);
    $("#expert").checked = on;
    try { localStorage.setItem("clstudio.expert", on ? "1" : "0"); } catch (e) { /* not remembered */ }
  }
  function savedExpert() {
    try { return localStorage.getItem("clstudio.expert") === "1"; } catch (e) { return false; }
  }

  async function toggleLanguage() {
    await Studio.setLanguage(Studio.language() === "en" ? "fr" : "en");
    $("#lang").textContent = Studio.language() === "en" ? "FR" : "EN";
    renderPath();
    renderProblems();
    renderForms();
    renderDiagnostics();
    renderResults();
    renderLibrary();
  }

  /* ── start ────────────────────────────────────────────────────── */

  document.addEventListener("DOMContentLoaded", async () => {
    await Studio.setLanguage(Studio.savedLanguage());
    $("#lang").textContent = Studio.language() === "en" ? "FR" : "EN";
    applyTheme(savedTheme());
    applyExpert(savedExpert());
    document.querySelectorAll("#steps button").forEach((b) => b.addEventListener("click", () => showTab(b.dataset.tab)));
    document.querySelectorAll("#result-tabs button").forEach((b) => b.addEventListener("click", () => showResult(b.dataset.result)));
    document.querySelectorAll("#modes button").forEach((b) => b.addEventListener("click", () => setMode(b.dataset.mode)));
    document.querySelectorAll("[data-copy]").forEach((b) => b.addEventListener("click", () => copyTable(b.dataset.copy)));
    $("#run").addEventListener("click", run);
    $("#stop").addEventListener("click", stop);
    $("#save").addEventListener("click", save);
    $("#saveas").addEventListener("click", saveAs);
    $("#open").addEventListener("click", open);
    $("#open2").addEventListener("click", open);
    $("#material-add").addEventListener("click", addMaterial);
    $("#title").addEventListener("input", (e) => { S.model.title = e.target.value; edited(); });
    $("#expert").addEventListener("change", (e) => applyExpert(e.target.checked));
    $("#copy").addEventListener("click", async () => {
      try { await navigator.clipboard.writeText(S.script); toast(t("copied")); }
      catch (e) { toast(t("copy_failed"), true); }
    });
    $("#theme").addEventListener("click", () => applyTheme(THEMES[(THEMES.indexOf(savedTheme()) + 1) % THEMES.length]));
    $("#lang").addEventListener("click", toggleLanguage);

    const [state, schema] = await Promise.all([Studio.get("/api/state"), Studio.get("/api/schema")]);
    if (state.status !== 200 || schema.status !== 200) {
      setBadge("bad", t("status_denied"));
      toast((state.data && state.data.error) || t("error"), true);
      return;
    }
    S.schema = schema.data;
    Forms.setSchema(schema.data);
    applyState(state.data, true);
    renderDiagnostics();
    loadLibrary();
    poll();
    setInterval(poll, 1500);
  });
})();
