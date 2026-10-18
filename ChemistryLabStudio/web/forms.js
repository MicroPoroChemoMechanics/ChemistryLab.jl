/* SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright © 2026 Jean-François Barthélémy (Cerema, UMR MCD)
 *
 * Forms drawn from the schema the server sends. A choice (a law of the degree
 * of reaction, a process, a rate law…) arrives as a family of kinds, each with
 * its fields: kind, default, unit, bounds, level. The page draws a list of the
 * kinds and the fields of the chosen one; nothing here knows any particular
 * choice, so a choice added in the server needs no change here.
 *
 * Labels: `k_<kind>` names a kind, `f_<field>` a field, `h_<field>` its help.
 */

const Forms = (() => {
  let schema = {};
  let catalog = null;

  const setSchema = (s) => { schema = s; };
  const setCatalog = (c) => { catalog = c; };

  function make(tag, attrs, ...children) {
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

  const label = (key, fallback) => { const s = t(key); return s === key ? fallback : s; };
  const kindLabel = (kind) => label("k_" + kind, kind);
  const fieldLabel = (name) => label("f_" + name, name);
  const help = (name) => { const s = t("h_" + name); return s === "h_" + name ? "" : s; };

  function membersOf(family) { return schema[family] || []; }
  function memberOf(family, kind) { return membersOf(family).find((m) => m.kind === kind); }

  /** The default object of a kind: its fields at their declared defaults. */
  function defaults(family, kind) {
    const m = memberOf(family, kind);
    const out = { kind };
    if (m) for (const f of m.fields) out[f.name] = structuredClone(f.default);
    return out;
  }

  /** The options of a `choice` field: fixed, or read from the catalog. */
  function options(f) {
    if (f.options && f.options.length) return f.options;
    if (f.source && catalog && catalog[f.source]) {
      return catalog[f.source].map((x) => (typeof x === "string" ? x : x.name));
    }
    return [];
  }

  function numberInput(f, value, onChange) {
    return make("input", {
      type: "number", step: "any", value: value === null || value === undefined ? "" : value,
      min: f.min === null ? undefined : f.min, max: f.max === null ? undefined : f.max,
      oninput: (e) => { if (e.target.value !== "") onChange(Number(e.target.value)); },
    });
  }

  /** One field: its label, its input, its unit, its help. */
  function field(f, value, onChange) {
    let input;
    if (f.kind === "number" || f.kind === "integer") {
      input = numberInput(f, value, (v) => onChange(f.kind === "integer" ? Math.round(v) : v));
    } else if (f.kind === "bool") {
      input = make("input", { type: "checkbox", checked: !!value, onchange: (e) => onChange(e.target.checked) });
    } else if (f.kind === "choice") {
      const opts = options(f);
      if (opts.length && !opts.includes(value)) opts.unshift(value);
      input = make("select", { onchange: (e) => onChange(e.target.value) },
        opts.map((o) => make("option", { value: o, text: label("o_" + o, o), selected: o === value })));
    } else if (f.kind === "numbers") {
      input = make("input", {
        type: "text", value: (value || []).join(" "), spellcheck: "false",
        title: t("numbers_hint"),
        onchange: (e) => onChange(e.target.value.split(/[\s;]+/).filter((x) => x !== "").map((x) => Number(x.replace(",", ".")))),
      });
    } else if (f.kind === "texts") {
      input = make("input", {
        type: "text", value: (value || []).join(", "), spellcheck: "false",
        onchange: (e) => onChange(e.target.value.split(",").map((x) => x.trim()).filter((x) => x)),
      });
    } else if (f.kind === "rows") {
      input = rowsTable(value || [], onChange, {});
    } else if (f.kind === "nested") {
      input = value ? kinded(f.source, value, onChange) : make("span", { class: "hint", text: "—" });
    } else {
      input = make("input", { type: "text", value: value || "", spellcheck: "false", oninput: (e) => onChange(e.target.value) });
    }
    const h = help(f.name);
    return make("div", { class: "field" + (f.level === "expert" ? " expert" : ""), title: h || null },
      make("label", {}, fieldLabel(f.name), f.unit ? make("span", { class: "unit", text: " (" + f.unit + ")" }) : null),
      input,
      h && f.kind !== "rows" ? make("p", { class: "hint", text: h }) : null);
  }

  /**
   * A choice of a family: a list of its kinds, then the fields of the chosen
   * one. Changing the kind starts from the defaults of the new one.
   */
  function kinded(family, value, onChange, opts) {
    const members = membersOf(family);
    const box = make("div", { class: "kinded" });
    function draw() {
      const m = memberOf(family, value.kind);
      const select = make("select", {
        onchange: (e) => { value = defaults(family, e.target.value); onChange(value); draw(); },
      }, members.map((x) => make("option", { value: x.kind, text: kindLabel(x.kind), selected: x.kind === value.kind })));
      const kindHelp = t("hk_" + value.kind);
      box.replaceChildren(
        (opts && opts.hideKind) ? null : make("div", { class: "field" },
          opts && opts.label ? make("label", { text: opts.label }) : null, select,
          kindHelp !== "hk_" + value.kind ? make("p", { class: "hint", text: kindHelp }) : null),
        ...(m ? m.fields : []).map((f) => field(f, value[f.name], (v) => { value[f.name] = v; onChange(value); }))
      );
    }
    draw();
    return box;
  }

  /* ── the input table ──────────────────────────────────────────── */

  const UNITS = ["%", "g", "kg", "mg", "mol", "mmol", "µmol", "mol/kg", "mmol/kg", "mol/L", "mmol/L"];

  /**
   * An editable table of rows {formula, amount, unit, notation}, with a paste
   * area: what is pasted is read by the server, shown, and kept on confirmation;
   * Ctrl-Z in the table undoes the last change.
   */
  function rowsTable(rows, onChange, opts) {
    const history = [];
    const box = make("div", { class: "rows" });
    const commit = (next) => { history.push(structuredClone(rows)); rows = next; onChange(rows); draw(); };

    function cellInput(r, i, key, type) {
      return make("input", {
        type, value: r[key], spellcheck: "false", step: type === "number" ? "any" : undefined,
        onchange: (e) => {
          const next = structuredClone(rows);
          next[i][key] = type === "number" ? Number(String(e.target.value).replace(",", ".")) : e.target.value;
          commit(next);
        },
      });
    }

    function draw() {
      const units = (opts && opts.units) || UNITS;
      const body = rows.map((r, i) => make("tr", {},
        make("td", {}, cellInput(r, i, "formula", "text")),
        make("td", {}, cellInput(r, i, "amount", "number")),
        make("td", {}, make("select", {
          onchange: (e) => { const next = structuredClone(rows); next[i].unit = e.target.value; commit(next); },
        }, units.map((u) => make("option", { value: u, text: u, selected: u === r.unit })))),
        make("td", {}, make("select", {
          title: t("notation_hint"),
          onchange: (e) => { const next = structuredClone(rows); next[i].notation = e.target.value; commit(next); },
        }, ["auto", "mendeleev", "cement"].map((n) => make("option", { value: n, text: t("notation_" + n), selected: n === r.notation })))),
        make("td", {}, make("button", { class: "small", text: "−", title: t("remove_row"),
          onclick: () => { const next = structuredClone(rows); next.splice(i, 1); commit(next); } }))
      ));
      const total = rows.filter((r) => r.unit === "%").reduce((s, r) => s + (Number(r.amount) || 0), 0);
      box.replaceChildren(
        make("table", { class: "data compact rows-table" },
          make("thead", {}, make("tr", {}, [t("formula"), t("amount"), t("unit"), t("notation"), ""].map((h) => make("th", { text: h })))),
          make("tbody", {}, body)),
        make("div", { class: "row" },
          make("button", { class: "small", text: t("add_row"),
            onclick: () => commit(rows.concat([{ formula: "", amount: 0, unit: (rows[rows.length - 1] || {}).unit || "%", notation: "auto" }])) }),
          make("button", { class: "small", text: t("paste_rows"), onclick: pasteArea }),
          make("button", { class: "small", text: t("copy_rows"), onclick: copyRows }),
          rows.some((r) => r.unit === "%") ? make("span", { class: "hint" + (Math.abs(total - 100) > 0.5 ? " problem" : ""), text: t("percent_total", { total: total.toFixed(2) }) }) : null)
      );
    }

    async function copyRows() {
      try {
        await navigator.clipboard.writeText(rows.map((r) => `${r.formula}\t${r.amount}\t${r.unit}`).join("\n"));
        toast(t("copied_rows"));
      } catch (e) { toast(t("copy_failed"), true); }
    }

    function pasteArea() {
      const area = make("textarea", { rows: 6, placeholder: t("paste_placeholder") });
      const preview = make("div", { class: "hint" });
      let parsed = null;
      area.addEventListener("input", async () => {
        const unit = (rows[rows.length - 1] || {}).unit || "%";
        const { status, data } = await Studio.post("/api/rows/parse", { text: area.value, unit });
        if (status !== 200) { preview.textContent = data.error || t("error"); parsed = null; return; }
        parsed = data.rows;
        preview.replaceChildren(
          make("span", { text: t("paste_read", { n: data.rows.length, skipped: data.skipped.length }) }),
          make("ul", {}, data.rows.slice(0, 12).map((r) => make("li", { text: `${r.formula} → ${r.amount} ${r.unit}` })))
        );
      });
      const dialog = make("div", { class: "card paste" },
        area, preview,
        make("div", { class: "row" },
          make("button", { class: "small primary", text: t("paste_keep"), onclick: () => {
            if (parsed && parsed.length) commit((opts && opts.replace) ? parsed : rows.concat(parsed));
            dialog.remove();
          } }),
          make("button", { class: "small", text: t("cancel"), onclick: () => dialog.remove() })));
      box.append(dialog);
      area.focus();
    }

    box.addEventListener("keydown", (e) => {
      if ((e.ctrlKey || e.metaKey) && e.key === "z" && history.length && e.target.tagName !== "TEXTAREA") {
        e.preventDefault();
        rows = history.pop();
        onChange(rows);
        draw();
      }
    });
    draw();
    return box;
  }

  return { setSchema, setCatalog, kinded, field, rowsTable, defaults, membersOf, make };
})();
