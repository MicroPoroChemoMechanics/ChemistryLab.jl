/* SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright © 2026 Jean-François Barthélémy (Cerema, UMR MCD)
 *
 * Forms drawn from the schema the server sends. A choice (a database, a law of
 * the degree of reaction, a calculation, a rate law…) arrives as a family of
 * kinds, each with its fields: kind, default, unit, bounds, level. The page
 * draws a list of the kinds and the fields of the chosen one; nothing here
 * knows any particular choice, so a choice added in the server needs no change
 * here. What a form cannot know (the constituents a law may apply to, the
 * units of a table, the fields another panel draws) comes from a context the
 * caller passes.
 *
 * Labels: `k_<kind>` names a kind and `hk_<kind>` explains it, `f_<field>`
 * names a field and `h_<field>` explains it, `o_<option>` names an option.
 */

const Forms = (() => {
  let schema = {};
  let catalog = null;

  const setSchema = (s) => { schema = s || {}; };
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
    for (const c of children.flat()) if (c !== null && c !== undefined && c !== false) e.append(c);
    return e;
  }

  const label = (key, fallback) => { const s = t(key); return s === key ? fallback : s; };
  const kindLabel = (kind) => label("k_" + kind, kind);
  const kindHelp = (kind) => label("hk_" + kind, "");
  const fieldLabel = (name) => label("f_" + name, name);
  const help = (name) => label("h_" + name, "");
  const optionLabel = (o) => label("o_" + o, o);

  const membersOf = (family) => schema[family] || [];
  const memberOf = (family, kind) => membersOf(family).find((m) => m.kind === kind);

  /** The default object of a kind: its fields at their declared defaults. */
  function defaults(family, kind) {
    const m = memberOf(family, kind);
    const out = { kind };
    if (m) for (const f of m.fields) out[f.name] = structuredClone(f.default);
    return out;
  }

  /** The options of a `choice` field: fixed, or read from the catalog. */
  function options(f) {
    if (f.options && f.options.length) return f.options.slice();
    if (f.source && catalog && catalog[f.source]) {
      return catalog[f.source].filter((x) => typeof x === "string" || !x.unavailable)
        .map((x) => (typeof x === "string" ? x : x.name));
    }
    return [];
  }

  const numberText = (v) => (v === null || v === undefined ? "" : String(v));
  const readNumber = (s) => Number(String(s).trim().replace(",", "."));

  function numberInput(f, value, onChange) {
    return make("input", {
      type: "text", inputmode: "decimal", value: numberText(value), spellcheck: "false",
      class: "num",
      oninput: (e) => {
        const x = readNumber(e.target.value);
        const bad = e.target.value.trim() === "" || !Number.isFinite(x) ||
          (f.min !== null && x < f.min) || (f.max !== null && x > f.max);
        e.target.classList.toggle("invalid", bad);
        if (e.target.value.trim() !== "" && Number.isFinite(x)) onChange(f.kind === "integer" ? Math.round(x) : x);
      },
    });
  }

  function numbersInput(value, onChange) {
    return make("input", {
      type: "text", value: (value || []).join("  "), spellcheck: "false", title: t("numbers_hint"),
      onchange: (e) => onChange(e.target.value.split(/[\s;]+/).filter((x) => x !== "").map(readNumber).filter(Number.isFinite)),
    });
  }

  function choiceInput(f, value, onChange) {
    const opts = options(f);
    if (value !== null && value !== undefined && !opts.includes(value)) opts.unshift(value);
    return make("select", { onchange: (e) => onChange(e.target.value) },
      opts.map((o) => make("option", { value: o, text: optionLabel(o), selected: o === value })));
  }

  /** A nested choice: "none" when the field may be empty, then the kinds. */
  function nestedInput(f, value, onChange, ctx) {
    const box = make("div", { class: "nested" });
    const optional = f.default === null;
    function draw() {
      const select = make("select", {
        onchange: (e) => {
          value = e.target.value === "" ? null : defaults(f.source, e.target.value);
          onChange(value);
          draw();
        },
      },
      optional ? make("option", { value: "", text: t("none"), selected: !value }) : null,
      membersOf(f.source).map((m) => make("option", { value: m.kind, text: kindLabel(m.kind), selected: !!value && m.kind === value.kind })));
      box.replaceChildren(...[select, value ? fieldsOf(f.source, value, (v) => { value = v; onChange(value); }, ctx) : null].filter(Boolean));
    }
    draw();
    return box;
  }

  /**
   * A choice per key (a law per constituent): the keys come from the context,
   * each with "none" or a kind of the family; the entries keep the order of
   * the keys.
   */
  function keyedInput(f, value, onChange, ctx) {
    const box = make("div", { class: "keyed" });
    const keys = ctx && ctx.keys ? ctx.keys(f) : [];
    let entries = (value || []).slice();
    const all = keys.concat(entries.map((e) => e.key).filter((k) => !keys.includes(k)));
    function set(key, v) {
      entries = entries.filter((e) => e.key !== key);
      if (v) entries.push({ key, value: v });
      entries.sort((a, b) => all.indexOf(a.key) - all.indexOf(b.key));
      onChange(entries);
    }
    function draw() {
      if (!all.length) {
        box.replaceChildren(make("p", { class: "hint", text: ctx && ctx.noKeys ? ctx.noKeys(f) : t("no_keys") }));
        return;
      }
      box.replaceChildren(...all.map((key) => {
        const entry = entries.find((e) => e.key === key);
        const note = ctx && ctx.note ? ctx.note(f, key) : "";
        const select = make("select", {
          onchange: (e) => { set(key, e.target.value === "" ? null : defaults(f.source, e.target.value)); draw(); },
        },
        make("option", { value: "", text: ctx && ctx.noneLabel ? ctx.noneLabel(f) : t("none"), selected: !entry }),
        membersOf(f.source).map((m) => make("option", { value: m.kind, text: kindLabel(m.kind), selected: !!entry && m.kind === entry.value.kind })));
        return make("div", { class: "keyed-row" },
          make("div", { class: "keyed-head" }, make("b", { text: key }), note ? make("span", { class: "hint", text: note }) : null, select),
          entry ? fieldsOf(f.source, entry.value, (v) => set(key, v), ctx) : null);
      }));
    }
    draw();
    return box;
  }

  /** One field: its label, its input, its unit, its help. */
  function field(f, value, onChange, ctx) {
    let input;
    if (f.kind === "number" || f.kind === "integer") input = numberInput(f, value, onChange);
    else if (f.kind === "bool") input = make("input", { type: "checkbox", checked: !!value, onchange: (e) => onChange(e.target.checked) });
    else if (f.kind === "choice") input = choiceInput(f, value, onChange);
    else if (f.kind === "numbers") input = numbersInput(value, onChange);
    else if (f.kind === "texts") {
      input = make("input", {
        type: "text", value: (value || []).join(", "), spellcheck: "false",
        onchange: (e) => onChange(e.target.value.split(",").map((x) => x.trim()).filter((x) => x)),
      });
    } else if (f.kind === "rows") input = rowsTable(value || [], onChange, { units: ctx && ctx.units ? ctx.units(f) : null });
    else if (f.kind === "nested") input = nestedInput(f, value, onChange, ctx);
    else if (f.kind === "keyed") input = keyedInput(f, value, onChange, ctx);
    else input = make("input", { type: "text", value: value || "", spellcheck: "false", oninput: (e) => onChange(e.target.value) });
    const h = help(f.name);
    const wide = ["rows", "keyed", "nested"].includes(f.kind);
    const inline = f.kind === "bool";
    return make("div", {
      class: "field" + (f.level === "expert" ? " expert" : "") + (wide ? " wide" : "") + (inline ? " inline" : ""),
      title: h || null, "data-field": f.name,
    },
    make("label", {}, inline ? input : null, fieldLabel(f.name), f.unit ? make("span", { class: "unit", text: " (" + f.unit + ")" }) : null),
    inline ? null : input);
  }

  /** The fields of an object of a family, without its kind list. */
  function fieldsOf(family, value, onChange, ctx) {
    const m = memberOf(family, value.kind);
    const skip = (ctx && ctx.skip) || [];
    return make("div", { class: "fields" },
      (m ? m.fields : []).filter((f) => !skip.includes(f.name))
        .map((f) => field(f, value[f.name], (v) => { value[f.name] = v; onChange(value); }, ctx)));
  }

  /**
   * A choice of a family: a list of its kinds (or of those `opts.only` keeps),
   * then the fields of the chosen one. Changing the kind starts from the
   * defaults of the new one.
   */
  function kinded(family, value, onChange, opts) {
    const o = opts || {};
    const box = make("div", { class: "kinded" });
    function draw() {
      const members = membersOf(family).filter((m) => !o.only || o.only(m) || m.kind === value.kind);
      const select = make("select", {
        onchange: (e) => { value = defaults(family, e.target.value); onChange(value); draw(); },
      }, members.map((x) => make("option", { value: x.kind, text: kindLabel(x.kind), selected: x.kind === value.kind })));
      const h = kindHelp(value.kind);
      // `replaceChildren` would write a null as the text "null"
      box.replaceChildren(...[
        o.hideKind ? null : make("div", { class: "field" },
          o.label ? make("label", { text: o.label }) : null, select, h ? make("p", { class: "hint", text: h }) : null),
        fieldsOf(family, value, (v) => { value = v; onChange(value); }, o.ctx),
      ].filter(Boolean));
    }
    draw();
    return box;
  }

  /* ── the input table ──────────────────────────────────────────── */

  const UNITS = ["%", "g", "kg", "mg", "mol", "mmol", "µmol", "mol/kg", "mmol/kg", "mol/L", "mmol/L"];

  /**
   * An editable table of rows {formula, amount, unit, notation}, with a paste
   * area: what is pasted is read by the server, shown, and kept on
   * confirmation. Ctrl-Z in the table undoes the last change; Copy puts the
   * rows on the clipboard as tab-separated text, for a spreadsheet.
   */
  function rowsTable(rows, onChange, opts) {
    const history = [];
    const units = (opts && opts.units) || UNITS;
    const box = make("div", { class: "rows" });
    const commit = (next) => { history.push(structuredClone(rows)); rows = next; onChange(rows); draw(); };

    function cellInput(r, i, key, number) {
      return make("input", {
        type: "text", value: r[key], spellcheck: "false", class: number ? "num" : null,
        onchange: (e) => {
          const next = structuredClone(rows);
          next[i][key] = number ? readNumber(e.target.value) : e.target.value.trim();
          commit(next);
        },
      });
    }

    function draw() {
      const body = rows.map((r, i) => make("tr", {},
        make("td", {}, cellInput(r, i, "formula", false)),
        make("td", {}, cellInput(r, i, "amount", true)),
        make("td", {}, make("select", {
          onchange: (e) => { const next = structuredClone(rows); next[i].unit = e.target.value; commit(next); },
        }, units.concat(units.includes(r.unit) ? [] : [r.unit]).map((u) => make("option", { value: u, text: u, selected: u === r.unit })))),
        make("td", {}, make("select", {
          title: t("notation_hint"),
          onchange: (e) => { const next = structuredClone(rows); next[i].notation = e.target.value; commit(next); },
        }, ["auto", "mendeleev", "cement"].map((n) => make("option", { value: n, text: t("notation_" + n), selected: n === r.notation })))),
        make("td", {}, make("button", {
          class: "small", text: "−", title: t("remove_row"),
          onclick: () => { const next = structuredClone(rows); next.splice(i, 1); commit(next); },
        }))));
      const percent = rows.filter((r) => r.unit === "%");
      const total = percent.reduce((s, r) => s + (Number(r.amount) || 0), 0);
      box.replaceChildren(
        make("table", { class: "data compact rows-table" },
          make("thead", {}, make("tr", {}, [t("formula"), t("amount"), t("unit"), t("notation"), ""].map((h) => make("th", { text: h })))),
          make("tbody", {}, body)),
        make("div", { class: "row" },
          make("button", {
            class: "small", text: t("add_row"),
            onclick: () => commit(rows.concat([{ formula: "", amount: 0, unit: (rows[rows.length - 1] || {}).unit || units[0], notation: "auto" }])),
          }),
          make("button", { class: "small", text: t("paste_rows"), onclick: pasteArea }),
          make("button", { class: "small", text: t("copy_rows"), onclick: copyRows, disabled: !rows.length }),
          percent.length ? make("span", { class: "hint" + (total > 100.5 ? " problem" : ""), text: t("percent_total", { total: total.toFixed(2) }) }) : null));
    }

    async function copyRows() {
      try {
        await navigator.clipboard.writeText(rows.map((r) => `${r.formula}\t${r.amount}\t${r.unit}`).join("\n"));
        toast(t("copied_rows"));
      } catch (e) { toast(t("copy_failed"), true); }
    }

    function pasteArea() {
      const area = make("textarea", { rows: 6, placeholder: t("paste_placeholder"), spellcheck: "false" });
      const preview = make("div", { class: "hint" });
      let parsed = null;
      area.addEventListener("input", async () => {
        const unit = (rows[rows.length - 1] || {}).unit || units[0];
        const { status, data } = await Studio.post("/api/rows/parse", { text: area.value, unit });
        if (status !== 200) { preview.textContent = data.error || t("error"); parsed = null; return; }
        parsed = data.rows;
        preview.replaceChildren(...[
          make("span", { text: t("paste_read", { n: data.rows.length, skipped: data.skipped.length }) }),
          make("ul", {}, data.rows.slice(0, 12).map((r) => make("li", { text: `${r.formula} → ${r.amount} ${r.unit}` }))),
          data.skipped.length ? make("p", { class: "problem", text: t("paste_skipped", { lines: data.skipped.slice(0, 5).join(" · ") }) }) : null,
        ].filter(Boolean));
      });
      const dialog = make("div", { class: "card paste" },
        area, preview,
        make("div", { class: "row" },
          make("button", {
            class: "small primary", text: t("paste_keep"),
            onclick: () => { if (parsed && parsed.length) commit(rows.concat(parsed)); dialog.remove(); },
          }),
          make("button", {
            class: "small", text: t("paste_replace"),
            onclick: () => { if (parsed && parsed.length) commit(parsed); dialog.remove(); },
          }),
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

  return {
    setSchema, setCatalog, kinded, field, fieldsOf, rowsTable, defaults, membersOf, memberOf, make,
    kindLabel, kindHelp, fieldLabel, optionLabel,
  };
})();
