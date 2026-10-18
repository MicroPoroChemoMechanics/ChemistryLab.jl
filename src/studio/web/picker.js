/* SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright © 2026 Jean-François Barthélémy (Cerema, UMR MCD)
 *
 * The file picker, after MFH Studio's. A browser hands a page the content of
 * a file but never its path, and saving back to the same file needs the path;
 * so the server lists its own folders. It also keeps working on a remote
 * machine, where a native dialog would show the wrong computer's disks.
 */

const Picker = (() => {
  const state = { path: "", selected: null, mode: "open", resolve: null };
  const q = (s) => document.querySelector(s);

  function row(icon, label, onActivate) {
    const d = document.createElement("div");
    d.className = "row";
    const i = document.createElement("span");
    i.className = "ico";
    i.textContent = icon;
    const s = document.createElement("span");
    s.textContent = label;
    d.append(i, s);
    d.addEventListener("click", () => onActivate(d, false));
    d.addEventListener("dblclick", () => onActivate(d, true));
    return d;
  }

  const join = (dir, name) => {
    const sep = dir.includes("\\") ? "\\" : "/";
    return dir.replace(/[\\/]+$/, "") + sep + name;
  };

  async function show(dir) {
    const { status, data } = await Studio.get("/api/browse?path=" + encodeURIComponent(dir || ""));
    if (status !== 200) {
      toast(data.error || t("error"), true);
      return;
    }
    state.path = data.path;
    state.selected = null;
    q("#picker-crumb").textContent = data.path;
    q("#picker-places").replaceChildren(
      ...data.places.map((p) => {
        const b = document.createElement("button");
        b.className = "small";
        b.textContent = t(p.label === "Home" ? "place_home" : "place_cwd");
        b.title = p.path;
        b.addEventListener("click", () => show(p.path));
        return b;
      })
    );
    const list = q("#picker-list");
    const rows = [];
    if (data.parent && data.parent !== data.path) rows.push(row("↑", "..", () => show(data.parent)));
    for (const name of data.folders) rows.push(row("📁", name, () => show(join(data.path, name))));
    for (const name of data.files) {
      rows.push(
        row("📄", name, (el, double) => {
          list.querySelectorAll(".row.on").forEach((n) => n.classList.remove("on"));
          el.classList.add("on");
          state.selected = join(data.path, name);
          q("#picker-name").value = name;
          if (double) accept();
        })
      );
    }
    if (!rows.length) {
      const e = document.createElement("div");
      e.className = "empty";
      e.textContent = t(state.mode === "save" ? "picker_empty_save" : "picker_empty_open");
      rows.push(e);
    }
    list.replaceChildren(...rows);
  }

  function accept() {
    const name = q("#picker-name").value.trim();
    let path = state.selected;
    if (state.mode === "save" || !path) {
      if (!name) {
        toast(t("picker_need_name"), true);
        return;
      }
      path = /[\\/]/.test(name) ? name : join(state.path, name);
    }
    close(path);
  }

  function close(result) {
    q("#picker").hidden = true;
    document.removeEventListener("keydown", onKey);
    if (state.resolve) state.resolve(result);
    state.resolve = null;
  }

  function onKey(e) {
    if (e.key === "Escape") close(null);
    else if (e.key === "Enter") accept();
  }

  /** Open the dialog; resolves with a path, or null when dismissed. */
  function pick(mode, startDir, suggestedName) {
    state.mode = mode;
    q("#picker-title").textContent = t(mode === "save" ? "picker_save" : "picker_open");
    q("#picker-ok").textContent = t(mode === "save" ? "save" : "open_short");
    q("#picker-name").value = suggestedName || "";
    q("#picker").hidden = false;
    document.addEventListener("keydown", onKey);
    show(startDir || "");
    return new Promise((resolve) => { state.resolve = resolve; });
  }

  document.addEventListener("DOMContentLoaded", () => {
    q("#picker-close").addEventListener("click", () => close(null));
    q("#picker-ok").addEventListener("click", accept);
    q("#picker").addEventListener("click", (e) => { if (e.target.id === "picker") close(null); });
  });

  return { pick };
})();
