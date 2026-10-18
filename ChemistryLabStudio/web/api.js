/* SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright © 2026 Jean-François Barthélémy (Cerema, UMR MCD)
 *
 * What every part of the page shares: the session token, the calls to the
 * server, the translations and the toasts.
 *
 * The token arrives in the address the terminal prints. It is kept for the tab
 * and removed from the address bar, so that a screenshot or a shared link does
 * not carry it; every call to `/api` sends it back in a header.
 */

const Studio = (() => {
  const params = new URLSearchParams(location.search);
  let token = params.get("token");
  try {
    if (token) sessionStorage.setItem("clstudio.token", token);
    else token = sessionStorage.getItem("clstudio.token");
  } catch (e) {
    /* storage unavailable: the token lives as long as the page */
  }
  if (params.has("token")) history.replaceState(null, "", location.pathname);

  async function call(method, path, body) {
    const init = { method, headers: { "X-Studio-Token": token || "" } };
    if (body !== undefined) {
      init.headers["Content-Type"] = "application/json";
      init.body = JSON.stringify(body);
    }
    const r = await fetch(path, init);
    let data = null;
    try {
      data = await r.json();
    } catch (e) {
      data = { error: r.statusText };
    }
    return { status: r.status, data };
  }

  const get = (path) => call("GET", path);
  const post = (path, body) => call("POST", path, body === undefined ? {} : body);

  /* ── translations ─────────────────────────────────────────────── */

  let lang = "en";
  let words = {};
  let english = {};

  async function loadWords(l) {
    const r = await fetch(`/i18n/${l}.json`);
    return r.ok ? r.json() : {};
  }

  /** The text of `key` in the current language, `{name}` replaced from `vars`. */
  function t(key, vars) {
    let s = words[key] !== undefined ? words[key] : english[key] !== undefined ? english[key] : key;
    if (vars) for (const [k, v] of Object.entries(vars)) s = s.split(`{${k}}`).join(String(v));
    return s;
  }

  function translatePage() {
    document.documentElement.lang = lang;
    document.querySelectorAll("[data-i18n]").forEach((el) => { el.textContent = t(el.dataset.i18n); });
    document.querySelectorAll("[data-i18n-title]").forEach((el) => { el.title = t(el.dataset.i18nTitle); });
  }

  async function setLanguage(l) {
    english = Object.keys(english).length ? english : await loadWords("en");
    words = l === "en" ? english : await loadWords(l);
    lang = l;
    try { localStorage.setItem("clstudio.lang", l); } catch (e) { /* not remembered */ }
    translatePage();
  }

  function savedLanguage() {
    try { return localStorage.getItem("clstudio.lang") || "en"; } catch (e) { return "en"; }
  }

  /* ── toasts ───────────────────────────────────────────────────── */

  function toast(text, bad) {
    const d = document.createElement("div");
    d.textContent = text;
    if (bad) d.className = "bad";
    document.getElementById("toast").append(d);
    setTimeout(() => d.remove(), bad ? 9000 : 3500);
  }

  return { get, post, t, setLanguage, savedLanguage, language: () => lang, translatePage, toast };
})();

const t = Studio.t;
const toast = Studio.toast;
