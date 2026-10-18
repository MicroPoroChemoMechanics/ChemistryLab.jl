#!/usr/bin/env node
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright © 2026 Jean-François Barthélémy (Cerema, UMR MCD)
//
// The screenshots of ChemistryLab Studio in the documentation and the README,
// taken from a running studio so that they can be taken again whenever the
// page changes.
//
//   julia -e 'using ChemistryLab; clstudio(no_browser = true)'   # in another terminal
//   node screenshots.mjs --url 'http://127.0.0.1:8766/?token=…' [--out ../../../docs/src/assets/studio]
//
// Each picture is taken in the light and in the dark theme:
//   studio_equilibrium_{light,dark}.png  the default cement paste at 28 days;
//   studio_species_{light,dark}.png      the System tab: phases, species, diagnosis;
//   studio_path_{light,dark}.png         a carbonation path;
//   studio_kinetics_{light,dark}.png     a kinetic run, its degrees of reaction.
// A kinetic run compiles for a few minutes the first time: the script waits.
// The pictures committed are then reduced to 256 colors (a third of the size,
// no visible loss), for instance with Python's PIL:
//   Image.open(f).convert("RGB").quantize(256, dither=Image.Dither.NONE).save(f, optimize=True)
//
// Needs Node and playwright-core with a Chromium build. playwright-core is found
// through $PLAYWRIGHT_CORE, then in the global @playwright/mcp install, then by
// plain resolution.

import { createRequire } from "node:module";
import { execSync } from "node:child_process";
import { mkdirSync } from "node:fs";
import { dirname, join, resolve } from "node:path";
import { fileURLToPath } from "node:url";

const here = dirname(fileURLToPath(import.meta.url));
const args = Object.fromEntries(process.argv.slice(2).reduce((acc, a, i, all) =>
  (a.startsWith("--") ? acc.concat([[a.slice(2), all[i + 1]]]) : acc), []));
if (!args.url) throw new Error("give the address of a running studio: --url 'http://127.0.0.1:8766/?token=…'");
const out = resolve(args.out || join(here, "..", "..", "docs", "src", "assets", "studio"));
const VIEW = { width: 1600, height: 950 };

function playwright() {
  const tries = [
    process.env.PLAYWRIGHT_CORE,
    (() => { try { return join(execSync("npm root -g").toString().trim(), "@playwright/mcp/node_modules/playwright-core"); } catch { return null; } })(),
    "playwright-core",
  ].filter(Boolean);
  const require = createRequire(import.meta.url);
  for (const t of tries) {
    try { return require(t); } catch { /* next */ }
  }
  throw new Error("playwright-core not found; set PLAYWRIGHT_CORE to its directory");
}

// The studio's own calls, with the token the page keeps.
const api = (page, path, body) => page.evaluate(([p, b]) => (b === undefined ? Studio.get(p) : Studio.post(p, b)), [path, body]);

async function ready(page, minutes) {
  await page.waitForFunction(() => document.querySelector("#status").classList.contains("ok"), null, { timeout: minutes * 60000 });
}

async function setModel(page, change) {
  const { data } = await api(page, "/api/state");
  const model = change(data.model);
  await api(page, "/api/model", { model });
  await page.reload();
  await ready(page, 2);
}

async function run(page, minutes) {
  await page.click("#run");
  // the badge turns busy once the model is sent (a short run may be over
  // before it is seen), then ready again with the results drawn
  await page.waitForFunction(() => document.querySelector("#status").classList.contains("busy"), null, { timeout: 30000 })
    .catch(() => { /* already over */ });
  await page.waitForFunction(() => !document.querySelector("#run").disabled && document.querySelector("#status").classList.contains("ok"),
    null, { timeout: minutes * 60000 });
  await page.waitForTimeout(1500);   // the plots settle
}

async function shoot(page, name) {
  for (const theme of ["light", "dark"]) {
    // the page reads its theme when it loads, the plots when they are drawn:
    // a reload in the theme, then the calculation again
    await page.evaluate((t) => { try { localStorage.setItem("clstudio.theme", t); } catch (e) { /* none */ } }, theme);
    await page.reload();
    await ready(page, 2);
    await page.waitForTimeout(300);
    await name.prepare(page);
    await page.waitForTimeout(800);
    const file = join(out, `studio_${name.id}_${theme}.png`);
    await page.screenshot({ path: file, clip: { x: 0, y: 0, ...VIEW } });
    console.log("wrote", file);
  }
}

const tab = (page, t) => page.click(`#steps button[data-tab=${t}]`);
const result = (page, r) => page.click(`#result-tabs button[data-result=${r}]`);

const { chromium } = playwright();
mkdirSync(out, { recursive: true });
const browser = await chromium.launch();
const page = await browser.newPage({ viewport: VIEW });
try {
  await page.goto(args.url);
  await page.evaluate(() => { try { localStorage.setItem("clstudio.lang", "en"); localStorage.setItem("clstudio.expert", "0"); } catch (e) { /* none */ } });
  await page.reload();
  await ready(page, 10);
  await api(page, "/api/library", { id: "cem1-52.5n-gruyaert2010" });
  await page.reload();
  await ready(page, 2);

  // the default paste, one equilibrium
  const equilibrium = { id: "equilibrium", prepare: async (p) => { await tab(p, "materials"); await run(p, 10); await result(p, "phases"); } };
  await shoot(page, equilibrium);
  const species = {
    id: "species", prepare: async (p) => {
      await tab(p, "system");
      await p.fill("#species-browser input[type=search]", "C");
      await p.waitForTimeout(2500);   // the diagnosis
      await run(p, 10);
      await result(p, "pore");
    },
  };
  await shoot(page, species);

  // a carbonation path
  await setModel(page, (m) => ({ ...m, calculation: { kind: "carbonation", amounts: [0, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0], age: 28 } }));
  const path = { id: "path", prepare: async (p) => { await tab(p, "calculation"); await run(p, 10); await result(p, "phases"); } };
  await shoot(page, path);

  // a kinetic run, the four clinker phases by the law of Parrott and Killoh
  const rates = ["C3S", "C2S", "C3A", "C4AF"].map((key) => ({ key, value: { kind: "parrott_killoh_avrami", blaine: 385, water_factor: true, powers_cap: true } }));
  await setModel(page, (m) => ({ ...m, calculation: { kind: "kinetics", horizon: 28, ages: [1, 2, 7, 28], rates, calorimeter: { kind: "isothermal" }, reltol: 1e-6 } }));
  const kinetics = { id: "kinetics", prepare: async (p) => { await tab(p, "calculation"); await run(p, 30); await result(p, "degrees"); } };
  await shoot(page, kinetics);
} finally {
  await browser.close();
}
