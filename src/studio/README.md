# ChemistryLab Studio

ChemistryLab Studio, a local web interface to ChemistryLab for cementitious
materials, shipped with the package as the submodule `ChemistryLab.Studio`.
Users start it with `using ChemistryLab; clstudio()`; the user guide is the
*Studio* chapter of the ChemistryLab documentation. This file is for those who
work on the studio.

## How it is built

- **The server** runs in the user's session: an HTTP.jl server on the loopback
  interface, which computes nothing itself. Every `/api` call carries the random
  token of the session and comes from the studio's own page (`Host` and
  `Origin` checked, JSON only); static files are guarded against `..`.
- **The compute process** is a second Julia process (Malt.jl) in the environment
  `@chemistrylab-studio`, which the first launch creates with the very
  ChemistryLab running (that version pinned, or the clone developed). It includes
  `worker/boot.jl`, which loads ChemistryLab and defines `StudioWorker`: the
  catalog, the run of a script, the extraction of its results. It can be
  interrupted, killed and restarted, and keeps its compilation between runs.
- **The model** (`model.jl`) holds the choices of a calculation; the page
  edits it, the server keeps it. `codegen.jl` writes the script it stands
  for; the studio runs the text of its sections 1 to 5, so the page shows what
  the script computes. The model is embedded at the end of the script with the
  fingerprint of the text (`readback.jl`).
- **One declaration per choice** (`schema.jl`, `features/`): a family of
  choices (database, system, activity model, material, law of reaction, mix,
  calculation, constraint, rate law, calorimeter) is an abstract type; each
  member declares its fields once (kind, default, unit, bounds, level) and
  writes its own lines of the script. The JSON of the model, its checks, and
  the forms of the page (`web/forms.js`, from `/api/schema`) all come from the
  declarations. Adding a choice is adding a type, its fields, its code and its
  words in `web/i18n/`; the tests refuse a member whose fields do not match its
  struct, a kind without its words, a script that does not parse.
- **Memory**: a run is refused below 1.5 GB available (the cgroups read on
  Linux, the inactive pages counted on macOS), and the compute process is
  killed if it outgrows what was left.
- **Stopping**: Ctrl-C ends a script at once, and the compute process follows
  (its connection closes, which ends Malt's loop); at the REPL the interrupt
  reaches the waiting loop, or is relayed to it by the task it landed in.
- **`core/`** holds what any studio needs and never refers to ChemistryLab:
  memory, browser, environment, compute process, server, launcher.
- **The page** (`web/`) is plain JavaScript without a build step; the vendored
  libraries are listed in `web/vendor/README.md`. Labels are translated in
  `web/i18n/{en,fr}.json`, which hold the same keys.
- **The library** of published formulations is `library/index.toml`: choices
  only, every number read from ChemistryLab's `data/literature` by the templates.

## Screenshots

The screenshots of the documentation and of the README are taken from a running
studio by `tools/screenshots.mjs` (Node and playwright-core), in the light and
the dark theme:

```bash
node src/studio/tools/screenshots.mjs --url 'http://127.0.0.1:8766/?token=…'
```

## Tests

The tests are `test/studio/`, part of the package's suite (`core` part); alone,
from the root of the repository:

```bash
julia --project=. -e 'include("test/studio/runtests.jl")'
```

The tests of the compute process need an environment holding ChemistryLab,
OptimaSolver, DynamicQuantities, OrdinaryDiffEq and Plots, named by
`CLSTUDIO_WORKER_PROJECT`;
without it they are skipped. A run with coverage measures the compute process
too, in the same mode. The studio's coverage is reported by the `Studio`
workflow, on the three systems, and left out of the main one.
