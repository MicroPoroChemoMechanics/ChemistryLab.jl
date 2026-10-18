# ChemistryLabStudio

ChemistryLab Studio, a local web interface to ChemistryLab for cementitious
materials. Users start it with `using ChemistryLabStudio; clstudio()`; the
user guide is the *Studio* chapter of the ChemistryLab documentation. This
file is for those who work on the studio.

## How it is built

- **The server** is this package: an HTTP.jl server on the loopback interface,
  light, that never loads ChemistryLab. Every `/api` call carries the random
  token of the session and comes from the studio's own page (`Host` and
  `Origin` checked, JSON only); static files are guarded against `..`.
- **The compute process** is a second Julia process (Malt.jl) in the environment
  `@chemistrylab-studio`, which the first launch creates. It includes
  `worker/boot.jl`, which loads ChemistryLab and defines `StudioWorker`: the
  catalog, the run of a script, the extraction of its results. It can be
  interrupted, killed and restarted, and keeps its compilation between runs.
- **The model** (`src/model.jl`) holds the choices of a calculation; the page
  edits it, the server keeps it. `src/codegen.jl` writes the script it stands
  for; the studio runs the text of its sections 1 to 5, so the page shows what
  the script computes. The model is embedded at the end of the script with the
  fingerprint of the text (`src/readback.jl`).
- **Memory**: a run is refused below 1.5 GB available (cgroups read, not only
  libuv's figures), and the compute process is killed if it outgrows what was
  left (Linux).
- **`src/core/`** holds what any studio needs and never refers to ChemistryLab:
  memory, browser, environment, compute process, server, launcher.
- **The page** (`web/`) is plain JavaScript without a build step; the vendored
  libraries are listed in `web/vendor/README.md`. Labels are translated in
  `web/i18n/{en,fr}.json`, which hold the same keys.
- **The library** of published formulations is `library/index.toml`: choices
  only, every number read from ChemistryLab's `data/literature` by the templates.

## Tests

```bash
julia --project=ChemistryLabStudio/test -e 'using Pkg; Pkg.instantiate(); include("ChemistryLabStudio/test/runtests.jl")'
```

The tests of the compute process need an environment holding ChemistryLab,
OptimaSolver, DynamicQuantities and Plots, named by `CLSTUDIO_WORKER_PROJECT`;
without it they are skipped. A run with coverage measures the compute process
too, in the same mode.
