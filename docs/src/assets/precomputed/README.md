# Precomputed trajectories

## If you changed the code, read this first

**The documentation build will refuse to run** — not warn, refuse — as soon as
these files are older than the code that produces them. That is deliberate: a
cache that can go stale silently is not a cache, it is a false claim, and the
site would keep showing the trajectory of a solver that no longer exists.

### What makes them stale

Two things, both checked by the guard at the top of `docs/make.jl`:

1. **A commit touching the code they depend on** — `src/`,
   `scripts/ionic_hydration.jl`, `scripts/hydration_calibration.jl`, or
   `scripts/precompute_docs.jl` itself. Documentation and tests are excluded:
   they cannot change a trajectory.
2. **A different resolved version of `OptimaSolver`.** The solver is source too,
   and it does not live in this repository, so no commit here can speak for it.
   This is not hypothetical: 0.5.3 changed what the certificate reports and 0.5.4
   changed it back, and `docs/Manifest.toml` once pinned 0.5.1 for months without
   anyone noticing.

### The procedure, in order

```bash
# 1. Commit your code change FIRST.
git add src/... && git commit -m "..."

# 2. Only then regenerate. The files record the commit they ran at, so a run
#    started before the commit records the OLD one and the guard still refuses.
julia --project=docs scripts/precompute_docs.jl

# 3. Check the headers actually moved.
head -8 docs/src/assets/precomputed/ionic_opc_phases.csv

# 4. Commit the regenerated files, ideally in the very next commit.
git add docs/src/assets/precomputed/ && git commit -m "docs: regenerate ..."

# 5. Now the build runs.
julia --project=docs docs/make.jl
```

### What it costs, and how not to hurt the machine

About **35 to 40 minutes** on two cores: two coupled hydration runs to 28 days
(225 s and 182 s of integration, plus the certified replay of 80 instants each),
then five coupled forward solves for the calibration pages at roughly 210 s
apiece, then three more for the clinker sensitivity.

Run it as **one heavy process at a time** — it is a long single-threaded job, and
starting a documentation build beside it helps nobody. `--only=ionic` and
`--only=calibration` exist so that a change affecting one group does not pay for
the other.

### Two traps, both met in practice

- **An interrupted run leaves a partial file.** The script writes each CSV as it
  goes, so a run killed in the middle of the sensitivity loop leaves
  `calibration_sensitivity.csv` with one row of three. The guard cannot see this:
  the header is already correct. **Check the row counts** against the table
  below, or simply re-run the group.
- **Regenerating before committing** records the previous commit in the header,
  and the guard refuses with a message that looks like the regeneration did not
  happen. It did; it was just early.

### If you are sure the change cannot move a trajectory

Regenerate anyway. A renamed variable cannot change a number, but proving that
costs more than the run does, and a stored result that no longer matches its
source is worse than a slow build. The guard is deliberately not configurable.

---

Every file here is **produced by a script in this repository**, not transcribed
from anywhere, and can be reproduced:

```bash
julia --project=docs scripts/precompute_docs.jl                    # all of it
julia --project=docs scripts/precompute_docs.jl --only=ionic       # one group
julia --project=docs scripts/precompute_docs.jl --only=calibration
julia --project=docs scripts/precompute_docs.jl --quick            # plumbing check only
```

## Why these are not computed at build time

A coupled hydration run costs one Gibbs minimization per accepted ODE step, and
replaying it certifies one equilibrium per reported instant. Measured on two
cores: **224 s** for one `run_ionic_hydration` to 28 days, **364 s** for one
coupled forward solve. The documentation calls for a dozen of those, which is
most of a three-hour build — against about ten minutes for everything else on
the site put together.

So they are computed once and the pages read the result.

**Nothing is coarsened to make that possible**, and the exchange runs the other
way: because the cost is paid once rather than at every build, the trajectories
are reported on **80** log-spaced instants instead of the 40 a build could
afford.

## What a file carries

Each `.csv` opens with a `#` block recording

- the package version and the **commit** it ran at;
- when it was generated;
- the sampling, the water/binder ratio, the temperature and the clinker
  composition;
- and **how many of its instants were proved optimal** against the KKT
  conditions, rather than merely converged.

That last line is the one to read before trusting a figure. A replayed instant
that fails to certify falls back to the interior-point composition, which is not
a proved optimum — on these runs it is a small minority, and the count says
exactly how small.

## The files

| file | what it is |
|:--|:--|
| `ionic_opc_phases.csv` | CEM I, w/b 0.50, 3.5 % limestone — volume fraction of each phase family, pore-solution pH, porosity |
| `ionic_opc_heat.csv` | the same run's isothermal calorimetry at 20 °C |
| `ionic_nolimestone_phases.csv` | the same paste with the limestone removed |
| `ionic_nolimestone_heat.csv` | its calorimetry |
| `calibration_target.csv` | measured heat against the published and the calibrated Parrott-Killoh parameters, on the calibration record |
| `calibration_holdout.csv` | the same on the holdout record, never fitted |
| `calibration_sensitivity.csv` | the fit with the alite content moved ±20 % |

The measured calorimetry in the last three comes from the CC-BY-4.0 deposit of
Šmilauer and Reiterman, [10.5281/zenodo.15212785](https://doi.org/10.5281/zenodo.15212785),
vendored under `data/experimental/`; only the computed columns are produced here.

## Reading them

`scripts/precomputed.jl` is the reader the pages use — a few lines of plain
Julia over `DelimitedFiles`, deliberately small enough to read:

```julia
include(joinpath(pkgdir(ChemistryLab), "scripts", "precomputed.jl"))
tbl = read_precomputed("ionic_opc_phases")
tbl.columns["pore_pH"]      # a column
tbl.provenance              # the header block, so a figure can say where it came from
```
