# Precomputed trajectories

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
