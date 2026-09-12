# Contributing to ChemistryLab.jl

## The one thing that will bite you first

**If you change anything under `src/`, the documentation build will refuse to
run until you regenerate the precomputed trajectories.**

It refuses rather than warns, on purpose. Four of the documentation's pages read
coupled kinetics/equilibrium runs that are computed **once** rather than at every
build — a single coupled forward solve costs about 210 s and the site calls for a
dozen, which is the difference between a 23-minute build and a three-hour one.
Those stored results are a cache, and a cache that can go stale without saying so
is not a cache but a false claim: the site would keep publishing the trajectory
of a solver that no longer exists.

So each stored file records the commit and the solver version it was produced
with, and `docs/make.jl` compares both against what the build resolves.

### What to do

```bash
# 1. Commit your code change FIRST -- the files record the commit they ran at.
git commit -m "..."

# 2. Regenerate (about 35-40 min; one heavy process at a time).
julia --project=docs scripts/precompute_docs.jl

# 3. Commit the regenerated CSVs.
git add docs/src/assets/precomputed/ && git commit -m "docs: regenerate ..."
```

`--only=ionic` and `--only=calibration` restrict it to one group when only one is
affected. `docs/src/assets/precomputed/README.md` carries the full procedure, the
cost of each part, and the two traps that have actually caught someone: an
interrupted run leaves a **partial** file the guard cannot detect, and
regenerating *before* committing records the previous commit.

**Do not try to decide that your change cannot have moved a trajectory.** Proving
it costs more than the run, and the guard is deliberately not configurable.

## Building the documentation

```bash
julia --project=docs docs/make.jl
```

Roughly 23 minutes. Every `@example` block runs, in **one shared process** — so a
global set on one page (`Plots.default(fontfamily = ...)` is the one that has
caused trouble) applies to every page built afterwards. `docs/make.jl` carries a
guard that refuses a page setting the plot font, and a timer that names any block
taking more than five seconds.

A `draft` build executes nothing and therefore proves nothing. Do not report a
page as working on the strength of one.

## Before you push

The pre-push hook runs these; running them yourself is faster than a rejected
push:

```bash
python3 .github/scripts/check_docrefs.py                       # @ref reachability
python3 .github/scripts/spelling_tool.py check <paths> --convention us
julia --project=@runic -e 'using Runic; exit(Runic.main(["--check", "--diff", "src", "test"]))'
```

## House rules worth knowing

- **US English everywhere** — source comments, docstrings, `docs/` prose, README
  and CHANGELOG entries, commit messages. `spelling_tool.py` scans prose without
  touching identifiers or fenced code. Database symbols are **not** ours to
  correct: CEMDATA18 spells them `monosulphate12`, and that is a key, not a word.
- **Never invent a reference.** Every DOI is resolved against Crossref before it
  is written, and every entry without one is checked against the publisher.
- **Synthetic or assumed data is labeled where it is used**, not in a preamble. A
  composition chosen at the midpoint of an EN 197-1 range is written `ASSUMED` at
  the point of use, because that is where a reader is least able to check it.
- **Molar masses come from the species**, never from a table typed by hand.
- The comparison with GEM-Selektor, Reaktoro and Optima is **constructive**:
  where this package differs, the difficulty is in the formulation, not in the
  other code.

## Tests

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Note that `Pkg.test()` resolves dependencies **from the registry**: a locally
modified dependency is *not* loaded, and a measurement that assumes otherwise is
silently invalid. Check `pkgdir(TheDependency)` before interpreting a result that
depends on one.
