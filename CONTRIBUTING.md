# Contributing to ChemistryLab.jl

## The documentation computes what it shows

Every number in the manual is computed by the build. Nothing is read from a
stored result, and that is a deliberate reversal: the heavy coupled
kinetics/equilibrium trajectories *were* precomputed into CSV files for a while,
because one coupled equilibrium cost 583 ms and the site needs thousands of them.

Two solver fixes later it costs 17 ms, and computing them is affordable again.
The exchange was worth making in both directions, and it is worth knowing why it
came back: a stored result is a claim about code that may since have changed, so
keeping the two in step needed a staleness guard, a documented procedure, and a
list of traps for the ways it could silently go wrong — an interrupted run
leaving a partial file whose header was already correct, a regeneration done
before the commit it was meant to follow. **None of that exists any more.** If
you change the solver, the next build simply shows the new answer.

What this costs is build time, and where it goes:

| | |
|:--|--:|
| one coupled hydration integration to 28 days | a few minutes |
| its certified replay on 40 reported instants | comparable |
| everything else on the site together | about ten minutes |

`scripts/precomputed.jl` holds those runs and **memoizes them per process**.
That matters more than it looks: Documenter runs every `@example` block of the
whole site in one process, so a page asking for a run's phase history and its
calorimetry as two tables gets one integration. The shared process is usually a
hazard — see the plot-font guard below — and here it is the thing that makes
this affordable.

If a change of yours makes the build much slower, the resolution of the reported
trajectories is the lever: `N_INSTANTS` in that script. Say so in the page rather
than quietly coarsening it.

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

## Deliberately deferred

Two pieces of work are identified, designed, and **not** in the current release.
They are recorded here rather than left as folklore, with enough detail to be
picked up without rediscovering the analysis.

### 1. Register the kinetic rate laws symbolically, as the thermodynamic ones are

`THERMO_MODELS` stores a thermodynamic model as a **symbolic expression plus the
unit of every parameter**, from which a compiled function is built. That is worth
having: the formula is inspectable by a user, the units are checked once at
registration rather than trusted at every call, and a new model is a table entry
rather than a hand-written closure.

`KINETICS_RATE_MODELS` has the same shape but holds only `:arrhenius`. The cement
laws — `parrot_killoh_avrami`, `waller`, `blaine_factor`, `humidity_factor` — are
hand-written closures.

**The obstacle, and why it is not the one it looks like.** A thermodynamic
function depends on scalars, `(T, P)`; a kinetic function receives the whole
state by name — `n[mineral_name]`, `lna[cat.species]` — and a named lookup is not
a symbolic variable. But reading the laws shows the lookup is confined to one
step: the state is used only to form the degree of reaction `ξ`, after which

```
r = n_init · A(T) · β_B · β_h · f(ξ, T, t; parameters)
```

and `f` **is** scalar. So the registrable object is the reduced law `f`, with
`:vars => [:T, :ξ, :t]` — `:vars` is already a list and takes more than one — and
the closure keeps the state lookup outside it.

`waller` needs one more decision: it branches on `ξ < 1e-12` for the initial
regime, and a branch does not go into a symbolic expression as it stands.

**Why it was not done in 0.18.0.** These closures are the hot path of every
coupled run, the one this release took from 594 s to 284 s, and `rate_models.jl`
is 1239 lines. The gain is real but presentational — a visible formula and units
checked at registration — not a correction. Refactoring a hot path without a
validation budget is the opposite of careful, whatever the destination.

### 2. A performance audit against Julia's own guidelines

Not a sweep for its own sake: a pass with a measurement behind each change, in
the spirit of the performance tips — no abstract types in struct fields, no
needless captures in closures (`R_gas = R_GAS` was one, and a `const` global is
folded where a captured local can be boxed), no arrays copied where a view or a
preallocated buffer serves, multiple dispatch rather than branching on type,
broadcast rather than temporaries.

Two candidates are already known and were left alone on purpose:

- `_one_speciation` allocates a fresh `Vector{Float64}` per speciation, twice
  (`[ustrip(us"mol", x) for x in eq_result.n]`). About 21 MB over a whole
  trajectory, against solves costing 18 ms each — and removing it means returning
  a shared buffer, which risks aliasing between the two candidates the caller
  compares. Measure before touching.
- `KineticFunc` holds `compiled::F` parametrically, which is right; check the
  same of every struct on the kinetic path before assuming it.

The rule that produced the gains in 0.18.0, and that this audit should keep:
**profile first, and let the measurement name the line**. Four wrong turns were
taken in that release by reasoning about where the time "must" be; the profiler
found it in one pass.
