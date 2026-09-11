# Cross-code checks

Scripts that run a second code on a problem this package has already solved, and
report how far the two answers are apart. They are **not** part of the test
suite: they need a conda environment that is not a dependency of anything here,
and they are run by hand when a result is worth a second opinion.

The comparisons that *are* checked on every CI run live in
`test/reference/` and `test/equilibrium_reference.jl`; see
[Validation against Reaktoro](../../docs/src/tutorials/reaktoro_comparison.md).

## Setting one up so that a disagreement means something

Three knobs have to match before a difference can be attributed to anything:

| knob | must match |
|:--|:--|
| thermodynamic database | the same file — both codes read `data/cemdata18-thermofun.json` |
| species list | the same species, and the same grouping into phases |
| activity model | the same convention, including its parameters |

A fourth is specific to a cement, and it is the one most easily missed: **the
element vector**. Handing each code the grams of oxide on the datasheet lets
each convert them with its own atomic-mass table, which separates the two
element vectors by about `1e-4` before any chemistry happens. The Julia side
therefore writes the anhydrous charge *in moles*, and the Python side reads it
and checks the element vectors against each other before equilibrating. On the
CEM I of `cem1_solid_solutions.md` they agree to eight decimal places.

## Running the CEM I solid-solution check

```bash
# once
conda create -n reaktoro-env -c conda-forge reaktoro thermofun

# the Julia side: solves, and writes the charge and the converged composition
julia --project=scripts scripts/crosscheck/cem1_solid_solutions.jl

# the Python side: replays them in Reaktoro and compares
conda run -n reaktoro-env python scripts/crosscheck/cem1_solid_solutions_reaktoro.py
```

Both write into `scripts/crosscheck/out/`, which is not versioned.

## What the CEM I check found, and one thing to know before running it

A full cement assemblage needs several hundred iterations. Reaktoro's default
cap is below that, and a run stopped at the cap returns an intermediate iterate
rather than an answer — with an impossible pH, which is at least easy to notice.
`EquilibriumOptions.optima.maxiters` is raised in the script for that reason.

With all eight solid solutions declared, neither a cold start nor a warm start
from the certified composition converged on the Reaktoro side within 4000
iterations, for any of the floors tried on the empty phases. That is a real
difficulty of the formulation rather than of one code: a solid solution whose
end-members are all at zero has no mole fractions, so its ideal-mixing term is
undefined there and its gradient depends on the direction of approach. Every
code regularizes it somehow. The comparison is therefore made against Reaktoro's
three-phase configuration, which is its documented setup for this system, and
the two answers agree to 0.043 units of pH and 0.08 % of total volume.
