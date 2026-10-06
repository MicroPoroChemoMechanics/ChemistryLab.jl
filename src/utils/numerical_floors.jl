# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

# ── Numerical floors, in one place ────────────────────────────────────────────
#
# Each small amount the package works with is defined here once, with the reason
# for its value, and every signature and formula refers to it by name. Written
# as a literal in each place, the default `ϵ = 1e-16` stood in some forty
# signatures, and the one role of it that had to change (the floor of the
# activities) could not be told from the others.

"""
    _AMOUNT_FLOOR

The default of the keyword `ϵ`, `1e-16` mol: the smallest amount the
interior-point back ends are bounded by, the amount of an absent product in a
cold state, and the size of the regularizations that keep a square root or a
ratio finite at zero (`sqrt(I + ϵ)`). It is no longer the floor of the
activities; see [`_ACTIVITY_FLOOR`](@ref).
"""
const _AMOUNT_FLOOR = 1.0e-16

"""
    _ACTIVITY_FLOOR

The amount, in moles, below which an activity model reads a species at this
amount: `1e-30`. It keeps `log n` finite, and nothing else.

It was `ϵ = 1e-16` until 0.28.0, the floor that also bounds the interior-point
back ends and fills the absent products of a cold state. As a floor on the
activities it is too high for a cement paste. With a few grams of water per 100 g
of binder, H⁺ at pH 13.5 to 14 is 1e-16 mol, at the floor: its activity no longer
followed its amount, a solve left it at 3e-100 mol, and a pH read from the amount
came out 0.09 high.

The other codes floor far lower, and one argument recurs. [Leal2017](@citet)
recall that in an interior-point minimization an unstable species ends at an
amount of the order of the barrier parameter, and that
[LealKulikKosakowski2016](@citet) recommended it below 1e-25 so that such a
species holds less than one molecule (1/N_A = 1.66e-24 mol) in a system of one
mole. GEMS3K takes that same molecule as
the least amount it considers (`lowPosNum`, 1.66e-24 mol) and eliminates a
solution species below `DcMin = 1e-30` mol (`ms_multi.h`). PHREEQC sets a molality
to zero below `MIN_LM = −30` in log (`global_structures.h`). Reaktoro bounds every
amount by `epsilon = 1e-16` mol and takes its barrier parameter from it
(`EquilibriumOptions.hpp`); its code states no reason for that value, and it is a
bound of its solver, as our `ϵ` is. So `ϵ` keeps its other two roles, and the
activity models floor at this value, below one molecule in any system up to a
million moles.
"""
const _ACTIVITY_FLOOR = 1.0e-30

"""
    _activity_floor(p) -> Real

The floor an activity model applies to the amounts of `p`: `p.ϵa`, set by
[`_build_params`](@ref) to `min(ϵ, _ACTIVITY_FLOOR)`, or `p.ϵ` for a parameter
tuple built by hand without it.
"""
_activity_floor(p) = hasproperty(p, :ϵa) ? p.ϵa : p.ϵ

"""
    _CERTIFICATE_FLOOR

The default `floor` of [`optimality_certificate`](@ref), `1e-25` mol: an amount
above it is interior, and held to the stationarity equality; below it a pure phase
is held to its inequality and a member of a present phase to the one-sided form of
the equality (OptimaSolver's `kkt_certificate`, whose default it is too).
"""
const _CERTIFICATE_FLOOR = 1.0e-25

"""
    _LOG_UNDERFLOW

The exponent below which an amount computed as `exp(w)` is taken as the
truncation `exp(-700) ≈ 1e-304`: the bound OptimaSolver's dual Newton clamps its
log amounts to (`W_FLOOR`), and the symmetric bound under which `exp` does not
overflow.
"""
const _LOG_UNDERFLOW = -700.0
