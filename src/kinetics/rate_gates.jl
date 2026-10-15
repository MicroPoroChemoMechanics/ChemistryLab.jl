# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

# ── rate_gates.jl ─────────────────────────────────────────────────────────────
#
# Laws that run one way only, and the smooth cuts they are built from. Both cuts
# below are zero on the forbidden side, exactly, and continuous with their slope
# at the threshold, so a stiff integrator meets no jump in a rate or in its
# Jacobian there.

# z²/(z + δ) above zero and nothing below: a rate proportional to the excess z
# over a threshold, continuous with its slope at z = 0, and within δ of z above
# it. Not symbolic: a branch on the sign of z.
_c1_cut(z::Real, δ::Real) = _primal(z) > 0 ? z^2 / (z + δ) : zero(z) * zero(δ)

# x²/(x² + δ²) above zero and nothing below: a factor that switches a law on
# across a width δ of its driving force, C¹ at x = 0 whatever the law it
# multiplies. Not symbolic: a branch on the sign of x.
_c1_step(x::Real, δ::Real) = _primal(x) > 0 ? x^2 / (x^2 + δ^2) : zero(x) * zero(δ)

# `ln Ω` of a reaction from its terms (`_saturation_terms`), written as a sum so
# that a large supersaturation does not overflow on its way through `exp`.
_log_saturation(terms, T, P, lna) =
    sum(ν * (lna[s] + g(; T = T, P = P, unit = false) / (R_GAS * T)) for (s, ν, g) in terms)

# The two directions a law may be restricted to, and the driving force of each:
# the opposite of the affinity's sign convention, positive where the law may run.
struct _Dissolution end
struct _Growth end
_driving_force(::_Dissolution, lnΩ) = -lnΩ
_driving_force(::_Growth, lnΩ) = lnΩ

# A dissolving phase stops when it is gone: n/(n + ε n₀), one down to a part in
# 10⁸ of the initial amount and zero at zero. A growing one has no such limit.
_depletion_gate(::_Dissolution, n, n₀) =
    (m = max(n, zero(n)); m / (m + SHRINK_FLOOR * max(n₀, m) + floatmin(Float64)))
_depletion_gate(::_Growth, n, n₀) = one(n)

_rate_refs(rate::KineticFunc) = rate.refs
_rate_refs(::Any) = NamedTuple()

"""
    dissolution_only(rate, cs, rxn; width = 1e-3) -> KineticFunc

The rate law `rate`, of the reaction `rxn` written as the dissolution of its
controlled species, restricted to one direction: dissolution only, while the
solution is undersaturated and the phase is present, or growth only, while it is
supersaturated.

```math
r_\\mathrm{diss} = r \\; s(-\\ln\\Omega)\\; \\frac{n}{n + \\varepsilon\\,n_0} ,
\\qquad
r_\\mathrm{growth} = r \\; s(\\ln\\Omega) ,
\\qquad
s(x) = \\begin{cases} \\dfrac{x^2}{x^2 + w^2} & x > 0 \\\\ 0 & x \\le 0 \\end{cases}
```

with `w = width`, `Ω` the saturation ratio of `rxn` ([`saturation_ratio`](@ref)),
`n` and `n₀` the current and initial amounts of the controlled species, and
`ε = SHRINK_FLOOR`. The factor `s` is zero on the forbidden side and continuous
with its slope at `Ω = 1`, whatever `rate` is, so neither the rate nor the
Jacobian of a stiff integrator jumps there. Above a few widths of `ln Ω` the
law is `rate` itself, to within `(w/ln Ω)²`.

This is the law of a phase that dissolves and never reprecipitates — an anhydrous
clinker phase, or a calcium sulfate whose solution stays supersaturated without
it growing back — and, with `growth_only`, of a phase that grows on its seeds
and is never redissolved. The Langmuir factor of an inhibitor is already a
mechanism of [`transition_state`](@ref), [`RateModelInhibitor`](@ref).

# Example

```julia
law = first_order_rate(k, cs, rxn, ShrinkingCoreArea(BlaineSurfaceArea(blaine)))
kr = KineticReaction(cs, rxn, dissolution_only(law, cs, rxn))
```

Not symbolic: the factor switches on the sign of `ln Ω`.
"""
dissolution_only(rate, cs::ChemicalSystem, rxn::AbstractReaction; width::Real = 1.0e-3) =
    _one_way(_Dissolution(), rate, cs, rxn, width, "dissolution_only")

"""
    growth_only(rate, cs, rxn; width = 1e-3) -> KineticFunc

[`dissolution_only`](@ref) in the other direction: `rate` while the solution is
supersaturated in the controlled species of `rxn`, zero otherwise, switched on
across a width `width` of `ln Ω`, with no depletion factor.
"""
growth_only(rate, cs::ChemicalSystem, rxn::AbstractReaction; width::Real = 1.0e-3) =
    _one_way(_Growth(), rate, cs, rxn, width, "growth_only")

function _one_way(dir, rate, cs::ChemicalSystem, rxn::AbstractReaction, width::Real, caller::AbstractString)
    width > 0 || throw(ArgumentError("$caller: `width` must be positive; got $width."))
    terms = _saturation_terms(cs, rxn, caller)
    i = _find_mineral_idx(cs, rxn)
    i === nothing && throw(ArgumentError("$caller: no controlled species found in reaction \"$(rxn.symbol)\"."))
    key = _rate_lookup_key(cs, cs.species[i])
    gated = (T, P, t, n, lna, n₀) -> rate(T, P, t, n, lna, n₀) *
        _c1_step(_driving_force(dir, _log_saturation(terms, T, P, lna)), width) *
        _depletion_gate(dir, n[key], n₀[key])
    return KineticFunc(gated, _rate_refs(rate), u"mol/s")
end
