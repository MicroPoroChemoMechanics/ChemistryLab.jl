# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

using SpecialFunctions: erf, erfc, erfinv

# ── size_classes.jl ───────────────────────────────────────────────────────────
#
# A ground mineral dissolves from the surface of its grains, and the fines carry
# most of that surface. Its particle size distribution is cut into classes of
# equal volume, each a copy of the species with its own amount and its own
# shrinking-core area, ready to receive one rate law.

"""
    AbstractPSD

A particle size distribution **by volume**: the fraction `F(d)` of the volume of
the powder held by particles of diameter below `d`. Its subtypes are
[`LogNormalPSD`](@ref) and [`TabulatedPSD`](@ref), and [`size_classes`](@ref)
cuts one into classes.
"""
abstract type AbstractPSD end

"""
    LogNormalPSD(; d50, d90) -> LogNormalPSD
    LogNormalPSD(; d50, σ) -> LogNormalPSD

A size distribution by volume whose diameter is lognormal,

```math
F(d) = \\Phi\\!\\left(\\frac{\\ln d - \\ln d_{50}}{\\sigma}\\right),
\\qquad
\\sigma = \\frac{\\ln(d_{90}/d_{50})}{\\Phi^{-1}(0.9)} ,
```

with `Φ` the standard normal distribution function. `d50` and `d90` are
lengths, a `Quantity` or a plain `Real` in m. Differentiable in both, and
generic in their number type.
"""
struct LogNormalPSD{T <: Real} <: AbstractPSD
    d50::T          # m
    σ::T
    function LogNormalPSD{T}(d50::Real, σ::Real) where {T <: Real}
        d50 > 0 || throw(ArgumentError("LogNormalPSD: d50 must be positive; got $d50 m."))
        σ > 0 || throw(ArgumentError("LogNormalPSD: σ must be positive; got $σ (d90 must exceed d50)."))
        return new{T}(d50, σ)
    end
end

function LogNormalPSD(; d50, d90 = nothing, σ = nothing)
    (d90 === nothing) == (σ === nothing) &&
        throw(ArgumentError("LogNormalPSD: give either `d90` or `σ`, one of the two."))
    m50 = _area_si(us"m", d50, "LogNormalPSD d50")
    s = σ === nothing ? log(_area_si(us"m", d90, "LogNormalPSD d90") / m50) / _probit(0.9) : float(σ)
    T = promote_type(typeof(m50), typeof(s))
    return LogNormalPSD{T}(m50, s)
end

"""
    TabulatedPSD(d, F) -> TabulatedPSD

A measured size distribution by volume: diameters `d` in ascending order (a
vector of lengths, `Quantity` or plain `Real` in m) and the cumulative volume
fraction `F` passing at each, from zero at the first diameter to its largest
value at the last, in percent or as a fraction (it is divided by its last
value). Between two points `F` is taken linear in `ln d`, which is how a laser
granulometer reports it.
"""
struct TabulatedPSD{T <: Real} <: AbstractPSD
    lnd::Vector{T}
    F::Vector{T}
end

function TabulatedPSD(d::AbstractVector, F::AbstractVector)
    length(d) == length(F) >= 2 ||
        throw(ArgumentError("TabulatedPSD: `d` and `F` must have the same length, two points at least."))
    m = [_area_si(us"m", x, "TabulatedPSD diameter") for x in d]
    all(>(0), m) && issorted(m; lt = <=) ||
        throw(ArgumentError("TabulatedPSD: the diameters must be positive and strictly increasing."))
    issorted(F) || throw(ArgumentError("TabulatedPSD: the cumulative fraction must not decrease."))
    last(F) > 0 || throw(ArgumentError("TabulatedPSD: the cumulative fraction never rises above zero."))
    iszero(first(F)) || throw(
        ArgumentError(
            "TabulatedPSD: the cumulative fraction must start at zero, at the first diameter " *
                "(got $(first(F))); the volume below it would be lost.",
        )
    )
    f = F ./ last(F)
    T = promote_type(eltype(m), eltype(f))
    return TabulatedPSD{T}(log.(T.(m)), T.(f))
end

# The standard normal distribution and its inverse.
_probit(p::Real) = sqrt(oftype(float(p), 2)) * erfinv(2p - 1)

# ∫ dF over the quantiles z ∈ (a, b) of the standard normal: Φ(b) − Φ(a),
# through erfc on the side where it does not cancel. Not symbolic: a branch on
# the side. `√2` is a plain number: a dual one would put −∞ × 0 into the
# derivative of a quantile at infinity.
function _normal_mass(a::Real, b::Real)
    r = sqrt(2 * one(float(_primal(a))))
    _primal(a) >= 0 && return (erfc(a / r) - erfc(b / r)) / 2
    _primal(b) <= 0 && return (erfc(-b / r) - erfc(-a / r)) / 2
    return (erf(b / r) - erf(a / r)) / 2
end

"""
    inverse_diameter_moment(psd, Fa, Fb) -> Real

`∫ (1/d) dF` over the particles between the volume fractions `Fa` and `Fb`, in
1/m: times `6/ρ`, the specific surface those particles add as spheres. Exact
for both distributions — a closed form for [`LogNormalPSD`](@ref),

```math
\\int_{F_a}^{F_b} \\frac{\\mathrm dF}{d}
  = e^{-\\ln d_{50} + \\sigma^2/2}\\,\\bigl[\\Phi(z_b + \\sigma) - \\Phi(z_a + \\sigma)\\bigr],
\\qquad z = \\Phi^{-1}(F) ,
```

and, for [`TabulatedPSD`](@ref), the integral of `e^{-u}` over `F` linear in
`u = ln d` between points — so that the classes of a distribution add up to its
whole surface to rounding.
"""
function inverse_diameter_moment(psd::LogNormalPSD, Fa::Real, Fb::Real)
    za, zb = _probit(Fa), _probit(Fb)
    return exp(-log(psd.d50) + psd.σ^2 / 2) * _normal_mass(za + psd.σ, zb + psd.σ)
end

function inverse_diameter_moment(psd::TabulatedPSD, Fa::Real, Fb::Real)
    acc = zero(promote_type(eltype(psd.F), typeof(float(Fa)), typeof(float(Fb))))
    for i in 1:(length(psd.F) - 1)
        F1, F2 = psd.F[i], psd.F[i + 1]
        lo, hi = max(F1, Fa), min(F2, Fb)
        (F2 > F1 && hi > lo) || continue
        u1, u2 = psd.lnd[i], psd.lnd[i + 1]
        slope = (u2 - u1) / (F2 - F1)
        ua, ub = u1 + slope * (lo - F1), u1 + slope * (hi - F1)
        # ∫ e^{−u} dF with dF = du/slope.
        acc += (exp(-ua) - exp(-ub)) / slope
    end
    return acc
end

"""
    SizeClass

One class of a ground mineral, as [`size_classes`](@ref) returns it.

# Fields

  - `species`: the species of the class — the mineral itself for the first,
    a copy under its own symbol ([`with_symbol`](@ref)) for the others;
  - `amount`: its initial amount, mol;
  - `surface`: its area model, `ShrinkingCoreArea(GeometricSurfaceArea(s))`;
  - `diameter`: its Sauter mean diameter, the diameter of the sphere with the
    class's ratio of volume to surface, m;
  - `fraction`: the fraction of the volume of the powder it holds.
"""
struct SizeClass{S <: AbstractSpecies, T <: Real, A <: AbstractSurfaceModel}
    species::S
    amount::T
    surface::A
    diameter::T
    fraction::T
end

"""
    size_classes(species, psd::AbstractPSD, K; amount | mass, density = nothing,
                 T = 298.15, P = 1e5, shape = 1, names = k -> …) -> Vector{SizeClass}

Cut a ground mineral into `K` classes of equal volume of its size distribution
`psd`, each ready to receive the same rate law.

Class `k` holds the volume fractions `((k−1)/K, k/K)`, a `K`-th of the powder,
and the particles in it, taken as spheres scaled by `shape`, have the specific
surface

```math
s_k = \\text{shape}\\;\\frac{6}{\\rho}\\;K \\int_{(k-1)/K}^{k/K} \\frac{\\mathrm dF}{d},
\\qquad d_k = \\frac{6\\,\\text{shape}}{\\rho\\, s_k} ,
```

`d_k` the Sauter diameter of the class ([`inverse_diameter_moment`](@ref)). The
classes add up to the amount given and to the specific surface of the whole
distribution, both to rounding; one class is the mineral itself with that whole
surface. Each class is a [`SizeClass`](@ref): its species, its amount, and the
area model `ShrinkingCoreArea(GeometricSurfaceArea(s_k))`, a grain shrinking
from the outside.

# Keywords

  - `amount` or `mass`: the whole powder, in mol or in kg (a `Quantity` or a
    plain `Real` in SI);
  - `density`: in kg/m³; by default the molar mass over the molar volume of the
    species, `M/V⁰(T, P)`, which it must then carry;
  - `shape`: the ratio of the particles' surface to that of spheres of the same
    volume, one by default;
  - `names`: the symbol of class `k`, by default the mineral's own for the first
    and `"\$(symbol)_\$k"` for the others.

# Example

```julia
classes = size_classes(subs["C3S"], LogNormalPSD(; d50 = 10u"μm", d90 = 30u"μm"), 4;
                       mass = 1u"kg")
species = vcat(aqueous, [c.species for c in classes])
# one reaction per class, the same law on each class's own surface:
laws = [transition_state(mechs, cs, dissolution(c.species), c.surface) for c in classes]
```

Differentiable in the parameters of the distribution and in the amount.
"""
function size_classes(
        species::AbstractSpecies, psd::AbstractPSD, K::Integer;
        amount = nothing, mass = nothing, density = nothing,
        T::Real = T_STANDARD, P::Real = P_STANDARD, shape::Real = 1,
        names = k -> k == 1 ? symbol(species) : "$(symbol(species))_$k",
    )
    K >= 1 || throw(ArgumentError("size_classes: K must be at least one; got $K."))
    (amount === nothing) == (mass === nothing) &&
        throw(ArgumentError("size_classes: give either the `amount` or the `mass` of the powder."))
    shape > 0 || throw(ArgumentError("size_classes: `shape` must be positive; got $shape."))
    M = _molar_mass_si(species)
    n = amount === nothing ? _area_si(us"kg", mass, "size_classes mass") / M :
        _area_si(us"mol", amount, "size_classes amount")
    ρ = density === nothing ? _species_density(species, T, P) :
        _area_si(us"kg/m^3", density, "size_classes density")
    bounds = range(0, 1; length = K + 1)
    return map(1:K) do k
        s = shape * 6 / ρ * K * inverse_diameter_moment(psd, bounds[k], bounds[k + 1])
        sym = names(k)
        sp = sym == symbol(species) ? species : with_symbol(species, sym)
        R = promote_type(typeof(n), typeof(s))
        SizeClass(sp, R(n / K), ShrinkingCoreArea(GeometricSurfaceArea(s)), R(6 * shape / (ρ * s)), R(1 / K))
    end
end

# The density of a species from its molar mass and molar volume, kg/m³.
function _species_density(species::AbstractSpecies, T::Real, P::Real)
    haskey(species, :V⁰) || throw(
        ArgumentError(
            "size_classes: species \"$(symbol(species))\" carries no molar volume `V⁰`; " *
                "give the `density` of the powder.",
        )
    )
    V = ustrip(us"m^3/mol", _molar_volume(species)(T = T * u"K", P = P * u"Pa"; unit = true))
    return _molar_mass_si(species) / V
end
