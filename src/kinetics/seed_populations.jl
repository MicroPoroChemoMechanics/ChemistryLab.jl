# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

# ── seed_populations.jl ───────────────────────────────────────────────────────
#
# A phase that precipitates from its own nuclei: a reservoir of seeds, a copy of
# the species fed by a nucleation rate, and growth on the surface those seeds
# offer, their number read from the volume of the reservoir and the size of one
# seed. Two shapes: spheres, and needles that lengthen up to a cap.

"""
    AbstractSeedGeometry

The shape of the seeds of a [`SeededGrowthArea`](@ref): [`SphericalSeeds`](@ref)
or [`NeedleSeeds`](@ref). A geometry gives the volume of one seed and the area
a population grows on.
"""
abstract type AbstractSeedGeometry end

"""
    SphericalSeeds(radius) -> SphericalSeeds

Seeds that are spheres of radius `radius` (a length, `Quantity` or `Real` in m),
growing as `N` equal spheres. With `N = V_res/v₀` seeds, `v₀ = 4πr₀³/3`, holding
the volume `V` of the phase and its seeds,

```math
\\mathcal{A} = (36\\pi N)^{1/3}\\, V^{2/3} = \\frac{3}{r_0}\\, V \\left(\\frac{V_\\text{res}}{V}\\right)^{1/3} ,
```

the area of the seeds themselves while nothing has grown on them.
"""
struct SphericalSeeds{T <: Real} <: AbstractSeedGeometry
    radius::T
    function SphericalSeeds{T}(r::Real) where {T <: Real}
        r > 0 || throw(ArgumentError("SphericalSeeds: the radius must be positive; got $r m."))
        return new{T}(r)
    end
end
SphericalSeeds(radius) = (r = _area_si(us"m", radius, "SphericalSeeds radius"); SphericalSeeds{typeof(r)}(r))

"""
    NeedleSeeds(; length, aspect, max_length, width = 1e-3) -> NeedleSeeds

Seeds that are needles: cylinders of length `l` and radius `aspect·l`, growing
homothetically from the length `length` up to `max_length`, where they meet one
another. With `N = V_res/v₀` seeds, `v₀ = π aspect² l₀³`, holding the volume `V`,
the current length is `l = l₀ (V/V_res)^{1/3}`, and the area the growth law
multiplies is the rate of change of the volume with the length,

```math
\\mathcal{A} = \\frac{3V}{l}\\; c\\!\\left(1 - \\frac{l}{l_\\max}\\right)
             = \\frac{3}{l_0}\\, V \\left(\\frac{V_\\text{res}}{V}\\right)^{1/3}
               c\\!\\left(1 - \\frac{l}{l_\\max}\\right) ,
```

so that a growth law `k(Ω)` per unit of that area is a lengthening rate
`dl/dt = V_m k(Ω)`. `c` stops the growth at `l_max`, continuously with its slope,
across a fraction `width` of the length: `c(z) = z²/(z + width)` above zero and
zero below. Lengths as `Quantity` or `Real` in m.
"""
struct NeedleSeeds{T <: Real} <: AbstractSeedGeometry
    length::T
    aspect::T
    max_length::T
    width::T
    function NeedleSeeds{T}(l0::Real, a::Real, lmax::Real, w::Real) where {T <: Real}
        (l0 > 0 && a > 0 && w > 0) || throw(ArgumentError("NeedleSeeds: length, aspect and width must be positive."))
        lmax > l0 || throw(ArgumentError("NeedleSeeds: max_length must exceed the seed length."))
        return new{T}(l0, a, lmax, w)
    end
end
function NeedleSeeds(; length, aspect, max_length, width = 1.0e-3)
    l0 = _area_si(us"m", length, "NeedleSeeds length")
    lm = _area_si(us"m", max_length, "NeedleSeeds max_length")
    v = promote(l0, float(aspect), lm, float(width))
    return NeedleSeeds{eltype(v)}(v...)
end

"""
    seed_volume(geometry) -> Real

The volume of one seed, m³.
"""
seed_volume(g::SphericalSeeds) = 4 * oftype(float(g.radius), π) * g.radius^3 / 3
seed_volume(g::NeedleSeeds) = oftype(g.length, π) * g.aspect^2 * g.length^3

# The growth area of a population holding the volume V, of which `f = V_res/V`
# is the seeds'. The cube root exactly: seeds are a small part of the volume
# (a floor on `f` misplaced the area by 2e-5 at f = 3e-4, and by a third at
# 1e-8), and the amounts of a run are floored above zero, so its slope stays
# finite there. Zero without seeds, its slope included. Not symbolic: that test.
_seed_cbrt(f) = _primal(f) > 0 ? cbrt(f) : zero(f)
_seed_growth_area(g::SphericalSeeds, V, f) = 3 / g.radius * V * _seed_cbrt(f)
function _seed_growth_area(g::NeedleSeeds, V, f)
    q = _seed_cbrt(f)                                   # (V_res/V)^{1/3} = l₀/l
    _primal(q) > 0 || return zero(V * q)
    l = g.length / q
    return 3 / g.length * V * q * _c1_cut(1 - l / g.max_length, g.width)
end

"""
    SeededGrowthArea(reservoir, geometry; T = 298.15, P = 1e5) -> SeededGrowthArea

The area a phase grows on when it grows on its own seeds: those of `reservoir`,
a copy of the phase ([`with_symbol`](@ref)) that a [`nucleation_rate`](@ref)
feeds, of the shape `geometry`. Give it as the surface of a
[`transition_state`](@ref) law of the phase, typically restricted to growth
([`growth_only`](@ref)):

```julia
seeds = SeededGrowthArea(cs["Portlandite_nuc"], SphericalSeeds(50u"nm"))
grow = growth_only(transition_state([RateMechanism(k, 1, 1)], cs, rxn, seeds), cs, rxn)
```

The volume of the population is `V = V_m (n_res + n)`, the molar volume `V_m`
of the reservoir (the phase's, the reservoir being a copy of it) at `T` and `P`.
The area is that of the geometry ([`SphericalSeeds`](@ref),
[`NeedleSeeds`](@ref)), zero without seeds.
"""
struct SeededGrowthArea{G <: AbstractSeedGeometry, T <: Real} <: AbstractSurfaceModel
    geometry::G
    reservoir::String
    V_m::T
end

function SeededGrowthArea(reservoir::AbstractSpecies, geometry::AbstractSeedGeometry; T::Real = T_STANDARD, P::Real = P_STANDARD)
    haskey(reservoir, :V⁰) || throw(
        ArgumentError("SeededGrowthArea: the reservoir \"$(symbol(reservoir))\" carries no molar volume `V⁰`."),
    )
    V = ustrip(us"m^3/mol", _molar_volume(reservoir)(T = T * u"K", P = P * u"Pa"; unit = true))
    return SeededGrowthArea(geometry, symbol(reservoir), V)
end

total_area(m::SeededGrowthArea, ::Real, ::Real, ::Real) = throw(
    ArgumentError(
        "SeededGrowthArea: the area of a growth on seeds depends on the seeds as well " *
            "as on the phase; it is the surface of a rate law (`transition_state`), " *
            "not of a site budget.",
    )
)

_check_reactive_area(cs::ChemicalSystem, m::SeededGrowthArea) =
    haskey(cs.dict_species, m.reservoir) || throw(
    ArgumentError("SeededGrowthArea: the reservoir \"$(m.reservoir)\" is not a species of the system."),
)

function _reactive_area(m::SeededGrowthArea, T, n, n_initial, key, M)
    nr = max(n[m.reservoir], zero(n[m.reservoir]))
    nb = max(n[key], zero(n[key]))
    V = m.V_m * (nr + nb)
    f = m.V_m * nr / max(V, floatmin(Float64))
    return _seed_growth_area(m.geometry, V, f)
end

# ── Nucleation ────────────────────────────────────────────────────────────────

"""
    AbstractNucleationLaw

The frequency of nucleation as a function of the saturation ratio, per unit of
the extent nuclei form on ([`nucleation_rate`](@ref)):
[`ThresholdNucleation`](@ref) or [`ClassicalNucleation`](@ref).
"""
abstract type AbstractNucleationLaw end

"""
    ThresholdNucleation(J, Ω_c; width = 1e-3) -> ThresholdNucleation

Nucleation above a critical saturation ratio `Ω_c`, at a frequency proportional
to the excess, `J c(Ω/Ω_c − 1)` with `c(z) = z²/(z + width)` above zero and zero
below: nothing forms below `Ω_c`, and the frequency is continuous with its slope
there.
"""
struct ThresholdNucleation{T <: Real} <: AbstractNucleationLaw
    J::T
    Ω_c::T
    width::T
    function ThresholdNucleation{T}(J::Real, Ω_c::Real, width::Real) where {T <: Real}
        (J >= 0 && Ω_c > 0 && width > 0) ||
            throw(ArgumentError("ThresholdNucleation: J ≥ 0, Ω_c > 0 and width > 0 are required."))
        return new{T}(J, Ω_c, width)
    end
end
function ThresholdNucleation(J::Real, Ω_c::Real; width::Real = 1.0e-3)
    v = promote(float(J), float(Ω_c), float(width))
    return ThresholdNucleation{eltype(v)}(v...)
end

"""
    ClassicalNucleation(J₂, B) -> ClassicalNucleation

The frequency of classical nucleation theory, written through its value `J₂`
at `Ω = 2`,

```math
J(\\Omega) = J_2 \\exp\\!\\left(\\frac{B}{\\ln^2 2} - \\frac{B}{\\ln^2\\Omega}\\right)
\\quad (\\Omega > 1), \\qquad J = 0 \\quad (\\Omega \\le 1) ,
```

`B` the barrier, `16πγ³v²/(3(kT)³)` for spherical nuclei. Every derivative of
`J` vanishes as `Ω → 1⁺`, so the frequency switches on smoothly.
"""
struct ClassicalNucleation{T <: Real} <: AbstractNucleationLaw
    J₂::T
    B::T
    function ClassicalNucleation{T}(J₂::Real, B::Real) where {T <: Real}
        (J₂ >= 0 && B > 0) || throw(ArgumentError("ClassicalNucleation: J₂ ≥ 0 and B > 0 are required."))
        return new{T}(J₂, B)
    end
end
function ClassicalNucleation(J₂::Real, B::Real)
    v = promote(float(J₂), float(B))
    return ClassicalNucleation{eltype(v)}(v...)
end

# Not symbolic: both laws switch on the sign of a value.
_nucleation_frequency(law::ThresholdNucleation, lnΩ) = law.J * _c1_cut(exp(lnΩ) / law.Ω_c - 1, law.width)
function _nucleation_frequency(law::ClassicalNucleation, lnΩ)
    x = _primal(lnΩ)
    # Beyond B/ln²Ω ≈ 700 the exponential is below the smallest normal number,
    # and its derivative would be 0 × ∞.
    (x > 0 && law.B / x^2 < 700) || return zero(lnΩ) * zero(law.J₂)
    return law.J₂ * exp(law.B / log(oftype(float(x), 2))^2 - law.B / lnΩ^2)
end

"""
    FixedExtent(x)

Nuclei form at a frequency per unit of a fixed extent `x` — an area in m², a
volume of solution, whatever `J` is given per.
"""
struct FixedExtent{T <: Real}
    value::T
end

"""
    SubstrateAmount(species)

Nuclei form on another solid in proportion to its amount, `J` per mole of it.
"""
struct SubstrateAmount
    species::String
end

"""
    SubstrateSurface(species, area)

Nuclei form on the surface of another solid, `J` per m² of it: the area model
`area` of `species` ([`total_area`](@ref)) at its current and initial amounts.
"""
struct SubstrateSurface{A <: AbstractSurfaceModel}
    species::String
    area::A
end

# The extent at the amounts of the run; `M` the substrate's molar mass when it
# has an area.
_nucleation_extent(e::FixedExtent, n, n0, M) = e.value
_nucleation_extent(e::SubstrateAmount, n, n0, M) = max(n[e.species], zero(n[e.species]))
_nucleation_extent(e::SubstrateSurface, n, n0, M) =
    total_area(e.area, max(n[e.species], zero(n[e.species])), n0[e.species], M)

_substrate_mass(::ChemicalSystem, ::FixedExtent) = 1.0
function _substrate_mass(cs::ChemicalSystem, e::Union{SubstrateAmount, SubstrateSurface})
    sp = get(cs.dict_species, e.species, nothing)
    sp === nothing && throw(ArgumentError("nucleation_rate: the substrate \"$(e.species)\" is not a species of the system."))
    return e isa SubstrateSurface ? _molar_mass_si(sp) : 1.0
end

"""
    nucleation_rate(law, cs, rxn, geometry; on, saturation = nothing,
                    T = 298.15, P = 1e5) -> KineticFunc

The rate at which the reservoir of seeds of a phase fills by nucleation: `rxn`
is the dissolution of the reservoir, a copy of the phase, and the rate is

```math
r = -\\,J(\\Omega)\\; E\\; \\frac{v_0}{V_m}
\\qquad [\\text{mol/s}],
```

negative, the reaction being written as a dissolution: `J` from `law`
([`ThresholdNucleation`](@ref), [`ClassicalNucleation`](@ref)), `E` the extent
nuclei form on (`on`: a [`FixedExtent`](@ref), the amount of another solid,
[`SubstrateAmount`](@ref), or its surface, [`SubstrateSurface`](@ref)), `v₀` the
volume of one seed ([`seed_volume`](@ref)) and `V_m` the molar volume of the
reservoir at `T` and `P`. `Ω` is the saturation ratio of `rxn`, the phase's,
unless `saturation`, a function `(T, P, lna) -> Ω`, says otherwise — the sum of
the ratios of the end-members of a solid solution, for instance.

The growth on these seeds is a separate law, of the phase, on a
[`SeededGrowthArea`](@ref).
"""
function nucleation_rate(
        law::AbstractNucleationLaw, cs::ChemicalSystem, rxn::AbstractReaction,
        geometry::AbstractSeedGeometry;
        on, saturation = nothing, T::Real = T_STANDARD, P::Real = P_STANDARD,
    )
    i = _find_mineral_idx(cs, rxn)
    i === nothing && throw(ArgumentError("nucleation_rate: no reservoir found in reaction \"$(rxn.symbol)\"."))
    reservoir = cs.species[i]
    haskey(reservoir, :V⁰) || throw(
        ArgumentError("nucleation_rate: the reservoir \"$(symbol(reservoir))\" carries no molar volume `V⁰`."),
    )
    V_m = ustrip(us"m^3/mol", _molar_volume(reservoir)(T = T * u"K", P = P * u"Pa"; unit = true))
    per_mole = seed_volume(geometry) / V_m
    M = _substrate_mass(cs, on)
    lnΩ = saturation === nothing ? _ReactionLogSaturation(_saturation_terms(cs, rxn, "nucleation_rate")) :
        _GivenLogSaturation(saturation)
    return _nucleation_law(law, lnΩ, on, M, per_mole)
end

# ln Ω of the reaction itself, or of a saturation ratio the caller gives.
struct _ReactionLogSaturation{V}
    terms::V
end
(f::_ReactionLogSaturation)(T, P, lna) = _log_saturation(f.terms, T, P, lna)
struct _GivenLogSaturation{F}
    Ω::F
end
(f::_GivenLogSaturation)(T, P, lna) = log(f.Ω(T, P, lna))

# The law behind a function barrier, specialized on each of its parts.
function _nucleation_law(law, lnΩ, on, M, per_mole)
    rate = (T, P, t, n, lna, n₀) -> -_nucleation_frequency(law, lnΩ(T, P, lna)) *
        _nucleation_extent(on, n, n₀, M) * per_mole
    return KineticFunc(rate, (T = T_STANDARD_Q, P = P_STANDARD_Q), u"mol/s")
end
