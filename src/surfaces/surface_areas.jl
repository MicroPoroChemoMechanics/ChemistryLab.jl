# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

using DynamicQuantities

# ── What this file is for ─────────────────────────────────────────────────────
#
# One abstraction for "how much area does this solid offer", serving both the
# reactive area of a dissolution rate law and the binding area that carries
# surface sites. The three areas are different *numbers* — a BET area is not the
# area accessible to every solute, and neither is a Blaine fineness — but they
# are the same *kind of object*, and writing them once is what lets a
# measurement method be checked rather than trusted.

"""
    abstract type AbstractSurfaceModel end

How much area a solid offers, in m², as a function of how much of it is left.

Concrete subtypes implement

```julia
total_area(model, n::Real, n₀::Real, M::Real) -> Real
```

where `n` is the current molar amount [mol], `n₀` the initial one [mol], and `M`
the molar mass [kg/mol]. The result is a **total** area in m², never a specific
one.

# Why `n₀` is in the signature from the start

Because an area that cannot see where it started cannot describe a
microstructure that evolves. `n` alone expresses "area proportional to the
remaining mass"; the ratio `n/n₀` is what expresses a shrinking grain
([`ShrinkingCoreArea`](@ref)). Passing both costs nothing and means no rate law
has to change its signature when a richer area model is added later.

[`surface_area`](@ref) is the three-argument form kept for callers that have no
initial amount to give; it is `total_area(model, n, n, M)`.

All methods are AD-compatible: no `Float64` conversion, and the element type of
the result follows its arguments.

See also: [`AbstractSpecificArea`](@ref), [`total_area`](@ref),
[`FixedSurfaceArea`](@ref), [`ShrinkingCoreArea`](@ref).
"""
abstract type AbstractSurfaceModel end

"""
    abstract type AbstractSpecificArea <: AbstractSurfaceModel end

A surface model defined by a **specific** area in m²/kg, together with the
measurement it comes from.

The measurement is carried by the *type*, not by a field, and that is the whole
point: [`area_ratio`](@ref) is defined only between two models of the same type,
so comparing a BET area with a Blaine fineness is a `MethodError`-free, named
refusal rather than a plausible wrong number.

Concrete subtypes implement [`specific_area`](@ref) and [`area_method`](@ref),
and inherit [`total_area`](@ref) as `specific_area(m) * max(n, 0) * M`.

See also: [`BETSurfaceArea`](@ref), [`BlaineSurfaceArea`](@ref),
[`GeometricSurfaceArea`](@ref).
"""
abstract type AbstractSpecificArea <: AbstractSurfaceModel end

# ── Unit handling ─────────────────────────────────────────────────────────────

"""
    _area_si(unit, x, what) -> Real

`x` as a bare number in `unit`, **refusing** a quantity that does not convert.

A plain `Real` is taken to be already in SI, which is the documented contract of
every constructor in this file. A `Quantity` is converted, and one whose
dimension does not match raises an `ArgumentError` naming the offending value.

This is deliberately *not* `safe_ustrip`, whose documented behaviour is to fall
back to stripping the quantity in its native unit when the dimensions differ.
That fallback is right where it is used — a dimensionless placeholder unit in
`ThermoFactory` — and wrong here, where it would silently turn `20_000u"m^2/kg"`
passed as a total area into the number 20 000 m².
"""
@inline _area_si(::UnionAbstractQuantity, x::Real, ::AbstractString) = float(x)

function _area_si(unit::UnionAbstractQuantity, x::UnionAbstractQuantity, what::AbstractString)
    # Two attempts rather than a dimension comparison: comparing
    # `SymbolicDimensions` with `Dimensions` raises inside DynamicQuantities, so
    # the conversion itself is the test. What changes here is the third branch,
    # which refuses instead of stripping.
    try
        return ustrip(unit, x)
    catch
        try
            return ustrip(uconvert(unit, x))
        catch
            throw(
                ArgumentError(
                    "$what expects a value convertible to $unit; got $x. " *
                        "A plain number is accepted and taken to be in $unit already.",
                )
            )
        end
    end
end

# ── FixedSurfaceArea ──────────────────────────────────────────────────────────

"""
    struct FixedSurfaceArea{T<:Real} <: AbstractSurfaceModel

A total reactive area in m², independent of how much solid is left.

Suitable for a short run, or when the area is controlled externally — a measured
powder in a reactor whose area is not expected to move over the experiment.

# Fields

  - `A`: total area [m²].

# Examples

```jldoctest
julia> FixedSurfaceArea(0.5)
FixedSurfaceArea{Float64}(0.5)

julia> FixedSurfaceArea(500.0u"cm^2")
FixedSurfaceArea{Float64}(0.05)
```

See also: [`BETSurfaceArea`](@ref), [`total_area`](@ref).
"""

# ── Why every struct here carries an inner constructor ────────────────────────
#
# Julia generates `Foo(x::T) where {T}` for a parametric struct, and that
# generated method is *more specific* than an outer `Foo(x)` written with an
# untyped argument. So the unit conversion and the refusal below were dead code
# for every plain `Real`: `BETSurfaceArea(90)` went straight to the generated
# constructor and came back as `BETSurfaceArea{Int64}`, unvalidated. Defining an
# inner constructor suppresses the generated ones, which makes the outer
# constructor the only way in — the point being that a check placed outside a
# constructor is a check that can be walked around.

struct FixedSurfaceArea{T <: Real} <: AbstractSurfaceModel
    A::T
    FixedSurfaceArea{T}(A::Real) where {T <: Real} = new{T}(convert(T, A))
end

"""
    FixedSurfaceArea(A) -> FixedSurfaceArea

Build a [`FixedSurfaceArea`](@ref) from a total area.

`A` is a plain `Real` in m², or a `Quantity` convertible to m². Unlike the
previous implementation the element type is **kept**, so a `ForwardDiff.Dual`
area differentiates through.
"""
function FixedSurfaceArea(A)
    a = _area_si(us"m^2", A, "FixedSurfaceArea")
    return FixedSurfaceArea{typeof(a)}(a)
end

"""
    total_area(model::FixedSurfaceArea, n, n₀, M) -> Real

The stored area, whatever the amounts. AD-compatible.
"""
total_area(model::FixedSurfaceArea, ::Real, ::Real, ::Real) = model.A

"""
    area_method(model) -> Symbol

The measurement a surface model's number comes from, for error messages.

`:fixed` for a [`FixedSurfaceArea`](@ref), `:BET`, `:Blaine` or `:geometric` for
the specific-area models. It exists so a refusal can name what it refused; the
*logic* that keeps two measurements apart is dispatch, not this symbol.
"""
area_method(::FixedSurfaceArea) = :fixed

# ── Specific-area models ──────────────────────────────────────────────────────

"""
    specific_area(model::AbstractSpecificArea) -> Real

The specific area of `model` in m²/kg.

Not defined for [`FixedSurfaceArea`](@ref), which carries a total area and no
mass to divide it by; asking is a `MethodError`, which is the honest answer.
"""
function specific_area end

"""
    struct BETSurfaceArea{T<:Real} <: AbstractSpecificArea

A specific area in m²/kg from a **BET** (Brunauer-Emmett-Teller) gas-adsorption
measurement.

```
A = A_specific × n × M        [m²]
```

the standard route in reactive-transport modeling (Palandri & Kharaka 2004).
Using it presumes the area accessible to the dissolving solute is the one the
adsorbed gas saw, which is an assumption about the pore structure, not a
measurement of it.

# Fields

  - `A_specific`: specific BET area [m²/kg].

# Examples

```jldoctest
julia> BETSurfaceArea(90.0)
BETSurfaceArea{Float64}(90.0)

julia> BETSurfaceArea(0.09u"m^2/g")
BETSurfaceArea{Float64}(90.0)
```

See also: [`BlaineSurfaceArea`](@ref), [`area_ratio`](@ref).
"""
struct BETSurfaceArea{T <: Real} <: AbstractSpecificArea
    A_specific::T
    BETSurfaceArea{T}(A_specific::Real) where {T <: Real} = new{T}(convert(T, A_specific))
end

"""
    struct BlaineSurfaceArea{T<:Real} <: AbstractSpecificArea

A specific area in m²/kg from a **Blaine** air-permeability fineness.

# Why this is a separate type from [`BETSurfaceArea`](@ref)

Because the two measure different things, and the package already says so in
prose three times over. Silica fume is about 20 000 m²/kg by BET, while the
effective *Blaine* fineness recommended for its Waller kinetics is about
2 000 m²/kg — a factor of ten. Passing the first where the second is expected
produces a hydration rate ten times too fast, with nothing to catch it.

Making the measurement a type is what turns that warning into a refusal:
[`area_ratio`](@ref) between a Blaine and a BET area raises, because no method
covers the pair.

# Fields

  - `A_specific`: Blaine fineness [m²/kg].

# Examples

```jldoctest
julia> BlaineSurfaceArea(385.0)
BlaineSurfaceArea{Float64}(385.0)

julia> BlaineSurfaceArea(380u"m^2/kg")
BlaineSurfaceArea{Float64}(380.0)
```

See also: [`blaine_factor`](@ref), [`area_ratio`](@ref).
"""
struct BlaineSurfaceArea{T <: Real} <: AbstractSpecificArea
    A_specific::T
    BlaineSurfaceArea{T}(A_specific::Real) where {T <: Real} = new{T}(convert(T, A_specific))
end

"""
    struct GeometricSurfaceArea{T<:Real} <: AbstractSpecificArea

A specific area in m²/kg computed from a geometric idealization — a particle
size distribution and an assumed shape — rather than measured by adsorption or
permeability.

Kept distinct from the two measured kinds for the same reason they are distinct
from each other: a geometric area is typically orders of magnitude below a BET
area on the same powder, because it counts no internal roughness.

# Fields

  - `A_specific`: geometric specific area [m²/kg].
"""
struct GeometricSurfaceArea{T <: Real} <: AbstractSpecificArea
    A_specific::T
    GeometricSurfaceArea{T}(A_specific::Real) where {T <: Real} = new{T}(convert(T, A_specific))
end

for (S, method, label) in (
        (:BETSurfaceArea, :(:BET), "BETSurfaceArea"),
        (:BlaineSurfaceArea, :(:Blaine), "BlaineSurfaceArea"),
        (:GeometricSurfaceArea, :(:geometric), "GeometricSurfaceArea"),
    )
    @eval begin
        function $S(A_specific)
            a = _area_si(us"m^2/kg", A_specific, $label)
            return $S{typeof(a)}(a)
        end
        area_method(::$S) = $method
        # Same measurement, so the ratio is meaningful. Defined per concrete
        # type rather than once on the abstract one, so that adding a fourth
        # measurement forces a decision about it instead of inheriting one.
        area_ratio(model::$S, reference::$S) =
            specific_area(model) / specific_area(reference)
    end
end

specific_area(model::AbstractSpecificArea) = model.A_specific

"""
    total_area(model::AbstractSpecificArea, n, n₀, M) -> Real

`specific_area(model) * max(n, 0) * M` [m²] — an area proportional to the mass
still present, which is the same as holding the specific area constant.

Clamped at zero so an amount that an ODE step has pushed slightly negative
returns no area rather than a negative one. The initial amount is unused here;
it matters only to a model that compares the two, such as
[`ShrinkingCoreArea`](@ref). AD-compatible.
"""
function total_area(model::AbstractSpecificArea, n::Real, ::Real, molar_mass::Real)
    return specific_area(model) * max(n, zero(n)) * molar_mass
end

# ── ShrinkingCoreArea ─────────────────────────────────────────────────────────

"""
    SHRINK_FLOOR

Fractional floor regularizing the power law of [`ShrinkingCoreArea`](@ref)
near exhaustion, `1e-8`.

# What it is for, and why a floor rather than a branch

`f^p` with a fractional `p` has an infinite derivative at `f = 0`, so a solid
that an integrator drives to exhaustion produces an infinite Jacobian entry and
an automatic derivative that is `Inf` or `NaN`. The regularized form used here,

```math
g(f) = \\frac{f\\,(f + f_c)^{p-1}}{(1 + f_c)^{p-1}}
```

is exactly zero at `f = 0`, exactly one at `f = 1`, has a finite derivative
everywhere, reduces to `g(f) = f` when `p = 1`, and differs from `f^p` by a
relative `f_c/(3f)` at `p = 2/3` — about 3e-9 at half depletion. It is written
without a comparison so that nothing branches on a value that may be a `Dual`.
"""
const SHRINK_FLOOR = 1.0e-8

"""
    struct ShrinkingCoreArea{A<:AbstractSpecificArea, T<:Real} <: AbstractSurfaceModel

An area that follows a population of grains being consumed, rather than the mass
that remains.

```math
\\mathcal{A}(n) = \\mathcal{A}_0 \\, g\\!\\left(\\frac{n}{n_0}\\right),
\\qquad g(f) \\simeq f^{\\,p}
```

with `𝒜₀` the initial area computed from `initial`, and `g` the regularized
power law described under [`SHRINK_FLOOR`](@ref).

# Choosing the exponent

`p = 2/3` is the geometric answer for spheres consumed from the outside: area
goes as the square of a radius, amount as its cube. `p = 1` recovers an area
proportional to the remaining mass, that is, exactly what a plain
[`BETSurfaceArea`](@ref) gives — so the two are one family, and the exponent is
the modeling choice that separates them. Values in between are common fits and
should be labeled as fits.

# What this does not do

It does not describe an area **blocked** by something growing on it — a hydrate
layer covering a clinker grain reduces the accessible area without consuming the
grain, so it depends on the whole assemblage and not on `n` alone. That is a
different model, and the precedent for writing it is [`PoreHumidity`](@ref),
which already reads the current composition to produce a geometric quantity.

# Fields

  - `initial`: the specific-area model giving the area at `n = n₀`, carrying its
    own measurement method.
  - `exponent`: `p`, dimensionless.

# Examples

```julia
ShrinkingCoreArea(BETSurfaceArea(90.0), 2 // 3)
ShrinkingCoreArea(GeometricSurfaceArea(0.4u"m^2/g"))      # p = 2/3 by default
```

See also: [`total_area`](@ref), [`SHRINK_FLOOR`](@ref).
"""
struct ShrinkingCoreArea{A <: AbstractSpecificArea, T <: Real} <: AbstractSurfaceModel
    initial::A
    exponent::T
    function ShrinkingCoreArea{A, T}(
            initial::AbstractSpecificArea, exponent::Real
        ) where {A <: AbstractSpecificArea, T <: Real}
        return new{A, T}(initial, convert(T, exponent))
    end
end

"""
    ShrinkingCoreArea(initial; exponent = 2//3) -> ShrinkingCoreArea

Build a [`ShrinkingCoreArea`](@ref). The default exponent is the spherical
`2/3`; state it explicitly whenever it is a fit rather than a geometry.
"""
function ShrinkingCoreArea(initial::AbstractSpecificArea; exponent = 2 // 3)
    p = float(exponent)
    return ShrinkingCoreArea{typeof(initial), typeof(p)}(initial, p)
end

# The positional form the docstring uses. It exists explicitly because defining
# an inner constructor suppressed the generated one, and a form that a docstring
# shows had better be a form that runs.
ShrinkingCoreArea(initial::AbstractSpecificArea, exponent::Real) =
    ShrinkingCoreArea(initial; exponent = exponent)

area_method(model::ShrinkingCoreArea) = area_method(model.initial)

"""
    _shrink_fraction(f, p) -> Real

The regularized power law `g` of [`SHRINK_FLOOR`](@ref), evaluated on an already
clamped fraction `f ≥ 0`. Exactly `f` when `p == 1`.
"""
@inline function _shrink_fraction(f::Real, p::Real)
    fc = oneunit(f) * SHRINK_FLOOR
    return f * (f + fc)^(p - 1) / (oneunit(f) + fc)^(p - 1)
end

"""
    total_area(model::ShrinkingCoreArea, n, n₀, M) -> Real

The initial area scaled by the regularized power law of `n/n₀`. AD-compatible,
and finite at `n = 0`.

Raises when `n₀` is zero, because the model is written as a fraction of an
initial amount and there is no sensible answer to "two thirds of nothing".
"""
function total_area(model::ShrinkingCoreArea, n::Real, n₀::Real, molar_mass::Real)
    iszero(n₀) && throw(
        ArgumentError(
            "ShrinkingCoreArea needs a non-zero initial amount n₀; got $n₀. " *
                "Use a specific-area model if the solid is absent at the start.",
        )
    )
    A₀ = total_area(model.initial, n₀, n₀, molar_mass)
    f = max(n, zero(n)) / n₀
    return A₀ * _shrink_fraction(f, model.exponent)
end

# ── The three-argument form, and the ratio ────────────────────────────────────

"""
    surface_area(model, n::Real, molar_mass::Real) -> Real

Total area in m², for a caller that has no initial amount to offer.

It is `total_area(model, n, n, M)` — exact for every model whose area does not
depend on where it started, and therefore **refused** by
[`ShrinkingCoreArea`](@ref), for which it would silently always return the
initial area.

Prefer [`total_area`](@ref) in new code.
"""
surface_area(model::AbstractSurfaceModel, n::Real, molar_mass::Real) =
    total_area(model, n, n, molar_mass)

function surface_area(::ShrinkingCoreArea, ::Real, ::Real)
    throw(
        ArgumentError(
            "ShrinkingCoreArea needs the initial amount: call " *
                "`total_area(model, n, n₀, M)` rather than `surface_area(model, n, M)`, " *
                "which would return the initial area at every step.",
        )
    )
end

"""
    area_ratio(model, reference) -> Real

The dimensionless ratio of two specific areas **of the same measurement kind**.

```jldoctest
julia> area_ratio(BlaineSurfaceArea(462.0), BlaineSurfaceArea(385.0))
1.2

julia> area_ratio(BETSurfaceArea(20_000.0), BlaineSurfaceArea(385.0))
ERROR: ArgumentError: cannot compare a BET area with a Blaine area: they measure different things, and their ratio is not a fineness factor. Convert deliberately, or use two models of the same kind.
```

# Why this refuses rather than divides

Because the two numbers are not commensurable, and the consequence of pretending
otherwise is quantitative: taking silica fume's BET area (about 20 000 m²/kg)
for its effective Blaine fineness (about 2 000 m²/kg) multiplies its hydration
rate by ten. The refusal is the mechanism behind
[`blaine_factor`](@ref) being safe to call.

See also: [`specific_area`](@ref), [`area_method`](@ref).
"""
function area_ratio(model::AbstractSpecificArea, reference::AbstractSpecificArea)
    throw(
        ArgumentError(
            "cannot compare a $(area_method(model)) area with a " *
                "$(area_method(reference)) area: they measure different things, and " *
                "their ratio is not a fineness factor. Convert deliberately, or use " *
                "two models of the same kind.",
        )
    )
end

# ── Surface: the support an area belongs to ───────────────────────────────────

"""
    struct Surface{M<:AbstractSurfaceModel}

A support and the area it offers: a named host solid together with its area
model.

# What this factorizes

The rate factories used to carry a surface model and rediscover the mineral it
belonged to from the reaction, through a helper that fell back to **0.1 kg/mol**
when the species had no molar mass — silently, and wrong by up to an order of
magnitude. A `Surface` names the host once, so the molar mass is looked up once
and its absence is an error that names the species.

The same object is what a family of surface sites will hang from, which is why
it lives here rather than in the kinetics: one support, one area, whether what
happens on it is dissolution or binding.

# Fields

  - `name`: a label for the surface, free-form.
  - `host`: the symbol of the species carrying it, or `nothing` when the support
    is externally prescribed and does not appear in the system.
  - `area`: the [`AbstractSurfaceModel`](@ref).

# Examples

```julia
Surface("calcite", "Cal", BETSurfaceArea(90.0))
Surface("inert sorbent", nothing, FixedSurfaceArea(0.5))
```
"""
struct Surface{M <: AbstractSurfaceModel}
    name::String
    host::Union{Nothing, String}
    area::M
end

"""
    Surface(name, host, area) -> Surface
    Surface(name, area) -> Surface

Build a [`Surface`](@ref). The two-argument form leaves the host unset, for a
support whose amount is prescribed rather than solved for.
"""
Surface(name::AbstractString, area::AbstractSurfaceModel) =
    Surface{typeof(area)}(String(name), nothing, area)

Surface(name::AbstractString, host, area::AbstractSurfaceModel) =
    Surface{typeof(area)}(String(name), host === nothing ? nothing : String(host), area)

area_method(s::Surface) = area_method(s.area)

function Base.show(io::IO, s::Surface)
    host = s.host === nothing ? "prescribed support" : "on $(s.host)"
    return print(io, "Surface(\"$(s.name)\", $host, $(nameof(typeof(s.area))))")
end
