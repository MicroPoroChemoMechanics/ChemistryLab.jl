# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

using DynamicQuantities
using ForwardDiff
using OrderedCollections

# ── Solid solution activity models ────────────────────────────────────────────

"""
    abstract type AbstractSolidSolutionModel end

Base type for activity models within a solid solution phase.
Each concrete subtype implements `_excess_ln_gamma(model, k, x, T)`.
"""
abstract type AbstractSolidSolutionModel end

"""
    struct IdealSolidSolutionModel <: AbstractSolidSolutionModel

Ideal (Temkin) solid solution mixing: `ln γᵢ = 0`, hence `ln aᵢ = ln xᵢ`.

Valid for any number of end-members.

# Example

```jldoctest
julia> m = IdealSolidSolutionModel()
IdealSolidSolutionModel()
```
"""
struct IdealSolidSolutionModel <: AbstractSolidSolutionModel end

"""
    struct RegularSolutionModel{T<:Real} <: AbstractSolidSolutionModel

Symmetric regular (multi-component Margules) solid solution.

The gap this fills is arity. [`RedlichKisterModel`](@ref) is the general binary
form and is restricted to two end-members, while [`IdealSolidSolutionModel`](@ref)
takes any number but no interaction at all. A C-S-H with six end-members, or the
CNASH and ECSH families of CEMDATA18, had no non-ideal option here.

Excess Gibbs energy, with one symmetric interaction parameter per pair:

```
G^ex = Σ_{i<j} W_ij x_i x_j
ln γ_k = (1/RT) [ Σ_{j≠k} W_kj x_j  −  Σ_{i<j} W_ij x_i x_j ]
```

`W` is in J/mol, symmetric, with a zero diagonal; only the off-diagonal entries
are read and `W[i,j]` is used for the pair `(i,j)`. `W_ij > 0` is repulsive —
it favors unmixing — and `W_ij = 0` recovers ideal mixing.

For two end-members this is exactly `RedlichKisterModel(a0 = W₁₂)`, which the
test suite checks; the point of it is `n > 2`.

AD-compatible: all computations propagate `ForwardDiff.Dual` numbers.

# Examples

```jldoctest
julia> m = RegularSolutionModel([0.0 4000.0; 4000.0 0.0]);

julia> m.W[1, 2]
4000.0
```

# References

  - Guggenheim, E.A. (1937). *Trans. Faraday Soc.* **33**, 151–159.
"""
struct RegularSolutionModel{T <: Real} <: AbstractSolidSolutionModel
    W::Matrix{T}

    # The validation lives in an *inner* constructor on purpose. As an outer
    # method on `AbstractMatrix` it was dead code for the commonest call: the
    # constructor Julia generates from the field declaration,
    # `RegularSolutionModel(::Matrix{T})`, is more specific than
    # `AbstractMatrix`, so a plain `Matrix` argument went straight into the
    # struct and a non-square or asymmetric `W` was accepted in silence.
    # Declaring an inner constructor suppresses the generated ones, which leaves
    # this the only way in.
    function RegularSolutionModel{T}(W::AbstractMatrix) where {T <: Real}
        n, m = size(W)
        n == m || throw(
            ArgumentError("RegularSolutionModel: W must be square, got $(size(W)).")
        )
        for i in 1:n, j in (i + 1):n
            isapprox(W[i, j], W[j, i]; rtol = 1.0e-12) || throw(
                ArgumentError(
                    "RegularSolutionModel: W must be symmetric; W[$i,$j] = $(W[i, j]) " *
                        "but W[$j,$i] = $(W[j, i]).",
                )
            )
        end
        return new{T}(Matrix{T}(W))
    end
end

"""
    RegularSolutionModel(W::AbstractMatrix) -> RegularSolutionModel

Construct a [`RegularSolutionModel`](@ref) from a symmetric matrix of
interaction parameters in J/mol. Raises if `W` is not square or not symmetric;
the diagonal is ignored. Integer entries are promoted to a floating-point type.
"""
RegularSolutionModel(W::AbstractMatrix) =
    RegularSolutionModel{float(eltype(W))}(W)

"""
    struct RedlichKisterModel{T<:Real} <: AbstractSolidSolutionModel

Binary Redlich-Kister (asymmetric Margules) model for non-ideal solid solutions.
Requires **exactly 2 end-members** per solid solution phase.

Parameters `a0`, `a1`, `a2` are in J/mol and divided by RT inside the activity closure.

Activity coefficients (Guggenheim / ThermoCalc convention):

```
ln γ₁ = (x₂²/RT)[a₀ + a₁(3x₁ − x₂) + a₂(x₁ − x₂)(5x₁ − x₂)]
ln γ₂ = (x₁²/RT)[a₀ − a₁(3x₂ − x₁) + a₂(x₂ − x₁)(5x₂ − x₁)]
```

AD-compatible: all computations propagate `ForwardDiff.Dual` numbers.

# Examples

```jldoctest
julia> m = RedlichKisterModel(a0 = 4000.0, a1 = 500.0)
RedlichKisterModel{Float64}(4000.0, 500.0, 0.0)

julia> m.a0
4000.0
```
"""
struct RedlichKisterModel{T <: Real} <: AbstractSolidSolutionModel
    a0::T   # J/mol — symmetric interaction parameter
    a1::T   # J/mol — asymmetry
    a2::T   # J/mol — higher-order correction
end

"""
    RedlichKisterModel(; a0=0.0, a1=0.0, a2=0.0)

Keyword constructor for [`RedlichKisterModel`](@ref). Parameters are promoted to a
common type.
"""
function RedlichKisterModel(; a0 = 0.0, a1 = 0.0, a2 = 0.0)
    vals = promote(a0, a1, a2)
    return RedlichKisterModel{eltype(vals)}(vals...)
end

# ── Solid solution phase ───────────────────────────────────────────────────────

"""
    abstract type AbstractSolidSolutionPhase end

Base type for solid solution phases. Concrete subtypes group a set of end-member
[`AbstractSpecies`](@ref) and associate them with an [`AbstractSolidSolutionModel`](@ref).
"""
abstract type AbstractSolidSolutionPhase end

"""
    struct SolidSolutionPhase{T<:AbstractSpecies, M<:AbstractSolidSolutionModel}
            <: AbstractSolidSolutionPhase

A solid-solution phase consisting of `end_members` (species with `AS_CRYSTAL` aggregate
state) mixing according to `model`.

End-members are automatically requalified to `SC_SSENDMEMBER` at construction time,
so database species with `SC_COMPONENT` can be passed directly.

# Construction

Use the keyword constructor:
```julia
SolidSolutionPhase(name, end_members; model = IdealSolidSolutionModel())
```

Validation at construction time:
- All end-members must have `aggregate_state == AS_CRYSTAL`.
- [`RedlichKisterModel`](@ref) requires exactly 2 end-members.
- The mixing energy must be convex, unless `instances > 1` or
  `check_convexity = false`.

# A miscibility gap: `instances`

`instances` is how many **coexisting compositions** the declaration may hold. It
is 1 for every phase the shipped data describes, and it is 1 because those models
are convex: a convex mixing energy has one minimum, so one composition describes
the phase.

Inside a spinodal it does not. Where `d²g/dx² < 0` the Gibbs minimum is the
**common-tangent pair** — two compositions of the same substance, coexisting —
and a formulation carrying one amount per species cannot write that down. So
`instances = 2` asks `ChemicalSystem` for a second copy of each end-member, under
a derived symbol (`monosulphate12#2`) sharing the same thermodynamic record, and
the minimization is free to put material in either lobe or in both.

This is how GEM-Selektor represents the same thing: CEMDATA18 ships the AFm and
AFt binaries under two names each, so that the user can declare them twice. The
difference here is only that the duplication is asked for by a keyword rather
than carried in the database.

`instances > 1` is **refused for a convex model**, and that is not a formality:
two instances of a convex phase are degenerate, every split of the amount between
them having the same energy, so the minimum becomes a flat manifold and the
optimizer is asked to choose a point on it for no reason. Inside a spinodal the
common-tangent pair is unique and the degeneracy does not arise.

```julia
# The published AFm sulfate/hydroxide parameters, whose spinodal is
# x in [0.631, 0.914] at 25 C. With one instance this is refused; with two it is
# the case the model was written for.
SolidSolutionPhase("AFm_SO4_OH", [c4ah13, monosulphate];
                   model = RedlichKisterModel(a0 = 20_000.0), instances = 2)
```

# Example

```jldoctest
julia> em1 = Species("Ca2SiO4"; aggregate_state=AS_CRYSTAL, class=SC_COMPONENT);

julia> em2 = Species("Ca3Si2O7"; aggregate_state=AS_CRYSTAL, class=SC_COMPONENT);

julia> ss = SolidSolutionPhase("CSH", [em1, em2])
SolidSolutionPhase{Species{Int64}, IdealSolidSolutionModel}
  name: CSH
  end-members (2): Ca2SiO4, Ca3Si2O7
  model: IdealSolidSolutionModel

julia> class(end_members(ss)[1])
SC_SSENDMEMBER::Class = 5
```
"""
struct SolidSolutionPhase{T <: AbstractSpecies, M <: AbstractSolidSolutionModel} <:
    AbstractSolidSolutionPhase
    name::String
    end_members::Vector{T}
    model::M
    # How many coexisting compositions this declaration is allowed to hold. One
    # for every convex phase, which is every phase the shipped data describes.
    # Greater than one only inside a miscibility gap; see the keyword
    # constructor, which refuses it otherwise.
    instances::Int
    # The name of the declaration this phase is an instance of. Equal to `name`
    # for an ordinary phase, and for the first instance of a multi-instance one;
    # the later instances are named `"$declared#k"`. `ChemicalSystem` uses it to
    # tell "the same substance declared twice on purpose" from "the same
    # substance declared twice by mistake", which it refuses.
    declared::String
end

"""
    SolidSolutionPhase(name, end_members; model=IdealSolidSolutionModel())

Construct and validate a [`SolidSolutionPhase`](@ref).

End-members whose `class` is not already `SC_SSENDMEMBER` are automatically
requalified via [`with_class`](@ref). Passing database species with
`SC_COMPONENT` therefore works directly, without a prior call to `with_class`.
"""
function SolidSolutionPhase(
        name::AbstractString,
        end_members::AbstractVector{<:AbstractSpecies};
        model::AbstractSolidSolutionModel = IdealSolidSolutionModel(),
        check_convexity::Bool = true, T::Real = 298.15,
        instances::Integer = 1, declared::AbstractString = name,
    )
    instances >= 1 || error(
        "SolidSolutionPhase \"$name\": `instances` is how many coexisting " *
            "compositions the phase may take, so it is at least 1; got $instances."
    )
    for sp in end_members
        aggregate_state(sp) == AS_CRYSTAL ||
            error(
            "SolidSolutionPhase: end-member \"$(symbol(sp))\" must have " *
                "aggregate_state = AS_CRYSTAL (got $(aggregate_state(sp)))",
        )
    end
    if model isa RedlichKisterModel
        length(end_members) == 2 ||
            error(
            "RedlichKisterModel requires exactly 2 end-members, " *
                "got $(length(end_members))",
        )
    end
    # A mixing energy that is concave somewhere does not describe one phase there.
    #
    # Refused rather than warned, and at construction rather than at the solve,
    # because what fails downstream fails obscurely: the minimization returns a
    # composition it cannot certify, with a small element-balance residual and no
    # phase reported missing, and nothing in that output names the cause. Measured
    # on the CEM II of `these_abcael`, whose two AFm/AFt Redlich-Kister sets are
    # concave over `x in [0.63, 0.91]` and `[0.56, 0.83]`: the solve stopped at
    # stationarity 9.2e-5 and element balance 2.6e-3, and the same eight phases
    # with ideal mixing certified to 1.1e-14.
    #
    # The physics is that the Gibbs minimum inside a spinodal is two coexisting
    # compositions, not one. Representing that needs the binary declared twice,
    # which a formulation with one entry per species cannot do — and neither
    # GEM-Selektor nor Reaktoro detects the condition either: both assume
    # convexity and leave the duplication to the user.
    #
    # `check_convexity = false` proceeds anyway, for a caller who knows the answer
    # stays outside the gap. The optimality certificate then loses its ground,
    # since its sufficiency rests on the problem being convex.
    #
    # `instances > 1` inverts the test rather than skipping it. Two instances of
    # a CONVEX phase are degenerate -- every way of splitting the amount between
    # them has the same energy, so the minimum is a flat manifold and the
    # optimizer is asked to pick a point on it for no reason. Inside a spinodal
    # they are not degenerate at all: the minimum is the common-tangent pair, and
    # it is unique. So a second instance is admitted exactly where it is needed
    # and refused where it would only add a null direction.
    if instances > 1
        gap = spinodal_interval(model, length(end_members); T = T)
        gap === nothing && error(
            "SolidSolutionPhase \"$name\": `instances = $instances` asks for " *
                "$instances coexisting compositions of this phase, but its mixing " *
                "energy is CONVEX at T = $(T) K, so it has one. The instances would " *
                "be degenerate -- every split of the amount between them has the " *
                "same energy -- and the minimization would be asked to choose a " *
                "point on a flat manifold. Use `instances = 1`, or a model whose " *
                "energy has a spinodal (`spinodal_interval` reports it)."
        )
    elseif check_convexity
        gap = spinodal_interval(model, length(end_members); T = T)
        gap === nothing || error(
            "SolidSolutionPhase \"$name\": the mixing energy of this model is " *
                "CONCAVE for x in [$(round(gap[1]; digits = 3)), " *
                "$(round(gap[2]; digits = 3))] at T = $(T) K, so the phase " *
                "unmixes there: the Gibbs minimum inside that interval is two " *
                "coexisting compositions, not one, and this formulation has a " *
                "single amount per species to describe it with. Declare it with " *
                "`instances = 2` to give it two, which is what a miscibility gap " *
                "needs; or use ideal mixing, or parameters that keep the energy " *
                "convex; or pass `check_convexity = false` to proceed with one " *
                "composition anyway — in which case the optimality certificate no " *
                "longer proves anything, its sufficiency resting on convexity."
        )
    end

    qualified = [
        class(sp) == SC_SSENDMEMBER ? sp : with_class(sp, SC_SSENDMEMBER)
            for sp in end_members
    ]
    T = eltype(qualified)
    return SolidSolutionPhase{T, typeof(model)}(
        String(name), collect(T, qualified), model, Int(instances), String(declared)
    )
end

# ── Convexity of the mixing energy ────────────────────────────────────────────

"""
    spinodal_interval(model, n_members; T = 298.15) -> Union{Nothing, Tuple{Float64, Float64}}

The interval of composition over which a binary mixing energy is **concave**, or
`nothing` when it is convex throughout.

A solid solution exists as one homogeneous phase only where its molar Gibbs
energy of mixing is convex. Where `d²g/dx² < 0` — inside the spinodal — the
minimum of `G` is not a single composition but **two coexisting ones**: the phase
unmixes, and the equilibrium is a miscibility gap.

For a binary with the package's Redlich-Kister convention,

```math
\\frac{g}{RT} = x\\ln x + (1-x)\\ln(1-x)
             + x(1-x)\\left[A_0 + A_1(2x-1) + A_2(2x-1)^2\\right] ,
\\qquad A_k = a_k/RT ,
```

and the second derivative is evaluated on a grid rather than in closed form, so
that the same routine covers `a₂` and the symmetric
[`RegularSolutionModel`](@ref) without a separate derivation. The classical
symmetric result is recovered as a check: `d²g/dx²` at `x = 1/2` is `4 − 2A₀`, so
a regular solution unmixes above `W = 2RT`.

`T` is the temperature the parameters are read at; they are stored in J/mol and
the criterion is `a/RT`, so a model that is convex at 25 °C may not be at 5 °C.

Returns `nothing` for an ideal model, whose second derivative is `1/x + 1/(1-x)`
and therefore positive everywhere, and for a phase with more than two
end-members, where a one-dimensional scan is not the right test — see the note in
[`SolidSolutionPhase`](@ref).

See also: [`RedlichKisterModel`](@ref), [`RegularSolutionModel`](@ref).
"""
function spinodal_interval(
        model::AbstractSolidSolutionModel, n_members::Int; T::Real = 298.15
    )
    n_members == 2 || return nothing
    RT = 8.31446261815324 * T          # J/(mol K), CODATA
    A0, A1, A2 = if model isa RedlichKisterModel
        (model.a0 / RT, model.a1 / RT, model.a2 / RT)
    elseif model isa RegularSolutionModel
        (model.W[1, 2] / RT, 0.0, 0.0)
    else
        return nothing                      # ideal: convex everywhere
    end
    (A0 == 0 && A1 == 0 && A2 == 0) && return nothing

    gx(x) = x * log(x) + (1 - x) * log(1 - x) +
        x * (1 - x) * (A0 + A1 * (2x - 1) + A2 * (2x - 1)^2)

    xs = range(1.0e-3, 1 - 1.0e-3; length = 2001)
    h = step(xs)
    lo, hi = Inf, -Inf
    for i in 2:(length(xs) - 1)
        d2 = (gx(xs[i + 1]) - 2gx(xs[i]) + gx(xs[i - 1])) / h^2
        if d2 < 0
            lo = min(lo, xs[i])
            hi = max(hi, xs[i])
        end
    end
    return isfinite(lo) ? (lo, hi) : nothing
end

# ── Accessors ─────────────────────────────────────────────────────────────────

"""
    name(ss::SolidSolutionPhase) -> String

Return the name of the solid solution phase.
"""
name(ss::SolidSolutionPhase) = ss.name

"""
    end_members(ss::SolidSolutionPhase) -> Vector{<:AbstractSpecies}

Return the end-member species of the solid solution phase.
"""
end_members(ss::SolidSolutionPhase) = ss.end_members

"""
    model(ss::SolidSolutionPhase) -> AbstractSolidSolutionModel

Return the activity model of the solid solution phase.
"""
model(ss::SolidSolutionPhase) = ss.model

# ── Display ───────────────────────────────────────────────────────────────────

function Base.show(io::IO, ss::SolidSolutionPhase{T, M}) where {T, M}
    em_names = join(symbol.(ss.end_members), ", ")
    println(io, "SolidSolutionPhase{$T, $M}")
    println(io, "  name: $(ss.name)")
    println(io, "  end-members ($(length(ss.end_members))): $em_names")
    print(io, "  model: $M")
    ss.instances > 1 && print(io, "\n  instances: $(ss.instances) (miscibility gap)")
    return nothing
end
