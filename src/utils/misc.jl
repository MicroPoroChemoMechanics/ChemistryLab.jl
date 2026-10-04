# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

using Crayons
using DynamicQuantities
using ForwardDiff

"""
    _adim(x) -> Real

`x` as a bare number, **checking** that it carries no dimension.

# What it replaces, and why that had to go

`thermo_factories.jl` used to extend some thirty `Base` math functions to
`DynamicQuantities.Quantity` by stripping the unit first:

```julia
Base.log(x::Quantity) = log(ustrip(x))     # removed
```

That is type piracy — the function and the type both belong elsewhere — so it
applied to **every** package loaded beside this one, whether or not it asked. And
it did not merely remove a dimension check; it made the result depend on an
invisible normalization, because `ustrip` returns the value in SI **base** units:

```
log(2u"m")       = 0.693      instead of raising
log(1u"mol/L")   = 6.908      the same quantity written two ways,
log(1u"mmol/L")  = 0.0        two different answers
```

# What this does instead

Dispatch rather than a branch, so `Float64`, `ForwardDiff.Dual` and
`Symbolics.Num` all take the first method untouched — none of them carries a
dimension — and only a genuine `Quantity` pays for the test.

A dimensionless `Quantity` passes, which is the case the extensions existed for:
a ratio whose units cancel is a number, and saying so explicitly is the point.
Anything else raises, which is the check being restored rather than a new
restriction.
"""
@inline _adim(x::Real) = x
@inline function _adim(x::DynamicQuantities.AbstractQuantity)
    iszero(dimension(x)) || throw(
        DimensionError(
            x, one(x)
        )
    )
    return ustrip(x)
end

"""Terminal color for charge display (cyan bold)."""
const COL_CHARGE = crayon"cyan bold"
"""Terminal color for parentheses (magenta bold)."""
const COL_PAR = crayon"magenta bold"
"""Terminal color for internal stoichiometric coefficients (red bold)."""
const COL_STOICH_INT = crayon"red bold"
"""Terminal color for external stoichiometric coefficients (yellow bold)."""
const COL_STOICH_EXT = crayon"yellow bold"

"""
    _primal(x::Real) -> Real

Extract the primal (non-derivative) value of `x` for use in control-flow decisions.

For plain `Real` numbers this is a no-op.  For forward-mode AD types such as
`ForwardDiff.Dual` (which are subtypes of `Real`) this returns `x` itself, but
since those types overload comparison operators to compare only the primal value,
writing `_primal(x) > threshold` makes the intent explicit and ensures that
branching decisions are always taken on physical values rather than derivatives.

Adding a method `_primal(x::YourDualType) = x.value` enables support for
non-standard AD frameworks that do not overload comparison operators.

# Examples
```julia
julia> _primal(3.14)
3.14
```
"""
@inline _primal(x::Real) = x
@inline _primal(x::ForwardDiff.Dual) = ForwardDiff.value(x)

"""
    print_title(title; crayon=:none, indent="", style=:none)

Print a formatted title with optional styling and indentation.

# Arguments

  - `title`: title string to display.
  - `crayon`: Crayon for colored output (default `:none`).
  - `indent`: left indentation string (default `""`).
  - `style`: formatting style - `:none`, `:underline`, or `:box` (default `:none`).

# Examples

```julia
julia> print_title("Test Title"; style=:none)
Test Title

julia> print_title("Section"; style=:underline, indent="  ")
  Section
  ───────
```
"""
function print_title(title; crayon = :none, indent = "", style = :none)
    draw = crayon != :none ? x -> println(crayon(x)) : println
    if style == :underline
        width = length(title)
        draw(indent * "$title")
        draw(indent * "─"^width)
    elseif style == :box
        width = length(title) + 4
        draw(indent * "┌" * "─"^width * "┐")
        draw(indent * "│  $title  │")
        draw(indent * "└" * "─"^width * "┘")
    else
        draw(indent * "$title")
    end
    return print(crayon"reset")
end

"""
    root_type(T) -> Type

Get the root type wrapper from a potentially parameterized type.
"""
root_type(T) = T isa UnionAll ? root_type(T.body) : T.name.wrapper

"""
    safe_ustrip(unit::UnionAbstractQuantity, q::UnionAbstractQuantity) -> Number
    safe_ustrip(unit::UnionAbstractQuantity, q) -> Number

Safely remove units from `q`, converting to `unit` when their dimensions match.
If `q` is a plain numeric value (no units) it is returned unchanged.

When the dimensions of `unit` and `q` differ (e.g. `unit` is dimensionless but
`q` carries a physical dimension) the conversion is skipped and `q` is stripped
in its native unit.  This prevents a `DimensionError` in contexts such as
`ThermoFactory` where a dimensionless placeholder unit `u"1"` may be passed
alongside a quantity with physical dimensions.

# Examples

```julia
julia> safe_ustrip(u"m", 5u"m")
5.0

julia> safe_ustrip(u"cm", 1u"m")
100.0

julia> safe_ustrip(u"m", 3.0)
3.0
```
"""
function safe_ustrip(unit::UnionAbstractQuantity, q::UnionAbstractQuantity)
    # Avoid comparing dimension types directly (SymbolicDimensions vs Dimensions
    # comparison raises an error in DynamicQuantities ≤ 1.11).  Instead, attempt
    # the direct conversion, then try uconvert + ustrip as fallback.
    try
        return ustrip(unit, q)
    catch
        try
            return ustrip(uconvert(unit, q))
        catch
            return ustrip(q)
        end
    end
end

safe_ustrip(::UnionAbstractQuantity, q) = ustrip(q)

"""
    safe_uconvert(qout::UnionAbstractQuantity{<:Any, <:AbstractSymbolicDimensions}, q::UnionAbstractQuantity{<:Any, <:Dimensions}) -> UnionAbstractQuantity
    safe_uconvert(qout::UnionAbstractQuantity{<:Any, <:AbstractSymbolicDimensions}, q) -> Any

Convert `q` to the unit type represented by `qout` using `uconvert` when `q`
is a quantity with compatible dimensions. If `q` is a plain numeric value,
return `q` unchanged.

# Arguments

  - `qout` : unitful quantity target (e.g., u"m", u"kg").
  - `q` : a quantity or a plain numeric value.

# Returns

  - Converted quantity (of type `qout`) or the original `q` when `q` has no units.

# Examples

```julia
julia> safe_uconvert(us"m", 5u"cm")
0.05 m

julia> safe_uconvert(us"m", 3.0)
3.0
```
"""
safe_uconvert(qout::UnionAbstractQuantity{<:Any, <:AbstractSymbolicDimensions}, q::UnionAbstractQuantity{<:Any, <:Dimensions}) = uconvert(qout, q)

safe_uconvert(::UnionAbstractQuantity{<:Any, <:AbstractSymbolicDimensions}, q) = q

"""
    _ensure_unit(unit, x::Real) -> Quantity
    _ensure_unit(unit, x) -> Quantity

Ensure a value is a `Quantity` with the given `unit`.

  - Plain `Real` → interpreted as SI, wrapped: `x * unit`.
  - `Quantity` → converted to `unit` via [`safe_uconvert`](@ref).
"""
_ensure_unit(unit, x::Real) = x * unit
_ensure_unit(unit, x) = safe_uconvert(unit, x)

"""
    force_uconvert(qout::UnionAbstractQuantity, q) -> AbstractQuantity

Strip units from `q` (converting to `qout` dimensions when possible) and
re-attach the unit of `qout`. Unlike `safe_uconvert`, this always returns
a quantity with the units of `qout`, even when the input is dimensionless.
"""
force_uconvert(qout::UnionAbstractQuantity, q) = safe_ustrip(qout, q) * qout

# A vector in the common number type of its entries: a function returning a
# constant beside one returning a dual number leaves neither a `Vector{Real}`
# nor a failed conversion to `Float64`.
_promoted(v) = convert(Vector{mapreduce(typeof, promote_type, v; init = Float64)}, v)

# ── the number type a value carries ──────────────────────────────────────────
#
# What an output must be typed by to keep a derivative: every input it is
# computed from, the parameters a closure captures included. Anything narrower
# either raises (a dual stored into a `Float64` container) or drops the
# derivative.

const _CAPTURE_DEPTH = 8

# The dual number type a value carries anywhere inside it (`Float64` when none):
# its fields, the variables a closure captures, the entries of a container.
_captured_number_type(x) = _captured_number_type(x, 0)
function _captured_number_type(x, depth::Int)
    depth > _CAPTURE_DEPTH && return Float64
    x isa ForwardDiff.Dual && return typeof(x)
    x isa Union{Number, AbstractString, Symbol, Module, Type, Nothing, Missing, Char, Function} &&
        !(x isa Function && nfields(x) > 0) && return Float64
    x isa DynamicQuantities.AbstractQuantity && return _captured_number_type(ustrip(x), depth + 1)
    if x isa AbstractArray
        E = eltype(x)
        E <: ForwardDiff.Dual && return E
        (E <: Number && isconcretetype(E)) && return Float64
        return mapreduce(e -> _captured_number_type(e, depth + 1), promote_type, x; init = Float64)
    end
    x isa AbstractDict && return mapreduce(e -> _captured_number_type(e, depth + 1), promote_type, values(x); init = Float64)
    x isa Union{Tuple, NamedTuple} && return mapreduce(e -> _captured_number_type(e, depth + 1), promote_type, values(x); init = Float64)
    T = Float64
    for k in 1:nfields(x)
        isdefined(x, k) || continue
        T = promote_type(T, _captured_number_type(getfield(x, k), depth + 1))
    end
    return T
end

# The same question asked of a parameter tuple at the level of TYPES, so that it
# costs nothing where it is asked at every evaluation of an activity model.
#
# It used to map over the VALUES of the tuple. An equilibrium hands a model a
# handful of parameters, but a kinetic run hands it the run's own parameters,
# some fifty fields, and a `map` over a tuple that long is not unrolled: each
# field was dispatched at run time, at every evaluation. Measured on the
# pore-humidity test, 44 % of the integration went there. The answer depends on
# the types alone, so it is computed from them, and `:foldable` lets the
# compiler fold it into a constant.
_number_type_of(x) = _number_type_t(typeof(x))

_number_type_t(::Type{T}) where {T <: ForwardDiff.Dual} = T
_number_type_t(::Type{<:DynamicQuantities.AbstractQuantity{T}}) where {T} = _number_type_t(T)
_number_type_t(::Type{<:AbstractArray{T}}) where {T <: ForwardDiff.Dual} = T
_number_type_t(::Type) = Float64
Base.@assume_effects :foldable function _number_type_t(::Type{T}) where {T <: Union{Tuple, NamedTuple}}
    R = Float64
    for F in fieldtypes(T)
        R = promote_type(R, _number_type_t(F))
    end
    return R
end
