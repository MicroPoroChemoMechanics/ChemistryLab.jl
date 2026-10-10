# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

using DynamicQuantities
using ForwardDiff
using OrderedCollections
using Symbolics

"""
    ADIM_MATH_FUNCTIONS

The mathematical functions that take a **dimensionless** argument.

Used to extend [`SymbolicFunc`](@ref) and [`NumericFunc`](@ref), which are this
package's own types, so that a compiled or symbolic expression can be written
with them. It is **not** used to extend them to `DynamicQuantities.Quantity`:
that was type piracy, it applied to every package loaded beside this one, and it
replaced a dimension error by an answer that depended on SI normalization. See
[`_adim`](@ref), which restores the check at the two places that needed one.
"""
const ADIM_MATH_FUNCTIONS = [
    :log,
    :log10,
    :log2,
    :log1p,
    :exp,
    :expm1,
    :exp2,
    :exp10,
    :sin,
    :cos,
    :tan,
    :csc,
    :sec,
    :cot,
    :asin,
    :acos,
    :atan,
    :acsc,
    :asec,
    :acot,
    :sinh,
    :cosh,
    :tanh,
    :csch,
    :sech,
    :coth,
    :asinh,
    :acosh,
    :atanh,
    :acsch,
    :asech,
    :acoth,
    :erf,
    :erfc,
    :gamma,
    :lgamma,
]

# ── Abstract type ──────────────────────────────────────────────────────────────

"""
    AbstractFunc

Abstract type for objects that can be called like functions.
"""
abstract type AbstractFunc end

# ── NumericFunc ────────────────────────────────────────────────────────────────

"""
    NumericFunc{N, F, R, Q} <: AbstractFunc

Closure-backed thermodynamic function for models that cannot be represented as symbolic
expressions (e.g. HKF [Helgeson1981](@cite), or any other numeric model).
Calling convention is identical to `SymbolicFunc`: `f(; T=..., P=..., unit=false)`.

Variable values are resolved in order: `kwarg > refs > _NF_DEFAULT_REFS`.
`refs` stores `Quantity` values so that unit information is preserved.

# Fields

  - `compiled`: closure `(vars...) → value` in SI units.
  - `vars`: names of the positional arguments (e.g. `(:T, :P)`).
  - `refs`: `NamedTuple` of default variable values as `Quantity` (e.g. `(T=298.15u"K", P=1e5u"Pa")`).
  - `unit`: output unit (DynamicQuantities `Quantity`).
"""
struct NumericFunc{N, F, R <: NamedTuple, Q} <: AbstractFunc
    compiled::F
    vars::NTuple{N, Symbol}
    refs::R
    unit::Q
end

# Global fallback for vars absent from refs
const _NF_DEFAULT_REFS = (T = T_STANDARD_Q, P = P_STANDARD_Q)

"""
    NumericFunc(f, vars, unit)

Construct a `NumericFunc` with no refs (fallback to `_NF_DEFAULT_REFS`).
"""
NumericFunc(f, vars::NTuple{N, Symbol}, unit) where {N} =
    NumericFunc(f, vars, NamedTuple(), unit)

"""
    NumericFunc(f, unit)

Backward-compatible constructor: assumes `vars = (:T, :P)` and empty refs.
"""
NumericFunc(f, unit) = NumericFunc(f, (:T, :P), NamedTuple(), unit)

@inline function (nf::NumericFunc{N})(; kwargs...) where {N}
    vals = ntuple(N) do i
        v = nf.vars[i]
        raw = haskey(kwargs, v) ? kwargs[v] :
            get(nf.refs, v, get(_NF_DEFAULT_REFS, v, nothing))
        ustrip(raw)   # Quantities must already be in SI; plain numbers are passed through
    end
    val = nf.compiled(vals...)
    return get(kwargs, :unit, false) ? val * nf.unit : val
end

function Base.show(io::IO, nf::NumericFunc)
    print(io, "NumericFunc [", dimension(nf.unit), "]")
    if !isempty(nf.vars)
        print(io, " ◆ vars=(", join(nf.vars, ", "), ")")
    end
    return if !isempty(nf.refs)
        print(io, " ◆ ", join(["$k=$v" for (k, v) in pairs(nf.refs)], ", "))
    end
end

function Base.show(io::IO, ::MIME"text/plain", nf::NumericFunc)
    println(io, "NumericFunc:")
    print(io, "  Unit: [", dimension(nf.unit), "]")
    print(io, "\n  Variables: ", join(nf.vars, ", "))
    return if !isempty(nf.refs)
        print(io, "\n  References: ", join(["$k=$v" for (k, v) in pairs(nf.refs)], ", "))
    end
end

# ---------------------------------------------------------------------------
# Arithmetic for NumericFunc
# ---------------------------------------------------------------------------

# Index-mapping helpers (used by binary NF-NF and cross-type operations)
@inline function _nf_call_from_combined(nf::NumericFunc, idx, args)
    vals = ntuple(i -> args[idx[i]], length(nf.vars))
    return nf.compiled(vals...)
end

_combined_nf_vars(nf1::NumericFunc, nf2::NumericFunc) =
    Tuple(union(nf1.vars, nf2.vars))

_var_indices(src_vars, combined_vars) =
    ntuple(i -> findfirst(==(src_vars[i]), combined_vars)::Int, length(src_vars))

# Unary
Base.:-(nf::NumericFunc) =
    NumericFunc((args...) -> -nf.compiled(args...), nf.vars, nf.refs, nf.unit)

# NF op NF — vars are united, refs are merged (second operand takes precedence)
function _combine_nf(op, nf1::NumericFunc, nf2::NumericFunc)
    new_vars = _combined_nf_vars(nf1, nf2)
    idx1 = _var_indices(nf1.vars, new_vars)
    idx2 = _var_indices(nf2.vars, new_vars)
    combined = (args...) -> op(
        _nf_call_from_combined(nf1, idx1, args),
        _nf_call_from_combined(nf2, idx2, args),
    )
    return NumericFunc(
        combined, new_vars, merge(nf1.refs, nf2.refs), oneunit(op(nf1.unit, nf2.unit))
    )
end

for op in (:+, :-, :*, :/)
    @eval Base.$op(nf1::NumericFunc, nf2::NumericFunc) =
        _combine_nf($op, nf1, nf2)
end

# NF op Number (and vice versa)
# For +/-: result unit = nf.unit (scalar assumed dimensionally compatible)
# For *//: result unit from unit algebra
function _combine_nf_scalar(op, nf::NumericFunc, x::Number, result_unit)
    return NumericFunc(
        (args...) -> op(nf.compiled(args...), ustrip(x)), nf.vars, nf.refs, result_unit
    )
end

function _combine_scalar_nf(op, x::Number, nf::NumericFunc, result_unit)
    return NumericFunc(
        (args...) -> op(ustrip(x), nf.compiled(args...)), nf.vars, nf.refs, result_unit
    )
end

for op in (:+, :-)
    @eval begin
        Base.$op(nf::NumericFunc, x::Number) =
            _combine_nf_scalar($op, nf, x, nf.unit)
        Base.$op(x::Number, nf::NumericFunc) =
            _combine_scalar_nf($op, x, nf, nf.unit)
    end
end

for op in (:*, :/)
    @eval begin
        Base.$op(nf::NumericFunc, x::Number) =
            _combine_nf_scalar($op, nf, x, oneunit(Base.$op(nf.unit, x)))
        Base.$op(x::Number, nf::NumericFunc) =
            _combine_scalar_nf($op, x, nf, oneunit(Base.$op(x, nf.unit)))
    end
end

# ── SymbolicFunc ───────────────────────────────────────────────────────────────

"""
    SymbolicFunc

Thermodynamic function with symbolic expression and compiled evaluation.
"""
struct SymbolicFunc{N, R <: NamedTuple, T, D} <: AbstractFunc
    symbolic::Num
    vars::NTuple{N, Symbol}
    refs::R
    compiled::RuntimeGeneratedFunction
    unit::Quantity{T, D}
end

# Call methods — specializations for N = 1, 2, 3 avoid runtime tuple dispatch overhead

@inline function (sf::SymbolicFunc{1})(; kwargs...)
    v = sf.vars[1]
    val = haskey(kwargs, v) ? kwargs[v] : sf.refs[v]
    return if get(kwargs, :unit, false)
        sf.compiled(ustrip(val)) * sf.unit
    else
        sf.compiled(ustrip(val))
    end
end

@inline function (sf::SymbolicFunc{2})(; kwargs...)
    v1, v2 = sf.vars
    val1 = haskey(kwargs, v1) ? kwargs[v1] : get(sf.refs, v1, nothing)
    val2 = haskey(kwargs, v2) ? kwargs[v2] : get(sf.refs, v2, nothing)
    return if get(kwargs, :unit, false)
        sf.compiled(ustrip(val1), ustrip(val2)) * sf.unit
    else
        sf.compiled(ustrip(val1), ustrip(val2))
    end
end

@inline function (sf::SymbolicFunc{3})(; kwargs...)
    v1, v2, v3 = sf.vars
    val1 = haskey(kwargs, v1) ? kwargs[v1] : get(sf.refs, v1, nothing)
    val2 = haskey(kwargs, v2) ? kwargs[v2] : get(sf.refs, v2, nothing)
    val3 = haskey(kwargs, v3) ? kwargs[v3] : get(sf.refs, v3, nothing)
    return if get(kwargs, :unit, false)
        sf.compiled(ustrip(val1), ustrip(val2), ustrip(val3)) * sf.unit
    else
        sf.compiled(ustrip(val1), ustrip(val2), ustrip(val3))
    end
end

@inline function (sf::SymbolicFunc{N})(; kwargs...) where {N}
    if isempty(kwargs)
        var_values = ntuple(i -> ustrip(sf.refs[sf.vars[i]]), N)
        return sf.compiled(var_values...)
    else
        merged = merge(sf.refs, kwargs)
        var_values = ntuple(i -> ustrip(merged[sf.vars[i]]), N)
        return if get(kwargs, :unit, false)
            sf.compiled(var_values...) * sf.unit
        else
            sf.compiled(var_values...)
        end
    end
end

"""
    Base.show(io::IO, sf::SymbolicFunc)

Compact string representation of a SymbolicFunc.
"""
function Base.show(io::IO, sf::SymbolicFunc)
    print(io, sf.symbolic, " [", dimension(sf.unit), "]")
    if !isempty(sf.vars)
        print(io, " ◆ vars=(", join(sf.vars, ", "), ")")
    end
    return if !isempty(sf.refs)
        print(io, " ◆ ", join(["$k=$v" for (k, v) in pairs(sf.refs)], ", "))
    end
end

"""
    Base.show(io::IO, ::MIME"text/plain", sf::SymbolicFunc)

Detailed string representation of a SymbolicFunc.
"""
function Base.show(io::IO, ::MIME"text/plain", sf::SymbolicFunc)
    println(io, "SymbolicFunc:")
    print(io, "  Expression: ")
    print(io, sf.symbolic, " [", dimension(sf.unit), "]")
    print(io, "\n  Variables: ", join(sf.vars, ", "))
    return if !isempty(sf.refs)
        print(io, "\n  References: ", join(["$k=$v" for (k, v) in pairs(sf.refs)], ", "))
    end
end

# ── ThermoFactory ──────────────────────────────────────────────────────────────

"""
    extract_vars_params(expr, vars) -> (Vector{Symbol}, Vector{Symbol})

Identify variables and parameters in an expression.

# Arguments

  - `expr`: symbolic expression to analyze.
  - `vars`: list of symbols considered as variables (others are parameters).

# Returns

  - Tuple of (variables, parameters) found in the expression.
"""
function extract_vars_params(expr, vars)
    params = Symbol[]
    vars_set = Set(vars)
    newvars = Symbol[]

    _scan_symbols!(newvars, params, vars_set, expr)
    unique!(newvars)
    unique!(params)
    return newvars, params
end

# The symbols of `ex`, to `newvars` when in `vars_set` and to `params` otherwise,
# leaving out the first argument of each expression: the function of a call.
function _scan_symbols!(newvars, params, vars_set, ex)
    if ex isa Symbol
        ex ∈ vars_set ? push!(newvars, ex) : push!(params, ex)
    elseif ex isa Expr
        for arg in @view ex.args[2:end]
            _scan_symbols!(newvars, params, vars_set, arg)
        end
    end
    return nothing
end

"""
    compile_symbolic(symbolic_expr, var_symbols)

Compile a symbolic expression to RuntimeGeneratedFunction.
"""
function compile_symbolic(symbolic_expr, var_symbols)
    return if length(var_symbols) == 1
        Symbolics.build_function(symbolic_expr, var_symbols[1]; expression = Val(false))
    else
        Symbolics.build_function(symbolic_expr, var_symbols...; expression = Val(false))
    end
end

# The units given to a factory, by symbol, and each read as a unit: a string is
# parsed, a quantity gives its unit, anything else is dimensionless.
_units_by_symbol(nt::NamedTuple) = Dict(pairs(nt))
_units_by_symbol(v::AbstractVector{<:Pair}) = Dict(v)
_as_unit(s::String) = uparse(s)
_as_unit(q::AbstractQuantity) = oneunit(q)
_as_unit(::Any) = u"1"

"""
    ThermoFactory{Q}

Factory for creating `SymbolicFunc` instances from expressions.
Units for each variable/parameter and the output unit are stored explicitly, so
no symbolic unit propagation is needed.
"""
struct ThermoFactory{Q}
    symbolic::Num
    vars::OrderedDict{Symbol, Num}
    params::OrderedDict{Symbol, Num}
    units::Dict{Symbol, Q}                                      # unit per var/param
    output_unit::Q                                              # explicit output unit
    cache::Dict{UInt64, Tuple{Num, RuntimeGeneratedFunction}}

    function ThermoFactory(
            symbolic::Num,
            vars::AbstractDict{Symbol, Num},
            params::AbstractDict{Symbol, Num},
            units::Dict{Symbol, Q},
            output_unit::Q,
        ) where {Q}
        return new{Q}(
            symbolic, vars, params, units, output_unit,
            Dict{UInt64, Tuple{Num, RuntimeGeneratedFunction}}(),
        )
    end
end

"""
    ThermoFactory(expr, vars=[:T, :P, :t, :x, :y, :z]; units=nothing, output_unit=nothing) -> ThermoFactory

Create a `ThermoFactory` from a symbolic expression.

# Arguments

  - `expr`: symbolic expression (Expr or Symbol).
  - `vars`: list of variable symbols (default: T, P, t, x, y, z).
  - `units`: dictionary mapping symbols to their units; a symbol it does not
    list is dimensionless.
  - `output_unit`: the unit of the expression's value (dimensionless when
    `nothing`).
"""
function ThermoFactory(
        expr,
        vars = [:T, :P, :t, :x, :y, :z];
        units = nothing,
        output_unit = nothing,
    )
    found_vars, params = extract_vars_params(expr, vars)
    var_sym_dict = OrderedDict{Symbol, Num}(v => Symbolics.variable(v) for v in found_vars)
    param_sym_dict = OrderedDict{Symbol, Num}(p => Symbolics.variable(p) for p in params)

    _fallback = u"1"
    if !isnothing(units)
        dict_units = _units_by_symbol(units)
        unit_dict = Dict{Symbol, typeof(_fallback)}(
            sym => (haskey(dict_units, sym) ? _as_unit(dict_units[sym]) : _fallback)
                for sym in Iterators.flatten((keys(var_sym_dict), keys(param_sym_dict)))
        )
    else
        unit_dict = Dict{Symbol, typeof(_fallback)}(
            sym => _fallback
                for sym in Iterators.flatten((keys(var_sym_dict), keys(param_sym_dict)))
        )
    end

    out_unit = isnothing(output_unit) ? _fallback : _as_unit(output_unit)

    all_symbols = merge(var_sym_dict, param_sym_dict)
    symbolic = Symbolics.wrap(Symbolics.parse_expr_to_symbolic(expr, all_symbols))
    return ThermoFactory(symbolic, var_sym_dict, param_sym_dict, unit_dict, out_unit)
end

"""
    get_unit(factory::ThermoFactory) -> Quantity

Get the output unit of the expression managed by the factory.
"""
function get_unit(factory::ThermoFactory)
    return factory.output_unit
end

"""
    get_unit(factory::ThermoFactory, sym::Symbol) -> Quantity

Get the unit of a specific variable or parameter in the factory.
"""
function get_unit(factory::ThermoFactory, sym::Symbol)
    return get(factory.units, sym, u"1")
end

# The compiled functions of every factory are memoized in its `cache`, and the
# factories are global (`THERMO_FACTORIES`): independent calculations running on
# threads, each building its species, insert into the same `Dict` at once.
# Concurrent insertions corrupt a `Dict`, and the documentation build, computing
# six coupled trajectories on four threads, crashed on it with a segmentation
# fault. The lock covers the lookup as well, since a read that races a resize is
# no safer; the compilation inside it runs once per key, so after the first call
# the lock guards a lookup only.
const _THERMO_FACTORY_LOCK = ReentrantLock()

# The cache entry of the expression compiled with its parameters as arguments.
const _GENERIC_CACHE_KEY = hash(:compiled_with_parameters_as_arguments)

_has_dual_value(x) = (v = x isa AbstractQuantity ? ustrip(x) : x; v isa ForwardDiff.Dual)

"""
    (factory::ThermoFactory)(; kwargs...)

Create a `SymbolicFunc` with caching for optimal performance. Safe to call from
several threads at once.
"""
function (factory::ThermoFactory)(; kwargs...)
    param_vals = Dict{Symbol, Any}(p => get(kwargs, p, 0.0) for p in keys(factory.params))
    refs = NamedTuple(
        [
            v => force_uconvert(get_unit(factory, v), kwargs[v]) for
                v in keys(factory.vars) if haskey(kwargs, v)
        ]
    )
    unit = get_unit(factory)

    # Parameters carrying dual numbers (thermodynamic data being differentiated)
    # are not substituted: the compiled function would be keyed on their
    # values, a dual hashing as its value, and the one compiled for the plain
    # data returned with every derivative zero. Nor would a new compilation per
    # value do, the cache growing by one function at every iteration of a fit.
    # The expression is compiled once with the parameters as arguments, and
    # evaluated at them.
    if any(_has_dual_value, values(param_vals))
        generic = lock(_THERMO_FACTORY_LOCK) do
            last(
                get!(factory.cache, _GENERIC_CACHE_KEY) do
                    args = vcat(collect(values(factory.vars)), collect(values(factory.params)))
                    (factory.symbolic, compile_symbolic(factory.symbolic, args))
                end,
            )
        end
        pv = Tuple(safe_ustrip(get_unit(factory, p), param_vals[p]) for p in keys(factory.params))
        return NumericFunc((args...) -> generic(args..., pv...), Tuple(keys(factory.vars)), refs, unit)
    end

    cache_key = hash(tuple(sort(collect(pairs(param_vals)); by = x -> x.first)...))

    simplified, compiled = lock(_THERMO_FACTORY_LOCK) do
        get!(factory.cache, cache_key) do
            substitutions = Dict(
                v => safe_ustrip(get_unit(factory, p), param_vals[p]) for
                    (p, v) in factory.params
            )
            substituted = Symbolics.substitute(factory.symbolic, substitutions)
            reduced = Symbolics.simplify(Symbolics.expand(substituted))
            (reduced, compile_symbolic(reduced, collect(keys(factory.vars))))
        end
    end

    return SymbolicFunc(simplified, Tuple(keys(factory.vars)), refs, compiled, unit)
end

"""
    Base.show(io::IO, factory::ThermoFactory)

Compact string representation of a ThermoFactory.
"""
function Base.show(io::IO, factory::ThermoFactory)
    println(io, factory.symbolic, " [", dimension(get_unit(factory)), "]")
    if !isempty(factory.params)
        print(io, " ◆ params = ", join(keys(factory.params), ", "))
    end
    return if !isempty(factory.vars)
        print(io, " ◆ vars = ", join(keys(factory.vars), ", "))
    end
end

"""
    Base.show(io::IO, ::MIME"text/plain", factory::ThermoFactory)

Detailed string representation of a ThermoFactory.
"""
function Base.show(io::IO, ::MIME"text/plain", factory::ThermoFactory)
    println(io, "ThermoFactory:")
    print(io, "  Expression: ")
    println(io, factory.symbolic, " [", dimension(get_unit(factory)), "]")

    if !isempty(factory.params)
        print(io, "  Parameters: ")
        println(io, join(sort(collect(keys(factory.params))), ", "))
    end

    print(io, "  Variables: ")
    return print(io, join(sort(collect(keys(factory.vars))), ", "))
end

# ── SymbolicFunc convenience constructors ──────────────────────────────────────
# (placed after ThermoFactory since they use it)

# The keywords of `ThermoFactory`; every other keyword of a `SymbolicFunc`
# constructor goes to the factory's call: parameter values and reference values.
const _FACTORY_KEYWORDS = (:units, :output_unit)
_factory_kwargs(kwargs) = (; (k => v for (k, v) in pairs(kwargs) if k in _FACTORY_KEYWORDS)...)
_call_kwargs(kwargs) = (; (k => v for (k, v) in pairs(kwargs) if !(k in _FACTORY_KEYWORDS))...)

"""
    SymbolicFunc(sym::Symbol; units, output_unit, kwargs...) -> SymbolicFunc

Create a `SymbolicFunc` from a single symbol. `units` and `output_unit` are those
of [`ThermoFactory`](@ref); the other keywords are the values of its parameters
and the reference values of its variable.
"""
function SymbolicFunc(sym::Symbol; kwargs...)
    factory = ThermoFactory(sym, [sym]; _factory_kwargs(kwargs)...)
    return factory(; _call_kwargs(kwargs)...)
end

"""
    SymbolicFunc(expr::Expr, vars=[:T, :P, :t, :x, :y, :z]; units, output_unit, kwargs...) -> SymbolicFunc

Create a `SymbolicFunc` from an expression. `units` and `output_unit` are those
of [`ThermoFactory`](@ref); the other keywords are the values of its parameters
and the reference values of its variables.
"""
function SymbolicFunc(expr::Expr, vars = [:T, :P, :t, :x, :y, :z]; kwargs...)
    factory = ThermoFactory(expr, vars; _factory_kwargs(kwargs)...)
    return factory(; _call_kwargs(kwargs)...)
end

"""
    SymbolicFunc(x::Quantity) -> SymbolicFunc

Create a constant `SymbolicFunc` from a quantity.
"""
function SymbolicFunc(x::Quantity)
    x = uexpand(x)
    # The unit in plain numbers, whatever the number type of the value: a
    # constant being differentiated (a fitted log K) is a dual quantity, and its
    # `oneunit` would be one too.
    u = Quantity(1.0, dimension(x))
    factory = ThermoFactory(:c; units = [:c => u], output_unit = u)
    return factory(; c = x)
end

"""
    SymbolicFunc(x::Number) -> SymbolicFunc

Create a constant `SymbolicFunc` from a number (unitless).
"""
function SymbolicFunc(x::Number)
    factory = ThermoFactory(:c)
    return factory(; c = x)
end

# ── SymbolicFunc arithmetic ────────────────────────────────────────────────────

"""
    combine_symbolic(op, sf1::SymbolicFunc, sf2::SymbolicFunc)

Combine two SymbolicFuncs.
"""
function combine_symbolic(op, sf1::SymbolicFunc, sf2::SymbolicFunc)
    all_vars = Tuple(union(sf1.vars, sf2.vars))
    refs = merge(sf1.refs, sf2.refs)

    combined = op(sf1.symbolic, sf2.symbolic)
    simplified = Symbolics.simplify(Symbolics.expand(combined))
    compiled = compile_symbolic(simplified, all_vars)

    return SymbolicFunc(
        simplified, all_vars, refs, compiled, oneunit(op(sf1.unit, sf2.unit))
    )
end

"""
    combine_symbolic(op, sf::SymbolicFunc, x::Number)

Combine SymbolicFunc with scalar.
"""
function combine_symbolic(op, sf::SymbolicFunc, x::Number)
    combined = op(sf.symbolic, ustrip(x))
    simplified = Symbolics.simplify(Symbolics.expand(combined))
    compiled = compile_symbolic(simplified, collect(sf.vars))

    return SymbolicFunc(simplified, sf.vars, sf.refs, compiled, oneunit(op(sf.unit, x)))
end

"""
    combine_symbolic(op, x::Number, sf::SymbolicFunc)

Combine scalar with SymbolicFunc.
"""
function combine_symbolic(op, x::Number, sf::SymbolicFunc)
    combined = op(ustrip(x), sf.symbolic)
    simplified = Symbolics.simplify(Symbolics.expand(combined))
    compiled = compile_symbolic(simplified, collect(sf.vars))

    return SymbolicFunc(simplified, sf.vars, sf.refs, compiled, oneunit(op(x, sf.unit)))
end

"""
    apply_symbolic(op, sf::SymbolicFunc)

Apply unary operation.
"""
function apply_symbolic(op, sf::SymbolicFunc)
    combined = op(sf.symbolic)
    simplified = Symbolics.simplify(Symbolics.expand(combined))
    compiled = compile_symbolic(simplified, collect(sf.vars))

    return SymbolicFunc(simplified, sf.vars, sf.refs, compiled, oneunit(op(sf.unit)))
end

# Binary operations (SymbolicFunc op SymbolicFunc, SymbolicFunc op Number)
for op in (:+, :-, :*, :/, :^)
    @eval begin
        Base.$op(sf1::SymbolicFunc, sf2::SymbolicFunc) = combine_symbolic($op, sf1, sf2)
        Base.$op(sf::SymbolicFunc, x::Number) = combine_symbolic($op, sf, x)
        Base.$op(x::Number, sf::SymbolicFunc) = combine_symbolic($op, x, sf)
    end
end

# Unary
Base.:-(sf::SymbolicFunc) = apply_symbolic(-, sf)

# Math functions (dimensionless)
for f in [ADIM_MATH_FUNCTIONS; :sqrt; :abs]
    if isdefined(Base, f)
        @eval Base.$f(sf::SymbolicFunc) = apply_symbolic($f, sf)
    end
end

# ── Cross-type: NumericFunc × SymbolicFunc ─────────────────────────────────────
# vars are the union of both sides; nf args are extracted by index mapping,
# sf is called via NamedTuple kwargs (it ignores keys it does not own).
# refs: nf.refs take precedence over sf.refs.

function _cross_type_nf_sf(op, nf::NumericFunc, sf::SymbolicFunc, result_unit)
    new_vars = Tuple(union(nf.vars, sf.vars))
    new_refs = merge(sf.refs, nf.refs)   # nf.refs take precedence
    nf_idx = _var_indices(nf.vars, new_vars)
    function combined(args...)
        nf_vals = ntuple(i -> args[nf_idx[i]], length(nf.vars))
        kw = NamedTuple{new_vars}(args)
        return op(nf.compiled(nf_vals...), sf(; pairs(kw)...))
    end
    return NumericFunc(combined, new_vars, new_refs, result_unit)
end

function _cross_type_sf_nf(op, sf::SymbolicFunc, nf::NumericFunc, result_unit)
    new_vars = Tuple(union(nf.vars, sf.vars))
    new_refs = merge(sf.refs, nf.refs)   # nf.refs take precedence
    nf_idx = _var_indices(nf.vars, new_vars)
    function combined(args...)
        nf_vals = ntuple(i -> args[nf_idx[i]], length(nf.vars))
        kw = NamedTuple{new_vars}(args)
        return op(sf(; pairs(kw)...), nf.compiled(nf_vals...))
    end
    return NumericFunc(combined, new_vars, new_refs, result_unit)
end

for op in (:+, :-, :*, :/)
    @eval begin
        Base.$op(nf::NumericFunc, sf::SymbolicFunc) =
            _cross_type_nf_sf($op, nf, sf, oneunit(Base.$op(nf.unit, sf.unit)))
        Base.$op(sf::SymbolicFunc, nf::NumericFunc) =
            _cross_type_sf_nf($op, sf, nf, oneunit(Base.$op(sf.unit, nf.unit)))
    end
end

# ── Symbolic differentiation ───────────────────────────────────────────────────

"""
    derivative(sf::SymbolicFunc, var::Symbol) -> SymbolicFunc

Compute the analytical derivative of `sf` with respect to variable `var` using
Symbolics.jl and return a new `SymbolicFunc` for the result.

The output unit is inferred as `sf.unit / unit_of_var`, where `unit_of_var` is
read from `sf.refs` (defaulting to dimensionless if `var` has no ref entry).

# Examples
```julia
Cp = SymbolicFunc(:(a + b*T + c/T^2); a=25.0, b=8e-3, c=-1.5e5)
dCp_dT = derivative(Cp, :T)   # ∂Cp/∂T
dCp_dT(T = 500.0)
```
"""
function derivative(sf::SymbolicFunc, var::Symbol)
    var_sym = Symbolics.variable(var)
    dsym = Symbolics.derivative(sf.symbolic, var_sym)
    dsym_simplified = Symbolics.simplify(Symbolics.expand(dsym))
    compiled = compile_symbolic(dsym_simplified, collect(sf.vars))
    # Output unit: [sf.unit] / [unit of var]
    var_unit = haskey(sf.refs, var) ? oneunit(sf.refs[var]) : u"1"
    deriv_unit = oneunit(sf.unit / var_unit)
    return SymbolicFunc(dsym_simplified, sf.vars, sf.refs, compiled, deriv_unit)
end
