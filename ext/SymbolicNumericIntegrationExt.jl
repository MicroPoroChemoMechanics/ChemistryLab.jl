# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

module SymbolicNumericIntegrationExt

using ChemistryLab
import ChemistryLab: add_thermo_model, extract_vars_params
using Symbolics: Num, variable, simplify, expand, parse_expr_to_symbolic, toexpr
using SymbolicNumericIntegration: integrate, terms

function ChemistryLab.add_thermo_model(model_name, Cpexpr::Expr, units = nothing)
    vars, params = extract_vars_params(Cpexpr, [:T])
    var_sym_dict = Dict{Symbol, Num}(v => variable(v) for v in vars)
    param_sym_dict = Dict{Symbol, Num}(p => variable(p) for p in params)
    all_symbols = merge(var_sym_dict, param_sym_dict)
    T = var_sym_dict[:T]

    Cp = simplify(expand(parse_expr_to_symbolic(Cpexpr, all_symbols)))
    H = integrate(Cp, T; symbolic = true, detailed = false)

    integS = sum(terms(Cp) ./ T)
    S = integrate(integS, T; symbolic = true, detailed = false)

    G = integrate(-S, T; symbolic = true, detailed = false)

    # Build a Dict{Symbol, Any} so we can add String default units without
    # triggering a type-mismatch push! when the caller passed Quantity values.
    units_dict = if isnothing(units)
        nothing
    else
        d = Dict{Symbol, Any}(units)
        get!(d, :Cp, "J/mol/K")
        get!(d, :S, "J/mol/K")
        get!(d, :H, "J/mol")
        get!(d, :G, "J/mol")
        collect(pairs(d))   # convert back to a vector of Pairs for ThermoFactory
    end

    dict_model = if isnothing(units_dict)
        Dict(:Cp => Cpexpr, :H => _model_expr(H), :S => _model_expr(S), :G => _model_expr(G))
    else
        Dict(
            :Cp => Cpexpr, :H => _model_expr(H), :S => _model_expr(S), :G => _model_expr(G),
            :units => units_dict,
        )
    end

    # Registered as any model is, under the lock its readers take.
    return add_thermo_model(model_name, dict_model)
end

# The expression of an integrated function, as the model reader reads it. From
# SymbolicUtils 4.53 on, `toexpr` writes an integer power `x^n` as
# `Base.literal_pow(^, x, Val(n))`, the form Julia lowers `x^2` to, which the
# reader takes for a call of unknown type: the substitution of the parameters
# then fails (`symtype(ex) <: Number`) on every function holding a power, H and
# G of a linear Cp. The power is written back as `x ^ n`.
_model_expr(f) = _plain_powers(toexpr(f))

_plain_powers(ex) = ex
function _plain_powers(ex::Expr)
    args = map(_plain_powers, ex.args)
    if ex.head === :call && length(args) == 4 && args[1] === Base.literal_pow && args[2] === (^) &&
            args[4] isa Expr && args[4].head === :call && length(args[4].args) == 2 && args[4].args[1] === Val
        return Expr(:call, ^, args[3], args[4].args[2])
    end
    return Expr(ex.head, args...)
end

end # module SymbolicNumericIntegrationExt
