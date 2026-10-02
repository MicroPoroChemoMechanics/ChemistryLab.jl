# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

using DataFrames: DataFrame

# ── Sequences of equilibria ───────────────────────────────────────────────────

"""
    ProcessResult

A process as a sequence of [`RecipeState`](@ref)s, one per value of its
parameter: `parameter` (the times of [`hydrate`](@ref), the replacement
fractions of [`blend`](@ref), the amounts of [`titrate`](@ref), the steps of
[`leach`](@ref)), `label` (what the parameter is) and `states`. Each state was
solved from the previous one. [`process_table`](@ref) tabulates it.
"""
struct ProcessResult
    label::Symbol
    parameter::Vector{Any}
    states::Vector{RecipeState}
end
Base.length(p::ProcessResult) = length(p.states)
Base.getindex(p::ProcessResult, i) = p.states[i]
Base.firstindex(p::ProcessResult) = firstindex(p.states)
Base.lastindex(p::ProcessResult) = lastindex(p.states)
Base.keys(p::ProcessResult) = keys(p.states)
Base.iterate(p::ProcessResult, s...) = iterate(p.states, s...)

"""
    hydrate(recipe, system, times; model = DiluteSolutionModel(), kwargs...) -> ProcessResult

The recipe at each of `times` (days, or `Quantity`s), its extents read at that
time, each equilibrium started from the previous answer. The keywords go to
[`equilibrate_certified`](@ref).
"""
function hydrate(r::Recipe, cs::ChemicalSystem, times; model::AbstractActivityModel = DiluteSolutionModel(), kwargs...)
    states = RecipeState[]
    prev = nothing
    for t in times
        rs, _ = equilibrate_certified(r, cs; t, model, start = prev, kwargs...)
        push!(states, rs)
        prev = rs.state
    end
    return ProcessResult(:time, collect(Any, times), states)
end

"""
    blend(base, scm, fractions, system; t = nothing, model = DiluteSolutionModel(), kwargs...)
        -> ProcessResult

`base` with the fraction `f` of its binder replaced by the material `scm`, for
each `f` of `fractions`: every other material of the binder keeps its share of
what remains, the water/binder ratio and the binder mass stay those of `base`.
"""
function blend(
        base::Recipe, scm::Material, fractions, cs::ChemicalSystem; t = nothing,
        model::AbstractActivityModel = DiluteSolutionModel(), kwargs...,
    )
    states = RecipeState[]
    prev = nothing
    for f in fractions
        0 <= f <= 1 || throw(ArgumentError("blend: a replacement fraction is in [0, 1]; got $f."))
        binder = [m => (1 - f) * x for (m, x) in base.binder if m !== scm]
        f > 0 && push!(binder, scm => f)
        r = Recipe(
            binder...; w_b = base.water_binder, binder_mass = base.binder_mass * u"g",
            additions = base.additions, T = base.T * u"K", P = base.P * u"Pa",
        )
        rs, _ = equilibrate_certified(r, cs; t, model, start = prev, kwargs...)
        push!(states, rs)
        prev = rs.state
    end
    return ProcessResult(:fraction, collect(Any, fractions), states)
end

"""
    titrate(rs, species, amounts; kwargs...) -> ProcessResult

The paste `rs` with each of `amounts` (mol, or `Quantity`s) of `species` added to
its element budget, the residue unchanged, each equilibrium started from the
previous one: the gas of a carbonation, the salt of a chloride ingress. `species`
is a symbol of the system, or a formula (`"NaCl"`) whose elements its primaries
hold, which need not be a species of it. [`carbonate`](@ref) and
[`add_salt`](@ref) are this with a chosen species.
"""
function titrate(rs::RecipeState, species::AbstractString, amounts; kwargs...)
    cs = rs.state.system
    j = findfirst(==(species), [symbol(s) for s in cs.species])
    col = j === nothing ? _formula_column(cs, species) : Float64.(cs.SM.A[:, j])
    states = RecipeState[]
    prev = rs.state
    for a in amounts
        n = _in_unit(us"mol", a)
        n >= 0 || throw(ArgumentError("titrate: a negative amount of $species."))
        b = rs.b .+ n .* col
        eq, cert = equilibrate_certified(prev; model = rs.model, b, kwargs...)
        push!(states, RecipeState(eq, cert, rs.recipe, rs.t, rs.initial, b, rs.residual, rs.model))
        prev = eq
    end
    return ProcessResult(Symbol(species), collect(Any, amounts), states)
end

# A formula that is not a species of the system, in its primaries.
function _formula_column(cs, formula)
    sp = try
        Species(formula)
    catch
        throw(ArgumentError("titrate: $formula is neither a species of the system nor a formula."))
    end
    return primary_decomposition(sp, cs.SM.primaries)
end

"""
    carbonate(rs, amounts; kwargs...) -> ProcessResult

[`titrate`](@ref) with dissolved carbon dioxide, `CO2@`.
"""
carbonate(rs::RecipeState, amounts; kwargs...) = titrate(rs, "CO2@", amounts; kwargs...)

"""
    add_salt(rs, salt, amounts; kwargs...) -> ProcessResult

[`titrate`](@ref) with the salt `salt`: its formula (`"NaCl"`, `"CaCl2"`), or a
solid salt the system declares. A salt enters as a whole, since one ion alone
would make the budget charged.
"""
add_salt(rs::RecipeState, salt::AbstractString, amounts; kwargs...) = titrate(rs, salt, amounts; kwargs...)

"""
    leach(rs, steps; renewal = nothing, kwargs...) -> ProcessResult

`steps` renewals of the pore solution: at each, the whole aqueous phase of the
last equilibrium is removed and replaced by `renewal` of pure water (g; by
default the water of the recipe), and the paste is equilibrated again with what
remains. The residue is unchanged. Each budget depends on the previous answer,
which is what makes this a sequence rather than a sweep.
"""
function leach(rs::RecipeState, steps::Integer; renewal = nothing, kwargs...)
    w = renewal === nothing ? rs.recipe.water_binder * rs.recipe.binder_mass : _in_unit(us"g", renewal)
    states = RecipeState[]
    prev = rs.state
    for _ in 1:steps
        start, b = _renewal(prev, w)
        eq, cert = equilibrate_certified(start; model = rs.model, b, kwargs...)
        push!(states, RecipeState(eq, cert, rs.recipe, rs.t, rs.initial, b, rs.residual, rs.model))
        prev = eq
    end
    return ProcessResult(:step, collect(Any, 1:steps), states)
end

# The state `prev` with its aqueous phase replaced by `w` grams of pure water, and
# its budget. Read in the system of `prev`, which holds a second instance of a
# phase declared `instances = :auto` once an earlier step has split it.
function _renewal(prev::ChemicalState, w)
    cs = prev.system
    iw = findfirst(==("H2O@"), [symbol(s) for s in cs.species])
    Mw = ustrip(us"g/mol", cs.species[iw][:M])
    n0 = ustrip.(us"mol", prev.n)
    # In the number type of the previous answer and of the renewal.
    n = collect(promote_type(eltype(n0), typeof(w / Mw)), n0)
    for i in cs.idx_aqueous
        n[i] = 0.0
    end
    n[iw] = w / Mw
    start = ChemicalState(cs; T = prev.T[1], P = prev.P[1], n = n .* u"mol")
    return start, Float64.(cs.SM.A) * n
end

"""
    process_table(p::ProcessResult; phases = nothing) -> DataFrame

One row per state of `p`: the parameter, whether it certified, the pH, the
porosity, the bound water (g/g of binder) and the mass (g per the recipe's binder
mass) of each of `phases`, by default every solid that appears in some state.
"""
function process_table(p::ProcessResult; phases = nothing)
    masses = [phase_masses(rs) for rs in p.states]
    names = phases === nothing ?
        unique(k for m in masses for k in keys(m) if !startswith(k, "unreacted ")) : collect(String, phases)
    df = DataFrame(p.label => p.parameter)
    df.certified = [rs.certificate.optimal for rs in p.states]
    df.pH = [pH(rs.state, rs.model) for rs in p.states]
    df.porosity = [porosity(rs).total for rs in p.states]
    df.bound_water = [bound_water(rs) for rs in p.states]
    for nm in names
        df[!, nm] = [get(m, nm, 0.0) for m in masses]
    end
    return df
end
