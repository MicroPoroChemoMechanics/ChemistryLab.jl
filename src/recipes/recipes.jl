# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

using OrderedCollections

# ── A recipe: materials, their masses, and the water ─────────────────────────

"""
    Recipe(binder...; w_b, binder_mass = 100u"g", additions = [],
           T = 293.15u"K", P = 1u"bar", balance = nothing)

A cementitious mix: the materials of the binder with their mass fractions
(`material => fraction` pairs, summing to one), the water/binder ratio `w_b`,
and optional `additions` by mass (`material => mass`, a salt or an admixture
outside the binder). `binder_mass` sets the scale: every amount the recipe
produces is for that mass of binder, 100 g by default, so that a phase mass read
off the answer is in grams per 100 g of binder.

Mass fractions that do not sum to one are refused, unless `balance` names the
material that takes the remainder; its given fraction is then ignored.

See [`budget`](@ref) for what enters the equilibrium, and
[`equilibrate_certified`](@ref), given the recipe and a system, to solve it.
"""
struct Recipe
    binder::Vector{Pair{Material, Float64}}
    additions::Vector{Pair{Material, Float64}}
    water_binder::Float64
    binder_mass::Float64
    T::Float64
    P::Float64
end
function Recipe(
        binder::Pair{Material, <:Real}...; w_b::Real,
        binder_mass = 100.0u"g", additions = Pair{Material, Float64}[],
        T = 293.15u"K", P = 1.0u"bar", balance = nothing,
    )
    isempty(binder) && throw(ArgumentError("Recipe: no binder material."))
    b = [m => float(f) for (m, f) in binder]
    if balance !== nothing
        k = findfirst(p -> first(p).name == balance, b)
        k === nothing && throw(ArgumentError("Recipe: `balance = \"$balance\"` names no material of the binder."))
        rest = sum(last(p) for (i, p) in enumerate(b) if i != k; init = 0.0)
        rest <= 1 || throw(ArgumentError("Recipe: the other materials already exceed the whole binder ($rest)."))
        b[k] = first(b[k]) => 1 - rest
    end
    tot = sum(last, b)
    isapprox(tot, 1; atol = 1.0e-6) || throw(
        ArgumentError(
            "Recipe: the binder mass fractions sum to $tot, not to one; give fractions of the " *
                "binder, or name the material that takes the remainder with `balance`."
        )
    )
    all(p -> last(p) >= 0, b) || throw(ArgumentError("Recipe: a negative mass fraction."))
    w_b >= 0 || throw(ArgumentError("Recipe: a negative water/binder ratio."))
    adds = [m => _in_unit(us"g", x) for (m, x) in additions]
    return Recipe(b, adds, float(w_b), _in_unit(us"g", binder_mass), _in_unit(us"K", T), _in_unit(us"Pa", P))
end

# Every material of a recipe with its mass in grams.
_material_masses(r::Recipe) = vcat([m => f * r.binder_mass for (m, f) in r.binder], r.additions)

# ── What enters the equilibrium, and what is kept aside ──────────────────────

"""
    budget(recipe, system; t = nothing) -> (; state, b, residual)

What the recipe puts into the equilibrium at time `t`, and what it keeps aside.

  - `state`: a [`ChemicalState`](@ref) of `system` holding the water and the
    reacted part of every mineral constituent, as species, at the recipe's
    temperature and pressure: its enthalpy is that of the reactants.
  - `b`: the element budget, that of `state` plus the reacted part of every
    oxide constituent through [`oxide_budget`](@ref).
  - `residual`: what is kept aside, one entry per constituent with its `mass`
    (g), `volume` (cm³) and `enthalpy` (J, of formation), `nothing` where a
    density or an enthalpy the constituent needs has no source, and `reason`:
    `:unreacted`, or `:not_in_system` for the reacted part of an oxide whose
    element the system has no primary for (the TiO2 or the P2O5 of a cement in a
    system without titanium or phosphorus).

The reacted part of a constituent is its mass times its effective extent (the
material's extent times its own, [`effective_extent`](@ref)). A mineral
constituent that reacts must be a species of `system`; one that does not need
not be.
"""
function budget(r::Recipe, cs::ChemicalSystem; t = nothing)
    st = ChemicalState(cs; T = r.T * u"K", P = r.P * u"Pa")
    haskey(cs.dict_species, "H2O@") || throw(ArgumentError("budget: the system has no H2O@ for the water of the recipe."))
    water = r.water_binder * r.binder_mass
    oxides = Pair{OrderedDict{String, Float64}, Float64}[]
    residual = NamedTuple[]
    for (m, mass) in _material_masses(r), c in m.constituents
        mc = mass * c.mass_fraction
        α = effective_extent(m, c, t)
        (0 <= α <= 1) || throw(ArgumentError("budget: the extent of $(c.name) in $(m.name) is $α at t = $t."))
        _add_reacted!(st, oxides, c, α * mc, m)
        (1 - α) * mc > 0 && push!(residual, _residue(c, m, (1 - α) * mc, r))
        _set_aside!(residual, c, m, α * mc, cs)
    end
    set_quantity!(st, "H2O@", moles(st, "H2O@") + water / ustrip(us"g/mol", cs.dict_species["H2O@"][:M]) * u"mol")
    b = Float64.(cs.SM.A) * ustrip.(us"mol", st.n)
    for (ox, mass) in oxides
        b .+= oxide_budget(ox, cs.SM.primaries; mass = mass * u"g")
    end
    return (; state = st, b, residual)
end

function _add_reacted!(st, oxides, c::MineralConstituent, mass, m)
    mass > 0 || return nothing
    sym = symbol(c.species)
    haskey(st.system.dict_species, sym) || throw(
        ArgumentError(
            "budget: $(c.name) of $(m.name) reacts, and $sym is not a species of the system. " *
                "Add it to the species list, or describe the constituent by its oxides."
        )
    )
    n = mass / ustrip(us"g/mol", c.species[:M])
    set_quantity!(st, sym, moles(st, sym) + n * u"mol")
    return nothing
end
function _add_reacted!(st, oxides, c::OxideConstituent, mass, m)
    mass > 0 || return nothing
    prim = st.system.SM.primaries
    held = OrderedDict(k => v for (k, v) in c.oxides if _representable(k, prim))
    isempty(held) || push!(oxides, held => mass)
    return nothing
end

# Whether the primaries of the system can hold the element of an oxide.
function _representable(ox, prim)
    return try
        primary_decomposition(Species(ox), prim)
        true
    catch
        false
    end
end

# The reacted part of an oxide the system has no element for is kept aside with
# the residue, said so by `reason`, rather than refused or dropped.
_set_aside!(residual, ::MineralConstituent, m, mass, cs) = nothing
function _set_aside!(residual, c::OxideConstituent, m, mass, cs)
    mass > 0 || return nothing
    prim = cs.SM.primaries
    for (ox, f) in c.oxides
        (f > 0 && !_representable(ox, prim)) || continue
        push!(
            residual,
            (;
                material = m.name, constituent = "$(c.name): $ox", mass = mass * f,
                volume = c.density === nothing ? nothing : mass * f / c.density,
                enthalpy = nothing, reason = :not_in_system,
            ),
        )
    end
    return nothing
end

function _residue(c::MineralConstituent, m, mass, r)
    n = mass / ustrip(us"g/mol", c.species[:M])
    V = _species_value(c.species, :V⁰, us"cm^3/mol", r)
    H = _species_value(c.species, :ΔₐH⁰, us"J/mol", r)
    return (;
        material = m.name, constituent = c.name, mass,
        volume = V === nothing ? nothing : n * V, enthalpy = H === nothing ? nothing : n * H,
        reason = :unreacted,
    )
end
function _residue(c::OxideConstituent, m, mass, r)
    return (;
        material = m.name, constituent = c.name, mass,
        volume = c.density === nothing ? nothing : mass / c.density,
        enthalpy = c.enthalpy === nothing ? nothing : mass * c.enthalpy,
        reason = :unreacted,
    )
end

# A thermodynamic property of a species at the recipe's T and P, or `nothing`
# when the species does not carry it; evaluated as the volumes of a state are.
function _species_value(sp, key, unit, r)
    haskey(sp, key) || return nothing
    return _in_unit(unit, sp[key](T = r.T * u"K", P = r.P * u"Pa"; unit = true))
end

# ── The answer, with its residue ─────────────────────────────────────────────

"""
    RecipeState

The certified equilibrium of a [`Recipe`](@ref) at a time, with what the recipe
kept aside. Fields: `state` (the equilibrated [`ChemicalState`](@ref)),
`certificate`, `recipe`, `t`, `initial` (the reactants, from [`budget`](@ref)),
`b`, `residual` and `model` (the activity model of the solve).

Read it with [`volume`](@ref), [`porosity`](@ref), [`phase_masses`](@ref),
[`bound_water`](@ref), [`pore_solution`](@ref) and [`enthalpy`](@ref), which all
count the residue where it belongs.
"""
struct RecipeState{S <: ChemicalState, C, M}
    state::S
    certificate::C
    recipe::Recipe
    t::Any
    initial::S
    b::Vector{Float64}
    residual::Vector{NamedTuple}
    model::M
end

"""
    equilibrate_certified(recipe::Recipe, system::ChemicalSystem; t = nothing,
                          model = DiluteSolutionModel(), start = nothing, kwargs...)
        -> (RecipeState, certificate)

The certified equilibrium of `recipe` at time `t` in `system`: its
[`budget`](@ref) solved by [`equilibrate_certified`](@ref), from `start` (a
state of `system`, typically the answer at an earlier time) when given, from the
reactants otherwise. The other keywords go to `equilibrate_certified`.
"""
function equilibrate_certified(
        r::Recipe, cs::ChemicalSystem; t = nothing, model::AbstractActivityModel = DiluteSolutionModel(),
        start::Union{Nothing, ChemicalState} = nothing, kwargs...,
    )
    bud = budget(r, cs; t)
    from = start === nothing ? bud.state : start
    eq, cert = equilibrate_certified(from; model, b = bud.b, kwargs...)
    return RecipeState(eq, cert, r, t, bud.state, bud.b, bud.residual, model), cert
end

"""
    residual_mass(rs::RecipeState) -> Float64

The mass (g) of what has not reacted.
"""
residual_mass(rs::RecipeState) = sum(x.mass for x in rs.residual; init = 0.0)

# The residue's volume and enthalpy, and the constituents that lack one.
function _residual_sum(rs, field)
    missing_ = [x.constituent for x in rs.residual if getfield(x, field) === nothing]
    return (; value = isempty(missing_) ? sum(getfield(x, field) for x in rs.residual; init = 0.0) : NaN, missing = missing_)
end

"""
    volume(rs::RecipeState) -> NamedTuple

The volumes (cm³) of the paste: `liquid`, `gas` and `solid` of the equilibrium,
`residual` of what has not reacted (`NaN` when a constituent has no density, see
`missing`), and `total`.
"""
function volume(rs::RecipeState)
    V = volume(rs.state)
    cm3(x) = _in_unit(us"cm^3", x)
    res = _residual_sum(rs, :volume)
    return (;
        liquid = cm3(V.liquid), gas = cm3(V.gas), solid = cm3(V.solid), residual = res.value,
        total = cm3(V.total) + res.value, missing = res.missing,
    )
end

"""
    porosity(rs::RecipeState) -> NamedTuple

The porosity of a sealed paste, as fractions of its initial volume (the
reactants, the water and the residue): `liquid` (filled by the pore solution and
any gas), `void` (the chemical shrinkage, what the reactions did not refill)
and `total`. The residue is counted in the solid of both volumes, which is what
makes the porosity of a paste with an unreacted part the right one.
"""
function porosity(rs::RecipeState)
    res = _residual_sum(rs, :volume).value
    V0 = _in_unit(us"cm^3", volume(rs.initial).total) + res
    V = volume(rs.state)
    liquid = _in_unit(us"cm^3", V.liquid + V.gas) / V0
    void = max((V0 - (_in_unit(us"cm^3", V.total) + res)) / V0, 0.0)
    return (; liquid, void, total = liquid + void)
end

"""
    volume_fractions(rs::RecipeState; void_key = "void") -> OrderedDict{String, Float64}

The volume fraction of every species of the equilibrium and of every unreacted
constituent (under `"unreacted <name>"`), relative to the initial volume of the
paste: the reactants, the water and the residue, the reference of
[`porosity`](@ref). The chemical-shrinkage void closes the sum under `void_key`,
so the fractions add up to one and `void` is `porosity(rs).void`.

A constituent whose residue has no sourced density leaves every fraction
undefined, and is refused by name rather than counted as nothing.
[`volume_fractions(state, groups)`](@ref) groups the species of the state alone.
"""
function volume_fractions(rs::RecipeState; void_key::AbstractString = "void")
    res = _residual_sum(rs, :volume)
    isempty(res.missing) || throw(
        ArgumentError(
            "volume_fractions: no sourced density for the unreacted part of " *
                join(res.missing, ", ") * ", so its volume, and every fraction, is unknown."
        )
    )
    V_initial = _in_unit(us"cm^3", volume(rs.initial).total)
    V0 = V_initial + res.value
    out = OrderedDict{String, Float64}()
    for (k, f) in volume_fractions(rs.state; reference = rs.initial, void_key)
        k == void_key && continue
        out[k] = f * V_initial / V0
    end
    for x in rs.residual
        key = "unreacted " * x.constituent
        out[key] = get(out, key, 0.0) + x.volume / V0
    end
    void = 1.0 - sum(values(out))
    void < -1.0e-10 && @warn "volume expanded beyond the initial volume; " *
        "the sealed-volume convention does not apply" excess = -void
    out[void_key] = max(void, 0.0)
    return out
end

"""
    phase_masses(rs::RecipeState; min_mass = 1e-6) -> OrderedDict{String, Float64}

The mass (g, for the recipe's binder mass, so g per 100 g of binder by default)
of every solid of the equilibrium above `min_mass`, largest first, then of each
unreacted constituent, under `"unreacted <name>"`.
"""
function phase_masses(rs::RecipeState; min_mass::Real = 1.0e-6)
    cs = rs.state.system
    out = OrderedDict{String, Float64}()
    solids = [(symbol(cs.species[i]), _in_unit(us"g", mass(rs.state, cs.species[i]))) for i in cs.idx_crystal]
    for (s, m) in sort(solids; by = last, rev = true)
        m > min_mass && (out[s] = m)
    end
    for x in rs.residual
        out["unreacted " * x.constituent] = get(out, "unreacted " * x.constituent, 0.0) + x.mass
    end
    return out
end

"""
    bound_water(rs::RecipeState; window = nothing, windows = nothing, min_mass = 1e-6)
        -> Float64

The bound water per gram of binder (g/g): the water the solids of the equilibrium
would lose on ignition ([`ignition_loss`](@ref)), plus that of an unreacted
mineral which holds some, over the binder mass.

With `window = (T₁, T₂)` (kelvin, as plain numbers or as quantities) and the
decomposition `windows` of the solids ([`DecompositionWindow`](@ref)), it is
instead the water released between the two temperatures, which is what a
thermogravimetric reading over that range weighs: a source that reports bound
water between 105 °C and 550 °C is compared with the same range. Every solid
holding more than `min_mass` grams of water, an unreacted mineral included,
needs a window releasing water; one without is refused by name, since its water
cannot be placed inside or outside the range.
"""
function bound_water(rs::RecipeState; window = nothing, windows = nothing, min_mass::Real = 1.0e-6)
    if window === nothing
        windows === nothing || throw(ArgumentError("bound_water: `windows` are read only with a temperature `window`."))
        w = _in_unit(us"g", ignition_loss(rs.state).water)
        for (_, g) in _unreacted_water(rs)
            w += g
        end
        return w / rs.recipe.binder_mass
    end
    held = vcat(
        Pair{String, Float64}[p.first => _in_unit(us"g", p.second) for p in bound_water_per_phase(rs.state)],
        _unreacted_water(rs),
    )
    windows === nothing && throw(
        ArgumentError("bound_water: a temperature window needs the decomposition `windows` of the solids.")
    )
    T1, T2 = _kelvin(first(window)), _kelvin(last(window))
    T1 < T2 || throw(ArgumentError("bound_water: the window ($T1 K, $T2 K) is empty."))
    _check_fractions(windows)
    water_windows = [w for w in windows if w.releases === :water]
    released = 0.0
    unplaced = String[]
    for (phase, g) in held
        ws = [w for w in water_windows if w.phase == phase]
        if isempty(ws)
            g > min_mass && push!(unplaced, phase)
            continue
        end
        released += g * sum(w.fraction * (released_fraction(w, T2) - released_fraction(w, T1)) for w in ws)
    end
    isempty(unplaced) || throw(
        ArgumentError(
            "bound_water: no decomposition window releasing water for " * join(unique(unplaced), ", ") *
                ", which hold water; give one per solid (see `DecompositionWindow`)."
        )
    )
    return released / rs.recipe.binder_mass
end

_kelvin(T::Real) = float(T)
_kelvin(T::DynamicQuantities.AbstractQuantity) = ustrip(us"K", T)

# The water (g) of each unreacted mineral that holds some, by species symbol.
function _unreacted_water(rs::RecipeState)
    out = Pair{String, Float64}[]
    Mw = ustrip(us"g/mol", Species("H2O")[:M])
    for (m, mass) in _material_masses(rs.recipe), c in m.constituents
        c isa MineralConstituent || continue
        α = effective_extent(m, c, rs.t)
        h = Float64(get(atoms(c.species), :H, 0))
        w = (1 - α) * mass * c.mass_fraction / ustrip(us"g/mol", c.species[:M]) * h / 2 * Mw
        w > 0 && push!(out, symbol(c.species) => w)
    end
    return out
end

"""
    pore_solution(rs::RecipeState) -> NamedTuple

The pore solution: `pH` in the activity convention of the solve's model, and
`elements`, the total molality (mol per kg of water) of each element dissolved,
all aqueous species counted.
"""
function pore_solution(rs::RecipeState)
    cs = rs.state.system
    iw = findfirst(==("H2O@"), [symbol(s) for s in cs.species])
    kgw = _in_unit(us"kg", mass(rs.state, cs.species[iw]))
    el = OrderedDict{Symbol, Float64}()
    for i in cs.idx_solutes
        n = ustrip(us"mol", rs.state.n[i])
        for (e, k) in atoms(cs.species[i])
            (e === :H || e === :O) && continue
            el[e] = get(el, e, 0.0) + k * n / kgw
        end
    end
    return (; pH = pH(rs.state, rs.model), elements = el)
end

"""
    enthalpy(rs::RecipeState) -> Float64

The enthalpy (J) of the paste, the equilibrium's plus the residue's; `NaN` when
an unreacted constituent has no sourced enthalpy of formation, so that a heat
computed from it cannot pass for complete.
"""
enthalpy(rs::RecipeState) = _in_unit(us"J", enthalpy(rs.state)) + _residual_sum(rs, :enthalpy).value

"""
    heat_release(rs1::RecipeState, rs2::RecipeState) -> Float64

The heat (J) a paste releases from the state `rs1` to the state `rs2` of the
same recipe, at the same temperature and pressure: the fall of its enthalpy,
`-(H₂ - H₁)`. The residue counts by what changed between the two. A constituent
set aside with the same mass in both (an inert crystal, an oxide the system has
no primary for) adds nothing, whether or not its enthalpy is sourced. One whose
unreacted mass changed needs a sourced enthalpy of formation, and the heat is
`NaN` without it.
"""
function heat_release(a::RecipeState, b::RecipeState)
    (a.state.T[1] == b.state.T[1] && a.state.P[1] == b.state.P[1]) || throw(
        ArgumentError("heat_release: the two states are at different temperatures or pressures; the heat is that of an isothermal, isobaric change."),
    )
    q = -(_in_unit(us"J", enthalpy(b.state)) - _in_unit(us"J", enthalpy(a.state)))
    ra = Dict(x.constituent => x for x in a.residual)
    rb = Dict(x.constituent => x for x in b.residual)
    for k in union(keys(ra), keys(rb))
        xa, xb = get(ra, k, nothing), get(rb, k, nothing)
        ma = xa === nothing ? 0.0 : xa.mass
        mb = xb === nothing ? 0.0 : xb.mass
        ma == mb && continue
        Ha = xa === nothing ? 0.0 : xa.enthalpy
        Hb = xb === nothing ? 0.0 : xb.enthalpy
        (Ha === nothing || Hb === nothing) && return NaN
        q -= Hb - Ha
    end
    return q
end

function Base.show(io::IO, rs::RecipeState)
    return print(
        io, "RecipeState(t = ", rs.t, ", certified ", rs.certificate.optimal, ", ",
        length(rs.residual), " unreacted parts, ", round(residual_mass(rs); sigdigits = 4), " g)",
    )
end
