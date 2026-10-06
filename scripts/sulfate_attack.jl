# =============================================================================
#  sulfate_attack.jl — a Portland cement paste and sodium sulfate, at equilibrium
#
#  One cement, two papers. Lothenbach et al. (2010) exposed mortars of a
#  laboratory CEM I to 4 and 44 g/L of Na2SO4 and measured the sulfur and the
#  calcium of the paste by SEM-EDS against the depth; Schmidt et al. (2008)
#  equilibrated crushed pastes of the same cement, with 0, 5 and 25 % of
#  limestone, with Na2SO4 solutions at 8 and 20 °C, and measured the solutions,
#  the phases and the thaumasite. This file computes both with Cemdata18: the
#  paste titrated by the solution, the axis of the 0D model of Lothenbach et al.,
#  and the closed batches of Schmidt et al. The page
#  docs/src/examples/sulfate_attack.md compares them. Nothing is fitted.
#
#  ASSUMED, as both papers' own calculations assume: the cement fully hydrated.
# =============================================================================

using ChemistryLab
using DynamicQuantities
using OptimaSolver

const SA_DB_SUBSTANCES = build_species(datapath("cemdata18-thermofun.json"); verbose = false)
const SA_DB = Dict(symbol(s) => s for s in SA_DB_SUBSTANCES)
sa_table(key, name) = literature_table(key, name)
sa_value(key, name) = ustrip(literature_value(key, name))

# The phases of a Portland paste, with what sulfate attack adds: thaumasite and
# gypsum, which the list holds as a reactant and which may precipitate.
const SA_PHASES = "Portland paste (Lothenbach and Winnefeld 2006)"
sa_system(; thaumasite = true) =
    phase_list_system(SA_PHASES, SA_DB_SUBSTANCES; add = thaumasite ? ["thaumasite"] : String[])

const SA_CLINKER = ("alite" => "C3S", "belite" => "C2S", "aluminate" => "C3A", "ferrite" => "C4AF")

"""
    sa_cement() -> Material

The laboratory CEM I of Lothenbach et al. (2010), Table 1, which Schmidt et al.
(2008) used as their OPC: its normative phases, the magnesium, alkalis and
sulfate the clinker carries, free lime, calcite, gypsum and the alkali
sulfates, all reacted.
"""
function sa_cement()
    t = sa_table("Lothenbach2010sulfate", "normative_phases")
    f(p) = ustrip(t.percent[findfirst(==(p), t.phase)]) / 100
    c = AbstractConstituent[MineralConstituent(SA_DB[r]; mass_fraction = f(p), name = p) for (p, r) in SA_CLINKER]
    minors = Dict(ox => f(ox) for ox in ("MgO", "K2O", "Na2O", "SO3"))
    total = sum(values(minors))
    push!(c, OxideConstituent("clinker minors", Dict(k => v / total for (k, v) in minors); mass_fraction = total))
    push!(c, OxideConstituent("free lime", Dict("CaO" => 1.0); mass_fraction = f("CaO_free")))
    push!(c, MineralConstituent(SA_DB["Cal"]; mass_fraction = f("CaCO3")))
    push!(c, MineralConstituent(SA_DB["Gp"]; mass_fraction = f("CaSO4.2H2O")))
    push!(c, MineralConstituent(SA_DB["K2SO4"]; mass_fraction = f("K2SO4")))
    push!(c, MineralConstituent(SA_DB["Na2SO4"]; mass_fraction = f("Na2SO4")))
    return Material("CEM I (Lothenbach et al. 2010)", :cement; constituents = c)
end

"""The limestone of Schmidt et al. (2008), taken as calcite."""
sa_limestone() = Material("limestone", :other; constituents = AbstractConstituent[MineralConstituent(SA_DB["Cal"]; mass_fraction = 1.0)])

"""
    sa_recipe(; limestone = 0.0, w_b, T) -> Recipe

The cement with `limestone` (fraction of the binder) replaced by limestone, at
the water/binder ratio `w_b` and the temperature `T` (K).
"""
function sa_recipe(; limestone = 0.0, w_b, T)
    binder = limestone > 0 ? (sa_cement() => 1 - limestone, sa_limestone() => limestone) : (sa_cement() => 1.0,)
    return Recipe(binder...; w_b, T = T * u"K")
end

const SA_MODEL = cemdata18_activity_model(:KOH)

"""
    sa_paste(; limestone = 0.0, w_b, T, thaumasite = true) -> RecipeState

The fully hydrated paste, certified.
"""
# The pastes already computed, by binder, w/b, temperature and phase set: the
# batches of one binder share theirs.
const _SA_PASTES = Dict{Tuple{Float64, Float64, Float64, Bool}, Any}()

function sa_paste(; limestone = 0.0, w_b, T, thaumasite = true)
    return get!(() -> _sa_paste(limestone, w_b, T, thaumasite), _SA_PASTES, (Float64(limestone), Float64(w_b), Float64(T), thaumasite))
end

function _sa_paste(limestone, w_b, T, thaumasite)
    cs = sa_system(; thaumasite)
    # At another temperature than 20 °C, from the answer at 20 °C: a cold start
    # on a cement does not always find its way.
    start = T == 293.15 ? nothing : sa_paste(; limestone, w_b, T = 293.15, thaumasite).state
    rs, cert = equilibrate_certified(sa_recipe(; limestone, w_b, T), cs; model = SA_MODEL, start)
    cert.optimal || error("sa_paste: the paste did not certify.")
    return rs
end

# The budget of `V` mL of a solution of `c` mol/L of Na2SO4: the salt and the
# water it is dissolved in, its density taken as one.
function _sa_solution_budget(cs, V, c)
    col(f) = (
        j = findfirst(==(f), [symbol(s) for s in cs.species]);
        j === nothing ? primary_decomposition(Species(f), cs.SM.primaries) : Float64.(cs.SM.A[:, j])
    )
    Mw = ustrip(us"g/mol", cs.dict_species["H2O@"][:M])
    return c * V / 1000 .* col("Na2SO4") .+ (V / Mw) .* col("H2O@")
end

"""
    sa_titration(rs, c, volumes) -> Vector{NamedTuple}

The paste `rs` mixed with each of `volumes` (mL per 100 g of binder, increasing)
of a solution of `c` mol/L of Na2SO4, each equilibrium started from the
previous one: the volume, whether it certified, the masses of the solids
(g per 100 g of binder), the SO3 and CaO of the solids as mass fractions of their
oxides, and the pH.
"""
function sa_titration(rs, c, volumes)
    cs = rs.state.system
    prev = rs.state
    out = NamedTuple[]
    for V in volumes
        b = rs.b .+ _sa_solution_budget(cs, V, c)
        eq, cert = equilibrate_certified(prev; model = rs.model, b)
        push!(out, (; V, certified = cert.optimal, _sa_solids(eq)..., pH = pH(eq, rs.model)))
        prev = eq
    end
    return out
end

# The solids of an equilibrium: the mass of each (g per 100 g of binder, the
# recipe's scale), and the SO3 and CaO of all of them over their oxides, water
# left out, the way an elemental analysis of the solid reports it.
function _sa_solids(eq)
    cs = eq.system
    n = ustrip.(us"mol", eq.n)
    solid = setdiff(eachindex(cs.species), cs.idx_aqueous)
    M(f) = ustrip(us"g/mol", Species(f)[:M])
    el(e) = sum(n[i] * Float64(get(atoms(cs.species[i]), e, 0)) for i in solid)
    oxides = Dict(
        "CaO" => el(:Ca) * M("CaO"), "SiO2" => el(:Si) * M("SiO2"), "Al2O3" => el(:Al) / 2 * M("Al2O3"),
        "Fe2O3" => el(:Fe) / 2 * M("Fe2O3"), "MgO" => el(:Mg) * M("MgO"), "SO3" => el(:S) * M("SO3"),
        "CO2" => el(:C) * M("CO2"), "Na2O" => el(:Na) / 2 * M("Na2O"), "K2O" => el(:K) / 2 * M("K2O"),
    )
    total = sum(values(oxides))
    mass(s) = (j = findfirst(x -> symbol(x) == s, cs.species); j === nothing ? 0.0 : n[j] * ustrip(us"g/mol", cs.species[j][:M]))
    ss(name) = (t = solid_solution_totals(eq, name); t === nothing ? 0.0 : 1000 * ustrip(t.mass))
    solids = sum(n[i] * ustrip(us"g/mol", cs.species[i][:M]) for i in solid)
    return (;
        SO3 = oxides["SO3"] / total, CaO = oxides["CaO"] / total, solids,
        gypsum = mass("Gp"), ettringite = mass("ettringite"), thaumasite = mass("thaumasite"),
        portlandite = mass("Portlandite"), calcite = mass("Cal"), monocarbonate = mass("monocarbonate"),
        afm_so4 = _sa_members(eq, ("monosulphate12", "C4AH13")), csh = ss("CSHQ"),
    )
end

# The mass (g) of the members `names` of a solid solution, every instance of
# it counted: a solution split in two names its second instance's members
# `name#2`.
function _sa_members(eq, names)
    cs = eq.system
    n = ustrip.(us"mol", eq.n)
    return sum(
        (
            n[j] * ustrip(us"g/mol", cs.species[j][:M]) for j in eachindex(n)
                if any(m -> symbol(cs.species[j]) == m || startswith(symbol(cs.species[j]), m * "#"), names)
        ); init = 0.0,
    )
end

"""
    sa_batch(; limestone, level, T, thaumasite = true) -> NamedTuple

A batch of Schmidt et al. (2008): the binder with `limestone` (fraction),
fully hydrated at w/b 0.35, mixed with seven times the mass of its hydrated
solids of the Na2SO4 solution of subsystem `level` ("A" or "B"), at `T` (K),
and equilibrated in one step. The solution's S, Ca and Na (mmol/L), the pH,
and the solids.

A paste of w/b 0.35 fully hydrated has used nearly all its water, which the
limewater of the curing replaced; its solids are therefore taken from the same
binder at w/b 0.5, from which the batch is also started.
"""
function sa_batch(; limestone, level, T, thaumasite = true)
    rs = sa_paste(; limestone, w_b = 0.5, T, thaumasite)
    sub = sa_table("Schmidt2008", "subsystems")
    c = ustrip(us"mol/L", sub.so4_concentration[findfirst(==(level), sub.subsystem)])
    V = sa_value("Schmidt2008", "liquid_solid_ratio") * _sa_solids(rs.state).solids
    cs = rs.state.system
    b0 = budget(sa_recipe(; limestone, w_b = sa_value("Schmidt2008", "w_b"), T), cs).b
    eq, cert = equilibrate_certified(rs.state; model = rs.model, b = b0 .+ _sa_solution_budget(cs, V, c))
    n = ustrip.(us"mol", eq.n)
    V_L = ustrip(uconvert(us"L", volume(eq).liquid))
    dissolved(e) = 1000 * sum(n[k] * Float64(get(atoms(cs.species[k]), e, 0)) for k in cs.idx_aqueous) / V_L
    return (;
        certified = cert.optimal, S = dissolved(:S), Ca = dissolved(:Ca), Na = dissolved(:Na),
        pH = pH(eq, rs.model), _sa_solids(eq)...,
    )
end

"""
    sa_events(rows) -> NamedTuple

What a titration of [`sa_titration`](@ref) passes through: the volume (mL per
100 g of binder) at which the AFm of sulfate is gone, at which gypsum first
forms (`NaN` if never), at which portlandite is gone, and the largest SO3/CaO of
the solids with the volume at which it is reached.
"""
function sa_events(rows)
    first_V(f) = (i = findfirst(f, rows); i === nothing ? NaN : rows[i].V)
    r = [x.SO3 / x.CaO for x in rows]
    k = argmax(r)
    return (;
        afm_gone = first_V(x -> x.afm_so4 < 1.0e-6), gypsum = first_V(x -> x.gypsum > 1.0e-6),
        portlandite_gone = first_V(x -> x.portlandite < 1.0e-6), ratio_max = r[k], V_ratio_max = rows[k].V,
        gypsum_max = maximum(x -> x.gypsum, rows),
    )
end
