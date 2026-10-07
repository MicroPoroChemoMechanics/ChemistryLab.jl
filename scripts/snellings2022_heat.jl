# =============================================================================
#  snellings2022_heat.jl — the heat of the slag-limestone cement pastes
#
#  The ternary cement of Snellings et al. (2022), 50 % CEM I 52.5 R, 40 % slag
#  and 10 % limestone, at the degrees of reaction the paper measured for its
#  clinker and its slag, and the heat it releases from the mixing, against the
#  isothermal calorimetry of the paper (Fig. 3, read in
#  data/literature/Snellings2022.json). Nothing is fitted: the heat is the fall
#  of the enthalpy of the paste between its state at mixing and its state at
#  each age, the hydrates and the solution from Cemdata18 and the unreacted
#  glass from glass_enthalpy, built from measured silicate glasses.
#
#  ASSUMED, each where it is made: the four clinker phases react at the one
#  degree the paper gives for the clinker; the limestone and the sulfates are
#  at equilibrium from the mixing in both states, so their heat is not counted;
#  the quartz is inert. The page docs/src/examples/slag_heat.md uses this.
# =============================================================================

isdefined(@__MODULE__, :sn22p_setup) || include(joinpath(@__DIR__, "snellings2022_pastes.jl"))

"""
    sn22h_system() -> ChemicalSystem

The phases of the pastes of `snellings2022_pastes.jl`, without the pseudo-species
of the glass: here the glass enters through its oxides, at a prescribed degree
of reaction.
"""
function sn22h_system()
    members = reduce(vcat, last.(SN22_SS))
    sp = speciation(
        collect(values(SN22_DB)), vcat(SN22_PURE, collect(keys(SN22_PK)), members);
        aggregate_state = [AS_AQUEOUS], exclude_species = split("H2@ O2@ CH4@"),
    )
    ss = [SolidSolutionPhase(n, [SN22_DB[m] for m in ms]) for (n, ms) in SN22_SS]
    return ChemicalSystem(sp, CEMDATA_PRIMARIES; solid_solutions = ss)
end

"""
    sn22h_glass(T_C) -> NamedTuple

The enthalpy of the glass of the slag at `T_C` °C, from [`glass_enthalpy`](@ref),
on the scale of Cemdata18 for the oxides it holds (lime, quartz, hematite and the
alkali oxides), aq17 and slop98 for the others. Phosphorus has no crystal in
these databases and is left out by name.
"""
const SN22H_REFERENCE = ([SN22_DB[s] for s in ("Lim", "Qtz", "Hem", "K2O", "Na2O")], ChemistryLab._default_glass_reference()...)

function sn22h_glass(T_C)
    glass = only(c for c in material_template("slag (Snellings 2022)", SN22_DB).constituents if c.name == "glass")
    ignore = Tuple(ox for ox in keys(glass.oxides) if ox == "P2O5")
    return glass_enthalpy(glass.oxides; T = _sn22_kelvin(T_C), reference = SN22H_REFERENCE, ignore)
end

# The enthalpy per gram of the oxides of the glass the system has no primary for,
# set aside as the glass reacts: their crystals', the enthalpy glass_enthalpy
# counted them at; phosphorus, left out of the glass, at zero as it was there.
function sn22h_set_aside()
    h(ox) = ustrip(us"J/mol", first(ChemistryLab._reference_oxide(ox, SN22H_REFERENCE))[:ΔₐH⁰](T = 298.15u"K", P = 1.0e5u"Pa"; unit = true)) /
        ustrip(us"g/mol", Species(ox)[:M])
    return Dict("TiO2" => h("TiO2"), "MnO" => h("MnO"), "P2O5" => 0.0)
end

# The degree of reaction the paper measured (fraction), at an age of its table.
function sn22h_degree(constituent, T_C, w_b, age)
    t = literature_table(SN22, "degree_of_reaction"; constituent, temperature_C = T_C, w_b, age = age * u"d")
    return ustrip(only(t.degree_percent)) / 100
end

"""
    sn22h_state(cs, T_C, w_b, αc, αs; h_glass, model) -> RecipeState

The paste, 100 g of binder, with the four clinker phases at the degree `αc` and
the glass of the slag at `αs`, its enthalpy `h_glass` (J/g).
"""
function sn22h_state(cs, T_C, w_b, αc, αs; h_glass, model)
    inert_quartz(m) = with_extents(m, Dict("Quartz" => 0.0))
    pc = with_extents(
        inert_quartz(material_template("PC (Snellings 2022)", SN22_DB)),
        Dict(p => αc for p in ("Alite", "Belite", "C3A", "C4AF"))
    )
    slag = with_extents(inert_quartz(material_template("slag (Snellings 2022)", SN22_DB)), Dict("glass" => αs))
    slag = with_enthalpy(slag, Dict("glass" => h_glass); source = "glass_enthalpy")
    ls = inert_quartz(material_template("limestone (Snellings 2022)", SN22_DB))
    r = Recipe(
        pc => sn22_value("pc_percent") / 100, slag => sn22_value("slag_percent") / 100,
        ls => sn22_value("limestone_percent") / 100; w_b, T = _sn22_kelvin(T_C)
    )
    rs, _ = equilibrate_certified(r, cs; model)
    return rs
end

"""
    sn22h_heat(cs, T_C, w_b; ages = (1, 2, 7, 28), model, vitrification = true) -> Vector

The heat (J per gram of Portland cement, as the paper reports it) released from
the mixing to each age, with the glass at its enthalpy from
[`glass_enthalpy`](@ref), or, with `vitrification = false`, at the enthalpy of
the crystals of its composition, the reference that shows what the enthalpy of
vitrification contributes.
"""
function sn22h_heat(cs, T_C, w_b; ages = (1, 2, 7, 28), model, vitrification = true)
    g = sn22h_glass(T_C)
    h = vitrification ? g.enthalpy : g.enthalpy - g.vitrification
    pc_g = sn22_value("pc_percent")                   # grams of Portland cement in 100 g of binder
    rs0 = sn22h_state(cs, T_C, w_b, 0.0, 0.0; h_glass = h, model)
    aside = sn22h_set_aside()
    return [
        heat_release(
            rs0, sn22h_state(
                cs, T_C, w_b, sn22h_degree("clinker", T_C, w_b, a),
                sn22h_degree("slag", T_C, w_b, a); h_glass = h, model
            ); set_aside = aside
        ) / pc_g
            for a in ages
    ]
end

# The cumulative heat the paper measured (J per gram of Portland cement), read off
# its Fig. 3, at the ages of the table nearest those asked.
function sn22h_measured(T_C, w_b, age)
    t = literature_table(SN22, "cumulative_heat"; temperature_C = T_C, w_b)
    days = ustrip.(us"d", t.age)
    i = argmin(abs.(days .- age))
    return (age = days[i], heat = ustrip(us"J/g", t.heat[i]))
end
