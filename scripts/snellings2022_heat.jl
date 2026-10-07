# =============================================================================
#  snellings2022_heat.jl — the heat of the slag-limestone cement pastes
#
#  The ternary cement of Snellings et al. (2022), 50 % CEM I 52.5 R, 40 % slag
#  and 10 % limestone, at the degrees of reaction the paper measured for its
#  four clinker phases and its slag, and the heat it releases from the mixing,
#  against the isothermal calorimetry of the paper (Fig. 3, read in
#  data/literature/Snellings2022.json). The heat is the fall of the enthalpy of
#  the paste between its state at mixing and its state at each age, the
#  hydrates and the solution from Cemdata18 and the unreacted glass from
#  glass_enthalpy, built from measured silicate glasses.
#
#  The same states give the bound water and the portlandite the paper measured
#  by thermogravimetry (Fig. 8, w/b 0.5). sn22h_progress fits one factor on the
#  degrees of reaction so that the bound water is the measured one; the heat at
#  that progress is then compared with the calorimetry, which nothing is fitted
#  on.
#
#  ASSUMED, each where it is made: the limestone and the sulfates are at
#  equilibrium from the mixing in both states, so their heat is not counted;
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

# The four clinker phases: their names in the template of the cement, and in
# Fig. 5 of the paper.
const SN22H_PHASES = ("Alite" => "C3S", "Belite" => "C2S", "C3A" => "C3A", "C4AF" => "C4AF")

"""
    sn22h_phase_degrees(T_C, w_b, age) -> Dict

The degree of reaction (fraction) of each clinker phase of the paste at `w_b`
cured at `T_C` °C, at an age of Fig. 5 (days): one minus its content over its
content at the mixing, the horizontal line of the figure, as the authors form
the degree of the clinker from the same contents. A content read a pixel below
zero is a phase that has reacted completely: the degree is kept within [0, 1].
"""
function sn22h_phase_degrees(T_C, w_b, age)
    initial = literature_table(SN22, "clinker_phase_initial")
    return Dict(
        name => begin
                c = literature_table(SN22, "clinker_phase_content"; phase, temperature_C = T_C, w_b, age = age * u"d")
                c0 = initial.percent[findfirst(==(phase), initial.phase)]
                clamp(1 - ustrip(only(c.percent)) / ustrip(c0), 0.0, 1.0)
            end
            for (name, phase) in SN22H_PHASES
    )
end

"""
    sn22h_state(cs, T_C, w_b, αc, αs; h_glass, model) -> RecipeState

The paste, 100 g of binder, with the clinker phases at the degrees `αc`, a
`Dict` by phase or one number for all four, and the glass of the slag at `αs`,
its enthalpy `h_glass` (J/g).
"""
function sn22h_state(cs, T_C, w_b, αc, αs; h_glass, model)
    inert_quartz(m) = with_extents(m, Dict("Quartz" => 0.0))
    extents = αc isa AbstractDict ? Dict{String, Any}(αc) : Dict{String, Any}(p => αc for p in first.(SN22H_PHASES))
    pc = with_extents(inert_quartz(material_template("PC (Snellings 2022)", SN22_DB)), extents)
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
    sn22h_origin(cs, T_C, w_b; model, vitrification = true) -> (; rs0, state, heat)

The paste at the mixing, `rs0`, nothing reacted; `state(αc, αs)`, the same
paste with its clinker phases at `αc` and its slag at `αs`; and `heat(rs)`,
the heat released from the mixing to the state `rs`, in J per gram of Portland
cement, as the paper reports its calorimetry. The glass is at its enthalpy from
[`glass_enthalpy`](@ref), or, with `vitrification = false`, at the enthalpy of
the crystals of its composition, the reference that shows what the enthalpy of
vitrification contributes.
"""
function sn22h_origin(cs, T_C, w_b; model, vitrification = true)
    g = sn22h_glass(T_C)
    h = vitrification ? g.enthalpy : g.enthalpy - g.vitrification
    rs0 = sn22h_state(cs, T_C, w_b, 0.0, 0.0; h_glass = h, model)
    aside = sn22h_set_aside()
    pc_g = sn22_value("pc_percent")                   # grams of Portland cement in 100 g of binder
    state(αc, αs) = sn22h_state(cs, T_C, w_b, αc, αs; h_glass = h, model)
    heat(rs) = heat_release(rs0, rs; set_aside = aside) / pc_g
    return (; rs0, state, heat)
end

"""
    sn22h_heat(cs, T_C, w_b; ages = (1, 2, 7, 28), model, vitrification = true) -> Vector

The heat (J per gram of Portland cement) released from the mixing to each age,
at the degrees of reaction of the diffraction: each clinker phase at its own
(Fig. 5), the slag at its (Fig. 6b).
"""
function sn22h_heat(cs, T_C, w_b; ages = (1, 2, 7, 28), model, vitrification = true)
    o = sn22h_origin(cs, T_C, w_b; model, vitrification)
    return [o.heat(o.state(sn22h_phase_degrees(T_C, w_b, a), sn22h_degree("slag", T_C, w_b, a))) for a in ages]
end

"""
    sn22h_tga(rs) -> (; bound_water, portlandite)

The bound water and the portlandite of the state `rs` as the thermogravimetry of
Fig. 8 reports them, in g per 100 g of the mass at 550 °C.

ASSUMED: every hydrate has lost all its water by 550 °C and nothing else has
been lost, the carbonates keeping their CO₂; the mass at 550 °C is every solid,
the unreacted part of the binder and the oxides set aside included, less that
water.
"""
function sn22h_tga(rs)
    g(q) = ustrip(uconvert(us"g", q))
    water = g(bound_water(rs.state))
    dry = g(mass(rs.state).solid) + sum(x.mass for x in rs.residual; init = 0.0) - water
    ch = g(mass(rs.state, rs.state.system.dict_species["Portlandite"]))
    return (bound_water = 100 * water / dry, portlandite = 100 * ch / dry)
end

# What the thermogravimetry of Fig. 8 measured on the paste at w/b 0.5 (g per
# 100 g of the mass at 550 °C), at an age of the figure (days).
sn22h_measured_tga(phase, T_C, age) = ustrip(
    only(literature_table(SN22, "hydrates_wb05"; phase, technique = "TGA", temperature_C = T_C, age = age * u"d").percent)
)

"""
    sn22h_progress(cs, T_C, age; model, origin) -> (; λ, state, tga)

The paste at w/b 0.5 at `age` (days), with the degrees of reaction of the
diffraction, the four clinker phases and the slag alike, multiplied by the
factor `λ` for which its bound water is the one the thermogravimetry measured
(Fig. 8). `origin` is the [`sn22h_origin`](@ref) of the paste.

ASSUMED: the diffraction tells how the reaction is shared among the five
constituents, the thermogravimetry how far it has gone. The bound water grows
with `λ`, which is found by regula falsi, in its Illinois form, from the bracket
0.1 to 1, to 0.1 % of the measured value.
"""
function sn22h_progress(cs, T_C, age; model, origin = sn22h_origin(cs, T_C, 0.5; model))
    αc = sn22h_phase_degrees(T_C, 0.5, age)
    αs = sn22h_degree("slag", T_C, 0.5, age)
    target = sn22h_measured_tga("bound water", T_C, age)
    at(λ) = origin.state(Dict(p => λ * a for (p, a) in αc), λ * αs)
    excess(rs) = sn22h_tga(rs).bound_water - target
    lo, hi = 0.1, 1.0
    flo, fhi = excess(at(lo)), excess(at(hi))
    flo < 0 < fhi || error("sn22h_progress: the measured bound water is not reached between λ = $lo and $hi")
    λ, rs, side = hi, at(hi), 0
    for _ in 1:30
        λ = (lo * fhi - hi * flo) / (fhi - flo)
        rs = at(λ)
        f = excess(rs)
        abs(f) <= 1.0e-3 * target && break
        if f > 0
            hi, fhi = λ, f
            side == -1 && (flo /= 2)
            side = -1
        else
            lo, flo = λ, f
            side == 1 && (fhi /= 2)
            side = 1
        end
    end
    return (; λ, state = rs, tga = sn22h_tga(rs))
end

# The cumulative heat the paper measured (J per gram of Portland cement), read off
# its Fig. 3, at the ages of the table nearest those asked.
function sn22h_measured(T_C, w_b, age)
    t = literature_table(SN22, "cumulative_heat"; temperature_C = T_C, w_b)
    days = ustrip.(us"d", t.age)
    i = argmin(abs.(days .- age))
    return (age = days[i], heat = ustrip(us"J/g", t.heat[i]))
end
