# =============================================================================
#  lothenbach_winnefeld_2006.jl — the CEM I 42.5 N of Lothenbach & Winnefeld (2006)
#
#  A Portland cement hydrated at w/c 0.5 and 20 °C, whose composition (their
#  Tables 1 and 2), clinker kinetics (Table 4, with the w/c factor of Section
#  4.1) and pore solution over a year (Table 3) are transcribed in
#  data/literature/LothenbachWinnefeld2006.json. The validation page, its test and
#  the generator of the GEMS3K replay (test/reference/xgems_lw2006.py) include
#  this file, so that the three compute the same paste.
#
#  Written once here so that they cannot drift apart. The assumptions are the
#  file's, and each is stated where it is made.
# =============================================================================

using ChemistryLab
using DynamicQuantities
using OptimaSolver
using OrderedCollections

const LW06 = "LothenbachWinnefeld2006"
const LW06_SUBSTANCES = build_species(datapath("cemdata18-thermofun.json"); verbose = false)
const LW06_DB = Dict(symbol(s) => s for s in LW06_SUBSTANCES)

lw06_value(name) = ustrip(literature_value(LW06, name))
lw06_table(name) = literature_table(LW06, name)

# The phases the paste may form. The clinker phases and the calcium sulfates are
# the reactants; the hydrates are those of the paper's assemblage, with the C-S-H
# as the CSHQ model of Cemdata18 and its alkali end-members. The siliceous
# hydrogarnet of later Cemdata versions is left out: the paper predates it, and
# the GEMS3K export the replay runs on does not carry it. The iron then goes to
# microcrystalline FeOOH.
const LW06_PURE = split(
    "C3S C2S C3A C4AF Gp Anh hemihydrate Cal Portlandite K2SO4 Na2SO4 syngenite " *
        "ettringite monosulphate12 monocarbonate hemicarbonate C4AH13 C3AH6 C3FH6 " *
        "hydrotalcite Brc FeOOHmic AlOHmic straetlingite Amor-Sl"
)
const LW06_CSH = ["CSHQ-JenD", "CSHQ-JenH", "CSHQ-TobD", "CSHQ-TobH", "KSiOH", "NaSiOH"]

"""The chemical system of the paste."""
function lw06_system()
    sp = speciation(
        LW06_SUBSTANCES, vcat(LW06_PURE, LW06_CSH);
        aggregate_state = [AS_AQUEOUS], exclude_species = split("H2@ O2@ CH4@"),
    )
    return ChemicalSystem(
        sp, CEMDATA_PRIMARIES;
        solid_solutions = [SolidSolutionPhase("CSHQ", [LW06_DB[m] for m in LW06_CSH])],
    )
end

# The clinker phases as the paper names them, with their database records.
const LW06_CLINKER = OrderedDict("alite" => "C3S", "belite" => "C2S", "aluminate" => "C3A", "ferrite" => "C4AF")

"""
    lw06_material() -> Material

The cement: the four clinker phases, each dissolving by the rate law of Parrott
and Killoh with the w/c factor and the fineness of the paper, and carrying its
own sodium, potassium, magnesium and sulfate (their Table 2), released with it;
then free lime, calcite, the three calcium sulfates and the readily soluble
alkali sulfates, available from the start.
"""
function lw06_material()
    T = lw06_value("curing_temperature") * u"K"
    w_c = lw06_value("water_cement_ratio")
    blaine = lw06_value("blaine") * u"m^2/kg"
    phases = lw06_table("normative_phases")
    fraction(p) = ustrip(phases.percent[findfirst(==(p), phases.phase)]) / 100
    minors = lw06_table("minor_elements_in_clinker")
    M(ox) = ustrip(us"g/mol", Species(ox)[:M])
    constituents = AbstractConstituent[]
    for (phase, record) in LW06_CLINKER
        ext = ParrottKillohExtent(record; T, w_c, blaine)
        push!(constituents, MineralConstituent(LW06_DB[record]; mass_fraction = fraction(phase), extent = ext, name = phase))
        # Table 2 gives mmol per 100 g of cement; as oxides, in g per 100 g.
        i = findfirst(==(phase), minors.phase)
        grams = OrderedDict(
            "Na2O" => ustrip(minors.Na[i]) / 2 * M("Na2O") / 1000, "K2O" => ustrip(minors.K[i]) / 2 * M("K2O") / 1000,
            "MgO" => ustrip(minors.Mg[i]) * M("MgO") / 1000, "SO3" => ustrip(minors.SO3[i]) * M("SO3") / 1000,
        )
        total = sum(values(grams))
        push!(
            constituents,
            OxideConstituent(
                "$phase minors", OrderedDict(k => v / total for (k, v) in grams);
                mass_fraction = total / 100, extent = ext,
            ),
        )
    end
    push!(constituents, OxideConstituent("free lime", Dict("CaO" => 1.0); mass_fraction = fraction("CaO")))
    push!(constituents, MineralConstituent(LW06_DB["Cal"]; mass_fraction = fraction("CaCO3")))
    push!(constituents, MineralConstituent(LW06_DB["Anh"]; mass_fraction = lw06_value("anhydrite_percent") / 100))
    push!(constituents, MineralConstituent(LW06_DB["hemihydrate"]; mass_fraction = lw06_value("hemihydrate_percent") / 100))
    push!(constituents, MineralConstituent(LW06_DB["Gp"]; mass_fraction = lw06_value("gypsum_percent") / 100))
    push!(constituents, MineralConstituent(LW06_DB["K2SO4"]; mass_fraction = fraction("K2SO4")))
    push!(constituents, MineralConstituent(LW06_DB["Na2SO4"]; mass_fraction = fraction("Na2SO4")))
    push!(constituents, OxideConstituent("SrO", Dict("SrO" => 1.0); mass_fraction = fraction("SrO")))
    return Material("CEM I 42.5 N (Lothenbach and Winnefeld 2006)", :cement; constituents)
end

"""The paste: the cement at w/c 0.5 and 20 °C."""
lw06_recipe() = Recipe(
    lw06_material() => 1.0; w_b = lw06_value("water_cement_ratio"),
    T = lw06_value("curing_temperature") * u"K",
)

"""The ages of Table 3, in hours."""
lw06_hours() = sort(unique(ustrip.(lw06_table("pore_solution").time_h)))

"""
    lw06_measured(hours, element) -> (; value, below_limit)

The concentration of `element` Table 3 gives at `hours`, in mmol/L, and whether
it is a detection limit rather than a measurement (the table prints those as
`<x`). `value` is `NaN` when the table gives nothing.
"""
function lw06_measured(hours, element)
    ps = lw06_table("pore_solution")
    i = findfirst(k -> ustrip(ps.time_h[k]) == hours && ps.element[k] == element, eachindex(ps.element))
    i === nothing && return (; value = NaN, below_limit = false)
    return (; value = ustrip(ps.concentration[i]), below_limit = ps.qualifier[i] == "below_detection_limit")
end

"""
    lw06_pore_solution(state) -> Dict

The dissolved elements of an equilibrium, and the hydroxide, in mmol per kg of
water, the elements summed over every aqueous species that holds them.
"""
function lw06_pore_solution(state)
    cs = state.system
    n = ustrip.(us"mol", state.n)
    iw = only(cs.idx_solvent)
    kg = n[iw] * ustrip(us"kg/mol", cs.species[iw][:M])
    out = Dict{String, Float64}()
    for i in cs.idx_solutes, (el, k) in atoms(cs.species[i])
        (el === :H || el === :O) && continue
        out[String(el)] = get(out, String(el), 0.0) + 1000 * k * n[i] / kg
    end
    out["OH-"] = 1000 * n[findfirst(s -> symbol(s) == "OH-", cs.species)] / kg
    return out
end

"""
    lw06_elements(cs, b) -> OrderedDict

A budget of `cs`, written in its primaries, as amounts of the elements (and of
the charge, `Zz`), in mol: what a second code is given to replay it.
"""
function lw06_elements(cs, b)
    content(p, el) = el === :Zz ? Float64(charge(p)) : Float64(get(atoms(p), el, 0))
    return OrderedDict(
        String(el) => sum(content(p, el) * x for (p, x) in zip(cs.SM.primaries, b)) for el in cs.CSM.primaries
    )
end
