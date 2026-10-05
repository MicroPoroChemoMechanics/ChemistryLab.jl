# =============================================================================
#  lothenbach2008_temperature.jl — the hydrates of two cements from 0 to 60 °C
#
#  Lothenbach, Matschei, Möschner and Glasser (2008) calculated the hydrates of
#  a sulfate-resisting CEM I 52.5 N HTS (SRPC) and of a CEM II/A-L 42.5 R (PLC)
#  between 0 and 60 °C, at the degree of hydration of 150 days kept the same at
#  every temperature (their Figs. 5 and 6), with the thermodynamic data of
#  cemdata2007, and found ettringite and monocarbonate giving way to
#  monosulfate above about 48 °C. This file computes the same pastes with
#  Cemdata18: the clinker phases at the degrees the law of Parrott and Killoh
#  gives them after 150 days, with the constants of the article (its Table 2),
#  the rest at equilibrium at each temperature. The page
#  docs/src/examples/hydrates_temperature.md compares the two calculations.
# =============================================================================

using ChemistryLab
using DynamicQuantities
using OptimaSolver
using OrderedCollections

const L08T = "Lothenbach2008"
l08t_value(name) = ustrip(literature_value(L08T, name))
l08t_table(name) = literature_table(L08T, name)

# The clinker phases as the article names them, with their database records.
const L08T_CLINKER = OrderedDict("alite" => "C3S", "belite" => "C2S", "aluminate" => "C3A", "ferrite" => "C4AF")
# The phases of the paste: those of the Portland paste of Lothenbach and
# Winnefeld (2006), whose list holds monocarbonate, hemicarbonate and the AFm
# solid solution where monosulfate is.
const L08T_PHASES = "Portland paste (Lothenbach and Winnefeld 2006)"
const L08T_DB = Dict(symbol(s) => s for s in build_species(datapath("cemdata18-thermofun.json"); verbose = false))

"""The chemical system of the pastes, on Cemdata18."""
l08t_system() = phase_list_system(L08T_PHASES, collect(values(L08T_DB)))

"""The water/cement ratio of `cement`: the paste of SRPC, the mortar of PLC."""
l08t_water_cement(cement) = cement == "SRPC" ? l08t_value("water_cement_ratio_srpc") : l08t_value("water_cement_ratio_mortar")

l08t_fraction(cement, phase) = (t = l08t_table("normative_phases"); ustrip(getproperty(t, Symbol(cement))[findfirst(==(phase), t.phase)]) / 100)

"""
    l08t_degrees(cement) -> OrderedDict

The degree of hydration of each clinker phase of `cement` after the hydration
time of Figs. 5 and 6 (150 days), under the law of Parrott and Killoh with the
constants of Table 2, the w/c factor of Section 3.1, the fineness of Table 1 and
the temperature the constants are given for, 20 °C.
"""
function l08t_degrees(cement)
    pk = l08t_table("parrott_killoh")
    w_c = l08t_water_cement(cement)
    blaine = literature_value(L08T, "blaine_" * lowercase(cement))
    days = ustrip(us"d", literature_value(L08T, "hydration_time_figs"))
    return OrderedDict(
        phase => begin
            i = findfirst(==(phase), pk.phase)
            parameters = (k₁ = ustrip(pk.K1[i]), n₁ = ustrip(pk.N1[i]), k₂ = ustrip(pk.K2[i]), k₃ = ustrip(pk.K3[i]), n₃ = ustrip(pk.N3[i]))
            extent(ParrottKillohExtent(record; T = PK84_PARAMS_C3S.T_ref, w_c, blaine, parameters), days)
        end for (phase, record) in L08T_CLINKER
    )
end

"""The overall degree of hydration of the clinker of `cement`, its phases weighted by their masses."""
function l08t_overall_degree(cement)
    α = l08t_degrees(cement)
    m = Dict(p => l08t_fraction(cement, p) for p in keys(L08T_CLINKER))
    return sum(m[p] * α[p] for p in keys(α)) / sum(values(m))
end

"""
    l08t_material(cement) -> Material

The cement of Table 1: the four clinker phases at their degrees after 150 days
(`l08t_degrees`), and the readily soluble part whole: free lime, calcite, the
calcium sulfate (as anhydrite, the formula the table prints) and the alkali
sulfates.

ASSUMED: the minor oxides in solid solution in the clinker phases (K2O, Na2O,
MgO, SO3) are released with the clinker as a whole, at its overall degree; the
article gives their totals and not how they are shared among the phases.
"""
function l08t_material(cement)
    α = l08t_degrees(cement)
    f(p) = l08t_fraction(cement, p)
    constituents = AbstractConstituent[
        MineralConstituent(L08T_DB[record]; mass_fraction = f(phase), extent = α[phase], name = phase)
            for (phase, record) in L08T_CLINKER
    ]
    minors = ("K2O", "Na2O", "MgO", "SO3")
    total = sum(f, minors)
    push!(
        constituents, OxideConstituent(
            "minor oxides", OrderedDict(ox => f(ox) / total for ox in minors);
            mass_fraction = total, extent = l08t_overall_degree(cement),
        ),
    )
    push!(constituents, OxideConstituent("free lime", Dict("CaO" => 1.0); mass_fraction = f("CaO (free)")))
    push!(constituents, MineralConstituent(L08T_DB["Cal"]; mass_fraction = f("CaCO3")))
    push!(constituents, MineralConstituent(L08T_DB["Anh"]; mass_fraction = f("CaSO4")))
    push!(constituents, MineralConstituent(L08T_DB["K2SO4"]; mass_fraction = f("K2SO4")))
    push!(constituents, MineralConstituent(L08T_DB["Na2SO4"]; mass_fraction = f("Na2SO4")))
    return Material("$cement (Lothenbach et al. 2008)", :cement; constituents)
end

"""The paste of `cement` at its w/c and at `T_C` °C."""
l08t_recipe(cement, T_C; material = l08t_material(cement)) =
    Recipe(material => 1.0; w_b = l08t_water_cement(cement), T = (T_C + 273.15) * u"K")

"""
    l08t_state(cement, T_C; cs = l08t_system(), model = cemdata18_activity_model(:KOH)) -> RecipeState

The certified equilibrium of the paste of `cement` at `T_C` °C.
"""
function l08t_state(cement, T_C; cs = l08t_system(), model = cemdata18_activity_model(:KOH), material = l08t_material(cement))
    rs, cert = equilibrate_certified(l08t_recipe(cement, T_C; material), cs; model)
    cert.optimal || error("l08t_state: the equilibrium of $cement at $T_C °C is not certified")
    return rs
end

"""
    l08t_volumes(rs) -> OrderedDict

The volume of each solid of the equilibrium `rs`, and of the unreacted clinker
phases, in cm³ per 100 g of cement, the phases above 0.01 cm³. The unreacted
minor oxides have no sourced density and are left out (0.2 cm³ at most).
"""
function l08t_volumes(rs)
    st = rs.state
    V = ustrip(uconvert(us"cm^3", volume(st).total))
    solids = Set(symbol(s) for s in st.system.species if aggregate_state(s) != AS_AQUEOUS)
    out = OrderedDict{String, Float64}()
    for (k, fr) in volume_fractions(st)
        k in solids && (out[k] = fr * V)
    end
    out["unreacted clinker"] = sum(x.volume for x in rs.residual if x.constituent in keys(L08T_CLINKER); init = 0.0)
    return OrderedDict(k => v for (k, v) in sort(collect(out); by = last, rev = true) if v > 0.01)
end

"""
    l08t_has_monosulfate(rs) -> Bool

Whether the equilibrium `rs` holds monosulfate (the `monosulphate12` member of
the AFm solid solution, above a micromole).
"""
function l08t_has_monosulfate(rs)
    st = rs.state
    i = findfirst(s -> symbol(s) == "monosulphate12", st.system.species)
    return ustrip(us"mol", st.n[i]) > 1.0e-6
end

"""
    l08t_transition(cement; lo = 20.0, hi = 60.0, tol = 0.05, cs = l08t_system(), model = ...) -> Float64

The temperature, °C, above which the equilibrium of `cement` holds monosulfate,
found by bisection between `lo` (without) and `hi` (with) to `tol`.
"""
function l08t_transition(
        cement; lo = 20.0, hi = 60.0, tol = 0.05, cs = l08t_system(),
        model = cemdata18_activity_model(:KOH), material = l08t_material(cement),
    )
    has(T) = l08t_has_monosulfate(l08t_state(cement, T; cs, model, material))
    has(lo) && error("l08t_transition: $cement already holds monosulfate at $lo °C")
    has(hi) || error("l08t_transition: $cement holds no monosulfate at $hi °C")
    while hi - lo > tol
        mid = (lo + hi) / 2
        has(mid) ? (hi = mid) : (lo = mid)
    end
    return (lo + hi) / 2
end

# The species of the equilibrium gathered as the figures of the article label
# their layers: the C-S-H with its end members, the AFm solid solution as
# monosulfate. `monosulphate12` is the Cemdata18 identifier.
const L08T_GROUPS = OrderedDict(
    "ettringite" => ("ettringite",), "monocarbonate" => ("monocarbonate",), "hemicarbonate" => ("hemicarbonate",),
    "monosulfate" => ("monosulphate12", "C4AH13"), "calcite" => ("Cal",), "portlandite" => ("Portlandite",),
    "hydrotalcite" => ("hydrotalcite",),
    "C-S-H" => ("CSHQ-JenD", "CSHQ-JenH", "CSHQ-TobD", "CSHQ-TobH", "KSiOH", "NaSiOH"),
    "iron hydroxide" => ("FeOOHmic",), "unhydrated clinker" => ("unreacted clinker",),
)

"""
    l08t_grouped(volumes) -> OrderedDict

The volumes of `l08t_volumes` gathered by `L08T_GROUPS`, anything else under
`"other"`.
"""
function l08t_grouped(v::AbstractDict)
    out = OrderedDict(g => sum(get(v, s, 0.0) for s in members) for (g, members) in L08T_GROUPS)
    rest = setdiff(keys(v), Iterators.flatten(values(L08T_GROUPS)))
    isempty(rest) || (out["other"] = sum(v[k] for k in rest))
    return out
end
