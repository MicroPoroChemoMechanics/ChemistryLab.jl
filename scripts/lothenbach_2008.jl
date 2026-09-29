# =============================================================================
#  lothenbach_2008.jl — the Portland cement with 4 % limestone of Lothenbach,
#  Le Saout, Gallucci and Scrivener (2008), with two models of its C-S-H
#
#  A Portland cement, alone (PC) and blended with 4 % limestone (PC4), hydrated
#  at w/c 0.4 and 20 °C, whose normative phases (their Table 1), clinker kinetics
#  (Table 3, with the critical degrees of Section 3.2) and pore solution over 400
#  days (Table 4) are transcribed in data/literature/LothenbachLeSaout2008.json.
#  Miron et al. (2022b) computed the PC4 paste with the CASH+NK model of its
#  C-S-H; the CASH+ page and its test include this file to compute it with CSHQ
#  and with CASH+NK, on the same recipe.
#
#  The assumptions are the file's, and each is stated where it is made.
# =============================================================================

using ChemistryLab
using DynamicQuantities
using OptimaSolver
using OrderedCollections

include(joinpath(@__DIR__, "validation_common.jl"))

const L08 = "LothenbachLeSaout2008"
l08_value(name) = ustrip(literature_value(L08, name))
l08_table(name) = literature_table(L08, name)

# The phases of the paste: those of the Portland paste of Lothenbach and
# Winnefeld (2006), a paste of the same laboratory and the same modeling, whose
# list already holds monocarbonate and hemicarbonate for a cement with calcite.
const L08_PHASES = "Portland paste (Lothenbach and Winnefeld 2006)"

"""
    l08_system(gel = :CSHQ) -> ChemicalSystem

The chemical system of the paste with its C-S-H as `CSHQ` (on Cemdata18, the
database it was fitted with) or as `CASH+NK` (on `cemdata18-cashplus.json`, with
the CaSiO3@ and Ca(OH)2@ complexes it was fitted with). For CASH+NK the neutral
complexes NaOH@ and KOH@ are left out, as Miron et al. (2022a) left them out when
they fitted the alkali end-members.
"""
function l08_system(gel::Symbol = :CSHQ)
    if gel === :CSHQ
        return phase_list_system(L08_PHASES, build_species(datapath("cemdata18-thermofun.json"); verbose = false))
    elseif gel === :CASHNK
        return phase_list_system(
            L08_PHASES, build_species(datapath("cemdata18-cashplus.json"); verbose = false);
            replace = Dict("CSHQ" => "CASH+NK"), exclude_aqueous = ["NaOH@", "KOH@"],
        )
    end
    throw(ArgumentError("l08_system: the gel is :CSHQ or :CASHNK"))
end

# The clinker phases as the paper names them, with their database records.
const L08_CLINKER = OrderedDict("alite" => "C3S", "belite" => "C2S", "aluminate" => "C3A", "ferrite" => "C4AF")

# How the minor oxides are shared among the clinker phases. The paper gives
# their totals (Table 1) and not their distribution; this file assumes the
# distribution of Lothenbach and Winnefeld (2006, Table 2, after Taylor), as a
# content per gram of each phase, scaled to the totals of this cement. That is an
# assumption of this file, stated here and on the page.
function l08_minor_distribution(cement::AbstractString)
    lw = literature_table("LothenbachWinnefeld2006", "minor_elements_in_clinker")
    lwp = literature_table("LothenbachWinnefeld2006", "normative_phases")
    phases = l08_table("normative_phases")
    col = Symbol(cement)
    own(p) = ustrip(getproperty(phases, col)[findfirst(==(p), phases.phase)])
    lwpct(p) = ustrip(lwp.percent[findfirst(==(p), lwp.phase)])
    # mmol per 100 g of cement in 2006, per gram of phase there, times this
    # cement's grams of phase: a share for each phase and element.
    share = Dict{Tuple{String, Symbol}, Float64}()
    for (el, col06) in ((:Na2O, :Na), (:K2O, :K), (:MgO, :Mg), (:SO3, :SO3))
        w = Dict(p => ustrip(getproperty(lw, col06)[findfirst(==(p), lw.phase)]) / lwpct(p) * own(p) for p in keys(L08_CLINKER))
        total = sum(values(w))
        for p in keys(L08_CLINKER)
            share[(p, el)] = w[p] / total
        end
    end
    minors = l08_table("clinker_minor_oxides")
    grams(ox) = ustrip(getproperty(minors, col)[findfirst(==(String(ox)), minors.oxide)])
    return OrderedDict(
        p => OrderedDict(String(ox) => share[(p, ox)] * grams(ox) for ox in (:Na2O, :K2O, :MgO, :SO3))
            for p in keys(L08_CLINKER)
    )
end

"""
    l08_material(cement = "PC4") -> Material

The cement: the four clinker phases, each dissolving by the rate law of Parrott
and Killoh with the constants, critical degrees, w/c and fineness of the paper,
and carrying its share of the minor oxides (`l08_minor_distribution`), released
with it; then periclase, free lime, calcite, gypsum and the readily soluble
alkali sulfates, available from the start.
"""
function l08_material(cement::AbstractString = "PC4")
    T = l08_value("curing_temperature") * u"K"
    w_c = l08_value("water_cement_ratio")
    blaine = l08_value("blaine_$cement") * u"m^2/kg"
    phases = l08_table("normative_phases")
    col = Symbol(cement)
    fraction(p) = ustrip(getproperty(phases, col)[findfirst(==(p), phases.phase)]) / 100
    pk = l08_table("parrott_killoh")
    minors = l08_minor_distribution(cement)
    db = Dict(symbol(s) => s for s in build_species(datapath("cemdata18-thermofun.json"); verbose = false))
    constituents = AbstractConstituent[]
    for (phase, record) in L08_CLINKER
        i = findfirst(==(phase), pk.phase)
        parameters = (k₁ = ustrip(pk.K1[i]), n₁ = ustrip(pk.N1[i]), k₂ = ustrip(pk.K2[i]), k₃ = ustrip(pk.K3[i]), n₃ = ustrip(pk.N3[i]))
        ext = ParrottKillohExtent(record; T, w_c, blaine, parameters, H = ustrip(pk.H[i]))
        grams = minors[phase]
        total = sum(values(grams))
        # The phase's mass holds its minor oxides (the normative phases sum to
        # 100 g with them inside), so the phase proper is the rest.
        push!(constituents, MineralConstituent(db[record]; mass_fraction = fraction(phase) - total / 100, extent = ext, name = phase))
        push!(
            constituents,
            OxideConstituent("$phase minors", OrderedDict(k => v / total for (k, v) in grams); mass_fraction = total / 100, extent = ext),
        )
    end
    push!(constituents, OxideConstituent("periclase", Dict("MgO" => 1.0); mass_fraction = fraction("MgO")))
    push!(constituents, OxideConstituent("free lime", Dict("CaO" => 1.0); mass_fraction = fraction("CaO")))
    push!(constituents, MineralConstituent(db["Cal"]; mass_fraction = fraction("CaCO3")))
    push!(constituents, MineralConstituent(db["Gp"]; mass_fraction = fraction("gypsum")))
    push!(constituents, MineralConstituent(db["K2SO4"]; mass_fraction = fraction("K2SO4")))
    push!(constituents, MineralConstituent(db["Na2SO4"]; mass_fraction = fraction("Na2SO4")))
    return Material("$cement (Lothenbach et al. 2008)", :cement; constituents)
end

"""The paste: the cement at w/c 0.4 and 20 °C."""
l08_recipe(cement::AbstractString = "PC4") = Recipe(
    l08_material(cement) => 1.0; w_b = l08_value("water_cement_ratio"),
    T = l08_value("curing_temperature") * u"K",
)

"""The ages of Table 4 at and after one day, in days."""
l08_days() = sort(unique(d for d in ustrip.(us"d", l08_table("pore_solution").time) if d >= 1))

"""
    l08_measured(cement, days, element) -> Float64

The concentration of `element` (a column of Table 4: Na, K, Ca, Al, Si, S, OH, or
pH) in the pore solution of `cement` at `days`, in mmol/L.
"""
function l08_measured(cement, days, element)
    ps = l08_table("pore_solution")
    i = findfirst(k -> ps.cement[k] == cement && ustrip(us"d", ps.time[k]) == days, eachindex(ps.cement))
    v = getproperty(ps, Symbol(element))[i]
    return element == "pH" ? Float64(v) : ustrip(us"mmol/L", v)
end
