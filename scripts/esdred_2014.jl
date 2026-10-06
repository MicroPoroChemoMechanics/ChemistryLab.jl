# =============================================================================
#  esdred_2014.jl — a Portland cement with 40 % silica fume, the ESDRED
#  shotcrete paste of Lothenbach et al. (2014)
#
#  60 % CEM I 42.5 N, 40 % silica fume, 4.8 g of an aluminum sulfate set
#  accelerator and 1.2 g of a polycarboxylate superplasticizer per 100 g of
#  binder, w/b 0.5, sealed at 20 °C for 3.5 years; the pore solutions analyzed
#  from one hour on (their Table 3), the clinker and the silica fume followed by
#  29Si NMR (their Table 2). Miron et al. (2022b) computed this paste to test the
#  alkali end-members of CASH+NK. This file computes it at the measured degrees
#  of reaction, with its C-S-H as CSHQ or as CASH+NK, everything else equal, and
#  with or without the formate of its set accelerator; the page
#  docs/src/tutorials/validation_silica_fume_paste.md compares the pore
#  solutions.
#
#  ASSUMED, each where it is made: alite and belite react at the degree the NMR
#  gives their sum (Table 2), aluminate and ferrite by the law of Parrott and
#  Killoh with the constants of this cement (Table 4 of Miron et al. 2022b, from
#  this paper); the minor oxides of the clinker are shared among its phases as
#  in Lothenbach and Winnefeld (2006) and released with them; the silica fume
#  reacts as a whole at the degree of Table 2, its oxides with it; the dolomite
#  is inert; the normative composition, which sums to 100.31 g as printed, is
#  scaled to 100 g; the accelerator is available from the start, its oxides dissolved
#  and its water (40 % of it, its solid content being 60 %) added to the mixing
#  water; its organic carbon, the formate, either left out, as Miron et al.
#  (2022b) left it, or kept in solution whole, where the authors estimate that
#  the solids take up to 60 % of it (Table 3): two bounds; the superplasticizer
#  is left out.
# =============================================================================

using ChemistryLab
using DynamicQuantities
using OptimaSolver
using OrderedCollections

isdefined(@__MODULE__, :l08_system) || include(joinpath(@__DIR__, "lothenbach_2008.jl"))

const ES14 = "Lothenbach2014"
es14_value(name) = ustrip(literature_value(ES14, name))
es14_table(name) = literature_table(ES14, name)
const ES14_DB = Dict(symbol(s) => s for s in build_species(datapath("cemdata18-thermofun.json"); verbose = false))

# g/100 g of `material` in Table 1, zero where the table gives none or a bound.
function es14_composition(material)
    t = es14_table("composition")
    return OrderedDict(
        String(t.component[i]) => (t.qualifier[i] == "measured" ? Float64(ustrip(t.value[i])) : 0.0)
            for i in eachindex(t.material) if t.material[i] == material
    )
end
# The normative composition of Table 1 sums to 100.31 g as printed; it is scaled
# to 100 g.
function es14_normative(phase)
    t = es14_table("normative_phases")
    return Float64(ustrip(only(literature_table(ES14, "normative_phases"; phase).percent))) * 100 / sum(Float64.(ustrip.(t.percent)))
end

"""
    es14_nmr_extents() -> (; silicates, silica_fume)

The degrees of reaction of Table 2: of alite and belite together, one minus
the clinker's share of the silicon over its share at 0 days; of the silica
fume, one minus its unreacted mass over its mass at 0 days.
"""
function es14_nmr_extents()
    t = es14_table("nmr")
    d = ustrip.(us"d", t.time_d)
    clinker, sf = Float64.(ustrip.(t.clinker)), Float64.(ustrip.(t.silica_fume_unreacted))
    return (; silicates = TabulatedExtent(d, 1 .- clinker ./ clinker[1]), silica_fume = TabulatedExtent(d, 1 .- sf ./ sf[1]))
end

# The minor oxides of the clinker (Table 1, 'present as solid solution in the
# clinker phases') shared among its four phases as Lothenbach and Winnefeld
# (2006) share theirs, per gram of each phase.
function es14_minor_distribution()
    lw = literature_table("LothenbachWinnefeld2006", "minor_elements_in_clinker")
    lwp = literature_table("LothenbachWinnefeld2006", "normative_phases")
    lwpct(p) = ustrip(lwp.percent[findfirst(==(p), lwp.phase)])
    minors = es14_table("alkalis_in_clinker")
    grams(ox) = Float64(ustrip(only(minors.percent[i] for i in eachindex(minors.oxide) if minors.oxide[i] == ox)))
    out = OrderedDict(p => OrderedDict{String, Float64}() for p in keys(L08_CLINKER))
    for (ox, col) in (("Na2O", :Na), ("K2O", :K), ("MgO", :Mg))
        w = Dict(p => ustrip(getproperty(lw, col)[findfirst(==(p), lw.phase)]) / lwpct(p) * es14_normative(p) for p in keys(L08_CLINKER))
        for p in keys(L08_CLINKER)
            out[p][ox] = w[p] / sum(values(w)) * grams(ox)
        end
    end
    return out
end

"""
    es14_cement() -> Material

The CEM I 42.5 N of Table 1: alite and belite at the degree of the NMR,
aluminate and ferrite by Parrott and Killoh, each with its share of the minor
oxides; periclase, free lime, calcite, the calcium sulfates, syngenite and the
alkali sulfates available from the start; the dolomite inert.
"""
function es14_cement()
    T = es14_value("temperature") * u"K"
    w_c = es14_value("water_binder_ratio") / (es14_value("cem_i_percent") / 100)
    blaine = es14_value("surface_cem_i") * u"m^2/g"
    pk = literature_table("Miron2022b", "parrott_killoh")
    nmr = es14_nmr_extents()
    minors = es14_minor_distribution()
    constituents = AbstractConstituent[]
    for (phase, record) in L08_CLINKER
        ext = if phase in ("alite", "belite")
            nmr.silicates
        else
            i = findfirst(==(phase), pk.phase)
            parameters = (k₁ = ustrip(pk.K1[i]), n₁ = ustrip(pk.N1[i]), k₂ = ustrip(pk.K2[i]), k₃ = ustrip(pk.K3[i]), n₃ = ustrip(pk.N3[i]))
            ParrottKillohExtent(record; T, w_c, blaine = uconvert(us"m^2/kg", blaine), parameters, H = ustrip(pk.H[i]))
        end
        grams = minors[phase]
        total = sum(values(grams))
        push!(constituents, MineralConstituent(ES14_DB[record]; mass_fraction = (es14_normative(phase) - total) / 100, extent = ext, name = phase))
        push!(constituents, OxideConstituent("$phase minors", OrderedDict(k => v / total for (k, v) in grams); mass_fraction = total / 100, extent = ext))
    end
    push!(constituents, OxideConstituent("periclase", Dict("MgO" => 1.0); mass_fraction = es14_normative("periclase") / 100))
    push!(constituents, OxideConstituent("free lime", Dict("CaO" => 1.0); mass_fraction = es14_normative("CaO_free") / 100))
    for (phase, record) in (
            ("calcite", "Cal"), ("gypsum", "Gp"), ("anhydrite", "Anh"), ("hemihydrate", "hemihydrate"),
            ("syngenite", "syngenite"), ("K2SO4", "K2SO4"), ("Na2SO4", "Na2SO4"),
        )
        push!(constituents, MineralConstituent(ES14_DB[record]; mass_fraction = es14_normative(phase) / 100, name = phase))
    end
    push!(constituents, MineralConstituent(ES14_DB["Dis-Dol"]; mass_fraction = es14_normative("dolomite") / 100, extent = ConstantExtent(0.0), name = "dolomite"))
    return Material("CEM I 42.5 N (Lothenbach et al. 2014)", :cement; constituents)
end

"""
    es14_silica_fume() -> Material

The silica fume of Table 1, its oxides reacting together at the degree of the
NMR; its loss on ignition left out.
"""
function es14_silica_fume()
    c = es14_composition("silica fume")
    oxides = OrderedDict(k => v for (k, v) in c if k in ("SiO2", "CaO", "Al2O3", "Fe2O3", "MgO", "K2O", "SO3") && v > 0)
    total = sum(values(oxides))
    glass = OxideConstituent("silica fume", OrderedDict(k => v / total for (k, v) in oxides); mass_fraction = 1.0, extent = es14_nmr_extents().silica_fume)
    return Material("silica fume (Lothenbach et al. 2014)", :scm; constituents = AbstractConstituent[glass])
end

"""
    es14_accelerator(; formate = false) -> (material, grams)

The set accelerator as an addition per 100 g of binder: its oxides of Table 1,
dissolved from the start, and the water that its solid content leaves; with
`formate`, its organic carbon too, as the formate it is (Table 3, footnote a:
202 mM at the start, which the carbon of Table 1 gives), entered as the oxide of
carbon of valence 2, CO, which the water turns into formic acid and formate.
"""
function es14_accelerator(; formate = false)
    c = es14_composition("set accelerator")
    oxides = OrderedDict(k => v for (k, v) in c if k in ("Al2O3", "SO3", "MgO", "K2O", "Na2O", "CaO", "SiO2", "Fe2O3") && v > 0)
    formate && (oxides["CO"] = c["DOC"] * ustrip(us"g/mol", calculate_molar_mass(composition(Formula("CO")))) / ustrip(us"g/mol", calculate_molar_mass(composition(Formula("C")))))
    water = 100 - c["solid_content"]
    total = sum(values(oxides)) + water
    parts = merge(OrderedDict(k => v / total for (k, v) in oxides), OrderedDict("H2O" => water / total))
    m = Material("set accelerator (Lothenbach et al. 2014)", :other; constituents = AbstractConstituent[OxideConstituent("accelerator", parts; mass_fraction = 1.0)])
    return m, es14_value("set_accelerator_dose") * total / 100
end

"""
    es14_recipe(; formate = false) -> Recipe

The paste: 60 % CEM I and 40 % silica fume at w/b 0.5 and 20 °C, with the
accelerator added, its formate with `formate`.
"""
function es14_recipe(; formate = false)
    acc, grams = es14_accelerator(; formate)
    return Recipe(
        es14_cement() => es14_value("cem_i_percent") / 100, es14_silica_fume() => es14_value("silica_fume_percent") / 100;
        w_b = es14_value("water_binder_ratio"), T = es14_value("temperature") * u"K", additions = [acc => grams * u"g"],
    )
end

# The aqueous species of sulfur below sulfate and of iron below iron(III): left
# out when the formate is there, so that nothing can oxidize it, as nothing
# does in the paste, where it stays as formate (Table 3).
const ES14_REDUCED = ["S-2", "HS-", "H2S@", "S2O3-2", "SO3-2", "HSO3-", "Fe+2", "FeOH+", "Fe(SO4)@", "Fe(HSO4)+", "Fe(CO3)@", "Fe(HCO3)+"]

"""
    es14_system(gel; formate = false) -> ChemicalSystem

The system of the Portland paste of Lothenbach and Winnefeld (2006)
(`l08_system`), its C-S-H as `:CSHQ` or `:CASHNK`; with `formate`, the formate
and formic acid of the SUPCRT organic database added and the reduced species of
sulfur and iron left out (`ES14_REDUCED`).

ASSUMED: the formate of `slop98-organic-thermofun.json` shares the reference
state of the Cemdata18 aqueous species, both built on the HKF model.
"""
function es14_system(gel; formate = false)
    formate || return l08_system(gel)
    excluded = vcat(ES14_REDUCED, gel === :CASHNK ? ["NaOH@", "KOH@"] : String[])
    db = gel === :CSHQ ? build_species(datapath("cemdata18-thermofun.json"); verbose = false) :
        build_species(datapath("cemdata18-cashplus.json"); verbose = false)
    base = phase_list_system(L08_PHASES, db; exclude_aqueous = excluded, replace = gel === :CSHQ ? Dict{String, String}() : Dict("CSHQ" => "CASH+NK"))
    org = Dict(symbol(s) => s for s in build_species(datapath("slop98-organic-thermofun.json"); verbose = false))
    return ChemicalSystem(vcat(base.species, [org["For-"], org["ForH@"]]), CEMDATA_PRIMARIES; solid_solutions = base.solid_solutions)
end

"""The ages of Table 3 at and after one day, in days."""
es14_days() = sort(unique(d for d in ustrip.(us"d", es14_table("pore_solution").time_d) if d >= 1))

"""
    es14_measured(days, element) -> Float64

The concentration of `element` in the pore solution at `days` (mmol/L; the pH
for `"pH"`), the mean of the replicates of that age; `NaN` below the
detection limit.
"""
function es14_measured(days, element)
    if element == "pH"
        t = es14_table("pH")
        v = [Float64(t.pH[i]) for i in eachindex(t.pH) if ustrip(us"d", t.time_d[i]) == days]
        return sum(v) / length(v)
    end
    t = es14_table("pore_solution")
    v = [
        t.qualifier[i] in ("measured", "estimated", "from_pH") ? ustrip(us"mmol/L", t.concentration[i]) : NaN
            for i in eachindex(t.element) if ustrip(us"d", t.time_d[i]) == days && t.element[i] == element
    ]
    return isempty(v) ? NaN : sum(v) / length(v)
end
