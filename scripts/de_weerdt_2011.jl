# =============================================================================
#  de_weerdt_2011.jl — the ternary cements of De Weerdt et al. (2011)
#
#  A clinker interground with gypsum (OPC), blended with 5 % limestone powder,
#  35 % siliceous fly ash, or 30 % fly ash and 5 % limestone, hydrated at w/b 0.5
#  and 20 °C. Their materials (Tables 1 to 4), the degrees of reaction of the
#  clinker phases (Table 7, XRD-Rietveld), the reaction of the fly ash (the fit
#  printed on Fig. 7), the portlandite (Table 7) and the pore solution (Table 8)
#  are transcribed in data/literature/DeWeerdt2011.json. The validation page, its
#  test and the generator of the GEMS3K replay include this file.
#
#  The assumptions are the file's, and each is stated where it is made.
# =============================================================================

using ChemistryLab
using DynamicQuantities
using OptimaSolver
using OrderedCollections

include(joinpath(@__DIR__, "validation_common.jl"))

const DW11 = "DeWeerdt2011"
const DW11_SUBSTANCES = build_species(datapath("cemdata18-thermofun.json"); verbose = false)
const DW11_DB = Dict(symbol(s) => s for s in DW11_SUBSTANCES)

dw11_value(name) = ustrip(literature_value(DW11, name))
dw11_table(name) = literature_table(DW11, name)

"""The four mixes of Table 4."""
const DW11_MIXES = ["OPC", "OPC-L", "OPC-FA", "OPC-FA-L"]

"""The ages of Table 8, in days."""
dw11_days() = sort(unique(ustrip.(dw11_table("pore_solution").time_d)))

# The phases the pastes may form: the phase list of the Lothenbach and Winnefeld
# page (data/phase_lists.toml), the fly ash and the limestone adding no new
# hydrate to it. The C-S-H is the CSHQ model of Cemdata18; the authors used a
# C-S-H with an Al/Si of 0.13 in the fly ash cements, which CSHQ, having no
# aluminum end-member, cannot take.
const DW11_PHASES = "Portland paste (Lothenbach and Winnefeld 2006)"

"""
The chemical system of the phase list, the AFm sulfate and hydroxide declared as
the Guggenheim binary Cemdata18 publishes (the GEMS3K export the replay runs
declares it so), one composition and a second only if the certificate finds the
phase wanting to split.
"""
dw11_system() = phase_list_system(DW11_PHASES, DW11_SUBSTANCES)

"""
    dw11_clinker(mix) -> Material

The clinker, by its Rietveld phases and its minor oxides, each phase dissolving
as Table 7 measured it in `mix`: its degree of reaction at an age is one minus
its content then over its content at 0 days, the content of a phase Table 7 no
longer detects taken as zero, as the column "%OPC reacted" of the table takes
it. The minor oxides (free lime, alkalis, magnesia, sulfate) are released with
the clinker as a whole, at that column's rate.
"""
function dw11_clinker(mix)
    days = sort(unique(ustrip(d) for (m, d) in zip(dw11_table("phase_contents").mix, dw11_table("phase_contents").time_d) if m == mix))
    content(quantity, d) = something(dw11_phase_content(mix, d, quantity), 0.0)
    extents = Dict{String, Any}()
    for phase in ("C3S", "C2S", "C3A", "C4AF")
        c0 = content(phase, first(days))
        extents[phase] = TabulatedExtent(days, [1 - content(phase, d) / c0 for d in days])
    end
    extents["minor oxides"] = TabulatedExtent(days, [content("opc_reacted", d) / 100 for d in days])
    return with_extents(material_template("clinker (De Weerdt 2011)", DW11_DB), extents)
end

"""What Table 7 gives for `quantity` in `mix` at `days` (wt.% of the dry content), or `nothing`."""
function dw11_phase_content(mix, days, quantity)
    pc = dw11_table("phase_contents")
    i = findfirst(k -> pc.mix[k] == mix && ustrip(pc.time_d[k]) == days && pc.quantity[k] == quantity, eachindex(pc.mix))
    return i === nothing ? nothing : ustrip(pc.percent[i])
end

"""The natural gypsum interground with the clinker: its CaSO4·2H2O, the rest left out."""
dw11_gypsum() = Material(
    "natural gypsum (De Weerdt 2011)", :other;
    constituents = [MineralConstituent(DW11_DB["Gp"]; mass_fraction = dw11_value("gypsum_purity_percent") / 100)],
)

"""
    dw11_fly_ash() -> Material

The fly ash: its crystals inert, its glass reacting at the rate of the fit Fig. 7
prints, y = a + b ln(t + c) percent of the fly ash at t days. The reacted fly ash
is taken from the glass alone, which the authors assume dissolves uniformly, so
the glass's own degree of reaction is y over its share of the fly ash.
"""
function dw11_fly_ash()
    fa = material_template("siliceous fly ash (De Weerdt 2011)", DW11_DB)
    share = only(c.mass_fraction for c in fa.constituents if c.name == "glass")
    a, b, c = dw11_value("fly_ash_fit_a"), dw11_value("fly_ash_fit_b"), dw11_value("fly_ash_fit_c")
    rate(t) = clamp((a + b * log(t + c)) / 100 / share, 0.0, 1.0)
    return with_extents(fa, Dict("glass" => rate))
end

"""The limestone powder: its calcite available, the rest of its analysis inert."""
dw11_limestone() = with_extents(material_template("limestone (De Weerdt 2011)", DW11_DB), Dict("minor oxides" => 0.0))

"""
    dw11_recipe(mix) -> Recipe

The paste of `mix` at w/b 0.5 and 20 °C, its binder in the proportions of Table 4,
the OPC being the clinker and 3.7 % of natural gypsum.
"""
function dw11_recipe(mix)
    m = dw11_table("mixes")
    i = findfirst(==(mix), m.mix)
    opc, fa, ls = ustrip(m.opc[i]) / 100, ustrip(m.fly_ash[i]) / 100, ustrip(m.limestone[i]) / 100
    g = dw11_value("gypsum_percent_of_opc") / 100
    parts = Pair{Material, Float64}[dw11_clinker(mix) => opc * (1 - g), dw11_gypsum() => opc * g]
    fa > 0 && push!(parts, dw11_fly_ash() => fa)
    ls > 0 && push!(parts, dw11_limestone() => ls)
    return Recipe(parts...; w_b = dw11_value("water_binder_ratio"), T = dw11_value("curing_temperature") * u"K")
end

"""The concentration of `element` Table 8 gives for `mix` at `days`, in mmol/L (the pH for `"pH"`)."""
function dw11_measured(mix, days, element)
    ps = dw11_table("pore_solution")
    i = findfirst(k -> ps.mix[k] == mix && ustrip(ps.time_d[k]) == days, eachindex(ps.mix))
    i === nothing && return NaN
    return ustrip(getproperty(ps, Symbol(element))[i])
end

"""The portlandite of `mix` at `days`, in wt.% of the dry paste (Table 7; `NaN` when not given)."""
dw11_measured_portlandite(mix, days) = something(dw11_phase_content(mix, days, "portlandite"), NaN)

"""The model's portlandite, in wt.% of the solids of the paste, the unreacted binder included."""
function dw11_portlandite(rs)
    masses = phase_masses(rs)
    return 100 * get(masses, "Portlandite", 0.0) / sum(values(masses))
end
