# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2026 Jean-François Barthélémy (Cerema, UMR MCD)

# Which species may form: a published phase list, with the departures the user
# makes from it (phases added or removed, another model of the C-S-H gel), or a
# selection of the user's own.

"""
    SystemChoice

The chemical system of a calculation.
"""
abstract type SystemChoice <: Choice end

"The models of the C-S-H gel; one at most may be declared."
const GEL_MODELS = ("CSHQ", "CNASH_ss", "CSH3T", "ECSH1", "ECSH2", "CASH+", "CASH+NK", "CASH+ext")

"A published phase list (`phase_lists()`), with the user's departures from it."
struct PhaseListSystem <: SystemChoice
    phase_list::String
    add::Vector{String}
    remove::Vector{String}
    gel::String
    exclude_aqueous::Vector{String}
end
kind(::Type{PhaseListSystem}) = "phase_list"
fields(::Type{PhaseListSystem}) = [
    Field("phase_list", :choice; default = "Portland paste (Lothenbach and Winnefeld 2006)", source = "phase_lists"),
    Field("add", :texts; default = String[]),
    Field("remove", :texts; default = String[]),
    Field("gel", :choice; default = "CSHQ", options = collect(GEL_MODELS)),
    Field("exclude_aqueous", :texts; default = String[], level = :expert),
]
function system_lines(x::PhaseListSystem)
    kws = String[]
    isempty(x.add) || push!(kws, "add = $(repr(x.add))")
    isempty(x.remove) || push!(kws, "remove = $(repr(x.remove))")
    x.gel == "CSHQ" || push!(kws, "replace = Dict(\"CSHQ\" => $(repr(x.gel)))")
    isempty(x.exclude_aqueous) || push!(kws, "exclude_aqueous = $(repr(x.exclude_aqueous))")
    tail = isempty(kws) ? "" : "; " * join(kws, ", ")
    lines = [
        "# The phases the paste may form, from a published phase list",
        "# (`phase_lists()` lists them)$(isempty(kws) ? "." : ", with the departures written after it.")",
        "cs = phase_list_system($(repr(x.phase_list)), substances$tail)",
    ]
    return lines
end

"A selection of the user's: pure phases, solid solutions, and the aqueous species their elements allow."
struct CustomSystem <: SystemChoice
    phases::Vector{String}
    solid_solutions::Vector{String}
    exclude_aqueous::Vector{String}
end
kind(::Type{CustomSystem}) = "custom"
fields(::Type{CustomSystem}) = [
    Field("phases", :texts; default = ["C3S", "C2S", "C3A", "C4AF", "Gp", "Cal", "Portlandite", "ettringite", "monosulphate12", "monocarbonate", "hydrotalcite", "FeOOHmic"]),
    Field("solid_solutions", :texts; default = ["CSHQ"]),
    Field("exclude_aqueous", :texts; default = ["H2@", "O2@", "CH4@"], level = :expert),
]
system_lines(x::CustomSystem) = [
    "# The phases chosen one by one: pure phases, solid solutions from ChemistryLab's",
    "# file of them, and the aqueous species their elements allow.",
    "db = Dict(symbol(s) => s for s in substances)",
    "solid_solutions = [p for p in build_solid_solutions(datapath(\"solid_solutions.toml\"), db) if ChemistryLab.name(p) in $(repr(x.solid_solutions))]",
    "members = unique(symbol(m) for p in solid_solutions for m in end_members(p))",
    "species = speciation(",
    "    substances, vcat($(repr(x.phases)), members);",
    "    aggregate_state = [AS_AQUEOUS], exclude_species = $(repr(x.exclude_aqueous)),",
    ")",
    "cs = ChemicalSystem(species, CEMDATA_PRIMARIES; solid_solutions)",
]

members(::Type{SystemChoice}) = (PhaseListSystem, CustomSystem)

system_problems(x::PhaseListSystem) = Problem[]
function system_problems(x::CustomSystem)
    out = Problem[]
    isempty(x.phases) && isempty(x.solid_solutions) && push!(out, Problem("empty_system", "The system has no phase."))
    gels = filter(in(GEL_MODELS), x.solid_solutions)
    length(gels) > 1 && push!(
        out, Problem(
            "two_gels", "$(join(gels, " and ")) are models of one C-S-H gel: declare one of them.",
            Dict{String, Any}("gels" => join(gels, ", ")),
        )
    )
    return out
end
