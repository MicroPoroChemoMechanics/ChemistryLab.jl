# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2026 Jean-François Barthélémy (Cerema, UMR MCD)

# What the page offers to choose from, read from ChemistryLab rather than
# copied into the studio: the published phase lists, the material templates
# with their constituents, the databases and their status.

"The system of a phase list, built once per list (the first build compiles for a while)."
const SYSTEMS = Dict{String, Any}()   # a `ChemicalSystem`'s type is huge; `Any` spares its compilation
const SYSTEMS_LOCK = ReentrantLock()

function system_of(list::AbstractString)
    return lock(SYSTEMS_LOCK) do
        get!(SYSTEMS, list) do
            pl = phase_list(list)
            phase_list_system(list, build_species(datapath(pl.database); verbose = false))
        end
    end
end

"A template's constituents, built in the system of `list`."
function template_entry(name::AbstractString, list::AbstractString)
    entry = Dict{String, Any}("name" => name)
    m = try
        material_template(name, system_of(list))
    catch err
        # a template whose species the list does not hold: listed, with the reason
        err isa KeyError || err isa ArgumentError || err isa ErrorException || rethrow()
        entry["unavailable"] = err isa KeyError ?
            "needs the species $(err.key), which the phase list \"$list\" does not declare" : sprint(showerror, err)
        return entry
    end
    entry["kind"] = String(m.kind)
    entry["source"] = m.source === nothing ? nothing : string(m.source)
    entry["constituents"] = [
        Dict{String, Any}("name" => c.name, "mass_fraction" => Float64(c.mass_fraction)) for c in m.constituents
    ]
    return entry
end

"""
    catalog() -> Dict{String, Any}

The choices ChemistryLab offers: phase lists with their source, material
templates with their constituents (built in the system of the first list),
and the databases with their status.
"""
function catalog()
    lists = phase_lists()
    first_list = first(lists)
    return Dict{String, Any}(
        "phase_lists" => [Dict{String, Any}("name" => n, "source" => string(something(phase_list(n).source, ""))) for n in lists],
        "templates" => [template_entry(n, first_list) for n in material_templates()],
        "databases" => [
            Dict{String, Any}("name" => d.name, "status" => String(d.status), "license" => string(d.license))
                for d in database_info(devnull)
        ],
        "activity_models" => ["cemdata18_KOH", "cemdata18_NaOH"],
        "versions" => versions(),
    )
end

"The versions of Julia and of the packages the calculations use."
versions() = Dict{String, Any}(
    "julia" => string(VERSION),
    "ChemistryLab" => string(pkgversion(ChemistryLab)),
    "OptimaSolver" => string(pkgversion(OptimaSolver)),
)

"""
    check_report() -> Vector{String}

The lines `clstudio(check = true)` prints: versions, the status of each
database, and the free memory as this process sees it.
"""
function check_report()
    v = versions()
    lines = ["Julia $(v["julia"]), ChemistryLab $(v["ChemistryLab"]), OptimaSolver $(v["OptimaSolver"])."]
    for d in database_info(devnull)
        push!(lines, "  database $(rpad(d.name, 40)) $(d.status)")
    end
    return lines
end
