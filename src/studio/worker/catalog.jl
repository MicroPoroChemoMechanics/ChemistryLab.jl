# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2026 Jean-François Barthélémy (Cerema, UMR MCD)

# What the page offers to choose from, read from ChemistryLab rather than
# copied into the studio: the databases and their status, the published phase
# lists, the material templates with their constituents, the solid solutions,
# and the species of a database.

# The species of a database, read once (the first read compiles for a while).
# A species vector's type is huge; `Any` spares compiling a Dict for it.
const SUBSTANCES = Dict{String, Any}()
const SUBSTANCES_LOCK = ReentrantLock()

function substances_of(database::AbstractString)
    return lock(SUBSTANCES_LOCK) do
        get!(() -> build_species(datapath(database); verbose = false), SUBSTANCES, database)
    end
end

"A template's constituents, built against the whole database (every species is there)."
function template_entry(name::AbstractString, db)
    m = material_template(name, db)
    return Dict{String, Any}(
        "name" => name, "kind" => String(m.kind), "source" => m.source === nothing ? nothing : string(m.source),
        "constituents" => [
            Dict{String, Any}(
                "name" => c.name, "mass_fraction" => Float64(c.mass_fraction),
                "species" => c isa MineralConstituent ? symbol(c.species) : nothing,
                # known by its oxides: in a supplementary material, a glass a
                # rate law can dissolve (`with_glass_species`)
                "oxides" => c isa OxideConstituent,
            ) for c in m.constituents
        ],
    )
end

"The solid solutions ChemistryLab ships, as its file declares them."
function solid_solution_entries()
    toml = TOML.parsefile(datapath("solid_solutions.toml"))
    return [
        Dict{String, Any}(
            "name" => e["name"], "end_members" => String.(e["end_members"]), "model" => get(e, "model", "ideal"),
            "database" => get(e, "database", nothing), "source" => string(get(e, "source", "")),
        ) for e in get(toml, "solid_solution", Any[])
    ]
end

"The published works ChemistryLab holds data from: the citation and the DOI of each, by key."
function references()
    return Dict{String, Any}(
        key => Dict{String, Any}("citation" => get(r.source, "citation", key), "doi" => get(r.source, "doi", nothing))
            for (key, r) in ((k, literature(k)) for k in available_literature())
    )
end

"""
    catalog() -> Dict{String, Any}

The choices ChemistryLab offers: the databases with their status, the phase
lists with their phases, the material templates with their constituents, the
solid solutions, the published works the data come from, and the versions.
"""
function catalog()
    db = Dict(symbol(s) => s for s in substances_of("cemdata18-thermofun.json"))
    lists = map(phase_lists()) do n
        pl = phase_list(n)
        Dict{String, Any}(
            "name" => n, "source" => string(something(pl.source, "")), "database" => pl.database,
            "reactants" => pl.reactants, "products" => pl.products, "solid_solutions" => pl.solid_solutions,
        )
    end
    return Dict{String, Any}(
        "databases" => [
            Dict{String, Any}("name" => d.name, "status" => String(d.status), "license" => string(d.license))
                for d in database_info(devnull)
        ],
        "phase_lists" => lists,
        "templates" => [template_entry(n, db) for n in material_templates()],
        "solid_solutions" => solid_solution_entries(),
        "references" => references(),
        "versions" => versions(),
    )
end

# The cement notation of a species, when it has one (not for ions).
function cement_notation(s)
    return try
        unicode(cemformula(CemSpecies(s)))
    catch err
        # a species whose elements are not all cement oxides, or an ion
        err isa ArgumentError || err isa ErrorException || err isa KeyError || rethrow()
        ""
    end
end

"""
    species_catalog(database) -> Vector{Dict}

The species of a database, for the page's picker: symbol, formula, cement
notation, aggregate state, elements, molar mass and temperature range.
"""
function species_catalog(database::AbstractString)
    return map(substances_of(database)) do s
        lo, hi = temperature_range(s)
        Dict{String, Any}(
            "symbol" => symbol(s), "formula" => unicode(formula(s)), "cement" => cement_notation(s),
            "state" => lowercase(replace(string(aggregate_state(s)), "AS_" => "")),
            "elements" => [string(e) for e in keys(atoms(s)) if e !== :Zz],
            "molar_mass" => finite(ustrip(us"g/mol", s[:M])),
            "T_min" => finite(lo), "T_max" => finite(hi),
        )
    end
end

"The versions of Julia and of the packages the calculations use."
versions() = Dict{String, Any}(
    "julia" => string(VERSION),
    "ChemistryLab" => string(pkgversion(ChemistryLab)),
    "OptimaSolver" => string(pkgversion(OptimaSolver)),
)

"""
    check_report() -> Vector{String}

The lines `clstudio(check = true)` prints: versions and the status of each
database.
"""
function check_report()
    v = versions()
    lines = ["Julia $(v["julia"]), ChemistryLab $(v["ChemistryLab"]), OptimaSolver $(v["OptimaSolver"])."]
    for d in database_info(devnull)
        push!(lines, "  database $(rpad(d.name, 40)) $(d.status)")
    end
    return lines
end
