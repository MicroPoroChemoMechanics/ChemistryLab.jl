# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

using TOML

# ── Phase lists of published pastes ───────────────────────────────────────────

"""
    phase_list(name) -> NamedTuple
    phase_lists() -> Vector{String}

A phase list of `data/phase_lists.toml`: the phases a paste may form, written
after the assemblage of the paper it names. Its fields are `name`, `database`
(the file its symbols are records of), `reactants` and `products` (pure phases),
`solid_solutions` (names in `data/solid_solutions.toml`), `instances` (what the
list declares a solid solution with, where it departs from that file),
`exclude_aqueous` and `source`. `phase_lists()` gives the names;
[`phase_list_system`](@ref) builds the system.
"""
function phase_list(name::AbstractString)
    e = _phase_list_entry(name)
    strings(key) = String.(get(e, key, String[]))
    return (;
        name = e["name"], database = e["database"],
        reactants = strings("reactants"), products = strings("products"),
        solid_solutions = strings("solid_solutions"),
        instances = Dict{String, Any}(k => (v == "auto" ? :auto : Int(v)) for (k, v) in get(e, "instances", Dict())),
        exclude_aqueous = strings("exclude_aqueous"), source = get(e, "source", nothing),
    )
end

phase_lists() = [e["name"] for e in _phase_list_entries()]

"""
    phase_list_system(name, substances; add = String[], remove = String[],
                      primaries = CEMDATA_PRIMARIES) -> ChemicalSystem

The chemical system of the phase list `name` ([`phase_list`](@ref)): its
reactants and products as pure phases, its solid solutions as
`data/solid_solutions.toml` declares them (with the `instances` the list gives),
and the aqueous species of `substances` its elements allow, those of
`exclude_aqueous` left out. `substances` are the species of the list's
database, `build_species(datapath(phase_list(name).database))`.

A calculation that departs from the list says so: `add` names pure phases to
declare beside the list's, `remove` pure phases of the list to leave out.

```julia
substances = build_species(datapath("cemdata18-thermofun.json"))
cs = phase_list_system("Portland paste (Lothenbach and Winnefeld 2006)", substances;
                       add = ["C4AClH10"])   # and Friedel's salt
```
"""
function phase_list_system(
        list::AbstractString, substances; add = String[], remove = String[],
        primaries = CEMDATA_PRIMARIES,
    )
    pl = phase_list(list)
    listed = vcat(pl.reactants, pl.products)
    stray = setdiff(remove, listed)
    isempty(stray) || throw(ArgumentError("phase_list_system: \"$list\" lists no pure phase $(join(stray, ", ")) to remove."))
    pure = vcat(setdiff(listed, remove), setdiff(add, listed))
    db = Dict(symbol(s) => s for s in substances)
    ss = [
        p for p in build_solid_solutions(datapath("solid_solutions.toml"), db; instances = pl.instances)
            if name(p) in pl.solid_solutions
    ]
    absent = setdiff(pl.solid_solutions, name.(ss))
    isempty(absent) || throw(
        ArgumentError(
            "phase_list_system: the solid solutions $(join(absent, ", ")) of \"$list\" could not be built " *
                "from these substances; the list is written for $(pl.database)."
        ),
    )
    members = unique(symbol(m) for p in ss for m in end_members(p))
    unknown = setdiff(vcat(pure, members), keys(db))
    isempty(unknown) || throw(
        ArgumentError(
            "phase_list_system: no species $(join(unknown, ", ")) among these substances; " *
                "\"$list\" is written for $(pl.database)."
        ),
    )
    sp = speciation(
        substances, vcat(pure, members);
        aggregate_state = [AS_AQUEOUS], exclude_species = pl.exclude_aqueous,
    )
    return ChemicalSystem(sp, primaries; solid_solutions = ss)
end

_phase_list_entries() = get(TOML.parsefile(resolve_data_path("phase_lists.toml")), "phase_list", Any[])
function _phase_list_entry(name)
    for e in _phase_list_entries()
        e["name"] == name && return e
    end
    throw(KeyError("no phase list \"$name\"; the lists are: " * join(phase_lists(), "; ")))
end
