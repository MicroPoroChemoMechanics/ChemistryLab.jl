# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2026 Jean-François Barthélémy (Cerema, UMR MCD)

# Where the species come from: a database ChemistryLab obtains by name (from
# its publisher, or built from Cemdata18 for the derived ones), or a file of
# the user's in any format `import_database` reads.

"""
    DatabaseChoice

The thermodynamic database of a calculation.
"""
abstract type DatabaseChoice <: Choice end

"A database ChemistryLab obtains by name (`database_info()` lists them)."
struct BuiltinDatabase <: DatabaseChoice
    name::String
end
kind(::Type{BuiltinDatabase}) = "builtin"
fields(::Type{BuiltinDatabase}) = [Field("name", :choice; default = "cemdata18-thermofun.json", source = "databases")]
database_lines(x::BuiltinDatabase) = [
    "# The species and their thermodynamic data, downloaded from the publisher of",
    "# the database at the first use (or built from Cemdata18) and kept afterwards.",
    "substances = build_species(datapath($(repr(x.name))); verbose = false)",
]

"A database file of the user's, in a format `import_database` reads."
struct ImportedDatabase <: DatabaseChoice
    path::String
    format::String
end
kind(::Type{ImportedDatabase}) = "imported"
fields(::Type{ImportedDatabase}) = [
    Field("path", :text; default = ""),
    Field("format", :choice; default = "auto", options = ["auto", "thermofun", "phreeqc", "gwb", "eq36", "reaktoro"]),
]
database_lines(x::ImportedDatabase) = [
    "# The species of a database file, read in its own format.",
    "_, table, _ = import_database($(repr(x.path)); format = :$(x.format))",
    "substances = build_species(table)",
]

members(::Type{DatabaseChoice}) = (BuiltinDatabase, ImportedDatabase)

database_problems(x::BuiltinDatabase) = bound_problems(x, "the database")
function database_problems(x::ImportedDatabase)
    isempty(strip(x.path)) || return Problem[]
    return [Problem("no_database_file", "Choose the database file to import.")]
end
