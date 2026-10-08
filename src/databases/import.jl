# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

# ── Importing a database of any format ───────────────────────────────────────

# The format of a file, from its name, then from its first lines.
function _database_format(path)
    file = lowercase(basename(path))
    ext = lowercase(splitext(file)[2])
    ext == ".json" && return :thermofun
    ext in (".yaml", ".yml") && return :reaktoro
    ext == ".tdat" && return :gwb
    startswith(file, "data0") && return :eq36
    head = lowercase(join(Iterators.take(eachline(path), 400), "\n"))
    occursin("dataset of thermodynamic data for gwb", head) && return :gwb
    startswith(strip(head), "data0") && return :eq36
    (occursin("solution_master_species", head) || occursin("solution_species", head)) && return :phreeqc
    throw(ArgumentError("$(basename(path)): the format of this file is not recognized; give it as `format`"))
end

"""
    import_database(path; format = :auto) -> (df_elements, df_substances, df_reactions)

Read a thermodynamic database of any of the formats the package reads, into the
three tables [`read_thermofun_database`](@ref) returns for a ThermoFun file:

| `format` | files | read by |
|:--|:--|:--|
| `:phreeqc` | PHREEQC databases (`.dat`) | [`read_phreeqc_database`](@ref) |
| `:gwb` | thermo datasets of The Geochemist's Workbench (`.tdat`) | [`read_gwb_database`](@ref) |
| `:eq36` | data0 files of EQ3/6 | [`read_eq36_database`](@ref) |
| `:thermofun` | ThermoFun databases (`.json`) | [`read_thermofun_database`](@ref) |
| `:reaktoro` | Reaktoro databases (`.yaml`) | [`read_reaktoro_database`](@ref) |

`:auto` decides by the name of the file, then by its first lines.
[`build_species`](@ref) builds the species of `df_substances` in every case.
`path` is a file, or the name of a database `datapath` knows.

# Example

```julia
_, substances, _ = import_database(datapath("phreeqc.dat"))
species = build_species(substances, ["H2O", "H+", "OH-", "Ca+2", "CO3-2", "HCO3-", "CO2", "Calcite"])
```
"""
function import_database(path::AbstractString; format::Symbol = :auto)
    p = resolve_data_path(path)
    f = format === :auto ? _database_format(p) : format
    f === :phreeqc && return read_phreeqc_database(p)
    f === :gwb && return read_gwb_database(p)
    f === :eq36 && return read_eq36_database(p)
    f === :reaktoro && return read_reaktoro_database(p)
    f === :thermofun && return read_thermofun_database(p)
    throw(ArgumentError("format is :auto, :phreeqc, :gwb, :eq36, :thermofun or :reaktoro; got :$f"))
end
