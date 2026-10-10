# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

using DynamicQuantities

"""
    build_sit_parameters(path; source = nothing) -> SITParameters

Read the `SIT` / `-epsilon` block of a PHREEQC-format database into a
[`SITParameters`](@ref) compilation.

Reading rather than transcribing is deliberate. A compilation of several hundred
interaction coefficients copied by hand is a transcription, which is the error
class every generator in `test/reference/` exists to remove, and this package
has already found two standard energies that had drifted that way.

**No database is stored in this package.** The compilation PHREEQC distributes
in `sit.dat`, the ANDRA/RWM *ThermoChimie* database, is obtained by
`datapath("sit.dat")` from the PHREEQC distribution; any other file the caller has
is read the same way. `source` defaults to the file's name and a truncated
SHA-256 of its contents, so a result can always say which compilation it came
from.

# Format

```
SIT
-epsilon
    Na+     Cl-     0.03
    ...
```

Three whitespace-separated fields per line: two species and the coefficient in
kg/mol. The block ends at the next keyword. The species are named as in the
database and keyed by the rule of every reader: `SO4--` is `SO4-2`, a neutral
species gets `@`. A line of another shape is not read, and a warning lists it.

# Example

```julia
p = build_sit_parameters("/path/to/sit.dat")
model = SITActivityModel(; parameters = p)
missing_epsilon_pairs(cs, model)     # what the compilation does not cover
```
"""
function build_sit_parameters(path::AbstractString; source = nothing)
    isfile(path) || throw(ArgumentError("no such database: $path"))
    src = source === nothing ? _source_tag(path) : String(source)
    pairs = Pair{Tuple{String, String}, Float64}[]
    unread = String[]
    for b in phreeqc_blocks(path)
        b.keyword == "SIT" || continue
        epsilon = false
        for (n, line) in b.lines
            if startswith(line, '-')
                # `-epsilon` opens the coefficients; another option closes them.
                epsilon = lowercase(line) == "-epsilon"
                continue
            end
            epsilon || continue
            fields = split(line)
            v = length(fields) == 3 ? tryparse(Float64, fields[3]) : nothing
            if v === nothing
                push!(unread, "line $n: `$line`")
                continue
            end
            a, c = _solute_symbol.(_phreeqc_name.(fields[1:2]))
            push!(pairs, (a, c) => v)
        end
    end
    isempty(unread) || @warn "build_sit_parameters: $(basename(path)) has ε lines of another shape than `species species ε`; not read: $(join(unread, "; "))"
    isempty(pairs) && throw(
        ArgumentError(
            "no SIT ε found in $path. A PHREEQC database carries them in a " *
                "`SIT` block under `-epsilon`; `phreeqc.dat` has none, and " *
                "`sit.dat` does."
        ),
    )
    # Read out of a named, hashed database: published, not estimated and not
    # this package's own measurement.
    return SITParameters(pairs; source = src, kind = PROV_PUBLISHED)
end
