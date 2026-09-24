# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

using DynamicQuantities
using SHA

"""
    build_sit_parameters(path; source = nothing) -> SITParameters

Read the `SIT` / `-epsilon` block of a PHREEQC-format database into a
[`SITParameters`](@ref) compilation.

Reading rather than transcribing is deliberate. A compilation of several hundred
interaction coefficients copied by hand is a transcription, which is the error
class every generator in `test/reference/` exists to remove, and this package
has already found two standard energies that had drifted that way.

**No such compilation ships with this package.** The one PHREEQC distributes in
`sit.dat` is the ANDRA/RWM *ThermoChimie* database, which is not USGS-authored
and whose terms are its own — so this reads a file the caller already has rather
than redistributing one. `source` defaults to the file's name and a truncated SHA-256 of its
contents, so a result can always say which compilation it came from.

# Format

```
SIT
-epsilon
    Na+     Cl-     0.03
    ...
```

Three whitespace-separated fields per line: two species symbols and the
coefficient in kg/mol. Blank lines and `#` comments are skipped; the block ends
at the first line that is neither.

# Example

```julia
p = build_sit_parameters("/path/to/sit.dat")
model = SITActivityModel(; parameters = p)
missing_epsilon_pairs(cs, model)     # what the compilation does not cover
```
"""
function build_sit_parameters(path::AbstractString; source = nothing)
    isfile(path) || throw(ArgumentError("no such database: $path"))
    text = read(path, String)
    digest = bytes2hex(sha256(text))
    src = source === nothing ?
        "$(basename(path)) sha256 $(first(digest, 12))" : String(source)

    pairs = Pair{Tuple{String, String}, Float64}[]
    in_block = false
    seen_epsilon = false
    for raw in split(text, '\n')
        line = strip(first(split(raw, '#')))
        if !in_block
            uppercase(line) == "SIT" && (in_block = true)
            continue
        end
        isempty(line) && continue
        if startswith(line, '-')
            seen_epsilon = lowercase(line) == "-epsilon"
            seen_epsilon && continue
            # another sub-keyword of the same block: stop reading coefficients
            seen_epsilon = false
            continue
        end
        fields = split(line)
        if length(fields) == 3 && seen_epsilon
            v = tryparse(Float64, fields[3])
            v === nothing && break
            push!(pairs, (String(fields[1]), String(fields[2])) => v)
        else
            break          # the next keyword block
        end
    end
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
