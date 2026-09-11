# SPDX-License-Identifier: LGPL-2.1-or-later
#
# Reading what `precompute_docs.jl` wrote.
#
# The coupled runs of the documentation cost tens of minutes each, so they are
# computed once by `scripts/precompute_docs.jl` and the pages plot the result.
# This is the seam between the two: a few lines of plain Julia, with no
# dependency beyond the standard library, so that a reader can see exactly what
# the page is reading and re-run the producer if they doubt it.

using DelimitedFiles
using Printf

const PRECOMPUTED_DIR =
    joinpath(@__DIR__, "..", "docs", "src", "assets", "precomputed")

"""
    read_precomputed(name) -> (; columns, provenance, header)

Read one precomputed table. `columns` maps a column name to its vector,
`provenance` is the `#` block the producer wrote — package version, commit,
sampling, and the composition the run was made at — and `header` is the column
order.

The provenance is returned rather than discarded because a stored number is
only worth what its traceability is: a figure drawn from one of these files can
say which commit produced it.
"""
function read_precomputed(name::AbstractString)
    path = joinpath(PRECOMPUTED_DIR, endswith(name, ".csv") ? name : name * ".csv")
    isfile(path) || error(
        "no precomputed file at $path. Run\n" *
            "    julia --project=docs scripts/precompute_docs.jl\n" *
            "to produce it. It takes a few minutes and writes every file the " *
            "documentation reads.",
    )
    provenance = String[]
    for line in eachline(path)
        startswith(line, "#") || break
        push!(provenance, strip(lstrip(line, ['#', ' '])))
    end
    raw = readdlm(path, ','; comments = true, comment_char = '#')
    header = String.(vec(raw[1, :]))
    data = Float64.(raw[2:end, :])
    columns = Dict(header[j] => data[:, j] for j in eachindex(header))
    return (; columns, provenance, header)
end

"""
    phase_families(tbl) -> Vector{String}

The phase-family columns of a phase-history table, in a fixed order chosen so a
stacked plot reads from the anhydrous grains at the bottom to the pore water at
the top — the order the classical hydration diagram uses.

Families absent from the table are skipped, so the same order serves a paste
with limestone and one without.
"""
function phase_families(tbl)
    order = [
        "anhydrous", "gypsum", "calcite",
        "C-S-H", "CH", "AFt", "AFm", "hydrogarnet", "FH3",
        "water", "void",
    ]
    return [f for f in order if haskey(tbl.columns, f)]
end

"""
    family_colors() -> Dict{String,Symbol}

One color per phase family, fixed here so that every figure in the
documentation uses the same one for the same phase — a reader should not have to
re-read the legend between two plots.
"""
function family_colors()
    return Dict(
        "anhydrous" => :grey40, "gypsum" => :khaki, "calcite" => :tan,
        "C-S-H" => :steelblue, "CH" => :seagreen, "AFt" => :orchid,
        "AFm" => :mediumpurple, "hydrogarnet" => :peru, "FH3" => :indianred,
        "water" => :lightskyblue, "void" => :white,
    )
end
