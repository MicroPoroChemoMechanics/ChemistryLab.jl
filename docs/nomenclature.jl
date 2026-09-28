# The nomenclature of the documentation, read from `nomenclature.toml`: the page
# `nomenclature.md` renders it, and `make.jl` writes it as the JSON the hints
# over the equations read (`src/.vitepress/nomenclature.json`, not versioned).
# The value of a physical constant is read from the library, never typed.

using ChemistryLab
using DynamicQuantities
using JSON
using Markdown
using Printf
using TOML

const NOMENCLATURE_FILE = joinpath(@__DIR__, "nomenclature.toml")

# The constants an entry may name, each with the unit it is shown in. They are
# those of ChemistryLab (`src/utils/constants.jl`), and for the two it does not
# name, the same CODATA values of DynamicQuantities it reads the others from.
const NOMENCLATURE_CONSTANTS = Dict(
    "R" => (ChemistryLab.R_GAS_Q, us"J/(mol*K)", "J/(mol K)"),
    "F" => (ChemistryLab.FARADAY_Q, us"C/mol", "C/mol"),
    "N_A" => (ChemistryLab.AVOGADRO_Q, us"1/mol", "1/mol"),
    "epsilon_0" => (ChemistryLab.VACUUM_PERMITTIVITY_Q, us"F/m", "F/m"),
    "k_B" => (DynamicQuantities.Constants.k_B, us"J/K", "J/K"),
    "e" => (DynamicQuantities.Constants.e, us"C", "C"),
)

const NOMENCLATURE_GROUPS = [
    "Physical constants", "State and thermodynamic functions", "Amounts and composition",
    "Aqueous solutions", "Surfaces", "Solid solutions", "Equilibrium and the solver",
    "Kinetics and hydration", "Cement chemist notation",
]

"""
    nomenclature_entries() -> Vector{NamedTuple}

Every entry of `nomenclature.toml`, with an `id` and, for a physical constant,
its `value` and `unit` read from the library.
"""
function nomenclature_entries()
    raw = TOML.parsefile(NOMENCLATURE_FILE)["symbol"]
    return map(enumerate(raw)) do (k, e)
        value, unit = nothing, get(e, "unit", "")
        if haskey(e, "constant")
            q, u, shown = NOMENCLATURE_CONSTANTS[e["constant"]]
            value, unit = @sprintf("%.10g", ustrip(u, q)), shown
        end
        (;
            id = "nomen-$k", group = e["group"], tex = e["tex"], label = e["label"], name = e["name"],
            unit, value, pages = String.(get(e, "pages", String[])), script = get(e, "script", false),
        )
    end
end

"""
    write_nomenclature_json(path)

The entries, as the JSON the plugin that typesets the formulas and the hints of
the theme read.
"""
function write_nomenclature_json(path)
    entries = [
        Dict(
            "id" => e.id, "tex" => e.tex, "label" => e.label, "name" => e.name, "unit" => e.unit,
            "value" => e.value, "pages" => e.pages, "script" => e.script,
        ) for e in nomenclature_entries()
    ]
    mkpath(dirname(path))
    open(io -> JSON.print(io, entries, 1), path, "w")
    return path
end

# The HTML the hints show (<sub>, <sup>, <b>), as the text of a Markdown table:
# Unicode subscripts and superscripts where every character has one, `_x` and
# `^x` otherwise.
const _SUBSCRIPTS = Dict(zip("0123456789+-−=()aehijklmnoprstuvx", "₀₁₂₃₄₅₆₇₈₉₊₋₋₌₍₎ₐₑₕᵢⱼₖₗₘₙₒₚᵣₛₜᵤᵥₓ"))
const _SUPERSCRIPTS = Dict(zip("0123456789+-−=()ni", "⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻⁻⁼⁽⁾ⁿⁱ"))
function _plain(s::AbstractString)
    shift(t, table, mark) = all(c -> haskey(table, c), t) ? join(table[c] for c in t) : mark * t
    s = replace(s, r"<sub>(.*?)</sub>" => m -> shift(match(r"<sub>(.*?)</sub>", m)[1], _SUBSCRIPTS, "_"))
    s = replace(s, r"<sup>(.*?)</sup>" => m -> shift(match(r"<sup>(.*?)</sup>", m)[1], _SUPERSCRIPTS, "^"))
    return replace(s, r"</?b>" => "")
end

# A page of the documentation as a link from `nomenclature.md`.
_nomenclature_page_link(p) = endswith(p, "/") ? "`$p`" : "[`$(first(splitext(p)))`]($p)"

"""
    nomenclature_markdown() -> Markdown.MD

The tables of the nomenclature page, one per group: the symbol, its meaning (with
the value of a constant), its unit, and the pages a meaning is restricted to.
"""
function nomenclature_markdown()
    entries = nomenclature_entries()
    io = IOBuffer()
    for g in NOMENCLATURE_GROUPS
        rows = [e for e in entries if e.group == g]
        isempty(rows) && continue
        println(io, "## ", g, "\n")
        println(io, "| Symbol | Meaning | Unit | Where |")
        println(io, "|:--|:--|:--|:--|")
        for e in rows
            meaning = _plain(e.value === nothing ? e.name : "$(e.name), $(e.value) $(e.unit)")
            where = isempty(e.pages) ? "" : join(_nomenclature_page_link.(e.pages), ", ")
            println(io, "| ", "``", e.tex, "``", " | ", meaning, " | ", _plain(e.unit), " | ", where, " |")
        end
        println(io)
    end
    return Markdown.parse(String(take!(io)))
end
