# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2026 Jean-François Barthélémy (Cerema, UMR MCD)

# The rows of an input table: a formula or a name, an amount and its unit, and
# the notation the formula is written in. Whatever the user pastes from a
# spreadsheet, a CSV file or a table copied from a PDF is read here, in the
# server, so that the reading is the same everywhere and can be tested; what a
# formula means is read by ChemistryLab in the compute process.

"""
    AmountRow(formula, amount, unit, notation)

One row of an input table: a formula (`CaO`, `C3S`) or a name (a species of the
database, a material), an amount, its unit (one of `UNITS`), and the notation
of the formula: `"auto"`, `"mendeleev"` or `"cement"`.
"""
struct AmountRow
    formula::String
    amount::Float64
    unit::String
    notation::String
end

"The units an amount may be given in, by what they measure."
const UNITS = (
    mass = ("g", "kg", "mg", "%"),
    amount = ("mol", "mmol", "µmol"),
    concentration = ("mol/kg", "mmol/kg", "mol/L", "mmol/L"),
)

const ALL_UNITS = Tuple(Iterators.flatten(values(UNITS)))

const NOTATIONS = ("auto", "mendeleev", "cement")

# The spellings a pasted unit may take, by the unit the studio uses.
const UNIT_ALIASES = Dict(
    "g" => "g", "kg" => "kg", "mg" => "mg",
    "%" => "%", "wt%" => "%", "wt.%" => "%", "mass%" => "%", "%m" => "%", "% m" => "%",
    "mol" => "mol", "mmol" => "mmol", "µmol" => "µmol", "μmol" => "µmol", "umol" => "µmol",
    "mol/kg" => "mol/kg", "mol/kgw" => "mol/kg", "mmol/kg" => "mmol/kg", "mmol/kgw" => "mmol/kg",
    "mol/l" => "mol/L", "m" => "mol/L", "mmol/l" => "mmol/L", "mm" => "mmol/L",
)

"""
    canonical_unit(s) -> Union{Nothing, String}

The unit a pasted spelling stands for (`"wt%"` is `"%"`, `"mmol/l"` is
`"mmol/L"`), or `nothing` when it names none.
"""
canonical_unit(s::AbstractString) = get(UNIT_ALIASES, lowercase(strip(s)), nothing)

function AmountRow(formula::AbstractString, amount::Real, unit::AbstractString, notation::AbstractString = "auto")
    unit in ALL_UNITS || throw(ArgumentError("unknown unit $(repr(unit)); the units are $(join(ALL_UNITS, ", "))"))
    notation in NOTATIONS || throw(ArgumentError("unknown notation $(repr(notation))"))
    return AmountRow(String(formula), Float64(amount), String(unit), String(notation))
end

to_dict(r::AmountRow) = Dict{String, Any}("formula" => r.formula, "amount" => r.amount, "unit" => r.unit, "notation" => r.notation)
plain(r::AmountRow) = to_dict(r)

"""
    amount_row(d) -> AmountRow

The row a JSON object describes.
"""
function amount_row(d)
    d isa AbstractDict || throw(ArgumentError("a row must be an object, got $(repr(d))"))
    formula, amount = get(d, "formula", nothing), get(d, "amount", nothing)
    formula isa AbstractString || throw(ArgumentError("a row needs a formula, got $(repr(formula))"))
    (amount isa Real && !(amount isa Bool)) || throw(ArgumentError("the amount of $formula must be a number, got $(repr(amount))"))
    unit = get(d, "unit", "%")
    unit isa AbstractString || throw(ArgumentError("the unit of $formula must be text"))
    notation = get(d, "notation", "auto")
    notation isa AbstractString || throw(ArgumentError("the notation of $formula must be text"))
    return AmountRow(formula, amount, unit, notation)
end

"""
    parse_number(s; decimal_comma) -> Union{Nothing, Float64}

A number as a spreadsheet or a French text writes it: `12.5`, `12,5` (when the
comma is the decimal mark), `1 234,5`, `1e-3`. `nothing` when it is not one.
"""
function parse_number(s::AbstractString; decimal_comma::Bool)
    t = replace(strip(s), r"[\s  ]" => "")
    decimal_comma && (t = replace(t, "," => "."))
    return tryparse(Float64, t)
end

"""
    split_cells(line, separator) -> Vector{String}

The cells of a line: split at the separator, or at runs of blanks for a table
copied from a PDF (`separator === nothing`).
"""
function split_cells(line::AbstractString, separator::Union{Nothing, Char})
    cells = separator === nothing ? split(line) : split(line, separator)
    return String.(strip.(cells))
end

"""
    detect_separator(lines) -> Union{Nothing, Char}

The separator of pasted lines: a tab (a spreadsheet), a semicolon (a French
CSV, whose decimal mark is the comma), a comma (an English CSV), or `nothing`
for cells separated by blanks.
"""
function detect_separator(lines)
    any(contains('\t'), lines) && return '\t'
    any(contains(';'), lines) && return ';'
    # a comma separates cells only if every line splits into at least two at it
    # and no line reads as "name  12,5" (a decimal comma between blanks)
    if all(l -> count(==(','), l) >= 1, lines) && !any(l -> occursin(r"^\S+\s+[-+]?\d+,\d+(\s|$)", l), lines)
        return ','
    end
    return nothing
end

"""
    parse_table(text; unit = "%") -> NamedTuple

The rows of a pasted block: `(; rows, skipped, separator, decimal_comma)`.
Each line gives a formula, an amount and, optionally, a unit; a line whose
amount is not a number is skipped (a header, a total) and reported in
`skipped`; the unit of a row without one is `unit`. A percent sign stuck to
the amount (`12.5%`) is its unit.
"""
function parse_table(text::AbstractString; unit::AbstractString = "%")
    unit in ALL_UNITS || throw(ArgumentError("unknown unit $(repr(unit))"))
    lines = filter(!isempty, String.(strip.(split(replace(text, "\r\n" => "\n", '\r' => '\n'), '\n'))))
    separator = detect_separator(lines)
    decimal_comma = separator !== ','
    rows, skipped = AmountRow[], String[]
    for line in lines
        cells = filter(!isempty, split_cells(line, separator))
        if length(cells) < 2
            push!(skipped, line)
            continue
        end
        formula, raw = cells[1], cells[2]
        row_unit = length(cells) >= 3 ? canonical_unit(cells[3]) : nothing
        if endswith(raw, '%')
            raw, row_unit = raw[1:prevind(raw, lastindex(raw))], "%"
        end
        amount = parse_number(raw; decimal_comma)
        if amount === nothing
            push!(skipped, line)
            continue
        end
        push!(rows, AmountRow(formula, amount, something(row_unit, unit), "auto"))
    end
    return (; rows, skipped, separator = separator === nothing ? "blanks" : string(separator), decimal_comma)
end

"""
    table_text(rows) -> String

Rows as tab-separated text, the way a spreadsheet takes a paste.
"""
table_text(rows) = join(("$(r.formula)\t$(r.amount)\t$(r.unit)" for r in rows), '\n')
