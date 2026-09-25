# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

using DynamicQuantities
using JSON
using OrderedCollections

# ── Published values, read from data files rather than typed into code ───────
#
# A value taken from an article is data, and it belongs where data is kept: in a
# file that names its source, says where in the source it was read, how it was
# checked and in what unit it is expressed. Typed into a source file it loses
# all of that, and a second copy typed elsewhere drifts from the first without
# anyone noticing. One file per source, `data/literature/<key>.json`, the key
# being the entry of `docs/src/refs.bib`: the bibliography says what the source
# is, the file says what was taken from it.
#
# Coefficients that define an equation of state or a published model (the
# water equation of state, the HKF constants, the Pitzer α and b) are not data
# of this kind: they are the model, and they stay in the code with their source
# in a comment. What is fitted or measured comes from here.

"""
    LITERATURE_SCHEMA

Identifier of the format of the files in `data/literature/`, written in each of
them as `"schema"`. See [`literature`](@ref) for the fields.
"""
const LITERATURE_SCHEMA = "chemistrylab-literature/1"

const _LITERATURE_KINDS = Dict(
    "published" => PROV_PUBLISHED,
    "measured" => PROV_MEASURED,
    "fitted" => PROV_FITTED,
    "estimated" => PROV_ESTIMATED,
    "placeholder" => PROV_PLACEHOLDER,
    "unstated" => PROV_UNSTATED,
)

"""
    struct LiteratureRecord

The content of one `data/literature/<key>.json` file, as returned by
[`literature`](@ref).

# Fields

  - `key`: the bibliography key, which is also the file name.
  - `path`: the file it was read from.
  - `source`: the `"source"` object of the file (key, DOI, citation, location).
  - `transcription`: how the values were transcribed and checked.
  - `quantities`: named scalars, each a [`Traced`](@ref) value carrying its
    provenance kind, its source and, when the source gives one, its uncertainty.
  - `tables`: named tables, each a `NamedTuple` of column vectors.
  - `notes`: free remarks recorded with the data.

`record[name]` returns the quantity `name`; [`literature_value`](@ref) and
[`literature_table`](@ref) are the shortcuts most callers need.
"""
struct LiteratureRecord
    key::String
    path::String
    source::Dict{String, Any}
    transcription::Dict{String, Any}
    quantities::OrderedDict{String, Traced}
    tables::OrderedDict{String, NamedTuple}
    notes::Vector{String}
end

function Base.getindex(r::LiteratureRecord, name::AbstractString)
    haskey(r.quantities, name) || throw(
        KeyError(
            "$(r.key) has no quantity \"$name\"; it has: " *
                join(keys(r.quantities), ", ")
        )
    )
    return r.quantities[name]
end

function Base.show(io::IO, r::LiteratureRecord)
    return print(
        io, "LiteratureRecord(\"", r.key, "\": ", length(r.quantities),
        " quantities, ", length(r.tables), " tables)"
    )
end

"""
    literature_path(key) -> String

The file that holds the values taken from the source `key`:
`data/literature/<key>.json` in the package directory.
"""
literature_path(key::AbstractString) =
    joinpath(pkgdir(@__MODULE__), "data", "literature", String(key) * ".json")

const _LITERATURE_CACHE = Dict{String, LiteratureRecord}()
const _LITERATURE_LOCK = ReentrantLock()

"""
    literature(key) -> LiteratureRecord

The values taken from the published source `key`, read from
`data/literature/<key>.json` and checked against the format
[`LITERATURE_SCHEMA`](@ref) as they are read. A malformed file is an error that
names the file and the field, never a silently different number.

# The format

```json
{
  "schema": "chemistrylab-literature/1",
  "source": {"key": "Powers1948", "doi": null,
             "citation": "…", "location": "p. 185"},
  "transcription": {"from": "…", "checked_against_source": false,
                    "note": "…"},
  "quantities": {
    "w_c_sealed": {"value": 0.42, "unit": "1", "kind": "published",
                   "uncertainty": null, "location": "…",
                   "description": "…"}
  },
  "tables": {
    "name": {"columns": ["m", "phi"], "units": ["mol/kg", "1"],
             "kind": "published", "location": "Table 16",
             "rows": [[0.1, 0.932], [0.2, 0.925]]}
  },
  "notes": ["…"]
}
```

  - `source.key` must equal the file name and be an entry of `refs.bib`, whose
    DOI, when it has one, `source.doi` repeats.
  - `kind` is one of `published`, `measured`, `fitted`, `estimated`,
    `placeholder`, `unstated`, the names of [`ProvenanceKind`](@ref).
  - `unit` is unit arithmetic read by `DynamicQuantities` (`"J/mol"`,
    `"m^2/g"`); `"1"` returns a plain number. A column whose unit is `null`
    holds text (a phase name, a label) and is returned as such.

Records are read once and cached. Reading one while the package precompiles
registers the file as a dependency, so that editing it recompiles the package
instead of leaving a stale value in the compiled image.

See also: [`literature_value`](@ref), [`literature_table`](@ref),
[`literature_row`](@ref), [`Traced`](@ref).
"""
function literature(key::AbstractString)
    k = String(key)
    return lock(_LITERATURE_LOCK) do
        get!(_LITERATURE_CACHE, k) do
            path = literature_path(k)
            isfile(path) || throw(
                ArgumentError(
                    "no published data for \"$k\": expected the file $path. " *
                        "Available: " * join(available_literature(), ", ")
                )
            )
            include_dependency(path)
            read_literature(path)
        end
    end
end

"""
    available_literature() -> Vector{String}

The keys of every source with a file in `data/literature/`.
"""
function available_literature()
    dir = joinpath(pkgdir(@__MODULE__), "data", "literature")
    isdir(dir) || return String[]
    return sort!([first(splitext(f)) for f in readdir(dir) if endswith(f, ".json")])
end

"""
    literature_value(key, name)

The value of the quantity `name` taken from the source `key`, with its unit
(a `DynamicQuantities` quantity) or as a plain number when it is dimensionless.
Its provenance is dropped here, which is where a value enters a calculation;
`literature(key)[name]` keeps it.
"""
literature_value(key::AbstractString, name::AbstractString) =
    value(literature(key)[name])

"""
    literature_table(key, name) -> NamedTuple

The table `name` taken from the source `key`, as a `NamedTuple` of columns.
Each column carries the unit the file declares for it, and a column declared
without a unit holds text.
"""
function literature_table(key::AbstractString, name::AbstractString)
    r = literature(key)
    haskey(r.tables, name) || throw(
        KeyError(
            "$(r.key) has no table \"$name\"; it has: " * join(keys(r.tables), ", ")
        )
    )
    return r.tables[name]
end

"""
    literature_row(key, table, label) -> NamedTuple

The row of the table `table` of the source `key` whose first column is `label`,
as a `NamedTuple` of its entries: the phase, mix or species a row is about is
what a caller knows, not its position.

```julia
co = literature_row("BaroghelBouny1999", "retention_fit", "CO")
VanGenuchten(; a = co.a, m = 1 / co.b)
```
"""
function literature_row(key::AbstractString, table::AbstractString, label::AbstractString)
    t = literature_table(key, table)
    first_col = first(t)
    i = findfirst(==(label), first_col)
    i === nothing && throw(
        KeyError(
            "$key, table \"$table\", has no row \"$label\"; its rows are: " *
                join(first_col, ", ")
        )
    )
    return map(col -> col[i], t)
end

# ── reading and checking one file ────────────────────────────────────────────

_literature_error(path, msg) = throw(ArgumentError("$(basename(path)): $msg"))

function _literature_field(d::AbstractDict, name, path, where)
    haskey(d, name) || _literature_error(path, "$where has no \"$name\"")
    return d[name]
end

# Unit arithmetic only, evaluated by walking the expression rather than by
# `uparse`: the parser of `uparse` evaluates into a module of its own, which is
# refused while this package precompiles (the constants of the rate laws are read
# at load time), and it can evaluate a call, which a data file must not be able
# to. Names resolve in the unit registry of `DynamicQuantities` and nowhere else.
const _UNIT_OPS = Dict{Symbol, Function}(:* => *, :/ => /, :^ => ^, :+ => +, :- => -)

# The registry of DynamicQuantities carries only the prefixes each unit commonly
# takes (`kPa` but not `MPa`), and extending it is a change of global state that
# every package in the session would see. A name the registry lacks is read here
# as an SI prefix followed by a unit it has, the factors being the exact powers
# of ten the SI defines.
const _SI_PREFIXES = Dict{Char, Float64}(
    'T' => 1.0e12, 'G' => 1.0e9, 'M' => 1.0e6, 'k' => 1.0e3, 'h' => 1.0e2,
    'd' => 1.0e-1, 'c' => 1.0e-2, 'm' => 1.0e-3, 'μ' => 1.0e-6, 'n' => 1.0e-9,
    'p' => 1.0e-12,
)

_registry_unit(s::Symbol) = isdefined(DynamicQuantities.Units, s) ?
    getfield(DynamicQuantities.Units, s) : nothing

_eval_unit(x::Real, path, where) = x
function _eval_unit(s::Symbol, path, where)
    x = _registry_unit(s)
    if x === nothing
        name = String(s)
        base = _registry_unit(Symbol(chop(name; head = 1, tail = 0)))
        (haskey(_SI_PREFIXES, first(name)) && base isa AbstractQuantity) ||
            _literature_error(path, "$where: unit \"$s\" is unknown")
        x = _SI_PREFIXES[first(name)] * base
    end
    x isa AbstractQuantity || _literature_error(path, "$where: \"$s\" is not a unit")
    return x
end
function _eval_unit(ex::Expr, path, where)
    (ex.head === :call && haskey(_UNIT_OPS, ex.args[1])) ||
        _literature_error(path, "$where: \"$ex\" is not unit arithmetic")
    return _UNIT_OPS[ex.args[1]]((_eval_unit(a, path, where) for a in ex.args[2:end])...)
end
_eval_unit(ex, path, where) = _literature_error(path, "$where: \"$ex\" is not unit arithmetic")

function _literature_unit(u, path, where)
    u == "1" && return 1
    u isa AbstractString || _literature_error(path, "$where: the unit must be a string")
    ex = try
        Meta.parse(u)
    catch
        _literature_error(path, "$where: \"$u\" is not unit arithmetic")
    end
    q = _eval_unit(ex, path, where)
    q isa AbstractQuantity || _literature_error(path, "$where: \"$u\" is not unit arithmetic")
    return q
end

function _literature_kind(k, path, where)
    haskey(_LITERATURE_KINDS, k) || _literature_error(
        path, "$where: kind \"$k\" is not one of " * join(sort!(collect(keys(_LITERATURE_KINDS))), ", ")
    )
    return _LITERATURE_KINDS[k]
end

_literature_number(x, path, where) =
    x isa Real ? Float64(x) : _literature_error(path, "$where: expected a number, got $(repr(x))")

_with_unit(x::Float64, ::Integer) = x
_with_unit(x::Float64, u) = x * u

"""
    read_literature(path) -> LiteratureRecord

Read and check one file of the format [`LITERATURE_SCHEMA`](@ref). Most callers
want [`literature`](@ref), which finds the file from the source key and caches
the result; this reads an arbitrary path and caches nothing.
"""
function read_literature(path::AbstractString)
    d = JSON.parsefile(path; dicttype = Dict{String, Any})
    get(d, "schema", nothing) == LITERATURE_SCHEMA ||
        _literature_error(path, "\"schema\" must be \"$LITERATURE_SCHEMA\"")

    source = _literature_field(d, "source", path, "the file")
    key = _literature_field(source, "key", path, "\"source\"")
    key == first(splitext(basename(path))) ||
        _literature_error(path, "\"source.key\" is \"$key\" but the file is named after another key")
    transcription = _literature_field(d, "transcription", path, "the file")
    label = string(key, (get(source, "location", nothing) === nothing ? "" : ", " * string(source["location"])))

    quantities = OrderedDict{String, Traced}()
    for (name, q) in sort!(collect(get(d, "quantities", Dict{String, Any}())); by = first)
        where = "quantity \"$name\""
        v = _literature_number(_literature_field(q, "value", path, where), path, where)
        u = _literature_unit(_literature_field(q, "unit", path, where), path, where)
        kind = _literature_kind(_literature_field(q, "kind", path, where), path, where)
        unc = get(q, "uncertainty", nothing)
        unc === nothing || (unc = _with_unit(_literature_number(unc, path, where), u))
        loc = get(q, "location", nothing)
        src = loc === nothing ? label : string(key, ", ", loc)
        quantities[name] = Traced(_with_unit(v, u), kind, src; uncertainty = unc)
    end

    tables = OrderedDict{String, NamedTuple}()
    for (name, t) in sort!(collect(get(d, "tables", Dict{String, Any}())); by = first)
        where = "table \"$name\""
        cols = String.(_literature_field(t, "columns", path, where))
        units = _literature_field(t, "units", path, where)
        length(units) == length(cols) ||
            _literature_error(path, "$where: $(length(cols)) columns but $(length(units)) units")
        _literature_kind(_literature_field(t, "kind", path, where), path, where)
        rows = _literature_field(t, "rows", path, where)
        for (i, row) in enumerate(rows)
            length(row) == length(cols) ||
                _literature_error(path, "$where: row $i has $(length(row)) entries for $(length(cols)) columns")
        end
        columns = map(eachindex(cols)) do j
            if units[j] === nothing
                String[string(row[j]) for row in rows]
            else
                u = _literature_unit(units[j], path, "$where, column \"$(cols[j])\"")
                [_with_unit(_literature_number(row[j], path, "$where, column \"$(cols[j])\""), u) for row in rows]
            end
        end
        tables[name] = NamedTuple{Tuple(Symbol.(cols))}(Tuple(columns))
    end

    notes = String[string(n) for n in get(d, "notes", String[])]
    return LiteratureRecord(key, String(path), source, transcription, quantities, tables, notes)
end
