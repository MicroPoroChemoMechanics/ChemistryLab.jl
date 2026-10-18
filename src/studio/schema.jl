# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2026 Jean-François Barthélémy (Cerema, UMR MCD)

# One declaration per choice. A choice the page offers (a law of the degree of
# reaction, an activity model, a process, a rate law…) is a Julia type that
# declares its fields once: their kind, default, unit, bounds and level. From
# that declaration come the reading and writing of the model's JSON, the bounds
# a model is checked against, the form the page draws and the keys of its
# labels. A type also writes its own Julia code, by a method of `code`. Adding a
# choice is adding a type, its fields and its code; nothing else is touched.

"""
    Field(name, kind; default, unit, min, max, options, source, level)

A field of a choice: its `name`, as in the JSON and in the struct; its `kind`
(`:number`, `:integer`, `:text`, `:bool`, `:choice`, `:numbers`, `:texts`,
`:rows`, `:nested`, `:keyed`); its `default`; the `unit` shown beside it; the bounds `min`
and `max` of a number; the fixed `options` of a `:choice`, or the catalog list
it takes them from (`source`), or the family of a `:nested` choice or of the
values of a `:keyed` one (a choice per key: a law per constituent); and the
`level`, `:novice` or `:expert` (hidden in the guided mode).
"""
struct Field
    name::String
    kind::Symbol
    default::Any
    unit::String
    min::Union{Nothing, Float64}
    max::Union{Nothing, Float64}
    options::Vector{String}
    source::String
    level::Symbol
end

const FIELD_KINDS = (:number, :integer, :text, :bool, :choice, :numbers, :texts, :rows, :nested, :keyed)

function Field(
        name::AbstractString, kind::Symbol; default = nothing, unit::AbstractString = "",
        min::Union{Nothing, Real} = nothing, max::Union{Nothing, Real} = nothing,
        options = String[], source::AbstractString = "", level::Symbol = :novice,
    )
    kind in FIELD_KINDS || throw(ArgumentError("unknown field kind $kind"))
    level in (:novice, :expert) || throw(ArgumentError("unknown field level $level"))
    return Field(
        String(name), kind, default, String(unit), min === nothing ? nothing : Float64(min),
        max === nothing ? nothing : Float64(max), String.(options), String(source), level,
    )
end

to_dict(f::Field) = Dict{String, Any}(
    "name" => f.name, "kind" => String(f.kind), "default" => plain(f.default), "unit" => f.unit,
    "min" => f.min, "max" => f.max, "options" => f.options, "source" => f.source, "level" => String(f.level),
)

"""
    Choice

A kind of thing the page lets the user choose, with fields of its own. Each
family (`ExtentChoice`, `ActivityChoice`…) is an abstract subtype listing its
members with `members`; each member gives its `kind` (the name in the JSON),
its `fields`, in the order of its struct, and its `code`.
"""
abstract type Choice end

"The concrete choices of a family, in the order the page lists them."
function members end

"The name of a choice in the JSON and in the labels of the page."
function kind end

"The fields of a choice, in the order of its struct."
function fields end

"""
    plain(x)

A value as JSON carries it: a choice as its object, a vector element by
element, anything else as it is.
"""
plain(x) = x
plain(x::AbstractVector) = Any[plain(v) for v in x]
plain(x::Choice) = to_dict(x)
plain(p::Pair{String, <:Choice}) = Dict{String, Any}("key" => first(p), "value" => to_dict(last(p)))

# Two choices are equal when their fields are: the default `==` of a struct
# compares its vectors by identity.
Base.:(==)(a::T, b::T) where {T <: Choice} = all(getfield(a, f) == getfield(b, f) for f in fieldnames(T))
Base.hash(x::T, h::UInt) where {T <: Choice} = foldr(hash, (getfield(x, f) for f in fieldnames(T)); init = hash(T, h))

"The JSON object of a choice: its kind, and its fields."
function to_dict(x::T) where {T <: Choice}
    d = Dict{String, Any}("kind" => kind(T))
    for f in fields(T)
        d[f.name] = plain(getfield(x, Symbol(f.name)))
    end
    return d
end

# The value of one field, from what the JSON holds, with its type checked: the
# page is ours, so a mismatch is a bug to name, not a value to guess.
function field_value(f::Field, v)
    v === nothing && f.default === nothing && f.kind === :nested && return nothing
    if f.kind === :number
        (v isa Real && !(v isa Bool)) || throw(ArgumentError("field \"$(f.name)\" must be a number, got $(repr(v))"))
        return Float64(v)
    elseif f.kind === :integer
        (v isa Real && !(v isa Bool) && isinteger(v)) ||
            throw(ArgumentError("field \"$(f.name)\" must be an integer, got $(repr(v))"))
        return Int(v)
    elseif f.kind in (:text, :choice)
        v isa AbstractString || throw(ArgumentError("field \"$(f.name)\" must be text, got $(repr(v))"))
        (f.kind === :choice && !isempty(f.options) && v ∉ f.options) &&
            throw(ArgumentError("field \"$(f.name)\" must be one of $(join(f.options, ", ")), got $(repr(v))"))
        return String(v)
    elseif f.kind === :bool
        v isa Bool || throw(ArgumentError("field \"$(f.name)\" must be true or false, got $(repr(v))"))
        return v
    elseif f.kind === :numbers
        (v isa AbstractVector && all(x -> x isa Real && !(x isa Bool), v)) ||
            throw(ArgumentError("field \"$(f.name)\" must be a list of numbers, got $(repr(v))"))
        return Float64.(v)
    elseif f.kind === :texts
        (v isa AbstractVector && all(x -> x isa AbstractString, v)) ||
            throw(ArgumentError("field \"$(f.name)\" must be a list of texts, got $(repr(v))"))
        return String.(v)
    elseif f.kind === :rows
        v isa AbstractVector || throw(ArgumentError("field \"$(f.name)\" must be a list of rows, got $(repr(v))"))
        return AmountRow[amount_row(r) for r in v]
    end
    F = FAMILIES[f.source]
    f.kind === :nested && return choice_from_dict(F, v)
    # :keyed — a choice of the family `source` for each key
    v isa AbstractVector || throw(ArgumentError("field \"$(f.name)\" must be a list, got $(repr(v))"))
    return Pair{String, F}[keyed_entry(f, F, e) for e in v]
end

function keyed_entry(f::Field, ::Type{F}, e) where {F <: Choice}
    (e isa AbstractDict && get(e, "key", nothing) isa AbstractString) ||
        throw(ArgumentError("each entry of \"$(f.name)\" must be an object with a key, got $(repr(e))"))
    return String(e["key"]) => choice_from_dict(F, get(e, "value", nothing))
end

"""
    build(T, d) -> T

The choice of type `T` an object describes, each field read with its type
checked, the declared default filling what the object omits.
"""
function build(::Type{T}, d::AbstractDict) where {T <: Choice}
    values = map(fields(T)) do f
        field_value(f, get(d, f.name, f.default))
    end
    return T(values...)
end

"""
    choice_from_dict(F, d) -> F

The member of family `F` an object names by its `kind`.
"""
function choice_from_dict(::Type{F}, d) where {F <: Choice}
    d isa AbstractDict || throw(ArgumentError("a $(family_name(F)) must be an object, got $(repr(d))"))
    k = get(d, "kind", nothing)
    for T in members(F)
        kind(T) == k && return build(T, d)
    end
    throw(ArgumentError("unknown $(family_name(F)) $(repr(k)); the known ones are $(join(kind.(members(F)), ", "))"))
end

"The families of choices, by the name their fields and the page use."
const FAMILIES = Dict{String, Type{<:Choice}}()

family_name(::Type{F}) where {F <: Choice} = only(k for (k, v) in FAMILIES if v === F)

"""
    register_family!(name, F)

Make family `F` known under `name`, to nested fields and to the schema route.
"""
register_family!(name::AbstractString, ::Type{F}) where {F <: Choice} = (FAMILIES[String(name)] = F; nothing)

"What the page knows of a member besides its fields; a family adds its own."
member_traits(::Type{<:Choice}) = Dict{String, Any}()

"""
    schema() -> Dict{String, Any}

Every family with its members and their fields, for the page to draw its forms.
"""
schema() = Dict{String, Any}(
    name => [merge(Dict{String, Any}("kind" => kind(T), "fields" => to_dict.(fields(T))), member_traits(T)) for T in members(F)]
        for (name, F) in FAMILIES
)

"""
    bound_problems(x, where_) -> Vector{Problem}

A problem for each number of choice `x` outside the bounds its field declares,
`where_` saying in which part of the model the choice sits.
"""
function bound_problems(x::T, where_::AbstractString) where {T <: Choice}
    out = Problem[]
    for f in fields(T)
        v = getfield(x, Symbol(f.name))
        if f.kind === :number
            out_of_bounds(f, v) && push!(out, range_problem(f, v, where_))
        elseif f.kind === :numbers
            any(w -> out_of_bounds(f, w), v) && push!(out, range_problem(f, first(w for w in v if out_of_bounds(f, w)), where_))
        elseif f.kind === :nested && v !== nothing
            append!(out, bound_problems(v, where_))
        elseif f.kind === :keyed
            for (k, c) in v
                append!(out, bound_problems(c, "$where_, $k"))
            end
        end
    end
    return out
end

out_of_bounds(f::Field, v::Real) = (f.min !== nothing && v < f.min) || (f.max !== nothing && v > f.max)

function range_problem(f::Field, v, where_)
    lo = f.min === nothing ? "" : "from $(f.min) "
    hi = f.max === nothing ? "" : "up to $(f.max)"
    return Problem(
        "out_of_range", "$(f.name) = $v in $where_ is outside its range ($(strip(lo * hi))).",
        Dict{String, Any}("field" => f.name, "value" => v, "where" => where_, "min" => f.min, "max" => f.max),
    )
end
