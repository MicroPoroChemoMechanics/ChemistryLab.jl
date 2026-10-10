# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)


# ── Reading a published sorption model ───────────────────────────────────────
#
# A sorption model is a body of fitted constants, and the thing that makes one
# usable by somebody else is not the numbers but what travels with them: which
# measurement each came from, and how well it is known. ClaySor 2023 writes both
# into every line —
#
#     Am+3 + Ilt_sOH = Ilt_sOAm+2 + H+   # ... error: 0.26 ref: Marinich_ea:2024:rep:
#
# — and a reader that kept the `-log_K` and dropped the rest would be throwing
# away the part that says how much to believe it. So this keeps all three, in a
# `Traced`.

"""
    struct SorptionReaction{T}

One reaction of a published sorption model: its equation as written, its
stoichiometry, and its `log K` with [`provenance`](@ref) and
[`uncertainty`](@ref).

# Fields

  - `equation`: the line as the database writes it, kept verbatim so a reader
    can check the parse.
  - `stoichiometry`: `species => coefficient`, negative for reactants and
    positive for products, as a rational number (`0.5 X` is `1//2`).
  - `log_K`: a [`Traced`](@ref) whose source is the `ref:` tag and whose
    uncertainty is the `error:` tag, when the entry carries them, in the number
    type `T` it is given in: a `log K` being fitted carries its derivative.
  - `comment`: the rest of the comment, which usually says what the reaction is
    in words.
"""
struct SorptionReaction{T <: Real}
    equation::String
    stoichiometry::Dict{String, Rational{Int}}
    log_K::Traced{T}
    comment::String
end

# Any string for the equation and the comment, as the struct took before it was
# parametric: the generated constructor of a parametric struct converts nothing.
SorptionReaction(equation::AbstractString, stoichiometry::AbstractDict, log_K::Traced{T}, comment::AbstractString) where {T} =
    SorptionReaction{T}(String(equation), Dict{String, Rational{Int}}(stoichiometry), log_K, String(comment))

"""
    struct SorptionSite

One site of a sorption model — a surface site family or an exchanger — with the
reactions written on it.

  - `master`: the master species, `"Mnt_s"`.
  - `reference`: the reference form, `"Mnt_sOH"` or `"Mntx-"`, which is what a
    site family's free member corresponds to.
  - `reactions`: every reaction whose products name this site.
  - `comment`: what the database calls it, for example `EdgeSite_S_mont`.
"""
struct SorptionSite
    master::String
    reference::String
    reactions::Vector{SorptionReaction}
    comment::String
end

"""
    struct SorptionModel

A published sorption model read from a PHREEQC-format database: its surface
site families, its exchangers, and where it came from.

Reading one is not the same as being able to solve it. A model is written
against a particular **aqueous** thermodynamic database — ClaySor 2023
[ClaySor2023data](@cite) says so in its own first lines, naming PSI/Nagra TDB
2020 — and its constants are that database's constants. Using them over a
different one is a different model, in exactly the way a surface constant fitted
with a diffuse layer is a different constant from one fitted without.

`header` keeps the file's own comment block, which is where such a statement
lives and where the reference list usually is.
"""
struct SorptionModel
    surfaces::Dict{String, SorptionSite}
    exchangers::Dict{String, SorptionSite}
    source::String
    header::String
end

function Base.show(io::IO, m::SorptionModel)
    ns = sum(length(s.reactions) for s in values(m.surfaces); init = 0)
    ne = sum(length(s.reactions) for s in values(m.exchangers); init = 0)
    return print(
        io, "SorptionModel(", length(m.surfaces), " surface sites/", ns,
        " reactions, ", length(m.exchangers), " exchangers/", ne,
        " reactions, ", m.source, ")",
    )
end

"""
    read_sorption_model(path) -> SorptionModel

Read the `SURFACE_MASTER_SPECIES`, `SURFACE_SPECIES`, `EXCHANGE_MASTER_SPECIES`
and `EXCHANGE_SPECIES` blocks of a PHREEQC-format database.

Reading rather than transcribing, for the reason every generator in
`test/reference/` exists: a hand-copied compilation is a transcription, and this
package has already found standard energies that had drifted that way.

**No sorption model ships with this package.** ClaySor 2023 is CC-BY-4.0 and
freely available from its Zenodo deposit [ClaySor2023data](@cite); point this at
your own copy.

# Example

```julia
m = read_sorption_model("claysor23_v0.7.dat")
keys(m.surfaces)                      # "Mnt_s", "Mnt_v", "Mnt_w", "Ilt_s", …
provenance_report(log_constants(m))   # how much of it is known how well
```

See also: [`log_constants`](@ref), [`reactions_involving`](@ref).
"""
function read_sorption_model(path::AbstractString)
    isfile(path) || throw(ArgumentError("no such database: $path"))
    text = read(path, String)
    src = _source_tag(path)

    header = join(
        Iterators.takewhile(
            l -> startswith(strip(l), "#") || isempty(strip(l)),
            split(text, '\n')
        ), '\n',
    )

    surfaces = Dict{String, SorptionSite}()
    exchangers = Dict{String, SorptionSite}()
    # The comments are kept: they hold the source and the uncertainty of a
    # constant (`ref:`, `error:`).
    blocks = _phreeqc_blocks(split(text, '\n'); comments = true)
    _read_masters!(surfaces, blocks, "SURFACE_MASTER_SPECIES")
    _read_masters!(exchangers, blocks, "EXCHANGE_MASTER_SPECIES")
    _read_reactions!(surfaces, blocks, "SURFACE_SPECIES", src)
    _read_reactions!(exchangers, blocks, "EXCHANGE_SPECIES", src)

    (isempty(surfaces) && isempty(exchangers)) && throw(
        ArgumentError(
            "no sorption model in $path: it declares neither " *
                "SURFACE_MASTER_SPECIES nor EXCHANGE_MASTER_SPECIES."
        ),
    )
    return SorptionModel(surfaces, exchangers, src, header)
end

# The lines of the blocks of `keyword`, in order.
_block_lines(blocks, keyword) = (line for b in blocks if b.keyword == keyword for (_, line) in b.lines)

function _read_masters!(into::Dict{String, SorptionSite}, blocks, keyword)
    for raw in _block_lines(blocks, keyword)
        body, comment = _split_comment(raw)
        fields = split(strip(body))
        length(fields) == 2 || continue
        into[fields[1]] = SorptionSite(
            String(fields[1]), _phreeqc_name(fields[2]), SorptionReaction[], comment,
        )
    end
    return into
end

function _read_reactions!(into::Dict{String, SorptionSite}, blocks, keyword, src)
    isempty(into) && return into
    pending = nothing
    pending_comment = ""
    for raw in _block_lines(blocks, keyword)
        body, comment = _split_comment(raw)
        s = strip(body)
        isempty(s) && continue
        # The options in every spelling PHREEQC reads (`-log_k`, `log_k`, `-l`),
        # as in the blocks of species of a database.
        opt = _species_option(first(split(s, r"[\s=]+")))
        if opt === nothing && occursin('=', s)
            pending, pending_comment = s, comment
        elseif opt === :log_k && pending !== nothing
            logk = _parse_log_k(s)
            logk === nothing && (pending = nothing; continue)
            stoich = _parse_sorption_stoichiometry(pending)
            rxn = SorptionReaction(
                pending, stoich,
                Traced(
                    logk, PROV_PUBLISHED, _ref_tag(pending_comment, src);
                    uncertainty = _error_tag(pending_comment),
                ),
                _plain_comment(pending_comment),
            )
            # A reaction belongs to the site whose master species prefixes a
            # species it produces — which is how PHREEQC itself resolves them.
            for (master, site) in into
                any(
                    sp -> startswith(sp, master) && stoich[sp] > 0,
                    keys(stoich),
                ) || continue
                push!(site.reactions, rxn)
                break
            end
            pending = nothing
        end
    end
    return into
end

# The constant of a `log_k` line, whatever the spelling of the option and
# whether an equal sign follows it (`-log_k = 0.0`).
function _parse_log_k(s)
    words = split(s, r"[\s=]+"; keepempty = false)
    return length(words) >= 2 ? tryparse(Float64, words[2]) : nothing
end

_ref_tag(comment, fallback) = let m = match(r"ref:\s*(\S+)", comment)
    m === nothing ? fallback : String(m.captures[1])
end

_error_tag(comment) = let m = match(r"error:\s*([0-9.eE+-]+)", comment)
    m === nothing ? nothing : tryparse(Float64, m.captures[1])
end

_plain_comment(comment) =
    strip(replace(comment, r"error:\s*[0-9.eE+-]+" => "", r"ref:\s*\S+" => ""))

"""
    _parse_sorption_stoichiometry(equation) -> Dict{String,Rational{Int}}

`species => coefficient` for a PHREEQC reaction line, negative on the left of
the `=` and positive on the right, a species on both sides counted once: the
equation is read by [`phreeqc_equation`](@ref), as every PHREEQC reaction is, and
its coefficients as written, exactly (`0.5` is `1//2`).
"""
function _parse_sorption_stoichiometry(equation::AbstractString)
    lhs, rhs = phreeqc_equation(equation)
    out = Dict{String, Rational{Int}}()
    for (side, sgn) in ((lhs, -1), (rhs, 1)), (c, name) in side
        out[name] = get(out, name, 0) + sgn * _exact_coefficient(c)
    end
    return out
end

# A coefficient as written, exactly: `2.0` is `2//1`, `0.5` is `1//2`.
_exact_coefficient(c::Float64) = isinteger(c) ? Int(c) // 1 : rationalize(Int, c)

"""
    log_constants(m::SorptionModel) -> Vector{<:Traced}

Every `log K` of the model, for [`provenance_report`](@ref).
"""
log_constants(m::SorptionModel) = vcat(
    [r.log_K for s in values(m.surfaces) for r in s.reactions],
    [r.log_K for s in values(m.exchangers) for r in s.reactions],
)

"""
    reactions_involving(m::SorptionModel, species) -> Vector{SorptionReaction}

Every reaction of `m` in which `species` appears, on either side.

The way to take a **subset** of a published model: a compilation covering thirty
elements is not something to import whole, and naming which part was taken is
part of saying what was reproduced.
"""
function reactions_involving(m::SorptionModel, species::AbstractString)
    out = SorptionReaction[]
    for site in Iterators.flatten((values(m.surfaces), values(m.exchangers)))
        for r in site.reactions
            haskey(r.stoichiometry, species) && push!(out, r)
        end
    end
    return out
end
