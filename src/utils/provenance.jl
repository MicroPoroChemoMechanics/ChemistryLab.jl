# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

# ── A number that says where it came from ────────────────────────────────────
#
# This repository already makes the distinction in data. `data/pitzer-reardon1990.toml`
# carries `origin = "estimated:<analog>"` on every coefficient it had to borrow,
# and says why in its own header: *"a caller reading beta0 for Al(OH)4- would
# otherwise have no way of knowing that no aluminate solution was ever measured
# to obtain it"*. What stops at the file should not: a value read out of it,
# combined with others and printed in a table has lost the one thing that said
# how much to trust it.

"""
    ProvenanceKind

How a number came to be, ordered from the weakest claim to the strongest.

  - `PROV_PLACEHOLDER` — a value standing in for one nobody has yet. It exists
    so a calculation can be *demonstrated*, and it is not evidence of anything.
  - `PROV_ESTIMATED` — borrowed from an analog or inferred from a correlation.
    Somebody reasoned, nobody measured.
  - `PROV_FITTED` — identified from data by an optimization. As sound as the
    data, the model and the identifiability of the parameter, which is three
    separate questions.
  - `PROV_PUBLISHED` — taken from a named source that reports it as a value.
  - `PROV_MEASURED` — determined experimentally in the work at hand.
  - `PROV_UNSTATED` — nothing is claimed. Deliberately the **weakest**, so a
    number that forgot to say where it came from never strengthens a result.

The order is what [`weakest`](@ref) uses, and the ordering is the point: a
result is only as sound as its least sound input.
"""
@enum ProvenanceKind begin
    PROV_UNSTATED
    PROV_PLACEHOLDER
    PROV_ESTIMATED
    PROV_FITTED
    PROV_PUBLISHED
    PROV_MEASURED
end

"""
    struct Traced{T}

A value carrying how it was obtained and from where.

`Traced` is **not** a `Real`. That is deliberate: a number that flowed silently
into arithmetic would arrive at the far end with its provenance gone, and the
whole point is that losing it should be an act rather than an accident. Unwrap
it with [`value`](@ref) at the boundary where it enters a calculation, which is
where a reader can see the claim being dropped.

# Fields

  - `value`: the number.
  - `kind`: a [`ProvenanceKind`](@ref).
  - `source`: free text — a database and version, a paper, the analog a value
    was borrowed from, the dataset it was fitted to.
  - `uncertainty`: what the source says about how well it is known, in the same
    unit as `value`, or `nothing` when it says nothing. Published sorption
    compilations report it per entry — ClaySor 2023 writes `error: 0.17` beside
    a `log K` — and a compilation where half the constants are known to 0.1 and
    half to 0.5 is not one number's worth of information.

# Example

```jldoctest
julia> using ChemistryLab

julia> ε = Traced(0.05, PROV_ESTIMATED, "analog: SO4-2");

julia> ChemistryLab.value(ε), provenance(ε)
(0.05, PROV_ESTIMATED)

julia> is_evidence(ε)
false
```

`value` and `source` are **not exported**, and reaching for them through the
module is not an oversight. Both names are taken: `Symbolics` and `SciMLBase`
each export a `value`, so exporting one here would make the bare name ambiguous
for anyone loading either beside this package — and a doc page that did would
fail on use rather than on import.

See also: [`weakest`](@ref), [`provenance_report`](@ref).
"""
struct Traced{T}
    value::T
    kind::ProvenanceKind
    source::String
    uncertainty::Union{Nothing, T}
end

Traced(v, kind::ProvenanceKind, source::AbstractString; uncertainty = nothing) =
    Traced(v, kind, String(source), uncertainty === nothing ? nothing : convert(typeof(v), uncertainty))
Traced(v) = Traced(v, PROV_UNSTATED, "")
Traced(v, kind::ProvenanceKind) = Traced(v, kind, "")

"""
    value(t::Traced)

The number, with its provenance dropped — explicitly, which is the only way it
should be dropped.
"""
value(t::Traced) = t.value
value(x) = x

"""
    provenance(t::Traced) -> ProvenanceKind

How `t` was obtained.
"""
provenance(t::Traced) = t.kind
provenance(::Any) = PROV_UNSTATED

"""
    source(t::Traced) -> String

Where `t` came from, as free text.
"""
source(t::Traced) = t.source
source(::Any) = ""

"""
    is_evidence(t) -> Bool

Whether `t` is a claim about the world somebody stands behind — `PROV_MEASURED`
or `PROV_PUBLISHED` — as opposed to something estimated, fitted, standing in, or
unstated.

A fitted value is deliberately **not** evidence here. It may be excellent, and
whether it is depends on the data, the model and the parameter's
identifiability, which is a judgement this predicate has no way to make.
"""
is_evidence(t) = provenance(t) in (PROV_MEASURED, PROV_PUBLISHED)

"""
    weakest(ts...) -> ProvenanceKind

The least sound provenance among its arguments — what a result computed from
them can honestly claim.

Nothing here propagates it automatically, because [`Traced`](@ref) is not a
number and does not do arithmetic. This is the function to call where a derived
quantity is assembled, so that the claim is written down rather than inherited
by accident.
"""
weakest(ts...) = isempty(ts) ? PROV_UNSTATED : minimum(provenance(t) for t in ts)

function Base.show(io::IO, t::Traced)
    print(io, t.value)
    t.uncertainty === nothing || print(io, " ± ", t.uncertainty)
    print(io, " [", _prov_label(t.kind))
    isempty(t.source) || print(io, ": ", t.source)
    return print(io, "]")
end

_prov_label(k::ProvenanceKind) = lowercase(replace(string(k), "PROV_" => ""))

"""
    uncertainty(t::Traced)

What the source says about how well `t` is known, or `nothing` when it says
nothing — which is not the same as saying it is exact.
"""
uncertainty(t::Traced) = t.uncertainty
uncertainty(::Any) = nothing

"""
    provenance_report(values) -> NamedTuple

How a collection of [`Traced`](@ref) values divides by provenance: a count per
kind, the weakest present, and whether every one of them is evidence.

Meant to be printed beside a result. A table of coefficients where nine are
published and one is a placeholder is a different object from one where all ten
are published, and the difference does not show in any of the numbers.

`values` may be any iterable of `Traced`, including the values of a `Dict`.
"""
function provenance_report(values)
    items = collect(values)
    counts = Dict(k => count(t -> provenance(t) === k, items) for k in instances(ProvenanceKind))
    return (;
        total = length(items),
        counts = counts,
        weakest = isempty(items) ? PROV_UNSTATED : weakest(items...),
        all_evidence = !isempty(items) && all(is_evidence, items),
        # How many say nothing about how well they are known — which is not the
        # same as saying they are exact.
        without_uncertainty = count(t -> uncertainty(t) === nothing, items),
    )
end
