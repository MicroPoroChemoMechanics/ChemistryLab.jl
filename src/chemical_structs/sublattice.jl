# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

using LinearAlgebra: rank

# ── Ideal mixing on sublattices ───────────────────────────────────────────────

"""
    struct SublatticeModel{T<:Real} <: AbstractSolidSolutionModel

Ideal mixing on several distinct sites (sublattices) of one formula unit: the
model of a solid solution whose end-members differ by what occupies a few
structural positions, such as the bridging tetrahedra and interlayer sites of a
C-S-H.

# The model

A formula unit has sites `s = 1 … S`, site `s` counted `mₛ` times in the
formula. Every end-member `k` puts exactly one species `σₛ(k)` on each site. The
fraction of site `s` held by species `i` is the sum of the mole fractions of
the end-members that put `i` there,

```math
y_{s,i} = \\sum_k [\\sigma_s(k) = i]\\, x_k ,
```

and the species mix at random on each site, independently of the other sites.
The configurational Gibbs energy per mole of formula units is therefore

```math
\\frac{G_\\text{mix}}{RT} = \\sum_s m_s \\sum_i y_{s,i} \\ln y_{s,i} ,
```

and the activity of an end-member is the product of its own site fractions:

```math
\\ln a_k = \\sum_s m_s \\ln y_{s,\\sigma_s(k)} .
```

That is exactly ``\\partial (n\\,G_\\text{mix}/RT)/\\partial n_k``: with
``N_{s,i} = n\\,y_{s,i}`` the amounts on each site, the derivative of
``\\sum_s m_s \\sum_i N_{s,i} \\ln(N_{s,i}/n)`` is
``\\sum_s m_s [\\ln y_{s,\\sigma_s(k)} + 1 - \\sum_i N_{s,i}/n]``, and the bracket
loses its last two terms because every end-member fills every site once. The
log activities are the gradient of one energy, so a certificate keeps its
`:global_minimum` scope.

# Two properties the solver has to know

**Convexity.** ``G_\\text{mix}`` is a sum of convex functions of the site
fractions, which are linear in ``x``, so it is convex. It is *strictly* convex
exactly when the occupancy matrix ``O`` — one row per (site, species), one
column per end-member, ``O_{(s,i),k} = [\\sigma_s(k) = i]`` — has full column
rank: a null vector of ``O`` is a direction along which no site fraction
changes. `rank` stores that rank.

**The dilute limit.** When ``x_k \\to 0``, the site fractions of the species
that ``k`` alone carries vanish with it, and the others do not, so

```math
\\ln a_k \\simeq e_k \\ln x_k + O(1), \\qquad
e_k = \\sum_{s \\,:\\, \\sigma_s(k) \\text{ unique to } k} m_s .
```

An end-member with ``e_k = 0`` has no species of its own: its activity stays
finite as it disappears, so it can be **exactly absent from a phase that is
present** — its stationarity is then an inequality, as a pure phase's is.
`exponents` stores the ``e_k``.

# Fields

  - `sites`: the name of each site;
  - `multiplicity`: ``m_s``, the coefficient of each site in the formula unit
    of the end-member records the model is used with;
  - `species`: for each site, the names of the species found there;
  - `occupancy`: `occupancy[s, k]` is the index in `species[s]` of the species
    end-member `k` puts on site `s`;
  - `exponents`: ``e_k``;
  - `rank`: the rank of ``O``.

A site that holds one species in every end-member contributes nothing and may
be left out.

# Construction

```julia
SublatticeModel(multiplicity, occupancy; sites)
```

with `occupancy` a matrix of species names, one row per site and one column per
end-member, in the order of the end-members of the phase. Two end-members with
the same occupancy on every site are refused: the model cannot tell them apart.
[`sublattice_model`](@ref) builds the model of a published phase from
`data/literature`.

# Example

The CSH3T model of Kulik (2011), Eq. (19), on its own formula unit: two sites
holding Si or Ca, and the pentameric end-member T5C carrying Ca on the first
and Si on the second.

```jldoctest
julia> m = SublatticeModel([1, 1], ["Si" "Ca" "Ca"; "Si" "Si" "Ca"]; sites = ["BTI1", "BTI2"]);

julia> m.exponents   # TobH and T2C each own a species on one site; T5C owns none
3-element Vector{Int64}:
 1
 0
 1

julia> m.rank
3
```

See also: [`site_fractions`](@ref), [`sublattice_model`](@ref),
[`IdealSolidSolutionModel`](@ref).

# References

  - [Kulik2011](@cite), Section 4.2 (CSH3T).
  - [Myers2014](@cite), Eqs. (15) and (18) (CNASH_ss).
"""
struct SublatticeModel{T <: Real} <: AbstractSolidSolutionModel
    sites::Vector{String}
    multiplicity::Vector{T}
    species::Vector{Vector{String}}
    occupancy::Matrix{Int}
    exponents::Vector{T}
    rank::Int

    # Inner, so that no generated constructor can bypass the checks and the
    # two derived fields.
    function SublatticeModel{T}(
            sites::AbstractVector{<:AbstractString}, multiplicity::AbstractVector,
            species::AbstractVector, occupancy::AbstractMatrix{<:Integer},
        ) where {T <: Real}
        S, K = size(occupancy)
        length(multiplicity) == S == length(sites) == length(species) || throw(
            ArgumentError(
                "SublatticeModel: $(length(sites)) site names, $(length(multiplicity)) " *
                    "multiplicities and $(length(species)) species lists for an occupancy " *
                    "of $S sites; they must agree."
            )
        )
        K >= 1 || throw(ArgumentError("SublatticeModel: no end-member."))
        for s in 1:S
            multiplicity[s] > 0 || throw(
                ArgumentError(
                    "SublatticeModel: site $(sites[s]) has multiplicity $(multiplicity[s]); " *
                        "a site counted zero or fewer times holds nothing."
                )
            )
            for k in 1:K
                1 <= occupancy[s, k] <= length(species[s]) || throw(
                    ArgumentError(
                        "SublatticeModel: end-member $k points at species " *
                            "$(occupancy[s, k]) of site $(sites[s]), which has " *
                            "$(length(species[s]))."
                    )
                )
            end
        end
        for k in 1:K, j in (k + 1):K
            view(occupancy, :, k) == view(occupancy, :, j) && throw(
                ArgumentError(
                    "SublatticeModel: end-members $k and $j put the same species on " *
                        "every site, so the model gives them the same activity and " *
                        "cannot tell them apart. Check the occupancy."
                )
            )
        end
        m = collect(T, multiplicity)
        occ = Matrix{Int}(occupancy)
        # e_k: the multiplicity of the sites on which k alone carries its species.
        e = [
            sum((m[s] for s in 1:S if count(==(occ[s, k]), view(occ, s, :)) == 1); init = zero(T))
                for k in 1:K
        ]
        O = zeros(sum(length, species), K)
        row = 0
        for s in 1:S
            for k in 1:K
                O[row + occ[s, k], k] = 1.0
            end
            row += length(species[s])
        end
        return new{T}(
            collect(String, sites), m, [collect(String, sp) for sp in species], occ, e, rank(O),
        )
    end
end

"""
    SublatticeModel(multiplicity, occupancy; sites) -> SublatticeModel

Build the model from the coefficient of each site and a matrix of species
names, `occupancy[s, k]` being what end-member `k` puts on site `s`. The species
of each site are numbered in the order they first appear along the row.
"""
function SublatticeModel(
        multiplicity::AbstractVector{<:Real}, occupancy::AbstractMatrix;
        sites::AbstractVector{<:AbstractString} = ["site $s" for s in 1:size(occupancy, 1)],
    )
    S, K = size(occupancy)
    labels = map(string, occupancy)
    species = [unique(view(labels, s, :)) for s in 1:S]
    idx = [findfirst(==(labels[s, k]), species[s]) for s in 1:S, k in 1:K]
    T = mapreduce(typeof, promote_type, multiplicity)
    return SublatticeModel{T}(sites, multiplicity, species, idx)
end

_n_members(m::SublatticeModel) = size(m.occupancy, 2)

"""
    site_fractions(model::SublatticeModel, x) -> Vector{Vector}

The site fractions ``y_{s,i}`` at the mole fractions `x` of the end-members:
`site_fractions(model, x)[s][i]` is the fraction of site `s` held by
`model.species[s][i]`. Each inner vector sums to `sum(x)`. The element type is
that of `x`, so the function carries dual and symbolic numbers.
"""
function site_fractions(m::SublatticeModel, x::AbstractVector)
    length(x) == _n_members(m) || throw(
        DimensionMismatch("site_fractions: $(length(x)) mole fractions for $(_n_members(m)) end-members")
    )
    y = [zeros(eltype(x), length(sp)) for sp in m.species]
    @inbounds for s in eachindex(y), k in eachindex(x)
        y[s][m.occupancy[s, k]] += x[k]
    end
    return y
end

function Base.show(io::IO, m::SublatticeModel)
    return print(
        io, "SublatticeModel(", length(m.sites), " sites, ", _n_members(m),
        " end-members, rank ", m.rank, ")",
    )
end

# ── The published models, from data/literature ───────────────────────────────

"""
    sublattice_model(ref, end_members) -> SublatticeModel

The sublattice model a published source gives for a phase, read from
`data/literature/<key>.json`, with its end-members in the order of
`end_members`.

`ref` is `"<key>:<model>"`, such as `"Myers2014:cnash"` or `"Kulik2011:csh3t"`;
the file holds the tables `<model>_sites` (the coefficient of each site in the
source's formula unit) and `<model>_occupancy` (the species each end-member puts
on each site, one row per end-member and site, with the database name of the
end-member in its first column).

# Formula units

A database may store an end-member on another formula unit than the source's:
the Cemdata18 records of CSH3T are half the units of Kulik (2011). The site
coefficients belong to the source's unit, and a model applied to records `f`
times as large has coefficients `f` times as large, because the configurational
entropy is extensive. When `end_members` are species, their formulas are
compared with the table `<model>_end_members`, the formula the source prints
for each end-member: every record must be the same multiple `f` of its
published formula, and the coefficients are multiplied by `f`. A record that is
not such a multiple is refused, since the model then describes another
substance. When `end_members` are names, `f = 1` and no check is made.

# Example

```julia
subs = Dict(symbol(s) => s for s in build_species(datapath("cemdata18-thermofun.json")))
m = sublattice_model("Kulik2011:csh3t", [subs[n] for n in ("CSH3T-TobH", "CSH3T-T5C", "CSH3T-T2C")])
m.multiplicity   # [0.5, 0.5]: the records are half the published unit
```
"""
function sublattice_model(ref::AbstractString, end_members::AbstractVector{<:AbstractString}; scale::Real = 1)
    key, prefix = _split_sublattice_ref(ref)
    sites = literature_table(key, "$(prefix)_sites")
    occ = literature_table(key, "$(prefix)_occupancy")
    site_names = String.(sites.site)
    m = Float64.(ustrip.(sites.multiplicity)) .* scale
    labels = Matrix{String}(undef, length(site_names), length(end_members))
    filled = falses(size(labels))
    for (em, site, sp) in zip(occ.end_member, occ.site, occ.species)
        k = findfirst(==(em), end_members)
        k === nothing && continue
        s = findfirst(==(site), site_names)
        s === nothing && error("sublattice_model: $ref places $em on site $site, which $(prefix)_sites does not list.")
        labels[s, k] = sp
        filled[s, k] = true
    end
    for k in eachindex(end_members), s in eachindex(site_names)
        filled[s, k] || error(
            "sublattice_model: $ref gives no species for $(end_members[k]) on site " *
                "$(site_names[s]); the end-members it describes are " *
                "$(join(unique(occ.end_member), ", "))."
        )
    end
    return SublatticeModel(m, labels; sites = site_names)
end

function sublattice_model(ref::AbstractString, end_members::AbstractVector{<:AbstractSpecies})
    key, prefix = _split_sublattice_ref(ref)
    names = [symbol(sp) for sp in end_members]
    published = literature_table(key, "$(prefix)_end_members")
    f = nothing
    for sp in end_members
        i = findfirst(==(symbol(sp)), published.end_member)
        i === nothing && error(
            "sublattice_model: $ref prints no formula for $(symbol(sp)); it describes " *
                "$(join(published.end_member, ", "))."
        )
        r = _formula_ratio(atoms(sp), composition(Formula(published.formula[i])))
        r === nothing && error(
            "sublattice_model: $(symbol(sp)) has the formula $(formula(sp)), which is " *
                "not a multiple of the formula $ref prints for it, $(published.formula[i]); " *
                "the model would describe another substance."
        )
        if f === nothing
            f = r
        elseif !isapprox(r, f; rtol = 1.0e-9)
            error(
                "sublattice_model: the records are not one multiple of the published " *
                    "formula units: $(symbol(sp)) is $(r) times its unit, the first " *
                    "end-member $(f) times. Site coefficients cannot be rescaled for both."
            )
        end
    end
    return sublattice_model(ref, names; scale = f)
end

function _split_sublattice_ref(ref::AbstractString)
    parts = split(ref, ":"; limit = 2)
    length(parts) == 2 || error(
        "sublattice_model: expected \"<literature key>:<model>\", such as " *
            "\"Myers2014:cnash\"; got \"$ref\"."
    )
    return String(parts[1]), String(parts[2])
end

# `a = f b` element by element, or `nothing` when no single `f` does.
function _formula_ratio(a::AbstractDict, b::AbstractDict)
    ka = Set(k for (k, v) in a if !iszero(v))
    kb = Set(k for (k, v) in b if !iszero(v))
    ka == kb || return nothing
    isempty(ka) && return nothing
    f = nothing
    for k in ka
        r = Float64(a[k]) / Float64(b[k])
        if f === nothing
            f = r
        elseif !isapprox(r, f; rtol = 1.0e-9)
            return nothing
        end
    end
    return f
end
