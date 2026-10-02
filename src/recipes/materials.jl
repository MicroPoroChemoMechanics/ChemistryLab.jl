# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

using OrderedCollections

# ── What a material is made of ────────────────────────────────────────────────

"""
    abstract type AbstractConstituent end

One part of a [`Material`](@ref): a database phase ([`MineralConstituent`](@ref))
or a part known only by its oxides ([`OxideConstituent`](@ref)), such as a glass.
Each carries its mass fraction of the material and its own degree of reaction.
"""
abstract type AbstractConstituent end

"""
    MineralConstituent(species; mass_fraction, extent = 1.0, name = symbol(species))

A crystalline phase of a material, given by its database species: alite in a
clinker, calcite in a limestone, quartz in a fly ash. What reacts enters the
equilibrium as that species; what does not keeps its mass, its molar volume and
its enthalpy of formation in the residue. `extent` is a number, a function of
time or an [`AbstractExtent`](@ref).
"""
struct MineralConstituent{S <: AbstractSpecies, E <: AbstractExtent, F <: Real} <: AbstractConstituent
    name::String
    species::S
    mass_fraction::F
    extent::E
end
function MineralConstituent(
        species::AbstractSpecies; mass_fraction::Real, extent = 1.0,
        name::AbstractString = symbol(species),
    )
    0 <= mass_fraction <= 1 || throw(ArgumentError("MineralConstituent $name: a mass fraction is in [0, 1]; got $mass_fraction."))
    aggregate_state(species) == AS_CRYSTAL || throw(
        ArgumentError("MineralConstituent $name: $(symbol(species)) is not a crystalline species.")
    )
    e = _as_extent(extent)
    mf = float(mass_fraction)
    return MineralConstituent{typeof(species), typeof(e), typeof(mf)}(String(name), species, mf, e)
end

"""
    OxideConstituent(name, oxides; mass_fraction = 1.0, extent = 1.0,
                     density = nothing, enthalpy = nothing, source = nothing)

A part of a material known by its oxide analysis only: the glass of a slag or of
a fly ash, the minor oxides of a clinker. `oxides` maps oxide formulas to their
mass fractions **in this constituent** (fractions, not percent), as published:
an analysis that closes a little above or below one is not renormalized, and
what it leaves unaccounted for (a loss on ignition, an undetermined part) is not
counted. What reacts enters the equilibrium through its elements
([`oxide_budget`](@ref)); what does not is residue.

A glass has no formula, so the density (g/cm³) and the standard enthalpy of
formation (J/g) its residue needs are given by the caller, with their `source`,
or left `nothing`: a volume or a heat that needs them is then reported missing,
never estimated.
"""
struct OxideConstituent{E <: AbstractExtent, F <: Real} <: AbstractConstituent
    name::String
    oxides::OrderedDict{String, F}
    mass_fraction::F
    extent::E
    density::Union{Nothing, F}
    enthalpy::Union{Nothing, F}
    source::Union{Nothing, String}
end
function OxideConstituent(
        name::AbstractString, oxides::AbstractDict; mass_fraction::Real = 1.0, extent = 1.0,
        density = nothing, enthalpy = nothing, source = nothing,
    )
    0 <= mass_fraction <= 1 || throw(ArgumentError("OxideConstituent $name: a mass fraction is in [0, 1]; got $mass_fraction."))
    d = density === nothing ? nothing : _in_unit(us"g/cm^3", density)
    h = enthalpy === nothing ? nothing : _in_unit(us"J/g", enthalpy)
    # In the number type of what it is given: an analysis or a density being
    # differentiated is dual.
    F = promote_type(
        typeof(float(mass_fraction)), (typeof(float(v)) for v in values(oxides))...,
        d === nothing ? Float64 : typeof(d), h === nothing ? Float64 : typeof(h),
    )
    ox = OrderedDict{String, F}(String(k) => float(v) for (k, v) in oxides)
    all(>=(0), values(ox)) || throw(ArgumentError("OxideConstituent $name: an oxide fraction is negative."))
    # An X-ray fluorescence analysis closes within a percent or so of 100 %, on
    # either side, and is kept as published; a sum far above one is percent
    # given for fractions.
    sum(values(ox); init = 0.0) <= 1.05 || throw(
        ArgumentError(
            "OxideConstituent $name: the oxide fractions sum to $(sum(values(ox))); " *
                "give fractions of the constituent's mass, not percent (see `literature_oxides`)."
        )
    )
    for k in keys(ox)
        _oxide_cation(k)   # refuses what is not an oxide of one element
    end
    e = _as_extent(extent)
    return OxideConstituent{typeof(e), F}(
        String(name), ox, F(mass_fraction), e, d === nothing ? nothing : F(d), h === nothing ? nothing : F(h),
        source === nothing ? nothing : String(source),
    )
end

# A number in the unit `u`: a bare number is taken as already in it, a quantity
# is converted first (a quantity may carry other symbolic units, m³ for cm³).
_in_unit(_, x::Real) = float(x)
_in_unit(u, x::DynamicQuantities.AbstractQuantity) = ustrip(u, uconvert(u, x))

const _MATERIAL_KINDS = (:cement, :scm, :water, :salt, :other)

"""
    Material(name, kind = :other; constituents, extent = 1.0, source = nothing)

A material of a recipe: a clinker or a cement (`:cement`), a supplementary
cementitious material (`:scm`), water (`:water`), a salt (`:salt`), anything
else (`:other`), made of `constituents`. `extent` is the degree of reaction of
the material as a whole, multiplied by that of each constituent. Water and salts
always react completely: an extent below one on them is refused.

The mass fractions of the constituents may sum to less than one (the rest is not
described and not counted), never to more.
"""
struct Material
    name::String
    kind::Symbol
    constituents::Vector{AbstractConstituent}
    extent::AbstractExtent
    source::Union{Nothing, String}
end
function Material(
        name::AbstractString, kind::Symbol = :other; constituents::AbstractVector,
        extent = 1.0, source = nothing,
    )
    kind in _MATERIAL_KINDS || throw(ArgumentError("Material $name: kind is one of $(_MATERIAL_KINDS); got :$kind."))
    isempty(constituents) && throw(ArgumentError("Material $name: no constituent."))
    tot = sum(c.mass_fraction for c in constituents)
    tot <= 1 + 1.0e-6 || throw(ArgumentError("Material $name: the mass fractions of its constituents sum to $tot, above one."))
    e = _as_extent(extent)
    if kind in (:water, :salt)
        all(c -> c.extent isa ConstantExtent && c.extent.value == 1, constituents) &&
            e isa ConstantExtent && e.value == 1 || throw(
            ArgumentError("Material $name: water and salts react completely; an extent below one is refused.")
        )
    end
    return Material(String(name), kind, collect(AbstractConstituent, constituents), e, source === nothing ? nothing : String(source))
end

"""
    effective_extent(material, constituent, t) -> Real

The degree of reaction of `constituent` in `material` at time `t`: the extent of
the material times the extent of the constituent.
"""
effective_extent(m::Material, c::AbstractConstituent, t) = extent(m.extent, t) * extent(c.extent, t)

function Base.show(io::IO, m::Material)
    return print(io, "Material(\"", m.name, "\", :", m.kind, ", ", length(m.constituents), " constituents)")
end

# ── Oxides, and what a formula is worth in them ──────────────────────────────

# The element an oxide carries besides oxygen, and how many of its atoms per
# formula: ("Al", 2) for Al2O3, ("H", 2) for H2O, ("C", 1) for CO2.
function _oxide_cation(ox::AbstractString)
    comp = composition(Formula(ox))
    others = [(k, v) for (k, v) in comp if k != :O && !iszero(v)]
    (length(others) == 1 && haskey(comp, :O)) || throw(
        ArgumentError("\"$ox\" is not an oxide of a single element (CaO, Al2O3, SO3, H2O …).")
    )
    return first(only(others)), Float64(last(only(others)))
end
_oxide_molar_mass(ox::AbstractString) = ustrip(us"g/mol", calculate_molar_mass(composition(Formula(ox))))

"""
    oxide_content(species, oxides) -> OrderedDict{String, Float64}

The mass fraction of each oxide of `oxides` in the phase `species`, from its
formula and the molar masses of the library: the oxide analysis of a pure phase.
"""
function oxide_content(sp::AbstractSpecies, oxides)
    a = atoms(sp)
    M = ustrip(us"g/mol", sp[:M])
    out = OrderedDict{String, Float64}()
    for ox in oxides
        el, k = _oxide_cation(ox)
        out[String(ox)] = Float64(get(a, el, 0)) / k * _oxide_molar_mass(ox) / M
    end
    return out
end

"""
    literature_oxides(key, table, material) -> OrderedDict{String, Float64}

The oxide analysis of `material` in the table `table` of
`data/literature/<key>.json`, as mass **fractions**. Two layouts are read: a long
table with the columns `material`, `oxide` and `percent`
(`literature_oxides("Durdzinski2017", "chemical_composition", "S1")`), and a wide
one with one column per oxide (`literature_oxides("Gruyaert2010", "oxides",
"BFS-CAL")`). What is not an oxide (a loss on ignition, a fineness) is left out.
"""
function literature_oxides(key::AbstractString, table::AbstractString, material::AbstractString)
    t = literature_table(key, table)
    out = OrderedDict{String, Float64}()
    if haskey(t, :oxide)
        for (m, ox, p) in zip(t.material, t.oxide, t.percent)
            (m == material && _is_oxide(String(ox))) || continue   # not a loss on ignition
            out[String(ox)] = ustrip(p) / 100
        end
        isempty(out) && throw(KeyError("$key, table \"$table\", has no material \"$material\""))
        return out
    end
    row = literature_row(key, table, material)
    for (k, v) in pairs(row)
        col = String(k)
        v isa Union{Real, DynamicQuantities.AbstractQuantity} || continue
        _is_oxide(col) || continue
        out[col] = ustrip(v) / 100
    end
    return out
end
_is_oxide(s::AbstractString) = try
    _oxide_cation(s)
    true
catch
    false
end

"""
    oxide_material(name, oxides; kind = :scm, extent = 1.0, density, enthalpy, source)
    oxide_material(name, key, table, material; kwargs...)

A material described by its oxide analysis alone, one [`OxideConstituent`](@ref)
that reacts as a whole: a slag, a fly ash, a calcined clay. The second form reads
the analysis with [`literature_oxides`](@ref) and records the source.
"""
function oxide_material(
        name::AbstractString, oxides::AbstractDict; kind::Symbol = :scm, extent = 1.0,
        density = nothing, enthalpy = nothing, source = nothing
    )
    c = OxideConstituent(name, oxides; mass_fraction = 1.0, extent = 1.0, density, enthalpy, source)
    return Material(name, kind; constituents = [c], extent, source)
end
oxide_material(name::AbstractString, key::AbstractString, table::AbstractString, material::AbstractString; kwargs...) =
    oxide_material(name, literature_oxides(key, table, material); source = "$key:$table:$material", kwargs...)

# ── From an oxide analysis to phases ──────────────────────────────────────────

"""
    bogue(oxides, species; sulfate = "Gp") -> (; phases, remainder)

Bogue's calculation, from the formulas of the library rather than from its
printed coefficients: the mass fractions of the four clinker phases (`C3S`,
`C2S`, `C3A`, `C4AF`), of the calcium sulfate `sulfate` (gypsum `Gp` by default,
anhydrite `Anh`, one sulfur per formula) carrying the `SO3` and of calcite
(`Cal`) carrying the `CO2`, that reproduce the oxide analysis `oxides` (mass
fractions). Gypsum brings its water, which an oxide analysis counts in the loss
on ignition: take the carrier the cement actually holds.
`species` maps those symbols to database species (a `Dict` or a
[`ChemicalSystem`](@ref)). The lime of the gypsum and of the calcite is taken off
first; the four clinker phases then reproduce CaO, SiO2, Al2O3 and Fe2O3 exactly.
`remainder` holds the oxides no phase takes (MgO, the alkalis, TiO2 …), for an
[`OxideConstituent`](@ref).

A negative fraction means the analysis is not a Portland clinker's, and is
refused.
"""
function bogue(oxides::AbstractDict, species; sulfate::AbstractString = "Gp")
    sp(s) = species[s]
    x(ox) = float(get(oxides, ox, 0.0))
    get(atoms(sp(sulfate)), :S, 0) == 1 && get(atoms(sp(sulfate)), :Ca, 0) == 1 || throw(
        ArgumentError("bogue: the sulfate carrier $sulfate is not a calcium sulfate with one sulfur.")
    )
    gyp = x("SO3") * ustrip(us"g/mol", sp(sulfate)[:M]) / _oxide_molar_mass("SO3")
    cal = x("CO2") * ustrip(us"g/mol", sp("Cal")[:M]) / _oxide_molar_mass("CO2")
    cao = x("CaO") - x("SO3") * _oxide_molar_mass("CaO") / _oxide_molar_mass("SO3") -
        x("CO2") * _oxide_molar_mass("CaO") / _oxide_molar_mass("CO2")
    clinker = ("C3S", "C2S", "C3A", "C4AF")
    main = ("CaO", "SiO2", "Al2O3", "Fe2O3")
    Mx = [oxide_content(sp(ph), main)[ox] for ox in main, ph in clinker]
    y = Mx \ [cao, x("SiO2"), x("Al2O3"), x("Fe2O3")]
    # In the number type of the analysis: a composition being differentiated
    # carries its derivatives into the phases.
    phases = OrderedDict{String, promote_type(eltype(y), typeof(gyp), typeof(cal))}()
    for (ph, v) in zip(clinker, y)
        phases[ph] = v
    end
    phases[String(sulfate)] = gyp
    phases["Cal"] = cal
    bad = [k for (k, v) in phases if v < -1.0e-9]
    isempty(bad) || throw(
        ArgumentError(
            "bogue: the analysis gives negative fractions of $(join(bad, ", ")); it is not " *
                "the analysis of a Portland clinker, which Bogue's calculation assumes."
        )
    )
    taken = Set(("CaO", "SiO2", "Al2O3", "Fe2O3", "SO3", "CO2"))
    remainder = OrderedDict(String(k) => float(v) for (k, v) in oxides if !(String(k) in taken) && v > 0)
    return (; phases, remainder)
end

"""
    decompose(oxides, phases, species) -> (; fractions, misfit, remainder)

The non-negative mass fractions of the phases `phases` (symbols, resolved in
`species`) whose oxides best reproduce the analysis `oxides`, in the least
squares sense over the oxides the analysis gives: the generalization of Bogue's
calculation to any set of phases, where a square system would give negative
fractions or have none. `misfit` is the largest absolute difference between the
analysis and the oxides of the decomposition (a mass fraction), `remainder` the
oxide fractions left over, clipped at zero.
"""
function decompose(oxides::AbstractDict, phases, species)
    ox = [String(k) for k in keys(oxides)]
    target = [float(oxides[k]) for k in keys(oxides)]
    A = reduce(hcat, [[oxide_content(species[p], ox)[o] for o in ox] for p in phases])
    f = _nnls(A, target)
    r = target .- A * f
    fractions = OrderedDict(String(p) => v for (p, v) in zip(phases, f))
    remainder = OrderedDict(o => max(v, 0.0) for (o, v) in zip(ox, r) if v > 1.0e-12)
    return (; fractions, misfit = maximum(abs, r; init = 0.0), remainder)
end

"""
    reactive_part(oxides, crystalline, species) -> (; oxides, mass_fraction, clipped)

The glass of a material by difference: its oxide analysis minus the oxides of
its crystalline phases (`crystalline` maps phase symbols to mass fractions of the
material, as a Rietveld analysis gives them). `oxides` of the result are
fractions **of the glass**, `mass_fraction` the glass's share of the material;
an oxide the crystals take more of than the analysis holds is clipped at zero and
listed in `clipped`.
"""
function reactive_part(oxides::AbstractDict, crystalline::AbstractDict, species)
    F = promote_type(Float64, (typeof(float(v)) for v in values(oxides))..., (typeof(float(w)) for w in values(crystalline))...)
    rest = OrderedDict{String, F}(String(k) => float(v) for (k, v) in oxides)
    for (ph, w) in crystalline
        c = oxide_content(species[String(ph)], keys(rest))
        for k in keys(rest)
            rest[k] -= w * c[k]
        end
    end
    clipped = [k for (k, v) in rest if v < -1.0e-9]
    frac = 1.0 - sum(values(crystalline); init = 0.0)
    frac > 0 || throw(ArgumentError("reactive_part: the crystalline phases take the whole material."))
    glass = OrderedDict(k => max(v, 0.0) / frac for (k, v) in rest)
    return (; oxides = glass, mass_fraction = frac, clipped)
end

# Lawson and Hanson's active-set method for min ‖Ax − b‖, x ≥ 0 (their Chapter
# 23), for the few columns a decomposition has.
function _nnls(A::AbstractMatrix, b::AbstractVector; tol = 1.0e-12, maxit = 30 * size(A, 2))
    n = size(A, 2)
    # In the number type of the data; the active set is chosen on the values,
    # and the answer is then the least-squares solution on it, derivatives and
    # all.
    T = promote_type(Float64, eltype(A), eltype(b))
    x = zeros(T, n)
    P = falses(n)
    w = transpose(A) * (b .- A * x)
    it = 0
    while any(.!P .& (w .> tol)) && it < maxit
        it += 1
        j = argmax([P[i] ? -Inf : _plain(w[i]) for i in 1:n])
        P[j] = true
        while true
            idx = findall(P)
            z = zeros(T, n)
            z[idx] = A[:, idx] \ b
            if all(z[idx] .> tol)
                x = z
                break
            end
            # Step back to the boundary of the feasible set along x → z.
            α = minimum(x[i] / max(x[i] - z[i], floatmin()) for i in idx if z[i] <= tol)
            x .+= α .* (z .- x)
            for i in idx
                x[i] <= tol && (P[i] = false; x[i] = 0.0)
            end
        end
        w = transpose(A) * (b .- A * x)
    end
    return x
end
