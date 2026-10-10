# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

# ── The enthalpy and heat capacity of a glass known by its oxides ─────────────
#
# A slag or a fly ash is mostly glass, and a glass has no formula and no entry in
# a thermodynamic database: what releases heat when it reacts in a paste has no
# enthalpy of formation to subtract. The silicate glasses of geochemistry have
# been measured, though: their heat capacity is additive in their oxides (Richet
# 1987), and the enthalpy by which a glass exceeds its crystal, its enthalpy of
# vitrification, is known by solution calorimetry for a few compositions near a
# slag's (Richet and Bottinga 1984, 1986; Navrotsky et al. 1980). This file
# builds a glass from those, and says what it had to assume.

using LinearAlgebra: I, det

# The glasses whose enthalpy of vitrification is measured at 298 K, each with the
# crystal it refers to (its symbol in the aq17 database) and where the value is
# read: the table, the column to match and its value.
const _MEASURED_GLASSES = (
    (name = "gehlenite", formula = "Ca2Al2SiO7", crystal = "Gehlenite", key = "RichetBottinga1986", match = (formula = "Ca2Al2SiO7",)),
    (name = "akermanite", formula = "Ca2MgSi2O7", crystal = "Akermanite", key = "RichetBottinga1986", match = (formula = "Ca2MgSi2O7",)),
    (name = "pseudowollastonite", formula = "CaSiO3", crystal = "Pseudowoll", key = "RichetBottinga1986", match = (formula = "CaSiO3",)),
    (name = "anorthite", formula = "CaAl2Si2O8", crystal = "Anorthite", key = "RichetBottinga1984", match = (crystal = "anorthite", reference = "Kracek and Neuvonen [26]")),
    (name = "diopside", formula = "CaMgSi2O6", crystal = "Diopside", key = "RichetBottinga1984", match = (crystal = "diopside", reference = "Kracek [33]")),
    (name = "silica", formula = "SiO2", crystal = "Quartz", key = "Navrotsky1980", match = (substance = "quartz", reference = "Kracek (1953)")),
)

# The four oxides the measured glasses are made of, and the crystal that stands
# for each in aq17.
const _GLASS_OXIDES = ("CaO", "MgO", "Al2O3", "SiO2")
const _OXIDE_CRYSTALS_AQ17 = Dict("CaO" => "Lime", "MgO" => "Periclase", "Al2O3" => "Corundum", "SiO2" => "Quartz")

# The vitrification enthalpy (J/mol) of each measured glass at 298 K and its
# uncertainty (J/mol, NaN where none is printed), read from the transcriptions.
function _vitrification_enthalpies()
    return map(_MEASURED_GLASSES) do g
        col = g.key == "Navrotsky1980" ? :T : :T_s
        t = literature_table(g.key, "vitrification_enthalpies"; g.match..., (col => 298.0u"K",)...)
        length(t.dH_v) == 1 || error("glass_enthalpy: $(g.key) holds $(length(t.dH_v)) rows for $(g.name) at 298 K, one expected.")
        h = ustrip(us"J/mol", only(t.dH_v))
        x = haskey(t, :dH_v_uncertainty) ? only(t.dH_v_uncertainty) : nothing
        u = (x === nothing || x === missing) ? NaN : ustrip(us"J/mol", x)
        (; g..., dH_v = h, uncertainty = u)
    end
end

# The oxides of a formula, as moles per formula unit in the order of _GLASS_OXIDES.
function _glass_oxide_counts(formula::AbstractString)
    a = atoms(Species(formula))
    return [Float64(get(a, :Ca, 0)), Float64(get(a, :Mg, 0)), Float64(get(a, :Al, 0)) / 2, Float64(get(a, :Si, 0))]
end

# The crystals of aq17 and the default reference, built once, under a lock as
# `_gel_models` is: two threads asking first would both build them and race on the
# reference.
const _GLASS_DATA_LOCK = ReentrantLock()
const _AQ17 = Ref{Any}(nothing)
function _aq17_crystals()
    return lock(_GLASS_DATA_LOCK) do
        if _AQ17[] === nothing
            _AQ17[] = Dict(symbol(s) => s for s in build_species(datapath("aq17-thermofun.json"); verbose = false))
        end
        _AQ17[]
    end
end

# The default reference for the enthalpies of formation of the oxides: aq17,
# then slop98, built once.
const _GLASS_REFERENCE = Ref{Any}(nothing)
function _default_glass_reference()
    return lock(_GLASS_DATA_LOCK) do
        if _GLASS_REFERENCE[] === nothing
            _GLASS_REFERENCE[] = (
                collect(values(_aq17_crystals())),
                build_species(datapath("slop98-inorganic-thermofun.json"); verbose = false),
            )
        end
        _GLASS_REFERENCE[]
    end
end

_formation_enthalpy(sp) = ustrip(us"J/mol", sp[:ΔₐH⁰](T = T_STANDARD_Q, P = P_STANDARD_Q; unit = true))

# The crystal of an oxide among `reference`, an ordered list of species
# collections: in the first collection that holds a crystal of that formula, the
# one of lowest standard Gibbs energy at 298.15 K, the stable polymorph (quartz,
# not coesite).
function _reference_oxide(formula::AbstractString, reference)
    target = atoms(Species(formula))
    for (k, db) in enumerate(reference)
        best, gbest = nothing, Inf
        for sp in db
            (aggregate_state(sp) == AS_CRYSTAL && charge(sp) == 0 && atoms(sp) == target) || continue
            g = ustrip(us"J/mol", sp[:ΔₐG⁰](T = T_STANDARD_Q, P = P_STANDARD_Q; unit = true))
            g < gbest && ((best, gbest) = (sp, g))
        end
        best === nothing || return (best, k)
    end
    return nothing
end

"""
    glass_enthalpy(oxides; T = 298.15u"K", reference = nothing, ignore = ()) -> NamedTuple

The standard enthalpy of formation of a silicate glass known by its oxide
analysis, per gram of material, built from measurements on silicate glasses and
stating what it assumes. `oxides` maps oxide formulas to mass fractions, as for
[`oxide_budget`](@ref) (not renormalized).

The CaO-MgO-Al2O3-SiO2 part of the glass is written as a combination of the six
glasses whose enthalpy of vitrification is measured at 298 K, gehlenite,
akermanite, pseudowollastonite, anorthite, diopside and silica
([RichetBottinga1986](@citet), [RichetBottinga1984](@citet),
[Navrotsky1980](@citet)), mixed ideally. The combination is chosen to leave as
little of the four oxides outside the measured glasses as possible; what is left
(a slag holds more lime than the measured glasses can take) is counted as the
crystalline oxide, its enthalpy of vitrification assumed zero. When several
combinations leave the same amount, their enthalpies differ by the measure of the
ideal-mixing assumption, and the midpoint is returned with the half-range as
`span`. The enthalpies of the crystals and of the four oxides are taken from the
aq17 database together, so that the enthalpy of forming the glass from its
oxides is consistent within one database; the enthalpies of formation of the
oxides themselves come from `reference`, an ordered list of species collections
searched in turn (by default aq17 then slop98 [Johnson1992](@cite)), so that the
glass can be put on the scale of the database of a calculation. The other oxides
of the analysis (Fe2O3, Na2O, K2O, TiO2, …) are counted as their crystals from
`reference`, with no enthalpy of vitrification (assumed).

At a temperature `T` other than 298.15 K the heat capacity of the glass,
[`glass_heat_capacity`](@ref), is integrated from 298.15 K.

An oxide with no crystal in `reference` is refused, unless it is listed in
`ignore`: it is then left out, its mass counted at zero enthalpy, and reported in
`ignored` with its mass fraction.

Returns `(; enthalpy, formation_298, from_oxides, vitrification, span,
uncertainty, norm, unassigned, others, ignored, oxide_sources)`: the enthalpy of formation
at `T` and at 298.15 K, the part from the oxides to the glass, the vitrification
part and its `span`, the `uncertainty` the printed uncertainties of the measured
enthalpies of vitrification give (those printed), all in J/g; `norm` the moles
of each measured glass per gram, `unassigned` the moles of the four oxides left
crystalline per gram, `others` the moles per gram of the oxides outside them,
and `oxide_sources` the database each oxide's enthalpy came from (its position
in `reference`). The construction is written out in
[The enthalpy of a glass](@ref sec-theory-glass).

# Examples

```julia
slag = Dict("CaO" => 0.416, "SiO2" => 0.352, "Al2O3" => 0.113, "MgO" => 0.060)
g = glass_enthalpy(slag)
OxideConstituent("glass", slag; enthalpy = g.enthalpy, source = "glass_enthalpy")
```
"""
function glass_enthalpy(oxides::AbstractDict{<:AbstractString, <:Real}; T = T_STANDARD_Q, reference = nothing, ignore = ())
    aq = _aq17_crystals()
    refs = reference === nothing ? _default_glass_reference() : reference
    glasses = _vitrification_enthalpies()

    # Moles of each oxide per gram of material, in the number type of the
    # fractions (a dual fraction is a composition being calibrated).
    R = mapreduce(f -> typeof(float(f)), promote_type, values(oxides); init = Float64)
    n = Dict{String, R}()
    ignored = Dict{String, R}()
    for (ox, f) in oxides
        f == 0 && continue
        f < 0 && throw(ArgumentError("glass_enthalpy: the mass fraction of `$ox` is negative ($f)."))
        if ox in ignore
            ignored[ox] = float(f)
            continue
        end
        n[ox] = float(f) / ustrip(us"g/mol", Species(ox)[:M])
    end
    isempty(n) && throw(ArgumentError("glass_enthalpy: the analysis carries no oxide with a positive mass fraction."))
    b = [get(n, ox, zero(R)) for ox in _GLASS_OXIDES]

    # The columns of the small linear program: the measured glasses, then one
    # crystalline oxide per oxide (the part no measured glass covers).
    A = hcat((_glass_oxide_counts(g.formula) for g in glasses)..., Matrix{Float64}(I, 4, 4))
    ng = length(glasses)
    cost = vcat(zeros(ng), ones(4))
    dHv = vcat([g.dH_v for g in glasses], zeros(4))
    sol = _norm_vertices(A, b)
    best = minimum(v -> cost' * v, sol)
    tol = 1.0e-12 * max(1.0, sum(b))
    optimal = filter(v -> cost' * v <= best + tol, sol)
    vit = [dHv' * v for v in optimal]
    lo, hi = extrema(vit)
    k = argmin(abs.(vit .- (lo + hi) / 2))
    ν = optimal[k]

    # The enthalpy of forming each measured crystal from its oxides, within aq17.
    hox = Dict(ox => _formation_enthalpy(aq[_OXIDE_CRYSTALS_AQ17[ox]]) for ox in _GLASS_OXIDES)
    from_ox_crystal = [
        _formation_enthalpy(aq[g.crystal]) - sum(_glass_oxide_counts(g.formula) .* [hox[ox] for ox in _GLASS_OXIDES])
            for g in glasses
    ]
    vitrification = (lo + hi) / 2
    from_oxides = sum(ν[1:ng] .* from_ox_crystal) + vitrification
    unc2 = sum((ν[i] * glasses[i].uncertainty)^2 for i in 1:ng if ν[i] > tol && !isnan(glasses[i].uncertainty); init = 0.0)

    # The enthalpies of formation of all the oxides of the analysis.
    sources = Dict{String, Int}()
    formation = from_oxides
    for (ox, m) in n
        found = _reference_oxide(ox, refs)
        found === nothing && throw(
            ArgumentError(
                "glass_enthalpy: no crystal of $ox in the reference databases; give `reference` a " *
                    "collection of species that holds one, or leave it out with `ignore = (\"$ox\",)`, " *
                    "which the result then reports."
            )
        )
        sp, k = found
        sources[ox] = k
        formation += m * _formation_enthalpy(sp)
    end

    TK = T isa Real ? float(T) : ustrip(us"K", T)
    enthalpy = formation
    if TK != T_STANDARD
        kept = Dict(ox => f for (ox, f) in oxides if !(ox in ignore))
        cp(θ) = ustrip(us"J/(g*K)", glass_heat_capacity(kept; T = θ * u"K", reference = refs))
        nodes, weights = _gauss_legendre_8()
        half = (TK - T_STANDARD) / 2
        mid = (TK + T_STANDARD) / 2
        enthalpy += half * sum(w * cp(mid + half * x) for (x, w) in zip(nodes, weights))
    end

    return (;
        enthalpy = enthalpy * u"J/g", formation_298 = formation * u"J/g", from_oxides = from_oxides * u"J/g",
        vitrification = vitrification * u"J/g", span = (hi - lo) / 2 * u"J/g", uncertainty = sqrt(unc2) * u"J/g",
        norm = Dict(glasses[i].name => ν[i] for i in 1:ng if ν[i] > tol),
        unassigned = Dict(_GLASS_OXIDES[i] => ν[ng + i] for i in 1:4 if ν[ng + i] > tol),
        others = Dict(ox => m for (ox, m) in n if !(ox in _GLASS_OXIDES)),
        ignored = ignored,
        oxide_sources = sources,
    )
end

# The vertices of {v ≥ 0 : A v = b}: every choice of as many columns as rows
# whose square system has a non-negative solution. A linear objective over the
# polytope reaches its extremes at vertices, and the problem here is small (four
# rows, ten columns, 210 choices).
function _norm_vertices(A::AbstractMatrix, b::AbstractVector)
    m, n = size(A)
    R = float(promote_type(eltype(A), eltype(b)))
    out = Vector{Vector{R}}()
    scale = max(1.0, maximum(abs, b))
    for cols in _combinations(n, m)
        B = A[:, cols]
        abs(det(B)) < 1.0e-12 && continue
        x = B \ b
        all(>=(-1.0e-12 * scale), x) || continue
        v = zeros(R, n)
        v[cols] .= max.(x, zero(R))
        any(w -> isapprox(w, v; atol = 1.0e-14 * scale), out) || push!(out, v)
    end
    isempty(out) && error("glass_enthalpy: no non-negative combination; this cannot happen with the crystalline oxides among the columns.")
    return out
end

function _combinations(n::Int, k::Int)
    k == 0 && return [Int[]]
    out = Vector{Vector{Int}}()
    for first in 1:(n - k + 1), rest in _combinations(n - first, k - 1)
        push!(out, vcat(first, rest .+ first))
    end
    return out
end

# Eight-point Gauss-Legendre on [-1, 1]: the integral of a heat capacity, a
# smooth function, over a few tens of kelvins is exact to rounding.
function _gauss_legendre_8()
    x = (
        -0.9602898564975363, -0.7966664774136267, -0.525532409916329, -0.1834346424956498,
        0.1834346424956498, 0.525532409916329, 0.7966664774136267, 0.9602898564975363,
    )
    w = (
        0.1012285362903763, 0.2223810344533745, 0.3137066458778873, 0.362683783378362,
        0.362683783378362, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763,
    )
    return x, w
end

"""
    glass_heat_capacity(oxides; T = 298.15u"K", reference = nothing) -> Quantity

The heat capacity of a silicate glass per gram, from its oxide analysis, as
[Richet1987](@citet) gives it: the sum over the oxides of their moles times
their partial molar heat capacity in the glass, `aᵢ + bᵢT + cᵢ/T² + dᵢ/√T` (his
Eqs. 1 and 3, Table I), to about 1 % below the glass transition and between
270 K and the upper temperature of each oxide's data. An oxide the paper does not
cover takes the heat capacity of its crystal (the paper's footnote), from
`reference` as in [`glass_enthalpy`](@ref).
"""
function glass_heat_capacity(oxides::AbstractDict{<:AbstractString, <:Real}; T = T_STANDARD_Q, reference = nothing)
    TK = T isa Real ? float(T) : ustrip(us"K", T)
    t = literature_table("Richet1987", "oxide_heat_capacity")
    cp = 0.0
    for (ox, f) in oxides
        f == 0 && continue
        f < 0 && throw(ArgumentError("glass_heat_capacity: the mass fraction of `$ox` is negative ($f)."))
        m = float(f) / ustrip(us"g/mol", Species(ox)[:M])
        i = findfirst(==(ox), t.oxide)
        if i === nothing
            found = _reference_oxide(ox, reference === nothing ? _default_glass_reference() : reference)
            found === nothing && throw(
                ArgumentError("glass_heat_capacity: $ox is not in Table I of Richet (1987) and no crystal of it is in the reference databases.")
            )
            cp += m * ustrip(us"J/(mol*K)", first(found)[:Cp⁰](T = TK * u"K", P = P_STANDARD_Q; unit = true))
        else
            a, bb, c, d = (ustrip(t.a[i]), ustrip(t.b[i]), ustrip(t.c[i]), ustrip(t.d[i]))
            cp += m * (a + bb * TK + c / TK^2 + d / sqrt(TK))
        end
    end
    return cp * u"J/(g*K)"
end
