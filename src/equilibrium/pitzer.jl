# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

# The Pitzer ion-interaction model. Unlike the extended Debye-Hückel models in
# `activities.jl`, this one is not a per-species formula: the excess Gibbs energy
# is a virial expansion, so a coefficient belongs to a *pair* of ions and another
# to a *triplet*, and the sums run over the whole solution. That is why it needs
# its own closure rather than a method of `_log10γ_ion`. Its parameters are
# caller input: the ThermoFun databases this package reads carry none, PHREEQC
# distributes a set in its `pitzer.dat`, and `data/pitzer-reardon1990.toml` ships
# the cement set of Reardon (1990).

using DynamicQuantities

"""
    PitzerParameters(; beta0, beta1, beta2, Cphi, theta, psi, lambda,
                       alpha1_22 = 1.4, alpha1 = 2.0, alpha2 = 12.0, b = 1.2,
                       origin = Dict{Tuple{String, String}, String}(),
                       temperature = Dict{Symbol, Any}())

The interaction parameters of a Pitzer model, keyed by species symbol.

Every keyword is **mandatory**, `origin` and the four shape constants excepted:
a Pitzer model is only as good as the set it is given, and a default would be a
number invented on the caller's behalf. Omitting one raises `UndefKeywordError`
before anything is computed.

# The tables

| keyword | keys | meaning |
|:--|:--|:--|
| `beta0`, `beta1`, `beta2` | `(cation, anion)` | the second virial coefficient of that pair, as three terms of one ionic-strength dependence |
| `Cphi` | `(cation, anion)` | the third virial coefficient of that pair |
| `theta` | `(ion, ion)` of **like** charge | interaction between two cations, or two anions |
| `psi` | `(ion, ion, ion)` | a triplet: two anions with a cation, or two cations with an anion |
| `lambda` | `(neutral, ion)` | a neutral solute with an ion — salting out, in Pitzer's form |

`theta`, `psi` and `lambda` entries are looked up in any order of their keys.
A **missing** `theta`, `psi` or `lambda` is taken as zero, which is the
convention of the literature these tables come from; a **missing** `beta0` for a
pair the system actually contains is an error, raised when the model is attached
to a [`ChemicalSystem`](@ref) — see [`PitzerActivityModel`](@ref).

# The shape constants

`alpha1 = 2.0`, `alpha1_22 = 1.4` (used when both ions carry a charge of
magnitude 2 or more), `alpha2 = 12.0` and `b = 1.2` kg^½·mol^−½ are **not**
fitted parameters: they fix the functional form, and the published `β` tables
were fitted assuming them. Changing one invalidates the table it is used with.
They are keywords so that a set fitted with a different convention can be used,
not so that they can be tuned.

# Temperature

The tables give each coefficient at 25 °C, `A₀`. `temperature` may give, for any
of them, the further coefficients of the form PHREEQC uses
[ParkhurstAppelo2013](@cite),

```math
P(T) = A_0 + A_1\\Big(\\frac1T - \\frac1{T_r}\\Big) + A_2 \\ln\\frac{T}{T_r}
     + A_3 (T - T_r) + A_4 (T^2 - T_r^2) + A_5\\Big(\\frac1{T^2} - \\frac1{T_r^2}\\Big),
\\qquad T_r = 298.15\\ \\mathrm{K},
```

keyed by the table and then by the entry's key: `temperature = Dict(:beta0 =>
Dict(("Na+", "Cl-") => (A₁, A₂, A₃, A₄, A₅)))`, fewer than five meaning the rest
are zero. They are used by a [`PitzerActivityModel`](@ref) built with
`temperature_dependent = true`, and vanish at ``T_r`` exactly.

# Provenance

`origin` maps a pair to a free-text note, and the loader
[`build_pitzer_parameters`](@ref) fills it from the data file. It exists because
a published set can contain values that were *estimated by analogy* rather than
measured, and a caller reading `beta0` has otherwise no way of telling the two
apart. [`pitzer_origin`](@ref) reads it back.

See also: [`PitzerActivityModel`](@ref), [`build_pitzer_parameters`](@ref).
"""
struct PitzerParameters{T <: Real}
    beta0::Dict{Tuple{String, String}, T}
    beta1::Dict{Tuple{String, String}, T}
    beta2::Dict{Tuple{String, String}, T}
    Cphi::Dict{Tuple{String, String}, T}
    theta::Dict{Tuple{String, String}, T}
    psi::Dict{Tuple{String, String, String}, T}
    lambda::Dict{Tuple{String, String}, T}
    alpha1::T
    alpha1_22::T
    alpha2::T
    b::T
    origin::Dict{Tuple{String, String}, String}
    temperature::Dict{Symbol, Dict{Any, NTuple{5, T}}}

    # Inner, so that the generated constructor cannot bypass the checks: a
    # `beta1` for a pair with no `beta0` is a typo in the caller's table, and a
    # `theta` between ions of opposite charge is a category error — that
    # interaction is what `beta` is for.
    function PitzerParameters{T}(
            beta0, beta1, beta2, Cphi, theta, psi, lambda,
            alpha1, alpha1_22, alpha2, b, origin, temperature,
        ) where {T <: Real}
        for kind in keys(temperature)
            kind in _PITZER_TABLES || throw(
                ArgumentError(
                    "PitzerParameters: temperature terms for :$kind, which is not one of " *
                        "the tables $(join(_PITZER_TABLES, ", ")).",
                )
            )
        end
        for (name, tbl) in (("beta1", beta1), ("beta2", beta2), ("Cphi", Cphi))
            for k in keys(tbl)
                haskey(beta0, k) || throw(
                    ArgumentError(
                        "PitzerParameters: $name has an entry for $k but beta0 does not. " *
                            "A pair is described by all four coefficients or by none.",
                    )
                )
            end
        end
        b > 0 || throw(ArgumentError("PitzerParameters: b must be positive, got $b."))
        alpha1 > 0 && alpha2 > 0 && alpha1_22 > 0 || throw(
            ArgumentError("PitzerParameters: the alpha constants must be positive.")
        )
        return new{T}(
            beta0, beta1, beta2, Cphi, theta, psi, lambda,
            alpha1, alpha1_22, alpha2, b, origin, temperature,
        )
    end
end

const _PITZER_TABLES = (:beta0, :beta1, :beta2, :Cphi, :theta, :psi, :lambda)

# The reference temperature of the temperature terms, PHREEQC's.
const _PITZER_TR = 298.15

# The five temperature coefficients A₁…A₅ of an entry, fewer meaning zeros.
function _five(::Type{T}, v) where {T}
    length(v) <= 5 || throw(
        ArgumentError("a Pitzer coefficient has at most five temperature terms; got $(length(v))."),
    )
    return ntuple(k -> k <= length(v) ? T(v[k]) : zero(T), 5)
end

function PitzerParameters(;
        beta0, beta1, beta2, Cphi, theta, psi, lambda,
        alpha1::Real = 2.0, alpha1_22::Real = 1.4, alpha2::Real = 12.0, b::Real = 1.2,
        origin::AbstractDict = Dict{Tuple{String, String}, String}(),
        temperature::AbstractDict = Dict{Symbol, Any}(),
    )
    Tterm = isempty(temperature) ? Float64 :
        mapreduce(d -> isempty(d) ? Float64 : mapreduce(v -> promote_type(map(typeof, v)...), promote_type, values(d)), promote_type, values(temperature))
    T = float(
        promote_type(
            eltype(values(beta0)), eltype(values(beta1)), eltype(values(beta2)),
            eltype(values(Cphi)), eltype(values(theta)), eltype(values(psi)),
            eltype(values(lambda)), typeof(alpha1), typeof(alpha1_22),
            typeof(alpha2), typeof(b), Tterm,
        )
    )
    conv(d, K) = Dict{K, T}(k => T(v) for (k, v) in d)
    P2 = Tuple{String, String}
    P3 = Tuple{String, String, String}
    return PitzerParameters{T}(
        conv(beta0, P2), conv(beta1, P2), conv(beta2, P2), conv(Cphi, P2),
        conv(theta, P2), conv(psi, P3), conv(lambda, P2),
        T(alpha1), T(alpha1_22), T(alpha2), T(b),
        Dict{P2, String}(k => String(v) for (k, v) in origin),
        Dict{Symbol, Dict{Any, NTuple{5, T}}}(
            Symbol(kind) => Dict{Any, NTuple{5, T}}(k => _five(T, v) for (k, v) in d)
                for (kind, d) in temperature
        ),
    )
end

"""
    _pitzer_tables_at(T, tables0, tablesT) -> tables

The coefficient tables at the temperature `T`: each entry `A₀` of `tables0`
moved by its temperature terms in `tablesT` (see [`PitzerParameters`](@ref)).
Each term vanishes at `T_r` exactly, so the tables are there their own values.
"""
function _pitzer_tables_at(T, tables0, tablesT)
    Tr = _PITZER_TR
    basis = (1 / T - 1 / Tr, log(T / Tr), T - Tr, T^2 - Tr^2, 1 / T^2 - 1 / Tr^2)
    at(a0, t) = a0 + (t[1] * basis[1] + t[2] * basis[2] + t[3] * basis[3] + t[4] * basis[4] + t[5] * basis[5])
    return map((a0, at_) -> map(at, a0, at_), tables0, tablesT)
end

"""
    pitzer_origin(p::PitzerParameters, cation, anion) -> String

What the parameters of that pair are: `"fitted"`, `"estimated:<analog>"`, or
`"unrecorded"` when the set carries no note.

A published Pitzer set can contain values obtained by analogy with a chemically
similar ion rather than by fitting a measurement — the cement set of
[Reardon1990](@citet) does so for every silicate, aluminate and ferrate pair,
which are exactly the ions a cement assemblage needs. This accessor is how that
shows up in a calculation instead of staying in a paper's footnote.
"""
pitzer_origin(p::PitzerParameters, cation::AbstractString, anion::AbstractString) =
    get(p.origin, (String(cation), String(anion)), "unrecorded")

# Symmetric lookups: the caller's table may key a pair or a triplet either way
# round, and the physics does not distinguish them.
@inline function _sym2(d, i, j)
    v = get(d, (i, j), nothing)
    v === nothing || return v
    return get(d, (j, i), zero(valtype(d)))
end
@inline function _sym3(d, i, j, k)
    for key in ((i, j, k), (i, k, j), (j, i, k), (j, k, i), (k, i, j), (k, j, i))
        v = get(d, key, nothing)
        v === nothing || return v
    end
    return zero(valtype(d))
end

"""
    PitzerActivityModel(; parameters, temperature_dependent = false, etheta = true,
                          missing_pairs = :refuse, A = $(_DH_A_25C))

The Pitzer ion-interaction activity model.

`parameters` is a [`PitzerParameters`](@ref) and has **no default**: this model
cannot be reached without one, which is the point. Whether the set is complete
is not a property of the set alone but of the set *relative to a species list*,
so the check happens when the model meets a [`ChemicalSystem`](@ref): every
cation-anion pair the system contains must have a `beta0` entry, and the error
names the pairs that do not.

`missing_pairs = :zero` takes the coefficients of such a pair as zero instead,
which is what PHREEQC does with a database whose `PITZER` block leaves pairs out
(`pitzer.dat` describes neither H⁺–OH⁻ nor Ca²⁺–CO₃²⁻): it is the convention of
the database, and [`database_activity_model`](@ref) builds the model of such a
database with it.

# Why this model rather than an extended Debye-Hückel one

`γ` and the osmotic coefficient are derived from **one** excess Gibbs energy — a
virial expansion in the molalities, with one coefficient per ion pair and one
per triplet — so the Gibbs-Duhem relation between the solutes and the solvent
holds by construction rather than approximately. That is what
[`DaviesActivityModel`](@ref) does not do and [`HKFActivityModel`](@ref) does
only up to a mean-radius approximation. See
[Activity models](@ref sec-theory-activity).

# The higher-order electrostatic terms

Two ions of like sign and *unsymmetrical* charge (Na⁺ with Ca²⁺, Cl⁻ with SO₄²⁻)
interact through ``{}^E\\theta_{ij}(I)`` besides ``\\theta_{ij}``, the terms of
[Pitzer1975](@citet): with ``x_{ij} = 6 z_i z_j A_\\varphi \\sqrt{I}``,

```math
{}^E\\theta_{ij} = \\frac{z_i z_j}{4I}\\Big[J(x_{ij}) - \\tfrac12 J(x_{ii}) - \\tfrac12 J(x_{jj})\\Big] ,
```

and ``{}^E\\theta'_{ij}`` its derivative with respect to ``I``. ``J`` is evaluated by
the Chebyshev approximation of Harvie, as PHREEQC evaluates it
[Plummer1988](@cite). They vanish for a symmetrical pair, so a single
electrolyte is unaffected, and they enter ``\\gamma`` and the osmotic coefficient
from the same excess Gibbs energy, so that the Gibbs–Duhem relation keeps holding
exactly. `etheta = false` leaves them out, for a set whose ``\\theta`` were fitted
without them (Pitzer and Mayorga's, for instance): ``\\theta`` and ``{}^E\\theta``
are fitted together, and a set is used with the convention it was fitted with.

# Temperature

`A` is the Debye–Hückel constant of the decimal logarithm, from which
``A_\\varphi = A \\ln 10 / 3``. It is used as given unless
`temperature_dependent = true`, which takes ``A_\\varphi`` from the water model at
the temperature of the state, and evaluates the temperature terms of the
coefficients when the set carries them (see [`PitzerParameters`](@ref)). A
published set fitted at one temperature carries none, and is then used at its
own values whatever the temperature, which is what the set can support.

# Example

```julia
p = build_pitzer_parameters(datapath("pitzer-reardon1990.toml"))
model = PitzerActivityModel(; parameters = p)
```

See also: [`PitzerParameters`](@ref), [`build_pitzer_parameters`](@ref),
[`pitzer_origin`](@ref).
"""
struct PitzerActivityModel{T <: Real, TA <: Real} <: AbstractActivityModel
    parameters::PitzerParameters{T}
    temperature_dependent::Bool
    etheta::Bool
    missing_pairs::Symbol
    A::TA
end

function PitzerActivityModel(;
        parameters, temperature_dependent::Bool = false, etheta::Bool = true,
        missing_pairs::Symbol = :refuse, A::Real = _DH_A_25C,
    )
    missing_pairs in (:refuse, :zero) ||
        throw(ArgumentError("missing_pairs is `:refuse` or `:zero`; got `:$missing_pairs`"))
    return PitzerActivityModel(parameters, temperature_dependent, etheta, missing_pairs, float(A))
end

concentration_scale(::PitzerActivityModel) = :molality

# ── The two ionic-strength functions of the model ────────────────────────────
#
# `g` and `g′` are 0/0 at zero ionic strength, so both carry a Taylor branch.
# The series are exact to the order written: with e^{-x} expanded,
# `1 - (1+x)e^{-x} = x²/2 - x³/3 + x⁴/8`, hence `g = 1 - 2x/3 + x²/4`, and
# `1 - (1+x+x²/2)e^{-x} = x³/6 - x⁴/8`, hence `g′ = -x/3 + x²/4`.
# Branching is on the primal so that a `ForwardDiff.Dual` keeps its derivative.
@inline function _pitzer_g(x::T) where {T <: Real}
    if abs(_primal(x)) < 1.0e-4
        return one(T) - (2 // 3) * x + x^2 / 4
    end
    return 2 * (one(T) - (one(T) + x) * exp(-x)) / x^2
end

@inline function _pitzer_gp(x::T) where {T <: Real}
    if abs(_primal(x)) < 1.0e-4
        return -x / 3 + x^2 / 4
    end
    return -2 * (one(T) - (one(T) + x + x^2 / 2) * exp(-x)) / x^2
end

# ── The higher-order electrostatic terms ────────────────────────────────────
#
# Pitzer's (1975) J(x), by the Chebyshev approximation of Harvie that PHREEQC
# uses (Plummer et al. 1988), with the coefficients Reaktoro carries in
# `computeJ0J1PHREEQC` (Reaktoro, LGPL-2.1-or-later; see NOTICE). The first 21
# serve x ≤ 1, the last 21 x > 1.
const _HARVIE_AK = (
    (
        1.925154014814667, -0.060076477753119, -0.029779077456514,
        -0.007299499690937, 0.000388260636404, 0.000636874599598,
        0.000036583601823, -0.000045036975204, -0.00000453789571,
        0.000002937706971, 0.000000396566462, -0.000000202099617,
        -0.000000025267769, 0.00000001352261, 0.000000001229405,
        -0.000000000821969, -0.000000000050847, 0.000000000046333,
        0.000000000001943, -0.000000000002563, -0.000000000010991,
    ),
    (
        0.628023320520852, 0.462762985338493, 0.150044637187895,
        -0.028796057604906, -0.036552745910311, -0.001668087945272,
        0.006519840398744, 0.001130378079086, -0.000887171310131,
        -0.000242107641309, 0.000087294451594, 0.000034682122751,
        -0.000004583768938, -0.000003548684306, -0.00000025045388,
        0.000000216991779, 0.00000008077957, 0.000000004558555,
        -0.000000006944757, -0.000000002849257, 0.000000000237816,
    ),
)

"""
    _pitzer_J(x) -> (J, xJ′)

Pitzer's function ``J(x)`` and ``x\\,J'(x)``, by Harvie's Chebyshev series in the
variable `z(x)`, evaluated by Clenshaw's recurrence with its derivative. The
derivative is that of the approximation itself, so that ``{}^E\\theta'`` is
exactly the derivative of ``{}^E\\theta`` and the model keeps deriving from one
energy. `x > 0`; a `ForwardDiff.Dual` passes through.
"""
function _pitzer_J(x::T) where {T <: Real}
    small = _primal(x) <= 1
    ak = small ? _HARVIE_AK[1] : _HARVIE_AK[2]
    # The series variable z(x), and x·dz/dx halved, as the derivative below
    # takes it.
    if small
        q = x^(1 // 5)
        z, hxdz = 4 * q - 2, (2 // 5) * q
    else
        q = x^(-1 // 10)
        z, hxdz = (40 * q - 22) / 9, -(2 // 9) * q
    end
    # Clenshaw: bₖ = z bₖ₊₁ − bₖ₊₂ + aₖ, and its derivative in z.
    b, b1, b2 = zero(z), zero(z), zero(z)
    d, d1, d2 = zero(z), zero(z), zero(z)
    bk2 = dk2 = zero(z)                 # b₂ and its derivative, kept on the way
    @inbounds for k in 21:-1:1
        b = z * b1 - b2 + ak[k]
        d = b1 + z * d1 - d2
        k == 3 && (bk2 = b; dk2 = d)
        b2, b1 = b1, b
        d2, d1 = d1, d
    end
    return (x / 4 - 1 + (b - bk2) / 2, x / 4 + hxdz * (d - dk2))
end

"""
    _etheta_pairs(z, Aφ, sI, Ie) -> (Eθ, Eθ′)

The higher-order electrostatic terms of every pair of a set of like-sign ions of
charges `z`, at the ionic strength `Ie` (regularized as the model regularizes
it) with `sI = √Ie`: symmetric matrices, zero for two ions of the same charge.
"""
function _etheta_pairs(z, Aφ, sI, Ie)
    n = length(z)
    TT = promote_type(typeof(Aφ), typeof(sI), typeof(Ie))
    E = zeros(TT, n, n)
    Ep = zeros(TT, n, n)
    @inbounds for i in 1:n, j in (i + 1):n
        z[i] == z[j] && continue
        zz = z[i] * z[j]
        Jij, xJij = _pitzer_J(6 * zz * Aφ * sI)
        Jii, xJii = _pitzer_J(6 * z[i]^2 * Aφ * sI)
        Jjj, xJjj = _pitzer_J(6 * z[j]^2 * Aφ * sI)
        e = zz / (4 * Ie) * (Jij - Jii / 2 - Jjj / 2)
        ep = zz / (8 * Ie^2) * (xJij - xJii / 2 - xJjj / 2) - e / Ie
        E[i, j] = E[j, i] = e
        Ep[i, j] = Ep[j, i] = ep
    end
    return E, Ep
end

"""
    activity_model(cs::ChemicalSystem, model::PitzerActivityModel) -> Function

Return the closure `lna(n, p)` of the Pitzer model for `cs`.

The completeness of `model.parameters` is checked **here**, against the species
`cs` actually contains, and a missing cation-anion pair raises rather than
defaulting to ideal behavior, unless the model takes it as zero
(`missing_pairs = :zero`).
"""
function activity_model(cs::ChemicalSystem, model::PitzerActivityModel)
    idx_solvent = only(cs.idx_solvent)
    idx_solutes = cs.idx_solutes
    idx_gas = cs.idx_gas
    has_gas = !isempty(idx_gas)
    # The mixing model of the gas phase: ideal, or an equation of state.
    gas_mix = _gas_mixing(cs)
    # The solid solutions and the site families, prepared once; see `_MixingTerms`.
    mix = _MixingTerms(cs)

    M_w = ustrip(us"kg/mol", cs.species[idx_solvent][:M])
    n_sp = lastindex(cs.species)
    par = model.parameters

    sym = [symbol(sp) for sp in cs.species]
    zv = Int8[charge(sp) for sp in cs.species]
    cats = [i for i in idx_solutes if zv[i] > 0]
    ans = [i for i in idx_solutes if zv[i] < 0]
    neus = [i for i in idx_solutes if iszero(zv[i])]

    # ── The refusal that makes "no half-parameterized Pitzer" true ───────────
    missing_pairs = Tuple{String, String}[]
    for c in cats, a in ans
        haskey(par.beta0, (sym[c], sym[a])) || push!(missing_pairs, (sym[c], sym[a]))
    end
    if !isempty(missing_pairs) && model.missing_pairs === :refuse
        listed = join(("$(c)/$(a)" for (c, a) in missing_pairs), ", ")
        throw(
            ArgumentError(
                "PitzerActivityModel: the parameter set has no beta0 for " *
                    "$(length(missing_pairs)) cation-anion pair(s) this system contains: " *
                    "$listed. A Pitzer model cannot fall back on ideal behavior for a " *
                    "pair it does not describe, so the model refuses rather than " *
                    "returning a number. Either supply the missing parameters or build " *
                    "a species list the set covers. Note that a set fitted for a " *
                    "dissociated speciation does not describe ion pairs such as " *
                    "Ca(SO4)@ or CaOH+: those associations are already inside its beta " *
                    "coefficients, and carrying them as species counts them twice.",
            )
        )
    end

    # Precomputed per-pair tables, Float64 and not differentiated.
    npair = (length(cats), length(ans))
    B0 = [get(par.beta0, (sym[c], sym[a]), zero(valtype(par.beta0))) for c in cats, a in ans]
    B1 = [get(par.beta1, (sym[c], sym[a]), 0.0) for c in cats, a in ans]
    B2 = [get(par.beta2, (sym[c], sym[a]), 0.0) for c in cats, a in ans]
    CC = [
        get(par.Cphi, (sym[c], sym[a]), 0.0) /
            (2 * sqrt(abs(Int(zv[c]) * Int(zv[a])))) for c in cats, a in ans
    ]
    A1 = [
        (abs(zv[c]) >= 2 && abs(zv[a]) >= 2) ? par.alpha1_22 : par.alpha1
            for c in cats, a in ans
    ]
    ΘCC = [_sym2(par.theta, sym[i], sym[j]) for i in cats, j in cats]
    ΘAA = [_sym2(par.theta, sym[i], sym[j]) for i in ans, j in ans]
    ΨCCA = [_sym3(par.psi, sym[i], sym[j], sym[a]) for i in cats, j in cats, a in ans]
    ΨAAC = [_sym3(par.psi, sym[i], sym[j], sym[c]) for i in ans, j in ans, c in cats]
    ΛNC = [_sym2(par.lambda, sym[nn], sym[i]) for nn in neus, i in cats]
    ΛNA = [_sym2(par.lambda, sym[nn], sym[i]) for nn in neus, i in ans]

    # The temperature terms of the same entries, looked up as the values are.
    tt = par.temperature
    Z5 = ntuple(_ -> zero(eltype(B0)), 5)
    tget(kind, key) = haskey(tt, kind) ? get(tt[kind], key, nothing) : nothing
    tsym2(kind, i, j) = something(tget(kind, (i, j)), tget(kind, (j, i)), Z5)
    tsym3(kind, i, j, k) = something(
        (tget(kind, key) for key in ((i, j, k), (i, k, j), (j, i, k), (j, k, i), (k, i, j), (k, j, i)))...,
        Z5,
    )
    tables0 = (B0, B1, B2, CC, ΘCC, ΘAA, ΨCCA, ΨAAC, ΛNC, ΛNA)
    tablesT = (
        [something(tget(:beta0, (sym[c], sym[a])), Z5) for c in cats, a in ans],
        [something(tget(:beta1, (sym[c], sym[a])), Z5) for c in cats, a in ans],
        [something(tget(:beta2, (sym[c], sym[a])), Z5) for c in cats, a in ans],
        [
            something(tget(:Cphi, (sym[c], sym[a])), Z5) ./ (2 * sqrt(abs(Int(zv[c]) * Int(zv[a]))))
                for c in cats, a in ans
        ],
        [tsym2(:theta, sym[i], sym[j]) for i in cats, j in cats],
        [tsym2(:theta, sym[i], sym[j]) for i in ans, j in ans],
        [tsym3(:psi, sym[i], sym[j], sym[a]) for i in cats, j in cats, a in ans],
        [tsym3(:psi, sym[i], sym[j], sym[c]) for i in ans, j in ans, c in cats],
        [tsym2(:lambda, sym[nn], sym[i]) for nn in neus, i in cats],
        [tsym2(:lambda, sym[nn], sym[i]) for nn in neus, i in ans],
    )

    α2 = par.alpha2
    bp = par.b
    A_fixed = model.A                    # log10-basis Debye-Hückel A
    temp_dep = model.temperature_dependent
    t_terms = temp_dep && any(!isempty, values(tt))
    # The charges of the cations and of the anions, for the higher-order
    # electrostatic terms of their unlike pairs.
    zC = [Int(zv[c]) for c in cats]
    zA = [Int(zv[a]) for a in ans]
    etheta = model.etheta
    MT = promote_type(_captured_number_type(mix), _captured_number_type(model))

    function lna(n::AbstractVector, p)
        ϵ = p.ϵ
        _n = max.(n, _activity_floor(p))
        # The number type of everything the output is computed from: the
        # amounts, the state's parameters and the model's coefficients.
        TT = promote_type(eltype(_n), _number_type_of(p), MT)
        out = zeros(TT, n_sp)

        # The coefficients at the temperature of the state, when the set
        # carries temperature terms and the model is asked to use them.
        local B0, B1, B2, CC, ΘCC, ΘAA, ΨCCA, ΨAAC, ΛNC, ΛNA
        B0, B1, B2, CC, ΘCC, ΘAA, ΨCCA, ΨAAC, ΛNC, ΛNA =
            (t_terms && hasproperty(p, :T)) ? _pitzer_tables_at(p.T, tables0, tablesT) : tables0

        A_log10 = if temp_dep && hasproperty(p, :T) && hasproperty(p, :P)
            hkf_debye_huckel_params(p.T, p.P).A
        else
            A_fixed
        end
        Aφ = A_log10 * log(10) / 3       # the osmotic-coefficient basis

        kg = _n[idx_solvent] * M_w
        # `m` is built from `max.(n, ϵ)` and is therefore strictly positive:
        # the log below takes it bare. Adding ϵ a second time returns
        # `log(2ϵ)` at the floor, an `ln 2` offset that propagates into the
        # saturation index of every phase built on a floored primary.
        m = [_n[i] / kg for i in eachindex(_n)]

        I = zero(TT)
        Σm = zero(TT)
        Z = zero(TT)
        @inbounds for i in idx_solutes
            I += m[i] * zv[i]^2
            Σm += m[i]
            Z += m[i] * abs(zv[i])
        end
        I /= 2
        sI = sqrt(I + ϵ)

        # The higher-order electrostatic terms, `Φ = θ + Eθ` and `Φ′ = Eθ′`;
        # zero tables when they are left out.
        EθC, EθpC = etheta ? _etheta_pairs(zC, Aφ, sI, I + ϵ) : (zeros(TT, 0, 0), zeros(TT, 0, 0))
        EθA, EθpA = etheta ? _etheta_pairs(zA, Aφ, sI, I + ϵ) : (zeros(TT, 0, 0), zeros(TT, 0, 0))
        eθC(i, j) = etheta ? EθC[i, j] : zero(TT)
        eθA(i, j) = etheta ? EθA[i, j] : zero(TT)
        eθpC(i, j) = etheta ? EθpC[i, j] : zero(TT)
        eθpA(i, j) = etheta ? EθpA[i, j] : zero(TT)

        # f^γ and the pair sums that make up F
        fγ = -Aφ * (sI / (1 + bp * sI) + 2 / bp * log1p(bp * sI))
        F = fγ
        @inbounds for (ic, c) in enumerate(cats), (ia, a) in enumerate(ans)
            Bp = (
                B1[ic, ia] * _pitzer_gp(A1[ic, ia] * sI) +
                    B2[ic, ia] * _pitzer_gp(α2 * sI)
            ) / (I + ϵ)
            F += m[c] * m[a] * Bp
        end
        if etheta
            @inbounds for (ic, c) in enumerate(cats), (jc, c2) in enumerate(cats)
                jc > ic && (F += m[c] * m[c2] * eθpC(ic, jc))
            end
            @inbounds for (ia, a) in enumerate(ans), (ja, a2) in enumerate(ans)
                ja > ia && (F += m[a] * m[a2] * eθpA(ia, ja))
            end
        end

        # Σ_c Σ_a m_c m_a C_ca, shared by every ion
        ΣmmC = zero(TT)
        @inbounds for (ic, c) in enumerate(cats), (ia, a) in enumerate(ans)
            ΣmmC += m[c] * m[a] * CC[ic, ia]
        end

        # ── cations ──────────────────────────────────────────────────────────
        @inbounds for (ic, c) in enumerate(cats)
            s = zv[c]^2 * F + abs(zv[c]) * ΣmmC
            for (ia, a) in enumerate(ans)
                B = B0[ic, ia] + B1[ic, ia] * _pitzer_g(A1[ic, ia] * sI) +
                    B2[ic, ia] * _pitzer_g(α2 * sI)
                s += m[a] * (2 * B + Z * CC[ic, ia])
            end
            for (jc, c2) in enumerate(cats)
                c2 == c && continue
                acc = 2 * (ΘCC[ic, jc] + eθC(ic, jc))
                for (ia, a) in enumerate(ans)
                    acc += m[a] * ΨCCA[ic, jc, ia]
                end
                s += m[c2] * acc
            end
            for (ia, a) in enumerate(ans), (ja, a2) in enumerate(ans)
                ja > ia || continue
                s += m[a] * m[a2] * ΨAAC[ia, ja, ic]
            end
            for (inn, nn) in enumerate(neus)
                s += 2 * m[nn] * ΛNC[inn, ic]
            end
            out[c] = s + log(m[c])
        end

        # ── anions ───────────────────────────────────────────────────────────
        @inbounds for (ia, a) in enumerate(ans)
            s = zv[a]^2 * F + abs(zv[a]) * ΣmmC
            for (ic, c) in enumerate(cats)
                B = B0[ic, ia] + B1[ic, ia] * _pitzer_g(A1[ic, ia] * sI) +
                    B2[ic, ia] * _pitzer_g(α2 * sI)
                s += m[c] * (2 * B + Z * CC[ic, ia])
            end
            for (ja, a2) in enumerate(ans)
                a2 == a && continue
                acc = 2 * (ΘAA[ia, ja] + eθA(ia, ja))
                for (ic, c) in enumerate(cats)
                    acc += m[c] * ΨAAC[ia, ja, ic]
                end
                s += m[a2] * acc
            end
            for (ic, c) in enumerate(cats), (jc, c2) in enumerate(cats)
                jc > ic || continue
                s += m[c] * m[c2] * ΨCCA[ic, jc, ia]
            end
            for (inn, nn) in enumerate(neus)
                s += 2 * m[nn] * ΛNA[inn, ia]
            end
            out[a] = s + log(m[a])
        end

        # ── neutral solutes ──────────────────────────────────────────────────
        @inbounds for (inn, nn) in enumerate(neus)
            s = zero(TT)
            for (ic, c) in enumerate(cats)
                s += 2 * m[c] * ΛNC[inn, ic]
            end
            for (ia, a) in enumerate(ans)
                s += 2 * m[a] * ΛNA[inn, ia]
            end
            out[nn] = s + log(m[nn])
        end

        # ── the solvent, from the osmotic coefficient ────────────────────────
        acc = -Aφ * I^(3 // 2) / (1 + bp * sI)
        @inbounds for (ic, c) in enumerate(cats), (ia, a) in enumerate(ans)
            Bφ = B0[ic, ia] + B1[ic, ia] * exp(-A1[ic, ia] * sI) +
                B2[ic, ia] * exp(-α2 * sI)
            acc += m[c] * m[a] * (Bφ + Z * CC[ic, ia])
        end
        @inbounds for (ic, c) in enumerate(cats), (jc, c2) in enumerate(cats)
            jc > ic || continue
            t = ΘCC[ic, jc] + eθC(ic, jc) + I * eθpC(ic, jc)
            for (ia, a) in enumerate(ans)
                t += m[a] * ΨCCA[ic, jc, ia]
            end
            acc += m[c] * m[c2] * t
        end
        @inbounds for (ia, a) in enumerate(ans), (ja, a2) in enumerate(ans)
            ja > ia || continue
            t = ΘAA[ia, ja] + eθA(ia, ja) + I * eθpA(ia, ja)
            for (ic, c) in enumerate(cats)
                t += m[c] * ΨAAC[ia, ja, ic]
            end
            acc += m[a] * m[a2] * t
        end
        @inbounds for (inn, nn) in enumerate(neus)
            for (ic, c) in enumerate(cats)
                acc += m[nn] * m[c] * ΛNC[inn, ic]
            end
            for (ia, a) in enumerate(ans)
                acc += m[nn] * m[a] * ΛNA[inn, ia]
            end
        end
        φ = 1 + 2 * acc / (Σm + ϵ)
        out[idx_solvent] = -M_w * Σm * φ

        has_gas && _gas_lna!(out, _n, idx_gas, p, gas_mix)
        # Solid solutions and surface sites mix on budgets of their own; leaving
        # either out would give its members unit activity, silently.
        _mixing_lna!(out, _n, mix, p, ϵ)
        return out
    end

    return lna
end
