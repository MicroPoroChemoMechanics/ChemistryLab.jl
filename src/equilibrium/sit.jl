# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

using DynamicQuantities

# ── Specific ion Interaction Theory ───────────────────────────────────────────

"""
    struct SITCoefficient{T<:Real}

One specific-ion interaction coefficient, `ε`, together with **where it came
from**.

A number in a thermodynamic model is not self-describing: `ε(Na⁺, Cl⁻)` measured
on a sodium chloride solution and `ε(Al(OH)₄⁻, Na⁺)` borrowed from a sulfate
analog because nobody ever measured an aluminate one are the same Julia
`Float64` and are not the same claim. This package already makes that
distinction in data — `data/pitzer-reardon1990.toml` carries
`origin = "estimated:<analog>"` on every entry it had to borrow — and this type
carries it into the code, so it survives into a printed table or a fitted
result instead of stopping at the file.

# Fields

  - `value`: the coefficient, in kg/mol.
  - `origin`: a short free-text provenance, for example the database and version
    it was read from, `"estimated:<analog>"` for a borrowed value, or
    `"fitted:<dataset>"` for one this package identified.

See also: [`SITParameters`](@ref), [`sit_epsilon`](@ref).
"""
struct SITCoefficient{T <: Real}
    value::T
    origin::String
end

SITCoefficient(value::Real) = SITCoefficient(value, "unstated")

Base.show(io::IO, c::SITCoefficient) = print(io, c.value, " (", c.origin, ")")

"""
    struct SITParameters{T<:Real}

A compilation of [`SITCoefficient`](@ref)s, keyed by the unordered pair of
species symbols, with the source of the compilation as a whole.

SIT sums `ε` over ions of **opposite** charge only, so the pair is unordered and
the lookup is symmetric: `sit_epsilon(p, "Na+", "Cl-")` and
`sit_epsilon(p, "Cl-", "Na+")` are the same entry.

# Fields

  - `epsilon`: `(a, b) => SITCoefficient`, with `a ≤ b` as strings so one pair
    has one key.
  - `source`: what the compilation is, as a whole — a database name and version,
    a paper, or a note saying it was assembled here.
"""
struct SITParameters{T <: Real}
    epsilon::Dict{Tuple{String, String}, SITCoefficient{T}}
    source::String
end

"""
    SITParameters(pairs; source = "unstated") -> SITParameters

Build a compilation from `(a, b) => ε` pairs, where `ε` is a number or a
[`SITCoefficient`](@ref). Keys are normalized so that a pair given either way
round is one entry, and giving the same pair twice with different values is an
error rather than a last-one-wins.
"""
function SITParameters(pairs; source::AbstractString = "unstated")
    # The element type FOLLOWS the values rather than being pinned to `Float64`.
    # Pinning it is the eltype-promotion trap this package has been caught by
    # before, and here it would be fatal to the point of the type: a coefficient
    # that cannot be a `ForwardDiff.Dual` cannot be IDENTIFIED from data, only
    # finite-differenced.
    collected = [
        (
            k,
            v isa SITCoefficient ? v : SITCoefficient(v, source),
        ) for (k, v) in pairs
    ]
    T = isempty(collected) ? Float64 :
        promote_type((typeof(c.value) for (_, c) in collected)...)
    out = Dict{Tuple{String, String}, SITCoefficient{T}}()
    for (k, v) in collected
        a, b = String(first(k)), String(last(k))
        key = a <= b ? (a, b) : (b, a)
        c = SITCoefficient(convert(T, v.value), v.origin)
        if haskey(out, key) && out[key].value != c.value
            throw(
                ArgumentError(
                    "two different ε for the pair ($a, $b): $(out[key].value) and " *
                        "$(c.value). SIT's ε is symmetric in the pair, so these are " *
                        "the same coefficient given twice."
                ),
            )
        end
        out[key] = c
    end
    return SITParameters(out, String(source))
end

SITParameters() = SITParameters(Pair{Tuple{String, String}, Float64}[]; source = "empty")

Base.length(p::SITParameters) = length(p.epsilon)
Base.isempty(p::SITParameters) = isempty(p.epsilon)

function Base.show(io::IO, p::SITParameters)
    return print(io, "SITParameters(", length(p.epsilon), " coefficients, ", p.source, ")")
end

"""
    sit_epsilon(p::SITParameters, a, b) -> SITCoefficient or nothing

The coefficient for the unordered pair `(a, b)`, or `nothing` when the
compilation does not carry it.

`nothing` and zero are different answers, and the distinction is the point: the
SIT literature's convention is to take an unlisted `ε` as zero, which is a
statement about what is customary and not about what was measured. See
[`missing_epsilon_pairs`](@ref) for asking a system which of its pairs are
resting on that convention.
"""
function sit_epsilon(p::SITParameters, a::AbstractString, b::AbstractString)
    key = a <= b ? (String(a), String(b)) : (String(b), String(a))
    return get(p.epsilon, key, nothing)
end

"""
    struct SITActivityModel{T<:Real, P} <: AbstractActivityModel

The Specific ion Interaction Theory of Brønsted, Guggenheim and Scatchard, in
the form the NEA thermodynamic reviews use.

```math
\\log_{10}\\gamma_i \\;=\\; -z_i^2\\,\\frac{A\\sqrt{I}}{1 + b\\sqrt{I}}
   \\;+\\; \\sum_k \\varepsilon(i,k)\\, m_k
```

with the sum over ions `k` of charge **opposite** to `i`, `m_k` in molality and
`b = 1.5 kg^{1/2} mol^{-1/2}`.

# What it is, against its neighbors in this package

Read the three together and the family is one idea with three answers to the
same question — what to do beyond Debye-Hückel's limiting law.

| | ion size | composition-dependent term |
|:--|:--|:--|
| [`DaviesActivityModel`](@ref) | none, absorbed into `1 + √I` | `−A z² b I`, one parameter for every ion |
| [`HKFActivityModel`](@ref) | per ion, `å` | `Ḃ I`, one parameter for every ion |
| **SIT** | none, absorbed into `1 + 1.5√I` | `Σ ε(i,k) m_k`, one parameter **per pair** |
| [`PitzerActivityModel`](@ref) | none | binary *and* ternary, several per pair |

So SIT sits between Davies and Pitzer exactly where its parameter count does. It
is the model the NEA reviews and the ANDRA/ThermoChimie database are calibrated
in, which is why a log K taken from those and used under Davies is not the
constant that was fitted.

# Validity, stated because it is narrower than the equation looks

SIT is meant for ionic strengths up to roughly 3 to 4 mol/kg, and it is a
*truncation*: the pairwise term is the first order of a virial expansion, so
above that range the terms it drops stop being small. It also says nothing about
ion pairs that are better described as species — where a complex forms, SIT
expects it in the speciation and not in `ε`.

# Fields

  - `A`: the Debye-Hückel slope, `0.509` at 25 °C, or recomputed per temperature
    when `temperature_dependent`.
  - `b`: the denominator coefficient, `1.5` by convention. It is a *convention*,
    not a fit: changing it makes every published `ε` inconsistent with the model
    that produced it, so it is a field only to make that visible.
  - `parameters`: a [`SITParameters`](@ref) compilation. **Empty by default**,
    which reduces the model to its Debye-Hückel term — a legitimate limiting
    case, and never silently a full SIT calculation.
  - `temperature_dependent`: whether `A` follows temperature.

# Water

The solvent activity uses the same mole-fraction approximation as
[`DaviesActivityModel`](@ref) rather than an osmotic coefficient. That is a
known departure from a full SIT treatment and is shared with this package's
other models; it does not enter a comparison made at prescribed proton
activity.

See also: [`SITParameters`](@ref), [`build_sit_parameters`](@ref),
[`missing_epsilon_pairs`](@ref).
"""
struct SITActivityModel{T <: Real, P} <: AbstractActivityModel
    A::T
    b::T
    parameters::P
    temperature_dependent::Bool
end

"""
    SITActivityModel(; A = 0.509, b = 1.5, parameters = SITParameters(),
                       temperature_dependent = false) -> SITActivityModel

Build a [`SITActivityModel`](@ref).
"""
function SITActivityModel(;
        A::Real = 0.509, b::Real = 1.5,
        parameters = SITParameters(),
        temperature_dependent::Bool = false,
    )
    b > 0 || throw(ArgumentError("the SIT denominator coefficient must be positive; got $b."))
    v = promote(float(A), float(b))
    return SITActivityModel{eltype(v), typeof(parameters)}(
        v[1], v[2], parameters, temperature_dependent,
    )
end

function Base.show(io::IO, m::SITActivityModel)
    return print(
        io, "SITActivityModel(A = ", m.A, ", b = ", m.b, ", ", m.parameters, ")",
    )
end

"""
    missing_epsilon_pairs(cs::ChemicalSystem, model::SITActivityModel) -> Vector

Every oppositely charged pair of aqueous solutes in `cs` for which `model` has
no `ε`, as `(a, b)` symbol tuples.

The SIT literature takes an unlisted coefficient as zero. That convention is
usually harmless and is never a measurement, so this reports what a calculation
is resting on rather than leaving it to be assumed. An empty result means every
pair the system can form is carried by the compilation.
"""
function missing_epsilon_pairs(cs::ChemicalSystem, model::SITActivityModel)
    syms = symbol.(cs.species)
    z = Int[charge(sp) for sp in cs.species]
    out = Tuple{String, String}[]
    for i in cs.idx_solutes, j in cs.idx_solutes
        i < j || continue
        (iszero(z[i]) || iszero(z[j])) && continue
        sign(z[i]) == sign(z[j]) && continue
        sit_epsilon(model.parameters, syms[i], syms[j]) === nothing &&
            push!(out, (syms[i], syms[j]))
    end
    return out
end

_sit_eltype(::SITParameters{T}) where {T} = T

"""
    _sit_epsilon_matrix(cs, model) -> Matrix

The `ε` a system needs, as a dense matrix over its species indices: zero for a
like-charged pair, for a neutral, and for a pair the compilation does not carry.

Its element type is the compilation's, so a `Dual` coefficient stays a `Dual`
all the way to the activity — which is what differentiating with respect to an
interaction coefficient requires.

Built once when the closure is, so the inner loop is a matrix read rather than a
dictionary lookup.
"""
function _sit_epsilon_matrix(cs::ChemicalSystem, model::SITActivityModel)
    n = length(cs.species)
    E = zeros(_sit_eltype(model.parameters), n, n)
    syms = symbol.(cs.species)
    z = Int[charge(sp) for sp in cs.species]
    for i in cs.idx_solutes, j in cs.idx_solutes
        (iszero(z[i]) || iszero(z[j])) && continue
        sign(z[i]) == sign(z[j]) && continue
        c = sit_epsilon(model.parameters, syms[i], syms[j])
        c === nothing && continue
        E[i, j] = c.value
    end
    return E
end

"""
    activity_model(cs::ChemicalSystem, model::SITActivityModel) -> Function

A closure `lna(n, p) -> Vector` of log-activities under SIT.

Unlike [`DaviesActivityModel`](@ref) and [`HKFActivityModel`](@ref), SIT has no
`_log10γ_ion`: its correction depends on the **whole composition** through
`Σ ε(i,k) m_k`, not on the ion's charge and the ionic strength alone, so there
is no per-ion formula to share with `activity_coefficients`. That is the same
reason [`PitzerActivityModel`](@ref) has none, and it is a property of the
model rather than an omission.

AD-compatible: every computation inside accepts `ForwardDiff.Dual`.
"""
function activity_model(cs::ChemicalSystem, model::SITActivityModel)

    idx_solvent = only(cs.idx_solvent)
    idx_solutes = cs.idx_solutes
    idx_gas = cs.idx_gas

    ss_groups = cs.ss_groups
    has_ss = !isempty(ss_groups)
    has_gas = !isempty(idx_gas)
    ss_models = has_ss ? map(ss -> ss.model, cs.solid_solutions) : nothing

    site_groups = cs.site_groups
    has_sites = !isempty(site_groups)
    site_models = has_sites ? map(f -> f.model, cs.site_families) : nothing
    site_denticity = has_sites ?
        [Int[denticity(f, sp) for sp in site_members(f)] for f in cs.site_families] :
        nothing
    site_charges = has_sites ?
        [Float64[charge(sp) for sp in site_members(f)] for f in cs.site_families] :
        nothing

    site_needs_I = has_sites && any(needs_ionic_strength, site_models)
    site_solvent = isempty(cs.idx_solvent) ? 0 : only(cs.idx_solvent)
    site_ions = site_needs_I ?
        [i for i in cs.idx_solutes if !iszero(charge(cs.species[i]))] : Int[]
    site_ion_z = Float64[charge(cs.species[i]) for i in site_ions]
    site_Mw = (site_needs_I && !iszero(site_solvent)) ?
        ustrip(us"kg/mol", cs.species[site_solvent][:M]) : 1.0
    site_support_idx, site_support_z = has_sites ? _support_members(cs) : (nothing, nothing)

    M_w = ustrip(us"kg/mol", cs.species[idx_solvent][:M])

    A_fixed = model.A
    b_sit = model.b
    temp_dep = model.temperature_dependent

    zv = Int8[charge(sp) for sp in cs.species]
    n_sp = lastindex(zv)
    idx_ions = [i for i in idx_solutes if !iszero(zv[i])]
    idx_neutrals = [i for i in idx_solutes if iszero(zv[i])]

    # The pairwise coefficients, resolved once. A pair the compilation does not
    # carry enters as zero, which is the SIT literature's convention;
    # `missing_epsilon_pairs` is how a caller finds out which those are.
    E = _sit_epsilon_matrix(cs, model)
    ln10 = log(10.0)

    function lna(n::AbstractVector, p)
        ϵ = p.ϵ
        _n = max.(n, ϵ)

        A = if temp_dep && hasproperty(p, :T) && hasproperty(p, :P)
            hkf_debye_huckel_params(p.T, p.P).A
        else
            A_fixed
        end

        # THE PROMOTION INCLUDES THE PARAMETERS, not only the composition.
        # Differentiating with respect to an interaction coefficient makes `ε` a
        # `Dual` while `n` stays a `Float64`, and an output vector typed on `n`
        # alone would refuse to hold the result — the trap that lets every piece
        # pass its own test and the chain break.
        ET = promote_type(eltype(_n), eltype(E), typeof(A))
        out = zeros(ET, n_sp)

        n_w = _n[idx_solvent]
        denom_mol = n_w * M_w

        I = zero(eltype(_n))
        @inbounds for i in idx_solutes
            I = I + (_n[i] / denom_mol) * zv[i]^2
        end
        I = I / 2
        sqrtI = sqrt(I + ϵ)
        D = A * sqrtI / (1 + b_sit * sqrtI)

        @inbounds for i in idx_ions
            pair = zero(ET)
            for k in idx_ions
                e = E[i, k]
                iszero(e) && continue
                pair += e * (_n[k] / denom_mol)
            end
            mᵢ = _n[i] / denom_mol
            out[i] = ln10 * (-zv[i]^2 * D + pair) + log(mᵢ + ϵ)
        end

        # A neutral solute has no Debye-Hückel term and, in SIT, an interaction
        # with ions only through an `ε` that is almost never tabulated — so this
        # is ideal unless the compilation says otherwise.
        @inbounds for i in idx_neutrals
            pair = zero(ET)
            for k in idx_ions
                e = E[i, k]
                iszero(e) && continue
                pair += e * (_n[k] / denom_mol)
            end
            mᵢ = _n[i] / denom_mol
            out[i] = ln10 * pair + log(mᵢ + ϵ)
        end

        n_aqueous = n_w + sum((_n[i] for i in idx_solutes); init = zero(eltype(_n)))
        out[idx_solvent] = log(n_w / n_aqueous)

        if has_gas
            n_gas = sum((_n[i] for i in idx_gas); init = zero(eltype(_n)))
            @inbounds for i in idx_gas
                out[i] = log(_n[i] / n_gas)
            end
        end

        if has_ss
            T_val = hasproperty(p, :T) ? p.T : 298.15
            _solid_solution_lna!(out, _n, ss_groups, ss_models, T_val, ϵ)
        end

        if has_sites
            T_val = hasproperty(p, :T) ? p.T : 298.15
            I_site = site_needs_I ?
                _aqueous_ionic_strength(_n, site_ions, site_ion_z, site_solvent, site_Mw) :
                zero(eltype(_n))
            ψ_site = hasproperty(p, :ψ_site) ? p.ψ_site : nothing
            _site_mixing_lna!(
                out, _n, site_groups, site_models, site_denticity, site_charges,
                I_site, T_val, ϵ, ψ_site, site_support_idx, site_support_z
            )
        end

        return out
    end

    return lna
end
