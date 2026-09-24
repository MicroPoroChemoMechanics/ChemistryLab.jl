# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)
# Portions of this file (the `transition_state` and `first_order_rate` rate
# factories implementing the Palandri-Kharaka / transition-state theory model
# of mineral dissolution-precipitation kinetics) are Julia ports adapted from
# the Reaktoro C++ library (https://github.com/reaktoro/reaktoro),
# Copyright © 2014-2024 Allan Leal, distributed under the LGPL-2.1-or-later.

using DynamicQuantities

# ── KineticReaction ───────────────────────────────────────────────────────────

"""
    struct KineticReaction{R<:AbstractReaction, F, H}

Associates a chemical [`Reaction`](@ref) with a compiled [`KineticFunc`](@ref).

Following Leal et al. (2017), **reactions** — not individual species — carry kinetics.
A single mineral can therefore appear as a reactant in multiple `KineticReaction` objects
(e.g. C₃A → ettringite and C₃A → monosulphate for multi-pathway cement hydration).
The ODE state is indexed by unique mineral species, and contributions from all reactions
that consume the same mineral are accumulated.

# Fields

  - `reaction`: the underlying [`Reaction`](@ref) / [`CemReaction`](@ref).
  - `rate_fn`: a [`KineticFunc`](@ref) (or any callable matching the six-argument signature
    `(T, P, t, n, lna, n_initial) -> Real`) computing r [mol/s].
  - `idx_mineral`: index of the primary (controlling) mineral species in the parent
    `ChemicalSystem`. Determined automatically as the first solid (AS_CRYSTAL) reactant.
  - `stoich`: stoichiometric coefficient vector for all species in the system.
    Sign convention: positive for products, negative for reactants.
  - `heat_per_mol`: enthalpy of reaction [J/mol], positive = exothermic (heat released).
    When `nothing` (default), the enthalpy is derived from the stoichiometric sum
    of species `:ΔₐH⁰` values.

# Constructors

**From a species name** (convenience, builds a minimal dissolution Reaction):

```julia
pk = parrot_killoh(PK_PARAMS_C3S, "C3S")
kr = KineticReaction(cs, "C3S", pk)
kr = KineticReaction(cs, "C3S", pk; heat_per_mol = 114_634.0)
```

**From an explicit Reaction** (multi-pathway):

```julia
pk_c3a = parrot_killoh(PK_PARAMS_C3A, "C3A")
kr_ett  = KineticReaction(cs, rxn_C3A_ettringite,   pk_c3a)
kr_mono = KineticReaction(cs, rxn_C3A_monosulphate, pk_c3a)
```

**Reaction-centric** (rate stored in `rxn.properties[:rate]`):

```julia
rxn[:rate] = parrot_killoh(PK_PARAMS_C3S, "C3S")
kr = KineticReaction(cs, rxn)
```
"""
struct KineticReaction{R <: AbstractReaction, F, H}
    reaction::R
    rate_fn::F           # KineticFunc or compatible callable
    idx_mineral::Int
    stoich::Vector{Float64}    # stoich coefficients for all species in system
    heat_per_mol::H            # Nothing or Float64: enthalpy [J/mol], positive = exothermic

    function KineticReaction{R, F, H}(
            reaction::R,
            rate_fn::F,
            idx_mineral::Int,
            stoich::Vector{Float64},
            heat_per_mol::H,
        ) where {R <: AbstractReaction, F, H}
        idx_mineral > 0 || throw(ArgumentError("idx_mineral must be a positive integer"))
        isempty(stoich) && throw(ArgumentError("stoich cannot be empty"))
        return new{R, F, H}(reaction, rate_fn, idx_mineral, stoich, heat_per_mol)
    end
end

"""
    KineticReaction(rxn, rate_fn, idx_mineral, stoich; heat_per_mol=nothing)

Low-level constructor: explicit `Reaction`, rate callable, index, and stoichiometry.
"""
function KineticReaction(
        rxn::R,
        rate_fn::F,
        idx_mineral::Integer,
        stoich::AbstractVector{<:Real};
        heat_per_mol = nothing,
    ) where {R <: AbstractReaction, F}
    hpm = _strip_heat_per_mol(heat_per_mol)
    return KineticReaction{R, F, typeof(hpm)}(
        rxn, rate_fn, Int(idx_mineral), Float64.(stoich), hpm,
    )
end

"""
    KineticReaction(cs::ChemicalSystem, species_name::AbstractString, rate_fn;
                    heat_per_mol=nothing) -> KineticReaction

Convenience constructor: look up `species_name` in `cs` and build a minimal dissolution
[`Reaction`](@ref) (species as sole reactant, no products) automatically.

The default stoichiometry places `-1.0` at the mineral index and `0.0` everywhere else.
"""
function KineticReaction(
        cs::ChemicalSystem,
        species_name::AbstractString,
        rate_fn;
        heat_per_mol = nothing,
    )
    idx = findfirst(
        sp -> phreeqc(formula(sp)) == species_name || string(symbol(sp)) == species_name,
        cs.species,
    )
    isnothing(idx) && throw(
        ArgumentError(
            "Species \"$species_name\" not found in ChemicalSystem. " *
                "Use phreeqc(formula(sp)) or symbol(sp) to check species names.",
        ),
    )

    sp = cs.species[idx]
    n_sp = length(cs.species)
    rxn = Reaction(
        OrderedDict(sp => 1),
        OrderedDict{typeof(sp), Int}();
        symbol = string(symbol(sp)),
        equal_sign = '→',
    )

    s = zeros(Float64, n_sp)
    s[idx] = -1.0
    hpm = _strip_heat_per_mol(heat_per_mol)
    return KineticReaction{typeof(rxn), typeof(rate_fn), typeof(hpm)}(
        rxn, rate_fn, Int(idx), s, hpm,
    )
end

"""
    KineticReaction(cs::ChemicalSystem, rxn::AbstractReaction, rate_fn;
                    heat_per_mol=nothing) -> KineticReaction

Construct a `KineticReaction` from an explicit [`Reaction`](@ref) object.

The controlling mineral index (`idx_mineral`) is determined automatically as the index
of the **first solid (AS_CRYSTAL) reactant** found in `rxn.reactants` that is present in
`cs`. The stoichiometric vector is derived from the reaction stoichiometry.

This constructor is the recommended entry point for multi-pathway kinetics:

```julia
pk_c3a = parrot_killoh(PK_PARAMS_C3A, "C3A")
kr_ett  = KineticReaction(cs, cs.dict_reactions["C3A_ettringite"],   pk_c3a)
kr_mono = KineticReaction(cs, cs.dict_reactions["C3A_monosulphate"], pk_c3a)
kp = KineticsProblem(cs, [kr_C3S, kr_ett, kr_mono], state0, tspan)
```
"""
function KineticReaction(
        cs::ChemicalSystem,
        rxn::R,
        rate_fn;
        heat_per_mol = nothing,
    ) where {R <: AbstractReaction}
    idx = _find_mineral_idx(cs, rxn)
    isnothing(idx) && throw(
        ArgumentError(
            "No reactant species of the given reaction found in the ChemicalSystem.",
        ),
    )
    stoich_vec = _normalise_to_mineral!(_stoich_from_reaction(cs, rxn), idx)
    hpm = _strip_heat_per_mol(heat_per_mol)
    return KineticReaction{R, typeof(rate_fn), typeof(hpm)}(
        rxn, rate_fn, Int(idx), stoich_vec, hpm,
    )
end

"""
    KineticReaction(cs::ChemicalSystem, rxn::AbstractReaction) -> KineticReaction

Reaction-centric constructor: build a `KineticReaction` from a [`Reaction`](@ref)
that carries its kinetics in `reaction.properties`.

Required property:
  - `rxn[:rate]` — a [`KineticFunc`](@ref) **or** any callable matching
    `(T, P, t, n, lna, n_initial) -> Real`. Non-`KineticFunc` callables are
    wrapped automatically in a `KineticFunc` with empty `refs`.

Optional property:
  - `rxn[:heat_per_mol]` — a `Number` giving the molar enthalpy [J/mol] for calorimetry.

# Examples

```julia
pk = parrot_killoh(PK_PARAMS_C3S, "C3S")
rxn[:rate]         = pk
rxn[:heat_per_mol] = 114_634.0
kr = KineticReaction(cs, rxn)

# Build problem directly from a list of annotated Reaction objects:
kp = KineticsProblem(cs, [rxn_C3S, rxn_C3A, rxn_C2S], state0, tspan)
```
"""
function KineticReaction(cs::ChemicalSystem, rxn::AbstractReaction)
    haskey(properties(rxn), :rate) || throw(
        ArgumentError(
            "Reaction \"$(rxn.symbol)\" must have a :rate entry in its properties. " *
                "Attach a KineticFunc via rxn[:rate] = parrot_killoh(...).",
        ),
    )

    rate_raw = properties(rxn)[:rate]
    rate_fn = rate_raw isa KineticFunc ? rate_raw :
        KineticFunc(rate_raw, NamedTuple(), u"mol/s")

    heat_raw = get(properties(rxn), :heat_per_mol, nothing)
    heat_val = _strip_heat_per_mol(heat_raw)

    idx = _find_mineral_idx(cs, rxn)
    isnothing(idx) && throw(
        ArgumentError(
            "No reactant species of reaction \"$(rxn.symbol)\" found in the ChemicalSystem.",
        ),
    )

    stoich_vec = _normalise_to_mineral!(_stoich_from_reaction(cs, rxn), idx)
    return KineticReaction{typeof(rxn), typeof(rate_fn), typeof(heat_val)}(
        rxn, rate_fn, Int(idx), stoich_vec, heat_val,
    )
end

# ── transition_state factory ──────────────────────────────────────────────────

"""
    transition_state(mechanisms, cs, rxn, surface; ϵ=1e-16) -> KineticFunc

Build a Transition-State Theory (TST) dissolution/precipitation rate function from a
list of [`RateMechanism`](@ref) objects, returning a [`KineticFunc`](@ref).

The compiled closure captures:
  - the host name and its molar mass, from the [`SurfaceSupport`](@ref) when one is
    given, otherwise rediscovered from `rxn` + `cs`
  - the area model, evaluated at **every** step from the current and initial
    amounts, so an area that follows the microstructure needs no change here
  - stoichiometry and `ΔₐG⁰` callables for all aqueous species (T-dependent Ω)

The net rate [mol/s] is:

```
r = A(n) × Σ_m [ k_m(T) × Π_cat(aᵢ^nᵢ) × (1 - Ω^p) × |1 - Ω^p|^(q-1) ]
```

where `Ω(T) = exp(Σ νᵢ ln aᵢ + Σ νᵢ ΔₐG°ᵢ(T)/(RT))` is re-evaluated at every ODE
step — correct for variable-temperature semi-adiabatic calorimetry.

# Arguments

  - `mechanisms`: vector of [`RateMechanism`](@ref) (acid, neutral, base, …).
  - `cs`: [`ChemicalSystem`](@ref) supplying `ΔₐG⁰` callables for aqueous species.
  - `rxn`: `AbstractReaction` defining stoichiometry and the mineral species.
  - `surface`: a [`SurfaceSupport`](@ref), which names the host solid and carries its
    area model, or an [`AbstractSurfaceModel`](@ref) alone, in which case the
    host is the first solid reactant of `rxn`.
  - `ϵ`: regularization floor near Ω = 1 (default `1e-16`).

!!! note "The molar mass is no longer guessed"
    A host species without an `:M` property used to fall back to 0.1 kg/mol,
    silently, which is wrong by up to an order of magnitude and scales the whole
    rate. It now raises, naming the species.

# Returns

A [`KineticFunc`](@ref) callable as
`f(T, P, t, n::StateView, lna::StateView, n_initial::StateView) -> Real [mol/s]`.

AD-compatible: all operations use generic Julia arithmetic; no `Float64` casts.

# References

  - Palandri, J.L. & Kharaka, Y.K. (2004). USGS Open-File Report 2004-1068.
  - Leal, A.M.M. et al. (2017). Pure Appl. Chem. 89, 597–643.
"""
function transition_state(
        mechanisms::AbstractVector{<:RateMechanism},
        cs::ChemicalSystem,
        rxn::AbstractReaction,
        surface::Union{SurfaceSupport, AbstractSurfaceModel};
        ϵ::Real = 1.0e-16,
    )
    mineral_name, M, area_model = _surface_context(cs, rxn, surface)
    stoich_species = _stoich_named(cs, rxn)   # Vector of (name, ν, ΔG°_fn)

    f = (T, _P, _t, n, lna, n_initial) -> begin
        n_m = max(n[mineral_name], oneunit(T) * 1.0e-30)
        n_m0 = max(n_initial[mineral_name], oneunit(T) * 1.0e-30)
        A = total_area(area_model, n_m, n_m0, M)
        ln_iap = sum(ν * lna[sp] for (sp, ν, _) in stoich_species)
        ln_K = -sum(ν * ΔG_fn(; T = T, unit = false) / (R_GAS * T) for (_, ν, ΔG_fn) in stoich_species)
        Ω = exp(ln_iap - ln_K)
        r = zero(promote_type(typeof(T), typeof(Ω), typeof(A)))
        for mech in mechanisms
            k_val = mech.k(; T = T)
            cat_term = one(r)
            for cat in mech.catalysts
                if haskey(lna, cat.species)
                    cat_term *= exp(cat.n * lna[cat.species])
                end
            end
            Ωp = Ω^mech.p
            diff = one(r) - Ωp
            sat = diff * (diff^2 + ϵ)^((mech.q - one(r)) / 2)
            r = r + k_val * cat_term * sat
        end
        return A * r
    end

    refs = (T = 298.15u"K", P = 1.0e5u"Pa")
    return KineticFunc(f, refs, u"mol/s")
end

# ── first_order_rate factory ──────────────────────────────────────────────────

"""
    first_order_rate(k, cs, rxn, surface_model; p=1.0, q=1.0, ϵ=1e-16) -> KineticFunc

Build a single-mechanism first-order TST rate as a [`KineticFunc`](@ref).

```
r = A(n) × k(T) × sign(1 - Ω) × |1 - Ω^p|^q
```

This is a convenience wrapper around [`transition_state`](@ref) with one no-catalyst
mechanism. Useful as a minimal test case or for empirical fits.

# Arguments

  - `k`: rate constant as an [`AbstractFunc`](@ref) (e.g. from
    [`arrhenius_rate_constant`](@ref)).
  - `cs`, `rxn`, `surface`: same as [`transition_state`](@ref).
  - `p`, `q`: saturation exponents (defaults `1.0`).
  - `ϵ`: regularization floor (default `1e-16`).

# Examples

```julia
k = arrhenius_rate_constant(1e-7, 40000.0)
rf = first_order_rate(k, cs, rxn, SurfaceSupport("calcite", "Cal", BETSurfaceArea(90.0)))
kr = KineticReaction(cs, rxn, rf)
```
"""
function first_order_rate(
        k::AbstractFunc,
        cs::ChemicalSystem,
        rxn::AbstractReaction,
        surface::Union{SurfaceSupport, AbstractSurfaceModel};
        p::Real = 1.0,
        q::Real = 1.0,
        ϵ::Real = 1.0e-16,
    )
    T_p = typeof(promote(p, q)[1])
    mech = RateMechanism{typeof(k), T_p}(k, T_p(p), T_p(q), RateModelCatalyst{T_p}[])
    return transition_state([mech], cs, rxn, surface; ϵ = ϵ)
end

# ── Internal helpers ──────────────────────────────────────────────────────────

# Convert heat_per_mol to Float64 SI [J/mol], or return nothing.
_strip_heat_per_mol(::Nothing) = nothing
_strip_heat_per_mol(h) = Float64(safe_ustrip(us"J/mol", h))

# Find the index in cs.species of the first solid (AS_CRYSTAL) reactant of rxn.
# Falls back to the first reactant present in cs if no crystal phase is found.
function _find_mineral_idx(cs::ChemicalSystem, rxn::AbstractReaction)
    # A species the system itself declares kinetic wins over any guess. The
    # fallbacks below infer the controlling mineral from the reaction alone,
    # and on a reaction generated by the `kinetic_species` API — where the
    # kinetic mineral may sit on either side — that lands on the first reactant
    # found, which is typically the solvent. The ODE then integrates water as
    # if it were the dissolving phase.
    if !isempty(cs.idx_kinetic)
        for (sp, _) in Iterators.flatten((rxn.reactants, rxn.products))
            i = findfirst(s -> s == sp, cs.species)
            !isnothing(i) && i in cs.idx_kinetic && return i
        end
    end
    for (sp, _) in rxn.reactants
        i = findfirst(s -> s == sp, cs.species)
        !isnothing(i) && aggregate_state(cs.species[i]) == AS_CRYSTAL && return i
    end
    for (sp, _) in rxn.reactants
        i = findfirst(s -> s == sp, cs.species)
        isnothing(i) || return i
    end
    return nothing
end

# Orient and scale a stoichiometry vector to Leal's convention for the kinetic
# partition: the controlling mineral carries ν = -1, so a positive rate is a
# dissolution and `dn_k/dt = νₖᵀ r` decreases it.
#
# Both halves matter. A reaction generated by the `kinetic_species` API comes out
# of a nullspace diagonalization with an arbitrary orientation, and when it lands
# with the mineral on the product side the ODE grows the clinker instead of
# consuming it. The scaling is the "corrects by 1/|νₖ|" the `ChemicalSystem`
# docstring promises. A hand-written dissolution reaction already satisfies both
# and is left untouched.
function _normalise_to_mineral!(stoich::AbstractVector{<:Real}, idx::Integer)
    ν = stoich[idx]
    iszero(ν) && throw(
        ArgumentError(
            "the controlling species has a zero stoichiometric coefficient; " *
                "the reaction cannot define its kinetics"
        )
    )
    stoich ./= -ν            # ν < 0 → scale only; ν > 0 → scale and flip
    return stoich
end

# Build the stoichiometric coefficient vector (length = length(cs.species)) from a
# reaction. Reactants get negative coefficients, products get positive ones.
function _stoich_from_reaction(cs::ChemicalSystem, rxn::AbstractReaction)
    s = zeros(Float64, length(cs.species))
    for (sp, ν) in rxn.reactants
        i = findfirst(s_ -> s_ == sp, cs.species)
        isnothing(i) || (s[i] -= Float64(ν))
    end
    for (sp, ν) in rxn.products
        i = findfirst(s_ -> s_ == sp, cs.species)
        isnothing(i) || (s[i] += Float64(ν))
    end
    return s
end

# Molar mass of a species in kg/mol, refusing to invent one.
#
# The previous default of 0.1 kg/mol was silent and multiplicative: it scales the
# reactive area, hence the whole rate. A species carrying no `:M` is a data
# problem, and naming it is the only useful answer.
function _molar_mass_si(sp::AbstractSpecies)
    haskey(properties(sp), :M) || throw(
        ArgumentError(
            "species \"$(symbol(sp))\" carries no molar mass `:M`, which the " *
                "reactive area needs. Supply it on the species, or use a " *
                "`FixedSurfaceArea`, whose area does not depend on a mass.",
        )
    )
    return Float64(ustrip(us"kg/mol", sp[:M]))
end

# The key a rate law uses to find its host in a `StateView`.
#
# THE SYMBOL, NOT THE FORMULA. `build_kinetics_params` registers both — the
# formula first, then the symbol — so a formula shared by two polymorphs is
# overwritten and resolves to whichever was declared last. Measured on calcite
# and aragonite, both `CaCO3`, holding 1 and 100 mol: asking for `Cal` came back
# with Arg's amount and a reactive area a hundred times too large. The host had
# been resolved correctly and the resolution was then thrown away.
#
# Nothing changes for the ordinary species: one built without an explicit symbol
# takes its formula as its symbol, so this returns the same string it always
# did. It differs exactly where the two differ, which is the case that was wrong.
#
# A species with no symbol at all falls back to its formula, and that formula is
# then required to be unique — a rate that silently follows another phase is
# worse than one that refuses to be built.
function _rate_lookup_key(cs::ChemicalSystem, sp::AbstractSpecies)
    sym = symbol(sp)
    isempty(sym) || return sym
    key = phreeqc(formula(sp))
    n = count(s -> phreeqc(formula(s)) == key, cs.species)
    n == 1 || throw(
        ArgumentError(
            "species with formula \"$key\" carries no symbol, and $n species of " *
                "this system share that formula, so it cannot name one of them. " *
                "Give the species a symbol — that is what distinguishes two " *
                "polymorphs.",
        )
    )
    return key
end

# Returns (mineral_name::String, M::Float64) for the controlling mineral in rxn.
function _mineral_name_and_mass(cs::ChemicalSystem, rxn::AbstractReaction)
    idx = _find_mineral_idx(cs, rxn)
    isnothing(idx) && throw(
        ArgumentError("No mineral reactant found in reaction \"$(rxn.symbol)\"."),
    )
    sp = cs.species[idx]
    return _rate_lookup_key(cs, sp), _molar_mass_si(sp)
end

# Returns (host_name::String, M::Float64, area_model) for a rate factory.
#
# Two entry points, one contract: a bare area model keeps the historical
# behavior of rediscovering the host from the reaction, while a `SurfaceSupport` names
# it once and is looked up by symbol.
_surface_context(cs::ChemicalSystem, rxn::AbstractReaction, m::AbstractSurfaceModel) =
    (_mineral_name_and_mass(cs, rxn)..., m)

function _surface_context(cs::ChemicalSystem, rxn::AbstractReaction, s::SurfaceSupport)
    s.host === nothing && return (_mineral_name_and_mass(cs, rxn)..., s.area)
    sp = get(cs.dict_species, s.host, nothing)
    sp === nothing && throw(
        ArgumentError(
            "SurfaceSupport \"$(s.name)\" names host \"$(s.host)\", which is not a species " *
                "of this system. Known symbols include " *
                "$(join(sort(collect(keys(cs.dict_species)))[1:min(end, 6)], ", ")), …",
        )
    )
    return (_rate_lookup_key(cs, sp), _molar_mass_si(sp), s.area)
end

# Returns Vector of (name::String, ν::Float64, ΔG_fn) for all species in rxn
# that are present in cs and have a :ΔₐG⁰ property.
function _stoich_named(cs::ChemicalSystem, rxn::AbstractReaction)
    result = Tuple{String, Float64, Any}[]
    for (sp, ν) in rxn.reactants
        i = findfirst(s -> s == sp, cs.species)
        isnothing(i) && continue
        sp_cs = cs.species[i]
        haskey(properties(sp_cs), :ΔₐG⁰) || continue
        push!(result, (phreeqc(formula(sp_cs)), -Float64(ν), sp_cs[:ΔₐG⁰]))
    end
    for (sp, ν) in rxn.products
        i = findfirst(s -> s == sp, cs.species)
        isnothing(i) && continue
        sp_cs = cs.species[i]
        haskey(properties(sp_cs), :ΔₐG⁰) || continue
        push!(result, (phreeqc(formula(sp_cs)), Float64(ν), sp_cs[:ΔₐG⁰]))
    end
    return result
end

# ── molar_mass ────────────────────────────────────────────────────────────────

"""
    molar_mass(kr::KineticReaction) -> Float64

Return the molar mass of the mineral species [kg/mol], used by every
specific-area model to turn an amount into an area.

Searches `kr.reaction.reactants` for a species with an `:M` property, and raises
when none has one: a default here is worse than a refusal, because it scales the
reactive area and therefore the whole rate.
"""
function molar_mass(kr::KineticReaction)
    for (sp_obj, _) in kr.reaction.reactants
        haskey(properties(sp_obj), :M) && return ustrip(us"kg/mol", sp_obj[:M])
    end
    throw(
        ArgumentError(
            "no reactant of \"$(kr.reaction.symbol)\" carries a molar mass `:M`, " *
                "which the reactive area needs.",
        )
    )
end
