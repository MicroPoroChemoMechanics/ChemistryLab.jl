# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

# ── The site budget, and whether anything checks it ──────────────────────────
#
# A `SiteFamily` declares a capacity. Until this file, nothing evaluated it:
# `site_moles` had no caller anywhere in `src`, and the site conservation row
# took its right-hand side from `b = A n₀` like every other row — that is, from
# the amounts the caller happened to put in the state. Declaring `1e-6 mol` of
# sites and initializing `1e-3 mol` of free sites gave an equilibrium holding
# `1e-3`, with `optimality_certificate` reporting `optimal = true`, because the
# certificate checks the budget it was given and has no way to know that the
# declaration said something else.
#
# The two numbers are not redundant: one is a physical statement about the
# sorbent, the other is a starting composition. They simply have to agree, and
# nothing made them.

"""
    InconsistentSiteBudget <: Exception

A state whose site amounts contradict the capacity its family declares.

Carries the family, the residual in mol, and the scale it was judged against,
so the message says how far off it is rather than only that it is off.
"""
struct InconsistentSiteBudget <: Exception
    family::String
    residual::Float64
    declared::Float64
    present::Float64
end

function Base.showerror(io::IO, e::InconsistentSiteBudget)
    return print(
        io,
        "InconsistentSiteBudget: SiteFamily \"", e.family, "\" declares ",
        e.declared, " mol of sites, but the state carries ", e.present,
        " mol across its members (residual ", e.residual, " mol).\n",
        "A capacity is a statement about the sorbent and the amounts are a ",
        "starting composition; they have to agree, and nothing but this check ",
        "makes them. Use `host_consistent_state` to derive the free-site ",
        "amount from the declaration, or correct whichever of the two is wrong.",
    )
end

"""
    declared_site_moles(state::ChemicalState, family::SiteFamily) -> Float64

The moles of sites `family` declares, evaluated on `state` — that is,
[`site_moles`](@ref) supplied with the host amount and molar mass the state
carries.

A [`TotalSiteAmount`](@ref) ignores all three and returns its number. The other
two do not, so a family that measures its capacity per unit mass or per unit
area **needs a host**, and one declared without naming a species is refused here
rather than quietly evaluated at zero — a capacity of zero is a different model,
not a missing input.
"""
function declared_site_moles(state::ChemicalState, family::SiteFamily)
    cap = site_capacity(family)
    support = surface_support(family)
    # A capacity stores a bare `Real` in SI — `_area_si` strips the unit at
    # construction — so there is nothing to unwrap here.
    cap isa TotalSiteAmount && return Float64(site_moles(family, 0.0, 0.0, 0.0))

    host = support.host
    host === nothing && throw(
        ArgumentError(
            "SiteFamily \"$(name(family))\" measures its capacity with a " *
                "$(nameof(typeof(cap))), which needs the amount of the solid carrying " *
                "the sites, but its SurfaceSupport names no host. Give the support a " *
                "host species, or declare the capacity as a `TotalSiteAmount`.",
        )
    )
    cs = state.system
    sp = get(cs.dict_species, host, nothing)
    sp === nothing && throw(
        ArgumentError(
            "SiteFamily \"$(name(family))\" names host \"$host\", which is not a " *
                "species of this system.",
        )
    )
    i = findfirst(s -> symbol(s) == host, cs.species)
    n_host = Float64(ustrip(us"mol", state.n[i]))
    M = _molar_mass_si(sp)
    return Float64(site_moles(family, n_host, n_host, M))
end

"""
    present_site_moles(state::ChemicalState, family::SiteFamily) -> Float64

The moles of sites the state actually carries for `family`: `Σ dᵢ nᵢ` over its
members, with `dᵢ` the denticity read from each member's formula.

This is exactly the left-hand side of the site conservation row, so comparing it
with [`declared_site_moles`](@ref) compares the declaration against what the
solver will conserve.
"""
function present_site_moles(state::ChemicalState, family::SiteFamily)
    cs = state.system
    total = 0.0
    for sp in site_members(family)
        i = findfirst(s -> symbol(s) == symbol(sp), cs.species)
        i === nothing && continue
        total += denticity(family, sp) * Float64(ustrip(us"mol", state.n[i]))
    end
    return total
end

"""
    site_budget_residual(state::ChemicalState) -> OrderedDict{String, Float64}

Per declared family, `present − declared` in mol.

Zero means the state's site amounts and the family's capacity say the same
thing. Anything else means the calculation will conserve the amounts and ignore
the declaration, which is what it did unconditionally before this existed.

Returns an empty dictionary for a system with no surface, so it is safe to call
on anything.

See also: [`host_consistent_state`](@ref), [`InconsistentSiteBudget`](@ref).
"""
function site_budget_residual(state::ChemicalState)
    out = OrderedDict{String, Float64}()
    fams = state.system.site_families
    fams === nothing && return out
    for f in fams
        out[name(f)] = present_site_moles(state, f) - declared_site_moles(state, f)
    end
    return out
end

"""
    check_site_budget(state::ChemicalState; rtol = 1e-6)

Throw [`InconsistentSiteBudget`](@ref) when a family's state amounts and its
declared capacity disagree by more than `rtol` of the declaration.

The tolerance is **relative to the declaration**, not absolute: a site budget is
`1e-6 mol` on an oxide and `1e-1 mol` on a clay, so an absolute threshold would
be meaningless on one of them. A declaration of exactly zero is compared
absolutely against the same `rtol`, since nothing else is available.

# Where `1e-6` comes from

From the floor, not from a round number. A state is normally initialized with
the occupied sites at the solver's `1e-12` floor rather than at zero, so a
family with two complexes and a `1e-3 mol` capacity is *correctly* initialized
and still `2e-9` off in relative terms — measured. A default of `1e-9` rejects
that, which is a check that fires on the ordinary case.

This exists to catch a **declaration** error, and those are orders of magnitude:
the case it was written for declared `1e-6 mol` against a state carrying
`1e-3`, a factor of a thousand. `1e-6` catches that with nine decades to spare
and leaves the floor alone.
"""
function check_site_budget(state::ChemicalState; rtol::Real = 1.0e-6)
    fams = state.system.site_families
    fams === nothing && return nothing
    for f in fams
        declared = declared_site_moles(state, f)
        present = present_site_moles(state, f)
        scale = max(abs(declared), abs(present))
        iszero(scale) && continue
        r = present - declared
        abs(r) <= rtol * scale ||
            throw(InconsistentSiteBudget(name(f), r, declared, present))
    end
    return nothing
end

"""
    host_consistent_state(state::ChemicalState) -> ChemicalState

`state` with each family's **free-site** amount derived from its declared
capacity instead of taken from the caller:

```math
n_{\\text{free}} = N_t - \\sum_{\\text{complexes}} d_i\\, n_i
```

so `site_budget_residual` comes back zero and the conservation row the solver
forms says what the family declared.

The free site is the one adjusted, and deliberately: it is the family's
reference member, the species whose amount carries no chemistry of its own. An
occupied site holds a sorbate whose element budget the rest of the system
accounts for, so moving it would silently move that element too.

Refused when the complexes already hold more sites than the capacity allows —
that is not an initialization to repair but a declaration to correct, and
scaling the complexes down to fit would invent a composition nobody asked for.

Everything outside the site families is copied unchanged.
"""
function host_consistent_state(state::ChemicalState)
    fams = state.system.site_families
    fams === nothing && return state
    cs = state.system
    n = copy(state.n)
    for f in fams
        declared = declared_site_moles(state, f)
        occupied = 0.0
        for sp in f.complexes
            i = findfirst(s -> symbol(s) == symbol(sp), cs.species)
            i === nothing && continue
            occupied += denticity(f, sp) * Float64(ustrip(us"mol", n[i]))
        end
        free = declared - occupied
        free >= 0 || throw(
            ArgumentError(
                "SiteFamily \"$(name(f))\" declares $declared mol of sites, but its " *
                    "complexes already occupy $occupied mol. There is no free-site " *
                    "amount that makes the two agree: raise the capacity or lower the " *
                    "occupied amounts. Scaling the complexes down to fit would invent " *
                    "a composition, and move the sorbates' elements with it.",
            )
        )
        j = findfirst(s -> symbol(s) == symbol(reference_member(f)), cs.species)
        j === nothing || (n[j] = free * u"mol")
    end
    return ChemicalState(cs, n; T = state.T[1], P = state.P[1])
end

# ── The site-density scale an adsorption constant refers to ──────────────────

"""
    REFERENCE_SITE_DENSITY_NM2

The conventional reference site density `Γ° = 12.05 nm⁻²` of
[Kulik2002](@cite), used to define the standard state of a monodentate surface
species.

# What it is for

An intrinsic adsorption constant is not a property of a surface alone: it is
fitted at some **total site density** `Γ_C`, and the value depends on that
choice. Dzombak and Morel fitted their hydrous-ferric-oxide constants at
`2.254 nm⁻²` for the weak sites and `0.056 nm⁻²` for the strong ones, a factor
of forty apart, so two of their own constants are not directly comparable with
each other, let alone with a constant from another compilation.

Fixing one conventional `Γ°` for every sorbent and every surface makes them
comparable, and [`convert_logk_site_density`](@ref) is the conversion. This is
a **choice of concentration scale**, not a measurement: it is one number, it is
stated here rather than assumed, and it can be changed.
"""
const REFERENCE_SITE_DENSITY_NM2 = 12.05

"""
    REFERENCE_SITE_DENSITY

[`REFERENCE_SITE_DENSITY_NM2`](@ref) in `mol/m²`, **derived** through
[`AVOGADRO`](@ref) rather than written again — `2.0009e-5 mol/m²`.

That it comes out so close to a round `2 × 10⁻⁵ mol/m²` is not a coincidence:
that is the number the convention was chosen to be, and `12.05 nm⁻²` is how it
reads in the units the surface literature uses.
"""
const REFERENCE_SITE_DENSITY = REFERENCE_SITE_DENSITY_NM2 * 1.0e18 / AVOGADRO

"""
    convert_logk_site_density(logK, Γ_C; Γ0 = REFERENCE_SITE_DENSITY_NM2,
                              free_site_side = :reactant) -> Float64

An intrinsic adsorption constant moved from the total site density it was
fitted at to a reference one — equation (21) of [Kulik2002](@cite):

```math
\\log K^\\circ = \\log K^{C} + \\log_{10}\\!\\frac{\\Gamma_C}{\\Gamma^\\circ}
```

with the sign of the last term **inverted** when the neutral functional group
sits on the product side of the reaction, which is what `free_site_side`
selects. Written `≡OH + H⁺ = ≡OH₂⁺` the free site is a reactant, the default.

Only the **ratio** of the two densities enters, so they may be given in any
unit as long as it is the same one. The default `Γ0` is in `nm⁻²`, which is how
the surface literature quotes a site density.

# Example

Dzombak and Morel's own two densities, which is the case Kulik uses to make
the point that their weak and strong constants are not on one scale:

```jldoctest
julia> using ChemistryLab

julia> round(convert_logk_site_density(0.0, 2.254), digits = 2)   # weak sites
-0.73

julia> round(convert_logk_site_density(0.0, 0.056), digits = 2)   # strong sites
-2.33
```

A constant fitted at a density **below** the reference moves down, and by
different amounts for the two site types — 1.6 log units apart here. Correlating
one against the other without this conversion compares two different scales.
"""
function convert_logk_site_density(
        logK::Real, Γ_C::Real;
        Γ0::Real = REFERENCE_SITE_DENSITY_NM2, free_site_side::Symbol = :reactant,
    )
    Γ_C > 0 || throw(ArgumentError("a site density must be positive; got $Γ_C."))
    Γ0 > 0 || throw(ArgumentError("the reference site density must be positive; got $Γ0."))
    free_site_side in (:reactant, :product) || throw(
        ArgumentError(
            "free_site_side is :reactant or :product; got :$free_site_side. It says " *
                "which side of the reaction the neutral functional group is written " *
                "on, and it flips the sign of the correction.",
        )
    )
    shift = log10(Γ_C / Γ0)
    return free_site_side === :reactant ? logK + shift : logK - shift
end

# ── ν, the moles of sites one mole of host carries ───────────────────────────

"""
    sites_per_host(family::SiteFamily, M_host) -> Float64

`ν`, the moles of sites one mole of the host carries, for a family whose
support is `SITES_FOLLOW_HOST`.

This is Kulik's equation (30) read as a coefficient: the moles of sites of a
surface type are `ψ · A_s · M · Γ` times the moles of sorbent, so the whole
dependence on the host is one number multiplying its amount.

# The capacity has to be homogeneous of degree one, and that is measured

`ν` only exists if the capacity really is proportional to the host's amount, and
which capacities are is not a list to maintain but a property to check. This
evaluates `site_moles` at two scaled amounts and requires the answer to scale
with them:

```math
\\text{site\\_moles}(\\lambda n, n_0, M) = \\lambda \\, \\text{site\\_moles}(n, n_0, M)
```

with `n₀` held **fixed**. Scaling `n₀` along with `n` would make a
[`ShrinkingCoreArea`](@ref)'s ratio `n/n₀` equal one at every `λ` and hide
exactly the nonlinearity this exists to find.

What the probe accepts and refuses follows without being enumerated:
[`MassSiteDensity`](@ref) and [`AreaSiteDensity`](@ref) over a specific area
pass, being `q M n` and `Γ a M n`; a [`TotalSiteAmount`](@ref) and an
`AreaSiteDensity` over a [`FixedSurfaceArea`](@ref) are refused, being constants
that do not follow anything; and a `ShrinkingCoreArea` passes only at `p = 1`,
where it genuinely is linear. A capacity type written later is judged on the
same evidence rather than on whether somebody remembered to add it here.

A capacity of zero is refused too: it passes the proportionality test for the
empty reason that `0 = λ·0`, and it is a family with no sites, which is a
declaration to correct.
"""
function sites_per_host(family::SiteFamily, M_host::Real)
    cap = site_capacity(family)
    support = surface_support(family)
    ν = Float64(site_moles(cap, support, 1.0, 1.0, M_host))

    ν > 0 || throw(
        ArgumentError(
            "SiteFamily \"$(name(family))\" evaluates to $ν mol of sites per mole of " *
                "host. A coupled family with no sites is a declaration to correct: " *
                "check the capacity, and the host's molar mass if the capacity is " *
                "measured per unit mass or area.",
        )
    )

    for λ in (0.37, 2.9)
        got = Float64(site_moles(cap, support, λ, 1.0, M_host))
        isapprox(got, λ * ν; rtol = 1.0e-10) || throw(
            ArgumentError(
                "SiteFamily \"$(name(family))\" is declared SITES_FOLLOW_HOST, but its " *
                    "$(nameof(typeof(cap))) is not proportional to the host's amount: " *
                    "scaling that amount by $λ changes the site budget by " *
                    "$(round(got / ν; sigdigits = 4)) instead. Measured, not assumed.\n" *
                    "A budget that follows the host has to be a coefficient times its " *
                    "amount, or the site row stops being linear and the problem stops " *
                    "being a polyhedron. Use a capacity measured per unit mass or per " *
                    "unit specific area, or leave the support at SITES_FIXED.",
            )
        )
    end
    return ν
end
