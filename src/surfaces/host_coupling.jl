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
    declared_site_moles(state::ChemicalState, family::SiteFamily) -> Real

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
    cap isa TotalSiteAmount && return float(site_moles(family, 0.0, 0.0, 0.0))

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
    # In the number type of the capacity and of the state: a site density or a
    # host amount being differentiated carries its derivative into the budget.
    n_host = ustrip(us"mol", state.n[i])
    M = _molar_mass_si(sp)
    return float(site_moles(family, n_host, n_host, M))
end

"""
    present_site_moles(state::ChemicalState, family::SiteFamily) -> Real

The moles of sites the state actually carries for `family`: `Σ dᵢ nᵢ` over its
members, with `dᵢ` the denticity read from each member's formula.

This is exactly the left-hand side of the site conservation row, so comparing it
with [`declared_site_moles`](@ref) compares the declaration against what the
solver will conserve.
"""
function present_site_moles(state::ChemicalState, family::SiteFamily)
    cs = state.system
    total = zero(_realtype(eltype(state.n)))
    for sp in site_members(family)
        i = findfirst(s -> symbol(s) == symbol(sp), cs.species)
        i === nothing && continue
        total += denticity(family, sp) * ustrip(us"mol", state.n[i])
    end
    return total
end

"""
    site_budget_residual(state::ChemicalState) -> OrderedDict{String, <:Real}

Per declared family, `present − declared` in mol.

Zero means the state's site amounts and the family's capacity say the same
thing. Anything else means the calculation will conserve the amounts and ignore
the declaration, which is what it did unconditionally before this existed.

Returns an empty dictionary for a system with no surface, so it is safe to call
on anything.

See also: [`host_consistent_state`](@ref), [`InconsistentSiteBudget`](@ref).
"""
function site_budget_residual(state::ChemicalState)
    fams = state.system.site_families
    fams === nothing && return OrderedDict{String, Float64}()
    r = [name(f) => present_site_moles(state, f) - declared_site_moles(state, f) for f in fams]
    return OrderedDict{String, mapreduce(typeof ∘ last, promote_type, r; init = Float64)}(r)
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
            throw(InconsistentSiteBudget(name(f), _plain(r), _plain(declared), _plain(present)))
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
        occupied = zero(_realtype(eltype(n)))
        for sp in f.complexes
            i = findfirst(s -> symbol(s) == symbol(sp), cs.species)
            i === nothing && continue
            occupied += denticity(f, sp) * ustrip(us"mol", n[i])
        end
        free = declared - occupied
        free >= 0 || throw(
            ArgumentError(
                "SiteFamily \"$(name(f))\" declares $(_plain(declared)) mol of sites, but its " *
                    "complexes already occupy $(_plain(occupied)) mol. There is no free-site " *
                    "amount that makes the two agree: raise the capacity or lower the " *
                    "occupied amounts. Scaling the complexes down to fit would invent " *
                    "a composition, and move the sorbates' elements with it.",
            )
        )
        j = findfirst(s -> symbol(s) == symbol(reference_member(f)), cs.species)
        if j !== nothing
            # A capacity being differentiated makes the free site, and so the
            # state, dual: the amounts are promoted to hold it.
            n = collect(promote_type(eltype(n), typeof(free * u"mol")), n)
            n[j] = free * u"mol"
        end
    end
    return ChemicalState(cs, n; T = state.T[1], P = state.P[1])
end

# ── The site-density scale an adsorption constant refers to ──────────────────

"""
    REFERENCE_SITE_DENSITY_NM2

The conventional reference site density `Γ° = $(literature_value("Kulik2002", "reference_site_density_nm2")) nm⁻²` of
[Kulik2002](@cite), used to define the standard state of a monodentate surface
species, read from `data/literature/Kulik2002.json`.

# What it is for

An intrinsic adsorption constant is not a property of a surface alone: it is
fitted at some **total site density** `Γ_C`, and the value depends on that
choice. [DzombakMorel1990](@citet) fitted their hydrous-ferric-oxide constants at two
site densities, one for the weak sites and one for the strong, a factor of forty
apart, so two of their own constants are not directly comparable with each
other, let alone with a constant from another compilation.

Fixing one conventional `Γ°` for every sorbent and every surface makes them
comparable, and [`convert_logk_site_density`](@ref) is the conversion. This is
a **choice of concentration scale**, not a measurement: it is one number, it is
stated here rather than assumed, and it can be changed.
"""
const REFERENCE_SITE_DENSITY_NM2 = literature_value("Kulik2002", "reference_site_density_nm2")

"""
    REFERENCE_SITE_DENSITY

[`REFERENCE_SITE_DENSITY_NM2`](@ref) in `mol/m²`, **derived** through
[`AVOGADRO`](@ref) rather than written again — `2.0009e-5 mol/m²`.

That it comes out so close to a round `2 × 10⁻⁵ mol/m²` is not a coincidence:
that is the number the convention was chosen to be, and the value in `nm⁻²` is
how it reads in the units the surface literature uses.
"""
const REFERENCE_SITE_DENSITY = REFERENCE_SITE_DENSITY_NM2 * 1.0e18 / AVOGADRO

# The spectral gap that decides whether a coupled matrix determines its site
# potential, used by `_refuse_unidentifiable_site`. It is `identifiable_rank`'s
# own default, kept rather than tuned, and the five systems it was checked on
# sit well clear of it on both sides: the tightest refusal has a ratio of 13
# across the gap (`1.39 / 0.105`, hydrous ferric oxide over a neutral component)
# and the tightest acceptance a ratio of 2.2 (`1.42 / 0.656`, the same oxide in
# a mixed-valence iron chloride solution). A threshold read off one machine's
# solve would not have that margin; this one is read off the matrix.
const _SITE_RANK_GAP = 5.0

"""
    convert_logk_site_density(logK, Γ_C; Γ0 = REFERENCE_SITE_DENSITY_NM2,
                              free_site_side = :reactant) -> Float64

An intrinsic adsorption constant moved from the total site density it was
fitted at to a reference one — equation (21) of [Kulik2002](@citet):

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

The two densities of [DzombakMorel1990](@citet), the case [Kulik2002](@citet)
uses to make the point that their weak and strong constants are not on one
scale:

```jldoctest
julia> using ChemistryLab

julia> weak = literature_value("Kulik2002", "dzombak_morel_weak_site_density_nm2");

julia> strong = literature_value("Kulik2002", "dzombak_morel_strong_site_density_nm2");

julia> round(convert_logk_site_density(0.0, weak), digits = 2)
-0.73

julia> round(convert_logk_site_density(0.0, strong), digits = 2)
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
    sites_per_host(family::SiteFamily, M_host) -> Real

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
    # In the number type of the capacity: a site density being differentiated
    # carries its derivative into the conservation matrix.
    ν = float(site_moles(cap, support, 1.0, 1.0, M_host))

    ν > 0 || throw(
        ArgumentError(
            "SiteFamily \"$(name(family))\" evaluates to $(_plain(ν)) mol of sites per mole of " *
                "host. A coupled family with no sites is a declaration to correct: " *
                "check the capacity, and the host's molar mass if the capacity is " *
                "measured per unit mass or area.",
        )
    )

    for λ in (0.37, 2.9)
        got = _plain(site_moles(cap, support, λ, 1.0, M_host))
        isapprox(got, λ * _plain(ν); rtol = 1.0e-10) || throw(
            ArgumentError(
                "SiteFamily \"$(name(family))\" is declared SITES_FOLLOW_HOST, but its " *
                    "$(nameof(typeof(cap))) is not proportional to the host's amount: " *
                    "scaling that amount by $λ changes the site budget by " *
                    "$(round(got / _plain(ν); sigdigits = 4)) instead. Measured, not assumed.\n" *
                    "A budget that follows the host has to be a coefficient times its " *
                    "amount, or the site row stops being linear and the problem stops " *
                    "being a polyhedron. Use a capacity measured per unit mass or per " *
                    "unit specific area, or leave the support at SITES_FIXED.",
            )
        )
    end
    return ν
end

# ── The coupling: the host's formula includes its sites ──────────────────────

"""
    site_coupling_rows(cs::ChemicalSystem) -> (Matrix, Vector{String})

One row per family whose support is `SITES_FOLLOW_HOST`, stating that the sites
in use equal `ν` times the host's amount:

```math
\\sum_k d_k\\, n_k \\;-\\; \\nu\\, n_{\\text{host}} \\;=\\; 0
```

with `d_k` the denticity of each member. Returns the rows and the family names
that label them; both are empty when nothing is coupled, so a system without a
coupled surface gets back exactly what it had.

These rows state the constraint so that it can be checked, and
[`conservation_matrix`](@ref) is how it is imposed: on the site row the two
coincide. Appending them to `SM.A` instead would not work, since `SM.A` already
carries a site row and two equations on one sum forbid the host to move.
"""
function site_coupling_rows(cs::ChemicalSystem)
    empty_rows = Matrix{Float64}(undef, 0, length(cs.species))
    fams = cs.site_families
    fams === nothing && return empty_rows, String[]
    coupled = [f for f in fams if surface_support(f).coupling === SITES_FOLLOW_HOST]
    isempty(coupled) && return empty_rows, String[]

    hosts = map(coupled) do f
        host = surface_support(f).host
        j = findfirst(s -> symbol(s) == host, cs.species)
        j === nothing && throw(
            ArgumentError(
                "SiteFamily \"$(name(f))\" follows host \"$host\", which is not a " *
                    "species of this system. A coupled family needs the solid whose " *
                    "amount its sites track.",
            )
        )
        (j, sites_per_host(f, _molar_mass_si(cs.species[j])))
    end
    # In the number type of the capacities.
    out = zeros(mapreduce(typeof ∘ last, promote_type, hosts; init = Float64), length(coupled), length(cs.species))
    labels = String[]
    for (r, f) in enumerate(coupled)
        j, ν = hosts[r]
        for sp in site_members(f)
            i = findfirst(s -> symbol(s) == symbol(sp), cs.species)
            i === nothing && continue
            out[r, i] += denticity(f, sp)
        end
        out[r, j] -= ν
        push!(labels, name(f))
    end
    return out, labels
end

"""
    _is_bare_site(sp, site::Symbol) -> Bool

Whether `sp` carries the site pseudo-element `site` **and nothing else** — no
real atom, possibly a charge.

Such a species is a *component* rather than a substance: the site with no
matter attached. A coupled family may take one as the primary of its site row,
or take its free site; [`conservation_matrix`](@ref) holds in either basis. The
test matters for `_refuse_unidentifiable_site`, which concerns the bare basis
only.
"""
function _is_bare_site(sp::AbstractSpecies, site::Symbol)
    a = atoms(sp)
    get(a, site, 0) == 1 || return false
    # A CHARGE IS ALLOWED, and usually required. The bare component carries no
    # matter, but it does carry whatever charge the free site carries with its
    # site symbol: `XsOH` is `Xs⁺ + OH⁻`, so the component is `Xs+`; an
    # exchanger `NaXc` is `Xc⁻ + Na⁺`, so its component is negative. Insisting
    # on neutrality here is what makes the basis degenerate — see
    # `_refuse_parasitic_charge`.
    return all(k === site || k === :Zz || iszero(v) for (k, v) in a)
end

"""
    _refuse_unidentifiable_site(cs, family, A, r)

Refuse a coupled family whose site potential the basis cannot determine.

# What goes wrong, measured

Every species carrying a site symbol carries it together with a fixed amount of
charge: `XsOH` decomposes as `H₂O − H⁺ + Xs + Zz` and `XsOCa⁺` as
`Ca²⁺ + H₂O − 2H⁺ + Xs + Zz`. Where that is the *only* way `Zz` enters, the site
row and the charge row of the coupled matrix are the same row up to the host
entry, so only `y_Xs + y_Zz` is identifiable.

The solver finds that out the hard way. Measured on portlandite carrying sites,
with a neutral bare component: the two multipliers ran to `+234 650` and
`−234 710` — four decades past any chemical potential — while their sum stayed
at `−60`; the dual Newton stalled at a KKT error of `2.3e-9`; and the host came
out at `0.0999 mol` where the uncoupled answer is `0.0883`.

Giving the component the charge the free site carries with its site symbol
removes the `Zz` primary altogether. The same run then converges at `4.8e-11`,
`max|y|` is `227`, the host lands on `0.08826` against `0.08827` uncoupled, and
the coupling holds to `1.8e-7`.

# Why this is measured and not inferred from `Zz` being a primary

Refusing on the *presence* of a charge component is what this did first, and it
refuses systems that are perfectly well posed. `Zz` survives as a primary
whenever charge is genuinely independent of the element rows — which is the
ordinary situation in a redox system. Measured, on hydrous ferric oxide over
`Fe(OH)₃(am)` with chloride and a mixed-valence iron speciation, singular
values of the coupled matrix:

| system | component | `σ` (smallest three) | identifiable rank |
|:--|:--|--:|:--|
| portlandite, no `Zz` row possible | `Xs+` | `1.19, 0.872` | 4 of 4 |
| portlandite | `Xs` | `1.21, 0.957, 1.7e-5` | **4 of 5** |
| HFO, Fe(III) only | `Xs` | `2.16, 1.39, 0.105` | **4 of 5** |
| HFO, Fe(II) and Fe(III), chloride | `Xs+` | `1.90, 1.35, 0.916` | 6 of 6 |
| HFO, Fe(II) and Fe(III), chloride | `Xs` | `2.14, 1.42, 0.656` | 6 of 6 |

The last two carry a `Zz` primary and are not degenerate at all — the charge row
is nonzero on an iron chloride complex, so it is not the site row. The old test
refused both of them, and the charge it then suggested was itself refused on the
next call: the two suggestions pointed at each other and the user went in a
circle.

[`identifiable_rank`](@ref) on the singular values separates the two groups with
the spectral gap it is built for — a factor of 13 at the tightest refusal
against 2.2 at the tightest acceptance. The charge suggested on a refusal is the
coefficient of `Zz` in the free site's own decomposition, which is correct
exactly where the refusal now fires.
"""
function _refuse_unidentifiable_site(
        cs::ChemicalSystem, family::SiteFamily, A::AbstractMatrix, r::Int,
    )
    zz = findfirst(p -> symbol(p) == "Zz", cs.SM.primaries)
    zz === nothing && return nothing
    # The decision is the identifiable rank of the matrix the solve will be
    # constrained with, NOT the presence of a charge component. See the
    # docstring for the five systems this was calibrated on.
    # On its values: a verdict, whatever the matrix carries.
    Av = _plain.(A)
    σ = svdvals(Av)
    identifiable_rank(σ; gap = _SITE_RANK_GAP) == size(A, 1) && return nothing

    site_row, charge_row = Av[r, :], Av[zz, :]
    nn = norm(site_row) * norm(charge_row)
    cosine = iszero(nn) ? 1.0 : abs(dot(site_row, charge_row)) / nn
    jf = findfirst(s -> symbol(s) == symbol(reference_member(family)), cs.species)
    z = jf === nothing ? 1 : Int(round(Float64(cs.SM.A[zz, jf])))
    sug = z == 0 ? "" : (z > 0 ? "+"^z : "-"^(-z))
    throw(
        ArgumentError(
            "SiteFamily \"$(name(family))\" follows its host over the bare component " *
                "\"$(symbol(cs.SM.primaries[r]))\", and the coupled matrix cannot " *
                "determine its site potential: the singular values end " *
                "$(join(string.(round.(σ[max(1, end - 2):end]; sigdigits = 3)), ", ")) " *
                "and the site row sits at |cos| = $(round(cosine; digits = 5)) from the " *
                "charge row, so only `y_$(family.site) + y_Zz` is identifiable.\n" *
                "Measured, the cost is not subtle: the two multipliers ran to " *
                "±2.3e5 while their sum stayed at −60, the solve stalled at a KKT " *
                "error of 2.3e-9, and the host came out thirteen percent off.\n" *
                "Declare the component with the charge the free site carries with its " *
                "site symbol — here `Species(\"$(family.site)$sug\")`.",
        )
    )
end


# Above this many log units on the host's saturation index, a coupled family is
# refused rather than solved. It is a statement about how much of an answer the
# free site's reference is allowed to be: 0.05 log units is 12 % on a
# solubility, below the spread between two databases for the same phase.
const _MAX_COUPLING_BIAS = 0.05

"""
    host_coupling_bias(cs::ChemicalSystem) -> OrderedDict{String, Float64}

Per coupled family, the shift in `log SI` that the reference energy of the free
site imposes on the host, in the units the answer is read in:

```math
\\text{bias} = \\frac{\\nu\\,\\left|\\Delta_a G^0_{\\text{free}}\\right|}{RT \\ln 10} .
```

# Why the free site of a coupled family has zero energy

With a fixed budget the free site's `ΔₐG⁰` cancels out of every surface
reaction, both sides carrying one site, and it is a gauge. With a budget that
follows its host it no longer cancels. [`conservation_matrix`](@ref) counts the
free sites as part of the host, so that an intact sorbent, every site free, has
the composition of the host's database formula; it has the database energy only
when the free site's energy is zero. That is the convention of
[Kulik2002](@cite) for surface groups that belong to their sorbent, and any
other value moves the host's solubility by the amount above. A family whose bias
exceeds `$(_MAX_COUPLING_BIAS)` log units is refused, and
[`site_family`](@ref) builds families that satisfy it, its free site at zero and
every complex written relative to it.

The energies are read at 298.15 K, where the surface constants are given.

See also: [`sites_per_host`](@ref), [`conservation_matrix`](@ref).
"""
function host_coupling_bias(cs::ChemicalSystem)
    out = OrderedDict{String, Real}()
    fams = cs.site_families
    fams === nothing && return out
    RT = R_GAS * 298.15
    for f in fams
        surface_support(f).coupling === SITES_FOLLOW_HOST || continue
        jf = findfirst(s -> symbol(s) == symbol(reference_member(f)), cs.species)
        jh = findfirst(s -> symbol(s) == surface_support(f).host, cs.species)
        (jf === nothing || jh === nothing) && continue
        g = _standard_gibbs(cs.species[jf])
        g === nothing && continue          # no standard energy: nothing to measure
        ν = sites_per_host(f, _molar_mass_si(cs.species[jh]))
        out[name(f)] = ν * abs(g) / (RT * log(10))
    end
    return OrderedDict{String, mapreduce(typeof, promote_type, values(out); init = Float64)}(out)
end

"""
    _standard_gibbs(sp) -> Union{Float64, Nothing}

`ΔₐG⁰` at 298.15 K and 1 bar, in J/mol, or `nothing` when the species carries
none. A missing property reads back as an `Int` zero rather than as an error
(`species.jl`), so the test is on the type and not on the value: a species whose
energy is genuinely zero must not be confused with one that has none.
"""
function _standard_gibbs(sp::AbstractSpecies)
    g = sp[:ΔₐG⁰]
    g isa Number && return nothing
    return ustrip(us"J/mol", g(T = 298.15u"K", P = 1.0e5u"Pa"; unit = true))
end

"""
    _refuse_biased_coupling(cs, family)

Refuse a coupled family whose free site is not at zero energy; see
[`host_coupling_bias`](@ref) for why a coupled family needs it and a fixed one
does not.
"""
function _refuse_biased_coupling(cs::ChemicalSystem, family::SiteFamily)
    bias = get(host_coupling_bias(cs), name(family), 0.0)
    bias ≤ _MAX_COUPLING_BIAS && return nothing
    jf = findfirst(s -> symbol(s) == symbol(reference_member(family)), cs.species)
    g = _standard_gibbs(cs.species[jf])
    throw(
        ArgumentError(
            "SiteFamily \"$(name(family))\" follows its host, and its free site " *
                "\"$(symbol(cs.species[jf]))\" has ΔₐG⁰ = $(round(g / 1000; digits = 1)) " *
                "kJ/mol, which moves the host's own saturation index by " *
                "$(round(_plain(bias); sigdigits = 3)) log units.\n" *
                "The coupling counts the free sites as part of the host, so an intact " *
                "sorbent has the database energy of the host only when the free site's " *
                "energy is zero. Set ΔₐG⁰ of the free site to 0 and write every complex " *
                "relative to it, as `site_family` does; see `host_coupling_bias`.",
        )
    )
end

"""
    _refuse_charged_free_site(family, free)

Refuse a coupled family whose free site carries a charge. The coupling moves `ν`
free sites with each mole of host, so a charged one would carry charge in and out
of the system as the host dissolves or grows.
"""
function _refuse_charged_free_site(family::SiteFamily, free::AbstractSpecies)
    z = charge(free)
    iszero(z) && return nothing
    throw(
        ArgumentError(
            "SiteFamily \"$(name(family))\" follows its host, but its free site " *
                "\"$(symbol(free))\" carries a charge of $z. The coupling counts " *
                "ν free sites as part of each mole of host, so a charged free site " *
                "would create or destroy charge as the host dissolves or grows.\n" *
                "Declare a neutral free site (for an exchanger, the site with its " *
                "compensating cation), or keep the support at SITES_FIXED.",
        )
    )
end

"""
    _refuse_host_short_of_site_matter(family, host, free, ν)

Refuse a coupled family whose host formula does not contain `ν` free sites' worth
of every element: the free sites are counted as part of the host, so their atoms
have to be among the host's own.
"""
function _refuse_host_short_of_site_matter(
        family::SiteFamily, host::AbstractSpecies, free::AbstractSpecies, ν::Real,
    )
    ah, af = atoms(host), atoms(free)
    for (e, k) in af
        (e === family.site || e === :Zz) && continue
        left = Float64(get(ah, e, 0)) - _plain(ν) * Float64(k)
        left ≥ -1.0e-12 && continue
        throw(
            ArgumentError(
                "SiteFamily \"$(name(family))\" follows host \"$(symbol(host))\" with " *
                    "ν = $(round(_plain(ν); sigdigits = 4)) sites per mole, but its free site " *
                    "\"$(symbol(free))\" carries $k $e per site and the host's formula " *
                    "holds only $(get(ah, e, 0)). The coupling counts the free sites " *
                    "as part of the host, so their atoms have to be among the host's " *
                    "own. Check the capacity, or the choice of host.",
            )
        )
    end
    return nothing
end

"""
    conservation_matrix(cs::ChemicalSystem) -> Matrix

The matrix the equilibrium is constrained with: `SM.A` when no family follows
its host, and otherwise `SM.A` with, for each coupled family, `ν` times the
column of its free site subtracted from the column of its host,

```math
A'_{:,\\,\\text{host}} = A_{:,\\,\\text{host}} - \\nu\\, A_{:,\\,\\text{free}} .
```

The host's database formula already contains the surface groups its sites are
made of, since a hydroxide carries the hydroxyls its surface exposes, so the free
sites are counted as part of the host rather than on top of it. On the site row
the subtraction reads `Σ dₖ nₖ − ν n_host = 0`, the statement of
[`site_coupling_rows`](@ref); on the element rows it removes from the host the
atoms the free sites carry, so that every atom is counted once; and, the free
site being neutral, it leaves the charge conserved while the host dissolves or
grows. This is the bookkeeping of [Kulik2002](@citet) for a sorbent whose surface
groups belong to it, and it holds in any basis of primaries: the free site may be
the primary of the site row, or a bare site component (`Species("Xs+")`) may be
declared in its place.

Refused, each by name: a free site that carries a charge, a host whose formula
does not contain `ν` free sites' worth of every element, a capacity that is not
proportional to the host's amount ([`sites_per_host`](@ref)), and a free site
whose energy is away from zero ([`host_coupling_bias`](@ref)).
"""
function conservation_matrix(cs::ChemicalSystem)
    A = Float64.(cs.SM.A)
    fams = cs.site_families
    fams === nothing && return A
    # In the number type of the capacities: a site density being differentiated
    # makes the host column dual. Under a solve on values (`_STRIP_TAGS`), on its
    # values.
    R = mapreduce(promote_type, fams; init = Float64) do f
        surface_support(f).coupling === SITES_FOLLOW_HOST || return Float64
        j = findfirst(s -> symbol(s) == surface_support(f).host, cs.species)
        j === nothing ? Float64 : typeof(_unscoped(sites_per_host(f, _molar_mass_si(cs.species[j]))))
    end
    A = Matrix{R}(A)
    A0 = copy(A)
    for f in fams
        surface_support(f).coupling === SITES_FOLLOW_HOST || continue
        r = findfirst(p -> get(atoms(p), f.site, 0) > 0, cs.SM.primaries)
        r === nothing && throw(
            ArgumentError(
                "SiteFamily \"$(name(f))\" follows its host, but no primary of this " *
                    "system carries :$(f.site), so there is no site row to couple.",
            )
        )
        host = surface_support(f).host
        j = findfirst(s -> symbol(s) == host, cs.species)
        j === nothing && throw(
            ArgumentError(
                "SiteFamily \"$(name(f))\" follows host \"$host\", which is not a " *
                    "species of this system.",
            )
        )
        free = reference_member(f)
        jf = findfirst(s -> symbol(s) == symbol(free), cs.species)
        jf === nothing && throw(
            ArgumentError(
                "SiteFamily \"$(name(f))\" follows its host, but its free site " *
                    "\"$(symbol(free))\" is not a species of this system.",
            )
        )
        _refuse_charged_free_site(f, cs.species[jf])
        ν = _unscoped(sites_per_host(f, _molar_mass_si(cs.species[j])))
        _refuse_host_short_of_site_matter(f, cs.species[j], cs.species[jf], ν)
        @views A[:, j] .-= ν .* A0[:, jf]
        # The bare basis is the one whose site row can coincide with the charge
        # row; the check reads the matrix the solve will use.
        _is_bare_site(cs.SM.primaries[r], f.site) && _refuse_unidentifiable_site(cs, f, A, r)
        _refuse_biased_coupling(cs, f)
    end
    return A
end

"""
    _constraint_matrix(cs::ChemicalSystem) -> AbstractMatrix

The matrix every solver path constrains the equilibrium with: `cs.SM.A` itself
when no family follows its host, and [`conservation_matrix`](@ref) otherwise.
Returning `SM.A` untouched keeps a system without a coupled family bit-identical
to what it was, and routing every path through here keeps a coupled one from
being solved on the uncoupled matrix by one of them.
"""
_constraint_matrix(cs::ChemicalSystem) =
    _has_coupled_family(cs) ? conservation_matrix(cs) : cs.SM.A

_has_coupled_family(cs::ChemicalSystem) =
    cs.site_families !== nothing &&
    any(f -> surface_support(f).coupling === SITES_FOLLOW_HOST, cs.site_families)
