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

# ── The coupling as an ADDED row, never as a changed composition ─────────────

"""
    site_coupling_rows(cs::ChemicalSystem) -> (Matrix{Float64}, Vector{String})

One row per family whose support is `SITES_FOLLOW_HOST`, stating that the sites
in use equal `ν` times the host's amount:

```math
\\sum_k d_k\\, n_k \\;-\\; \\nu\\, n_{\\text{host}} \\;=\\; 0
```

with `d_k` the denticity of each member. Returns the rows and the family names
that label them; both are empty when nothing is coupled, so a system without a
coupled surface gets back exactly what it had.

# Why an added row, and not a coefficient in the composition matrix

Because the obvious alternative is **impossible**, and that is worth stating
once rather than rediscovering.

Writing the coupling as `A[site, host] -= ν` in the projected matrix changes
what the system conserves: the rows of `A` are indexed by primary species, so
subtracting there subtracts the free site's whole composition — and a free site
carries real atoms, an oxygen and a hydrogen for `XsOH`. The system then
conserves `M n − ν m_free n_host`, which creates `ν` moles of oxygen and `ν` of
hydrogen per mole of host. For Dzombak and Morel's weak sites, `ν = 0.2`: seven
percent of the oxygen of `Fe(OH)₃`, invented.

The repair would be to subtract `ν` times a preimage of the **pure** site
pseudo-element instead. No such preimage exists. Measured on two systems of
different structure — an amphoteric oxide over `[:H, :O, :Xs, :Zz]` and a
cation exchanger over `[:Na, :K, :H, :O, :Xc, :Zz]` — the least-squares residual
`‖M_indep v − Xs_unit‖` comes out `0.378` and `0.500`, not zero. The reason is
structural rather than a quirk of a basis: a site symbol never appears alone.
Every species carrying it carries it attached to matter, and only one site
species can be primary, so no combination of primaries yields a bare site with
every real atom and the charge at zero.

An added row has none of this to answer for. `SM.A` **is** the encoding of
element conservation, and leaving it alone leaves that conservation exact: the
surface species carry their own oxygen and hydrogen in their own formulas, the
host carries its own, and growing the site population draws them from the water
through the ordinary element rows. Automatically, with nothing to correct.

What it does cost is that `saturation_indices` reads `SM.A` and does not see
this row; the host's reported index has to be taught about it separately, or it
would disagree with the stationarity the solver actually reached.

# These rows state the constraint; they are not how it is imposed

[`conservation_matrix`](@ref) imposes it, as a single `−ν` in the site row's
host column. This returns the same statement in row form, which is what lets the
two be checked against each other rather than believed.

Appending these rows *instead* does not work, and the measurement is worth
keeping: `SM.A` already carries a site row, so a second one tying the same sum
to the host leaves two equations on one quantity, and together they say
`n_host = n_host,0` — the host may not dissolve at all. Run on portlandite, the
solve returned `MaxIters`, the host moved from 0.1 to 0.0883 mol regardless, and
the site total stayed at exactly `ν × 0.1`: the solver satisfied the old row and
violated the new one. The constraint has to **replace** the site row, not join
it.
"""
function site_coupling_rows(cs::ChemicalSystem)
    empty_rows = Matrix{Float64}(undef, 0, length(cs.species))
    fams = cs.site_families
    fams === nothing && return empty_rows, String[]
    coupled = [f for f in fams if surface_support(f).coupling === SITES_FOLLOW_HOST]
    isempty(coupled) && return empty_rows, String[]

    out = zeros(Float64, length(coupled), length(cs.species))
    labels = String[]
    for (r, f) in enumerate(coupled)
        host = surface_support(f).host
        j = findfirst(s -> symbol(s) == host, cs.species)
        j === nothing && throw(
            ArgumentError(
                "SiteFamily \"$(name(f))\" follows host \"$host\", which is not a " *
                    "species of this system. A coupled family needs the solid whose " *
                    "amount its sites track.",
            )
        )
        ν = sites_per_host(f, _molar_mass_si(cs.species[j]))
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
real atom, no charge.

Such a species is not a chemical species at all. It is a *component*: the pure
site, with no matter attached. That is exactly what a coupled family needs its
primary to be, and why one is allowed to sit among the primaries without being
among the species.
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
    _refuse_parasitic_charge(cs, family, r)

Refuse a coupled family whose bare component leaves `Zz` among the primaries.

# What goes wrong, measured

Every species carrying a site symbol carries it together with a fixed amount of
charge: `XsOH` decomposes as `H₂O − H⁺ + Xs + Zz` and `XsOCa⁺` as
`Ca²⁺ + H₂O − 2H⁺ + Xs + Zz`. The site component and the charge component
therefore appear in the **same ratio** everywhere, so only their sum is
identifiable and the basis is degenerate along `y_Xs + y_Zz`.

The solver finds that out the hard way. Measured on portlandite carrying sites,
with a neutral bare component: the two multipliers ran to `+234 650` and
`−234 710` — four decades past any chemical potential — while their sum stayed
at `−60`; the dual Newton stalled at a KKT error of `2.3e-9`; and the host came
out at `0.0999 mol` where the uncoupled answer is `0.0883`.

Giving the component the charge the free site carries with its site symbol
removes the `Zz` primary altogether. The same run then converges at `4.8e-11`,
`max|y|` is `227`, the host lands on `0.08826` against `0.08827` uncoupled, and
the coupling holds to `1.8e-7`.

The charge to use is not guessed: it is the coefficient of `Zz` in the free
site's own decomposition, which is what this reports.
"""
function _refuse_parasitic_charge(cs::ChemicalSystem, family::SiteFamily, r::Int)
    zz = findfirst(p -> symbol(p) == "Zz", cs.SM.primaries)
    zz === nothing && return nothing
    jf = findfirst(s -> symbol(s) == symbol(reference_member(family)), cs.species)
    z = jf === nothing ? 1 : Int(round(Float64(cs.SM.A[zz, jf])))
    sug = z == 0 ? "" : (z > 0 ? "+"^z : "-"^(-z))
    throw(
        ArgumentError(
            "SiteFamily \"$(name(family))\" follows its host over the bare component " *
                "\"$(symbol(cs.SM.primaries[r]))\", and that leaves :Zz among the " *
                "primaries. Every species carrying :$(family.site) carries it with a " *
                "fixed amount of charge, so the site and the charge components appear " *
                "in one ratio everywhere and only their sum is identifiable.\n" *
                "Measured, the cost is not subtle: the two multipliers ran to " *
                "±2.3e5 while their sum stayed at −60, the solve stalled at a KKT " *
                "error of 2.3e-9, and the host came out thirteen percent off.\n" *
                "Declare the component with the charge the free site carries with its " *
                "site symbol — here `Species(\"$(family.site)$sug\")`.",
        )
    )
end

"""
    conservation_matrix(cs::ChemicalSystem) -> Matrix{Float64}

The matrix the equilibrium is constrained with: `SM.A` as it stands when no
family follows its host, and `SM.A` with `ν` subtracted from each coupled
family's `(site row, host column)` entry when one does.

That single entry is the whole coupling. It states

```math
\\sum_k d_k n_k - \\nu n_{\\text{host}} = 0
```

and [`site_coupling_rows`](@ref) returns the same statement in row form, which
is how the two are checked against each other.

# Why one entry is enough, and why it was not before

Subtracting from a row of the projected matrix subtracts the **primary's whole
composition**, not the site symbol alone: `M = M_indep A`, so a correction `v`
in primary coordinates removes `M_indep v` of matter. Getting the pure site out
of it needs `M_indep v = Xs_unit`, a preimage of the bare pseudo-element.

With the free site for primary that preimage does **not exist**. Measured, as
the least-squares residual `‖M_indep v − Xs_unit‖`: `0.378` on an amphoteric
oxide over `[:H, :O, :Xs, :Zz]`, `0.500` on a cation exchanger over
`[:Na, :K, :H, :O, :Xc, :Zz]`. The reason is structural — a site symbol never
appears alone, every species carrying it carries it attached to matter, and only
one site species can be primary. Subtracting anyway invents `ν` moles of oxygen
and `ν` of hydrogen per mole of host: seven percent of the oxygen of `Fe(OH)₃`
at Dzombak and Morel's weak-site density.

Declaring the **bare** site as the primary removes the obstruction rather than
working around it. The residual is then `0.0` exactly and the preimage is the
unit vector, so subtracting `ν` from that one entry subtracts `ν` times a
component carrying no atom and no charge. Element conservation is exact by
construction, and nothing else in the matrix moves.

The component carries the **charge** the free site carries with its site symbol
— `XsOH` is `Xs⁺ + OH⁻`, so the component is `Xs+`. A neutral one leaves `Zz`
among the primaries and the basis is then degenerate; `_refuse_parasitic_charge`
says so, with what it costs.

# Measured

Uncoupled, the substitution costs nothing: the same system over a bare-site
component and over the free site returns the same host amount to eight digits,
converges in both, elements to `5.6e-16` against `1.8e-15`.

Coupled, on portlandite carrying sites at `Γ = 1e-5 mol/m²` over `90 m²/kg`, at
three host amounts:

| host | converged | `Σdn / n_host / ν − 1` | elements | host, coupled vs fixed |
|:--|:--|--:|--:|:--|
| 0.05 | yes, kkt `1.4e-14` | `3.9e-7` | `4.2e-15` | 0.0382614 / 0.0382737 |
| 0.10 | yes, kkt `5.7e-14` | `1.7e-7` | `1.5e-15` | 0.0882592 / 0.0882715 |
| 0.20 | yes, kkt `4.5e-11` | `8.0e-8` | `2.1e-12` | 0.1882548 / 0.188267 |

The fixed budget is off by 13 %, 6 % and 31 % on the same three, which is not a
defect of it: a budget that does not follow its host cannot track one.
"""
function conservation_matrix(cs::ChemicalSystem)
    A = Float64.(cs.SM.A)
    fams = cs.site_families
    fams === nothing && return A
    for f in fams
        surface_support(f).coupling === SITES_FOLLOW_HOST || continue
        r = findfirst(p -> get(atoms(p), f.site, 0) > 0, cs.SM.primaries)
        r === nothing && throw(
            ArgumentError(
                "SiteFamily \"$(name(f))\" follows its host, but no primary of this " *
                    "system carries :$(f.site), so there is no site row to couple.",
            )
        )
        prim = cs.SM.primaries[r]
        _is_bare_site(prim, f.site) || throw(
            ArgumentError(
                "SiteFamily \"$(name(f))\" follows its host, but its site primary is " *
                    "\"$(symbol(prim))\", which carries more than :$(f.site).\n" *
                    "A coupled family needs the BARE site as its component. " *
                    "Subtracting from the row of a primary that carries matter " *
                    "subtracts that matter too — measured, ν moles of oxygen and ν of " *
                    "hydrogen invented per mole of host, seven percent of the oxygen " *
                    "of Fe(OH)₃ at Dzombak and Morel's weak-site density.\n" *
                    "Declare `Species(\"$(f.site)\")` among the primaries. It need not " *
                    "be among the species: it is a component, not a substance.",
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
        _refuse_parasitic_charge(cs, f, r)
        A[r, j] -= sites_per_host(f, _molar_mass_si(cs.species[j]))
    end
    return A
end
