# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

# ── A site is a conserved quantity that is not an element ─────────────────────
#
# The package already has one of those: electric charge, carried as the
# pseudo-element `:Zz` in a species' formula, turned into a conservation row by
# the ordinary matrix assembly, and kept or dropped by an exact rank test. A
# surface site family works the same way, with its own pseudo-element from
# `SITE_SYMBOLS`. Nothing in `StoichMatrix` had to change for that.

"""
    abstract type AbstractSiteMixingModel end

How the species of a [`SiteFamily`](@ref) mix on their shared budget of sites.

Concrete subtypes implement

```julia
_site_excess_ln_gamma(model, k, x, T) -> Real
```

the departure from ideality of the `k`-th member, given the site fractions `x`
of the whole family. The ideal part, `ln x_k`, is added by the caller, exactly
as it is for a solid solution.

Only [`IdealSiteMixing`](@ref) exists in this release. The others the literature
uses — Frumkin's interaction term, the quasi-chemical approximation for a
multidentate adsorbate — are separate models with their own parameters and their
own validation, and this is where they will attach.

See also: [`IdealSiteMixing`](@ref), [`SiteFamily`](@ref).
"""
abstract type AbstractSiteMixingModel end

"""
    struct IdealSiteMixing <: AbstractSiteMixingModel

Occupied and free sites mix ideally: the activity of a member is its site
fraction, `a_j = n_j / N_t`.

# What this is, and what it is not

It is the production form of the Langmuir model, and Langmuir *falls out* of it
rather than being imposed: eliminating the free site from the equilibrium of
`n_free + Σ n_j = N_t` gives `n_j/N_t = β_j /(1 + Σ_k β_k)` with
`β_j = K_j a_j`, saturation and competition included, with no isotherm written
anywhere.

It is **not** a place to also apply a surface activity coefficient of the
`1/(1 − θ)` kind. That factor is what an eliminated-free-site formulation needs
in order to recover what this one already produces; applying both counts the
same physics twice. The ratio is available as a diagnostic — see the tests — and
the identity between the two forms is asserted rather than assumed.
"""
struct IdealSiteMixing <: AbstractSiteMixingModel end

"""
    supports_multidentate(model::AbstractSiteMixingModel) -> Bool

Whether `model` describes a species occupying more than one site.

`false` by default, and deliberately: the site **balance** holds for any
denticity — a species carrying two of the family's pseudo-elements contributes
twice to it — but a *mixing* law does not follow from the balance. Counting the
ways a molecule can straddle two neighboring sites is a combinatorial problem
with its own answer, and returning a number without having solved it would be
worse than refusing.

The ion-exchange conventions are the exception, and not because they solved that
problem: on a permanent-charge exchanger the "sites" being counted are **units
of charge**, which a divalent cation neutralizes two of without straddling
anything. Their multidentate case is bookkeeping, not combinatorics.

See also: [`VanselowMixing`](@ref), [`GainesThomasMixing`](@ref).
"""
supports_multidentate(::AbstractSiteMixingModel) = false

"""
    struct VanselowMixing <: AbstractSiteMixingModel

Cation exchange in the **Vanselow** convention: the activity of an exchanger
species is its **mole fraction**, counting molecules.

```math
a_i = \\frac{n_i}{\\sum_j n_j}
```

A calcium and a sodium on the exchanger count as one particle each, whatever
charge they neutralize.

See also: [`GainesThomasMixing`](@ref), and the note on conversion there.
"""
struct VanselowMixing <: AbstractSiteMixingModel end
supports_multidentate(::VanselowMixing) = true

"""
    struct GainesThomasMixing <: AbstractSiteMixingModel

Cation exchange in the **Gaines-Thomas** convention: the activity of an
exchanger species is its **equivalent fraction**, counting the charge it
compensates.

```math
a_i = \\frac{z_i \\, n_i}{\\sum_j z_j \\, n_j}
```

with `z_i` read from the formula as the number of the family's pseudo-elements
the species carries — which on a permanent-charge exchanger *is* the charge it
neutralizes.

# The two conventions are not interchangeable, and the conversion is not a factor

For a homovalent exchange, `Na⁺/K⁺`, the two fractions are proportional and the
selectivity coefficients agree. For a heterovalent one, `Na⁺/Ca²⁺`, they do not:
one calcium and one sodium are one particle each but one and two charges, so
`x_i` and `E_i` diverge, and so do the constants fitted under each.

The literature sometimes quotes conversion factors of 2, 3 or 4. Those are
**trace-composition limits**, not constant offsets: the exact relation depends
on the exchanger's composition, which is what the calculation is solving for.
This package therefore **declares** the convention and converts nothing
implicitly. A constant fitted under one convention used under the other is a
different model.

See also: [`VanselowMixing`](@ref), [`supports_multidentate`](@ref).
"""
struct GainesThomasMixing <: AbstractSiteMixingModel end
supports_multidentate(::GainesThomasMixing) = true

"""
    struct ConstantCapacitance{M<:AbstractSiteMixingModel, T<:Real} <: AbstractSiteMixingModel

A charged surface, in the **constant-capacitance** model: the potential is
proportional to the charge the surface carries.

```math
\\sigma = C\\, \\Psi, \\qquad
\\sigma = \\frac{F}{\\mathcal{A}} \\sum_k z_k n_k
```

with `C` the capacitance in F/m², `𝒜` the surface area in m², and `z_k` the
formal charge of each member. It **decorates** another mixing model rather than
replacing one: the site fractions are still whatever `base` says they are, and
this adds the electrical work of putting a charge on a charged surface.

# Why this needs no unknown of its own

The literature presents an electrostatic surface model as one extra unknown per
surface, `Ψ`, with one extra equation to close it. That is true of the diffuse
layer, where `Ψ` depends on the ionic strength and the closure cannot be
inverted. It is **not** true here: `σ = CΨ` makes `Ψ` an explicit function of
the composition,

```math
\\tilde\\psi \\equiv \\frac{F\\Psi}{RT}
 = \\frac{F^2}{C\\,\\mathcal{A}\\,RT} \\sum_k z_k n_k
```

so the whole model is a composition-dependent term in the chemical potential —
which is what an activity coefficient is. It belongs with the mixing, not with
the constraints, and putting it there is what keeps the solver unchanged.

# The convexity, written rather than assumed

The electrical work of charging the surface is
`G_el = ∫₀^σ Ψ(s)\\,ds \\cdot 𝒜`, and with `Ψ = σ/C` that integrates to

```math
G_{\\mathrm{el}}(n) = \\frac{F^2}{2\\,C\\,\\mathcal{A}} \\left(\\sum_k z_k n_k\\right)^2
```

a quadratic form in `n` with Hessian `\\frac{F^2}{C\\mathcal{A}} z z^{\\mathsf T}`,
positive semi-definite for any positive capacitance. So this term is **convex**,
the site mixing it decorates is convex, and the equilibrium certificate covers
the sum unchanged. That is a proof, not an expectation, and it is the reason
this model could be added without reopening the certificate.

Its gradient is `∂G_el/∂n_j = z_j F Ψ`, the electrochemical work the physics
asks for — which is the check that the energy above is the right one.

# How far the solve reaches, measured

Convexity makes the minimum unique, so any failure to find it is numerical and
not a second answer. And there is one: eliminating `Ψ` puts the whole
electrostatic stiffness into the composition dependence of an activity, and the
Newton loses it when that stiffness grows.

The scale of it is dimensionless and worth computing before a run:

```math
\\tilde\\psi_{\\max} = \\frac{F^2 N}{C\\,\\mathcal{A}\\,RT}
```

`N` being the site budget — the potential the surface would reach with every
site charged. Measured on hydrous ferric oxide, `N = 2·10⁻⁴` mol on 53.4 m²:

| `C` [F/m²] | `ψ̃_max` | stationarity of the solve |
|---:|---:|---:|
| 10 | 1.4 | 3e-16 |
| 5 | 2.8 | 2e-16 |
| 3 | 4.7 | 2e-16 |
| 2 | 7.0 | 0.04 — lost |
| 1.2 | 11.7 | 0.07 — lost |

So `ψ̃_max ≲ 5` is reached directly, and stepping `C` down from a large value
while reusing the previous composition extends the range but does not remove the
limit. Oxide capacitances of 1–3 F/m² therefore sit at the edge of it.

**What lifts it is the formulation, not the tolerance.** Carrying `Ψ` as an
unknown of the solve with `σ = CΨ` as its closing equation is mathematically the
same problem — that is what eliminating it proved — but the Newton then controls
the potential directly instead of meeting it through a stiff exponential. The
literature's extra unknown is a preconditioner, and it is what this model is
missing rather than a correction to it.

# Fields

  - `base`: the site mixing this decorates, usually [`IdealSiteMixing`](@ref).
  - `C`: capacitance [F/m²]. Values of 1–3 F/m² are the usual range for an oxide.
  - `area`: the charged area [m²].

# What it does to a titration

A surface that has already taken protons resists taking more, because the work
of adding a charge to an object that is already charged grows with the charge.
The visible effect is a **flattened** titration curve: the transitions spread
over more pH units than the constants alone would give. A set of constants
fitted *with* an electrostatic term and used *without* one — or the reverse —
therefore describes a different surface.

See also: [`IdealSiteMixing`](@ref), [`SiteFamily`](@ref).
"""
struct ConstantCapacitance{M <: AbstractSiteMixingModel, T <: Real} <:
    AbstractSiteMixingModel
    base::M
    C::T
    area::T
    function ConstantCapacitance{M, T}(base::AbstractSiteMixingModel, C::Real, area::Real) where {M <: AbstractSiteMixingModel, T <: Real}
        C > 0 || throw(ArgumentError("capacitance must be positive; got $C F/m²."))
        area > 0 || throw(ArgumentError("area must be positive; got $area m²."))
        return new{M, T}(base, convert(T, C), convert(T, area))
    end
end

"""
    ConstantCapacitance(base, C, area) -> ConstantCapacitance
    ConstantCapacitance(; C, area, base = IdealSiteMixing()) -> ConstantCapacitance

Build a [`ConstantCapacitance`](@ref). `C` is a capacitance in F/m² and `area`
an area in m², each a plain `Real` in SI or a `Quantity`.
"""
function ConstantCapacitance(
        base::AbstractSiteMixingModel, C, area
    )
    c = _area_si(us"F/m^2", C, "ConstantCapacitance capacitance")
    a = _area_si(us"m^2", area, "ConstantCapacitance area")
    v = promote(c, a)
    return ConstantCapacitance{typeof(base), eltype(v)}(base, v[1], v[2])
end

ConstantCapacitance(; C, area, base::AbstractSiteMixingModel = IdealSiteMixing()) =
    ConstantCapacitance(base, C, area)

# The decoration is transparent to everything the base model decides.
supports_multidentate(m::ConstantCapacitance) = supports_multidentate(m.base)

"""
    abstract type AbstractSiteCapacity end

How many moles of sites a family offers.

Concrete subtypes implement

```julia
site_moles(capacity, support, n_host::Real, n_host₀::Real, M_host::Real) -> Real
```

returning moles of sites.

# Why three of them, and not one

Because the published data comes in three shapes, and converting between them
needs a number nobody measured. A specific area with a site density per square
meter ([`AreaSiteDensity`](@ref)) is the oxide literature's form; a capacity per
kilogram of dry solid ([`MassSiteDensity`](@ref)) is the clay literature's, and
turning the second into the first would mean **inventing a BET area** to divide
by. A prescribed total ([`TotalSiteAmount`](@ref)) is what a fixed sorbent in a
batch experiment actually gives.

See also: [`SiteFamily`](@ref), [`site_moles`](@ref).
"""
abstract type AbstractSiteCapacity end

"""
    struct AreaSiteDensity{T<:Real} <: AbstractSiteCapacity

Sites per unit area, in mol/m². The area comes from the family's
[`SurfaceSupport`](@ref), so the measurement it was made with travels with it.

# Fields

  - `Γ_C`: site density [mol/m²].
"""
struct AreaSiteDensity{T <: Real} <: AbstractSiteCapacity
    Γ_C::T
    AreaSiteDensity{T}(Γ_C::Real) where {T <: Real} = new{T}(convert(T, Γ_C))
end

"""
    struct MassSiteDensity{T<:Real} <: AbstractSiteCapacity

Sites per kilogram of dry support, in mol/kg — the form a clay exchange capacity
is published in, and the one that needs no area at all.

# Fields

  - `q`: site capacity [mol/kg].
"""
struct MassSiteDensity{T <: Real} <: AbstractSiteCapacity
    q::T
    MassSiteDensity{T}(q::Real) where {T <: Real} = new{T}(convert(T, q))
end

"""
    struct TotalSiteAmount{T<:Real} <: AbstractSiteCapacity

A prescribed total number of moles of sites, independent of how much support is
present. The honest description of a batch experiment on a fixed sorbent, and
the form the first milestone's fixed-support scope uses.

# Fields

  - `N`: total sites [mol].
"""
struct TotalSiteAmount{T <: Real} <: AbstractSiteCapacity
    N::T
    TotalSiteAmount{T}(N::Real) where {T <: Real} = new{T}(convert(T, N))
end

# Outer constructors, with the unit check in one place. The inner constructors
# above exist so these are reachable at all — Julia's generated `Foo(x::T)` is
# more specific than an untyped outer one, and would bypass them.
for (S, unit, label) in (
        (:AreaSiteDensity, :(us"mol/m^2"), "AreaSiteDensity"),
        (:MassSiteDensity, :(us"mol/kg"), "MassSiteDensity"),
        (:TotalSiteAmount, :(us"mol"), "TotalSiteAmount"),
    )
    @eval function $S(value)
        v = _area_si($unit, value, $label)
        v >= zero(v) || throw(ArgumentError($label * " must be non-negative; got $value."))
        return $S{typeof(v)}(v)
    end
end

"""
    site_moles(capacity, support, n_host, n_host₀, M_host) -> Real

Moles of sites the family offers, given how much host solid is present.

`n_host` and `n_host₀` are the current and initial molar amounts of the support
species [mol] and `M_host` its molar mass [kg/mol]. A [`TotalSiteAmount`](@ref)
ignores all three; the other two do not, which is what will let a support that
precipitates carry its sites with it.

The first milestone holds the support fixed, so this is evaluated once. Writing
it as a function of the state rather than as a number is what makes an evolving
support a change of *when* it is called, not of the data model.
"""
function site_moles end

site_moles(c::TotalSiteAmount, ::SurfaceSupport, ::Real, ::Real, ::Real) = c.N

site_moles(c::MassSiteDensity, ::SurfaceSupport, n::Real, ::Real, M::Real) =
    c.q * max(n, zero(n)) * M

site_moles(c::AreaSiteDensity, s::SurfaceSupport, n::Real, n₀::Real, M::Real) =
    c.Γ_C * total_area(s.area, n, n₀, M)

# ── The family itself ─────────────────────────────────────────────────────────

"""
    struct SiteFamily{S<:AbstractSpecies, C<:AbstractSiteCapacity, U<:SurfaceSupport}

One family of surface sites: a group of species that share a finite site budget
and mix on it.

A family is to a surface what a [`SolidSolutionPhase`](@ref) is to a solid: a
named group of species with its own mixing. The one thing it has that a solid
solution does not is a **conservation row** — its total is fixed by a budget
rather than free — and that row is produced by the ordinary matrix assembly,
because every member carries the family's pseudo-element in its formula.

# The free site is a species, and it is the reference

`XsOH` is as much a member as `XsOH2+` is. It occupies a site, it takes part in
the mixing, and its amount is what makes saturation happen: as free sites run
out, adding another bound molecule costs more. Eliminating it gives the familiar
Langmuir expression, but the production form keeps it explicit — one species per
state of a site, with reaction-consistent standard potentials.

It is also the **primary** the system should be built on, and that is not a
detail of taste. A bare `Species("Xs")` declared as the site primary is
inconsistent in charge with a neutral free site, which flips the exact rank test
in `StoichMatrix` and silently adds a spurious charge component; and it never
enters `cs.species`, so `saturation_indices` would zero its potential and report
a wrong index for every surface species without a word.

# Denticity is read from the formula, not declared beside it

How many sites a molecule occupies is the coefficient of the family's symbol in
its formula — `Xs2OCa` occupies two. Declaring it separately would create a way
for the two to disagree. The first milestone **refuses** anything but one,
because ideal mixing of occupied and free sites is exact only for a monodentate
species; the quasi-chemical treatment of the rest is a later, separate model.

# Fields

  - `name`: the family's label, free-form (`"Hfo_s"`, `"Hfo_w"`).
  - `site`: its pseudo-element, one of [`SITE_SYMBOLS`](@ref).
  - `free_site`: the unoccupied site species.
  - `complexes`: the occupied ones.
  - `capacity`: how many moles of sites, as an [`AbstractSiteCapacity`](@ref).
  - `support`: the [`SurfaceSupport`](@ref) carrying it.
  - `model`: how the members mix on the budget, an
    [`AbstractSiteMixingModel`](@ref).

# Examples

```julia
support = SurfaceSupport("hydrous ferric oxide", "Fe(OH)3", BETSurfaceArea(600.0u"m^2/g"))
free    = Species("XsOH"; aggregate_state = AS_SURFACE, class = SC_SURFCOMPLEX)
bound   = Species("XsOH2+"; aggregate_state = AS_SURFACE, class = SC_SURFCOMPLEX)
SiteFamily("Hfo_s", free, [bound]; capacity = TotalSiteAmount(5.0e-6u"mol"), support)
```

See also: [`site_moles`](@ref), [`denticity`](@ref), [`SITE_SYMBOLS`](@ref).
"""
struct SiteFamily{
        S <: AbstractSpecies, C <: AbstractSiteCapacity, U <: SurfaceSupport,
        M <: AbstractSiteMixingModel,
    }
    name::String
    site::Symbol
    free_site::S
    complexes::Vector{S}
    capacity::C
    support::U
    model::M
end

"""
    denticity(family::SiteFamily, sp::AbstractSpecies) -> Int

How many sites of `family` one molecule of `sp` occupies, read from the
coefficient of the family's pseudo-element in its formula. Zero when `sp` does
not belong to the family.
"""
denticity(family::SiteFamily, sp::AbstractSpecies) =
    Int(get(atoms(sp), family.site, 0))

"""
    SiteFamily(name, free_site, complexes; capacity, support) -> SiteFamily

Build and validate a [`SiteFamily`](@ref).

Every member is requalified to `AS_SURFACE` / `SC_SURFCOMPLEX`, so a species
built from a database record can be passed as it is.

Refused at construction, each because the alternative is a wrong number rather
than an error:

  - a free site whose formula carries no site symbol, or two members carrying
    different ones — the family would have no budget to share;
  - a free site occupying anything other than one site;
  - a complex of denticity zero (it is not a member) or above one (ideal site
    mixing does not describe it; see the note on denticity above);
  - a duplicate member.
"""
function SiteFamily(
        name::AbstractString,
        free_site::AbstractSpecies,
        complexes::AbstractVector{<:AbstractSpecies} = AbstractSpecies[];
        capacity::AbstractSiteCapacity,
        support::SurfaceSupport,
        model::AbstractSiteMixingModel = IdealSiteMixing(),
    )
    site = _family_site_symbol(name, free_site, complexes)

    members = vcat([free_site], collect(complexes))
    syms = symbol.(members)
    if length(unique(syms)) != length(syms)
        dup = [s for s in unique(syms) if count(==(s), syms) > 1]
        throw(
            ArgumentError(
                "SiteFamily \"$name\": $(join(dup, ", ")) appears more than once. " *
                    "One species is one state of a site, and the site balance would " *
                    "count it twice.",
            )
        )
    end

    d_free = Int(get(atoms(free_site), site, 0))
    (d_free == 1 || (d_free > 1 && supports_multidentate(model))) || throw(
        ArgumentError(
            "SiteFamily \"$name\": the reference member \"$(symbol(free_site))\" " *
                "occupies $d_free sites; it must occupy at least one, and more than " *
                "one only under a model that describes multidentate occupancy.",
        )
    )

    for sp in complexes
        d = Int(get(atoms(sp), site, 0))
        d == 0 && throw(
            ArgumentError(
                "SiteFamily \"$name\": \"$(symbol(sp))\" carries no :$site, so it is " *
                    "not a member of this family.",
            )
        )
        (d == 1 || supports_multidentate(model)) || throw(
            ArgumentError(
                "SiteFamily \"$name\": \"$(symbol(sp))\" occupies $d sites, and " *
                    "$(nameof(typeof(model))) does not describe that. The site balance " *
                    "holds for any denticity, but a mixing law does not follow from " *
                    "the balance: ideal mixing of occupied and free sites is exact " *
                    "only for one site per molecule. An exchange convention " *
                    "(`VanselowMixing`, `GainesThomasMixing`) does handle it, because " *
                    "there the sites counted are units of charge; for a surface " *
                    "complex that genuinely straddles two sites, the quasi-chemical " *
                    "model this release does not provide is the one needed.",
            )
        )
    end

    qualified = [_as_surface_species(sp) for sp in members]
    return SiteFamily{
        eltype(qualified), typeof(capacity), typeof(support), typeof(model),
    }(
        String(name), site, first(qualified), qualified[2:end],
        capacity, support, model,
    )
end

# The one site symbol every member must agree on.
function _family_site_symbol(name, free_site, complexes)
    sites = Symbol[]
    for sp in vcat([free_site], collect(complexes))
        found = [k for k in keys(atoms(sp)) if is_site_symbol(k)]
        length(found) <= 1 || throw(
            ArgumentError(
                "SiteFamily \"$name\": \"$(symbol(sp))\" carries several site " *
                    "symbols ($(join(found, ", "))). A species occupies sites of one " *
                    "family; a bridge across two families is a different model.",
            )
        )
        append!(sites, found)
    end
    isempty(sites) && throw(
        ArgumentError(
            "SiteFamily \"$name\": no member carries a site symbol. A surface " *
                "species declares which family it belongs to through its formula, " *
                "for instance \"XsOH\"; see `SITE_SYMBOLS`.",
        )
    )
    length(unique(sites)) == 1 || throw(
        ArgumentError(
            "SiteFamily \"$name\": members carry different site symbols " *
                "($(join(unique(sites), ", "))). One family, one symbol.",
        )
    )
    return first(sites)
end

_as_surface_species(sp::AbstractSpecies) =
    with_class(with_aggregate_state(sp, AS_SURFACE), SC_SURFCOMPLEX)

"""
    name(family::SiteFamily) -> String
    site_members(family::SiteFamily) -> Vector{<:AbstractSpecies}
    site_capacity(family::SiteFamily) -> AbstractSiteCapacity
    surface_support(family::SiteFamily) -> SurfaceSupport

Accessors of a [`SiteFamily`](@ref). `members` lists the free site first, which
is the order the site mixing and the solver's reference member both expect.
"""
name(family::SiteFamily) = family.name

"""
    reference_member(family::SiteFamily) -> AbstractSpecies

The member the mixing is written against, and the one the solver carries instead
of inverting — `site_members(family)` puts it first for that reason.

On an oxide it is the **unoccupied** site, and `family.free_site` names it
literally: as the surface fills, its amount falls and every occupied state costs
more to form, which is where saturation comes from.

On a permanent-charge **exchanger** there is no unoccupied site at all: every
charge is compensated by some cation, and what this returns is the form chosen
as the reference of the exchange — usually the abundant monovalent one, `Na-X`.
The field keeps the name `free_site` from the case it was written for; this
accessor exists to say what it means in the case it was not.
"""
reference_member(family::SiteFamily) = family.free_site
site_members(family::SiteFamily) = vcat([family.free_site], family.complexes)
site_capacity(family::SiteFamily) = family.capacity
surface_support(family::SiteFamily) = family.support

"""
    site_mixing_model(family::SiteFamily) -> AbstractSiteMixingModel

How the members of `family` mix on their shared site budget.
"""
site_mixing_model(family::SiteFamily) = family.model

"""
    site_moles(family::SiteFamily, n_host, n_host₀, M_host) -> Real

The family's site budget, from its capacity and its support.
"""
site_moles(family::SiteFamily, n::Real, n₀::Real, M::Real) =
    site_moles(family.capacity, family.support, n, n₀, M)

function Base.show(io::IO, f::SiteFamily)
    return print(
        io,
        "SiteFamily(\"$(f.name)\", :$(f.site), $(length(f.complexes) + 1) members, ",
        "$(nameof(typeof(f.capacity))) on $(f.support.name))",
    )
end
