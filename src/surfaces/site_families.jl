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

# Examples

```julia
support = SurfaceSupport("hydrous ferric oxide", "Fe(OH)3", BETSurfaceArea(600.0u"m^2/g"))
free    = Species("XsOH"; aggregate_state = AS_SURFACE, class = SC_SURFCOMPLEX)
bound   = Species("XsOH2+"; aggregate_state = AS_SURFACE, class = SC_SURFCOMPLEX)
SiteFamily("Hfo_s", free, [bound]; capacity = TotalSiteAmount(5.0e-6u"mol"), support)
```

See also: [`site_moles`](@ref), [`denticity`](@ref), [`SITE_SYMBOLS`](@ref).
"""
struct SiteFamily{S <: AbstractSpecies, C <: AbstractSiteCapacity, U <: SurfaceSupport}
    name::String
    site::Symbol
    free_site::S
    complexes::Vector{S}
    capacity::C
    support::U
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
    d_free == 1 || throw(
        ArgumentError(
            "SiteFamily \"$name\": the free site \"$(symbol(free_site))\" occupies " *
                "$d_free sites; it must occupy exactly one.",
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
        d == 1 || throw(
            ArgumentError(
                "SiteFamily \"$name\": \"$(symbol(sp))\" occupies $d sites. The site " *
                    "balance holds for any denticity, but ideal mixing of occupied " *
                    "and free sites is exact only for one, so a multidentate species " *
                    "needs a quasi-chemical activity model this release does not " *
                    "provide. Declare it as monodentate, or wait for that model.",
            )
        )
    end

    qualified = [_as_surface_species(sp) for sp in members]
    return SiteFamily{eltype(qualified), typeof(capacity), typeof(support)}(
        String(name), site, first(qualified), qualified[2:end], capacity, support
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
site_members(family::SiteFamily) = vcat([family.free_site], family.complexes)
site_capacity(family::SiteFamily) = family.capacity
surface_support(family::SiteFamily) = family.support

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
