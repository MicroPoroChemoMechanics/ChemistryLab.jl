# [Surface areas](@id sec-manual-surfaces)

A reaction on a solid needs to know how much of that solid the liquid can reach.
This page is the syntax for saying so: the area models, what each one assumes,
and why the measurement an area comes from is part of its type.

## An area is a model, not a number

```@example surfaces
using ChemistryLab, DynamicQuantities

BETSurfaceArea(90.0u"m^2/kg")
```

Every model answers one question — [`total_area`](@ref), in m² — from the current
molar amount, the initial one, and the molar mass:

```@example surfaces
bet = BETSurfaceArea(90.0)           # a plain number is read as m²/kg
total_area(bet, 0.01, 0.05, 0.1)     # n = 0.01 mol, n₀ = 0.05 mol, M = 0.1 kg/mol
```

The initial amount is in the signature even where it is unused, so that a model
which does depend on it needs no change anywhere else. [`surface_area`](@ref) is
the two-amount-free shorthand, `total_area(model, n, n, M)`.

## The four models

| model | area | what it assumes |
|:--|:--|:--|
| [`FixedSurfaceArea`](@ref) | a total m², constant | the area is set externally and does not move over the run |
| [`BETSurfaceArea`](@ref) | `A_spec × n × M` | the specific area is constant, so area follows the remaining **mass** |
| [`BlaineSurfaceArea`](@ref) | same law, Blaine fineness | as above, for a binder characterized by air permeability |
| [`ShrinkingCoreArea`](@ref) | `A₀ (n/n₀)^p` | grains consumed from the outside, `p = 2/3` for spheres |

`BETSurfaceArea` and `ShrinkingCoreArea` with `p = 1` are the same law, which is
the point of writing the second as a family rather than a special case:

```@example surfaces
linear = ShrinkingCoreArea(bet; exponent = 1)
total_area(linear, 0.01, 0.05, 0.1) ≈ total_area(bet, 0.01, 0.05, 0.1)
```

At `p = 2/3` the area outlives the mass, because area goes as the square of a
radius and amount as its cube:

```@example surfaces
sphere = ShrinkingCoreArea(bet)                     # p = 2/3 by default
n₀ = 0.05
[(f, total_area(sphere, f * n₀, n₀, 0.1)) for f in (1.0, 0.5, 0.1)]
```

At exhaustion the power law has an infinite slope, which an integrator does
reach. The regularization is documented under [`SHRINK_FLOOR`](@ref): the area is
exactly zero at zero, exactly `A₀` at `n₀`, and differentiable throughout.

!!! note "What no model here does"
    None of them describes an area **blocked** by something growing on it. A
    hydrate layer covering a clinker grain reduces the accessible area without
    consuming the grain, so it depends on the whole assemblage rather than on one
    amount. That is a different object, and [`PoreHumidity`](@ref) is the
    precedent for writing one.

## Two areas of different measurements do not compare

This is the part worth reading even if the rest looks obvious. A BET area and a
Blaine fineness are different measurements of different things, and the package
now refuses to mix them:

```@example surfaces
try
    area_ratio(BETSurfaceArea(20_000.0), BlaineSurfaceArea(385.0))
catch err
    showerror(stdout, err)
end
```

The numbers in that example are silica fume's: about 20 000 m²/kg by BET, about
2 000 m²/kg as the effective Blaine fineness its kinetics were fitted with. Using
the first where the second belongs multiplies its hydration rate by ten, and
before this refusal nothing stopped it.

Same measurement, and the ratio is just a ratio:

```@example surfaces
area_ratio(BlaineSurfaceArea(462.0), BlaineSurfaceArea(385.0))
```

which is exactly what [`blaine_factor`](@ref) computes. A bare number handed to
it is still read as a Blaine fineness — that is the historical contract and it
has not changed — so type the argument when you want the check:

```@example surfaces
blaine_factor(462u"m^2/kg"), blaine_factor(BlaineSurfaceArea(462.0))
```

## A surface names its support

[`SurfaceSupport`](@ref) ties an area model to the solid carrying it:

```@example surfaces
SurfaceSupport("calcite", "Cal", BETSurfaceArea(90.0))
```

The rate factories take one in place of a bare area model, and that is what lets
the molar mass be looked up once from the named species. It matters: a species
with no `:M` used to fall back to 0.1 kg/mol silently, which scales the reactive
area and therefore the whole rate. It now raises, naming the species.

For a sorbent whose amount is prescribed rather than solved for, leave the host
out:

```@example surfaces
SurfaceSupport("inert sorbent", FixedSurfaceArea(0.5u"m^2"))
```

## See also

  - [Rate laws, and every parameter in them](@ref sec-theory-kinetics) — where the reactive area enters
    a rate law, and why the fineness correction of the empirical cement laws is a
    constant rather than an area.
  - [Surfaces API](@ref) — the docstrings.

## Site families: a conserved quantity that is not an element

A surface that only has an area does nothing. What makes it chemistry is a
finite number of **sites**, and species that occupy them.

The package already had a conserved quantity that is not a chemical element:
electric charge, carried as the pseudo-element `:Zz` in a formula and turned
into a conservation row by the ordinary matrix assembly. A site family works the
same way, with a pseudo-element of its own from [`SITE_SYMBOLS`](@ref):

```@example sites
using ChemistryLab, DynamicQuantities

surf(sym) = Species(sym; aggregate_state = AS_SURFACE, class = SC_SURFCOMPLEX)

free = surf("XsOH")          # the free site — a species like any other
prot = surf("XsOH2+")        # protonated
depr = surf("XsO-")          # deprotonated
atoms(free)
```

`:Xs` is the family; `O` and `H` are the real atoms. That matters for the mass:

```@example sites
free[:M]                     # oxygen and hydrogen only
```

**The site weighs nothing, by declaration.** The mass of the support is carried
once, by its own mineral species; adding it again on every occupied site would
be double counting. The same applies to volume — a surface complex has no
standard molar volume, because the volume it occupies is the host's.

### Declaring the family

```@example sites
support = SurfaceSupport("hydrous ferric oxide", nothing, FixedSurfaceArea(600.0u"m^2"))
family = SiteFamily("Hfo_s", free, [prot, depr];
                    capacity = TotalSiteAmount(5.0e-6u"mol"), support)
```

The capacity comes in three shapes because the published data does, and
converting between them needs a number nobody measured:
[`TotalSiteAmount`](@ref) for a prescribed sorbent, [`AreaSiteDensity`](@ref) in
mol/m² for the oxide literature, [`MassSiteDensity`](@ref) in mol/kg for a clay
exchange capacity — the last needing no area at all.

```@example sites
site_moles(family, 0.0, 0.0, 0.0)     # this one ignores the host entirely
```

### The site balance falls out of the matrix

Build the system with the **free site species as the primary**, and the
conservation row appears on its own:

```@example sites
aq(sym, cl = SC_AQSOLUTE) = Species(sym; aggregate_state = AS_AQUEOUS, class = cl)
species = [aq("H2O@", SC_AQSOLVENT), aq("H+"), aq("OH-"), aq("Ca+2"), free, prot, depr]
cs = ChemicalSystem(species, [species[1], species[2], species[4], free];
                    site_families = [family])
row = cs.SM.A[findfirst(p -> symbol(p) == "XsOH", cs.SM.primaries), :]
Int.(row)
```

One per surface species, zero on every aqueous one: that row **is**
`n_free + n_protonated + n_deprotonated = N_t`. And the rest of the matrix is
the protolysis, `XsOH2+ = XsOH + H+` and `XsO- = XsOH − H+`, written by the same
assembly that writes every other reaction.

!!! warning "With a fixed support, the primary is the free site"
    Declaring a **neutral** `Species("Xs")` as the site primary looks equivalent
    and is not. It is inconsistent in charge with a neutral free site, which
    flips the exact rank test inside [`StoichMatrix`](@ref) and adds a spurious
    charge component — silently.

    A support whose sites **follow their host** is the one case that needs a
    bare component, and it needs a *charged* one. See below: the package refuses
    the other two declarations by name rather than letting either through.

### What is refused, and why

Each of these produces a wrong number rather than an error if it is allowed
through, so each raises at construction: a member carrying no site symbol or two
of them; a free site occupying anything but one site; a species in two families;
two families sharing one pseudo-element; and an `AS_SURFACE` species belonging
to no declared family at all — it would carry a pseudo-element into the matrix,
and so a row, with nothing mixing on it.

Multidentate species are refused too, for now. The site balance holds for any
denticity — `Xs2OCa` contributes two — but ideal mixing of occupied and free
sites is exact only for one, and the quasi-chemical treatment of the rest is a
separate model this release does not provide. Saying so is better than returning
a number that looks like an isotherm.

### Where the species end up

`AS_SURFACE` is its own aggregate state, deliberately. Surface complexes are
**not** in `crystal(cs)`, because `idx_crystal` is what the solver and the start
repair read as "a pure mineral phase", which a site occupancy is not:

```@example sites
symbol.(surface(cs))
```

They are nonetheless counted in the **solid** compartment of a state, because
that is where their matter is.

## What mixes on the sites, and what Langmuir has to do with it

Declaring the family is half of it. The other half is the activity of its
members, and it is the simplest one there is: a **site fraction**.

```@example sites
grp = cs.site_groups[1]              # the family's members, free site first
n = zeros(length(cs.species))
n[1] = 55.5                          # water
n[grp] .= [5.0e-4, 3.0e-4, 2.0e-4]   # free, protonated, deprotonated
lna = activity_model(cs, DiluteSolutionModel())(n, (ϵ = 1.0e-30, T = 298.15))
[exp(lna[i]) for i in grp]           # each one is n_i / Σ n
```

Half the sites free, three tenths protonated, two tenths deprotonated — the
composition put in, read back as activities.

That is [`IdealSiteMixing`](@ref), and every activity model the package ships
applies it — dilute, HKF, Davies and Pitzer alike — because a site fraction owes
nothing to the ionic strength of the solution beside it.

### Langmuir is the consequence, not the premise

No isotherm is written anywhere in this package. Put the site balance
`n_free + Σ n_j = N` next to the mass-action law of each binding reaction, with
the free site's activity in it, and eliminate the free site:

```math
\frac{n_j}{N} = \frac{\beta_j}{1 + \sum_k \beta_k},
\qquad \beta_j = K_j\,a_j
```

— competitive Langmuir, with the **shared denominator** that is the signature of
a finite capacity. `test/surface_complexation.jl` checks it against that derived
form on an amphoteric surface across four pH values, to a relative 1e-6.

!!! warning "Do not apply a surface activity coefficient on top of this"
    The literature's surface activity coefficient, `γ_L = 1/(1 − θ)`, exists to
    make an *eliminated* free site reproduce what an explicit one already does.
    This formulation keeps the free site as a species, so the correction is
    already in the mixing. Applying it again counts saturation twice.

    The identity is worth knowing as a **diagnostic**, and the test suite asserts
    it rather than assuming it: at equilibrium `1/(1 − θ)` equals `N/n_free` to
    1e-8.

### In the solver, a family is a phase

A site family reaches the dual solver as a mixing phase, like the aqueous
solution and like a solid solution — never as a set of pure phases, which is
what it would be if it kept its default classification. Two settings differ from
a solid solution's, and neither is cosmetic:

  - the reference member is the **free site**, not the most abundant one. It is
    the state a fresh surface is mostly in, and a family cannot run out of it
    without the occupied states taking its place;
  - the phase is **always present**. A family's total is pinned by a
    conservation row, so a presence test on the size of its members would refuse
    to activate it from a cold start and leave `Σ n = N > 0` unsatisfiable from
    the first Newton iteration.

## Cation exchange: declaring the convention

A permanent-charge exchanger uses the same machinery with the sites counted as
**units of charge**. Write the divalent form with two of the family's
pseudo-elements and the conservation row becomes the cation exchange capacity:

```@example exchange
using ChemistryLab, DynamicQuantities

ex(sym) = Species(sym; aggregate_state = AS_SURFACE, class = SC_SURFCOMPLEX)
aq(sym, cl = SC_AQSOLUTE) = Species(sym; aggregate_state = AS_AQUEOUS, class = cl)

naX, kX, caX = ex("NaXc"), ex("KXc"), ex("CaXc2")
atoms(caX)          # two units of charge, from the formula
```

The convention is part of the declaration. There is no default, and nothing is
converted behind your back:

```@example exchange
support = SurfaceSupport("clay", nothing, FixedSurfaceArea(1.0))
family = SiteFamily("X", naX, [kX, caX];
                    capacity = TotalSiteAmount(1.0e-3u"mol"),   # moles of charge
                    support, model = GainesThomasMixing())
site_mixing_model(family)
```

```@example exchange
species = [aq("H2O@", SC_AQSOLVENT), aq("H+"), aq("Na+"), aq("K+"), aq("Ca+2"),
           naX, kX, caX]
cs = ChemicalSystem(species, species[[1, 2, 3, 4, 5, 6]]; site_families = [family])
Int.(cs.SM.A[findfirst(p -> symbol(p) == "NaXc", cs.SM.primaries), :])
```

One for each monovalent form, **two** for the calcium: the row counts
equivalents, and it counts them because the formula said so.

!!! warning "The two conventions are not a rescaling of each other"
    [`VanselowMixing`](@ref) takes the activity to be a mole fraction,
    [`GainesThomasMixing`](@ref) an equivalent fraction. On a homovalent
    exchange they are the same number. On a heterovalent one they are not, and
    the conversion between the constants fitted under each depends on the
    exchanger's composition — which is what you are solving for.

    Factors of 2, 3 or 4 appear in the literature; they are trace-composition
    limits. A constant fitted under one convention and used under the other is
    a different model.

There is no *free* site on an exchanger — every charge is compensated — so the
first member is a reference form rather than a vacancy.
[`reference_member`](@ref) says so where `free_site` would mislead.

Multidentate occupancy is refused under [`IdealSiteMixing`](@ref) and accepted
under either exchange convention. That is not an inconsistency: on an exchanger
the sites counted are charges, and a divalent cation neutralizes two of them
without straddling two surface groups. A surface complex that genuinely
straddles two sites is a combinatorial problem this release does not solve.


## When the support precipitates

Everything above poses the site budget once. A sorbent that appears during a
calculation — a C-S-H forming as a paste hydrates — carries its sites with it,
so the budget has to follow its amount instead.

Ask for it on the support, and give the family a capacity measured per unit mass
or per unit specific area, since those are the ones proportional to how much
host there is:

```julia
support = SurfaceSupport(
    "C-S-H", "CSHQ-JenD", BETSurfaceArea(90.0u"m^2/kg");
    coupling = SITES_FOLLOW_HOST,
)
family = SiteFamily("Xs", free, [bound]; capacity = AreaSiteDensity(1.0e-5), support)
```

Naming a host does **not** couple anything on its own: the kinetics has named
one since long before, to find the amount a rate law scales with. The coupling
is asked for.

`ν`, the moles of sites one mole of host carries, is then `Γ·a·M` — or `q·M` for
a mass density — and [`sites_per_host`](@ref) returns it. Which capacities
qualify is **measured**, not listed: it evaluates the capacity at two scaled
host amounts and requires the budget to scale with them. A
[`TotalSiteAmount`](@ref) and an area density over a [`FixedSurfaceArea`](@ref)
are constants and are refused; a [`ShrinkingCoreArea`](@ref) passes at `p = 1`
and nowhere else.

### The component is the bare site, carrying a charge

A coupled family needs its **component** — not its free site — as the primary,
and that component carries the charge the free site carries with its site
symbol. `XsOH` is `Xs⁺ + OH⁻`, so the component is `Species("Xs+")`; an
exchanger `NaXc` is `Xc⁻ + Na⁺`, so its component is negative.

```julia
bare = Species("Xs+"; aggregate_state = AS_SURFACE, class = SC_SURFCOMPLEX)
cs = ChemicalSystem(species, [h2o, hp, ca, bare]; site_families = [family])
```

It need not be among the species: it is a component, not a substance.

Both ways of getting this wrong are refused, by name and with what they cost.
Over the **free site**, subtracting the coupling would subtract that site's real
atoms too — seven percent of the oxygen of `Fe(OH)₃` invented, at Dzombak and
Morel's weak-site density. Over a **neutral** bare site, only the sum of the
site and charge potentials is identifiable, and the solver runs the two to
`±2.3e5` before stalling.

### The free site's reference energy stops being a gauge

With a **fixed** budget, `ΔₐG⁰` of the free site cancels out of every surface
reaction — both sides carry a site — so setting it to zero costs nothing.
Measured: shifting a whole family by up to 20 kJ/mol moves the host amount and
its saturation index by `5e-14`, which is the solver's own noise.

With a budget that **follows its host** it does not cancel. The host carries
`−ν` of the site component, so the site potential enters the host's own chemical
potential and its saturation index moves by

```math
\Delta \log \mathrm{SI}_{\text{host}} = -\,\frac{\nu\,\Delta}{RT \ln 10}
```

verified to five decimals over 40 kJ/mol of shift. At Dzombak and Morel's
weak-site density, `ν = 0.2`, that is **0.35 log units of solubility per
10 kJ/mol** of reference energy.

!!! warning "At a realistic site density this is a parameter, not a convention"
    Measured on amorphous ferric hydroxide carrying `ν = 0.2`: with
    `ΔₐG⁰(free site) = 0` the coupled host comes out **2.3 log units
    undersaturated and dissolves completely**, where the same system with a
    fixed budget keeps its solid at equilibrium. The coupling arithmetic is not
    at fault — the constraint holds to `4e-9` and the elements to `1e-14` — the
    reference energy is.

    The package does not choose that energy for you, and zero is not a safe
    default here. Until a standard state is fixed — Kulik's `Γ°` convention is
    the candidate, and [`convert_logk_site_density`](@ref) is the part of it
    that is implemented — treat a coupled family on a phase whose stability
    matters as **a model with one more parameter in it**, and check the host's
    saturation index against the same system uncoupled.

### What to check afterwards

[`site_budget_residual`](@ref) reports, per family, the moles of sites the state
carries minus the moles its capacity declares. Zero means the two say the same
thing. [`host_consistent_state`](@ref) derives the free-site amount from the
declaration so they do, and [`check_site_budget`](@ref) turns a disagreement into
an error instead of a silently different calculation.
