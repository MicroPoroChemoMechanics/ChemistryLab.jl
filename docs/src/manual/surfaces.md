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

!!! warning "The primary is the free site, never a bare site species"
    Declaring `Species("Xs")` as the site primary looks equivalent and is not. A
    neutral bare site is inconsistent in charge with a neutral free site, which
    flips the exact rank test inside [`StoichMatrix`](@ref) and adds a spurious
    charge component — silently. It also never enters `cs.species`, so
    [`saturation_indices`](@ref) would zero its potential and report a wrong
    index for every surface species, again without a word.

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
