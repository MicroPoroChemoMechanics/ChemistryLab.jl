# [Chemistry that happens on a surface](@id sec-theory-surface)

Most of this package treats a solid as a *bulk* phase: so many moles of calcite,
with a composition and a chemical potential, and nothing about its shape. That
is enough for dissolution and precipitation. It is not enough for the class of
phenomena where **what matters is the interface**: a proton that attaches to an
oxide's outer layer, an alkali that a C-S-H holds without ever forming a
separate phase, a chloride that a cement binds without becoming Friedel's salt.

This page explains what is added to describe those, in the order the ideas
depend on one another. No system is built and nothing is solved here; the
[syntax](@ref sec-manual-surfaces) and a
[worked case](@ref sec-example-surface-langmuir) are elsewhere.

## 1. The vocabulary, because the words are not interchangeable

A grain of mineral is the **support**. Some of the atomic groups exposed on its
outside can react with dissolved species; those are **binding sites**. Sites
that behave the same way and share one capacity form a **site family**. A
dissolved species that attaches to one becomes a **surface complex**.

The canonical reaction is not an exotic one — it is an acid-base reaction:

```math
\equiv\!\mathrm{XOH} + \mathrm{H}^{+} \;\rightleftharpoons\; \equiv\!\mathrm{XOH_2^{+}}
```

The `≡X` is not an element to look up in a periodic table: it means *attached to
the support*. Read the reaction carefully, because two things happen and one
does not:

  - the **charge** of the surface changes, and so does the distribution of
    hydrogen;
  - the **number of sites** does not. One site on the left, one on the right.

That second fact is the whole reason a surface needs a formalism of its own.
Together with its mirror image, `≡XOH ⇌ ≡XO⁻ + H⁺`, it is what makes an oxide
*amphoteric*: positively charged in acid, negatively charged in base, and
neutral somewhere between.

Two more words, because they are easy to conflate. **Adsorption** is
accumulation at an interface; **precipitation** is the appearance of a new
phase with its own stability. They compete for the same dissolved matter, which
is exactly why a surface belongs inside the general equilibrium calculation
rather than beside it. And **monodentate** means one site per adsorbed molecule,
**multidentate** several — a distinction that follows from how a site is
defined, and cannot be settled by a good curve fit.

## 2. A site is a conserved quantity that is not an element

Every calculation in this package rests on conservation: whatever redistributes,
the atoms are still there. That is the matrix equation

```math
A\,n = b
```

with one row per component and one column per species. Sites obey a law of
exactly the same shape — *every site is either free or occupied, and the total
is fixed while the support is* — but a site is not a chemical element, and its
row is not an element row.

The package already had one conserved quantity that is not an element: **electric
charge**. It rides in a species' formula as a placeholder symbol, and the matrix
assembly turns it into a row without knowing it is special. Sites are carried
the same way, each family with its own placeholder. For a single family with a
free site and two occupied states, the row reads

```math
n_{\equiv\mathrm{XOH}} + n_{\equiv\mathrm{XOH_2^{+}}} + n_{\equiv\mathrm{XO^{-}}} = N
```

with `N` the family's site budget in moles. Nothing was added to the matrix
machinery to obtain it.

!!! note "Where `N` comes from"
    Three ways, because published data comes in three shapes and converting
    between them needs a number nobody measured: a site density per unit area
    times an area, a capacity per kilogram of dry support times its mass, or a
    prescribed total. The second needs no area at all, which matters — turning a
    clay's exchange capacity into the first would mean **inventing** a BET area
    to divide by.

## 3. The activity of something bound to a surface

A chemical potential is always written the same way,

```math
\mu_j = \mu_j^\circ + RT \ln a_j
```

and the modeling is entirely in what `a_j` means. For a dissolved species it is
built from a concentration; for a pure solid it is one; for an end-member of a
solid solution it is a mole fraction.

For a species occupying a site it is a **site fraction**:

```math
a_j = \frac{n_j}{N}
```

its share of the family's budget. The free site is a species like any other and
has an activity of its own — and that is the part worth pausing on, because it
is where saturation comes from. As sites fill, `n` of the free site falls, its
activity falls with it, and the reaction that consumes it becomes harder to
drive. Nothing in the code says "stop at saturation"; saturation is what running
out of a reactant looks like.

!!! warning "Reference state and capacity are two different things"
    The literature distinguishes a **reference density** `Γ°`, which fixes where
    the activity scale has its zero, from a **capacity density** `Γ_C`, which is
    how many sites physically exist [Kulik2002](@cite). Changing the first means
    changing the standard potentials to match, and changes no physics. Changing
    the second changes the physics. Conflating them is the most common way a
    set of published constants stops meaning what it meant.

## 4. Langmuir is a consequence, not an ingredient

Here is the part that surprises people who have met adsorption through
isotherms. **No isotherm is written anywhere in this package.**

Take one family, a free site `≡X`, and adsorbates `j` binding one site each.
Each binding reaction has its law of mass action,

```math
K_j = \frac{a_{\equiv\mathrm{X}j}}{a_{\equiv\mathrm{X}}\, a_j}
      = \frac{n_{\equiv\mathrm{X}j}}{n_{\equiv\mathrm{X}}\, a_j}
```

— the two `N` cancel, which is why the result does not depend on how the budget
was specified. Write `β_j = K_j a_j`, the *drive* of adsorbate `j` at the
dissolved activity it actually has. Then `n_{≡Xj} = β_j n_{≡X}`, and the site
balance `n_{≡X} + Σ_j n_{≡Xj} = N` closes it:

```math
\boxed{\;\frac{n_{\equiv\mathrm{X}j}}{N} = \frac{\beta_j}{1 + \sum_k \beta_k},
\qquad
\frac{n_{\equiv\mathrm{X}}}{N} = \frac{1}{1 + \sum_k \beta_k}\;}
```

That is the competitive Langmuir isotherm. It was not assumed; it fell out of a
conservation law and a mass-action law, which is all the package ever applies.

Read the denominator, because it carries the physics:

  - a weak drive, `β ≪ 1`, leaves the sites mostly free and gives
    `n_j/N ≈ β_j` — the linear regime, where sorption looks like a partition
    coefficient;
  - a strong drive, `β ≫ 1`, gives `n_j/N → 1`: saturation;
  - **the sum is shared**. Every adsorbate appears in every other's denominator,
    because they are drawing on one finite capacity. That shared denominator
    *is* competition, and it is the thing no single-species isotherm can express.

## 5. The surface activity coefficient, and why it is absent here

The literature often writes the same physics with the free site **eliminated**,
carrying only the occupied species. That formulation then needs a correction
factor to reproduce what the free site was doing, and [Kulik2006](@cite) derives
it: for monodentate species on one family,

```math
\gamma_L = \frac{1}{1 - \theta_\Sigma},
\qquad \theta_\Sigma = \frac{\sum_j n_j}{N}
```

with `θ_Σ` the occupied fraction. It grows without bound as the sites fill,
which is exactly the resistance to further binding that an explicit free site
provides by simply becoming scarce.

This package keeps the free site, so **`γ_L` is not applied**. Applying it on top
would count saturation twice. The identity between the two formulations is worth
knowing, and it is verified rather than asserted: at equilibrium, `1/(1 − θ_Σ)`
equals `N` divided by the free-site amount, which is the same statement written
twice.

!!! note "What about several sites per molecule?"
    The site balance holds for any denticity — a bidentate species contributes
    twice to it. Ideal *mixing* of occupied and free sites does not: counting the
    ways a molecule can straddle two neighboring sites is a different
    combinatorial problem, treated by the quasi-chemical approximation
    [Kulik2006](@cite), whose own derivation excludes mixed adsorbates. This
    release therefore **refuses** a multidentate species rather than returning a
    number that looks like an isotherm.

## 6. What the certificate still means

This package does not only solve; it [proves](@ref sec-theory-certificate) that
what it returns is a constrained minimum. A new kind of species could break that
proof, so it is worth saying why this one does not.

Ideal site mixing contributes, per family,

```math
G_{\text{sites}} = \sum_j n_j \bar\mu_j^\circ
   + RT \sum_j n_j \ln \frac{n_j}{N}
```

and the entropy term `Σ n ln(n/Σn)` is convex — the same function, and the same
argument, as an ideal solid solution. The site balance is linear. So the
surface adds a convex term under a linear constraint, and the certificate covers
it unchanged.

That is true of *this* model and is not a general statement about surfaces. Two
extensions on the roadmap will need the question reopened rather than inherited:
an attractive lateral-interaction term can make the mixing energy non-convex
near half coverage, and an electrostatic term depends on the ionic strength,
hence on the aqueous composition, so convexity at fixed ionic strength is not
convexity in the composition. Where that fails, the honest outcome is a
certificate that **refuses**, as it already does for a concave solid solution —
not one that quietly means less.

## 7. What this page does not cover

Saying what is absent is part of describing what is present.

  - **No surface potential.** Real surfaces are charged, and that charge
    attracts counter-ions and repels co-ions, changing what binds. The models
    that describe it — constant capacitance, the diffuse double layer, the
    Donnan approximation — add one unknown per surface and a closure relating
    charge to potential. None is in this release, so a set of constants fitted
    *with* an electrostatic term must not be used here: it would be a different
    model wearing the same numbers.
  - **No evolving support.** The site budget is fixed. In a hydrating cement the
    support is a phase that precipitates, so its sites appear with it and the
    site row becomes bilinear — the one thing the linear budget `A n = b` has
    never had to carry.
  - **No multidentate species**, per the note above.
  - **No ion exchange conventions.** Cation exchange on a permanent-charge
    surface is close kin to what is here, but its capacity is counted in moles
    of *charge* and its selectivity coefficients come in two conventions,
    Vanselow and Gaines-Thomas, that are not interchangeable
    [Marinich2025](@cite).

## See also

  - [Surface areas](@ref sec-manual-surfaces) — the syntax, and how a site
    family is declared.
  - [Adsorption on a single site family](@ref sec-example-surface-langmuir) —
    the smallest complete case, computed, against the closed form above.
  - [Proving that an answer is the answer](@ref sec-theory-certificate) — what
    the certificate checks, and what it does not.
  - [Solid solutions](@ref sec-theory-solid-solutions) — the other mixing phase, whose arithmetic is the
    same and whose closure is not.
