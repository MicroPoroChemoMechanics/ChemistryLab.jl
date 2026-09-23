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

## 7. Cation exchange: the same machinery, counting charge

A permanent-charge clay is not an oxide. Its charge does not come and go with
pH — it is built into the mineral, from substitutions in the lattice — and
there is no such thing as an unoccupied site: every unit of charge is
compensated by some cation, and what varies is *which*.

That sounds like a different model. It is the same one, with the sites counted
differently.

### The budget is in charge, and it says so in the formula

Write the exchanger species as `Na-X` and `Ca-X₂`. The second carries **two** of
the family's pseudo-elements because a calcium neutralizes two units of charge,
and the conservation row that comes out of the matrix assembly is therefore

```math
n_{\mathrm{Na\text{-}X}} + n_{\mathrm{K\text{-}X}} + 2\,n_{\mathrm{Ca\text{-}X_2}} = \mathrm{CEC}
```

— the **cation exchange capacity**, in moles of charge. Nothing was declared to
make that happen: the coefficient is in the formula, and the row counts what
the formula says.

There is also no free site to serve as the reference of the mixing, so one of
the forms is chosen as the reference — the abundant monovalent one, usually.

### Two conventions, and no factor between them

Here the two literatures part company, and the package refuses to guess which
one a number came from.

| convention | activity of an exchanger species | counts |
|:--|:--|:--|
| **Vanselow** | ``x_i = n_i / \sum_j n_j`` | particles |
| **Gaines-Thomas** | ``E_i = z_i n_i / \sum_j z_j n_j`` | charges |

One calcium and one sodium are **one particle each** and **two charges against
one**. For a homovalent exchange, `Na⁺/K⁺`, every ``z`` is 1 and the two
fractions are the same number — which is why an oxide surface never has to
choose. For a heterovalent one they diverge, and so do the selectivity
coefficients fitted under each.

The literature quotes conversion factors of 2, 3 or 4 between them. Those are
**trace-composition limits**, not constant offsets: the exact relation depends
on the exchanger's composition, which is what the calculation is solving for
[Marinich2025](@cite). So this package converts nothing implicitly. The
convention is part of the declaration, and a constant fitted under one and used
under the other is a different model — not a rescaled one.

Both conventions are checked against Reaktoro, which carries both, on a
heterovalent Na/K/Ca exchange, and agree to 4 parts in 10⁹. On that composition
they differ from **each other** by 48 % on the sodium, which is the point.

## 8. A charged surface, and why it needs no new unknown

Everything above ignored one thing: a surface that binds protons **becomes
charged**, and the work of putting one more charge on an already charged object
is not zero. The chemical potential of a member acquires an electrical term,

```math
\mu_j = \mu_j^\circ + RT\ln a_j + z_j F \Psi
```

with `Ψ` the potential of the surface plane and `z_j` the species' formal
charge. What relates `Ψ` to the charge the surface carries is a **closure**, and
the simplest is a capacitor: ``\sigma = C\,\Psi``.

### The elimination, and what it buys

The literature presents an electrostatic surface model as one extra unknown per
surface with one extra equation. For the diffuse layer that is unavoidable: `Ψ`
there depends on the ionic strength through a relation with no closed inverse.
For a **constant capacitance** it is not. The closure inverts, so

```math
\tilde\psi \equiv \frac{F\Psi}{RT}
 = \frac{F^2}{C\,\mathcal{A}\,RT}\sum_k z_k n_k
```

is an explicit function of the composition — and a composition-dependent term in
a chemical potential is exactly what an activity coefficient is. So the model
belongs with the mixing, and the solver needs no new machinery at all.

### The certificate survives it, and that is a proof

The electrical work of charging the surface is ``\int_0^\sigma \Psi\,ds``
over the area, which with ``\Psi = \sigma/C`` integrates to

```math
G_{\mathrm{el}}(n) = \frac{F^2}{2\,C\,\mathcal{A}}\left(\sum_k z_k n_k\right)^2
```

a quadratic form whose Hessian ``(F^2/C\mathcal{A})\,z z^{\mathsf T}`` is
positive semi-definite for any positive capacitance. Convex term, convex mixing,
linear constraint: the certificate of §6 covers the sum unchanged. Its gradient
is ``z_j F \Psi``, which is how one knows the energy is the right one.

That is *this* model. It does not transfer: the diffuse layer's `Ψ` depends on
the ionic strength, hence on the aqueous composition, so convexity at fixed
ionic strength is not convexity in the composition, and that question has to be
reopened rather than inherited.

### What it does, and where it stops

A surface already charged resists charging further, so a titration curve
**flattens**: the transitions spread over more pH units than the constants alone
would give. On hydrous ferric oxide at pH 5, protonation falls from 0.995 without
the term to 0.81 at `C = 3 F/m²`.

Convexity makes the minimum unique, so any failure to find it is numerical.
There is one, and its scale is

```math
\tilde\psi_{\max} = \frac{F^2 N}{C\,\mathcal{A}\,RT}
```

the potential the surface would reach fully charged. Below about 5 the solve is
found to machine precision; above it the Newton loses it. On that same oxide
that is `C ≳ 3 F/m²`, which puts the usual oxide range of 1–3 at the edge.

**The remedy is the formulation, not the tolerance.** Carrying `Ψ` as an unknown
with ``\sigma = C\Psi`` as its equation is the same problem — the elimination
proved that — but the Newton then controls the potential directly instead of
meeting it through a stiff exponential. The extra unknown of the textbooks is a
preconditioner.

## 9. The diffuse layer, and the price of eliminating a potential

§8 removed an unknown by inverting a closure. The diffuse layer is where that
trick is usually said to stop, and it is worth being precise about what stops
and what does not.

### The closure inverts too

Gouy and Chapman's solution of the Poisson-Boltzmann equation beside a flat
surface, for a symmetric 1:1 electrolyte, relates the charge the surface carries
to the potential it raises and to how well the solution screens it:

```math
\sigma \;=\; \kappa\sqrt{I}\;\sinh\!\left(\frac{F\Psi}{2RT}\right),
\qquad \kappa = \sqrt{8\,\varepsilon_r\varepsilon_0 RT\rho}
```

with ``I`` the ionic strength of the bulk solution and ``\rho = 1000`` kg/m³ the
factor that makes it a volumetric concentration. This is transcendental in
``\Psi`` — and *monotone* in it, which is a different thing. A monotone relation
inverts, and this one inverts in closed form:

```math
\tilde\psi \;=\; 2\,\operatorname{asinh}\!\left(\frac{\sigma}{\kappa\sqrt{I}}\right)
```

`asinh` is smooth, its derivative is bounded by one, and it costs the solver a
single call. So the diffuse layer, like the constant capacitance, is an
*activity coefficient*: `ln a_k = ln x_k + z_k ψ̃`, with the same convention and
a different closure. No unknown is added.

### What is actually lost is not the unknown

Two different questions are usually merged into one word, "consistency", and
this is the model that separates them.

| | is it a gradient? | is that gradient extensive? |
|:--|:--|:--|
| ideal site mixing | yes | yes |
| constant capacitance | **yes** | no |
| diffuse layer | **no** | no |

The first column is the one the certificate needs: a minimization has to be
minimizing *something*. Symmetry of the activity Jacobian,
``\partial \ln a_i/\partial n_j = \partial \ln a_j/\partial n_i``, is exactly
the statement that such a something exists.

A constant capacitance passes it — its charging work
``F^2(z\cdot n)^2/(2C\mathcal{A})`` is a genuine potential — and fails the
second column, because that expression is homogeneous of degree two in the
amounts while a Gibbs energy is homogeneous of degree one. That failure is not a
defect. The area is a *parameter* of the problem, fixed from outside like a
volume; scaling the amounts without scaling the surface is not the extensive
scaling Gibbs-Duhem is about.

A diffuse layer fails the first column, and that one is a defect. ``\tilde\psi``
depends on ``I``, which depends on the aqueous composition, while the activity
of an aqueous ion does not depend in return on what is bound to the surface. The
coupling is one-way, so the Jacobian is asymmetric, and an asymmetric Jacobian
is the Hessian of nothing.

This is a property of the model, not of any implementation of it. Treating the
bulk as a reservoir whose ionic strength is a parameter is precisely the
approximation that lets a diffuse layer be written without carrying its own
inventory of counter-ions — the approximation Dzombak and Morel make, and the
one PHREEQC's default `SURFACE` block makes. What comes back from such a solve
is a **self-consistent speciation**, mass action and conservation satisfied
together. It is not a certified minimum, and this package says which it is
rather than letting the word "certificate" cover both.

### And a second price, which is the one that bites

Writing a potential as an activity coefficient means the solver reaches it by
successive substitution: a composition implies activities, which imply a
composition. That iteration converges while the activity moves less than the
composition does — while the dimensionless sensitivity

```math
N\left|\frac{\partial\tilde\psi}{\partial n}\right|
\;=\;\frac{2NF}{\mathcal{A}\,\kappa\sqrt{I}\,\sqrt{1+u^2}},
\qquad u = \frac{\sigma}{\kappa\sqrt{I}}
```

is small. Above roughly five it is not, and what the solve returns is not a
second root but a point that violates mass action outright.

Three things follow, and all three are measured rather than argued:

  - **The physics is never the problem.** Charging a surface always opposes
    further charging, so the equilibrium is unique and stable at every
    composition. Only the elimination fails.
  - **A diffuse layer is hard in the middle and easy at the edges.** Far from
    the point of zero charge ``u`` is large, `asinh` flattens, and the same
    system that will not solve at pH 7 solves at pH 4. It is also harder in a
    *dilute* background, as ``1/\sqrt{I}``, which is the opposite of the usual
    intuition about difficult chemistry.
  - **A constant capacitance is judged once.** Its sensitivity,
    ``NF^2/(C\mathcal{A}RT)``, does not depend on the composition at all.

The forecast is [`electrostatic_stiffness`](@ref) and the threshold
[`ELECTROSTATIC_STIFFNESS_LIMIT`](@ref); the verdict, as always, is the
stationarity residual of the certificate, which separates the two regimes by
fourteen orders of magnitude. The honest resolution — carrying ``\Psi`` as an
unknown with its own equation, so the coupling is linearized instead of
iterated — is not in this package yet, and the page says so rather than leaving
a reader to discover it at pH 7.

## 10. What this page does not cover

Saying what is absent is part of describing what is present.

  - **No Donnan approximation, and no diffuse-layer inventory.** The diffuse
    layer is here (§9), in the form that carries a potential and not an ion
    census: the counter-ions accumulated in the layer are not tracked as a
    separate reservoir, which is PHREEQC's default and not its `-diffuse_layer`
    option. The Donnan approximation, which a compacted clay needs, is absent.
  - **No surface potential beyond the point where eliminating it works.** §9
    measures where that is. Above it the answer is refused by the certificate
    rather than returned.
  - **No evolving support.** The site budget is fixed. In a hydrating cement the
    support is a phase that precipitates, so its sites appear with it and the
    site row becomes bilinear — the one thing the linear budget `A n = b` has
    never had to carry.
  - **No multidentate species**, per the note above.
  - **No lateral interactions.** Neighbors on a surface affect each other's
    binding energy, and the models that describe it — Frumkin's interaction
    term, the quasi-chemical approximation — can make the mixing energy
    non-convex. Neither is here.

## See also

  - [Surface areas](@ref sec-manual-surfaces) — the syntax, and how a site
    family is declared.
  - [Adsorption on a single site family](@ref sec-example-surface-langmuir) —
    the smallest complete case, computed, against the closed form above.
  - [Proving that an answer is the answer](@ref sec-theory-certificate) — what
    the certificate checks, and what it does not.
  - [Solid solutions](@ref sec-theory-solid-solutions) — the other mixing phase, whose arithmetic is the
    same and whose closure is not.
