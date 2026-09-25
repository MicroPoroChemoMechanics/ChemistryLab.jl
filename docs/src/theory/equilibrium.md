# [Proving that an answer is the answer](@id sec-theory-certificate)

!!! info "Before this page"
    [Thermochemistry](@ref sec-theory-thermo) §4, where the minimization and its
    multipliers are set up, and [Standard states](@ref sec-theory-standard-states).

An interior-point method minimizes ``G`` by walking the interior of the feasible
set, and on a cement equilibrium it stops on `MaxIters` — at any tolerance.
Whether the point it returns is the minimum is then an open question. This page
is why the question has an answer at all, what a certificate checks, and how a
solver can aim at the optimality conditions directly instead of at the objective.

The formulation being certified is the one set out in
[Thermochemistry](@ref sec-theory-thermo) §4: minimize ``\sum_i n_i g_i(n)``
subject to ``\mathbf{A}n = b`` and ``n \ge 0``.

## Which equilibrium is computed

A system is said to be at equilibrium when none of its properties changes with
time and when, after a small disturbance that is later removed, it returns to the
same state [AndersonCrerar1993](@cite) (§3.3). The definition is operational,
and it is met by states that are not the minimum of the Gibbs energy: aragonite
kept at room temperature, or a mixture of hydrogen and oxygen, persists
indefinitely although calcite and water lie lower. Such a state is metastable,
separated from the stable one by an energy barrier that the conditions do not
allow the system to cross, whereas the stable equilibrium is the lowest state
compatible with the constraints.

A minimization knows nothing about barriers, and the state it returns is the
stable equilibrium of the system that was posed. Metastability therefore enters a
calculation only through the way the system is posed, and it does so in two ways.
A species absent from the list cannot form, so that leaving it out computes the
equilibrium that is metastable with respect to it. The first calculation of
[Getting started](@ref sec-quickstart) leaves out dissolved hydrogen, oxygen and
methane for this reason: with them, the minimization would also settle the
oxidation state of carbon, that is a redox equilibrium between carbonate and
methane which water at room temperature does not reach without a biological
catalyst, and which the example has no intention of describing.
[Choosing the species list](@ref man-choosing-species) treats the same decision
for cement phases, where a missing phase is a far more common error than a
deliberate exclusion. The second way is to withdraw some amounts from the budget
altogether, which is what the coupling with kinetics does with the phases it
integrates in time.

The latter is a partial equilibrium, in which only a subset of the possible
transformations has reached its balance while the others proceed at their own
pace. It is the whole principle of
the coupling between kinetics and equilibrium
([Coupling kinetics and equilibrium](@ref sec-coupling)), where the aqueous
speciation and the precipitation of hydrates are taken as instantaneous while
the dissolution of the clinker phases follows rate laws. Each instant of such a
trajectory is then the stable equilibrium of a smaller system, the one left once
the kinetic amounts have been withdrawn from the budget, and it is this
equilibrium that [`speciated_states`](@ref) recomputes and certifies. A last distinction concerns space.
A system out of equilibrium as a whole may be made of regions each at
equilibrium, which is called local equilibrium; every calculation of this package
is zero-dimensional and describes one such region, transport between regions
being outside its scope ([The water budget of a hydrating paste](@ref sec-theory-water-budget)
argues what this excludes for a paste).

The certificate below proves that a composition is the stable equilibrium of the
posed system. Whether the posed system, with its species list and its fixed
amounts, is the right one is a modeling question on which it has nothing to say.

## Why the question has an answer

Write the Gibbs energy in `RT` units as `G(n) = Σᵢ nᵢ μᵢ(n)`. Its ideal part

```math
\varphi(n) = \sum_i n_i \ln\frac{n_i}{N}, \qquad N = \sum_j n_j
```

has Hessian ``\operatorname{diag}(1/n_i) - \tfrac1N \mathbf{1}\mathbf{1}^\top``,
and for any ``v``

```math
v^\top \nabla^2\varphi\, v = \sum_i \frac{v_i^2}{n_i} - \frac{1}{N}\Bigl(\sum_i v_i\Bigr)^2 \;\ge\; 0
```

by Cauchy–Schwarz applied to ``v_i = (v_i/\sqrt{n_i})\sqrt{n_i}``. A pure phase
has unit activity, so it contributes a term **linear** in its amount. Hence `G`
is convex, the feasible set ``\{An=b,\ n\ge 0\}`` is a polyhedron, and — the
constraints being affine, so that the linearity constraint qualification holds
everywhere — the KKT conditions are **necessary and sufficient**.

Two consequences follow. The minimizer is unique, so a solver returning different
answers from different starting points is not finding local minima but stopping
short of stationarity. And optimality can be *checked*: a composition satisfying
the KKT conditions is proved globally optimal.

## The certificate

[`optimality_certificate`](@ref) checks the three conditions, on any composition
and whatever produced it. Writing ``u = -A^\top y`` for the element potentials:

| condition | on which species | meaning |
|:--|:--|:--|
| ``\mu_i + (A^\top y)_i = 0`` | interior (`n > floor`) | stationarity |
| ``An = b`` | — | conservation of matter |
| ``u_i \le g_i`` | a **pure** phase at its bound | that phase undersaturated |
| ``\ln \sum_i \exp(u_i - g_i - \ln\gamma_i) \le 0`` | a **mixing** phase held entirely absent | that solution cannot form |

The last row is Michelsen's tangent-plane measure, and it is a separate test
because a mixing phase needs one: its members are never exactly zero while it
exists, so they are neither interior nor at a bound, and a solid solution left out
of the assemblage used to pass the certificate **unexamined**. The trial
composition is refined against the phase's own activity model, so the test is not
the ideal approximation.

Two subtleties decide whether the check is meaningful.

A species **at its bound** obeys the inequality, not the equality. Imposing the
equality on an amount held at `1e-16` whose mass-action value is `e⁻³⁰⁰`
misstates its log-activity by 263 `RT` units, and the check then reports a
residual of 74 for a composition solved to `5e-12`.

A species carrying a **vanished component** is absent by the *constraint*, not by
thermodynamics, and its saturation index is meaningless — the element potential
of a component nobody supplies is determined by nothing. The test for that is not
`bₖ ≈ 0` but `bₖ ≈ 0` **with the non-zero entries of row `k` sharing a sign**:
only then does ``\sum_i A_{ki} n_i = 0`` with ``n \ge 0`` force each term to
vanish. The `H⁺` row carries `+1` for `H⁺` and `−1` for `OH⁻`, so its zero total
is the ordinary state of pure water; treating it as degenerate kills the entire
acid–base system and returns pH 7.000 with the calcite undissolved.

## Mass action or minimization

Programs computing a speciation fall into two families, according to the data
they take and the unknowns they solve for [AndersonCrerar1993](@cite) (§19.2).
The first writes, for every species outside a chosen basis, the law of mass
action of its formation reaction with a tabulated equilibrium constant,
completes these laws with the mass balances and the charge balance, and solves
the resulting nonlinear system by a Newton method on the activities of the basis
species; the phases allowed to precipitate are declared, and whether each of them
is present is decided from its saturation index. PHREEQC
[ParkhurstAppelo2013](@cite) is built on this formulation. The second minimizes
the Gibbs energy of the whole system under the mass balances, with the standard
potentials of all species as the only thermodynamic data, so that no reaction is
written and the assemblage is part of the result rather than of the input; it is
the formulation of GEM-Selektor [Kulik2013](@cite), of Reaktoro
[Leal2017](@cite) and of this package.

From the same data, both families compute the same state whenever both
converge, which follows from
the stationarity conditions of [Thermochemistry](@ref sec-theory-thermo) §4.
Taking as basis the primary species, whose potentials are the multipliers
``y_c``, the condition written for a present species ``s`` reads

```math
\ln a_s - \sum_c A_{cs}\ln a_c
  = -\frac{1}{RT}\Bigl(\Delta_a G_s^\circ - \sum_c A_{cs}\,\Delta_a G_c^\circ\Bigr)
  = \ln K_s ,
```

which is the law of mass action of the reaction forming ``s`` from the basis,
with the equilibrium constant implied by the standard potentials, while the
inequality written for an absent phase is the condition that its saturation index
be negative. What differs is what each family requires and what it guesses. A set
of equilibrium constants may be gathered reaction by reaction from separate
sources, whereas a minimization requires standard potentials consistent across
all species; conversely, a mass-action solver decides the presence of each
declared phase by a procedure added to its Newton iteration, whereas a
minimization decides the assemblage from the same conditions that define the
answer. The
convexity established above is a property of the minimization problem, and it is
what makes a certificate of global optimality possible at all.

## The certifying solver

[`DualEquilibriumSolver`](@ref) solves the KKT system directly, in element
potentials. From ``\mu_i + (A^\top y)_i = 0`` an aqueous species obeys the
mass-action law ``a_i = \exp(u_i - g_i)``, and a pure phase is present exactly
when ``u_i = g_i``, absent when undersaturated — the classical phase-stability
criterion.

Two levels. The inner one inverts the **solutes'** mass-action laws at fixed
potentials and fixed solvent amount; the outer is a Newton on ``1 + m + |P|``
unknowns — the solvent, the `m` element potentials, and the amounts of the
active phases. Parameterizing the solutes by ``\ln n`` makes their positivity
automatic, which is what removes the fraction-to-boundary limit that caps the
interior-point step at every iteration.

The solvent is deliberately **not** inverted through its own mass-action law: its
activity is a mole fraction, so ``\ln a_w \le 0`` always, and an arbitrary `y` can
demand more, for which no finite composition exists. It belongs to the outer
system, where the balance determines it.

!!! note "What it buys, measured"
    On calcite in pure water the certified pH is **9.90** against an
    interior-point 6.96 — not an imprecision but a wrong answer, and one nothing
    in that solver's output reveals. On the Reaktoro reference (calcite, CO₂ and
    water) both routes now agree with Reaktoro on every species: above `10⁻⁵` mol
    to `10⁻³` relative, the trace ions to 5 %, the worst being `CaOH⁺` at ×1.032
    on 1.6 nmol. That reference used to carry a `@test_broken` for `CaOH⁺` at
    ×2.47; what closed it was the convergence test moving to the true KKT error at
    `μ = 0` (`OptimaSolver` 0.4.1), and `test/equilibrium_reference.jl` is now 26
    plain assertions.

    [`speciated_states`](@ref) certifies every instant it replays and names any it
    cannot. On a full ordinary Portland cement over 28 days, all forty replayed
    instants are certified, with element balances between `1e-11` and `1e-13` mol.

## Where to go next

The chapter continues with [Activity models](@ref sec-theory-activity), the
first of the places where a mixture stops being ideal and the one on which every
aqueous equilibrium depends. The solver described here is driven in practice
from the tutorial [Chemical Equilibrium](@ref sec-equilibrium), and the API
entries are [`equilibrate_certified`](@ref), [`optimality_certificate`](@ref),
[`DualEquilibriumSolver`](@ref) and [`saturation_indices`](@ref).
