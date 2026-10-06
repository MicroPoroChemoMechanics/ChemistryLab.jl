# [Proving that an answer is the answer](@id sec-theory-certificate)

!!! info "Before this page"
    [Thermochemistry](@ref sec-theory-thermo) §4, where the minimization and its
    multipliers are set up, and [Standard states](@ref sec-theory-standard-states).

An interior-point method minimizes ``G`` by walking the interior of the feasible
set, and on a cement equilibrium it stops on `MaxIters` — at any tolerance.
Whether the point it returns is the minimum is then an open question. This page
is why the question has an answer at all, what a certificate checks, and how a
solver can aim at the optimality conditions directly instead of at the objective.

For a consistent Gibbs potential, the formulation being certified is the one
set out in
[Thermochemistry](@ref sec-theory-thermo) §4: minimize ``\sum_i n_i\, g_i(\mathbf{n})``
subject to ``\mathbf{A}\mathbf{n} = \mathbf{b}`` and ``\mathbf{n} \ge 0``.

## Which equilibrium is computed

A system is said to be at equilibrium when none of its properties changes with
time and when, after a small disturbance that is later removed, it returns to the
same state [AndersonCrerar1993](@cite) (§3.3). The definition is operational,
and it is met by states that are not the minimum of the Gibbs energy: aragonite
kept at room temperature, or a mixture of hydrogen and oxygen, persists
on the observation timescale although calcite and water lie lower. Such a state
is metastable,
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
([Kinetics under partial equilibrium](@ref sec-theory-pe-partition)), where the aqueous
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

The certificate below checks stationarity, conservation and phase conditions
for the posed system. These prove a stable global equilibrium when the model
defines a convex Gibbs potential, as qualified below. Whether the species list,
fixed amounts and boundary conditions describe the experiment is a separate
modeling question.

## Why the question has an answer

Write the Gibbs energy in `RT` units as `G(n) = Σᵢ nᵢ μᵢ(n)`. Its ideal part

```math
\varphi(\mathbf{n}) = \sum_i n_i \ln\frac{n_i}{N}, \qquad N = \sum_j n_j
```

has Hessian ``\operatorname{diag}(1/n_i) - \tfrac1N \mathbf{1}\mathbf{1}^\mathsf{T}``,
and for any ``\mathbf{v}``

```math
\mathbf{v}^\mathsf{T} \nabla^2\varphi\, \mathbf{v} = \sum_i \frac{v_i^2}{n_i} - \frac{1}{N}\Bigl(\sum_i v_i\Bigr)^2 \;\ge\; 0
```

by Cauchy–Schwarz applied to ``v_i = (v_i/\sqrt{n_i})\sqrt{n_i}``. A pure phase
has unit activity, so it contributes a term **linear** in its amount. For this
ideal formulation, `G` is convex, the feasible set
``\{\mathbf{A}\mathbf{n}=\mathbf{b},\ \mathbf{n}\ge 0\}`` is a polyhedron, and — the
constraints being affine, so that the linearity constraint qualification holds
everywhere — the KKT conditions are **necessary and sufficient**.

Optimality can then be *checked*: a composition satisfying the KKT conditions
is globally optimal. Convexity alone does not imply a unique composition: the
ideal Hessian has a scaling null direction, and pure-phase terms are linear.
Uniqueness requires strict convexity on the feasible directions or additional
conditions excluding degeneracy. Different starting points in a convex problem
can reveal incomplete convergence or degenerate global minima; they cannot
produce distinct isolated local minima of different energies.

### [What changes with a non-ideal model](@id sec-theory-convexity)

Two requirements must be checked separately. First, the supplied potentials
must be the gradient of one ``G``; the Maxwell and Gibbs-Duhem checks in
[Activity models](@ref sec-theory-potential) show why extended Debye-Hückel
formulas do not guarantee this. If no such potential exists, the certificate
checks equilibrium residuals and phase inequalities, not a minimum of a Gibbs
energy. Second, when a potential does exist, its **total** mixing energy must
be convex on the feasible set. Non-ideal solid solutions can violate this,
requiring phase-stability tests and phase separation
([Solid solutions](@ref sec-theory-solid-solutions)). A KKT point of a
non-convex potential alone need not be a minimum. The composition-dependent
water shift of [`CapillaryWater`](@ref) likewise lies outside the ideal proof.

These mathematical issues differ from kinetic metastability. An activation
barrier is a feature of a molecular or nucleation pathway, not a second minimum
that must appear in the bulk composition objective
([Metastable does not mean a local minimum of G](@ref sec-theory-metastable)
draws the two side by side). A thermodynamic solve does
not model barrier crossing. Suppressing a phase or holding kinetic amounts
fixed changes the feasible set, and a constrained minimum can be globally
optimal on that set while remaining metastable relative to an excluded
transformation. Waiting longer changes the constraints through kinetics; it
does not ask the equilibrium optimizer to climb an activation barrier.

## The certificate

[`optimality_certificate`](@ref) checks the conditions below, on any composition
and whatever produced it: stationarity, conservation, and a phase condition for
each kind of absent or vanishing species. Writing ``\mathbf{u} = -\mathbf{A}^\mathsf{T}\mathbf{y}`` for the potentials the element potentials give each species:

| condition | on which species | meaning |
|:--|:--|:--|
| ``\mu_i + (\mathbf{A}^\mathsf{T}\mathbf{y})_i = 0`` | interior (`n > floor`) | stationarity |
| ``\mathbf{A}\mathbf{n} = \mathbf{b}`` | — | conservation of matter |
| ``u_i \le g_i`` | a **pure** phase at its bound | that phase undersaturated |
| ``\mu_i + (\mathbf{A}^\mathsf{T}\mathbf{y})_i \ge 0`` | a member of a **present** phase below the floor, whose potential the interior determines | it holds no less than its potentials give it |
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
residual of 74 for a composition solved to `5e-12`. A member of a present phase
below the floor is read the same way, as the truncation of a smaller exact
amount: it may hold more than that amount, not less. Excluded from every test
until 0.28.0, such a member left behind by the search passed, as H⁺ at 3e-100 mol
did in a cement paste whose potentials gave it 1.2e-16.

A species carrying a **vanished component** is absent by the *constraint*, not by
thermodynamics, and its saturation index is meaningless — the element potential
of a component nobody supplies is determined by nothing. The test for that is not
`bₖ ≈ 0` but `bₖ ≈ 0` **with the non-zero entries of row `k` sharing a sign**:
only then does ``\sum_i A_{ki} n_i = 0`` with ``\mathbf{n} \ge 0`` force each term to
vanish. The `H⁺` row carries `+1` for `H⁺` and `−1` for `OH⁻`, so its zero total
is the ordinary state of pure water; treating it as degenerate kills the entire
acid–base system and returns pH 7.000 with the calcite undissolved.

### [What the certificate proves, and when](@id sec-theory-certificate-scope)

These conditions are sufficient for a global minimum only when the chemical
potentials are the gradient of one Gibbs energy and that energy is convex over
the feasible set. The first property can be measured at the audited composition,
since the second derivatives of one energy commute: the Jacobian of the log
activities has to be symmetric. It is for the ideal dilute model, whose solvent
row is the partner of the solutes' ``\ln m_i``, for Davies on ions, for the
Debye-Hückel form with a common ion size and no linear term, and for Pitzer's
equations; it is not for the extended forms in general use, B-dot with its ion
sizes and its linear term and Davies with a neutral solute, nor for SIT, nor for
a diffuse layer. With these, a certified equilibrium is a composition consistent
with its own activities rather than the minimum of an energy.

The second property is not a property of the point: a Hessian positive at the
answer proves a minimum there and nothing elsewhere. A global scope therefore
rests on a convexity established for every composition, phase by phase, the
Gibbs energy being the sum of the energies of the phases. A pure phase
contributes linearly; an ideal mixture (a gas, an ideal solid solution, ideal
site mixing, a constant capacitance, whose charging energy is a convex quadratic
of the charge) a convex function; a non-ideal solid solution what
[`mixing_convexity`](@ref) decides; and the aqueous phase what follows.

#### A convexity bound for the Debye-Hückel form

Let every ion carry ``\ln\gamma_i = -A' z_i^2 f(I)``, ``A' = \ln 10\,A``, with
one function ``f`` for all of them, ``f(I) = \sqrt{I}/(1+\kappa\sqrt{I})`` and
``\kappa = B\mathring{a}`` for the Debye-Hückel form with a common ion size, and
let the neutral solutes carry none. With ``w = n_w M_w`` the mass of solvent and
``I = \sum_s z_s^2 n_s/(2w)``, the sums over ``s`` running over the solutes,
these coefficients and the solvent row the Gibbs-Duhem relation pairs with them
derive from

```math
\frac{G}{RT} = \sum_i n_i \frac{\mu_i^\circ}{RT}
  + \sum_s n_s\Big(\ln\frac{n_s}{w} - 1\Big) + w\,h(I) ,
\qquad h' = -2A' f .
```

The ideal part has the quadratic form

```math
Q(\mathbf{v}) = \sum_s n_s\, u_s^2 ,
\qquad u_s = \frac{v_s}{n_s} - \frac{v_w}{n_w} ,
```

which is positive. The excess part ``E = w\,h(I)`` has the differential
``\mathrm{d}E = h M_w\,\mathrm{d}n_w + h'\,\mathbf{c}\cdot\mathrm{d}\mathbf{n}``,
with ``c_s = z_s^2/2`` and ``c_w = -I M_w``, so that
``w\,\mathrm{d}I = \mathbf{c}\cdot\mathrm{d}\mathbf{n}``. Differentiating once more,
the two terms in ``\mathrm{d}n_w\,\mathrm{d}I`` cancel and

```math
\mathrm{d}^2 E = \frac{h''}{w}\,(\mathbf{c}\cdot\mathbf{v})^2 ,
```

a term of rank one, negative wherever ``f`` increases. Substituting
``v_s = n_s(u_s + v_w/n_w)``, the solvent terms cancel,
``\mathbf{c}\cdot\mathbf{v} = \sum_s \tfrac{1}{2} z_s^2 n_s u_s``, and the
Cauchy-Schwarz inequality bounds it by the ideal part,

```math
(\mathbf{c}\cdot\mathbf{v})^2
\le Q(\mathbf{v}) \sum_s \frac{z_s^4 n_s}{4}
\le Q(\mathbf{v})\,\frac{z_{\max}^2}{2}\, w I .
```

The Hessian is therefore positive as soon as ``-h''(I)\,z_{\max}^2 I/2 \le 1``,
that is ``A' z_{\max}^2\, I f'(I) \le 1``. For the Debye-Hückel form,
``I f'(I) = \sqrt{I}/\big(2(1+\kappa\sqrt{I})^2\big)``, whose supremum, reached
at ``\sqrt{I} = 1/\kappa``, is ``1/(8\kappa)``, so that the condition

```math
\frac{\ln 10\; A\, z_{\max}^2}{8\, B\mathring{a}} \le 1
```

proves the aqueous phase convex whatever its composition. At 25 °C, with
``A = 0.5114`` and ``B = 0.3288``, it is 0.45 for divalent ions and an ion size
of 4 Å, and 1.01 for trivalent ones, which are therefore not proved; the limiting
law, without an ion size, never is. Davies' form has ``\kappa = 1`` and adds
``A' b\, z_i^2 I`` to ``\ln\gamma_i``, the derivative of
``w A' b I^2 = A' b\,(\sum_s z_s^2 n_s)^2/(4w)``, a convex function (a square
over a positive linear one), so that the same argument proves it convex when
``\ln 10\,A z_{\max}^2/8 \le 1``: 0.59 for divalent ions.

#### Where convexity is not proved

The certificate then asks what the answer is. The Hessian of ``G/RT`` is the
Jacobian of the log activities, symmetric here; restricted to the species
present, to the directions ``\mathbf{v}`` that conserve matter,
``\mathbf{A}\mathbf{v} = 0``, and scaled by ``\sqrt{n_i}`` so that an ideal solute
contributes exactly one whatever its amount, its smallest eigenvalue is reported
as `reduced_curvature`. Positive, the answer is a strict local minimum on its
active set, which a non-convex model may hold beside others. Negative, a
direction that conserves matter lowers the energy and the answer is a saddle: a
Pitzer set whose ``\beta^{(0)}`` is made strongly negative places a sodium chloride
solution saturated with halite there. Zero, the minimum is not strict.

The certificate reports which case applies as `scope`, with the reasons in
`scope_reasons`: `:global_minimum` when both properties hold and the constraint
is variational; `:local_minimum` when the energy exists, its convexity is not
proved and the reduced curvature is positive; `:kkt_point` when that curvature
is not positive, or when the constraint is not variational (a retention law that
shifts the activity of water, a temperature or a pressure made an unknown); and
`:self_consistent` when the activities are not a gradient. `optimal` keeps its
meaning in all four, the conditions being met to tolerance.

## [Mass action or minimization](@id sec-theory-mass-action)

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
Taking as basis the primary species, whose potentials are the negated
multipliers ``-y_c``, the condition written for a present species ``s`` reads

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
answer. When the model defines a convex Gibbs potential, the stationarity and
phase conditions also certify global optimality. Equality of mass-action
residuals alone does not establish the existence or convexity of that potential.

## The certifying solver

[`DualEquilibriumSolver`](@ref) solves the KKT system directly, in element
potentials. From ``\mu_i + (\mathbf{A}^\mathsf{T}\mathbf{y})_i = 0`` an aqueous species obeys the
mass-action law ``a_i = \exp(u_i - g_i)``, and a pure phase is present exactly
when ``u_i = g_i``, absent when undersaturated — the classical phase-stability
criterion.

Two levels. The inner one inverts the **solutes'** mass-action laws at fixed
potentials and fixed solvent amount; the outer is a Newton on
``|\Phi| + m + |P|`` unknowns — the total of each mixing phase present
(``\Phi``: the solvent of the aqueous phase, and each solid solution present),
the `m` element potentials, and the amounts of the active pure phases ``P``; a
constraint adds its own parameters. Parameterizing the solutes by ``\ln n_i``
makes their positivity automatic, which is what removes the fraction-to-boundary
limit that caps the interior-point step at every iteration.

A solid solution is a mixing phase without a solvent. At given potentials its
composition is explicit,

```math
x_i = \frac{\exp(u_i - g_i - \ln\gamma_i)}{\sum_j \exp(u_j - g_j - \ln\gamma_j)} ,
```

each end-member taking the larger share the cheaper it is against the potentials
of its elements, and the phase is present exactly when
``\sum_i \exp(u_i - g_i - \ln\gamma_i) = 1``, the equation its total adds to the
outer level. Which solid solutions are present is therefore not an input: an
absent one enters when that sum, measured with its activity coefficients at a
trial composition — the tangent-plane distance of [Michelsen1982](@citet) — says
it would lower the Gibbs energy, and leaves when its total vanishes, the number
of present phases bounded by the phase rule. Where ``\gamma_i`` depends on the
composition (non-ideal and sublattice models) the composition is a root, found
by Newton's method at each step.

For the Debye–Hückel, Davies and Truesdell–Jones models the inner level is one
equation. Their activity coefficients depend on the composition through the
ionic strength alone, so at a given ``I`` every solute is explicit,
``\ln m_i = u_i - g_i - \ln\gamma_i(I)``, and ``I`` solves
``I = \tfrac12\sum_i z_i^2 m_i(I)``, which is how PHREEQC carries the ionic
strength, as an unknown of its own [ParkhurstAppelo2013](@cite). The root taken is
the first one above the dilute limit, the branch connected to it, unless the
composition the solve holds lies above it on another branch. An ion of high
valence makes the equation cross zero three times: its activity coefficient
falls so fast with ``I`` that its amount, and ``I`` with it, outgrow ``I``. Under
the balances the equilibrium can then sit on the middle root, where an open
solution would be unstable — silica in a potassium hydroxide solution at 80 °C,
with the tetramer ``\mathrm{Si_4O_{10}^{4-}}`` carrying most of the dissolved
silicon, sits at 1.38 mol/kg between roots at 0.45 and 4.3 — and the first root
alone never gives it. The inversion therefore takes the root nearest the ionic
strength of the iterate when that lies above the first one, up to four times
the range the model states; with one root, that is the first one. The limiting
law keeps the first root: past its range it can have none there, the equation
dipping toward zero without reaching it and crossing again only through the
``\dot B I`` term, at thousands of mol/kg. The inversion then says that the
potentials hold no composition instead of iterating on them, and the outer
level moves on the sweeps until they hold one again. Recovering the solutes one
by one instead cycles where multivalent ions couple strongly through ``I``.

SIT and Pitzer add terms in the molalities themselves, the ``\varepsilon(i,k)\,m_k``
of the specific ion interaction and the pair and triplet sums of Pitzer, and no
single equation gives the solutes back. Their inner level is Newton's method on
the log-amounts of the solutes, ``h_i(\ln n) = u_i - g_i``, with the Jacobian of
the model itself, exact by forward differentiation. It starts from the better of
the composition the solve holds and the one the model's Debye–Hückel part gives
through the ionic strength, and steps at most 30 in any log-amount, halving the
step until the squared residual falls.

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
