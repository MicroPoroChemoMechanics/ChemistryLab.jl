# [Solid solutions](@id sec-theory-solid-solutions)

!!! info "Before this page"
    [Thermochemistry](@ref sec-theory-thermo) and [Standard states](@ref
    sec-theory-standard-states) §2, where the standard state of an end-member is
    set.

A pure solid has ``a = 1``: it is its own standard state, its chemical potential
does not depend on how much of it there is, and it either precipitates or it does
not. A **solid solution** is a single crystalline phase whose composition varies
continuously between end-members, and that changes the thermodynamics
qualitatively: the activity of an end-member depends on the composition of the
phase, so the phase can absorb a component gradually instead of appearing all at
once.

This matters in cement more than anywhere else. C-S-H is not a compound with a
formula but a phase of variable Ca/Si; the AFm phases exchange sulfate, carbonate
and hydroxide on the same interlayer site; alkalis partition into hydrates rather
than staying in solution. A model that only knows pure phases predicts sharp
appearances and disappearances that are not observed.

## 1. Ideal mixing: the activity *is* the mole fraction

Put ``n_k`` moles of each end-member into one phase, with mole fractions
``x_k = n_k/\sum_j n_j``. Even with no interaction at all, mixing lowers the
Gibbs energy, because there are more ways to arrange a mixture than a pure
solid. Per mole of phase,

```math
G^{\text{mix}} = RT\sum_k x_k \ln x_k \;\; (<0) ,
```

which is entropy alone — the ``-T\Delta S`` of mixing. Differentiating for one
end-member gives its chemical potential and hence its activity:

```math
\mu_k = \mu_k^\circ + RT\ln x_k
\qquad\Longrightarrow\qquad
\boxed{\;\ln a_k = \ln x_k\;}
```

That is [`IdealSolidSolutionModel`](@ref), and it is the default. Note what it
already buys: an end-member at ``x_k = 10^{-3}`` has ``\ln a_k = -6.9``, so its
saturation index is shifted by ``-3`` decades relative to the pure phase. A trace
component is *stabilized* by being diluted in a host, which is why a solid
solution can take up an ion that would never precipitate as its own phase.

### The curve, its tangent, and the plane of the rest of the system

![Left: the Gibbs energy of an ideal binary solid solution below the line of its two end-members unmixed, and the tangent at one composition meeting the two edges at the chemical potentials of the end-members. Right: the plane the rest of the system imposes, against three curves of a solid solution: above it, absent; touching it, present; below it, the phase would form.](../assets/theory/solid_solution_tangent.svg)

Per mole of phase, the Gibbs energy of a binary of end-members A and B is the
straight line of the two end-members unmixed, ``(1-x)\,\mu_A^\circ + x\,\mu_B^\circ``,
plus ``G^{\text{mix}}``, negative between the two: the curve lies below the line,
and the same matter costs less mixed than apart (left panel). The tangent to the
curve at a composition ``x^\ast`` meets the edges ``x = 0`` and ``x = 1`` at the
chemical potentials of the two end-members in the phase, and it is through these
potentials that the phase exchanges matter with the rest of the system.

At equilibrium each element has one potential, shared by every phase, and those
potentials give each end-member a potential ``u_i`` in units of ``RT``, which
defines a plane over the compositions of the phase (right panel). A curve lying
above that plane everywhere is the Gibbs energy of an absent phase, none of its
compositions being as cheap as what the rest of the system already offers; a
curve touching it is that of a present phase, at the composition of the contact;
a curve dipping below it is that of a phase which would lower the Gibbs energy by
forming, so that the state is not yet an equilibrium. With ``g_i = \mu_i^\circ/RT``,
the phase is present exactly when ``\sum_i \exp(u_i - g_i - \ln\gamma_i) = 1``,
which [the certifying solver](@ref sec-theory-certificate) uses as an equation.
For ideal mixing ``\gamma_i = 1``, and ``\exp(u_i - g_i)`` is then the saturation
ratio ``\Omega_i`` end-member ``i`` would have as a pure phase: an ideal solid
solution forms as soon as ``\sum_i \Omega_i`` reaches one, while every
``\Omega_i`` is still below one and no end-member would precipitate on its own.
[Solid solutions in a calculation](@ref sec-tutorial-solid-solutions) reads these
ratios on a C-S-H gel.

## 2. Non-ideal mixing: the excess Gibbs energy

Real end-members interact. Everything beyond the ideal term is collected into
the **excess** Gibbs energy ``G^{\text{ex}}``, and the activity coefficient is
its partial molar derivative:

```math
G = \underbrace{\sum_k x_k\mu_k^\circ}_{\text{mechanical mixture}}
  + \underbrace{RT\sum_k x_k\ln x_k}_{\text{ideal mixing}}
  + \underbrace{G^{\text{ex}}(x,T)}_{\text{interaction}} ,
\qquad
\ln\gamma_k = \frac{\partial}{\partial n_k}\!\left(\frac{n\,G^{\text{ex}}}{RT}\right) ,
```

and then ``\ln a_k = \ln x_k + \ln\gamma_k``. This is exactly what the code
does: for these models `_ss_log_activities!` writes
`log(x[k] + ϵ) + _excess_ln_gamma(model, k, x, T)` into the log-activity vector,
so the two terms are visibly separate and a model supplies only the second. The
model of section 7 is written differently, and has its own method.

## 3. `RegularSolutionModel` — one interaction energy per pair

The simplest non-ideal form keeps only pairwise contacts, with one energy per
pair of end-members:

```math
\frac{G^{\text{ex}}}{RT} = \sum_{i<j}\frac{W_{ij}}{RT}\,x_i x_j ,
\qquad
\ln\gamma_k = \sum_{j\neq k}\frac{W_{kj}}{RT}x_j \;-\; \sum_{i<j}\frac{W_{ij}}{RT}x_i x_j ,
```

which is the symmetric multicomponent Margules expression, and it is
`_excess_ln_gamma(::RegularSolutionModel, …)` line for line. For a **binary** it
collapses to the form usually quoted:

```math
\ln\gamma_1 = \frac{W_{12}}{RT}x_2^2 ,
\qquad
\ln\gamma_2 = \frac{W_{12}}{RT}x_1^2 .
```

In the picture of nearest neighbors behind the model, ``W_{ij}`` (J/mol) is the
energy of an ``i``–``j`` contact less the mean of the ``i``–``i`` and ``j``–``j``
contacts, counted over the contacts of a mole. Fitted to data, it often has to
depend on temperature or pressure, and it then mixes enthalpy and entropy, but
its sign keeps the sign of the departure from ideal mixing
[Richet2001; Sec. 11.4b, pp. 262–263](@cite), and that sign is physical:

  - ``W_{ij} > 0`` — unlike neighbors are unfavorable. ``\gamma > 1``, the phase
    resists mixing, and above a threshold it unmixes;
  - ``W_{ij} < 0`` — unlike neighbors are favorable. ``\gamma < 1``, mixing is
    stabilized beyond ideal, and the phase may order.

`W` must be square and symmetric, which a **validating inner constructor**
enforces; the reason it is inner rather than outer is recorded in the source, and
it is a mistake worth not repeating.

### The threshold is ``W = 2RT``, and it is checkable

For a symmetric binary, the molar Gibbs energy of mixing is
``G^{\text{mix}}/RT = x\ln x + (1-x)\ln(1-x) + (W/RT)x(1-x)``. A phase is
unstable to unmixing where that function is concave, and

```math
\frac{\partial^2}{\partial x^2}\frac{G^{\text{mix}}}{RT}
 = \frac{1}{x} + \frac{1}{1-x} - \frac{2W}{RT} ,
```

which at ``x = 1/2`` is ``4 - 2W/RT``: negative as soon as ``W > 2RT``. So
``W/RT = 2`` — about **5 kJ/mol at 25 °C** — is the critical point of a symmetric
regular solution: below it one homogeneous phase, above it a miscibility gap
that widens as ``W`` grows. The threshold is evaluated, and the sign of the
second derivative tabulated on either side of it, in
[Solid solution models, in numbers](@ref sec-app-solid-solutions).

This matters in practice because **nothing in the activity expression detects
it**. Past the threshold it goes on returning values, and those values describe
a single phase that should have split. The check is therefore made at
construction: [`SolidSolutionPhase`](@ref) evaluates the second derivative and
refuses a model that unmixes at the temperature given, naming the interval,
unless the phase is declared with two instances (section 6) or the check is
switched off with `check_convexity = false`.

## 4. `RedlichKisterModel` — a binary that need not be symmetric

A regular solution is symmetric by construction: swapping the end-members leaves
``G^{\text{ex}}`` unchanged. Real binaries often are not, so the excess is
expanded in a polynomial of the composition difference. As implemented, for two
end-members,

```math
\ln\gamma_1 = \frac{x_2^2}{RT}\Big[a_0 + a_1(3x_1 - x_2) + a_2(x_1 - x_2)(5x_1 - x_2)\Big] ,
```
```math
\ln\gamma_2 = \frac{x_1^2}{RT}\Big[a_0 - a_1(3x_2 - x_1) + a_2(x_2 - x_1)(5x_2 - x_1)\Big] ,
```

with ``a_0, a_1, a_2`` in J/mol. Read the roles off the formulas: ``a_0`` is the
symmetric term and is exactly a regular solution's ``W_{12}``; ``a_1`` and
``a_2`` are the asymmetric corrections, so setting them to zero must reproduce
the regular model. That is an identity between two independently written
methods, and it is checked rather than trusted in
[Solid solution models, in numbers](@ref sec-app-solid-solutions), where the two
agree to the last digit.

`RedlichKisterModel` **requires exactly two end-members**, and
[`SolidSolutionPhase`](@ref) refuses the combination at construction rather than
letting a ternary reach an expression written for a binary. For three or more
end-members the choices are the ideal model, `RegularSolutionModel` with a full
``\mathbf{W}`` matrix, or the asymmetric models of section 9: the subregular
model, whose binaries are this model to the first order, the series of every
order, and the van Laar model.

### Published dimensionless parameters, and the temperature

Cemdata18 gives its non-ideal AFt and AFm binaries dimensionless Guggenheim
parameters ``\alpha_0`` and ``\alpha_1`` [Lothenbach2019](@cite) (its Table 1,
notes a, b, g and h, in `data/literature/Lothenbach2019.json`), and no
temperature dependence for them. The package turns them into
``a_k = \alpha_k R T_0`` with ``T_0 = 298.15`` K and keeps ``a_k`` in J/mol at every
temperature: the excess Gibbs energy is then all enthalpy, and the interaction
it puts into ``\ln\gamma`` at another temperature is ``\alpha_k T_0/T``, 7 % larger at
5 °C and 16 % smaller at 80 °C. Holding ``\alpha_k`` itself constant would make
the excess all entropy and leave ``\ln\gamma`` the same at every temperature.
The two limits have classical names, a *regular* solution for the first and an
*athermal* one for the second [Richet2001; Secs. 11.4b and 11.5e, pp. 262 and 270](@cite).
Neither choice is in the source, and no data in the package decides between
them; the edges of a gap measured at two temperatures would. What the choice
does to the gap of the AFm SO₄/OH binary, monosulfate and C₄AH₁₃, whose
printed edges are the 25 °C ones:

```@example guggenheim
using ChemistryLab, DynamicQuantities, Printf
p = literature_row("Lothenbach2019", "guggenheim_parameters", "AFm SO4/OH")
α0, α1 = ustrip(p.alpha0), ustrip(p.alpha1)
rk(T) = RedlichKisterModel(a0 = α0 * R_GAS * T, a1 = α1 * R_GAS * T)
@printf("printed gap: %.2f to %.2f\n", ustrip(p.gap_from), ustrip(p.gap_to))
println(" °C   a in J/mol held (the package)   α held")
for T_C in (5, 25, 50, 80)
    T = T_C + 273.15
    held_a = common_tangent(rk(298.15); T)
    held_α = common_tangent(rk(T); T)
    @printf("%3d       %.3f to %.3f            %.3f to %.3f\n", T_C, held_a..., held_α...)
end
```

Both give the printed gap at 25 °C. The package's choice narrows it as the
temperature rises, from 0.48–0.98 at 5 °C to 0.57–0.94 at 80 °C (the mole
fraction of C₄AH₁₃, the OH end, as the edges are printed); the other keeps it
where Cemdata18 prints it.

## 5. What the three models do to an activity

Evaluated over the whole composition range in
[Solid solution models, in numbers](@ref sec-app-solid-solutions), with
``W = \pm 4`` kJ/mol — inside the stable range of §3, so that the numbers
describe a phase that stays homogeneous. Two things to read off it.

A positive ``W`` pushes the activity of a dilute end-member **above** its mole
fraction: the host is rejecting it. A negative ``W`` pulls it below: the host is
stabilizing it. And at ``x_k \to 1`` all three models converge, because
``\gamma_k \to 1`` as the phase becomes pure — the standard state of an
end-member is the pure end-member, so they agree there by construction.

The ideal case is the one to keep in mind for cement. An end-member at
``x_k = 10^{-3}`` has ``a_k = 10^{-3}``, so its saturation index sits three
decades below the pure phase: **a trace component is stabilized simply by being
diluted in a host.** That is how a solid solution takes up an ion which would
never precipitate as a phase of its own, and it is the whole reason C-S-H and
the AFm phases can absorb alkalis, sulfate and carbonate continuously instead of
in jumps.

## 6. The spinodal, the common tangent, and why one amount is not enough

Sections 3 and 4 gave models whose excess term can be strong enough to make the
molar Gibbs energy of mixing **concave** over an interval. This section is about
what that means for the equilibrium, because it is not a detail of the model: it
changes what the answer *is*.

### Convexity is the hypothesis under everything else

For a binary write the molar Gibbs energy of mixing, in units of ``RT``, as

```math
\frac{g(x)}{RT} = x\ln x + (1-x)\ln(1-x) + \frac{g^{\mathrm{ex}}(x)}{RT} .
```

The ideal part is convex everywhere — its second derivative is
``1/x + 1/(1-x) > 0`` — so unmixing is always the excess term's doing. Where

```math
\frac{\mathrm{d}^2 g}{\mathrm{d}x^2} < 0
```

the phase is **inside its spinodal**, and there the straight line joining two
compositions lies *below* the curve between them. A system at an intermediate
overall composition therefore lowers its energy by separating into those two, and
the minimum of ``G`` is not a point but a **pair**.

For the symmetric regular solution of section 3 the criterion is exactly
``\mathrm{d}^2g/\mathrm{d}x^2 = 4 - 2W/RT`` at ``x = \tfrac12``, so unmixing
begins at ``W = 2RT`` — the threshold section 3 already quoted, here derived from
the same inequality. [`spinodal_interval`](@ref) evaluates the second derivative
on a grid rather than in closed form, which is why the same routine covers the
asymmetric Redlich-Kister case without a separate derivation.

### The answer inside a gap is the common tangent

Which pair? Write ``g`` now for the molar Gibbs energy of the phase with its
standard terms included,

```math
g(x) = (1-x)\,\mu_A^\circ + x\,\mu_B^\circ + \Delta_{\text{mix}}g(x) ,
```

``x`` the mole fraction of the end-member B and ``\Delta_{\text{mix}}g`` the
energy of mixing above; the linear terms change neither the second derivative
nor, as shown below, the pair. Let ``n`` moles of the phase, of overall composition
``\bar{x}``, be shared between two parts ``\alpha`` and ``\beta`` holding
``n_\alpha`` and ``n_\beta`` moles at compositions ``x_\alpha`` and ``x_\beta``.
At fixed ``T`` and ``P`` the equilibrium minimizes

```math
G = n_\alpha\, g(x_\alpha) + n_\beta\, g(x_\beta)
\quad\text{subject to}\quad
n_\alpha + n_\beta = n , \qquad n_\alpha x_\alpha + n_\beta x_\beta = n\,\bar{x} ,
```

the conservation of the two end-members. With a multiplier ``\lambda_0`` for
the first constraint and ``\lambda_1`` for the second, the Lagrangian is
``\mathcal{L} = G - \lambda_0\,(n_\alpha + n_\beta - n) - \lambda_1\,(n_\alpha x_\alpha + n_\beta x_\beta - n\bar{x})``,
and its stationarity with respect to the compositions and to the amounts, both
parts being present (``n_\alpha, n_\beta > 0``), reads

```math
\frac{\partial\mathcal{L}}{\partial x_\varphi} = n_\varphi\,\bigl[g'(x_\varphi) - \lambda_1\bigr] = 0 ,
\qquad
\frac{\partial\mathcal{L}}{\partial n_\varphi} = g(x_\varphi) - \lambda_0 - \lambda_1 x_\varphi = 0 ,
\qquad \varphi = \alpha, \beta .
```

The four equations say one thing: the straight line
``\ell(x) = \lambda_0 + \lambda_1 x``
passes through the curve at ``x_\alpha`` and at ``x_\beta`` and has the slope of
the curve at both, so that it is tangent to ``g`` at the two compositions, the
**common tangent**. Eliminating ``\lambda_0`` gives the two equations in two
unknowns that [`common_tangent`](@ref) solves:

```math
g'(x_\alpha) = g'(x_\beta) = \frac{g(x_\beta) - g(x_\alpha)}{x_\beta - x_\alpha} .
```

The multipliers are chemical potentials. The tangent to ``g`` at a composition
``x`` meets the vertical ``x = 0`` at ``g - x\,g'``, the chemical potential of
A at that composition, and the vertical ``x = 1`` at ``g + (1-x)\,g'``, that of
B, partial molar quantities being read as the intercepts of the tangent
[Richet2001; Sec. 3.2c, Fig. 3.2, pp. 51–52](@cite). The stationarity
therefore states that ``\mu_A = \lambda_0`` and ``\mu_B = \lambda_0 + \lambda_1``
in both parts: each end-member has the same chemical potential in the two
compositions, which is the condition of equilibrium between any two phases
[Richet2001; Sec. 7.3, Fig. 7.5, pp. 167–169](@cite).

Nothing requires the two points to be at the same height:
``g(x_\beta) - g(x_\alpha) = \lambda_1 (x_\beta - x_\alpha)``, which vanishes
only for a horizontal tangent, ``\lambda_1 = 0``, as for a symmetric mixing energy
between end-members of equal standard potentials. Nor are they the minima of
``g``, where ``g' = 0`` while here ``g' = \lambda_1``. A term linear in ``x``, as the
standard potentials add, tilts the curve and its tangent alike and moves
neither point: the pair depends on the mixing model only.

What the split gains follows from the lever rule,
``n_\alpha = n\,(x_\beta - \bar{x})/(x_\beta - x_\alpha)``: the Gibbs energy of
the two parts is ``G = n\,\ell(\bar{x})``, the height of the tangent at
``\bar{x}``, below the ``n\,g(\bar{x})`` of a single composition wherever the
curve rises above its tangent. Between ``x_\alpha`` and ``x_\beta`` the
equilibrium is this mixture, and the energy follows the tangent rather than the
curve. Stationarity alone does not select it, however: a single composition,
``x_\alpha = x_\beta = \bar{x}``, satisfies the same equations with the tangent
drawn at ``\bar{x}``. What separates the two is a condition on the whole curve,
that it lie nowhere below the line, which is the tangent-plane test of
[Michelsen1982](@citet), the subject of a subsection below.

The common tangent construction is wider than the spinodal: the *binodal*
``[x_\alpha, x_\beta]`` contains the spinodal, and between the two the phase is
metastable rather than unstable. Inside the spinodal any small fluctuation of
composition lowers ``G`` and grows by itself; between the spinodal and the
binodal every small fluctuation raises it, and the phase splits only through a
nucleus large enough to pay for its interface
[Richet2001; Sec. 7.3c, p. 170](@cite). A minimization sees only the tangent.

![Left: the molar Gibbs energy of an asymmetric binary, tilted by the standard potentials of its end-members, with two wells, the spinodal shaded between its inflection points, and the common tangent touching the curve at two compositions and two heights; extended to x = 0 and x = 1, the tangent marks the chemical potentials of A and B. Right: a composition inside the gap, which costs more as one phase than the same matter split into the two compositions of the tangent; the phase is then held twice, one instance at each.](../assets/theory/miscibility_gap_instances.svg)

The right panel states the difficulty the rest of this section resolves: inside
the gap the equilibrium holds the phase at two compositions at once, in the
proportions of the lever rule, which a formulation carrying one amount per
end-member has no means to express.

### Computing the pair from the model alone

The pair is not found by minimizing. It is **computed from the model alone**, by
[`common_tangent`](@ref): two equations in two unknowns, solved by Newton with
`ForwardDiff` supplying the derivatives. Nothing about the rest of the chemical
system enters — not the aqueous solution, not the other phases, not the element
budget — which is why it costs microseconds and can be used as the starting point
of a full equilibrium rather than as its result.

The construction is that of [GlynnReardon1990](@citet) for binary solid
solutions, which PHREEQC [ParkhurstAppelo2013](@cite) uses as well. It was validated here against a case
with a closed form: for a symmetric model the pair must be symmetric about
``x = 1/2``, and it comes out at ``(0.070720,\ 0.929280)`` with a residual of
``2.3\times10^{-13}`` and the symmetry exact to the last bit. For a symmetric
regular model the pair has a closed form,
``W/RT = \ln[(1-x)/x]/(1-2x)`` [Richet2001; Sec. 10.4b, Eq. (10.29), p. 230](@cite),
which these two compositions satisfy with ``W = 3RT`` to the six figures given.

[`miscibility_split`](@ref) then applies the lever rule. Inside the gap the two
**compositions are fixed** and only their proportions move with the overall
composition — which is what makes a miscibility gap flat in a phase diagram, and
what a single-composition answer cannot reproduce at any resolution.

### Why a formulation with one amount per species cannot hold it

The composition vector of a `ChemicalSystem` has one entry per species, so a
species belongs to exactly one phase and a phase has exactly one composition. A
common-tangent pair is *two compositions of one substance*, and there is no way
to write that down with one amount each.

So the formulation has to be given the substance twice. That is what
`SolidSolutionPhase(...; instances = 2)` does: `ChemicalSystem` builds a second
copy of each end-member under a derived symbol (`monosulphate12#2`) sharing the
same thermodynamic record, and each copy belongs to its own phase. The groups
stay disjoint, the duplicated columns of the conservation matrix are copies of
columns already there — so the row rank, and with it conservation, is untouched —
and the minimization is free to populate either lobe or both.

It is refused for a convex model, and the reason is worth stating because it is
not conservatism. Two instances of a convex phase are **degenerate**: the energy
is the same however the amount is split between them, so the minimum becomes a
flat manifold and the optimizer is asked to pick a point on it arbitrarily.
Inside a spinodal the common-tangent pair is unique and no such direction exists.

### Michelsen's tangent-plane distance, which is what detects it

Representation is one problem; knowing that it is needed is another. A phase
sitting inside its own spinodal satisfies every first-order condition — its
end-members are stationary, the element balance closes, nothing absent is
supersaturated — so stationarity alone certifies a non-minimum.

The test that separates them is Michelsen's [Michelsen1982](@cite): for a trial composition
``\hat{\mathbf{x}}`` of the phase, the **tangent-plane distance**

```math
D(\hat{\mathbf{x}}) = \sum_k \hat{x}_k\,
   \bigl[\,\mu_k(\hat{\mathbf{x}}) - \mu_k(\mathbf{x}^\star)\,\bigr]
```

measures how far the Gibbs surface at ``\hat{\mathbf{x}}`` lies above the tangent
plane drawn at the reported composition ``\mathbf{x}^\star``. If ``D`` is
negative anywhere, some other composition lies *below* that plane and the
reported state is not a minimum. At ``\hat{\mathbf{x}} = \mathbf{x}^\star`` the
distance is exactly zero and is a stationary point of the search, which is why
the measure has to be started from the **end-member corners** — the other lobe is
where the negative value lives.

`OptimaSolver`'s `phase_split_measure` computes this for every present
mole-fraction phase and folds the result into the certificate. Concretely:
`check_convexity = false` stops being a silent loss of the proof — a concave
declaration that does unmix now fails to certify, with the phase named.

### A second composition, declared or added when needed

The certificate's measure of a phase is the sum over its members of their
activities from the **dual** solution divided by their activity coefficients
from the **primal** one, an estimate of the mole fraction the phase would take:

```math
\Omega_k = \sum_{j\in l_k}\pi_j , \qquad
\pi_j = \frac{\omega_j(\hat{\mathbf{u}})}{\gamma_j(\hat{\mathbf{n}})} .
```

It is the phase stability index ``\Lambda_k = \log_{10}\Omega_k`` of
[Kulik2013](@citet), who derive it from the same KKT conditions and call
``\Omega_k`` *"a generalization of the saturation index"*. The package
evaluates it in two places: as ``\Omega = \sum_i 10^{SI_i}`` (the ideal case,
``\gamma = 1``) when it repairs a start, and through the log-sum-exp of
``d_j = u_i - g_i - \ln\gamma_i`` in `phase_split_measure`.

Holding a gap takes the phase twice. Cemdata18 ships the AFm and AFt binaries
under two names each, so that the database itself declares each binary twice;
`instances = 2` declares the second composition with a keyword.
`instances = :auto` asks for less: the phase is solved with one composition,
and a second instance is added only when the certificate finds it wanting to
split. On a binary whose overall composition the element budget holds inside
the gap, calcite and magnesite in equal amounts, that finds the common-tangent
pair:

```@example split
using ChemistryLab, DynamicQuantities, Printf
sp = Dict(symbol(s) => s for s in build_species(datapath("slop98-inorganic-thermofun.json"); verbose = false))
gap = RedlichKisterModel(a0 = 14_000.0)
carbonate = SolidSolutionPhase("carbonate", [sp["Cal"], sp["Mgs"]]; model = gap, instances = :auto)
cs = ChemicalSystem([sp[s] for s in split("H2O@ H+ OH- CO2@ HCO3- CO3-2 Ca+2 Mg+2 Cal Mgs")],
                    ["H2O@", "H+", "Ca+2", "Mg+2", "CO3-2", "Zz"]; solid_solutions = [carbonate])
st = ChemicalState(cs)
set_quantity!(st, "H2O@", 1.0u"kg"); set_quantity!(st, "Cal", 0.025u"mol"); set_quantity!(st, "Mgs", 0.025u"mol")
eq, cert = equilibrate_certified(st; b = Float64.(cs.SM.A) * ustrip.(us"mol", st.n))
n = ustrip.(us"mol", eq.n)
@printf("certified %s, instances added to %s\n", cert.optimal, cert.instances_added)
for (g, ph) in zip(eq.system.ss_groups, eq.system.solid_solutions)
    @printf("  %-12s %.5f mol, x(calcite) = %.4f\n", name(ph), sum(n[g]), n[g[1]] / sum(n[g]))
end
@printf("common tangent of the model: %.4f and %.4f\n", common_tangent(gap)...)
```

It does not make a nonconvex problem easy: the certificate names a phase to split
only once the one-composition solve has converged.

The executed counterpart of this section is
[the miscibility-gap page](@ref ex-miscibility-gap), which runs one cement three
ways and reports the certificate each time.

## 7. Mixing on several sites: `SublatticeModel`

Sections 1 to 4 mix **end-members**. In a C-S-H the entities that actually mix
are smaller: a few structural positions of the silicate chain and of the
interlayer, each held by one species or another. An end-member is then one
particular filling of those positions, and mixing is random on each position
separately. That is the model of [Kulik2011](@citet) for CSH3T and of
[Myers2014](@citet) for the CNASH gel, and it is what [`SublatticeModel`](@ref) implements.

**The model.** A formula unit has sites ``s``, site ``s`` counted ``m_s`` times,
and end-member ``k`` puts the species ``\sigma_s(k)`` on site ``s``. The fraction
of site ``s`` held by species ``i`` is the total mole fraction of the end-members
that put it there, ``y_{s,i} = \sum_k [\sigma_s(k) = i]\,x_k``, and the
configurational Gibbs energy is the ideal one of each site:

```math
\frac{G^{\text{mix}}}{RT} = \sum_s m_s \sum_i y_{s,i}\ln y_{s,i}
\qquad\Longrightarrow\qquad
\boxed{\;\ln a_k = \sum_s m_s \ln y_{s,\sigma_s(k)}\;}
```

The implication is ``\partial(n\,G^{\text{mix}}/RT)/\partial n_k``, exact because
every end-member fills every site once. One site counted once is section 1 again.
[Kulik2011](@citet) and [Myers2014](@citet) write the same model as a
*fictive activity coefficient*
``\lambda_k = a_k/x_k``, which is what [`excess_ln_gamma_expression`](@ref)
returns for this model.

**CSH3T.** Two bridging-tetrahedral sites hold either Si or Ca. TobH puts Si on
both, T2C puts Ca on both, and the ordered T5C puts Ca on the first and Si on the
second [Kulik2011; Eq. 19](@cite). The Cemdata18 records are half of its
formula units, so on them each site counts one half:

```@example sublattice
using ChemistryLab
subs = Dict(symbol(s) => s for s in build_species(datapath("cemdata18-thermofun.json"); verbose = false))
csh3t = sublattice_model("Kulik2011:csh3t", [subs[n] for n in ("CSH3T-TobH", "CSH3T-T5C", "CSH3T-T2C")])
csh3t.multiplicity
```

`sublattice_model` finds that factor itself, by comparing each record's formula
with the one [Kulik2011](@citet) prints, and refuses a record that is not such
a multiple.
At ``\mathbf{x} = (0.2, 0.5, 0.3)`` the first site holds Si only through TobH, the second
through TobH and T5C:

```@example sublattice
site_fractions(csh3t, [0.2, 0.5, 0.3])
```

and the fictive activity coefficient of T5C is Eq. (20) of [Kulik2011](@citet),
halved:

```math
\ln\lambda_{\mathrm{T5C}} = \tfrac12\ln(x_{\mathrm{TobH}} + x_{\mathrm{T5C}})
  + \tfrac12\ln(x_{\mathrm{T2C}} + x_{\mathrm{T5C}}) - \ln x_{\mathrm{T5C}}
```

The block below builds that expression from the model and subtracts the formula
above; the difference simplifies to zero:

```@example sublattice
using Symbolics
x = [Symbolics.variable(:x, i) for i in 1:3]   # x₁ = TobH, x₂ = T5C, x₃ = T2C
λ = excess_ln_gamma_expression(csh3t, 2, 3)
simplify(expand(λ - (0.5log(x[1] + x[2]) + 0.5log(x[3] + x[2]) - log(x[2]))))
```

**Two properties the solver relies on.**

  - *Convexity.* ``G^{\text{mix}}`` is a sum of convex functions of the site
    fractions, which are linear in ``\mathbf{x}``, so it is convex, and strictly so when
    the occupancy matrix (one row per site and species, one column per
    end-member) has full column rank, stored as `rank`. The activities are the
    gradient of that energy, so a certificate keeps its `:global_minimum` scope.
  - *The dilute limit.* As ``x_k \to 0`` only the site fractions of the species
    ``k`` alone carries vanish with it: ``\ln a_k \simeq e_k \ln x_k``, with
    ``e_k`` the multiplicity of those sites, stored as `exponents`. T5C owns no
    species of its own, ``e_k = 0``: its activity stays finite as it disappears,
    and it can be absent from a CSH3T that is present.

```@example sublattice
csh3t.exponents
```

**What the solver does with it.** The composition of a mixing phase is found by
successive substitution, and on a sublattice model that iteration diverges: its
rate is set by the sum of the multiplicities, nine for the CNASH gel. Such a phase
is therefore inverted by Newton's method, and the members with ``e_k = 0`` are
declared to the certificate as able to leave the phase, which then tests them by
the inequality a pure phase obeys (OptimaSolver's `SolutionPhase`, keywords
`newton` and `bounded_members`). A system with a sublattice phase is first solved
with that phase mixing its end-members ideally, and the answer is the start of
the sublattice solve.

`data/solid_solutions.toml` ships `CNASH_ss` as the model of [Myers2014](@citet) and `CSH3T` as
Cemdata18 ships it, ideal between end-members; `sublattice_model("Kulik2011:csh3t",
members)` gives the site form.

## 8. Sites with the energies of the compounds: `CompoundEnergyModel`

Section 7 mixes ideally on the sites and keeps the end-members' own standard
energies, ``\sum_k x_k G^\circ_k``. That is exact only when the end-members are
independent in the site fractions. The CASH+ model of C-S-H [Kulik2022](@cite)
has two mixing sites: the bridging tetrahedron holds a silicate, a vacancy or a
calcium (S, v, C), and the interlayer a vacancy or a calcium (v, C). Its six
end-members are the six ways of filling them, TSvh to TCCh, but three site
fractions fix the composition, and the six members are not independent. The
reciprocal reaction

```math
\mathrm{TSvh} + \mathrm{TCCh} = \mathrm{TSCh} + \mathrm{TCvh}
```

leaves every site fraction unchanged, yet its standard Gibbs energy is not zero.
The compound energy formalism, in which [Kulik2022](@citet) write the model, gives such
a reaction its energy. [`CompoundEnergyModel`](@ref) implements it.

**The Gibbs energy.** The end-members must be every *compound* of the sites, each
once. Compound ``j`` puts species ``j_s`` on site ``s`` (the ``\sigma_s(j)`` of
section 7). Per formula unit,

```math
G = \underbrace{\sum_j G^\circ_j \prod_s y_{s,j_s}}_{G_\text{ref}}
  + RT \sum_s m_s \sum_i y_{s,i}\ln y_{s,i}
  + \sum_s \sum_{i<l} W_{s,il}\, y_{s,i}\, y_{s,l} .
```

The first term, the *reference surface*, equals ``G^\circ_j`` at compound ``j``
and is linear in the fractions of each site, so a reaction between compounds that
changes no site fraction keeps its energy. The last term is a regular interaction
between two species of one site, symmetric, with the parameters ``W`` of
[Kulik2022](@citet). ``G`` depends on the site fractions alone.

**The activities.** With ``y_{s,i} = n_{s,i}/n``, the chemical potential of
member ``k`` is ``\mu_k = \partial(nG)/\partial n_k``. For a function of the site
fractions, ``\partial y_{s,i}/\partial n_k = ([k_s = i] - y_{s,i})/n``, so

```math
\mu_k = G + \sum_s\Big(\frac{\partial G}{\partial y_{s,k_s}} - \sum_i y_{s,i}\frac{\partial G}{\partial y_{s,i}}\Big).
```

On the reference surface ``\sum_i y_{s,i}\,\partial G_\text{ref}/\partial y_{s,i} = G_\text{ref}``
for every site, because ``G_\text{ref}`` is linear in each site's fractions. With
``n_\text{site}`` mixing sites, and writing ``\ln a_k = (\mu_k - G^\circ_k)/RT``,

```math
\ln a_k = \sum_s m_s \ln y_{s,k_s}
  + \frac{1}{RT}\Big(\sum_s \frac{\partial G_\text{ref}}{\partial y_{s,k_s}} - (n_\text{site}-1)\,G_\text{ref} - G^\circ_k\Big)
  + \frac{1}{RT}\sum_s\Big(\sum_l W_{s,k_s l}\, y_{s,l} - \sum_{i<l} W_{s,il}\, y_{s,i}\, y_{s,l}\Big).
```

The middle term vanishes at every compound. It also does not change when each
``G^\circ_j`` is shifted by the energies of its elements, since those are a sum
over the sites of the species' shares, which the reference surface interpolates
exactly. So the convention of the database does not matter, but the energies of
the members do: the solver hands them to the model at the temperature of the
solve.

**Three consequences.**

  - *The amounts are not unique.* Every split of the same site fractions between
    the members has the same ``G`` and the same element content: in CASH+, six
    amounts for three independent site fractions and one total. Left to choose, the solver returns one
    of them, which one depending on where the search started. It is therefore given
    one split, the product of the site fractions, ``x_j = \prod_s y_{s,j_s}``, by adding to the energy it
    minimizes ``RT\,D(\mathbf{x})``, with

    ```math
    D(\mathbf{x}) = \sum_j x_j \ln x_j - \sum_s \sum_i y_{s,i}\ln y_{s,i} .
    ```

    ``D`` is the divergence of ``\mathbf{x}`` from the product of its own site fractions:
    it is never negative, and it vanishes only at that product, where its gradient
    ``\ln x_k - \sum_s \ln y_{s,k_s}`` vanishes too. The minimum over the splits
    of one set of site fractions is therefore the product, with the energy and
    the chemical potentials of the model. The activity of a vanishing member now
    goes as its mole fraction, as in a mixture of end-members. Measured on twelve
    pastes of CASH+NK, the equilibria certify with or without ``D``, at the same
    cost; what ``D`` adds is an answer that does not depend on the start.
  - *Convexity must be decided.* The reference surface is multilinear, and a
    reciprocal energy can make ``G`` concave somewhere. In the site fractions, each
    site contributes a curvature of at least ``2m_s + \mu_s``, the configurational
    bound of section 9 plus the smallest tangent eigenvalue ``\mu_s`` of its
    interactions ``\mathbf{W}/RT``, and with two sites the reference surface couples them
    through a constant matrix ``\mathbf{C}``, the ``G^\circ_j/RT`` of the compounds projected
    on the two tangent spaces. So ``G`` is convex when
    ``(2m_1 + \mu_1)(2m_2 + \mu_2) > \lVert \mathbf{C} \rVert^2``. That proves the CASH+
    core convex (34.3 against 13.2 at 25 °C), and a certificate on it keeps its
    `:global_minimum` scope as far as the gel goes; the bound does not decide
    CASH+NK, whose sampled Hessian shows no concave point, and its scope is then
    `:local_minimum` or `:kkt_point`, according to the curvature at the answer
    ([`mixing_convexity`](@ref)).
  - *Euler's relation holds.* ``\sum_k x_k\,\mu_k = G`` for any amounts with the
    site fractions ``\mathbf{y}``, as it must for a Gibbs energy.

The block below checks the gradient and the last point on the core model, at a
composition of the six members. The activities are compared with the gradient of
``n(G/RT + D)`` coded from the formulas above and differentiated by ForwardDiff;
at the product split, ``D`` and its gradient vanish:

```@example cef
using ChemistryLab, DynamicQuantities, ForwardDiff, LinearAlgebra
subs = Dict(symbol(s) => s for s in build_species(datapath("cemdata18-cashplus.json"); verbose = false))
members = ["TSvh", "TSCh", "Tvvh", "TCvh", "TvCh", "TCCh"]
cash = compound_energy_model("Kulik2022:cashplus", [subs[m] for m in members])
T = 298.15
g = [ustrip(us"J/mol", subs[m][:ΔₐG⁰](T = T * u"K", P = 1.0e5u"Pa"; unit = true)) for m in members] ./ (ChemistryLab.R_GAS * T)
o = cash.lattice.occupancy
function nG(n; split = true)
    x = n ./ sum(n)
    y = site_fractions(cash, x)
    G = sum(g[j] * y[1][o[1, j]] * y[2][o[2, j]] for j in 1:6) + sum(v * log(v) for ys in y for v in ys)
    G += sum(W / (ChemistryLab.R_GAS * T) * y[s][i] * y[s][l] for (s, i, l, W) in cash.interactions)
    split && (G += sum(v * log(v) for v in x) - sum(v * log(v) for ys in y for v in ys))   # D(x)
    return sum(n) * G
end
n = [0.10, 0.25, 0.05, 0.20, 0.15, 0.25]
lna = ChemistryLab._ss_log_activities!(zeros(6), 1:6, n ./ sum(n), cash, T, 0.0, g)
μ = ForwardDiff.gradient(nG, n)
y = site_fractions(cash, n ./ sum(n))
xp = [y[1][o[1, j]] * y[2][o[2, j]] for j in 1:6]          # the product split
(gradient = maximum(abs, lna .+ g .- μ), euler = dot(n, μ) - nG(n),
 D_at_product = nG(xp) - nG(xp; split = false))
```

All three are at the rounding of the arithmetic. The model and its data are
in `data/literature/Kulik2022.json`, the sodium and potassium of
[Miron2022a, Miron2022b](@citet) in `Miron2022a.json` and `Miron2022b.json`, and the
twelve end-members in the database `cemdata18-cashplus.json`.
[The CASH+ page](@ref ex-cashplus-csh) computes the C-S-H in water and in alkali
solutions with it.

## 9. More than two end-members

### Asymmetric mixing of any number of end-members: `SubregularSolutionModel`

The regular model of section 3 is symmetric on every edge, and the
Redlich–Kister model of section 4 is asymmetric but binary.
[HelffrichWood1989](@citet) write the asymmetric Margules model for any number of
end-members (their Eq. 5′):

```math
G^{\text{ex}} = \sum_{i<j} x_i x_j \Bigl\{ W_{ij}\Bigl[x_j + \tfrac12 \sum_{k\neq i,j} x_k\Bigr]
              + W_{ji}\Bigl[x_i + \tfrac12 \sum_{k\neq i,j} x_k\Bigr] \Bigr\}
              + \sum_{i<j<k} W_{ijk}\, x_i x_j x_k .
```

On the edge ``i``–``j`` the bracket sums vanish and the excess is
``x_i x_j (W_{ij} x_j + W_{ji} x_i)``: ``W_{ij}`` is ``RT\ln\gamma_i`` at infinite
dilution of ``i`` in ``j``, and ``W_{ji}`` the converse. The ternary coefficients
``W_{ijk}`` are a separate measurement. [HelffrichWood1989](@citet) show that the binaries
do not determine them, even when every binary is symmetric, and that no
coefficient of higher order appears.

[HelffrichWood1989](@citet) give the activity coefficients (their Eq. 6′)
without the intermediate
steps. They are short once the excess is written on the simplex, where
``\sum_{k\neq i,j} x_k = 1 - x_i - x_j``. The bracket of a pair is then

```math
W_{ij}\Bigl[x_j + \tfrac12(1 - x_i - x_j)\Bigr] + W_{ji}\Bigl[x_i + \tfrac12(1 - x_i - x_j)\Bigr]
= a_{ij} + b_{ij}(x_i - x_j),
\qquad
a_{ij} = \frac{W_{ij} + W_{ji}}{2},\quad b_{ij} = \frac{W_{ji} - W_{ij}}{2},
```

so each pair is a Redlich–Kister term of the first order, ``a_0 = a_{ij}`` and
``a_1 = b_{ij}``, and two end-members are exactly the model of section 4 with
those two coefficients. Write ``g = G^{\text{ex}}/RT`` in this form, as a
function of ``x_1, \dots, x_n``. The excess of ``n_k`` moles is ``N g(\mathbf{n}/N)``,
``N = \sum_k n_k``, and since ``\partial x_l/\partial n_k = (\delta_{kl} - x_l)/N``,

```math
\ln\gamma_k = \frac{\partial (N g)}{\partial n_k}
            = g + \frac{\partial g}{\partial x_k} - \sum_l x_l \frac{\partial g}{\partial x_l} .
```

The formula holds whatever ``g`` is off the simplex, since ``N g(\mathbf{n}/N)``
reads ``g`` only on it. For a pair, ``t = x_i x_j[a + b(x_i - x_j)]`` has

```math
\frac{\partial t}{\partial x_i} = x_j\bigl[a + b(2x_i - x_j)\bigr],
\qquad
\frac{\partial t}{\partial x_j} = x_i\bigl[a + b(x_i - 2x_j)\bigr],
\qquad
\sum_l x_l \frac{\partial t}{\partial x_l} = x_i x_j\bigl[2a + 3b(x_i - x_j)\bigr],
```

and a triple, ``W_{ijk} x_i x_j x_k``, has ``\sum_l x_l\,\partial/\partial x_l`` equal to
three times itself. Collecting, with every coefficient divided by ``RT``,

```math
\ln\gamma_k = \sum_{j\neq k} x_j\bigl[a_{kj} + b_{kj}(2x_k - x_j)\bigr]
            - \sum_{i<j} x_i x_j\bigl[a_{ij} + 2b_{ij}(x_i - x_j)\bigr]
            + \sum_{\substack{i<j\\ i,j\neq k}} W_{kij}\, x_i x_j
            - 2\sum_{i<j<l} W_{ijl}\, x_i x_j x_l ,
```

with ``b_{kj} = (W_{jk} - W_{kj})/2``, which is
`_excess_ln_gamma(::SubregularSolutionModel, …)` term by term. Setting
``W_{ij} = W_{ji}`` and every ``W_{ijk}`` to zero leaves the regular model of
section 3, ``b = 0``. The test suite checks this expression against
[HelffrichWood1989; Eq. 6′](@cite), for each end-member of a quaternary with its twelve binary and four
ternary coefficients, and Gibbs–Duhem and ``\sum_k x_k \ln\gamma_k = g`` as
identities.

### A series for each pair: `MulticomponentRedlichKisterModel`

[RedlichKister1948](@citet) write the excess of a binary as a series in the
difference of the two mole fractions, which changes sign when the components
are exchanged, and that of a ternary as the sum of its three binaries, each
evaluated at the mole fractions of the ternary as they are, plus a term in the
product of the three (their Eqs. 17 and 18):

```math
\frac{G^{\text{ex}}}{RT} = \sum_{(i,j)} x_i x_j \sum_{k\ge 0} L^{ij}_k (x_i - x_j)^k
 + \sum_{(i,j,l)} x_i x_j x_l \bigl[C + D_1 (x_j - x_l) + D_2 (x_l - x_i)\bigr] .
```

The order of a pair matters for its odd terms, which change sign with it, and
the order of a triple for ``D_1`` and ``D_2``. Beyond three components
[RedlichKister1948](@citet) add the further pairs and triples and a term in four mole fractions, which it
does not write out and which is not implemented. Its coefficients are in units
of ``2.303\,RT``, the paper working with decimal logarithms.

The activity coefficients follow from Eq. (14) of [RedlichKister1948](@citet),
the partial molar
derivative written in mole fractions, which is the formula of the subregular
model above. For a pair, with ``d = x_i - x_j``, ``S = \sum_k L_k d^k`` and
``S' = \sum_k k L_k d^{k-1}``, the term ``t = x_i x_j S`` has

```math
\frac{\partial t}{\partial x_i} = x_j S + x_i x_j S',
\qquad
\frac{\partial t}{\partial x_j} = x_i S - x_i x_j S',
\qquad
\sum_m x_m \frac{\partial t}{\partial x_m} = 2 x_i x_j S + x_i x_j\, d\, S' ,
```

and for a triple, with ``U = C + D_1(x_j - x_l) + D_2(x_l - x_i)``, the term
``t = x_i x_j x_l U`` has

```math
\frac{\partial t}{\partial x_i} = x_j x_l U - x_i x_j x_l D_2,
\quad
\frac{\partial t}{\partial x_j} = x_i x_l U + x_i x_j x_l D_1,
\quad
\frac{\partial t}{\partial x_l} = x_i x_j U + x_i x_j x_l (D_2 - D_1),
\quad
\sum_m x_m \frac{\partial t}{\partial x_m} = 3t + x_i x_j x_l (U - C) .
```

``\ln\gamma_r`` is the sum over every pair and triple of ``t - \sum_m x_m
\partial t/\partial x_m``, plus the derivatives with respect to ``x_r`` of the
terms that contain it, which is `_excess_ln_gamma(::MulticomponentRedlichKisterModel, …)`.
Two end-members give the model of section 4 with ``a_k = L_k``. To the first
order the series is the subregular model above, with
``L_0 = (W_{ij} + W_{ji})/2`` and ``L_1 = (W_{ji} - W_{ij})/2``: [RedlichKister1948](@citet)
call the first-order binary "the equation of Margules", and the
subregular section shows why the ternary extensions coincide.

[RedlichKister1948](@citet) also give the ratio of two activity coefficients in
a ternary (their Eq. 19). Its term in ``C_{12}``, ``C_{12}[3(x_1 - x_2)^2 - 1]/2``, is the value
on the binary: in the ternary, Eq. (14) gives ``C_{12}[(x_1 - x_2)^2 - 2x_1x_2]``,
which is that plus ``C_{12}\,x_3(2 - x_3)/2``. Every other term of Eq. (19)
agrees with Eq. (14). The implementation follows Eq. (14), and the test suite
checks it against the paper's worked ternary, heptane, methanol and toluene with
the coefficients of its Eq. (21) and the association of methanol left out,
through its Eqs. (22) and (23).

### Unequal sizes: `VanLaarModel`

[HollandPowell2003](@citet) make the regular model asymmetric by giving each
end-member a size ``\alpha_i`` and weighting the mole fractions by it,
``\varphi_i = \alpha_i x_i / \sum_l \alpha_l x_l``:

```math
G^{\text{ex}} = \sum_{i<j} \varphi_i \varphi_j B_{ij},
\qquad
B_{ij} = \frac{2 \sum_l \alpha_l x_l}{\alpha_i + \alpha_j}\, W_{ij} .
```

Only the ratios of the sizes matter. Equal sizes give back the regular model of
section 3, which [HollandPowell2003](@citet) call the symmetric
formalism, and two end-members
give the van Laar binary, asymmetric as soon as the sizes differ. The papers
that fit ``W_{ij}`` and ``\alpha_i`` often make them depend on temperature and
pressure; the model takes their values at the conditions of the calculation.

[HollandPowell2003](@citet) write the activity coefficients (their Eq. 2) with
interactions rescaled by the size of the end-member considered,
``W^*_{ij} = 2\alpha_k W_{ij}/(\alpha_i + \alpha_j)``, their double sum running over
``j > i``, the range the excess energy implies. A shorter
route goes through the excess itself. Substituting ``\varphi``,

```math
G^{\text{ex}} = \frac{Q(\mathbf{x})}{A(\mathbf{x})},
\qquad
Q = \sum_{i<j} x_i x_j w_{ij},
\quad
A = \sum_l \alpha_l x_l,
\quad
w_{ij} = \frac{2\alpha_i\alpha_j W_{ij}}{\alpha_i + \alpha_j},
```

a quadratic form over a linear one. The excess of ``n_k`` moles is then
``N G^{\text{ex}}(\mathbf{n}/N) = Q(\mathbf{n})/A(\mathbf{n})``, homogeneous of
degree one in the amounts, and its derivative is

```math
RT \ln\gamma_k = \frac{\partial}{\partial n_k}\frac{Q(\mathbf{n})}{A(\mathbf{n})}
 = \frac{\sum_{j\neq k} x_j w_{kj}}{A} - \frac{\alpha_k\, Q}{A^2} ,
```

which is `_excess_ln_gamma(::VanLaarModel, …)`. With equal sizes ``\alpha``,
``w_{ij} = \alpha W_{ij}`` and ``A = \alpha``, and the regular expression of
section 3 comes back; for a binary it is Eqs. (4) and (5) of
[HollandPowell2003](@citet),
``RT\ln\gamma_1 = 2\alpha_1/(\alpha_1 + \alpha_2)\,\varphi_2^2 W_{12}``. The test
suite checks it against those equations for the alkali feldspar of the paper,
against its Eq. (13) for calcite, magnesite and dolomite, with the paper's
parameters, and ``\sum_k x_k RT\ln\gamma_k`` against ``G^{\text{ex}}``.

Neither this model nor the series above has the Redlich-Kister form of section
4 with three coefficients, so their binary is scanned for a spinodal and its
common tangent found from the excess itself, ``\sum_k x_k \ln\gamma_k``, in the
mole fraction of the first end-member.

### Convexity with more than two end-members

A one-dimensional scan decides a binary. With ``n`` end-members the question is
whether the Hessian of ``g/RT`` is positive on the tangent space of the simplex,
``\{\mathbf{d} : \sum_i d_i = 0\}``, at every composition, and
[`mixing_convexity`](@ref) answers it by model.

For a regular model, ``g/RT = \sum_i x_i \ln x_i + \sum_{i<j} w_{ij} x_i x_j`` with
``w_{ij} = W_{ij}/RT``, and the Hessian is ``\operatorname{diag}(1/\mathbf{x}) + \mathbf{W}/RT``. The first term
is at least 2 on the tangent space: for a unit ``\mathbf{d}`` with ``\sum_i d_i = 0``, the
Cauchy–Schwarz inequality gives
``\sum_i d_i^2/x_i \ge (\sum_i |d_i|)^2 / \sum_i x_i = (\sum_i |d_i|)^2``, and the
ℓ1 norm of such a ``\mathbf{d}`` is at least ``\sqrt 2`` (its positive and negative parts
have equal sums). Hence

```math
\lambda_{\min}\!\left(\mathbf{Q}^\mathsf{T}\,\frac{\mathbf{W}}{RT}\,\mathbf{Q}\right) \ge -2
\quad\Longrightarrow\quad \text{convex},
```

``\mathbf{Q}`` an orthonormal basis of the tangent space. For two end-members this is
``W \le 2RT``, the threshold of section 3; for more it is a sufficient
condition. In the other direction, a pair with ``W_{ij} > 2RT`` is a witness:
at the middle of that edge the second derivative along it is ``4 - 2w_{ij} < 0``.
Between the two, the smallest eigenvalue of the projected Hessian is searched on
a lattice of compositions, and a negative one is a witness too; finding none
proves nothing, and the verdict is then `:undecided`. The subregular, the
Redlich-Kister and the van Laar models are decided by that search alone beyond
two end-members: their Hessian depends on the composition through the
asymmetric and the ternary terms, and no bound of the kind above is written for
them.

```@example ss_convexity
using ChemistryLab, DynamicQuantities, Printf
RT = ChemistryLab.R_GAS * 298.15
W(a, b, c) = [0.0 a b; a 0.0 c; b c 0.0] .* RT
[mixing_convexity(RegularSolutionModel(W(1.9, 1.9, 1.9)), 3).verdict,
 mixing_convexity(RegularSolutionModel(W(3.0, 0.0, 0.0)), 3).verdict]
```

A phase whose verdict is `:nonconvex` is refused at construction, as a concave
binary is, and admitted with two instances. The certificate of an answer is not
scoped `:global_minimum` unless every mixing phase is proved convex; it is then
`:local_minimum` where the Hessian over the directions that conserve matter is
positive definite at the answer, and `:kkt_point` otherwise.

### A member that is a mixture of two others

Two different things can make one end-member the mixture of two others in
composition. In an **ordered** member the difference of Gibbs energy is the point:
the pentameric T5C of CSH3T is the average of TobH and T2C less an ordering
energy [Kulik2011; Eq. 18](@cite), and the siliceous hydrogarnet
C3AS0.41H5.18 lies between C3AH6 and C3AS0.84H4.32. When the Gibbs energy is
the average as well, the phase holds one substance twice, once as a member and
once as a mixture of two, and ideal mixing counts its configurations twice.
[`SolidSolutionPhase`](@ref) warns when that holds within ``0.1\,RT``; the
MgAl-OH-LDH ternary of Cemdata18 is the case, its M6A member the average of M4A
and M8A:

```@example ss_convexity
subs = Dict(symbol(s) => s for s in build_species(datapath("cemdata18-thermofun.json"); verbose = false))
g(s) = ustrip(us"kJ/mol", subs[s][:ΔₐG⁰](T = 298.15u"K", P = 1.0e5u"Pa"; unit = true))
@printf("MgAl-OH-LDH  M6A - (M4A + M8A)/2  = %+.4f kJ/mol\n",
        g("M6A-OH-LDH") - (g("M4A-OH-LDH") + g("M8A-OH-LDH")) / 2)
@printf("CSH3T        T5C - (TobH + T2C)/2 = %+.2f kJ/mol\n",
        g("CSH3T-T5C") - (g("CSH3T-TobH") + g("CSH3T-T2C")) / 2)
```

The LDH difference, 0.3 J/mol, is about 10⁻⁴ RT, where the −4.35 kJ/mol of T5C
is an ordering energy. The shipped LDH entry keeps the published model with
`acknowledge_degenerate = true`.

## 10. How a solid solution is declared

A [`SolidSolutionPhase`](@ref) names its end-members and carries a model:

```julia
phase = SolidSolutionPhase("CSHQ", [sp["CSHQ-JenD"], sp["CSHQ-JenH"],
                                    sp["CSHQ-TobD"], sp["CSHQ-TobH"]];
                           model = IdealSolidSolutionModel())
```

Every end-member must be `AS_CRYSTAL`, and the phase is checked against its
model's arity **and against the convexity of its mixing energy** at construction;
`instances` is how a declaration asks for more than one coexisting composition,
which section 6 is about. Sets of end-members with their interaction
parameters are read from a TOML through
[`build_solid_solutions`](@ref), which is also the only interaction-parameter
file format in the package — `data/solid_solutions.toml` is the shipped example.

Two consequences to keep in mind when reading results:

  - a `ChemicalSystem` groups end-members into `ss_groups`, and **every** aqueous
    activity model must call the solid-solution branch or the end-members come
    back with ``\ln a = 0``, i.e. treated as pure phases;
  - an end-member's [`saturation_indices`](@ref) entry is relative to its
    **current** mole fraction, since its activity is ``x_k\gamma_k``. For an
    end-member sitting at the solver's lower bound that is a statement about a
    vanishing phase, not about whether the solid solution would form.

## See also

  - [Thermochemistry](@ref sec-theory-thermo) — where ``\ln a`` enters ``\mu``
  - [Activity models](@ref sec-theory-activity) — the aqueous half
  - [Solving an equilibrium](@ref sec-solving) — declaring and solving with
    solid solutions
