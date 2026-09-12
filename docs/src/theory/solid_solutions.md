# [Solid solutions](@id sec-theory-solid-solutions)

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
does: `_solid_solution_lna!` writes `log(x[k] + ϵ) + _excess_ln_gamma(model, k, x, T)`
into the log-activity vector, so the two terms are visibly separate and a model
supplies only the second.

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

``W_{ij}`` (J/mol) is the energy of an ``i``–``j`` contact *relative to the
average* of ``i``–``i`` and ``j``–``j``, so its sign is physical:

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

This matters in practice because **nothing detects it**. A solid solution is
entered as one phase, the activity expression goes on returning values past the
threshold, and those values describe a metastable single phase. Checking
``W/RT`` against 2 is the caller's job.

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
end-members the choices are the ideal model or `RegularSolutionModel` with a
full ``W`` matrix.

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

Which pair? The two compositions ``x_\alpha < x_\beta`` at which the chemical
potentials of *both* end-members agree:

```math
\left.\frac{\mathrm{d}g}{\mathrm{d}x}\right|_{x_\alpha}
  = \left.\frac{\mathrm{d}g}{\mathrm{d}x}\right|_{x_\beta}
  = \frac{g(x_\beta) - g(x_\alpha)}{x_\beta - x_\alpha} ,
```

which says geometrically that one straight line is tangent to the curve at both
points — the **common tangent**. Between ``x_\alpha`` and ``x_\beta`` the
equilibrium state is a mixture of the two, in the proportions the lever rule
gives, and the energy follows the tangent line rather than the curve.

The common tangent construction is wider than the spinodal: the *binodal*
``[x_\alpha, x_\beta]`` contains the spinodal, and between the two the phase is
metastable rather than unstable. A minimization sees only the tangent.

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

The test that separates them is Michelsen's: for a trial composition
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

### The same criterion, in another code

GEM-Selektor's `PhaseSelection` computes, for every phase, a stability index

```math
\Lambda_k = \log_{10}\Omega_k , \qquad
\Omega_k = \sum_{j\in l_k}\pi_j , \qquad
\pi_j = \frac{\omega_j(\hat{\mathbf{u}})}{\gamma_j(\hat{\mathbf{n}})} ,
```

the activity from the **dual** solution divided by the activity coefficient from
the **primal** one — an estimate of the mole fraction. That is term for term what
this package computes, in `_repair_start`'s ``\Omega = \sum_i 10^{SI_i}`` (the
ideal case ``\gamma = 1``) and in `phase_split_measure`'s log-sum-exp of
``d_j = u_i - g_i - \ln\gamma_i``. The GEMS3K paper calls ``\Omega_k``
*"a generalization of the saturation index"* and notes that it derives from the
KKT conditions [Kulik2013](@cite).

So the criterion is one object, arrived at independently. What differs is only
the declaration: CEMDATA18 ships the AFm and AFt binaries under two names each,
so a GEMS user represents a gap by declaring the binary twice in the database;
`instances = 2` asks for the same thing with a keyword. **Neither code splits a
phase by itself.**

The executed counterpart of this section is
[the miscibility-gap page](@ref ex-miscibility-gap), which runs one cement three
ways and reports the certificate each time.

## 7. How a solid solution is declared

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
  - [Chemical Equilibrium](@ref sec-equilibrium) — declaring and solving with
    solid solutions
