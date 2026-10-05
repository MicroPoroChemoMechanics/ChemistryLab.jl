# [Kinetics under partial equilibrium: the right-hand side and the calorimeters](@id sec-theory-pe-kinetics)

!!! info "Before this page"
    [Coupling kinetics and equilibrium](@ref sec-coupling), which introduces the
    partition and the state ``(\mathbf{b}_e, \mathbf{n}_k)``, and
    [Energies, enthalpies and the chemical potential](@ref sec-theory-basics),
    where the enthalpy of a state is defined.

A kinetic run under partial equilibrium integrates the slow reactions and solves
the fast ones by a Gibbs energy minimization. The formulation of
[Leal2015](@citet) makes the right-hand side of the differential system a
function of its state alone, the minimization being carried out inside that
function, and gives its Jacobian exactly by the implicit-function theorem; any
integrator can then be used, the stiff ones included, and its error control
sees the coupled system as a whole. This page writes that formulation as the
package implements it, and extends it to the energy balance of a calorimeter,
where the temperature, which the minimization depends on, is itself an unknown.

## 1. The state and the right-hand side

The species are split into a kinetic partition, of amounts ``\mathbf{n}_k``,
and an equilibrium partition, of amounts ``\mathbf{n}_e``; the kinetic
reactions, of rates ``\mathbf{r}`` and stoichiometric matrix
``\boldsymbol{\nu} = [\boldsymbol{\nu}_e\;\boldsymbol{\nu}_k]``, are the only
ones that move the budget ``\mathbf{b}_e = \mathbf{A}_e\mathbf{n}_e`` of the
equilibrium partition. The state is

```math
\mathbf{u} = (\mathbf{b}_e,\ \mathbf{n}_k,\ \boldsymbol{\xi}),
```

``\boldsymbol{\xi}`` being the extents of the kinetic reactions, and it evolves
according to [Leal2015; Eqs. 2.25–2.26](@cite)

```math
\frac{\mathrm{d}\mathbf{n}_k}{\mathrm{d}t} = \boldsymbol{\nu}_k^\mathsf{T}\mathbf{r},
\qquad
\frac{\mathrm{d}\mathbf{b}_e}{\mathrm{d}t} = \mathbf{A}_e\,\boldsymbol{\nu}_e^\mathsf{T}\mathbf{r},
\qquad
\frac{\mathrm{d}\boldsymbol{\xi}}{\mathrm{d}t} = \mathbf{r},
\qquad
\mathbf{r} = \mathbf{r}(\mathbf{n}_e, \mathbf{n}_k, T, t),
```

the amounts of the equilibrium partition being those of the minimum of the Gibbs
energy at the budget and the temperature of the state,

```math
\mathbf{n}_e = \varphi(\mathbf{b}_e, T)
= \arg\min_{\mathbf{n}\,\ge\,0}\ G(\mathbf{n}; T, P)
\quad\text{subject to}\quad \mathbf{A}_e\,\mathbf{n} = \mathbf{b}_e .
```

The right-hand side ``\mathbf{f}(\mathbf{u}, t)`` is thereby a function of the
state, since ``\varphi`` is: the same state always gives the same derivative,
whichever point of the trajectory the integrator evaluates it at, a stage, a
rejected step or a Jacobian probe.

It is evaluated as written whenever a rate law reads the equilibrium partition,
through an activity, a saturation ratio or the amount of an equilibrium species,
and whenever a semi-adiabatic calorimeter is coupled (Section 3): the
minimization is solved by the certified solver, warm-started from the last
accepted step (`speciation = :rhs`). When no rate law reads it, as for a
Parrott–Killoh or a Waller law, ``\mathbf{f}`` does not depend on ``\varphi`` at
all, so that freezing the partition within a step changes nothing, and the
minimization is solved once per accepted step only to report it
(`speciation = :frozen`), by the interior point; [`speciated_states`](@ref)
certifies it at any instant afterwards. `integrate` tells the two cases apart
from the rate laws themselves.

## 2. The Jacobian

A stiff method needs ``\mathbf{J} = \partial\mathbf{f}/\partial\mathbf{u}``.
Applying the chain rule to the rates gives [Leal2015; Eqs. 2.34–2.35](@cite)

```math
\frac{\partial\mathbf{r}}{\partial\mathbf{u}}
= \frac{\partial\mathbf{r}}{\partial\mathbf{n}_e}\,\frac{\partial\varphi}{\partial\mathbf{b}_e}\,\frac{\partial\mathbf{b}_e}{\partial\mathbf{u}}
+ \frac{\partial\mathbf{r}}{\partial\mathbf{n}_k}\,\frac{\partial\mathbf{n}_k}{\partial\mathbf{u}},
```

in which only ``\partial\varphi/\partial\mathbf{b}_e`` is not explicit. It
follows from the optimality conditions of the minimization. Writing them
``\mathbf{F}(\mathbf{n}_e, \mathbf{y}; \mathbf{b}_e, T) = \mathbf{0}``, the
stationarity of the Lagrangian of the species present and the conservation
``\mathbf{A}_e\mathbf{n}_e = \mathbf{b}_e``, with ``\mathbf{y}`` the multipliers
of the conservation rows, and differentiating them at the answer with its active
set held yields

```math
\frac{\partial\mathbf{F}}{\partial(\mathbf{n}_e, \mathbf{y})}
\begin{bmatrix} \partial\mathbf{n}_e/\partial\mathbf{b}_e \\ \partial\mathbf{y}/\partial\mathbf{b}_e \end{bmatrix}
= -\frac{\partial\mathbf{F}}{\partial\mathbf{b}_e},
```

and the same with ``T`` in place of ``\mathbf{b}_e``. The package does not form
these matrices column by column. The minimization is solved on plain numbers,
and when the state carries the dual numbers of automatic differentiation the
answer is lifted by one Newton correction with the Jacobian of ``\mathbf{F}``
at the answer,

```math
\mathbf{n}_e^\star \;\leftarrow\; \mathbf{n}_e^\star
- \left[\frac{\partial\mathbf{F}}{\partial(\mathbf{n}_e,\mathbf{y})}\right]^{-1}
\mathbf{F}(\mathbf{n}_e^\star, \mathbf{y}^\star; \mathbf{b}_e, T),
```

whose value part vanishes and whose dual part is exactly the derivative above,
in every direction the dual numbers carry. The Jacobian of ``\mathbf{f}`` is thus
exact to first order, which is all a Rosenbrock method or the Newton iteration of
a BDF method requires, and it holds the derivative of the partition rather than
the zero a partition frozen within a step would contribute.

## 3. The calorimeters

The enthalpy of the paste is the sum of the standard enthalpies of its species,

```math
H(\mathbf{n}, T) = \sum_i n_i\, h_i^\circ(T),
```

the excess enthalpies of the activity models being left out, as
[Energies, enthalpies and the chemical potential](@ref sec-theory-basics) §3
states. Under partial equilibrium the kinetic reactions only release ions, the
hydrates being precipitated by the minimization, so that the heat of those
reactions is not the heat of the paste; the heat is the change of ``H`` over the
whole composition, ``\mathbf{n}_e = \varphi(\mathbf{b}_e, T)`` included. The
reference is the partition at the start, ``H_0 = H(\varphi(\mathbf{b}_e(0),
T_0), \mathbf{n}_k(0), T_0)``, so that the first equilibrium releases nothing.

### The isothermal cell

The bath holds the paste at ``T_0`` and takes the heat ``Q`` the paste releases,
so that the enthalpy of the paste and of what the bath has received does not
change, ``H + Q = H_0``. The heat is therefore a function of the state,

```math
Q(\mathbf{u}) = H_0 - H\bigl(\varphi(\mathbf{b}_e, T_0), \mathbf{n}_k, T_0\bigr),
```

and requires no integration. Its rate along the trajectory is the derivative of
the same function,

```math
\dot q = -\frac{\partial H}{\partial\mathbf{n}_e}\,\frac{\partial\varphi}{\partial\mathbf{b}_e}\,\frac{\mathrm{d}\mathbf{b}_e}{\mathrm{d}t}
- \frac{\partial H}{\partial\mathbf{n}_k}\,\frac{\mathrm{d}\mathbf{n}_k}{\mathrm{d}t},
```

evaluated by lifting the partition along ``\mathrm{d}\mathbf{u}/\mathrm{d}t``
([`heat_flow`](@ref)), and the certified replay
([`heat_release`](@ref)) computes both from the same partition.

### The semi-adiabatic cell

The vessel and what it holds besides the paste, of heat capacity ``C_v``, take
the temperature of the paste, and heat leaves through the walls at the rate
``\mathcal{L}(T - T_{\rm env})`` ([`SemiAdiabaticCalorimeter`](@ref)). The
enthalpy of the cell, ``H + C_v(T - T_0)``, changes by the losses alone, so that
the state is extended by its change since the start,

```math
\mathbf{u} = (\mathbf{b}_e,\ \mathbf{n}_k,\ \boldsymbol{\xi},\ \Delta H),
\qquad
\frac{\mathrm{d}\,\Delta H}{\mathrm{d}t} = -\mathcal{L}(T - T_{\rm env}),
```

and the temperature is the root of the energy balance

```math
\mathcal{R}(\mathbf{u}, T) = H\bigl(\varphi(\mathbf{b}_e, T), \mathbf{n}_k, T\bigr) - H_0 + C_v\,(T - T_0) - \Delta H = 0 .
```

The temperature is solved for, not integrated. Its derivative with respect to
``T`` at fixed state,

```math
C_{\rm eq} = \frac{\partial\mathcal{R}}{\partial T}
= C_v + \sum_i n_i\,\frac{\mathrm{d}h_i^\circ}{\mathrm{d}T}
+ \sum_{i \in e} h_i^\circ\,\frac{\partial\varphi_i}{\partial T},
```

is the heat capacity of the cell at equilibrium, that of the vessel and of the
paste at fixed composition, plus the heat the partition takes up as it shifts
with the temperature. It is positive wherever the balance has a unique root, and
the temperature is found by Newton's method,
``T \leftarrow T - \mathcal{R}(\mathbf{u}, T)/C_{\rm eq}``, from the temperature
of the last accepted step, each iteration solving the minimization at the
current ``T`` and lifting it in temperature to obtain ``C_{\rm eq}``; a balance
whose ``C_{\rm eq}`` is not positive makes the evaluation fail, and the
integrator rejects the step.

The implicit-function theorem applied to ``\mathcal{R}`` then gives the
derivatives of the temperature the Jacobian needs,

```math
\frac{\partial T}{\partial\mathbf{u}} = -\frac{1}{C_{\rm eq}}\,\frac{\partial\mathcal{R}}{\partial\mathbf{u}},
\qquad
\frac{\partial\mathcal{R}}{\partial\mathbf{b}_e} = \sum_{i \in e} h_i^\circ\,\frac{\partial\varphi_i}{\partial\mathbf{b}_e},
\quad
\frac{\partial\mathcal{R}}{\partial\mathbf{n}_k} = \mathbf{h}_k^\circ,
\quad
\frac{\partial\mathcal{R}}{\partial\,\Delta H} = -1,
```

and the package obtains them, as in Section 2, by a single correction lifted
into the dual numbers of the state, ``T \leftarrow T^\star - \bigl(\mathcal{R}(\mathbf{u},
T^\star) - \mathcal{R}^\star\bigr)/C_{\rm eq}``, which keeps the value ``T^\star``
at which the partition was solved. The partition at the state is then lifted in
``\mathbf{b}_e`` and in ``T`` together, so that every entry of the Jacobian holds
the first derivatives of ``\varphi`` and of the root, and nothing else.

That is the reason for carrying ``\Delta H`` rather than ``T``. Written on the
temperature, the balance would read
``C_{\rm eq}\,\mathrm{d}T/\mathrm{d}t = -\partial\mathcal{R}/\partial\mathbf{b}_e \cdot
\mathrm{d}\mathbf{b}_e/\mathrm{d}t - \dots - \mathcal{L}``, a right-hand side that
holds ``\partial\varphi/\partial\mathbf{b}_e`` and ``\partial\varphi/\partial T``
themselves, so that its Jacobian would need the second derivatives of the
minimization. Following the derivatives of a partition frozen at the last
accepted step instead, and correcting the temperature by the heat of the next
re-speciation, makes the right-hand side depend on the history of the
integration and puts part of the energy balance outside its error control.

The heat and the temperature the accessors report are functions of the state:
[`temperature_profile`](@ref) solves the root at each instant,
[`cumulative_heat`](@ref) returns ``Q = C_v(T - T_0) - \Delta H``, the enthalpy
the paste has lost, and [`heat_flow`](@ref) its rate,
``\dot q = C_v\,\mathrm{d}T/\mathrm{d}t + \mathcal{L}(T - T_{\rm env})``.

### Without an equilibrium partition

When the kinetic reactions produce the hydrates themselves and no minimization is
solved, every amount follows from the extents,
``\mathbf{n} = \mathbf{n}(0) + \boldsymbol{\nu}^\mathsf{T}\boldsymbol{\xi}``, independently of the temperature, and
the same balance reduces to an equation the state can carry directly,

```math
\Bigl(C_v + \sum_i n_i\,c_{p,i}^\circ(T)\Bigr)\frac{\mathrm{d}T}{\mathrm{d}t}
= \sum_j r_j\,\bigl(-\Delta_r H_j^\circ(T)\bigr) - \mathcal{L}(T - T_{\rm env}),
```

whose right-hand side holds no derivative of a minimization; this formulation
integrates the temperature, and an isothermal cell the heat
``\mathrm{d}Q/\mathrm{d}t = \sum_j r_j(-\Delta_r H_j^\circ)``.

## 4. What is checked

The tests of the calorimeters (`test/kinetics/test_calorimetry.jl`) hold the
formulation to identities that a route sharing nothing with it computes:

- the enthalpy of an adiabatic cell, ``H + C_v(T - T_0)`` evaluated on the
  certified replay at the temperature the run solved, stays ``H_0`` to within
  ``10^{-7}`` of the heat released;
- ``1/(\partial T/\partial\,\Delta H)``, taken through the lifted root, equals
  ``C_{\rm eq}`` computed with the certified equilibrium differentiated with
  respect to its temperature, to every printed digit, and
  ``\partial T/\partial\mathbf{b}_e`` the same with respect to the budget;
- the heat of an isothermal cell equals the enthalpy difference of the certified
  replay, and its rate the rate of the replay;
- a Rosenbrock method, a BDF method and an explicit Runge–Kutta method integrate
  the same temperature, to ``10^{-5}`` K, with losses through the walls equal to
  ``-\Delta H``.

## 5. Limits

- The excess enthalpies of the activity models are not in ``H``, whereas the
  minimization carries the temperature dependence of the activity coefficients;
  ``C_{\rm eq}`` then holds a term that the excess enthalpies would have made a
  quadratic form, and its positivity is not guaranteed by the stability of the
  equilibrium alone. A balance with a non-positive ``C_{\rm eq}`` is refused at
  the evaluation rather than corrected.
- ``\varphi`` is smooth only while the assemblage does not change. Where a phase
  appears or vanishes it has a kink, which the integrator meets as a loss of
  accuracy and answers by shortening its step; no event is located there.
- The partition is the answer of the certified solver to its tolerance, and the
  temperature is solved to ``10^{-9}`` K.

## See also

- [Coupling kinetics and equilibrium](@ref sec-coupling), the partition and why
  the state carries ``\mathbf{b}_e``.
- [A CEM I 52.5 N mortar in a semi-adiabatic calorimeter, inside the kinetics](@ref ex-semiadiabatic),
  the semi-adiabatic cell on a cement.
- [Writing a kinetic model](@ref sec-kinetics-syntax), the syntax of the
  calorimeters.
