# [Thermochemistry: the quantities and how they connect](@id sec-theory-thermo)

[Energies, enthalpies and the chemical potential](@ref sec-theory-basics) ends
with the chemical potential of a species written through the apparent Gibbs
energy a database stores. This page follows that quantity into the code: the
potential the solver minimizes, the models by which its standard part is
evaluated at any temperature and pressure, the minimization and the multipliers
it returns, and the quantities read from them.

## 1. The potential the solver minimizes

The terms of the elements and of the charge that the chemical potential contains
besides the apparent Gibbs energy change neither an equilibrium composition nor
a heat computed between two states
([The chemical potential and the apparent Gibbs energy](@ref sec-theory-mu-apparent)),
and the package leaves them out. It works with the remaining potential divided by
``RT``,

```math
g_i \;=\; \frac{\Delta_a G_i^\circ(T,P)}{RT} + \ln a_i
      \;=\; \underbrace{\texttt{ΔₐG⁰overRT[i]}}_{\text{database}}
      \;+\; \underbrace{\texttt{lna(n,p)[i]}}_{\text{activity model}} ,
```

which the rest of the chapter calls the chemical potential divided by ``RT``,
the difference shifting the component potentials of §4 by constants and leaving
compositions and heats unchanged. The two parts come from different
places: the standard part belongs to the species and is read from the database
(§3), whereas the logarithm of the activity depends on the composition of the
phase and is computed by the [activity model](@ref sec-theory-activity).
[`build_potentials`](@ref) returns their sum as a closure,

```julia
μ(n, p) = p.ΔₐG⁰overRT .+ lna(n, p)      #  μᵢ/RT  =  ΔₐG⁰ᵢ/RT  +  ln aᵢ
```

The division by ``RT`` makes the potential dimensionless, which is the form the
minimizer differentiates, and brings its standard part and the logarithm of the
activity to the same scale.

## 2. Activities, standard states, and the conventions in use

An activity is a ratio to a standard state, and the choice of the standard state
is a convention that must be stated for the numbers to have a meaning. The
conventions of the package are summarized below and argued, class by class, in
[Standard states](@ref sec-theory-standard-states):

| class | activity | standard state | convention |
|:--|:--|:--|:--|
| aqueous solute | ``a_i = \gamma_i\, m_i/m^\circ`` | ``m^\circ = 1`` mol/kg of solvent | molality |
| aqueous solute (ideal model) | ``a_i = c_i/c^\circ`` | ``c^\circ = 1`` mol/L | molarity |
| solvent (water) | ``a_w`` | pure water | mole fraction |
| pure solid, pure phase | ``a_i = 1`` ⇒ ``\ln a_i = 0`` | the pure substance | — |
| gas in an ideal mixture | ``a_i = x_i\,P/P_r`` | the pure ideal gas at ``P_r = 1`` bar | mole fraction |
| solid-solution end-member | ``a_k = \gamma_k x_k`` | the pure end-member | mole fraction |

[`concentration_scale`](@ref) is how a model declares which of the first two it
uses, and the aqueous accessors read it rather than guessing:

```@example thermo
using ChemistryLab
concentration_scale(DiluteSolutionModel()), concentration_scale(HKFActivityModel())
```

An activity coefficient is thus a property of a species in a given solution
rather than of the species alone, and it tends to one in the limit that defines
the standard state, infinite dilution for a solute and purity for a solid.

## 3. How ``\Delta_a G^\circ(T,P)`` is evaluated

The apparent Gibbs energy of a species is its value of formation at
``(T_r, P_r)`` plus integrals of its heat capacity and of its volume
([Apparent and formation Gibbs energies](@ref sec-theory-apparent)). A database
supplies, for each species, the reference values and a model for these
integrals, and three models cover the species of the databases the package
reads.

Minerals and gases carry a heat-capacity polynomial, the model `:cp_ft_equation`,
of up to eleven terms in powers of ``T``, ``\sqrt{T}`` and ``\ln T``.
`THERMO_MODELS` (`src/thermodynamics/thermo_models.jl`) stores ``C_p(T)``
together with the analytic antiderivatives of the same coefficients for ``S``,
``H`` and ``G``, so that the relations

```math
C_p = \left(\frac{\partial H}{\partial T}\right)_P ,
\qquad
S(T) = S(T_r) + \int_{T_r}^{T}\frac{C_p}{T'}\,\mathrm{d}T' ,
\qquad
G = H - TS
```

hold by construction rather than numerically. The model has no pressure term of
its own; a record that declares a constant molar volume adds ``V^\circ (P - P_r)``
to its Gibbs energy and enthalpy, as [Standard states](@ref sec-theory-standard-states)
§1 explains. A heat capacity given on several temperature intervals separated by
phase transitions is followed into each of them: the functions are anchored on
the interval that contains ``T_r``, carried continuously to the next, and a
transition the record places at a boundary adds its enthalpy ``\Delta H_t`` and
its entropy ``\Delta S_t``, the Gibbs energy staying continuous; the molar
volume stays the record's at ``T_r``. An entry that gives a
single heat capacity at ``T_r`` is extrapolated with the same model reduced to its
constant term; in CEMDATA18, the solvent, the zeolites and the magnesium silicate
hydrates are among them. An entry that gives none, but gives the entropy, is
extrapolated with a zero heat capacity, so that its apparent Gibbs energy still
decreases as ``-S^\circ`` with temperature.

Aqueous solutes use the model `:solute_hkf88_reaktoro`, the
Helgeson-Kirkham-Flowers (HKF) equation of state [Helgeson1981](@cite) in its
revised form [TangerHelgeson1988](@cite), whose standard Gibbs energy adds two
contributions of different origins to the reference value
[AndersonCrerar1993](@cite) (§17.9). The first, called nonsolvation, describes
the species itself through a heat capacity (coefficients ``c_1``, ``c_2``) and a
volume (coefficients ``a_1`` to ``a_4``) that depend on ``T`` and ``P``, the
pressure entering through ``P - P_r`` and ``\ln[(\Psi + P)/(\Psi + P_r)]`` with
``\Psi = 2600`` bar. The second is the energy of solvation of a charge in a
dielectric continuum, given by the Born equation

```math
\Delta G_{\text{s}} = \omega\left(\frac{1}{\varepsilon_r} - 1\right) ,
```

where ``\varepsilon_r`` is the relative permittivity of water at ``T`` and ``P`` and
``\omega`` the Born coefficient of the species, itself a function of ``T`` and
``P`` for an ion. The code evaluates it through the Born function
``Z = -1/\varepsilon_r`` as ``-\omega(Z+1)``, measured from its value at the
reference conditions, a term that follows the permittivity of water, which falls
from about 78 at 25 °C to about 55 at 100 °C; the same model of water supplies the
parameters of the activity models. The standard state of a solute is thus taken
at the pressure of the system, and the parameters of an ion are conventional,
those of H⁺ being zero at every temperature and pressure.

For the Maier-Kelley heat capacity, ``C_p^\circ = a_0 + a_1 T + a_2 T^{-2}``,
both temperature integrals are elementary, and the closed form can be set
against the function the package builds from the same data. With the parameters
of calcite used in [Thermodynamic Functions](@ref sec-thermodynamics), at 500 K:

```@example thermo
using DynamicQuantities
Tr, T = 298.15, 500.0
# Thermoddem's calcite, as transcribed in data/literature/Blanc2012.json (SI units)
cal = literature_row("Blanc2012", "individual_properties", "Calcite")
mk = literature_row("Blanc2012", "maier_kelley", "Calcite")
ΔfG, S = ustrip(cal.dfG), ustrip(cal.S)
a₀, a₁, a₂ = ustrip(mk.a), ustrip(mk.b), ustrip(mk.c)

∫Cp = a₀ * (T - Tr) + a₁ / 2 * (T^2 - Tr^2) - a₂ * (1 / T - 1 / Tr)        # ∫ Cp dT
∫Cp_T = a₀ * log(T / Tr) + a₁ * (T - Tr) - a₂ / 2 * (1 / T^2 - 1 / Tr^2)   # ∫ Cp/T dT
closed_form = ΔfG - S * (T - Tr) + ∫Cp - T * ∫Cp_T

dtf = build_thermo_functions(
    :cp_ft_equation,
    Dict(
        :S⁰ => cal.S, :ΔₐH⁰ => cal.dfH, :ΔₐG⁰ => cal.dfG,
        :a₀ => mk.a, :a₁ => mk.b, :a₂ => mk.c, :T => Tr * u"K",
    ),
)
(code = dtf[:ΔₐG⁰](T = T), closed_form = closed_form)
```

This is the quantity every species carries as `ΔₐG⁰`, and the one the minimizer
reads, divided by ``RT``, as `ΔₐG⁰overRT`; the subscript ``a`` is not a variant
spelling of ``f``. A value of ``\Delta_f G^\circ(T)`` read from a table built in
the traditional convention cannot be combined with the apparent energies of a
database, since the two differ by ``\sum_e \alpha_{ei}\,[G_e^\circ(T) -
G_e^\circ(T_r)]`` and this difference no longer cancels between species taken
from different sources. A second apparent convention, due to Berman and Brown,
also removes the elemental entropies at ``T_r``
[AndersonCrerar1993](@cite) (§7.4.2); its values differ from the former by the
constant ``T_r\sum_e \alpha_{ei}\, S_e^\circ(T_r)``, which is consistent within
one database and inconsistent across two. Which convention the code implements
can be read on the formula itself: the anchor ``\Delta_a G_i^\circ(T_r) =
\Delta_f G_i^\circ(T_r)`` rules out the values of Berman and Brown, and the
absolute entropy in the linear term rules out the traditional ones.

## 4. Equilibrium is a constrained minimization, and its dual is the useful part

At fixed ``T`` and ``P``, equilibrium is the composition that minimizes ``G``
subject to the matter available, which is a linear constraint, since the
elements and the charge are conserved whatever the reactions:

```math
\min_{\mathbf{n} \ge 0} \; \sum_i n_i\, g_i(\mathbf{n})
\qquad\text{subject to}\qquad
\mathbf{A}\,\mathbf{n} = \mathbf{b} ,
```

where ``\mathbf{A}`` is the conservation matrix `cs.SM.A`, with one row per
component and one column per species, and ``\mathbf{b}`` the budget of the components.
The components are the primary species of the system, or the elements with the
charge, and the two choices express the same constraints
([Formation from primary species](@ref sec-theory-primaries)). This is the
formulation of [Leal2017](@citet), in which no list of reactions is needed: the
reactions are the moves of ``\mathbf{n}`` within the null space of ``\mathbf{A}``.

The Lagrange multipliers of the equality constraints are the useful output.
Writing ``y_c`` for the multiplier of row ``c``, and
``u_s = -\sum_c A_{cs}\, y_c`` for the potential they give species ``s`` (the
components of ``\mathbf{u} = -\mathbf{A}^\mathsf{T}\mathbf{y}``), the first-order
conditions are

```math
g_s \;=\; u_s \quad\text{for every species present},
\qquad
g_s \;\ge\; u_s \quad\text{for every species absent.}
```

The ``-y_c`` are the component potentials, or element potentials when the
components are elements: one number per component, from which the chemical
potential of any species follows as a scalar product,
``\mathbf{u} = -\mathbf{A}^\mathsf{T}\mathbf{y}``. The two lines above are
the complementarity conditions of the Karush-Kuhn-Tucker (KKT) system, and they
are what [`optimality_certificate`](@ref) checks; the certificate itself, and the
conditions under which it proves a global minimum, are the subject of
[Proving that an answer is the answer](@ref sec-theory-certificate).

## 5. The saturation index

The column of ``\mathbf{A}`` for a species ``s`` is its reaction of formation
from the primary species, ``\sum_c A_{cs}\,\mathrm{B}_c \rightarrow s``, whose Gibbs
energy is ``\Delta_r G = \Delta_r G^\circ + RT\ln Q_r`` and vanishes at equilibrium
([Reactions and equilibrium constants](@ref sec-theory-reactions)). For a solid, the
reverse reaction is its dissolution into the primaries, ``Q_r`` of that dissolution
is the ion activity product IAP, ``K`` the solubility product, and the logarithm
of their ratio the saturation index. Since the rows of ``\mathbf{A}`` are
labeled by primary species, a component potential is the potential ``g`` of
the corresponding primary, and the index becomes a difference of potentials,
with no equilibrium constant to look up,

```math
\mathrm{LogSI}_s \;=\; \frac{u_s - g_s}{\ln 10}
\;=\; \log_{10}\frac{\mathrm{IAP}}{K_{sp}}
\;=\; -\frac{\Delta_r G}{RT\ln 10} ,
```

``\Delta_r G`` being that of the reaction of formation. A negative index denotes
an undersaturated phase, a zero index a phase in equilibrium with the solution,
and a positive index a phase that should have precipitated; this is what
[`saturation_indices`](@ref) returns.

Two consequences of the identity make the index a check rather than a
convention. Since ``K_{sp}`` never enters the computation, being implied by the
standard potentials, an index cannot disagree with the ``\Delta_a G^\circ`` used
by the rest of the calculation. And since every phase present at an equilibrium
satisfies the first line of the conditions of §4, its index must vanish to the
tolerance of the solve, failing which the state is not an equilibrium and no
other index of the result has a meaning. The test suite
(`test/aqueous_properties.jl`) recomputes the index by hand from this formula and
compares it with the function.

## 6. Volume, porosity and chemical shrinkage

Volumes are treated as ideal: the volume of a phase is the sum of the standard
molar volumes of its species,

```math
V = \sum_i n_i V_i^\circ(T,P) ,
```

with no excess volume of mixing. This assumption underlies every porosity the
package reports. [`volume`](@ref) returns the split by aggregate state,
[`porosity`](@ref) the void fraction relative to a reference state, and
[`chemical_shrinkage`](@ref) the volume the reaction itself consumes, the
hydrates occupying less than the water and the clinker they were made from,
which is the mechanism of self-desiccation
([Self-desiccation](@ref sec-self-desiccation)).

A species carrying no ``V^\circ`` contributes zero and would corrupt a porosity
in silence, which is why the volume machinery reports
`missing_molar_volumes(state)` and why [`CapillaryWater`](@ref) refuses to be
constructed when that list is not empty.

## 7. What is assumed, and what is not

| assumed | not assumed |
|:--|:--|
| one well-mixed phase per aggregate state | a list of reactions, a reaction path, a sequence |
| ideal molar volumes, no excess volume | ideal activities, which are the business of the activity model |
| the domain of validity of the activity model | the phases present, which are a result |
| that ``\mathbf{A}\mathbf{n} = \mathbf{b}`` is the whole of the conservation | that a point satisfying the optimality conditions is the minimum |

The last cell rests on the certificate. The optimality conditions prove a global
minimum when the potentials derive from one convex Gibbs energy, which holds for
some activity models and not for others, and fails as well for a solid solution
declared inside a miscibility gap or for a constraint that shifts the activity of
water; the certificate then reports a KKT point or a composition consistent with
its own activities, as stated in
[What the certificate proves, and when](@ref sec-theory-certificate-scope).

## Where to go next

[Standard states](@ref sec-theory-standard-states) is the next page of the
chapter and states, class by class, what each activity of §2 is measured from.
The other pages this one leads to are:

  - [Proving that an answer is the answer](@ref sec-theory-certificate), the
    certificate built on the conditions of §4;
  - [Activity models](@ref sec-theory-activity), the ``\ln a_i`` part of ``g_i``;
  - [Solid solutions](@ref sec-theory-solid-solutions), the part due to mixing
    in a solid;
  - [Chemical Equilibrium](@ref sec-equilibrium), which drives the solver.
