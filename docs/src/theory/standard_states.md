# [Standard states](@id sec-theory-standard-states)

!!! info "Before this page"
    [Thermochemistry](@ref sec-theory-thermo), whose notation is used
    throughout, in particular the split ``g_i = \mu_i/RT = \mu_i^\circ/RT + \ln a_i``.

The chemical potential of a species is split in
[Thermochemistry](@ref sec-theory-thermo) into a standard potential, read from a
database, and an activity, computed by a model. The split is not unique: it rests
on a reference state of the species, its standard state, which the database and
the activity model must share and which is fixed by convention rather than by
physics. This page states the convention adopted by the package for each class
of species, the way to pass from one convention to another, and the one case in
which the reference energy ceases to be a convention and becomes a statement
about matter, namely a surface site whose budget follows a dissolving solid.

## 1. One potential, two ways of writing it

Once ``\mu_i^\circ`` is chosen, the relation

```math
\mu_i = \mu_i^\circ + RT\ln a_i
```

defines the activity, and conversely. The chemical potential itself is a property
of the system that does not depend on any convention, and at equilibrium it takes
the same value in every phase between which the species can be exchanged;
neither the activity nor the standard potential shares that property, only their
combination does. Carbon dioxide distributed between a gas phase and water thus
has one chemical potential and two activities, a fraction of the pressure in the
gas and a molality in the solution. The ratio of the two activities is the Henry
constant ``K_H``, and equating the two expressions of the potential shows that
``RT\ln K_H = \mu^\circ_{\text{gas}} - \mu^\circ_{\text{aq}}`` is nothing but
the difference between the two standard potentials
[AndersonCrerar1993](@cite) (§12.6).

A standard state is a state of the pure substance, or of a reference solution,
completely specified by its temperature, its pressure, its composition and its
state of aggregation. The phrase "25 °C and 1 bar", often met in place of a
standard state, fixes two of these four attributes and leaves the other two open
[AndersonCrerar1993](@cite) (§12.2). Nothing requires the state to be
realizable either: the standard state of a solute used below is a hypothetical
solution, whose interest lies in the fact that its properties are known, not in
the possibility of preparing it.

The temperature of the standard state is always that of the system. The relation
above results from integrating ``\mathrm{d}\mu_i = RT\,\mathrm{d}\ln a_i`` at
constant temperature, from the standard state up to the state considered, so
that both must be taken at the same ``T``; this is why ``\Delta_a G_i^\circ`` is
stored as a function of temperature and evaluated at the temperature of each
state (see [Apparent and formation Gibbs energies](@ref sec-theory-apparent))
rather than as a number.

Pressure is treated according to the model attached to the species. The
Helgeson-Kirkham-Flowers equation of state of aqueous solutes depends on ``P``,
so that the standard state of a solute is at the pressure of the system. A
record that declares a molar volume independent of temperature and pressure
(`mv_constant` in the ThermoFun databases: the solids, the solutes described by
a heat-capacity polynomial, and the solvent, whose equation of state is not
implemented) has its standard state at the pressure of the system as well, its
standard Gibbs energy and enthalpy carrying

```math
\int_{P_r}^{P} V_i^\circ\,\mathrm{d}P = V_i^\circ\,(P - P_r) ,
```

which vanishes at ``P_r = 1`` bar ([`P_STANDARD`](@ref)). For portlandite, with
``V^\circ \simeq 33\ \mathrm{cm^3/mol}``, it amounts to ``0.012\,RT`` at 10 bar
and to ``1.3\,RT`` at 1 kbar: negligible for a laboratory sample or a
structure, not for a deep reservoir. The compressibility of the solvent is
neglected in this term, which misses about 0.7 % of it at 300 bar. A gas is
referred to the pure ideal gas at ``P_r``, and the pressure enters its activity
(§2). A species built by hand keeps the functions it is given: a heat-capacity
polynomial alone carries no pressure term.

## 2. The conventions in use, class by class

### Solutes: the hypothetical ideal one-molal solution

The standard state of a solute is a solution at the molality ``m^\circ = 1``
mol/kg in which the solute would behave as it does at infinite dilution, that is
a solution obeying Henry's law up to one molal. Such a solution does not exist,
the interactions between ions being significant well below that concentration,
and its properties are obtained by extrapolating measurements made on dilute
solutions. The two parts of the choice have separate reasons
[AndersonCrerar1993](@cite) (§12.4.4 and §12.4.6). Referring the ideal behavior
to infinite dilution makes the activity coefficient tend to one in the very limit
where the measurements are made, so that every departure from ideality is carried
by ``\gamma_i`` alone; the value of one molal makes the activity numerically
equal to the molality in that limit, whereas any other value would add a
constant to every standard potential without any benefit. The activity then
reads

```math
\ln a_i = \ln\gamma_i + \ln\frac{m_i}{m^\circ},
\qquad
m_i = \frac{n_i}{n_w M_w} ,
```

``n_w`` being the amount of water and ``M_w`` its molar mass.
Here ``M_w`` is expressed in kg/mol, so ``m_i`` is in mol/kg of **solvent**,
not mol/g or mol/kg of solution. Molarity instead uses solution volume, in
mol/L. At fixed composition, thermal expansion changes molarity but not
molality; reaction or water exchange can change either.

[`HKFActivityModel`](@ref), [`DaviesActivityModel`](@ref) and the Pitzer model
compute ``\gamma_i`` from the composition of the solution
([Activity models](@ref sec-theory-activity)). The dilute model takes
``\gamma_i = 1`` and forms the same ratio ``n_i/(n_w M_w)``, which it declares to
the aqueous accessors as a molarity, on the ground that a dilute solution has a
density close to 1 kg/L; [`concentration_scale`](@ref) says which reading a model
adopts.

### The solvent: pure water

Water is referred to the pure liquid at the temperature of the system, with
Raoult's law as the ideal limit, so that ``a_w \to 1`` as the solution tends to
pure water. The models obtain it from the osmotic coefficient ``\varphi``, which
is one for the ideal dilute model,

```math
\ln a_w = -M_w\,\varphi \sum_j m_j ,
```

which is the form imposed by the Gibbs-Duhem relation once the activity
coefficients of the solutes are given, as argued in
[Activity models](@ref sec-theory-activity) §3.

This is a consistency requirement for a thermodynamic model. Integrating a
chosen path does not by itself prove that arbitrary solute activity laws come
from one potential; the integrability checks in
[Activity models](@ref sec-theory-potential) establish the limits of that claim.

### Pure solids

A pure solid is its own standard state, its activity is one and its chemical
potential reduces to ``\mu_i^\circ``. It follows that a pure phase contributes to
the Gibbs energy a term linear in its amount, and that the minimization decides
whether it is present by comparing ``\mu_i^\circ`` with the combination of
component potentials its formula implies ([Thermochemistry](@ref sec-theory-thermo)
§4 and §5).

### Gases

A gas is referred to the pure ideal gas at ``P_r = 1`` bar, and its activity in
an ideal mixture is the ratio of its fugacity ``x_i P`` to that pressure,

```math
a_i = x_i\,\frac{P}{P_r} ,
```

so that its chemical potential grows with pressure as ``(\partial\mu_i/\partial
P)_T = RT/P``, the molar volume of an ideal gas. That is also the volume the
package gives a gas, whether its record declares the ideal gas (`mv_pvnrt`) or it
is built without a molar volume, so that the two are consistent at every
pressure.

### End-members of a solid solution

An end-member is referred to its own pure phase, as a pure solid is, and its
activity in an ideal solution is its mole fraction in the phase, ``a_k = x_k``;
the excess models add ``\ln\gamma_k``, with ``\gamma_k \to 1`` as the phase
becomes pure ([Solid solutions](@ref sec-theory-solid-solutions)). The mole
fraction depends on the way the formula of each end-member is written. Doubling a
formula halves the number of moles for a given mass, changes every mole fraction
of the phase and doubles the standard potential, so that database values cannot
be carried over to another formula unit. The mixing implemented here counts one
site per formula unit, which ties the size of the mole to the formula as written
in the database [AndersonCrerar1993](@cite) (§12.7).

### Species bound to a surface

A species occupying a site is referred to a surface entirely covered by it, and
its activity is its site fraction ``a_j = n_j/N``, ``N`` being the budget of the
site family ([Chemistry that happens on a surface](@ref sec-theory-surface) §3).
The reference density ``\Gamma^\circ`` fixing the zero of that scale is distinct
from the capacity ``\Gamma_C`` fixing how many sites exist [Kulik2002](@cite):
the former is a convention of the kind discussed on this page, the latter a
physical property of the sorbent.

## 3. Changing from one convention to another

Since ``\mu_i`` does not depend on the convention, passing from an old standard
state to a new one amounts to transferring a constant from one term to the other,

```math
\mu_i^{\circ,\text{new}} - \mu_i^{\circ,\text{old}}
  = RT\ln\frac{a_i^{\text{old}}}{a_i^{\text{new}}} ,
```

the right-hand side being evaluated in any state, and most conveniently in a
limit where both activity coefficients are known [AndersonCrerar1993](@cite)
(§12.6.1). The most frequent instance in aqueous chemistry is the passage between
the mole-fraction scale of a solute, used by formulations that treat the aqueous
phase as a mixture like any other, and the molality scale used here. With
``x_i = n_i/(n_w + \sum_j n_j)`` and ``m_i = n_i/(n_w M_w)``, the ratio
``x_i/(m_i M_w)`` equals ``x_w`` and tends to one at infinite dilution, where both
activity coefficients do, which yields

```math
\mu_i^\circ(\text{molality}) = \mu_i^\circ(\text{mole fraction}) + RT\ln(m^\circ M_w),
\qquad
\ln\gamma_i^{(x)} = \ln\gamma_i^{(m)} - \ln x_w .
```

The constant is not small, as its evaluation with the molar mass of water
computed from the atomic masses shows:

```@example stdstates
using ChemistryLab, DynamicQuantities
Mw = Species("H2O")[:M]
m° = 1.0u"mol/kg"
T = 298.15u"K"
uconvert(us"kJ/mol", R_GAS_Q * T * log(ustrip(m° * Mw)))
```

A standard potential transcribed from a mole-fraction database without this
correction shifts the solubility of every phase formed from that solute by
``\log_{10}(1/(m^\circ M_w)) \simeq 1.74`` per unit of stoichiometric
coefficient.

## 4. When the reference energy stops being a convention

A convention can be changed at no cost as long as the constant it introduces
cancels from every computed quantity. For a surface site with a fixed budget it
does: every surface reaction carries one site on each side, so that adding the
same constant to the ``\Delta_a G^\circ`` of all the members of a family changes
no reaction energy. The test suite checks it by shifting a whole family by as
much as 20 kJ/mol, after which the amount of the sorbent is unchanged to
``10^{-6}`` in relative terms, which is the tolerance of the solve itself.

The cancellation fails as soon as the budget follows the amount of the solid that
carries the sites, the case treated in
[Chemistry that happens on a surface](@ref sec-theory-surface) §10. The free
sites are then counted as part of the host, whose database formula already
contains the surface groups they are made of, so that a grain whose sites are all
free has the composition of that formula; it has the database energy only if the
free site carries no energy of its own. The reference is therefore fixed at
zero, following [Kulik2002](@citet), and it is no longer a convention: any other
value moves the host's solubility by ``\nu\,|\Delta_a G^\circ_{\text{free}}|/(RT
\ln 10)``, which [`host_coupling_bias`](@ref) evaluates. Giving the free site
`XsOH` the energy of the oxygen and hydrogen it holds, ``\mu^\circ(\mathrm{H_2O})
- \mu^\circ(\mathrm{H^+}) = -237.2`` kJ/mol, would count that energy twice, and
at the weak-site density of [DzombakMorel1990](@citet), ``\nu = 0.2``, it is worth
8.3 log units; a coupled family whose free site is more than 0.05 log units away
from zero is refused at construction.

## Where to go next

[Proving that an answer is the answer](@ref sec-theory-certificate) is the next
page of the chapter: it uses the potentials defined here to state what makes a
computed equilibrium provably the minimum. The two places where the activity
departs from its ideal form are treated in
[Activity models](@ref sec-theory-activity) and
[Solid solutions](@ref sec-theory-solid-solutions), and the surface conventions
in [Chemistry that happens on a surface](@ref sec-theory-surface).
