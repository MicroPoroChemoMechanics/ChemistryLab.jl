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
so that the standard state of a solute is at the pressure of the system. The
heat-capacity polynomial used for solids and gases carries no pressure term, and
their standard state is at ``P_r = 1`` bar whatever the pressure of the state.
This is exact at 1 bar. Above it, a condensed phase misses the contribution
``\int_{P_r}^{P} V_i^\circ\,\mathrm{d}P \simeq V_i^\circ (P - P_r)`` and an ideal
gas the contribution ``RT\ln(P/P_r)``. For portlandite, with ``V^\circ \simeq
33\ \mathrm{cm^3/mol}``, the first amounts to ``0.012\,RT`` at 10 bar and to
``1.3\,RT`` at 1 kbar: negligible for a laboratory sample or a structure, not
for a deep reservoir.

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
pure water. The dilute model takes ``a_w = x_w``, the mole fraction of water in
the aqueous phase. The other models obtain it from the osmotic coefficient
``\varphi``,

```math
\ln a_w = -M_w\,\varphi \sum_j m_j ,
```

which is the form imposed by the Gibbs-Duhem relation once the activity
coefficients of the solutes are given, as argued in
[Activity models](@ref sec-theory-activity) §3.

### Pure solids

A pure solid is its own standard state, its activity is one and its chemical
potential reduces to ``\mu_i^\circ``. It follows that a pure phase contributes to
the Gibbs energy a term linear in its amount, and that the minimization decides
whether it is present by comparing ``\mu_i^\circ`` with the combination of
component potentials its formula implies ([Thermochemistry](@ref sec-theory-thermo)
§4 and §5).

### Gases

A gas is referred to the pure ideal gas at ``P_r = 1`` bar, and its activity in
an ideal mixture is then ``a_i = x_i P/P_r``. The package retains the mole
fraction ``x_i`` alone, which coincides with it at 1 bar, and the remark on
pressure made in §1 applies.

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
[Chemistry that happens on a surface](@ref sec-theory-surface) §10. The host then
carries ``-\nu`` of the site component, the reference energy of the free site
enters its chemical potential, and that energy becomes a statement about matter.
A free site `XsOH` holds an oxygen and a hydrogen, and assigning it
``\Delta_a G^\circ = 0`` states that a surface hydroxyl forms from its elements at
no cost. The error is the energy of that matter, which [`host_coupling_bias`](@ref)
evaluates from the conservation matrix. At the weak-site density of
[DzombakMorel1990](@cite) for hydrous ferric oxide, where ``\nu = 0.2``, it
amounts to 8.3 log units on the solubility of the host, enough to dissolve it
completely where the same system with a fixed budget keeps its solid. The value
consistent with the rest of the database is the standard energy of the matter
the site carries, ``\mu^\circ(\mathrm{H_2O}) - \mu^\circ(\mathrm{H^+}) =
-237.2`` kJ/mol for a hydroxylated oxide, and a coupled family left at zero is
refused at construction with that value in the error message.
[Kulik2002](@cite) avoids the question altogether by keeping the free site out of
the balance, as a surface solvent of fixed activity; the surface page compares
the two routes.

## Where to go next

[Proving that an answer is the answer](@ref sec-theory-certificate) is the next
page of the chapter: it uses the potentials defined here to state what makes a
computed equilibrium provably the minimum. The two places where the activity
departs from its ideal form are treated in
[Activity models](@ref sec-theory-activity) and
[Solid solutions](@ref sec-theory-solid-solutions), and the surface conventions
in [Chemistry that happens on a surface](@ref sec-theory-surface).
