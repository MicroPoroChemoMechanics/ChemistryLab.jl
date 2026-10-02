# [Energy and entropy balances](@id sec-theory-energy-entropy)

[Energies, enthalpies and the chemical potential](@ref sec-theory-basics)
derives the first and second laws, the enthalpy a calorimeter measures and the
Gibbs energy the solver minimizes. This page does not repeat those
definitions. It writes the same two laws as **balances** on a chosen boundary,
which is where the conditions behind ``q_P=\Delta H`` and behind the minimum of
``G`` become explicit. It also explains why a reacting system can heat up while
remaining closed, and why an equilibrium calculation alone cannot predict a
calorimetric curve.

## 1. Choose the boundary before writing the balance

A system is the part of the experiment under consideration. Its boundary
determines which transfers belong in its balance.

| description | matter crosses the boundary | heat crosses the boundary | work crosses the boundary |
|:--|:--|:--|:--|
| closed | no | allowed | allowed |
| open | allowed | allowed | allowed |
| adiabatic | depends on the boundary | no | allowed |
| isolated | no | no | no |

These descriptions can be combined: a sealed, insulated piston-cylinder is
closed and adiabatic, but can do expansion work. A sealed, insulated, rigid
vessel with no other work exchange is isolated. A sealed sample can also be
isothermal if a thermostat removes the heat of reaction.

In a **closed reacting system**, the species amounts change even though no
matter crosses the boundary. For example,

```math
\mathrm{CaO + H_2O \longrightarrow Ca(OH)_2}
```

consumes lime and water and produces portlandite while conserving calcium,
hydrogen, oxygen and charge. In the package's notation the fixed quantity is
``\mathbf{A}\mathbf{n}=\mathbf{b}``, not each ``n_i`` or the sum of the species
amounts.

This is classical chemical thermodynamics: nuclear transformations and the
relativistic contribution of exchanged energy to a sample's mass are outside
the model. Closure concerns material transfer; it does not require the energy
of the sample to remain constant. In a relativistic description, exchanged
energy contributes to the system's mass; classical element budgets and molar
masses neglect that correction while retaining the chemical energy balance.

## 2. What the first law needs before ``q_P=\Delta H``

The heat ``q`` and the work ``w`` are counted positive when received by the
system, as in [Energies, enthalpies and the chemical potential](@ref sec-theory-basics) §1.
Neglecting changes in bulk kinetic and gravitational potential energy, the
first law for a closed system reads, in differential form,

```math
\mathrm{d}U = \delta q - P_{\mathrm{ext}}\,\mathrm{d}V + \delta w_{\mathrm{other}}.
```

Heat and work describe transfers along a process: neither is a stored property
of the sample. This is why their infinitesimal transfers use ``\delta`` whereas
a state function uses ``\mathrm{d}``. See
[MIT's discussion of heat and work](https://www.ocw.mit.edu/ans7870/16/16.unified/thermoF03/chapter_3.htm).

The pressure in expansion work is the **external** pressure. Replacing it by
the system's pressure requires mechanical equilibrium. Constant pressure alone
does not rule out electrical work, stirring or other forms of work. With
``H=U+PV``, a mechanically equilibrated closed system obeys

```math
\mathrm{d}H
= \mathrm{d}U + P\,\mathrm{d}V + V\,\mathrm{d}P
= \delta q + \delta w_{\mathrm{other}} + V\,\mathrm{d}P.
```

At constant pressure, with only pressure-volume work, this gives
``\Delta H=q_P``. With other work, the balance is instead
``\Delta H=q_P+w_{\mathrm{other}}``. A rigid closed calorimeter with no other
work measures ``\Delta U=q_V``; conversion to ``\Delta H`` requires the change
in ``PV`` and any corrections needed to refer the measurement to the stated
temperature and standard states.

An enthalpy change is an energy difference, in joules. A heat rate is an energy
transfer per unit time, in watts. Enthalpy changes can include reaction,
heating, phase change and mixing; they are not restricted to bond energies.

## [3. Why an adiabatic reaction can raise the temperature](@id sec-theory-adiabatic)

For a closed adiabatic system at constant pressure, with only pressure-volume
work, ``H`` remains constant. For a rigid closed adiabatic system with no work,
it is ``U`` that remains constant. Neither condition requires constant
temperature.

Consider the constant-pressure case. Let the initial state have composition
``\mathbf{n}_0`` at ``T_0`` and the final state composition ``\mathbf{n}_f`` at
``T_f``. Evaluating ``H(T_f,P,\mathbf{n}_f)-H(T_0,P,\mathbf{n}_0)=0`` through an
imaginary path, first changing the composition at ``T_0`` and then heating the
final composition, gives the balance already stated in
[Energies, enthalpies and the chemical potential](@ref sec-theory-basics) §1,

```math
\underbrace{H(T_0,P,\mathbf{n}_f)-H(T_0,P,\mathbf{n}_0)}_{\Delta H_{\mathrm{composition}}(T_0)}
+ \int_{T_0}^{T_f} C_p(T,P,\mathbf{n}_f)\,\mathrm{d}T = 0.
```

The first term is a **difference between two compositions**, not the enthalpy
of either state alone. Both terms are total energies: ``C_p`` here is the heat
capacity of the final sample, in J/K. If a molar reaction enthalpy
``\Delta_r H`` is used, it must be integrated over the reaction extent, in moles.
Writing ``\Delta_r H+\int C_p\,\mathrm{d}T=0`` without a common amount basis
would mix J/mol and J.

An exothermic composition change makes the first term negative, and heating
can compensate it. The integral does not presume that ``T_f`` is known: this
equation determines it. If the final composition itself depends on ``T_f``,
composition and temperature must be solved together. Phase transitions along
the imaginary heating path require their enthalpy jumps as well. A vessel's
heat capacity belongs in this balance if the vessel is part of the system.

## 4. Entropy exchange and entropy production

The inequality of Clausius in
[Energies, enthalpies and the chemical potential](@ref sec-theory-basics) §2
becomes a balance once the entropy created inside the system is named. For a
closed system exchanging heat at a boundary temperature ``T_b``,

```math
\mathrm{d}S = \frac{\delta q}{T_b} + \mathrm{d}S_{\mathrm{gen}},
\qquad \mathrm{d}S_{\mathrm{gen}} \ge 0.
```

With several heat contacts, sum ``\delta q_k/T_{b,k}``. The denominator belongs
to the heat contact in the entropy balance; it cannot in general be replaced
by a single temperature of a sample containing thermal gradients. An
adiabatic closed system has ``\Delta S=S_{\mathrm{gen}}\ge0`` even while its
temperature changes. A closed system that exports enough entropy as heat can
have ``\Delta S<0``; the nonnegative quantity is the entropy production.

To evaluate the difference between two thermodynamic states, one may use any
**reversible reference path** connecting them:

```math
\Delta S = \int_{\mathrm{initial}}^{\mathrm{final}}
              \frac{\delta q_{\mathrm{rev}}}{T}.
```

This is a path integral of heat transfer divided by the temperature at which
each transfer occurs. For reversible heating at fixed pressure and composition,
``\delta q_{\mathrm{rev}}=C_p\,\mathrm{d}T``, so it becomes
``\int_{T_0}^{T_f} C_p/T\,\mathrm{d}T``. For an isothermal reversible path it
becomes ``q_{\mathrm{rev}}/T``. These are different parameterizations of the
same definition, not a general replacement of heat by temperature.

Entropy describes the thermodynamic availability of energy and the possible
microscopic states. Reading it as disorder is a useful picture, but the picture
alone does not establish a balance or the sign of a reaction entropy. The
distinction between conservation and direction is developed in
[AndersonCrerar1993](@cite), Chapter 5.

## [5. What the Gibbs energy tells us](@id sec-theory-gibbs-criterion)

``G=H-TS`` is defined at each thermodynamic state, and its existence does not
require the temperature to remain constant along a process. The **minimum
criterion** derived in
[Energies, enthalpies and the chemical potential](@ref sec-theory-basics) §2
does require specified constraints. For a closed system whose initial and final
states share the temperature ``T`` of one heat reservoir and the same imposed
pressure, with only pressure-volume work, the entropy balance of §4 applied to
the system and the reservoir reads

```math
\Delta G=\Delta H-T\Delta S_{\mathrm{system}}
=-T\Delta S_{\mathrm{universe}}\le0,
\qquad q_P=\Delta H .
```

Equilibrium minimizes ``G`` over the allowed compositions at fixed ``T``,
``P`` and conservation budget. Other boundaries give other criteria: an
isolated rigid system maximizes ``S`` at fixed ``U``, ``V`` and budget, and an
adiabatic constant-pressure system with only pressure-volume work maximizes
``S`` at fixed ``H``, ``P`` and budget.

``G`` at equilibrium need not be zero. Its first variation vanishes for allowed
reaction directions that can proceed both ways; at a boundary, such as an
absent pure phase, only the feasible direction is tested and an inequality
replaces that equality. A finite difference between an initial state and an
equilibrium state need not vanish either.

### Why ``\Delta H`` and ``T\Delta S`` are not two heat sources

At fixed ``T`` and ``P``, a reversible process can deliver useful work in
addition to expansion work. Then ``q_{\mathrm{rev}}=T\Delta S`` and the first
law gives

```math
w_{\mathrm{other,rev}}=\Delta H-q_{\mathrm{rev}}=\Delta G.
```

With the received-work convention, the maximum useful work **delivered** is
``-\Delta G``. In a process with no useful work extraction, ``q_P=\Delta H``
instead, and irreversibility accounts for the difference from
``T\Delta S``. These statements concern different processes between the same
states. They do not split the actual heat into a "chemical" contribution
``\Delta H`` and a separate "thermal" contribution ``T\Delta S``.
The conditions for the minimum criterion and for useful work extraction are
also derived in [Oxford's chemical thermodynamics notes](https://manolopoulos.chem.ox.ac.uk/downloads/thermo2.pdf),
Section A.

## 6. State functions and kinetic paths

A state includes temperature, pressure, composition and phase information.
Two histories reaching exactly the same initial and final states have the same
``\Delta H`` and ``\Delta G``. If one history's heating triggers a different
reaction and a different final composition, its final state is different and
so can be its state-function differences.

At an observation time, a hydrating paste retains slowly reacting clinker and
may contain phases whose formation or conversion is inhibited. The coupling
in [Coupling kinetics and equilibrium](@ref sec-coupling) integrates the slow
amounts and re-equilibrates the remaining budget. It computes successive
**partial equilibria**, not thermal jumps over molecular activation barriers.
Rate laws supply the time dependence, and a thermal balance supplies a changing
temperature when that is modeled.

### [What the package counts as heat](@id sec-theory-heat-output)

For a closed isothermal constant-pressure sample with only pressure-volume
work, the heat released up to time ``t``, counted positive as in
[Energies, enthalpies and the chemical potential](@ref sec-theory-basics) §1, is

```math
Q(t)=H(T,P,\mathbf{n}_0)-H(T,P,\mathbf{n}(t)),
\qquad
\dot Q=-\frac{\mathrm{d}H}{\mathrm{d}t}.
```

The rate depends on the trajectory; the cumulative value depends on the two
states actually reached. Changing temperature, matter transfer or other work
requires the corresponding balance terms. In particular, an adiabatic
temperature rise is not an outward heat transfer.

The package evaluates ``H=\sum_i n_i\,\Delta_a H_i^\circ`` with
[`enthalpy`](@ref): the standard enthalpies of the species, without the excess
enthalpies of the aqueous or solid-solution activity models, as
[Energies, enthalpies and the chemical potential](@ref sec-theory-basics) §3
explains. Using a non-ideal activity model for equilibrium therefore does not
by itself add its heat of mixing to the reported calorimetry. For a
concentration-based activity, the temperature derivative behind an excess
enthalpy also includes the thermal expansion of the solution, which changes
molarity but not molality.

[`heat_release`](@ref) differences this enthalpy between two certified
speciations on the same conserved budget, so it counts the precipitation of
hydrates as well as the dissolution of the clinker. In a kinetic run, the heat
of [`cumulative_heat`](@ref) depends on the formulation. When the kinetic
reactions themselves produce the hydrates, it is [`heat_rate`](@ref), summed
over those reactions. Under **partial equilibrium** the kinetic reactions only
release ions, so the run instead follows ``-\mathrm{d}H/\mathrm{d}t`` over the
whole composition, the equilibrium part included, and refuses a system in
which a species has no enthalpy of formation; [`missing_enthalpy`](@ref) lists
such species. Its integral then matches [`heat_release`](@ref) to within the
accuracy of the in-run partition.

Waiting longer does not by itself guarantee a unique observed assemblage.
Accessibility of transformations, boundary conditions, the species list and
possible degeneracy of equilibrium all matter. The distinction between
physical metastability, a constrained equilibrium and numerical stationarity
is developed in [Proving that an answer is the answer](@ref sec-theory-certificate).

## Where to go next

[Thermochemistry](@ref sec-theory-thermo) takes up these quantities in the
notation of the code, with the potential the solver minimizes and the models by
which ``\Delta_a G^\circ(T,P)`` is evaluated.
