# [Energies, enthalpies and the chemical potential](@id sec-theory-basics)

The pages of this chapter use the chemical potential, the standard Gibbs energy
of a species and its enthalpy of formation as familiar objects, although the
relations between them are rarely obvious on a first reading, since a
calorimeter measures a change of enthalpy whereas the solver minimizes a Gibbs
energy assembled from chemical potentials. This page rebuilds them from the two
laws of thermodynamics, mostly in the order in which [AndersonCrerar1993](@cite)
introduce them, up to the relation between the chemical potential of a species
and the apparent Gibbs energy a database stores, which is where
[Thermochemistry](@ref sec-theory-thermo) begins.

## 1. Heat, work, and why a calorimeter measures an enthalpy

A system is the part of the world under consideration, closed when it exchanges
energy but no matter with its surroundings, and its state is fixed by a few
variables, the temperature ``T``, the pressure ``P`` and the amounts ``n_i`` of
its constituents. A state function is a quantity whose value depends on the
state only, so that its change between two states does not depend on the way
followed from one to the other. The internal energy ``U`` is such a function,
whereas the heat ``q`` and the work ``w`` received by the system are not: the
same change of state can be brought about with more work and less heat, or the
reverse. The first law states that their sum is nevertheless fixed by the two
end states [AndersonCrerar1993](@cite) (§4.6),

```math
\Delta U = q + w .
```

In a chemical system held at constant pressure, the only work exchanged is
usually that of the pressure against the change of volume, ``w = -P\,\Delta V``,
and the heat received is then ``q_P = \Delta U + P\,\Delta V``. It is therefore
natural to introduce the enthalpy

```math
H = U + PV ,
```

another state function, whose change at constant pressure is exactly the heat
received, ``q_P = \Delta H``. A reaction that gives off heat at constant pressure
is one whose enthalpy decreases, and the heat recorded by an isothermal
calorimeter, counted positive when released, is ``Q = -\Delta H``. The
difference between ``\Delta H`` and ``\Delta U`` is the work ``P\,\Delta V``,
small for condensed phases and considerable as soon as a gas is produced or
consumed.

The heat capacity at constant pressure measures how the enthalpy grows with
temperature at fixed composition,

```math
C_p = \left(\frac{\partial H}{\partial T}\right)_{P,n} ,
```

and it governs the other kind of calorimeter. In a vessel that exchanges no
heat, the enthalpy of the contents is constant, so that the enthalpy released by
the reaction at the initial temperature is found again as the heating of the
contents after the reaction, ``\Delta H(T_0) + \int_{T_0}^{T} C_p\,\mathrm{d}T' = 0``. A
semi-adiabatic vessel adds the heat it loses to this balance, and it is the
subject of [A CEM I 52.5 N mortar in a semi-adiabatic calorimeter, inside the kinetics](@ref ex-semiadiabatic).

## 2. Entropy, and the potential of a system held at fixed temperature and pressure

The first law says nothing about the direction of a change. The second law
supplies it through a further state function, the entropy ``S``, which obeys, for
any transformation of a closed system at temperature ``T``, the inequality of
Clausius [AndersonCrerar1993](@cite) (§5.2, §5.8)

```math
\mathrm{d}S \;\ge\; \frac{\delta q}{T} ,
```

the equality holding for a reversible transformation only. Unlike ``U`` and
``H``, the entropy can be given an absolute value, the third law fixing its zero;
§5 defines it, together with the entropy of formation it must not be confused
with.

At constant temperature and pressure, with no work other than that of the
pressure, ``\delta q = \mathrm{d}H`` and the inequality becomes
``\mathrm{d}H - T\,\mathrm{d}S \le 0``. It is thus the function

```math
G = H - TS ,
```

the Gibbs energy, that can only decrease in a spontaneous transformation at fixed
``T`` and ``P``, and reaches its minimum at equilibrium. This is the principle
the solver of this package implements. The two terms of ``G`` weigh two
tendencies against each other: a change releasing heat lowers ``H``, a change
creating disorder raises ``S``, and a reaction that absorbs heat can still
proceed when the entropy it creates outweighs it, as the dissolution of many
salts does. The Gibbs energy is also the Legendre transform of the internal
energy whose natural variables are ``T`` and ``P`` (§5.4), so that for a closed
system of fixed composition

```math
\mathrm{d}G = -S\,\mathrm{d}T + V\,\mathrm{d}P ,
\qquad
\left(\frac{\partial G}{\partial T}\right)_P = -S .
```

Combining the second relation with ``G = H - TS`` eliminates the entropy and
yields the Gibbs-Helmholtz relation,

```math
\left(\frac{\partial (G/T)}{\partial T}\right)_P = -\frac{H}{T^2} ,
```

by which the enthalpy, the quantity calorimetry measures, is read on the
temperature dependence of the Gibbs energy, the quantity equilibrium depends on.

## 3. Open systems: the chemical potential

When matter is exchanged or transformed, ``G`` depends on the amounts as well, and
its differential gains one term per constituent,

```math
\mathrm{d}G = -S\,\mathrm{d}T + V\,\mathrm{d}P + \sum_i \mu_i\,\mathrm{d}n_i ,
\qquad
\mu_i = \left(\frac{\partial G}{\partial n_i}\right)_{T,P,n_{j\neq i}} .
```

The chemical potential ``\mu_i`` is thus the partial molar Gibbs energy of
constituent ``i``: the change of ``G`` per mole of ``i`` added to a system so
large that its composition does not change [AndersonCrerar1993](@cite) (§9.2).
Every extensive quantity has partial molar counterparts defined in the same
way, among which the partial molar enthalpy ``h_i = (\partial H/\partial
n_i)_{T,P,n_{j\neq i}}`` and entropy ``s_i``. Differentiating ``G = H - TS``
with respect to ``n_i`` gives ``\mu_i = h_i - T s_i``, and differentiating the
Gibbs-Helmholtz relation of §2 in the same way gives

```math
\left(\frac{\partial (\mu_i/T)}{\partial T}\right)_{P,n} = -\frac{h_i}{T^2} .
```

This is the link between the two families of quantities: the chemical potential
carries the Gibbs energy of a constituent, its temperature dependence carries
the enthalpy.

At fixed temperature and pressure, multiplying every amount by ``\lambda``
multiplies ``G`` and ``H`` by ``\lambda``. Euler's theorem on functions
homogeneous of degree one then expresses each of them as the sum of the amounts
weighted by the partial molar quantities (§9.2),

```math
G = \sum_i n_i\,\mu_i ,
\qquad
H = \sum_i n_i\,h_i ,
```

so that the enthalpy of a system, and hence the heat of a change at constant
pressure, is obtained by weighting the partial molar enthalpies by the amounts.
Differentiating the first expression and subtracting ``\mathrm{d}G = \sum_i
\mu_i\,\mathrm{d}n_i``, valid at fixed ``T`` and ``P``, leaves the Gibbs-Duhem
relation, which holds within each phase,

```math
\sum_i n_i\,\mathrm{d}\mu_i = 0 \qquad (T,\ P \text{ fixed}) .
```

The chemical potentials of the constituents of one phase are thus bound
together, and a model cannot prescribe the activity of the solvent
independently of the activity coefficients it assigns to the solutes; it is by
this relation that
[Activity models](@ref sec-theory-activity) §3 obtains the activity of water.

The chemical potential of a constituent of a mixture is written as a standard
part, the potential ``\mu_i^\circ(T,P)`` of a reference state, plus a
contribution of the composition through the activity ``a_i``,

```math
\mu_i = \mu_i^\circ(T,P) + RT\ln a_i ,
```

the reference state being the subject of [Standard states](@ref sec-theory-standard-states).
Inserted in the relation above, the logarithm of an ideal activity does not
depend on temperature at fixed composition and drops out, so that the partial
molar enthalpy of a constituent of an ideal mixture is its standard enthalpy
``h_i^\circ``; a non-ideal mixture adds an excess enthalpy, carried by the
temperature dependence of the activity coefficients. The enthalpy of a state
evaluated by [`enthalpy`](@ref), and with it the heat of the calorimeters, is
``\sum_i n_i h_i^\circ``, which leaves the excess enthalpies out.

At equilibrium, the minimum of ``G`` under the conservation of the elements
implies two conditions, independent of the path by which the system reached
it (§14.2): a constituent present in two phases has the same chemical potential in
both, and for any reaction ``\sum_i \nu_i\,\mathrm{A}_i = 0`` among the constituents,
counted positive for the products,

```math
\Delta_r G \;\equiv\; \sum_i \nu_i\,\mu_i = 0 .
```

Away from equilibrium, ``\mathcal{A} = -\Delta_r G`` is the affinity of the
reaction, positive when it proceeds forward (§14.5), which is the quantity the
rate laws of [Rate laws](@ref sec-theory-kinetics) are written with.

## [4. Reactions and equilibrium constants](@id sec-theory-reactions)

Inserting ``\mu_i = \mu_i^\circ + RT\ln a_i`` into ``\Delta_r G`` separates the
part of the standard states from that of the composition,

```math
\Delta_r G = \Delta_r G^\circ + RT\ln Q_r ,
\qquad
\Delta_r G^\circ = \sum_i \nu_i\,\mu_i^\circ ,
\qquad
Q_r = \prod_i a_i^{\nu_i} ,
```

and the condition of equilibrium ``\Delta_r G = 0`` fixes the value ``Q_r`` takes
there, the equilibrium constant [AndersonCrerar1993](@cite) (§13.1),

```math
\ln K = -\frac{\Delta_r G^\circ}{RT} .
```

The standard Gibbs energy of reaction decides which side a reaction favors
between reactants and products in their standard states, and the equilibrium
constant expresses the same thing in terms of activities. Written as above,
``\Delta_r G^\circ`` requires the standard potentials of the species, which no
measurement gives; §5 shows that the energies of formation a database tabulates
are sufficient, and §6 how they depend on temperature and pressure.

The package itself writes no reaction to find an equilibrium: it minimizes
``G`` over all the constituents at once, and the conditions ``\Delta_r G = 0`` of
every reaction that can be written among them follow from the optimality
conditions of the minimization ([Proving that an answer is the answer](@ref sec-theory-certificate)).

## 5. Energies known up to a constant, and the reactions of formation

The internal energy, the enthalpy and the Gibbs energy are known up to an
additive constant only, since a measurement gives their changes and never their
values. A database cannot tabulate the change of every reaction of interest,
which are too many, and tabulates instead, for each substance, the change of the
one reaction that forms it from its elements, each taken in its stable form at
the conditions of reference, ``T_r = 298.15`` K and ``P_r = 1`` bar
[AndersonCrerar1993](@cite) (§7.1, §7.2): the standard enthalpy of formation
``\Delta_f H^\circ`` and the standard Gibbs energy of formation
``\Delta_f G^\circ``. That these values determine ``\Delta_r G^\circ`` follows
from the conservation of the elements and of the charge in a reaction, and the
argument below also delimits the freedom left in their definition, of which §6
makes use.

### Balanced reactions

Let species ``\mathrm{A}_i`` contain ``\alpha_{ei}`` atoms of element ``e`` and carry the
charge ``z_i``. The reaction ``\sum_i \nu_i\,\mathrm{A}_i = 0`` is balanced when it
conserves every element and the charge,

```math
\sum_i \nu_i\,\alpha_{ei} = 0 \quad\text{for every element } e ,
\qquad
\sum_i \nu_i\,z_i = 0 .
```

The vector ``\boldsymbol{\nu}`` of the stoichiometric coefficients thus belongs to
the null space of the matrix whose rows hold the ``\alpha_{ei}`` of each element
and the charges ``z_i``, which the package builds as the canonical stoichiometric
matrix, the charge being its row `Zz` ([Stoichiometric Matrix](@ref sec-stoich-matrices)).
The slaking of lime, CaO + H₂O → Ca(OH)₂, has ``\boldsymbol{\nu} = (-1, -1, 1)``
and conserves one calcium, two hydrogens and two oxygens.

### Formation from the elements

The reaction of formation of ``\mathrm{A}_i`` takes ``\alpha_{ei}`` atoms of each element in
its reference form, half a mole of O₂ gas for each atom of oxygen for instance.
An ion cannot be formed that way without a counter-charge, and the convention is
to complete its formation with ``z_i`` hydrogen ions turned into hydrogen gas
[AndersonCrerar1993](@cite) (§17.3), as in Ca + 2 H⁺ → Ca²⁺ + H₂ for the calcium
ion. In both cases, the standard Gibbs energy of formation reads

```math
\Delta_f G_i^\circ = \mu_i^\circ - \sum_e \alpha_{ei}\,G_e^\circ - z_i\,G_Z^\circ ,
\qquad
G_Z^\circ = \mu_{\mathrm{H^+}}^\circ - G_{\mathrm{H}}^\circ ,
```

where ``G_e^\circ`` is the standard molar Gibbs energy of element ``e`` in its
reference form, counted per atom, and ``G_Z^\circ`` the Gibbs energy of formation
of the hydrogen ion from hydrogen gas, ``G_{\mathrm{H}}^\circ`` being half the
potential of H₂. Neither is known, and neither needs to be.
The definition makes the energy of formation of an element in its reference form
zero, and that of H⁺ zero as well, which is the convention on ions: the energies
of formation of all other ions follow from it. Inserted in ``\Delta_r G^\circ =
\sum_i \nu_i\,\mu_i^\circ``, the energies of formation give

```math
\Delta_r G^\circ = \sum_i \nu_i\,\Delta_f G_i^\circ
  + \sum_e G_e^\circ \sum_i \nu_i\,\alpha_{ei}
  + G_Z^\circ \sum_i \nu_i\,z_i ,
```

in which the last two sums vanish for a balanced reaction, so that

```math
\Delta_r G^\circ = \sum_i \nu_i\,\Delta_f G_i^\circ .
```

The proof requires nothing of ``G_e^\circ`` and ``G_Z^\circ`` besides being the
same for all the species of the reaction: any values assigned to them on that
condition leave ``\Delta_r G^\circ`` unchanged. The same argument applied to the
enthalpy yields ``\Delta_r H^\circ = \sum_i \nu_i\,\Delta_f H_i^\circ``, which is
Hess's law: the change in a reaction is the same whether the reaction is carried
out directly or by decomposing the reactants into their elements and then
forming the products from them.

### [Absolute entropy and entropy of formation](@id sec-theory-absolute-entropy)

The entropy escapes the indeterminacy of the energies, because the third law
fixes its zero: the entropy of a pure substance in a perfect crystalline form
tends to zero as the temperature tends to 0 K [AndersonCrerar1993](@cite)
(§6.5). The entropy at any other temperature then follows from the second law,
``\mathrm{d}S = \delta q_{\rm rev}/T = C_p\,\mathrm{d}T/T`` at constant pressure,
integrated from 0 K, with the entropy of each phase transition met on the way,
the ratio of its enthalpy to the temperature at which it occurs,

```math
S^\circ(T) = \int_0^{T} \frac{C_p^\circ(T')}{T'}\,\mathrm{d}T'
  + \sum_k \frac{\Delta_{\rm trs} H_k}{T_k} .
```

This is the absolute, or third-law, entropy, tabulated as ``S^\circ`` at
``T_r``: the heat capacity is measured by calorimetry from a few kelvins upward
and extrapolated to 0 K (§7.3, §7.5), and ``S^\circ`` is a positive quantity. A
substance that is not a perfect crystal at 0 K, a glass or a crystal whose
disorder freezes in on cooling, retains there a residual entropy, which adds to
the integral (§6.5.3); the glass of a blast-furnace slag is of this kind.

The entropy of formation is of another nature: it is the change of entropy in
the reaction of formation,

```math
\Delta_f S_i^\circ = S_i^\circ - \sum_e \alpha_{ei}\, S_e^\circ ,
```

where ``S_e^\circ`` is the absolute entropy of element ``e`` in its reference
form, counted per atom, half that of O₂ gas for oxygen. It is the entropy of
formation that is bound to the enthalpy and the Gibbs energy of formation by the
relation of §2 applied to the reaction of formation,

```math
\Delta_f G_i^\circ = \Delta_f H_i^\circ - T\,\Delta_f S_i^\circ ,
```

whereas the entropy of a reaction is computed from the absolute entropies,
``\Delta_r S^\circ = \sum_i \nu_i S_i^\circ``, the entropies of the elements
canceling for the reason given above. The two are easily confused, and they
differ even in sign: the absolute entropy of lime is positive, its entropy of
formation from calcium metal and oxygen gas negative, since the gas it consumes
is far more disordered than the solid it forms.

The aqueous ions carry a third kind of entropy. The convention that sets the
formation properties of H⁺ to zero sets its tabulated entropy to zero as well, and
the entropy tabulated for an ion of charge ``z_i`` is then a conventional
entropy, ``S_i^\circ - z_i\,S_{\mathrm{H^+}}^\circ``, its absolute entropy less
``z_i`` times that of the hydrogen ion, which thermodynamics alone does not
determine (§17.3). A conventional entropy can be negative, as that of Ca²⁺ is,
and the unknown term cancels from any reaction, whose charges balance,
``\sum_i \nu_i z_i = 0``. Heat capacities and volumes of ions are conventional in
the same way.

The values of CEMDATA18 at ``T_r`` show each of these quantities, the element
entropies being read from the same file; the identity between the three
formation quantities of lime holds to the internal consistency of the tabulated
values:

```@example basics
using ChemistryLab, Printf
db = Dict(symbol(s) => s for s in build_species(datapath("cemdata18-thermofun.json"); verbose = false))
S_e = Dict(e["symbol"] => only(e["entropy"]["values"])
           for e in ChemistryLab.JSON.parsefile(datapath("cemdata18-thermofun.json"))["elements"])
Tr = 298.15
G(s, T) = db[s][:ΔₐG⁰](T = T, unit = false)
H(s, T) = db[s][:ΔₐH⁰](T = T, unit = false)
S(s, T) = db[s][:S⁰](T = T, unit = false)
for s in ("O2", "H2", "H+", "Ca+2")
    @printf("%-4s  ΔfG = %9.0f, ΔfH = %9.0f J/mol, S = %6.2f J/(mol K)\n", s, G(s, Tr), H(s, Tr), S(s, Tr))
end
@printf("oxygen per atom: %.3f J/(mol K), half that of O2\n", S_e["O"])
ΔfS = S("Lim", Tr) - sum(k * S_e[string(e)] for (e, k) in atoms(db["Lim"]))
@printf("lime: S = %.2f, ΔfS = %.2f J/(mol K); ΔfG - (ΔfH - Tr ΔfS) = %.0f J/mol\n",
        S("Lim", Tr), ΔfS, G("Lim", Tr) - (H("Lim", Tr) - Tr * ΔfS))
```

With the entropy of reaction defined above, the three quantities of a reaction
are bound by ``\Delta_r G^\circ = \Delta_r H^\circ - T\Delta_r S^\circ``, which the
slaking of lime verifies with the same values:

```@example basics
slaking = ("Portlandite" => 1, "Lim" => -1, "H2O@" => -1)
Δr(f, T) = sum(ν * f(s, T) for (s, ν) in slaking)
@printf("slaking: ΔrH = %.1f kJ/mol, ΔrG = %.1f kJ/mol, ΔrS = %.2f J/(mol K); ΔrG - (ΔrH - Tr ΔrS) = %.0f J/mol\n",
        Δr(H, Tr) / 1000, Δr(G, Tr) / 1000, Δr(S, Tr), Δr(G, Tr) - (Δr(H, Tr) - Tr * Δr(S, Tr)))
```

### [Formation from primary species](@id sec-theory-primaries)

The elements are not the only possible reference. Let ``B_1, \dots, B_C`` be
species of the system, called primary, whose columns of ``\alpha_{ec}`` and
``z_c`` are linearly independent and generate those of all the other species.
Each species then decomposes in a unique way on the primaries,

```math
\alpha_{ei} = \sum_c A_{ci}\,\alpha_{ec} \quad\text{for every element } e ,
\qquad
z_i = \sum_c A_{ci}\,z_c ,
```

which is the balanced reaction ``\sum_c A_{ci}\,\mathrm{B}_c \rightarrow \mathrm{A}_i`` forming
``\mathrm{A}_i`` from the primaries, with ``A_{ci} = 1`` if ``\mathrm{A}_i`` is the primary
``\mathrm{B}_c`` and zero otherwise. The matrix of the ``A_{ci}`` is the
conservation matrix ``\mathbf{A}``, `SM.A` in the package, whose rows are labeled by the
primaries. Substituting this decomposition in the balance of the elements and of
the charge shows that a reaction conserves them if and only if it conserves the
primaries,

```math
\sum_i \nu_i\,A_{ci} = 0 \quad\text{for every primary } c ,
```

the converse resting on the independence of the columns of the primaries. The
reaction forming ``\mathrm{A}_i`` from the primaries being balanced, its standard Gibbs
energy follows from the Gibbs energies of formation by the result above, and it
defines
the equilibrium constant ``K_i`` of that reaction,

```math
\Delta_r G_i^\circ = \mu_i^\circ - \sum_c A_{ci}\,\mu_c^\circ
  = \Delta_f G_i^\circ - \sum_c A_{ci}\,\Delta_f G_c^\circ
  = -RT\ln K_i .
```

The argument used for the elements applies with the primaries in their place,
the unknown ``\mu_c^\circ`` taking the part of ``G_e^\circ`` and ``G_Z^\circ``,
and any balanced reaction satisfies

```math
\Delta_r G^\circ = \sum_i \nu_i\,\Delta_r G_i^\circ = -RT\sum_i \nu_i \ln K_i .
```

A table of formation constants ``K_i`` from a basis of primary species, which is
the form of the data of mass-action programs
([Mass action or minimization](@ref sec-theory-mass-action)), thus carries the
same information as a table of Gibbs energies of formation, for the species it
contains. For lime and portlandite, with Ca²⁺, H₂O and H⁺ as
primaries, the matrix and the two routes to the Gibbs energy of slaking read:

```@example basics
basis = [db["Ca+2"], db["H2O@"], db["H+"]]
SM = StoichMatrix([db["Lim"], db["Portlandite"], basis...], basis)
names = symbol.(SM.species)
ν = [get(Dict(slaking), s, 0) for s in names]     # slaking, in the order of the columns
# Gibbs energy of formation of each species from the primaries (rows of SM.A)
ΔrG_B = [G(s, Tr) - sum(SM.A[c, j] * G(symbol(b), Tr) for (c, b) in enumerate(SM.primaries))
         for (j, s) in enumerate(names)]
logK(s) = -ΔrG_B[findfirst(==(s), names)] / (R_GAS * Tr * log(10))
@printf("log K of formation from the primaries: lime %.2f, portlandite %.2f\n", logK("Lim"), logK("Portlandite"))
@printf("primaries conserved by slaking: %s; ΔrG = %.3f kJ/mol from the primaries, %.3f from the elements\n",
        SM.A * ν, sum(ν .* ΔrG_B) / 1000, Δr(G, Tr) / 1000)
SM
```

The oxides of a silicate, a carbonate or a sulfate form such a basis, and the
enthalpy of formation from the oxides is in use for these compounds (§7.4.4).
Its values are much smaller in magnitude, and so are their uncertainties, which
no longer include those of forming the oxides from the elements, but they cannot
be combined with values of formation from the elements in one calculation. It is
the natural scale for a glass, and
[CEM I 52.5 N and slag in an isothermal calorimeter, read off the states](@ref sec-example-isothermal)
places the slag of a blended cement on it.

## 6. Temperature, pressure, and the apparent quantities of a database

### The standard potential at any temperature and pressure

The standard state of a species being taken at the temperature and the pressure
of the system, its potential obeys the relation of §2 for a substance of fixed
composition, with the standard molar entropy and volume,

```math
\mathrm{d}\mu_i^\circ = -S_i^\circ\,\mathrm{d}T + V_i^\circ\,\mathrm{d}P ,
```

whereas a gas, whose standard state stays at ``P_r`` whatever the pressure,
keeps the first term only and receives the pressure through its activity
([Standard states](@ref sec-theory-standard-states)). The relation is integrated
from ``(T_r, P_r)`` in temperature at ``P_r``, where
``S_i^\circ(T') = S_i^\circ(T_r) + \int_{T_r}^{T'} C_{P,i}^\circ\,\mathrm{d}T''/T''``,
and then in pressure at ``T``. Exchanging the order of the two temperature
integrals gives [AndersonCrerar1993](@cite) (§7.4, §7.6)

```math
\mu_i^\circ(T,P) - \mu_i^\circ(T_r,P_r) =
  - S_i^\circ(T_r,P_r)\,(T - T_r)
  + \int_{T_r}^{T} C_{P,i}^\circ\,\mathrm{d}T'
  - T\int_{T_r}^{T} \frac{C_{P,i}^\circ}{T'}\,\mathrm{d}T'
  + \int_{P_r}^{P} V_i^\circ(T,P')\,\mathrm{d}P' ,
```

the heat capacity being taken at ``P_r``. The enthalpy follows in the same way
from ``\mathrm{d}H = C_p\,\mathrm{d}T + [V - T(\partial V/\partial T)_P]\,\mathrm{d}P``,

```math
H_i^\circ(T,P) - H_i^\circ(T_r,P_r) =
  \int_{T_r}^{T} C_{P,i}^\circ\,\mathrm{d}T'
  + \int_{P_r}^{P} \left[V_i^\circ - T\left(\frac{\partial V_i^\circ}{\partial T}\right)_{P}\right]\mathrm{d}P' .
```

The entropy that multiplies ``T - T_r`` is the absolute entropy of §5, and
both expressions involve properties of the species alone.

### [Apparent and formation Gibbs energies](@id sec-theory-apparent)

Away from ``(T_r, P_r)``, an energy of formation depends on the conditions at
which the elements are taken, and two choices are in use
[AndersonCrerar1993](@cite) (§7.4.1, §7.4.2). The traditional quantities of
formation ``\Delta_f G^\circ(T,P)`` take the elements at the temperature ``T``
and at 1 bar, which requires their heat capacities, and their phase changes,
over the whole range. The apparent quantities of formation leave them at
``(T_r, P_r)``: in the definition of §5, ``G_e^\circ`` becomes the constant
``G_e^\circ(T_r,P_r)`` and ``G_Z^\circ`` becomes
``\mu_{\mathrm{H^+}}^\circ(T,P) - G_{\mathrm{H}}^\circ(T_r,P_r)``,

```math
\Delta_a G_i^\circ(T,P) = \mu_i^\circ(T,P) - \sum_e \alpha_{ei}\,G_e^\circ(T_r,P_r)
  - z_i\left[\mu_{\mathrm{H^+}}^\circ(T,P) - G_{\mathrm{H}}^\circ(T_r,P_r)\right] ,
```

so that the apparent Gibbs energy of H⁺ is zero at every temperature and
pressure. This convention, due to Benson and Helgeson, is that of the databases
the package reads, which stores it as `ΔₐG⁰`, and as `ΔₐH⁰` the apparent
enthalpy defined in the same way. Both choices assign to the elements and to the
charge energies common to all species, and the argument of §5 cancels them from
any balanced reaction,

```math
\Delta_r G^\circ(T,P) = \sum_i \nu_i\,\Delta_f G_i^\circ(T,P)
  = \sum_i \nu_i\,\Delta_a G_i^\circ(T,P) .
```

The two definitions coincide at ``(T_r, P_r)``, where the tabulated value is
``\Delta_f G_i^\circ(T_r,P_r)``. Since the terms of the elements are constant,
the apparent energy of a neutral species changes with ``T`` and ``P`` as its
standard potential does, which the integral above gives explicitly,

```math
\Delta_a G_i^\circ(T,P) = \Delta_f G_i^\circ(T_r,P_r)
  - S_i^\circ(T_r,P_r)\,(T - T_r)
  + \int_{T_r}^{T} C_{P,i}^\circ\,\mathrm{d}T'
  - T\int_{T_r}^{T} \frac{C_{P,i}^\circ}{T'}\,\mathrm{d}T'
  + \int_{P_r}^{P} V_i^\circ\,\mathrm{d}P' ,
```

and the apparent enthalpy is ``\Delta_f H_i^\circ(T_r,P_r)`` plus the
increment given by the enthalpy relation above. For an ion, the term in
``\mu_{\mathrm{H^+}}^\circ(T,P)`` subtracts ``z_i`` times the variation of the
hydrogen ion, and the same expression holds with the conventional entropy, heat
capacity and volume of §5.

### [The chemical potential and the apparent Gibbs energy](@id sec-theory-mu-apparent)

Combining ``\mu_i = \mu_i^\circ + RT\ln a_i`` with the definition of the apparent
Gibbs energy relates the chemical potential of a species to the quantity a
database stores,

```math
\mu_i = \Delta_a G_i^\circ(T,P) + RT\ln a_i
  + \sum_e \alpha_{ei}\,G_e^\circ(T_r,P_r)
  + z_i\left[\mu_{\mathrm{H^+}}^\circ(T,P) - G_{\mathrm{H}}^\circ(T_r,P_r)\right] .
```

The last two terms are unknown, and the package leaves them out: the potential
it minimizes is ``\Delta_a G_i^\circ + RT\ln a_i``. The omission changes no
equilibrium. Summed over the species with their amounts, the omitted terms add
to the Gibbs energy of the system ``\sum_e b_e\,G_e^\circ(T_r,P_r) +
b_Z\,[\mu_{\mathrm{H^+}}^\circ(T,P) - G_{\mathrm{H}}^\circ(T_r,P_r)]``, where
``b_e = \sum_i \alpha_{ei}\,n_i`` is the amount of element ``e`` and ``b_Z =
\sum_i z_i\,n_i`` the charge. This sum takes the same value for every
composition that satisfies the conservation constraints of the minimization, so
that the two Gibbs energies reach their minimum at the same composition, and
only the multipliers of the constraints, the component potentials of
[Thermochemistry](@ref sec-theory-thermo) §4, are shifted, by
``G_e^\circ(T_r,P_r)/RT`` for element ``e`` and by the bracket divided by ``RT``
for the charge. With primary species as components, the omitted terms are a
combination of the ``A_{ci}`` by the decomposition of §5, and the conclusion
is the same.

The enthalpy is treated likewise. In a closed system without net charge, the
apparent enthalpy differs from the absolute one by ``\sum_e b_e\,H_e^\circ(T_r,P_r)``,
a constant, and heats computed as differences of apparent enthalpies are
therefore exact; the temperature derivative of an apparent enthalpy is the heat
capacity of the species itself, the elements contributing nothing.

One relation does not survive for a species on its own. Since the elements are
frozen at ``T_r``, the apparent Gibbs energy of a species divided by ``T`` does
not obey the Gibbs-Helmholtz relation of §2 with its apparent enthalpy: for a
neutral species the difference is ``\sum_e \alpha_{ei}\,[G_e^\circ(T_r) -
H_e^\circ(T_r)]/T^2 = -T_r \sum_e \alpha_{ei} S_e^\circ(T_r)/T^2``, the absolute
entropies of the elements at ``T_r`` appearing where their enthalpies and Gibbs
energies cancel. The relation is recovered for any balanced reaction, where these
terms cancel in turn. The slaking of lime shows both at 60 °C, a centered
difference giving the left-hand side, up to the misfit of the tabulated values
met in §5, which divides by ``T^2`` in the same way:

```@example basics
T, δ = 333.15, 0.01
# ∂(g/T)/∂T + h/T², which Gibbs-Helmholtz sets to zero, by a centered difference
gibbs_helmholtz(g, h) = (g(T + δ) / (T + δ) - g(T - δ) / (T - δ)) / 2δ + h(T) / T^2
# what the tabulated values miss of G = H - TS at Tr, divided by T²
misfit(g, h, s) = (h(Tr) - g(Tr) - Tr * s(Tr)) / T^2
@printf("slaking: %8.5f J/(mol K²), all of it the misfit %8.5f of the tabulated values\n",
        gibbs_helmholtz(t -> Δr(G, t), t -> Δr(H, t)), misfit(t -> Δr(G, t), t -> Δr(H, t), t -> Δr(S, t)))
@printf("lime:    %8.5f J/(mol K²) = %8.5f from the elements + %8.5f of misfit\n",
        gibbs_helmholtz(t -> G("Lim", t), t -> H("Lim", t)),
        -Tr * sum(k * S_e[string(e)] for (e, k) in atoms(db["Lim"])) / T^2,
        misfit(t -> G("Lim", t), t -> H("Lim", t), t -> ΔfS))
```

### The heat of reaction

The standard enthalpy of reaction is the heat given off at constant pressure when
the reaction proceeds by one mole, in the standard states, and it decides how
``K`` changes with temperature. Applying the relation of §3 to each constituent
and summing over a balanced reaction, where it holds for apparent quantities,
yields the van 't Hoff relation

```math
\frac{\mathrm{d}\ln K}{\mathrm{d}T} = \frac{\Delta_r H^\circ}{RT^2} ,
```

by which an exothermic reaction is displaced towards its reactants on heating.
The enthalpy of reaction varies with temperature in turn through the heat
capacities, ``\mathrm{d}\Delta_r H^\circ/\mathrm{d}T = \Delta_r C_p^\circ``
(Kirchhoff's relation), and a database that stores ``C_p^\circ(T)`` for every
species therefore determines ``K(T)`` completely (§13.3). The pressure acts
through the volume of reaction: differentiating ``\ln K = -\Delta_r G^\circ/RT``
at fixed temperature, with ``(\partial\mu_i^\circ/\partial P)_T = V_i^\circ``,
gives ``(\partial\ln K/\partial P)_T = -\Delta_r V^\circ/RT``.

The package computes heat without writing a reaction either. The enthalpy of a
state being ``\sum_i n_i h_i^\circ``, the heat released between two states at the
same temperature is the difference of their enthalpies, which is how
[`heat_release`](@ref) and the calorimeters compute it. The heat of the
calorimeters under partial equilibrium, the temperature shift of an equilibrium
in a semi-adiabatic cell and the sensitivities of the minimization all rest on
the two remarks made above on the omitted terms.

## Where to go next

[Thermochemistry](@ref sec-theory-thermo) takes up these quantities in the
notation of the code, with the chemical potential divided by ``RT`` that the
solver differentiates and the models by which ``\Delta_a G^\circ(T,P)`` is
evaluated, and [Standard states](@ref sec-theory-standard-states) states what
each activity is measured from.

  - [CEM I 52.5 N and slag in an isothermal calorimeter, read off the states](@ref sec-example-isothermal)
    and [A CEM I 52.5 N mortar in a semi-adiabatic calorimeter, inside the kinetics](@ref ex-semiadiabatic),
    the two calorimeters of §1 on a cement.
  - [Proving that an answer is the answer](@ref sec-theory-certificate), the
    optimality conditions from which the conditions of §3 follow.
