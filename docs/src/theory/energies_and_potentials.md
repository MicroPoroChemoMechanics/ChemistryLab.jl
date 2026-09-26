# [Energies, enthalpies and the chemical potential](@id sec-theory-basics)

The pages of this chapter use the chemical potential, the standard Gibbs energy
of a species and its enthalpy of formation as familiar objects, although the
relations between them are rarely obvious on a first reading, since a
calorimeter measures a change of enthalpy whereas the solver minimizes a Gibbs
energy assembled from chemical potentials. This page rebuilds them from the two laws of thermodynamics, in the
order in which [AndersonCrerar1993](@cite) introduce them (chapters 4, 5, 7, 9
and 14), and ends where [Thermochemistry](@ref sec-theory-thermo) begins, with
the quantities the code evaluates.

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
C_P = \left(\frac{\partial H}{\partial T}\right)_{P,n} ,
```

and it governs the other kind of calorimeter. In a vessel that exchanges no
heat, the enthalpy of the contents is constant, so that the enthalpy released by
the reaction at the initial temperature is found again as the heating of the
contents after the reaction, ``\Delta H(T_0) + \int_{T_0}^{T} C_P\,\mathrm{d}T' = 0``. A
semi-adiabatic vessel adds the heat it loses to this balance, and it is the
subject of [A semi-adiabatic calorimeter, inside the kinetics](@ref ex-semiadiabatic).

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
§4 defines it, together with the entropy of formation it must not be confused
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
the enthalpy. Extensivity then gives ``G = \sum_i n_i\mu_i`` and
``H = \sum_i n_i h_i`` (§9.2), so that the enthalpy of a system, and hence the
heat of a change at constant pressure, is obtained by weighting the partial
molar enthalpies by the amounts.

The chemical potential of a constituent of a mixture is written as a standard
part, the potential ``\mu_i^\circ(T,P)`` of a reference state, plus a
contribution of the composition through the activity ``a_i``,
``\mu_i = \mu_i^\circ + RT\ln a_i`` ([Standard states](@ref sec-theory-standard-states)).
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
both, and for any reaction ``\sum_i \nu_i A_i = 0`` among the constituents,
counted positive for the products,

```math
\Delta_r G \;\equiv\; \sum_i \nu_i\,\mu_i = 0 .
```

Away from equilibrium, ``\mathcal{A} = -\Delta_r G`` is the affinity of the
reaction, positive when it proceeds forward (§14.5), which is the quantity the
rate laws of [Rate laws](@ref sec-theory-kinetics) are written with.

## 4. Energies known up to a constant, and the reaction of formation

The internal energy, the enthalpy and the Gibbs energy are known up to an
additive constant only, since a measurement gives their changes and never their
values. A database cannot tabulate the change of every reaction of interest,
which are too many, and tabulates instead, for each substance, the change of the
one reaction that forms it from its elements, each taken in its stable form at
the conditions of reference, ``T_r = 298.15`` K and ``P_r = 1`` bar
[AndersonCrerar1993](@cite) (§7.1, §7.2): the standard enthalpy of formation
``\Delta_f H^\circ`` and the standard Gibbs energy of formation
``\Delta_f G^\circ``. For an element in its reference form, that reaction forms
the element from itself, and its enthalpy and Gibbs energy of formation are
zero. This follows from the definition, the absolute energies of the elements
being finite and unknown, whereas the absolute entropy of the same element is
not zero. The aqueous ions require
a genuine convention besides, because an ion cannot be formed without a
counter-charge: all the formation properties of H⁺ are set to zero, and those of
every other ion follow (§17.3).

### [Absolute entropy and entropy of formation](@id sec-theory-absolute-entropy)

The entropy escapes the indeterminacy of the energies, because the third law
fixes its zero: the entropy of a pure substance in a perfect crystalline form
tends to zero as the temperature tends to 0 K [AndersonCrerar1993](@cite)
(§6.5). The entropy at any other temperature then follows from the second law,
``\mathrm{d}S = \delta q_{\rm rev}/T = C_P\,\mathrm{d}T/T`` at constant pressure,
integrated from 0 K, with the entropy of each phase transition met on the way,
the ratio of its enthalpy to the temperature at which it occurs,

```math
S^\circ(T) = \int_0^{T} \frac{C_P^\circ(T')}{T'}\,\mathrm{d}T'
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
\Delta_f S_i^\circ = S_i^\circ - \sum_e \nu_{ei}\, S_e^\circ ,
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
canceling as in Hess's law below. The two are easily confused, and they differ
even in sign: the absolute entropy of lime is positive, its entropy of formation
from calcium metal and oxygen gas negative, since the gas it consumes is far more
disordered than the solid it forms.

The aqueous ions carry a third kind of entropy. The convention that sets the
formation properties of H⁺ to zero sets its tabulated entropy to zero as well, and
the entropy tabulated for an ion of charge ``z_i`` is then a conventional
entropy, ``S_i^\circ - z_i\,S_{\mathrm{H^+}}^\circ``, its absolute entropy less
``z_i`` times that of the hydrogen ion, which thermodynamics alone does not
determine (§17.3). A conventional entropy can be negative, as that of Ca²⁺ is,
and the unknown term cancels from any reaction, whose charges balance,
``\sum_i \nu_i z_i = 0``.

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

Because ``H`` and ``G`` are state functions, the change in any balanced reaction
is the same whether the reaction is carried out directly or by first
decomposing the reactants into their elements and then forming the products
from them. This is Hess's law, and it yields

```math
\Delta_r H^\circ = \sum_i \nu_i\,\Delta_f H_i^\circ ,
\qquad
\Delta_r G^\circ = \sum_i \nu_i\,\Delta_f G_i^\circ ,
```

the elements appearing on both sides of the balanced reaction and canceling
between the two steps. With the entropy of reaction defined above, the three
quantities of a reaction are bound by ``\Delta_r G^\circ = \Delta_r H^\circ -
T\Delta_r S^\circ``, which the slaking of lime, CaO + H₂O → Ca(OH)₂, verifies
with the same values:

```@example basics
slaking = ("Portlandite" => 1, "Lim" => -1, "H2O@" => -1)
Δr(f, T) = sum(ν * f(s, T) for (s, ν) in slaking)
@printf("slaking: ΔrH = %.1f kJ/mol, ΔrG = %.1f kJ/mol, ΔrS = %.2f J/(mol K); ΔrG - (ΔrH - Tr ΔrS) = %.0f J/mol\n",
        Δr(H, Tr) / 1000, Δr(G, Tr) / 1000, Δr(S, Tr), Δr(G, Tr) - (Δr(H, Tr) - Tr * Δr(S, Tr)))
```

The reference need not be the elements. For silicates, oxides, carbonates and
sulfates, the enthalpy of formation from the oxides, the change of the reaction
forming the compound from its constituent oxides, is also in use (§7.4.4). Its
values are much smaller in magnitude, and so are their uncertainties, which no
longer include those of forming the oxides from the elements; the two references
cannot be mixed within one calculation. It is the natural scale for a glass, and
[An isothermal calorimeter, read off the states](@ref sec-example-isothermal)
places the slag of a blended cement on it.

## 5. Reactions, equilibrium constants and the heat of reaction

Inserting ``\mu_i = \mu_i^\circ + RT\ln a_i`` into ``\Delta_r G`` separates a
standard part from the composition,

```math
\Delta_r G = \Delta_r G^\circ + RT\ln Q ,
\qquad
Q = \prod_i a_i^{\nu_i} ,
```

and the condition of equilibrium ``\Delta_r G = 0`` fixes the value ``Q`` takes
there, the equilibrium constant [AndersonCrerar1993](@cite) (§13.1),

```math
\ln K = -\frac{\Delta_r G^\circ}{RT} .
```

The standard Gibbs energy of reaction decides which side a reaction favors
between reactants and products in their standard states, and the equilibrium
constant expresses the same thing in terms of activities. The standard enthalpy
of reaction plays a different part: it is the heat given off at constant
pressure when the reaction proceeds by one mole, in the standard states, and it
decides how ``K`` changes with temperature. Applying the relation of §3 to each
constituent and summing yields the van 't Hoff relation

```math
\frac{\mathrm{d}\ln K}{\mathrm{d}T} = \frac{\Delta_r H^\circ}{RT^2} ,
```

by which an exothermic reaction is displaced towards its reactants on heating.
The enthalpy of reaction varies with temperature in turn through the heat
capacities, ``\mathrm{d}\Delta_r H^\circ/\mathrm{d}T = \Delta_r C_P^\circ``
(Kirchhoff's relation), and a database that stores ``C_P^\circ(T)`` for every
species therefore determines ``K(T)`` completely (§13.3).

The package itself writes no reaction to find an equilibrium: it minimizes
``G`` over all the constituents at once, and the conditions ``\Delta_r G = 0`` of
every reaction that can be written among them follow from the optimality
conditions of the minimization ([Proving that an answer is the answer](@ref sec-theory-certificate)).
The same holds for heat. The enthalpy of a state being ``\sum_i n_i h_i^\circ``,
the heat released between two states at the same temperature is the difference
of their enthalpies, with no reaction to write down, which is how
[`heat_release`](@ref) and the calorimeters compute it.

## 6. Temperature, and the apparent quantities of a database

The enthalpy and the Gibbs energy of formation at a temperature ``T`` other than
``T_r`` can be defined in two ways (§7.4). The standard quantities of formation
``\Delta_f H^\circ(T)`` and ``\Delta_f G^\circ(T)`` compare the substance with
its elements at the same temperature ``T``, which requires the heat capacities of
the elements, and their phase changes, over the whole range. The apparent
quantities of formation leave the elements at ``T_r``,

```math
\Delta_a H_i^\circ(T) = H_i^\circ(T) - \sum_e \nu_{ei}\,H_e^\circ(T_r) ,
\qquad
\Delta_a G_i^\circ(T) = G_i^\circ(T) - \sum_e \nu_{ei}\,G_e^\circ(T_r) ,
```

where ``\nu_{ei}`` is the number of atoms of element ``e`` in species ``i``. The
two definitions coincide at ``T_r``, and since the terms of the elements are
constants of each species, they cancel from every balanced reaction: Hess's law
and the relations of §5 hold with apparent quantities exactly as with absolute
ones. This is the convention of the databases the package reads, the one it
stores as `ΔₐH⁰` and `ΔₐG⁰`, and [Apparent and formation Gibbs energies](@ref sec-theory-apparent)
gives its evaluation from the heat capacity.

One relation does not survive for a species on its own. Since the elements are
frozen at ``T_r``, the apparent Gibbs energy of a species divided by ``T`` does
not obey the Gibbs-Helmholtz relation of §2 with its apparent enthalpy: the
difference is ``\sum_e \nu_{ei}\,[G_e^\circ(T_r) - H_e^\circ(T_r)]/T^2 = -T_r
\sum_e \nu_{ei} S_e^\circ(T_r)/T^2``, the absolute entropies of the elements at
``T_r`` appearing where their enthalpies and Gibbs energies cancel. The relation
is recovered for any balanced reaction, where these terms cancel in turn. The
slaking of lime shows both at 60 °C, a centered difference giving the left-hand
side, up to the misfit of the tabulated values met above, which divides by
``T^2`` in the same way:

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

The same terms enter the chemical potential of every species as a sum over its
elements, ``\sum_e \nu_{ei} c_e``, with one constant ``c_e`` per element. At
equilibrium the chemical potentials are combinations of the potentials of the
elements, ``\mu_i = \sum_e \nu_{ei}\, y_e`` ([Thermochemistry](@ref sec-theory-thermo)
§4), so that such terms are absorbed by the element potentials and change
neither the equilibrium composition nor its derivatives. For the enthalpy, the
apparent enthalpy of a closed system differs from its absolute value by
``\sum_e b_e H_e^\circ(T_r)``, the amounts ``b_e`` of the elements being
conserved, and heats computed as differences of apparent enthalpies are exact;
the temperature derivative of an apparent enthalpy is the heat capacity of the
species itself, the elements contributing nothing. The heat of the calorimeters
under partial equilibrium, the temperature shift of an equilibrium in a
semi-adiabatic cell and the sensitivities of the minimization all rest on these
two remarks.

## Where to go next

[Thermochemistry](@ref sec-theory-thermo) takes up these quantities in the
notation of the code, with the chemical potential divided by ``RT`` that the
solver differentiates, and [Standard states](@ref sec-theory-standard-states)
states what each activity is measured from.

  - [An isothermal calorimeter, read off the states](@ref sec-example-isothermal)
    and [A semi-adiabatic calorimeter, inside the kinetics](@ref ex-semiadiabatic),
    the two calorimeters of §1 on a cement.
  - [Proving that an answer is the answer](@ref sec-theory-certificate), the
    optimality conditions of §5.
