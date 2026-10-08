# [Formation quantities and the database](@id sec-theory-formation)

[The chemical potential and reactions](@ref sec-theory-basics) writes every
equilibrium in terms of standard potentials ``\mu_i^\circ``, which no
measurement gives. This page shows what a database tabulates instead, how those
numbers are measured, and how they are carried to any temperature and pressure,
up to the apparent Gibbs energy the package reads, which is where
[Thermochemistry](@ref sec-theory-thermo) begins. It is the most technical of the
foundations, and a first reading can come back to it later.

!!! tip "Questions this page answers"
      - Energies are known only up to a constant: what can a database tabulate,
        and why is it enough?
      - How are enthalpies, entropies and Gibbs energies of formation measured?
      - What separates an absolute entropy, an entropy of formation and the
        conventional entropy of an ion?
      - How does a standard potential change with temperature and pressure, and
        what is the apparent Gibbs energy a database stores?

## 1. Energies known up to a constant, and the reactions of formation

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
argument below also delimits the freedom left in their definition, of which §2
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
No value is assigned to them: the zeros follow from the definition itself.
For an element in its reference form, ``\mu_i^\circ`` is ``G_e^\circ`` times
its number of atoms, and the formula gives zero; for H⁺, one atom of hydrogen and
one charge give ``\mu_{\mathrm{H^+}}^\circ - G_{\mathrm{H}}^\circ - G_Z^\circ = 0``,
whatever the values of ``G_{\mathrm{H}}^\circ`` and ``G_Z^\circ``. That the
energy of formation of H⁺ is zero is the convention on ions, and it is built
into the reaction of formation chosen for an ion rather than into a value
given to an unknown: the energies of formation of all other ions follow from
it. Inserted in ``\Delta_r G^\circ =
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

The entropy is in a different position from the energies. Calorimetry gives
only its changes, but a zero common to all substances can be chosen: the third
law takes the entropy of a pure substance in a perfect crystalline form, in
internal equilibrium, to tend to zero as the temperature tends to 0 K
[AndersonCrerar1993](@cite) (§6.5). The choice leaves out what no chemical
reaction changes, down to the constituents of the nuclei, and is a convention
rather than a measurement [Richet2001; Sec. 4.4b, p. 75](@cite). The entropy at any other temperature then follows from the second law,
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
the integral (§6.5.3); the glass of a blast-furnace slag is of this kind. That
entropy is not a property of the substance alone: it is the configurational
entropy the phase had when its structure stopped relaxing on cooling, and it
depends on the temperature at which that happened, so on the thermal history of
the sample [Richet2001; Sec. 6.3a, Eq. (6.38), p. 144](@cite).

The residual entropy escapes the calorimetry from 0 K, which starts from the
glass as it already is, disorder included. It is obtained by closing a cycle
with the crystal of the same composition. Both are heated to the melting
temperature ``T_m``: the crystal melts there, while the glass has turned into a
supercooled liquid and then into the same liquid, so the two paths end in one
state. Equating the entropies of that liquid gives

```math
S_{\rm res} = \int_0^{T_m} \frac{C_{p,\text{crystal}}}{T}\,\mathrm{d}T
  + \frac{\Delta_{\rm fus} H}{T_m}
  - \int_0^{T_m} \frac{C_{p,\text{glass}}}{T}\,\mathrm{d}T ,
```

the last heat capacity being that of the glass and, above its glass
transition, of the supercooled liquid [Richet2001; Sec. 4.4c, p. 76](@cite).
The cycle crosses the glass transition, which is not reversible, and integrating
``C_p/T`` across it treats it as if it were; the entropy produced there is small
against the residual entropy, which is why the cycle is accepted
[Richet2001; Sec. 6.4d, p. 152](@cite); [RichetBottinga1986](@citet) put it at
least an order of magnitude below. None of this enters a calculation of
this package. No database carries the Gibbs energy of a slag glass: the glass
enters the budget through its oxides, and its reaction is prescribed, by a rate
law or by an imposed degree of reaction, rather than decided by its energy.
What can be read off a measurement is its enthalpy, which
[CEM I 52.5 N and slag in an isothermal calorimeter, read off the states](@ref sec-example-isothermal)
derives from the heat the pastes release.

The entropy of formation is of another nature: it is the change of entropy in
the reaction of formation,

```math
\Delta_f S_i^\circ = S_i^\circ - \sum_e \alpha_{ei}\, S_e^\circ ,
```

where ``S_e^\circ`` is the absolute entropy of element ``e`` in its reference
form, counted per atom, half that of O₂ gas for oxygen. It is the entropy of
formation that is bound to the enthalpy and the Gibbs energy of formation by the
relation of [The two laws, and what the Gibbs energy measures](@ref sec-theory-laws)
§5 applied to the reaction of formation,

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
the same way, and they describe the ion at infinite dilution: those of OH⁻, its
entropy, heat capacity and volume, are all negative, which reflects the water
around the ion rather than an error of the data
[Richet2001; Sec. 12.3b, p. 296](@cite).

!!! note "How these numbers are measured"
    No instrument reads a Gibbs energy. Every tabulated value combines a few
    kinds of measurement [AndersonCrerar1993](@cite).

      - **Enthalpies of formation** come from calorimetry. When a substance
        burns cleanly from its elements, the heat of combustion gives its
        enthalpy of formation directly; a bomb calorimeter works at constant
        volume, so it measures ``\Delta U``, which the change of ``PV`` converts
        to ``\Delta H`` ([The two laws, and what the Gibbs energy measures](@ref sec-theory-laws)
        §2). Most substances cannot be formed that way, and their enthalpy is
        obtained by Hess's law from reactions that can be measured, typically
        the dissolution of the compound and of its oxides in the same acid: the
        heats combine as the reactions do.
      - **Absolute entropies** come from the heat capacity, measured from a few
        kelvins upward and integrated as above.
      - **Gibbs energies of formation** follow from the two, as
        ``\Delta_f G^\circ = \Delta_f H^\circ - T\Delta_f S^\circ``, or from an
        equilibrium measured directly. A solubility gives ``\ln K`` and hence
        ``\Delta_r G^\circ = -RT\ln K``; the reversible voltage ``E^\circ`` of an
        electrochemical cell gives ``\Delta_r G^\circ = -\nu_e F E^\circ``,
        ``\nu_e`` being the number of electrons the reaction transfers. Either
        fixes one unknown energy of formation once those of the other
        participants are known.

    The hydrates of a cement can be neither burned nor formed from their
    elements, and most of the Gibbs energies CEMDATA18 gives them are derived
    from measured solubilities [Lothenbach2019](@cite). A solubility measured at
    several temperatures yields the enthalpy of reaction as well, through the
    van 't Hoff relation of §2.

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

### [Energies counted from the primaries: the gauge of a database of reactions](@id sec-theory-gauge)

A database of reactions, of which PHREEQC's are the common example, does not
give the formation of its species from the elements but only from its master
species: the log K of each reaction, ``L_i(T) = \log_{10} K_i(T)``, as a function
of temperature, with the master species ``\mathrm{B}_c`` defined by identity
reactions (`Ca+2 = Ca+2`). By the relation above, it determines

```math
\mu_i^\circ(T) - \sum_c A_{ci}\,\mu_c^\circ(T) = -RT\ln 10\; L_i(T)
```

and nothing else: the ``\mu_c^\circ`` of the master species are not in it. The
package reads such a database by setting them to zero at every temperature,

```math
\mu_c^\circ(T) = 0 , \qquad \mu_i^\circ(T) = -RT\ln 10\; L_i(T) ,
```

which is a choice of zero, a gauge, and not a measurement. It does not change
an equilibrium. The energies so defined differ from any others consistent with
the database by

```math
\delta\mu_i^\circ = \sum_c A_{ci}\,\lambda_c ,
\qquad \lambda_c = -\mu_c^\circ(T) ,
```

a combination of the rows of the conservation matrix, the same for every
species, solute or mineral, since each decomposes on the same master species.
On a composition ``\mathbf{n}`` that meets the balances ``\mathbf{A}\mathbf{n} =
\mathbf{b}``, the Gibbs energy changes by

```math
\sum_i n_i\,\delta\mu_i^\circ = \sum_c \lambda_c \sum_i A_{ci}\,n_i = \sum_c \lambda_c\,b_c ,
```

a constant of the feasible set, and the composition that minimizes it is the
same. The master species of a database span the elements and the charge (one
per element, H⁺ for hydrogen and H₂O for oxygen, with the electron for the
charge), so that every species of the database decomposes on them and the
argument applies to all.

It applies within one database only. Species taken from two databases, or from a
database of reactions and a database of formation properties, carry energies
whose differences from the formation energies are two different combinations
``\sum_c A_{ci}\lambda_c`` and ``\sum_c A_{ci}\lambda'_c``. Their sum over a
feasible composition is no longer a constant, and the equilibrium is changed by
an amount that has no physical meaning: a master species of PHREEQC is at zero
where its formation energy from the elements is several hundred kilojoules per
mole. Each species read from a database therefore records its gauge under
`:gauge` (the database file and its digest for a database of reactions, the
formation from the elements for ThermoFun), and [`ChemicalSystem`](@ref) refuses
a system that mixes two.

The other standard quantities follow from the temperature dependence of
``L_i``, by the Gibbs–Helmholtz relation, and are counted in the same gauge:

```math
\Delta_a H_i^\circ = RT^2 \ln 10\,\frac{\mathrm{d}L_i}{\mathrm{d}T} ,
\qquad
S_i^\circ = \frac{\Delta_a H_i^\circ - \mu_i^\circ}{T} ,
\qquad
C_{p,i}^\circ = R\ln 10\left(2T\,\frac{\mathrm{d}L_i}{\mathrm{d}T} + T^2\,\frac{\mathrm{d}^2L_i}{\mathrm{d}T^2}\right) .
```

With phreeqc.dat, the master species Ca²⁺ and CO₃²⁻ are at zero, and the Gibbs
energy of calcite is that of its dissolution constant, ``\log_{10} K = -8.48``:

```@example basics
pdb = read_phreeqc_database(datapath("phreeqc.dat"))
ps = Dict(symbol(s) => s for s in build_species(pdb, ["Ca+2", "CO3-2", "CaCO3", "Calcite"]))
for s in ("Ca+2", "CO3-2", "CaCO3@", "Calcite")
    @printf("%-8s %9.3f kJ/mol\n", s, ps[s][:ΔₐG⁰](T = Tr) / 1000)
end
@printf("log10 K of the dissolution of calcite: %.2f\n", ps["Calcite"][:ΔₐG⁰](T = Tr) / (R_GAS * Tr * log(10)))
```

## 2. Temperature, pressure, and the apparent quantities of a database

### The standard potential at any temperature and pressure

The standard state of a species being taken at the temperature and the pressure
of the system, its potential obeys the relation of
[The two laws, and what the Gibbs energy measures](@ref sec-theory-laws) §5
for a substance of fixed
composition, with the standard molar entropy and volume,

```math
\mathrm{d}\mu_i^\circ = -S_i^\circ\,\mathrm{d}T + V_i^\circ\,\mathrm{d}P ,
```

whereas a gas, whose standard state stays at ``P_r`` whatever the pressure,
keeps the first term only and receives the pressure through its activity
([Standard states](@ref sec-theory-standard-states)). The relation is integrated
from ``(T_r, P_r)`` in temperature at ``P_r``, where
``S_i^\circ(T') = S_i^\circ(T_r) + \int_{T_r}^{T'} C_{P,i}^\circ\,\mathrm{d}T''/T''``,
and then in pressure at ``T``. The order follows the data: heat capacities are
measured near 1 bar, and integrating in pressure first would need ``C_p`` at
``P``, which only the equation of state gives, through
``(\partial C_p/\partial P)_T = -T\,(\partial^2 V/\partial T^2)_P``
[Richet2001; Sec. 8.1a, Eq. (8.4), p. 174](@cite). The temperature step is a double integral,

```math
-\int_{T_r}^{T} S_i^\circ(T')\,\mathrm{d}T'
  = -S_i^\circ(T_r)\,(T - T_r)
  - \int_{T_r}^{T}\!\!\int_{T_r}^{T'} \frac{C_{P,i}^\circ(T'')}{T''}\,\mathrm{d}T''\,\mathrm{d}T' ,
```

and exchanging the order of its two integrations, ``T''`` running from ``T_r``
to ``T`` and ``T'`` from ``T''`` to ``T``, turns the last term into a single
integral,

```math
\int_{T_r}^{T} \frac{C_{P,i}^\circ(T'')}{T''}\,(T - T'')\,\mathrm{d}T''
  = T\int_{T_r}^{T} \frac{C_{P,i}^\circ}{T''}\,\mathrm{d}T''
  - \int_{T_r}^{T} C_{P,i}^\circ\,\mathrm{d}T'' .
```

With the pressure step added, this gives [AndersonCrerar1993](@cite) (§7.4, §7.6)

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

The entropy that multiplies ``T - T_r`` is the absolute entropy of §1, and
both expressions involve properties of the species alone.

### [Apparent and formation Gibbs energies](@id sec-theory-apparent)

Away from ``(T_r, P_r)``, an energy of formation depends on the conditions at
which the elements are taken, and two choices are in use
[AndersonCrerar1993](@cite) (§7.4.1, §7.4.2). The traditional quantities of
formation ``\Delta_f G^\circ(T,P)`` take the elements at the temperature ``T``
and at 1 bar, which requires their heat capacities, and their phase changes,
over the whole range. The apparent quantities of formation leave them at
``(T_r, P_r)``: in the definition of §1, ``G_e^\circ`` becomes the constant
``G_e^\circ(T_r,P_r)`` and ``G_Z^\circ`` becomes
``\mu_{\mathrm{H^+}}^\circ(T,P) - G_{\mathrm{H}}^\circ(T_r,P_r)``,

```math
\Delta_a G_i^\circ(T,P) = \mu_i^\circ(T,P) - \sum_e \alpha_{ei}\,G_e^\circ(T_r,P_r)
  - z_i\left[\mu_{\mathrm{H^+}}^\circ(T,P) - G_{\mathrm{H}}^\circ(T_r,P_r)\right] ,
```

so that the apparent Gibbs energy of H⁺ is zero at every temperature and
pressure. This convention is that of the databases the package reads, which stores it as `ΔₐG⁰`, and as `ΔₐH⁰` the apparent
enthalpy defined in the same way. Both choices assign to the elements and to the
charge energies common to all species, and the argument of §1 cancels them from
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
capacity and volume of §1.

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
combination of the ``A_{ci}`` by the decomposition of §1, and the conclusion
is the same.

The enthalpy is treated likewise. In a closed system without net charge, the
apparent enthalpy differs from the absolute one by ``\sum_e b_e\,H_e^\circ(T_r,P_r)``,
a constant, and heats computed as differences of apparent enthalpies are
therefore exact; the temperature derivative of an apparent enthalpy is the heat
capacity of the species itself, the elements contributing nothing.

One relation does not survive for a species on its own. Since the elements are
frozen at ``T_r``, the apparent Gibbs energy of a species divided by ``T`` does
not obey the Gibbs-Helmholtz relation
([The two laws, and what the Gibbs energy measures](@ref sec-theory-laws) §5) with its apparent enthalpy: for a
neutral species the difference is ``\sum_e \alpha_{ei}\,[G_e^\circ(T_r) -
H_e^\circ(T_r)]/T^2 = -T_r \sum_e \alpha_{ei} S_e^\circ(T_r)/T^2``, the absolute
entropies of the elements at ``T_r`` appearing where their enthalpies and Gibbs
energies cancel. The relation is recovered for any balanced reaction, where these
terms cancel in turn. The slaking of lime shows both at 60 °C, a centered
difference giving the left-hand side, up to the misfit of the tabulated values
met in §1, which divides by ``T^2`` in the same way:

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

### [The heat of reaction](@id sec-theory-heat-of-reaction)

The standard enthalpy of reaction is the heat the system receives at constant
pressure when the reaction proceeds by one mole, in the standard states; the heat
it releases is ``-\Delta_r H^\circ``, and a reaction is exothermic when
``\Delta_r H^\circ < 0`` [Richet2001; Sec. 4.5b, p. 80](@cite). It decides how
``K`` changes with temperature. Applying the relation of
[The chemical potential and reactions](@ref sec-theory-basics) §1 to each constituent
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
calorimeters under partial equilibrium, the heat capacity of a semi-adiabatic
cell at equilibrium and the sensitivities of the minimization all rest on the
two remarks made above on the omitted terms
([Kinetics under partial equilibrium](@ref sec-theory-pe-kinetics)).

## Where to go next

[Thermochemistry](@ref sec-theory-thermo) takes up these quantities in the
notation of the code, with the chemical potential divided by ``RT`` that the
solver differentiates and the models by which ``\Delta_a G^\circ(T,P)`` is
evaluated, and [Standard states](@ref sec-theory-standard-states) states what
each activity is measured from.

  - [CEM I 52.5 N and slag in an isothermal calorimeter, read off the states](@ref sec-example-isothermal),
    which places the glass of a slag on the scale of formation from the oxides.
