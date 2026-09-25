# [Getting started](@id sec-quickstart)

This page takes a first calculation from installation to a solved equilibrium in
three steps: the thermodynamic data of a few species are read from a database
distributed with the package, the reaction between them is written and its
equilibrium constant evaluated, and the equilibrium of calcite with water is
finally computed. No knowledge of chemical thermodynamics is assumed. Each
quantity is defined in the chapter [Theory](@ref sec-theory), to which the text
refers at the point where the quantity is first needed.

## Installation

ChemistryLab.jl is registered in the General registry and is installed with the
Julia package manager, either from the Pkg mode of the REPL (entered by typing
`]`)

```julia
pkg> add ChemistryLab
```

or, equivalently, through the `Pkg` API:

```julia
julia> import Pkg; Pkg.add("ChemistryLab")
```

## A first reaction: calcite in water

The dissolution of calcite is written

$\ce{CaCO3 <=> Ca^2+ + CO3^2-}$

and its equilibrium constant ``K``, here a solubility product, is related to the
standard Gibbs energy of reaction by

```math
\Delta_r G^\circ = -RT \ln K ,
\qquad
\Delta_r G^\circ = \sum_i \nu_i\, \Delta_a G_i^\circ
```

where ``\nu_i`` are the stoichiometric coefficients of the reaction, negative
for the reactants, and ``\Delta_a G_i^\circ`` the Gibbs energies of the species
at the temperature ``T`` and the reference pressure of 1 bar. The subscript
``a`` stands for *apparent*: the package works with apparent Gibbs energies of
formation, which coincide with the ordinary Gibbs energies of formation at the
reference temperature of 298.15 K and differ from them elsewhere. The
distinction, and the reason for it, are the subject of
[Apparent and formation Gibbs energies](@ref sec-theory-apparent).

The species are read from CEMDATA18 [Lothenbach2019](@cite), a database for
cement systems whose file is distributed with the package as a copy of the one
maintained on [ThermoHub](https://github.com/thermohub). `build_species` reads
the file and returns a vector of species whose thermodynamic functions are
already compiled; `speciation` then keeps the species whose elements are all
found among those of a few seed species, here calcium, carbon, hydrogen and
oxygen from `Cal`, `H2O@` and `CO2`. The three aqueous gases listed in
`exclude_species` are left out because their presence would open a redox
equilibrium that this example does not need.

```@example from_scratch
using ChemistryLab
using DynamicQuantities # for unit management

all_species = build_species(datapath("cemdata18-thermofun.json"); verbose = false)
species_calcite = speciation(all_species, split("Cal H2O@ CO2");
                             aggregate_state=[AS_AQUEOUS],
                             exclude_species=split("H2@ O2@ CH4@"))
dict_species_calcite = Dict(symbol(s) => s for s in species_calcite)
```

Each species carries its molar mass and its standard thermodynamic functions of
temperature: heat capacity, entropy, enthalpy and Gibbs energy. They can be
inspected and plotted directly.

```@example from_scratch
dict_species_calcite["Cal"]
```

```@example from_scratch
using Plots

p1 = plot(xlabel="Temperature [°C]", ylabel="Cp⁰ [J/mol/K]", title="Heat capacity of calcite \nas a function of temperature")
plot!(p1, θ -> dict_species_calcite["Cal"].Cp⁰(T = θ*ua"degC"), 0:0.1:100, label="Cp⁰")
```

Writing the reactions requires a choice of independent species, the
*primaries*, in terms of which every other species is expressed. The
stoichiometric matrix collects those decompositions, one column per species and
one row per primary; its construction is detailed in
[Stoichiometric matrices](@ref ex-stoich-matrix).

```@example from_scratch
primaries = [dict_species_calcite[s] for s in split("H2O@ H+ CO3-2 Ca+2")]
SM = StoichMatrix(collect(values(dict_species_calcite)), primaries)
pprint(SM)
```

The reactions follow from the matrix, each carrying its own thermodynamic
functions of temperature, among which the equilibrium constant.

```@example from_scratch
list_reactions = reactions(SM)
dict_reactions_calcite = Dict(r.symbol => r for r in list_reactions)
```

```@example from_scratch
dict_reactions_calcite["Cal"].logK⁰
```

```@example from_scratch
p2 = plot(xlabel="Temperature [°C]", ylabel="pKs", title="Solubility product (pKs) of calcite \nas a function of temperature")
plot!(p2, θ -> dict_reactions_calcite["Cal"].logK⁰(T = θ*ua"degC"), 0:0.1:100, label="pKs")
```

## A first equilibrium

The previous section evaluated a property of one reaction. Solving the
equilibrium of a system is a different question: the amounts of all species are
sought that minimize the Gibbs energy of the system under the conservation of
its elements. Three objects are involved.

| object | role |
|:--|:--|
| [`ChemicalSystem`](@ref) | the immutable description of the system: species, primaries, stoichiometric matrices; built once and reused |
| [`ChemicalState`](@ref) | the mutable state: amounts in mol, temperature and pressure |
| [`equilibrate`](@ref) | the solver, which tries every available route and returns the state that [`optimality_certificate`](@ref) proves to be the minimum; [`equilibrate_certified`](@ref) returns the certificate as well |

```@example from_scratch
using Optimization, OptimizationIpopt
using DynamicQuantities

# ChemicalSystem: declare the species and which four are the independent basis
primaries_eq = [dict_species_calcite[s] for s in split("H2O@ H+ CO3-2 Ca+2")]
cs = ChemicalSystem(collect(values(dict_species_calcite)), primaries_eq)

# ChemicalState: set the initial amounts
state = ChemicalState(cs)
set_quantity!(state, "Cal",  1e-3u"mol")   # 1 mmol of calcite
set_quantity!(state, "H2O@", 1.0u"kg")     # 1 kg of water

# Seed H⁺ and OH⁻ at pH 4 (trace amounts to help convergence)
V = volume(state)
set_quantity!(state, "H+",  1e-4u"mol/L" * V.liquid)
set_quantity!(state, "OH-", 1e-10u"mol/L" * V.liquid)

# Solve: find the Gibbs-energy minimum
state_eq = equilibrate(state)
nothing # hide
```

```@raw html
<details><summary>The solved state in full — every species, with its amount</summary>
```

```@example from_scratch
state_eq
```

```@raw html
</details>
```

```@example from_scratch
println("pH = ", round(pH(state_eq), digits = 2))
```

The returned state gives access to every derived quantity: pH, the volumes of
the phases, the amount of each species. Why the minimum exists, why it is
unique, and how the solver proves that the state it returns is that minimum are
explained in [Proving that an answer is the answer](@ref sec-theory-certificate);
the options of the solver are described in
[Chemical Equilibrium](@ref sec-equilibrium).

## Where to go next

The documentation can be entered from three directions, depending on what the
reader already knows.

- A reader new to chemical thermodynamics is best served by the chapter
  [Theory](@ref sec-theory) read in its stated order, starting from
  [Thermochemistry](@ref sec-theory-thermo), before the Manual, which describes
  one kind of object per page from [Species](@ref sec-species) onwards.
- A reader who knows what is to be computed can go to the tutorial
  [Chemical Equilibrium](@ref sec-equilibrium), then to the aqueous examples,
  whose answers can be checked against closed forms, starting with
  [CO₂ dissolution and the carbonate system](@ref sec-co2-carbonate).
- A reader holding the oxide analysis of a cement will find the route from an
  analysis to a chemical system in [Bogue calculation](@ref ex-bogue), and the
  worked binders from [A CEM I from its clinker phases](@ref sec-cem1-from-clinker)
  onwards.

## Citing ChemistryLab

When ChemistryLab is used in published work, it is to be cited as follows.

```bibtex
@software{chemistrylab_jl,
  author       = {Barthélémy, Jean-François and
                  Soive, Anthony},
  title        = {ChemistryLab.jl: Numerical laboratory for
                   computational chemistry},
  doi          = {10.5281/zenodo.17756074},
  url          = {https://doi.org/10.5281/zenodo.17756074},
}
```
