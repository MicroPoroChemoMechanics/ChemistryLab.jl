# Chemical systems and states

```@index
Pages = ["chemical_systems.md"]
```

## ChemicalSystem

```@autodocs
Modules = [ChemistryLab]
Pages = ["chemical_structs/chemical_systems.jl"]
```

## ChemicalState

```@autodocs
Modules = [ChemistryLab]
Pages = ["chemical_structs/chemical_states.jl"]
```

## Volume fractions

```@autodocs
Modules = [ChemistryLab]
Pages = ["chemical_structs/volume_fractions.jl"]
```

## An oxide analysis as an element budget

```@autodocs
Modules = [ChemistryLab]
Pages = ["chemical_structs/oxide_budget.jl"]
```

## Loss on ignition

What a solid assemblage would lose on heating, from the formulas its phases
carry — the total a thermogram integrates to. See
[The water budget of a hydrating paste](@ref sec-theory-water-budget) for what
this gives and what a full thermogravimetric operator would still need.

```@autodocs
Modules = [ChemistryLab]
Pages   = ["chemical_structs/ignition_loss.jl"]
```

## Thermogravimetry

Turning that total into a curve: a decomposition window per phase, carrying its
own provenance, and the operator that makes them **identifiable from a measured
thermogram** rather than only supplied. See
[A thermogram, and the windows it takes to have one](@ref sec-example-tga).

```@autodocs
Modules = [ChemistryLab]
Pages   = ["chemical_structs/thermogram.jl"]
```
