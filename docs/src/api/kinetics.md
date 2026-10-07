# Kinetics API

Chemical kinetics module: rate models, kinetic reactions, ODE problem setup, solvers, and calorimetry.

## Rate constants and models

```@autodocs
Modules = [ChemistryLab]
Pages   = ["kinetics/rate_models.jl"]
```

## Kinetic reactions

The area models the rate factories consume live on their own page, since they
are shared with surface chemistry: see [Surfaces API](@ref).

```@autodocs
Modules = [ChemistryLab]
Pages   = ["kinetics/kinetics_reactions.jl"]
```

## A slow surface

A reaction moving a site from one state to another, at a rate that vanishes
where the equilibrium puts the surface. See
[A slow surface](@ref sec-theory-pe-slow-surface).

```@autodocs
Modules = [ChemistryLab]
Pages   = ["kinetics/sorption_rates.jl"]
```

## Kinetics problem

```@autodocs
Modules = [ChemistryLab]
Pages   = ["kinetics/kinetics_problems.jl"]
```

## Kinetics solver

```@autodocs
Modules = [ChemistryLab]
Pages   = ["kinetics/kinetics_solver.jl"]
```

## Calorimetry

```@autodocs
Modules = [ChemistryLab]
Pages   = ["kinetics/calorimetry.jl"]
```

## Post-processing a solution

```@autodocs
Modules = [ChemistryLab]
Pages   = ["kinetics/kinetics_postprocessing.jl"]
```

## The second law along a run

The affinity of each kinetic reaction and the power it dissipates, on the
certified compositions of a run. See
[What a rate law may be](@ref sec-theory-kinetics-admissible).

```@autodocs
Modules = [ChemistryLab]
Pages   = ["kinetics/dissipation.jl"]
```

## The implicit kinetic step

One fully implicit problem per step, with the reaction extents as unknowns of the
same Gibbs minimization. See
[Which route: two ways to advance in time](@ref sec-kinetics-routes).

```@autodocs
Modules = [ChemistryLab]
Pages   = ["kinetics/implicit_step.jl"]
```
