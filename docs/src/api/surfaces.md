# Surfaces API

How much area a solid offers, and to what. One abstraction serves the reactive
area of a dissolution rate law and the binding area that carries surface sites:
the three areas of a problem — reactive, accessible, binding — are different
numbers, but they are the same kind of object, and the measurement each comes
from is carried by its type.

## Areas

```@autodocs
Modules = [ChemistryLab]
Pages   = ["surfaces/surface_areas.jl"]
```

## Site families

A site is a conserved quantity that is not a chemical element, and it is carried
the way electric charge already was: as a pseudo-element in a species' formula,
turned into a conservation row by the ordinary matrix assembly.

```@autodocs
Modules = [ChemistryLab]
Pages   = ["surfaces/site_families.jl"]
```

## The surface potential as an unknown

A [`DiffuseLayer`](@ref) written as an activity coefficient is exact and not
always reachable: the solver's inner loop recovers a mixing phase with `lnγ`
read at the previous iterate, which contracts only while
[`electrostatic_stiffness`](@ref) stays below
[`ELECTROSTATIC_STIFFNESS_LIMIT`](@ref). Carrying the potential as an unknown of
the outer Newton removes that fixed point rather than taming it.

```@autodocs
Modules = [ChemistryLab]
Pages   = ["surfaces/surface_potential.jl"]
```
