# Utilities

```@index
Pages = ["utils.md"]
```

## Module

```@autodocs
Modules = [ChemistryLab]
Pages = ["ChemistryLab.jl"]
```

## Physical constants

The gas constant, the Faraday constant and the electric constant, taken from
`DynamicQuantities` rather than written down, each in two forms: the bare SI
number for arithmetic and the dimensioned quantity for anything that must carry
its unit. See [Where the numbers come from](@ref sec-manual-numbers) for how to
reach these, and everything else the package already knows, from a script.

```@autodocs
Modules = [ChemistryLab]
Pages = ["utils/constants.jl"]
```

## Miscellaneous helpers

```@autodocs
Modules = [ChemistryLab]
Pages = ["utils/misc.jl"]
```

## Sub/superscripts

```@autodocs
Modules = [ChemistryLab]
Pages = ["utils/subsuperscripts.jl"]
```

## Provenance

A number that says how it was obtained and from where, so the claim survives out
of a data file and into a table, a figure or a fitted result. See
[Where the numbers come from](@ref sec-manual-numbers) for when to reach for it.

```@autodocs
Modules = [ChemistryLab]
Pages   = ["utils/provenance.jl"]
```
