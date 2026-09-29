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

## Numerical floors

The small amounts the solvers and the activity models work with, each defined
once with the reason for its value: the default `ϵ` of the solvers, the floor of
the activities, the floor of the certificate and the exponent below which an
amount is taken as zero.

```@autodocs
Modules = [ChemistryLab]
Pages = ["utils/numerical_floors.jl"]
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

## Identifiability

Which parameters a measurement can actually determine, and which it only
appears to. See [Where the numbers come from](@ref sec-manual-numbers) for how
this closes onto [`Traced`](@ref): a parameter the data did not constrain comes
back as a placeholder rather than as a fitted value.

```@autodocs
Modules = [ChemistryLab]
Pages   = ["utils/identifiability.jl"]
```
