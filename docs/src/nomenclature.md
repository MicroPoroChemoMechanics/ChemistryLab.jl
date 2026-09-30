# [Nomenclature](@id nomenclature)

The symbols of the formulas of this documentation, grouped by subject. A symbol
that means different things in different chapters is listed once for each
meaning, with the pages where that meaning holds; a meaning listed without pages
holds everywhere else. The value of a physical constant is the one ChemistryLab
computes with, read from the library when this page is built.

Hovering an equation on any page shows the symbols it holds, with their meaning
on that page.

## Typography

A scalar, and a component of a vector or a matrix, is set in italic: ``n_i``,
``A_{ci}``, ``y_c``. A vector is set in bold upright lower case, ``\mathbf{n}``,
``\mathbf{b}``, ``\mathbf{y}``, a matrix in bold upright capitals,
``\mathbf{A}``, ``\mathbf{H}``, and a Greek vector or matrix in bold,
``\boldsymbol{\nu}``, ``\boldsymbol{\xi}``; the transpose is
``\mathbf{A}^\mathsf{T}``. So ``\mathbf{A}\mathbf{n} = \mathbf{b}`` reads, one row
at a time, ``\sum_i A_{ci}\, n_i = b_c``. A species written with a letter in a
general reaction is set upright, as a chemical symbol is:
``\sum_i \nu_i\,\mathrm{A}_i = 0``. Each quantity has one symbol throughout,
and a symbol has one meaning on a page.

```@eval
using ChemistryLab
include(joinpath(pkgdir(ChemistryLab), "docs", "nomenclature.jl"))
nomenclature_markdown()
```
