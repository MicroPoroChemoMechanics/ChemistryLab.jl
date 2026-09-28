# [Nomenclature](@id nomenclature)

The symbols of the formulas of this documentation, grouped by subject. A symbol
that means different things in different chapters is listed once for each
meaning, with the pages where that meaning holds; a meaning listed without pages
holds everywhere else. The value of a physical constant is the one ChemistryLab
computes with, read from the library when this page is built.

Hovering an equation on any page shows the symbols it holds, with their meaning
on that page.

```@eval
using ChemistryLab
include(joinpath(pkgdir(ChemistryLab), "docs", "nomenclature.jl"))
nomenclature_markdown()
```
