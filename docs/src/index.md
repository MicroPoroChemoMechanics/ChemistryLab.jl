```@raw html
---
# https://vitepress.dev/reference/default-theme-home-page
layout: home

hero:
  name: "ChemistryLab.jl"
  text: "Chemistry you can script"
  tagline: Formulas, species, reactions, thermodynamic databases and equilibrium — for aqueous geochemistry and cement chemistry.
  image:
    src: /logo.png
    alt: ChemistryLab
  actions:
    - theme: brand
      text: Get started
      link: /quickstart
    - theme: alt
      text: Theory
      link: /theory/
    - theme: alt
      text: Applications
      link: /examples/co2_carbonate_system
    - theme: alt
      text: View on GitHub
      link: https://github.com/MicroPoroChemoMechanics/ChemistryLab.jl

features:
  - icon: 🚀
    title: Getting started
    details: Install the package, evaluate the solubility constant of a reaction, then solve and certify a first equilibrium.
    link: /quickstart
  - icon: 🎓
    title: Theory
    details: Why each calculation is the right one, written for a reader who has never studied chemical thermodynamics, with every equation the one the code evaluates.
    link: /theory/
  - icon: 🧰
    title: Manual
    details: One kind of object per page — formulas, species, databases, stoichiometric matrices, systems and states, surfaces, cement notation.
    link: /manual/where_the_numbers_come_from
  - icon: 🧭
    title: Tutorials
    details: A calculation driven from start to finish — an equilibrium, a kinetic trajectory, their coupling, the self-desiccation of a paste.
    link: /tutorials/equilibrium
  - icon: 🧪
    title: Applications
    details: Complete cases and what the modeling choices cost in numbers — aqueous equilibria, surfaces, Portland and blended cements, hydration in time.
    link: /examples/co2_carbonate_system
  - icon: 📖
    title: API reference
    details: Every exported function and type, grouped by topic.
    link: /api/formulas
---
```

## What it does

`ChemistryLab` handles chemical formulas, species and reactions as first-class
objects, reads thermodynamic data from ThermoFun and Cemdata, and solves
equilibrium by Gibbs-energy minimization. Kinetics and chemo-mechanical coupling
build on the same objects.

It is written for work that has to be reproducible and scripted: aqueous
geochemistry, cement chemistry, and any problem where speciation, a database and
a solver have to be driven from code rather than from a dialog box.

## A first calculation

Calcite in water. The species come from a thermodynamic database, four of them
are declared as the independent basis, and the equilibrium state follows from a
Gibbs-energy minimization:

```@example home
using ChemistryLab, DynamicQuantities
using OptimaSolver          # the default equilibrium back-end

species = speciation(build_species(datapath("cemdata18-thermofun.json"); verbose = false),
                     split("Cal H2O@ CO2");
                     aggregate_state = [AS_AQUEOUS],
                     exclude_species = split("H2@ O2@ CH4@"))
byname = Dict(symbol(s) => s for s in species)

system = ChemicalSystem(collect(values(byname)),
                        [byname[s] for s in split("H2O@ H+ CO3-2 Ca+2")])

state = ChemicalState(system)
set_quantity!(state, "Cal", 1e-3u"mol")     # 1 mmol of calcite
set_quantity!(state, "H2O@", 1.0u"kg")      # in 1 kg of water
V = volume(state)
set_quantity!(state, "H+",  1e-4u"mol/L" * V.liquid)
set_quantity!(state, "OH-", 1e-10u"mol/L" * V.liquid)

equilibrated = equilibrate(state)
nothing # hide
```

The state carries everything derived from it — pH, phase volumes, individual
amounts:

```@example home
using Printf
@printf("pH                = %.2f\n", pH(equilibrated))
@printf("dissolved Ca(2+)  = %.3e mol\n", ustrip(moles(equilibrated, "Ca+2")))
@printf("remaining calcite = %.3e mol\n", ustrip(moles(equilibrated, "Cal")))
```

[Getting started](@ref sec-quickstart) takes the same problem more slowly, and
shows how the solubility constant is obtained analytically before any solver is
involved.

## Reading paths

The documentation is organized by the question each chapter answers.
[Theory](@ref sec-theory) explains why a calculation is the right one, the
Manual describes how each kind of object is written, the Tutorials drive a
calculation from start to finish, and the Applications show complete cases and
what the modeling choices cost in numbers. Three entry points follow from it,
depending on what the reader already knows.

**New to chemical thermodynamics.** The chapter [Theory](@ref sec-theory) is
written for this reader and is best read in the order it states, beginning with
[Thermochemistry](@ref sec-theory-thermo) and
[Standard states](@ref sec-theory-standard-states), which fix the notation and
the conventions every other page relies on. The Manual then introduces the
objects one at a time, from [Species](@ref sec-species) onwards, and the tutorial
[Chemical Equilibrium](@ref sec-equilibrium) puts them together.

**Knowing what is to be computed.** [Getting started](@ref sec-quickstart) and
the tutorial [Chemical Equilibrium](@ref sec-equilibrium) are enough to write a
first calculation. The aqueous applications, beginning with
[CO₂ dissolution and the carbonate system](@ref sec-co2-carbonate), are small
enough to check by hand and are the place to test one's understanding before a
larger system.

**Holding the analysis of a cement.** The route from an oxide analysis to a
chemical system is laid out in [Bogue calculation](@ref ex-bogue) and
[Choosing the species list](@ref man-choosing-species); the binders of EN 197-1
are then worked in increasing order of difficulty, from
[A CEM I from its clinker phases](@ref sec-cem1-from-clinker) to the composite
cements, and [The binders, and what distinguishes them](@ref man-binder-families)
says which of the package's models each family requires.

## Where this comes from

This package would not exist without two bodies of work that came before it, and
it is worth saying so on the first page rather than in a footnote.

**GEM-Selektor and GEMS3K**, from the Paul Scherrer Institute and Empa
[Kulik2013](@cite), established the Gibbs energy minimization approach used here,
and much of the vocabulary with it — the phase stability index this package
computes as ``\Omega`` is the same quantity as their ``\Lambda_k``, reached from
the same KKT conditions. **CEMDATA18** [Lothenbach2019](@cite), the thermodynamic
database behind every cement calculation in this manual, is their laboratory's
work and ships here unchanged; the zeolite extension is transcribed from two
further papers by the same group. Nothing here would produce a number without it.

**Reaktoro** [Leal2017](@cite), by Allan Leal and contributors, is both an
ancestor and a reference: parts of the thermodynamics and kinetics here are Julia
ports of its C++ implementation, and it is the code this package checks itself
against throughout — see [the validation page](@ref Validation-against-Reaktoro).
Where a result differs, the burden of proof has been on us.

Both are mature, carefully built and widely used, and both address a wider range
of problems than this package attempts. What ChemistryLab tries to add is
narrower: a Julia-native formulation in which an equilibrium comes with a
**proof** of its optimality rather than a converged iterate, differentiable end
to end so that a calibration can be posed as an optimization, and with the
cementitious special cases — cement chemist notation, Bogue, the oxide-budget
entry route for a glass, the binder families of EN 197-1 — treated as first-class
rather than as an application layer.

That is an addition to their work, not a substitute for it.
