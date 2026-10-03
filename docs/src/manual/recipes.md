# [Recipes: materials, extents and what has not reacted](@id man-recipes)

!!! info "Before this page"
    [The binders, and what distinguishes them](@ref man-binder-families) and
    [Solving an equilibrium](@ref sec-solving).

A cement calculation starts from **materials** (a clinker, a slag, a limestone,
water), their **masses**, and **how much of each has reacted** at the age of
interest. What has reacted enters the equilibrium; what has not stays in the
paste as a solid that takes up volume and holds matter, without taking part in
the chemistry. Written by hand, this bookkeeping is repeated on every page and
each copy can drift from the others. This page describes the layer that does it
once: [`Material`](@ref), [`Recipe`](@ref), [`budget`](@ref) and the processes
built on them.

## Materials

A material is a list of **constituents**, each with its mass fraction of the
material and its own degree of reaction:

  - a [`MineralConstituent`](@ref) is a database phase (alite, calcite, quartz);
  - an [`OxideConstituent`](@ref) is known by its oxide analysis only, like the
    glass of a slag, and enters the equilibrium through its elements.

Published materials are built from `data/recipe_templates.toml`, which names
the literature file each number is read from. The system below is a small paste:
the clinker phases, the usual hydrates and CSHQ.

```@example recipes
using ChemistryLab, OptimaSolver, DynamicQuantities
substances = build_species(datapath("cemdata18-thermofun.json"); verbose = false)
db = Dict(symbol(s) => s for s in substances)
pure = split("C3S C2S C3A C4AF Gp Anh Cal Portlandite ettringite monosulphate12 " *
             "monocarbonate hemicarbonate C4AH13 C3AH6 C3FH6 straetlingite " *
             "hydrotalcite Mg2AlC0.5OH Brc FeOOHmic AlOHmic Amor-Sl K2SO4")
csh = ["CSHQ-JenD", "CSHQ-JenH", "CSHQ-TobD", "CSHQ-TobH", "KSiOH", "NaSiOH"]
species = speciation(substances, vcat(pure, csh); aggregate_state = [AS_AQUEOUS],
                     exclude_species = split("H2@ O2@ CH4@"))
cs = ChemicalSystem(species, CEMDATA_PRIMARIES;
                    solid_solutions = [SolidSolutionPhase("CSHQ", [db[m] for m in csh])])
material_templates()
```

This list is written out by hand. The phases a paste may form are a modeling
choice, and `data/phase_lists.toml` records that choice for published pastes,
each list after the paper it names. [`phase_list_system`](@ref) builds a system
from one, and [the processes page](@ref ex-cement-processes) uses it.

The Portland cement of the RILEM round robin by Bogue's calculation on its
oxide analysis ([`bogue`](@ref)): four clinker phases and gypsum, and the oxides
no phase takes as one oxide constituent.

```@example recipes
pc = material_template("PC (Durdzinski 2017), Bogue", cs)
[(c.name, round(c.mass_fraction; digits = 3)) for c in pc.constituents]
```

## How much has reacted

An extent is a number between 0 and 1, a function of time, or one of
[`ConstantExtent`](@ref), [`TabulatedExtent`](@ref) (measured degrees of
reaction, interpolated), [`LogisticExtent`](@ref), [`ParrottKillohExtent`](@ref)
(the rate law of Parrott and Killoh integrated) and [`CappedExtent`](@ref)
(any of them under Powers' water limit). [`with_extents`](@ref) sets them by
constituent, and `material_extent` sets one for the material as a whole; the
effective extent of a constituent is the product of the two
([`effective_extent`](@ref)). A glass described by its oxides alone has one
extent, the material's.

```@example recipes
pc = with_extents(pc, Dict("C3S" => 0.80, "C2S" => 0.45, "C3A" => 0.90, "C4AF" => 0.60))
slag = with_extents(material_template("S1 slag (Durdzinski 2017)", cs), Dict{String, Any}();
                    material_extent = 0.40)
nothing # hide
```

## A recipe, and what it puts into the equilibrium

A [`Recipe`](@ref) gives the binder's mass fractions and the water/binder ratio.
Every amount is for 100 g of binder unless `binder_mass` says otherwise, so a
mass read off the answer is in grams per 100 g of binder.

```@example recipes
recipe = Recipe(pc => 0.7, slag => 0.3; w_b = 0.45)
bud = budget(recipe, cs)
[(x.constituent, round(x.mass; digits = 2), x.reason) for x in bud.residual]
```

The rules [`budget`](@ref) applies:

| what | where it goes |
|:--|:--|
| the reacted part of a mineral constituent | the equilibrium, as that species |
| the reacted part of an oxide constituent | the equilibrium, as its elements ([`oxide_budget`](@ref)) |
| water, salts | the equilibrium, always completely |
| the unreacted part of a constituent | the residue: its mass, its volume, its enthalpy |
| an oxide whose element the system has no primary for | the residue, with `reason = :not_in_system` |

The residue keeps its volume from the molar volume of its species or from the
density of an oxide constituent, and its enthalpy from the enthalpy of
formation. A glass has no formula, so no database gives either: unless its
source does,
they are `nothing`, and a volume or a heat that needs them is reported missing
(`NaN`), never estimated.

## The answer

[`equilibrate_certified`](@ref), given a recipe and a system, solves the budget and returns a [`RecipeState`](@ref).

```@example recipes
model = cemdata18_activity_model(:KOH)
rs, cert = equilibrate_certified(recipe, cs; model)
cert.optimal
```

What is read off it counts the residue where it belongs.
[`phase_masses`](@ref) lists the solids, the unreacted parts last:

```@example recipes
first(phase_masses(rs), 6)
```

[`porosity`](@ref) is referred to the initial volume, that of the reactants, the
water and the residue:

```math
\phi = \frac{V_\text{liquid} + V_\text{gas} + V_\text{void}}{V_0},
\qquad V_\text{void} = V_0 - (V_\text{equilibrium} + V_\text{residue}) ,
```

the void being the chemical shrinkage the reactions did not refill. With a slag
glass of unknown density, the residue's volume is missing and the porosity
cannot be computed: `NaN`, which is the honest answer. [`volume_fractions`](@ref)
of a recipe state divides the same initial volume among the phases, the
unreacted parts and the void, and is refused for the same reason. [`bound_water`](@ref) is
per gram of binder and [`pore_solution`](@ref) gives the pH and the dissolved
elements, in mol per kg of water.

```@example recipes
bound_water(rs), pore_solution(rs).pH
```

## Bound water and heat

[`bound_water`](@ref) is the water the solids would lose on ignition
([`ignition_loss`](@ref)), with that of an unreacted mineral such as gypsum,
over the binder mass. It is the whole of what a thermogram integrates to. A
thermogravimetric reading taken between 105 °C and 550 °C, a common choice,
leaves part of it out, so the two are compared only over the same temperature
window: `bound_water(rs; window = (T₁, T₂), windows)` counts what the solids
release between the two temperatures, from a [`DecompositionWindow`](@ref) for
each solid that holds water, and refuses a solid without one. [The thermogram
page](@ref sec-example-tga) shows where such windows come from.

At constant temperature and pressure, the heat a paste releases between two
states is the fall of its enthalpy:

```math
Q = -\left[H_2 - H_1\right] .
```

[`heat_release`](@ref), given two states of one paste, computes it. The residue
counts by what changed. The unreacted part of a clinker phase has the enthalpy
of formation of its database record, so its change is counted. An oxide the
system cannot hold, such as the titanium of this cement, keeps the same mass in
both states and adds nothing, although no enthalpy is known for it. Where the
unreacted part of a constituent without a sourced enthalpy changes, the heat is
`NaN` rather than a number that leaves it out. The slag glass of this recipe is
such a constituent, and [`enthalpy`](@ref) of the whole paste is `NaN` for the
same reason.

Below, the cement alone at two sets of extents, both chosen for the
illustration, gives the heat released from the first to the second, in joules
per gram of cement:

```@example recipes
early = with_extents(pc, Dict("C3S" => 0.40, "C2S" => 0.10, "C3A" => 0.50, "C4AF" => 0.30))
rs1, _ = equilibrate_certified(Recipe(early => 1.0; w_b = 0.45), cs; model)
rs2, _ = equilibrate_certified(Recipe(pc => 1.0; w_b = 0.45), cs; model)
(slag_paste = enthalpy(rs), heat = heat_release(rs1, rs2) / 100)
```

## Processes

A process is a sequence of equilibria, each started from the previous answer,
returned as a [`ProcessResult`](@ref) that [`process_table`](@ref) turns into a
table:

  - [`hydrate`](@ref): the recipe at several ages, its extents read at each;
  - [`blend`](@ref): part of the binder replaced by a material, at several
    fractions;
  - [`titrate`](@ref), [`carbonate`](@ref), [`add_salt`](@ref): a species added
    to the budget in steps;
  - [`leach`](@ref): the pore solution replaced by water, step after step.

```@example recipes
co = carbonate(rs, [0.0, 0.2, 0.5])
process_table(co; phases = ["Portlandite", "Cal"])
```

## From a recipe to a kinetic run

[`hydrate`](@ref) imposes the extents. Rate laws can decide them instead: a
[`KineticsProblem`](@ref) built from the recipe takes each constituent named in
`rates` whole and unreacted, as a kinetic species that dissolves into the
primaries of the system, and every other constituent as the recipe has it at
the start of the run. The same recipe then drives both.

```@example recipes
rates = Dict(
    "C3S" => parrott_killoh_avrami(PK84_PARAMS_C3S, "C3S"),
    "C2S" => parrott_killoh_avrami(PK84_PARAMS_C2S, "C2S"),
)
kp = KineticsProblem(Recipe(pc => 1.0; w_b = 0.45), cs, rates, (0.0, 28 * 86400.0))
for kr in kp.kinetic_reactions
    println(kr.reaction)
end
```

The dissolutions are those the conservation matrix of the system gives, so they
conserve every element by construction. The reacted part of a constituent known
by its oxides (the alkalis and the minor oxides of a clinker) enters as the
primaries that carry its elements; its heat of dissolution is not part of the
heat of the run. A glass has no formula to dissolve and cannot be given a rate.
[`integrate`](@ref) then runs the problem, with an equilibrium solver for the
partial equilibrium of the products, as in the [kinetics tutorial](@ref
sec-kinetics).

## Where to go next

The rate laws a recipe attaches to its constituents are written as described in
[Writing a kinetic model](@ref sec-kinetics-syntax), and the run itself is the one
of the tutorial [Chemical Kinetics](@ref sec-kinetics). The blended binders use
recipes throughout, from [CEM II/A-LL and CEM II/B-S](@ref ex-cem2-blended)
onwards.
