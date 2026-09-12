# [The binders, and what distinguishes them](@id man-binder-families)

A "cement" is not one material. EN 197-1 recognizes five main types and some
twenty-seven products, and the difference between them is **how much of the
clinker has been replaced, and by what**. That single question decides the
chemistry: which elements enter the paste, which hydrates can form, and which of
this package's models the calculation needs.

This page is the map. The calculations themselves are on the pages it points to.

## The families

| type | name | clinker | main other constituent |
|:--|:--|--:|:--|
| **CEM I** | Portland | 95–100 % | — |
| **CEM II/A-S** | Portland-slag | 80–94 % | blastfurnace slag 6–20 % |
| **CEM II/B-S** | Portland-slag | 65–79 % | blastfurnace slag 21–35 % |
| **CEM II/A-D** | Portland-silica fume | 90–94 % | silica fume 6–10 % |
| **CEM II/A-P, /A-Q** | Portland-pozzolana | 80–94 % | natural (P) or calcined (Q) pozzolana 6–20 % |
| **CEM II/B-P, /B-Q** | Portland-pozzolana | 65–79 % | pozzolana 21–35 % |
| **CEM II/A-V, /A-W** | Portland-fly ash | 80–94 % | siliceous (V) or calcareous (W) fly ash 6–20 % |
| **CEM II/B-V, /B-W** | Portland-fly ash | 65–79 % | fly ash 21–35 % |
| **CEM II/A-T, /B-T** | Portland-burnt shale | 80–94 / 65–79 % | burnt shale 6–20 / 21–35 % |
| **CEM II/A-L, /A-LL** | Portland-limestone | 80–94 % | limestone 6–20 % |
| **CEM II/B-L, /B-LL** | Portland-limestone | 65–79 % | limestone 21–35 % |
| **CEM II/A-M, /B-M** | Portland-composite | 80–94 / 65–79 % | a mixture of the above |
| **CEM III/A** | Blastfurnace | 35–64 % | slag 36–65 % |
| **CEM III/B** | Blastfurnace | 20–34 % | slag 66–80 % |
| **CEM III/C** | Blastfurnace | 5–19 % | slag 81–95 % |
| **CEM IV/A** | Pozzolanic | 65–89 % | D + P + Q + V + W, 11–35 % |
| **CEM IV/B** | Pozzolanic | 45–64 % | D + P + Q + V + W, 36–55 % |
| **CEM V/A** | Composite | 40–64 % | slag 18–30 % **and** pozzolana/fly ash 18–30 % |
| **CEM V/B** | Composite | 20–38 % | slag 31–49 % **and** pozzolana/fly ash 31–49 % |

!!! warning "The standard is the authority, not this table"
    The ranges above summarize **EN 197-1**, *Cement — Part 1: Composition,
    specifications and conformity criteria for common cements*. They are given
    as context for choosing a composition to calculate; percentages are by mass
    of the main constituents, excluding calcium sulfate and minor additional
    constituents. Anything with consequences — a conformity claim, a
    specification — must be read from the standard itself, which this package
    does not ship and cannot substitute for.

    **EN 197-5** adds two families this table does not carry: CEM II/C-M and
    CEM VI, which push the replacement further than EN 197-1 allows.

## What each constituent brings to the chemistry

The percentages matter because of what they *do*. Each constituent changes the
element budget, and the element budget decides which hydrates the minimization
can form.

| constituent | brings | consequence for the calculation |
|:--|:--|:--|
| **clinker** | Ca, Si, Al, Fe, and all of the heat | the reference case; see [CEM I from the clinker up](@ref) |
| **blastfurnace slag (S)** | Ca, Si, Al, **Mg**, and sulfur as **S(-II)** | needs a [redox](@ref theory-redox) treatment: the slag's sulfide meets the pore solution's sulfate, and both must be held at once. The Mg forms hydrotalcite |
| **fly ash (V, W)** | Si, Al, alkalis; W also Ca | the Al goes into the C-S-H, which then needs `CNASH_ss` rather than `CSHQ` |
| **natural/calcined pozzolana (P, Q)** | Si, Al, alkalis | as fly ash, and at high alkali the zeolites become stable |
| **silica fume (D)** | Si, and nothing else | lowers the Ca/Si of the C-S-H; no new phase family |
| **limestone (L, LL)** | **CO₃** | changes the aluminate sequence: monocarboaluminate forms instead of monosulphate, which stabilizes the ettringite |
| **burnt shale (T)** | Ca, Si, Al, sulfate | behaves as a weak clinker plus a pozzolana |

Read the table as a list of **model requirements**, because that is what it is.
A CEM I needs none of this package's later machinery; a CEM III needs the
oxidation state to be a conserved quantity; a CEM IV or V needs a C-S-H that can
take aluminum and alkalis.

## What the heat says about all of it

Every joule of hydration heat comes from the clinker. Replacing clinker with
anything lowers it, and the ordering of the measured records this package ships
is exactly the ordering of the replacement:

| record | family | final *Q* |
|:--|:--|--:|
| `122-cemI-52.5R-cizkovice` | CEM I | **376 J/g** |
| `116-cemI-52.5R-ladce` | CEM I | 355 J/g |
| `165-cemII-A-LL-42.5R-hranice` | CEM II/A-LL | 329 J/g |
| `149-cemII-B-S-32.5R-mokra` | CEM II/B-S | 296 J/g |
| `184-cemIII-A-42.5N-hranice` | CEM III/A | 261 J/g |
| `200-cemV-A-S-V-32.5R-prachovice` | CEM V-A (S-V) | 259 J/g |
| `121-cemIII-B-32.5N-mokra` | CEM III/B | **234 J/g** |

Measured on one instrument at 20 °C, from the CC-BY-4.0 deposit of Šmilauer and
Reiterman [Smilauer2025data](@cite); see `data/experimental/README.md` for the
full provenance and for one inconsistency found in the source metadata.

!!! danger "What these records do **not** report"
    The deposit gives the calorimetry, the Blaine fineness and the water/binder
    ratio. It gives **neither the clinker phase composition nor the actual
    replacement level** of any blend.

    So a calculation of one of these cements has to *assume* a composition
    inside the EN 197-1 range of its family. Every page that does so says which
    number is measured and which is assumed, at the point of use. A page that
    presented "40 % clinker" as a property of the specimen would be inventing
    it, and the distinction matters most exactly where a reader is least able to
    check.

## Where the calculations are

| page | what it does |
|:--|:--|
| [Bogue calculation](@ref) | an oxide analysis to clinker phases — the entry point for a CEM I |
| [CEM I from the clinker up](@ref) | the reference paste, hydrated and certified |
| [A CEM I at equilibrium, with every solid solution declared](@ref) | the phase list as a modeling decision |
| [The full Portland cement, through its pore solution](@ref ex-ionic-opc) | the coupled run, and its calorimetry |
| [Calibrating hydration kinetics](@ref) | the inverse problem, against measured calorimetry |

And one page per blended family, each on the same shape — element budget in,
certified assemblage out, measured calorimetry beside it:

| page | family | what it is about |
|:--|:--|:--|
| [Two CEM II, and the two different things a replacement can do](@ref ex-cem2-blended) | CEM II | a carbonate that rewrites the aluminate sequence, against a glass that brings magnesium |
| [A blastfurnace cement, and the oxidation state it needs](@ref ex-cem3-slag) | CEM III | the glass entry route, hydrotalcite, and sulfur at two oxidation states |
| [A pozzolanic binder, and the C-S-H that has to carry the aluminum](@ref ex-cem4-pozzolanic) | CEM IV | `CSHQ` against `CNASH_ss`, and portlandite as the limiting reagent |
| [A composite binder: two glasses at once](@ref ex-cem5-composite) | CEM V | all four difficulties simultaneously, on one additive budget |

The CEM IV page is the one without a measured specimen behind it — the deposit
carries no record of that family — and it says so at its head rather than in a
footnote.
