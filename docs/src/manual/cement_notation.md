# [Cement chemist notation](@id man-cement-notation)

Cement chemistry writes its phases in a shorthand that is opaque on first sight
and indispensable once learned: alite is `C3S`, portlandite is `CH`, ettringite
is `C6AS̄3H32`. This page is the key to it, and to how ChemistryLab reads it.

The idea is simple. A cement phase is an assembly of **oxides**, so each oxide
gets a one-letter symbol and a phase is written as the oxides it contains, with
their molar multiplicities. `C3S` is three parts lime to one of silica —
``3\,\ce{CaO}\cdot\ce{SiO2}`` — which is ``\ce{Ca3SiO5}`` written the way a
cement chemist thinks about it: not as an arbitrary silicate, but as *this much
lime combined with that much silica*.

## The alphabet

| Symbol | Oxide | Name | Molar mass |
|:------:|:------|:-----|-----------:|
| `C`  | ``\ce{CaO}``    | lime                  |  56.08 g/mol |
| `S`  | ``\ce{SiO2}``   | silica                |  60.08 g/mol |
| `A`  | ``\ce{Al2O3}``  | alumina               | 101.96 g/mol |
| `F`  | ``\ce{Fe2O3}``  | ferric oxide          | 159.69 g/mol |
| `M`  | ``\ce{MgO}``    | magnesia              |  40.30 g/mol |
| `K`  | ``\ce{K2O}``    | potassium oxide       |  94.20 g/mol |
| `N`  | ``\ce{Na2O}``   | sodium oxide          |  61.98 g/mol |
| `T`  | ``\ce{TiO2}``   | titania               |  79.86 g/mol |
| `P`  | ``\ce{P2O5}``   | phosphorus pentoxide  | 141.94 g/mol |
| `H`  | ``\ce{H2O}``    | water                 |  18.01 g/mol |
| `C̄`  | ``\ce{CO2}``    | carbon dioxide        |  44.01 g/mol |
| `S̄`  | ``\ce{SO3}``    | sulfur trioxide       |  80.06 g/mol |
| `N̄`  | ``\ce{NO3}``    | nitrate               |  62.00 g/mol |

The first four are the oxides a clinker is made of, and the four a datasheet
always reports; the next five are the minor oxides; the last three are the
acidic ones that take a bar.

Two conventions decide everything, and they are where a newcomer stumbles:

- **A letter stands for an oxide, not an element.** `C` is lime ``\ce{CaO}``, not
  carbon; `S` is silica ``\ce{SiO2}``, not sulfur; `N` is ``\ce{Na2O}``, not
  nitrogen.
- **The bar marks the acidic oxides** whose letter is already taken. Carbon
  dioxide is `C̄` because `C` is lime, and sulfur trioxide is `S̄` because `S` is
  silica. The bar is a combining macron (`U+0304`): type the letter, then that
  character.

The table above is the mapping the parser itself uses,
[`CEMENT_TO_MENDELEEV`](@ref). The block below prints it from the code and
checks the two agree, so the table cannot quietly drift from what ChemistryLab
accepts:

```@example cemnot
using ChemistryLab
using DynamicQuantities

# (formula, molar mass) exactly as the table above states them
documented = Dict(
    :C => ("CaO", 56.08), :S => ("SiO2", 60.08), :A => ("Al2O3", 101.96),
    :F => ("Fe2O3", 159.69), :M => ("MgO", 40.30), :K => ("K2O", 94.20),
    :N => ("Na2O", 61.98), :T => ("TiO2", 79.86), :P => ("P2O5", 141.94),
    :H => ("H2O", 18.01), :C̄ => ("CO2", 44.01), :S̄ => ("SO3", 80.06),
    :N̄ => ("NO3", 62.00),
)

for (letter, oxide) in CEMENT_TO_MENDELEEV
    sp = Species(oxide)
    formula, mass = documented[letter]
    @assert atoms(sp) == atoms(Species(formula)) "table disagrees on $letter"
    computed = round(ustrip(us"g/mol", sp.M), digits = 2)
    @assert abs(computed - mass) < 0.005 "table gives $mass for $letter, the data give $computed"
    println(rpad(string(letter), 3), " = ", rpad(unicode(sp), 8),
            lpad(computed, 8), " g/mol")
end
```

`H` is water, which is why hydrates carry a large `H` count: `C4AH13` is
``4\,\ce{CaO}\cdot\ce{Al2O3}\cdot 13\,\ce{H2O}``.

## Reading and building a phase

[`CemSpecies`](@ref) parses the shorthand. Ordinary ASCII digits and Unicode
subscripts are both accepted, and non-integer multiplicities are allowed — which
matters for C-S-H, whose lime-to-silica ratio is a composition, not a constant:

```@example cemnot
for name in ["C3S", "C2S", "C3A", "C4AF", "CH", "C₄AH₁₃", "C1.7SH4"]
    sp = CemSpecies(name)
    println(rpad(name, 9), " -> ", rpad(unicode(sp), 12),
            "  M = ", round(ustrip(us"g/mol", sp.M), digits = 2), " g/mol")
end
```

The composition in ordinary elements is always available, so the shorthand is a
way of *writing* a species and never a different kind of object:

```@example cemnot
atoms(CemSpecies("C3S"))
```

## [The trap: `Species` and `CemSpecies` read the same string differently](@id man-cement-notation-trap)

This is the single most expensive mistake to make, and it is silent:

```@example cemnot
cem = CemSpecies("C3S")
ord = Species("C3S")

println("CemSpecies(\"C3S\") = ", unicode(cem), "   M = ",
        round(ustrip(us"g/mol", cem.M), digits = 2), " g/mol   (3 CaO + SiO2)")
println("Species(\"C3S\")    = ", unicode(ord), "   M = ",
        round(ustrip(us"g/mol", ord.M), digits = 2), " g/mol   (3 carbons + 1 sulfur)")
```

Both **display as `C₃S`**, and they differ by more than a factor of three in
molar mass. [`Species`](@ref) applies the ordinary rules of chemical formulas, in
which `C` is carbon and `S` is sulfur; [`CemSpecies`](@ref) applies the cement
convention. Neither is wrong — they answer different questions — but a recipe
that reaches for the wrong one is wrong everywhere downstream and raises no
error.

!!! warning "Read a database phase by name, never by re-parsing its symbol"
    A phase read from a thermodynamic database already carries its composition
    and its molar mass. Look it up by name — `byname["C3S"]` — rather than
    rebuilding it from its symbol with `Species`, which would read `C3S` as
    three carbons and a sulfur. This is the reason the worked examples in this
    documentation never write a molar mass by hand.

## The common phases, and what they weigh

Nothing below is typed from a table: each mass is computed from the formula and
the element data.

```@example cemnot
phases = [
    ("C3S",       "alite"),
    ("C2S",       "belite"),
    ("C3A",       "aluminate"),
    ("C4AF",      "ferrite"),
    ("CS̄H2",      "gypsum"),
    ("CH",        "portlandite"),
    ("CC̄",        "calcite"),
    ("C6AS̄3H32",  "ettringite (AFt)"),
    ("C4AS̄H12",   "monosulfoaluminate (AFm)"),
    ("C4AH13",    "hydroxy-AFm"),
    ("C2ASH8",    "strätlingite"),
    ("C3AH6",     "hydrogarnet"),
    ("M5AH13",    "hydrotalcite"),
]

for (name, english) in phases
    sp = CemSpecies(name)
    println(rpad(unicode(sp), 13), rpad(english, 28),
            lpad(round(ustrip(us"g/mol", sp.M), digits = 2), 8), " g/mol")
end
```

!!! note "Ettringite carries alumina, and the shorthand must say so"
    Ettringite is ``\ce{C6AS̄3H32}`` — six lime, **one alumina**, three sulfate,
    thirty-two water — at 1255 g/mol. Dropping the `A` gives `C6S̄3H32`, which
    parses without complaint, weighs 1153 g/mol and is not a cement phase. The
    shorthand is compact enough that an omission looks like a typo and behaves
    like a different substance.

## Where the notation is used

- [`CemSpecies`](@ref) for the species themselves;
- the oxide components of a stoichiometric decomposition — see
  [Stoichiometric matrices](@ref ex-stoich-matrix), where a clinker phase is
  expressed over `C`, `S`, `A`, `F`;
- the [Bogue calculation](@ref), which is that decomposition inverted and
  converted to mass.
