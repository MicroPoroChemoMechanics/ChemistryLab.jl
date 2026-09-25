# [Cement chemist notation](@id man-cement-notation)

!!! info "Before this page"
    [Chemical Formula Manipulation](@ref sec-formulas).

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

```@example cemnot
using ChemistryLab
using DynamicQuantities
using Markdown, Printf

# The names are this page's. The formulas are the parser's own mapping,
# CEMENT_TO_MENDELEEV, and the molar masses are computed from them.
oxide_names = [
    :C => "lime", :S => "silica", :A => "alumina", :F => "ferric oxide",
    :M => "magnesia", :K => "potassium oxide", :N => "sodium oxide",
    :T => "titania", :P => "phosphorus pentoxide", :H => "water",
    :C̄ => "carbon dioxide", :S̄ => "sulfur trioxide", :N̄ => "nitrate",
]
mapping = Dict(CEMENT_TO_MENDELEEV)
@assert Set(first.(oxide_names)) == Set(keys(mapping))   # every letter, and no other

rows = ["| Symbol | Oxide | Name | Molar mass |", "|:------:|:------|:-----|-----------:|"]
for (letter, oxide_name) in oxide_names
    formula = join(string(el, n == 1 ? "" : n) for (el, n) in mapping[letter])
    M = ustrip(us"g/mol", Species(formula).M)
    push!(rows, @sprintf("| `%s` | %s | %s | %.2f g/mol |",
                         letter, phreeqc_to_unicode(formula), oxide_name, M))
end
Markdown.parse(join(rows, "\n"))
```

The table is produced from the mapping the parser itself uses,
[`CEMENT_TO_MENDELEEV`](@ref), and its molar masses are computed from the
formulas, so it cannot drift from what ChemistryLab accepts.

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
- the [Bogue Calculation](@ref ex-bogue), which is that decomposition inverted and
  converted to mass.

## Where to go next

The notation is attached to species in [Cement Species](@ref man-cement-species),
and its inversion from an oxide analysis to phase masses is
[Bogue Calculation](@ref ex-bogue).
