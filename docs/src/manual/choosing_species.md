# [Choosing the species list](@id man-choosing-species)

A Gibbs energy minimization answers the question *"of the phases I was told
about, which assemblage has the lowest energy?"*. It cannot form a phase it was
not given, and it will not tell you that the one it needed was missing. It will
report a poor element balance, or a pH of 7, or an assemblage that looks almost
right — and none of those says "you forgot a phase".

This page is the check-list that turns that silent failure into a decision you
make on purpose. **Every trap below was met while writing this manual**, and each
is stated with what it looked like when it happened.

## The rule underneath all of it

> Every element in the budget must have somewhere to go.

That is not a figure of speech. The constraint is ``\mathbf{A}\mathbf{n} =
\mathbf{b}``: if a component of ``\mathbf{b}`` is non-zero and no declared species
carries it in a form the paste can reach, the feasible set is empty or nearly so,
and the solve fails in ways that look like numerical trouble rather than like a
modeling omission.

So the first thing to do with a new formulation is not to solve it. It is to ask,
element by element, **which declared phase receives it**.

```@example choosing
using ChemistryLab
using DynamicQuantities
using Printf

substances = build_species(datapath("cemdata18-thermofun.json"); verbose = false)
byname = Dict(symbol(s) => s for s in substances)

"""Which declared crystalline species can carry a given component, and how many."""
function carriers(names, element)
    out = [n for n in names if haskey(atoms(byname[n]), Symbol(element))]
    return out
end

nothing # hide
```

## The check-list, by what the formulation contains

| the formulation brings | it must be able to form | why, and what happens otherwise |
|:--|:--|:--|
| **clinker alone** (CEM I) | `Portlandite`, a C-S-H family, `ettringite`, an AFm, a hydrogarnet, `FeOOHmic` | the reference case; anything missing here shows up immediately |
| **limestone** (CEM II/L, LL) | `Cal`, **`monocarbonate`, `hemicarbonate`** | without the carbonate AFm the carbonate has only calcite, and the aluminate sequence is wrong rather than merely incomplete |
| **slag** (CEM II/S, CEM III) | **`hydrotalcite`**, and the sulfur ladder | the magnesium has no calcium-aluminate host at all; and slag sulfur arrives as S(-II) while the pore solution carries S(+VI) |
| **fly ash, pozzolana** (CEM II/V, CEM IV) | **`CNASH_ss`**, **the siliceous hydrogarnet**, `straetlingite`, `AlOHmic` | the ash brings nearly as much Al as Si, and `CSHQ` has **no aluminum end-member at all** |
| **high replacement** (> ~35 %) | **zeolites** — [the extension](@ref sec-zeolites) | past that the alkalis exceed what the C-A-S-H and the sulfates can hold |
| **alkalis from any source** | `K2SO4`, `syngenite`, `Na2SO4`, and an alkali-bearing gel | otherwise they stay in solution and the pH comes out too high |

## Seven traps, each one met in practice

### 1. Declaring one end-member of a family and not the other

The worst of them, because the calculation runs and the answer is wrong.

`C3AFS0.84H4.32` is the **aluminum–iron** member of the siliceous hydrogarnet;
`C3AS0.84H4.32` is the **aluminum** one. Declaring the first and not the second
was enough to break a CEM V: of 0.228 mol of aluminum in the budget, only 0.095
found a phase, the element balance stopped at 2.7e-01, and the worst
supersaturation was *negative* — so nothing was asking to form. Nothing in that
output points at a missing species.

```@example choosing
family = ["C3AS0.84H4.32", "C3AFS0.84H4.32", "C3FS0.84H4.32"]
for n in family
    @printf("  %-18s %s\n", n, unicode(byname[n]))
end
```

Read a family by its composition, not by its name.

### 2. Declaring a gel model that cannot hold the element you are adding

`CSHQ` has four calcium-silicate end-members and two alkali ones — and **no
aluminum end-member**. On a pozzolanic binder, every atom of aluminum the ash
brings is then forced into the AFm/AFt phases, and the aluminum balance of the
paste is wrong by construction.

```@example choosing
CSHQ = ["CSHQ-JenD", "CSHQ-JenH", "CSHQ-TobD", "CSHQ-TobH", "KSiOH", "NaSiOH"]
CNASH = ["T2C-CNASHss", "T5C-CNASHss", "TobH-CNASHss",
         "5CA", "5CNA", "INFCA", "INFCN", "INFCNA"]
@printf("CSHQ     end-members carrying Al : %d of %d\n",
        length(carriers(CSHQ, "Al")), length(CSHQ))
@printf("CNASH_ss end-members carrying Al : %d of %d\n",
        length(carriers(CNASH, "Al")), length(CNASH))
```

### 3. Declaring two models of the same phase

`CSHQ`, `CNASH_ss` and the `ECSH` family are three *models of one gel*, not three
phases. Declaring two counts the same hydrate twice. `ChemicalSystem` refuses the
pair when their end-members share a composition — `KSiOH`, `ECSH1-KSH` and
`ECSH2-KSH` are all `((KOH)2.5SiO2H2O)0.2` — but **choose deliberately** rather
than relying on the refusal.

### 4. Declaring a phase whose element you did not put in the budget

The mirror image of trap 1, and it is not harmless. Declaring `KSiOH` and
`NaSiOH` on a paste with no potassium and no sodium leaves them pinned at the
solver's floor; their saturation indices then run to +7 and the certificate has
to decide what that means.

The expensive version of this trap is quieter, and it was met on these very
pages. A **Bogue calculation returns four phases and no alkalis** — Na₂O and K₂O
are minor oxides, outside the four-phase decomposition — so a binder entered
through Bogue has no sodium and no potassium in its budget however carefully its
C-S-H is declared. Nothing fails. Every solve certifies. What comes out is a pH
of **12.51 on every paste**, because portlandite is then the only thing setting
it, and a portlandite buffer is by construction insensitive to everything else.

!!! tip "The diagnostic is the invariance, not the value"
    Three CEM II pastes and a CEM III, differing in replacement level, in w/b and
    in assemblage, all returned **12.510 to three decimals**. A quantity that
    does not move when the inputs move is either buffered or not computed from
    them — here both. The corroboration was a fifth calculation: the CEM V, whose
    fly ash carries 2.5 % K₂O, was the only one with alkalis in its budget and
    the only one that did not return 12.510.

    A real cement pore solution sits above 13 for exactly this reason: the
    alkalis dissolve almost completely and stay in solution, while the calcium is
    held at the portlandite floor. If your pH comes out at 12.5 and will not
    move, look at the `K+` and `Na+` rows of the budget before looking at the
    solver.

Two of the twenty-eight zeolites are a chloride and a nitrate sodalite. On a
binder carrying neither element, declaring them pulls every aqueous chloride and
nitrate species in the database into the system on a budget of exactly zero —
a larger, slower problem for no phase that can form.

```@example choosing
zeo = build_species(datapath("cemdata18-zeolites.json"); verbose = false)
zeo_by = Dict(symbol(s) => s for s in zeo)
added = sort(collect(setdiff(Set(keys(zeo_by)), Set(keys(byname)))))
carries(n, el) = haskey(atoms(zeo_by[n]), Symbol(el))
drop = [n for n in added if carries(n, "Cl") || carries(n, "N")]
@printf("%d zeolites added, %d of them carry Cl or N: %s\n",
        length(added), length(drop), join(drop, ", "))
```

### 5. Forgetting that a glass has no formula

Slag, fly ash and natural pozzolana have no formula unit, so they cannot be
entered as a species. They enter as an **oxide analysis** converted to component
totals by [`oxide_budget`](@ref), added to the clinker's budget. See
[the CEM III page](@ref ex-cem3-slag).

And a budget says what a material *contains*, never what it does: a slag and a
quartz sand of the same analysis give the same ``\mathbf{b}``.

### 6. Reading a failed solve as a numerical problem

The failure modes of a missing phase, as they actually appeared:

| symptom | what it usually means |
|:--|:--|
| element balance stuck at 1e-1 with **negative** worst supersaturation | an element has nowhere to go — nothing is asking to form because nothing *can* |
| pH exactly 6.999 | the solve failed and returned neutral water |
| pH exactly 12.51, and the same on unrelated pastes | portlandite is the only buffer, because the budget carries no alkalis |
| an assemblage still containing anhydrous clinker | the solve never reached hydration |
| a sweep whose pH jumps around non-monotonically | isolated failed solves inside an otherwise fine scan |
| **one** point failing between two that certify | a starting point, not an infeasibility — walk to it by continuation from its neighbor |
| a saturation index of +7 on a phase at 1e-305 mol | a declared phase whose element is absent from the budget |

Most of these say "add a phase", and all of them did mean exactly that — except
the last, which says the opposite and is worth separating out. A configuration
that has **no** admissible assemblage fails everywhere near itself; a point that
fails while both of its neighbors certify has an answer the search did not
reach. The remedy there is not a species but a **start**: solve the easy neighbor
first and continue from it. On a convex problem the minimum is unique, so a
continuation cannot change what is found — only whether it is found — and the
certificate still decides every point. Both blended-binder sweeps in this
documentation are written that way, and both had a point that needed it.

### 7. Declaring a solid solution whose range cannot reach where the answer is

Trap 1's family has a second edge to it, and it is sharper. The siliceous
hydrogarnet is a substitution on **two sites**, so its three members are not
three points on one axis but

| member | occupancy | ``x(\mathrm{Al})`` |
|:--|:--|--:|
| `C3AS0.84H4.32` | (AlAl)O₃ | 1.0 |
| `C3AFS0.84H4.32` | (AlFe³⁺)O₃ | 0.5 |
| `C3FS0.84H4.32` | (Fe³⁺Fe³⁺)O₃ | 0.0 |

and CEMDATA18 declares the binary between the **middle and the iron end**
[Lothenbach2019](@cite) — spanning ``x(\mathrm{Al}) \in [0, 0.5]`` and no
further. A CEM I is iron-rich through its ferrite phase and never needs more. A
binder whose pozzolana brings twice as much aluminum as iron does, and the
declared phase **cannot go there**: the minimization will drive the solution to
its aluminum-richest admissible composition and stop, with aluminum left over and
nothing in the output saying the range was the binding constraint.

The remedy used on [the CEM IV](@ref ex-cem4-pozzolanic) and
[CEM V](@ref ex-cem5-composite) pages is to declare `C3AS0.84H4.32` as a separate
**pure phase** beside the binary, which is how the aluminum-rich half of the
series becomes reachable at all. It is an approximation and it is named as one:
as a pure phase it carries no mixing entropy, where a site-fraction model over
``x(\mathrm{Al}) \in [0,1]`` would.

!!! warning "Extending the solid solution to three end-members would be worse"
    The tempting fix — declare all three as an ideal ternary — is wrong, and
    wrong in a way that certifies. Three compositions of a two-site substitution
    are not three independent end-members; an ideal ternary over them counts
    configurations that do not exist and gets the mixing entropy wrong, so it
    returns a confident answer to a model nobody published. A solid solution is
    only ever as good as the model that was fitted for it: declare the range that
    was fitted, and handle what lies outside it explicitly.

## A habit worth adopting

Before trusting a new formulation, print the budget and ask where each component
goes:

```@example choosing
pure = split("C3S C2S C3A C4AF Gp Portlandite ettringite monosulphate12 " *
             "C3AH6 C3AS0.84H4.32 straetlingite hydrotalcite AlOHmic Brc")
for el in ("Al", "Mg", "S", "K")
    c = carriers(pure, el)
    @printf("  %-3s can go into %d declared phase(s): %s\n",
            el, length(c), isempty(c) ? "NONE — check the budget!" : join(c, ", "))
end
```

Potassium comes back empty in that list, which is the point: on a paste whose
clinker carries alkalis, that line is the warning the solver will not give you.
