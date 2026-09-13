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

## Six traps, each one met in practice

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
| an assemblage still containing anhydrous clinker | the solve never reached hydration |
| a sweep whose pH jumps around non-monotonically | isolated failed solves inside an otherwise fine scan |
| a saturation index of +7 on a phase at 1e-305 mol | a declared phase whose element is absent from the budget |

None of these says "add a phase", and all of them did mean exactly that.

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
