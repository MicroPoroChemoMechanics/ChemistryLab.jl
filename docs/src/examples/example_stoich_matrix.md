# [Stoichiometric matrices](@id ex-stoich-matrix)

Every calculation in this package rests on one object: a matrix that says how
each species is built out of a chosen set of components. Equilibrium is the
minimization of Gibbs energy *subject to* `A n = b`, and `A` is that matrix. This
page shows it — over atoms, over primary species, and in cement notation — and
ends with a matrix whose entries are **symbols rather than numbers**, which is
what lets a phase of variable composition be carried through the algebra.

## What the matrix is

Take four clinker phases and ask what each is made of. Answer it in **atoms**
and you get the canonical stoichiometric matrix: one row per element, one column
per species, the entry being how many of that element the species contains.

```@example stoich
using ChemistryLab
using DynamicQuantities
using Symbolics
using LinearAlgebra

clinker = CemSpecies.(split("C3S C2S C3A C4AF"))
csm_atoms = CanonicalStoichMatrix(clinker)
pprint(csm_atoms.A, csm_atoms.primaries, clinker; label = :name)
```

Read a column: `C3S` is three calcium, one silicon, five oxygen — which is
``\ce{Ca3SiO5}``. Read a row: calcium appears in every one of the four phases,
with the multiplicities on that line.

## The same phases over oxide components

Atoms are not the only possible components. A cement chemist decomposes the same
four phases over **oxides**, and the matrix changes shape accordingly — four rows
instead of five, because there are four oxides and five elements:

```@example stoich
oxides = CemSpecies.(split("C S A F"))
cs_ox = ChemicalSystem(clinker, oxides)
pprint(cs_ox.SM.A, oxides, clinker; label = :name)
```

That is the same chemistry in a different basis, and the column for `C3S` now
reads what its name says: three lime, one silica. Neither basis is more correct;
the choice is made by what the problem is posed in. A datasheet reports oxides,
so the oxide basis is the one the [Bogue calculation](@ref) inverts.

!!! note "`C` is lime here, not carbon"
    The components above are [`CemSpecies`](@ref), so `C` is ``\ce{CaO}`` and `S`
    is ``\ce{SiO2}``. Written as ordinary [`Species`](@ref) the same letters mean
    carbon and sulfur, and the matrix would be a different one — see
    [the trap](@ref man-cement-notation-trap).

## Over primary species, from a database

For an aqueous system the natural components are not atoms but a set of
**primary species** — the ions a reaction is conventionally written over. A
`ChemicalSystem` computes that matrix at construction:

```@example stoich
substances = build_species(datapath("cemdata18-thermofun.json"); verbose = false)
species = speciation(substances, split("C3S Portlandite H2O@");
                     aggregate_state = [AS_AQUEOUS])
byname = Dict(symbol(s) => s for s in species)
primaries = [byname[s] for s in CEMDATA_PRIMARIES if haskey(byname, s)]

cs = ChemicalSystem(species, primaries)
length(cs.species), length(cs.SM.primaries)
```

The matrix itself, over the first few species so the table stays readable:

```@example stoich
pprint(cs.SM.A[:, 1:8], cs.SM.primaries, cs.species[1:8]; label = :symbol)
```

Each column is a **formation reaction**: the species in the header, written over
the primaries. `reactions` turns those columns back into reactions:

```@raw html
<details><summary>Every independent reaction of this system, reconstructed from the matrix</summary>
```

```@example stoich
reactions(cs.SM)
```

```@raw html
</details>
```

## A matrix with symbols in it

C-S-H has no fixed formula. Its lime-to-silica ratio, its aluminum uptake and its
water content all vary with the cement and with age, so a matrix that describes a
hydrated paste has to carry those as **parameters** rather than numbers.

This is the decomposition of [Chen2007a](@cite) and [Chen2007b](@cite), who write
a blended-cement paste over its oxides with a parameterized C-S-H. Take
``\ce{C_â S A_b̂ H_ĝ}`` alongside five hydrates of fixed composition:

```@example stoich
@variables â b̂ ĝ

CSH = CemSpecies(Dict(:C => â, :S => 1, :A => b̂, :H => ĝ))
hydrates = [
    CSH,
    CemSpecies("M5AH13"),      # hydrotalcite
    CemSpecies("C6AFS2H8"),    # siliceous hydrogarnet
    CemSpecies("C6AS̄3H32"),    # ettringite
    CemSpecies("C2ASH8"),      # strätlingite
    CemSpecies("C4AH13"),      # hydroxy-AFm
]

csm_h = CanonicalStoichMatrix(hydrates)
pprint(csm_h.A, csm_h.primaries, hydrates; label = :name)
```

The `â`, `b̂` and `ĝ` sit in the matrix exactly where numbers sit for the other
phases, and every operation below carries them through.

### Inverting it

Seven oxides and six hydrates is not square, so one row must go. Dropping `H`
expresses the hydrates over the **anhydrous** oxides, water following as a
consequence rather than as an input — which is how a hydration calculation is
posed:

```@example stoich
keep = findall(!=(:H), csm_h.primaries)
A_sq = Symbolics.Num.(csm_h.A[keep, :])
anhydrous = CemSpecies.(string.(csm_h.primaries[keep]))

A_inv = simplify.(inv(A_sq))
pprint(A_inv, hydrates, anhydrous; label = :name)
```

Each row is one hydrate written over the anhydrous oxides, **as a function of the
C-S-H composition**. Substituting a composition collapses it to numbers — here a
C-S-H at C/S = 1.8 taking a little aluminum:

```@example stoich
at = Dict(â => 1.8, b̂ => 0.1, ĝ => 4.0)
A_num = Float64.(Symbolics.value.(substitute.(A_inv, (at,))))
pprint(round.(A_num; digits = 4), hydrates, anhydrous; label = :name)
```

### Why carry the symbols at all

Because the derivative is then free. How much the ettringite content responds to
the C-S-H lime ratio is a partial derivative of an entry of this matrix, and with
the symbol still in place it is obtained by differentiating rather than by
re-solving at two nearby compositions:

```@example stoich
∂ = Symbolics.derivative(A_inv[1, 1], â)
simplify(∂)
```

## What this matrix does *not* decide

It is worth being explicit, because the matrix is necessary and nowhere near
sufficient. `A` fixes what is **conservable**: any composition `n` with
`A n = b` respects the element balance. It says nothing about which of those
compositions is the stable one — that is the Gibbs energy's business, and it
needs the thermodynamic data the matrix knows nothing about. A stoichiometric
matrix built from formulas alone, as on this page, can balance a reaction and
cannot predict whether it proceeds.
