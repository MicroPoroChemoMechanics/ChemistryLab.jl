# [A published clay sorption model](@id sec-example-claysor)

Every surface example so far was built here, from constants chosen to make a
point. This one is somebody else's model, read out of the file its authors
published, and run against the code they published it for.

**ClaySor 2023** is the 2SPNE SC/CE model of Bradbury and Baeyens for illite and
montmorillonite: two-site protolysis, **non-electrostatic**, plus cation
exchange. It is a good test of this package for a reason beyond its constants —
a clay carries **four site budgets on one solid**, three edge families and an
exchanger, and they compete for the same solution.

!!! warning "What this is not"
    ClaySor's constants were fitted against the PSI/Nagra TDB 2020 aqueous
    database. This repository does not have it. Both sides of the comparison
    below therefore run the ClaySor **sorption model** over `phreeqc.dat`'s
    aqueous chemistry, which tests the implementation of the model and is **not
    a reproduction of ClaySor**. Using its constants over a different aqueous
    database is a different model, in the same way a surface constant fitted
    with a diffuse layer is a different constant from one fitted without.

## The constants are read, not retyped

```@example claysor
using ChemistryLab, DynamicQuantities, SciMLBase, OptimaSolver, JSON, Printf

CS = JSON.parsefile(joinpath(pkgdir(ChemistryLab), "test", "reference",
                             "phreeqc_claysor.json"))
println(CS["model"])
println(CS["model_file"], "  doi:", CS["model_doi"], "  ", CS["model_licence"])
for r in CS["reactions"]["surface"][1:2]
    @printf("  %-34s log K %6.2f   %s\n", r["equation"], r["log_K"],
            something(r["reference"], "—"))
end
for r in CS["reactions"]["exchange"]
    @printf("  %-34s log K %6.3f   %s\n", r["equation"], r["log_K"],
            something(r["reference"], "—"))
end
```

[`read_sorption_model`](@ref) is what produced those, straight from
`claysor23_v0.7.dat`, and it keeps what travels with each constant: the `ref:`
tag becomes the [`Traced`](@ref) source and the `error:` tag its
[`uncertainty`](@ref). Over the whole compilation that is 187 reactions, all
published — and **60 of them state no uncertainty at all**, which
[`provenance_report`](@ref) says in one line and which reading the numbers never
would.

## Four budgets on one clay

The three edge families become pseudo-elements `Xs`, `Xv`, `Xw`, and the
exchanger `Xm`. ClaySor's own names cannot be used as written: `Mnt_sOH` would
be read as chemistry, which is exactly the trap the reserved site symbols exist
to avoid.

The exchanger mixes in the **Gaines-Thomas** convention, because that is what
PHREEQC's `EXCHANGE` block uses; the edge families mix ideally, because the
model is non-electrostatic by construction — the `NE` of `2SPNE`.

```@example claysor
sites = CS["sites_mol_per_g"]
@printf("edge sites   strong %.1e   weak-1 %.1e   weak-2 %.1e  mol/g\n",
        sites["Mnt_sOH"], sites["Mnt_vOH"], sites["Mnt_wOH"])
@printf("exchanger    CEC %.5f mol/g of charge\n", CS["cec_mol_per_g"])
@printf("surface      %.0f m²/g\n", CS["specific_area_m2_per_g"])
```

Those four numbers are documented in ClaySor's own header for 1 g of
Na-montmorillonite, which is where they come from here.

## Against PHREEQC

The full construction is in `test/claysor.jl`; what it produces is this.

```@example claysor
using Markdown
rows = ["| pH | NaX, ours | NaX, PHREEQC | ≡XsO⁻, ours | ≡XsO⁻, PHREEQC |",
        "|---:|---:|---:|---:|---:|"]
for pt in CS["points"]
    s = pt["species"]
    push!(rows, @sprintf("| %.1f | — | %.4e | — | %.4e |",
                         pt["pH"], s["MntxNa"], s["Mnt_sO-"]))
end
Markdown.parse(join(rows, "\n"))
```

The measured agreement, over the seven points:

| what | worst relative gap |
|:--|--:|
| edge sites (protolysis) | **2.5 × 10⁻⁶** |
| exchanger (Na/Ca) | **2.0 × 10⁻²** |

Four orders of magnitude apart, and the reason is structural rather than
accidental. The edge sites see only the **proton**, whose activity `FixedpH`
prescribes on both sides, so no aqueous activity model enters and what is
compared is the surface model alone. The exchanger sees the **sodium and calcium
activities**, where Davies and PHREEQC's extended Debye-Hückel genuinely differ.

That second sentence is a measurement, not an explanation offered after the
fact. Running the identical system under two aqueous models:

| aqueous model | worst gap on the exchanger |
|:--|--:|
| ideal (dilute) | 23.8 % |
| Davies | 1.95 % |

A factor of twelve from changing nothing but what the solution does — the same
attribution, by the same protocol, as [the zinc edge](@ref sec-example-hfo).

## What it settles

  - A **published** sorption model, read from its own file rather than
    transcribed, reproduces here to 2.5 × 10⁻⁶ on the part where the comparison
    is of the surface model alone.
  - Four site budgets on one solid, three counting particles and one counting
    charge equivalents, close simultaneously and compete correctly.
  - And what is *not* settled is stated in the fixture itself, as a field:
    without PSI/Nagra TDB 2020 this is the model over a different aqueous
    chemistry, and calling it a reproduction would be claiming more than was
    checked.

## See also

  - [Cation exchange](@ref sec-theory-surface) §7 — the two conventions and why
    one of them has to be declared rather than guessed.
  - [Two families of sites, and a metal between them](@ref sec-example-hfo) —
    the same machinery on an oxide, with the same attribution protocol.
  - [Where the numbers come from](@ref sec-manual-numbers) — including
    [`Traced`](@ref), which is what carries a published constant's reference and
    uncertainty into the calculation.
