# The zeolite extension of CEMDATA18

`data/cemdata18-zeolites.json` is CEMDATA18 with 28 zeolites appended. It is
**generated**, and `regenerate.jl` in this directory is what generates it:

    julia --project=docs data/zeolites/regenerate.jl

Run it to reproduce the shipped file. It is deterministic, and it refuses to
write anything if any of its three checks fails.

## Why a separate database, and why these phases

A pozzolanic or an alkali-activated binder at high alkalinity precipitates
zeolites. Without them in the species list the alkalis have nowhere to go but
the pore solution, and the calculated pH comes out too high — the error is in
the phase list, not in the solver. CEMDATA18 carries five zeolites —
`chabazite`, `natrolite`, `zeoliteP_Ca`, `zeoliteX` and `zeoliteY` — and the families
a blended cement actually needs — clinoptilolite, heulandite, mordenite,
phillipsite, analcime, stilbite, the gismondine/faujasite/LTA series, in both
their Na and their K forms — are not there.

The data are transcribed from two open-access papers by the laboratory that
produced CEMDATA18 itself:

| series | phases | source |
|:--|--:|:--|
| Na | 14 | Ma, B. & Lothenbach, B. (2020), *Synthesis, characterization, and thermodynamic study of selected Na-based zeolites*, **Cement and Concrete Research 135**, 106111, [10.1016/j.cemconres.2020.106111](https://doi.org/10.1016/j.cemconres.2020.106111) |
| K | 14 | Ma, B. & Lothenbach, B. (2021), *Synthesis, characterization, and thermodynamic study of selected K-based zeolites*, **Cement and Concrete Research 148**, 106537, [10.1016/j.cemconres.2021.106537](https://doi.org/10.1016/j.cemconres.2021.106537) |

Both DOIs were resolved against the Crossref REST API; the titles, authors,
journal, volumes, pages and years above are as Crossref returns them.

**Every number in `zeolite_data.jl` is copied from the published table.** None is
estimated, interpolated, averaged or adjusted. The `origin` field records the
paper's own footnote for `S⁰` and `Cp⁰` — `:measured`, `:additivity` or
`:literature` — and is metadata that no calculation reads.

### Sources that were examined and rejected

Merging two thermodynamic datasets is only defensible when they share a
reference state, and the usual failure is silent: an offset of a few kJ/mol on
`Na+` moves every dissolution equilibrium by an order of magnitude with no
solver complaining. Three other candidates were compared against CEMDATA18 on
every phase they share:

- **MINES19** offsets `Na+` by +5.71 kJ/mol and `K+` by +5.57 — roughly a factor
  of ten on activity, through the ions every dissolution reaction passes.
  Rejected.
- **SLOP16** carries literal placeholders (`G = 1000 J/mol`) for chabazite and
  natrolite. Rejected.
- **PSI/Nagra** agrees exactly (median |Δ| = 0.000 kJ/mol over 108 shared phases)
  but adds nothing here: of its 578 entries absent from the shipped file, 570 are
  radionuclides.

## What the generator refuses to do

1. **Overwrite.** A zeolite whose symbol already exists in CEMDATA18 aborts the
   build. Nothing in the base dataset is modified, so `natrolite` and `NAT-Na`
   coexist and the caller chooses between them.
2. **Emit an unbalanced dissolution.** Each phase's congruent dissolution is
   re-derived from its formula string — parsed by the same parser
   `build_species` will use — and must balance element by element and in charge.
   A formula this script accepted but the package could not read would otherwise
   fail much later and much less legibly.
3. **Emit a phase whose `log Ksp` does not close.** Both papers publish `log Ksp`
   *and* `ΔfG⁰` for the same phase, referred to the CEMDATA18 primary species.
   Recomputing one from the other through CEMDATA18's own aqueous Gibbs energies
   closes the loop: it can only agree if the reference states coincide **and**
   the transcription is exact. The tolerance is 0.05 log units, and all 28 agree
   to within 0.026.

That third check is the reason this merge can be offered at all, and it runs at
every regeneration rather than being asserted in a comment.

## A second, independent consistency check

On the phases both CEMDATA18 and the papers carry, the papers' `S⁰`, `Cp⁰` and
`V⁰` reproduce the shipped values to rounding — natrolite: 360 against
359.73 J/(mol·K), 359 against 359.23, 169.36 against 169.2 cm³/mol. Only `ΔfG⁰`
differs, and deliberately: it is a revision from new solubility measurements
(−5305.15 against −5325.70 kJ/mol). Nothing is overwritten, so both values remain
available under different symbols.

## Using it

```julia
using ChemistryLab
substances = build_species(datapath("cemdata18-zeolites.json"); verbose = false)
```

The file is a superset of `cemdata18-thermofun.json`: every symbol the base file
carries is present with the same values, so a script that loads the extension
instead of the base gets the same answer unless it names one of the new phases.
