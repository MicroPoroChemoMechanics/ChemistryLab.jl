# [Database extensions, filters and solid solutions](@id sec-databases)

!!! info "Before this page"
    [Importing thermodynamic databases](@ref sec-importing-databases).

Once a database is read ([Importing thermodynamic databases](@ref
sec-importing-databases)), three things remain before a calculation: the phases
a cement needs and the database lacks, which ChemistryLab adds in databases of
its own built on the downloaded one; the species of the problem, selected from
the hundreds the database holds; and the solid solutions, which group species
into phases with a mixing model.

## [The zeolite extension of CEMDATA18](@id sec-zeolites)

`cemdata18-zeolites.json` is a database ChemistryLab **builds** from Cemdata18
on first use: CEMDATA18 unchanged, with 28 zeolites appended.

The reason it exists is a gap in the phase list. A pozzolanic or an
alkali-activated binder at high alkalinity precipitates zeolites; without them in
the species list the alkalis have nowhere to go but the pore solution, and the
calculated pH comes out too high. A blended binder taken to full reaction can
also stop certifying without them, and that is a different matter: there it is
the Debye–Hückel limiting law (`å = 0`) that no longer closes, and a model with an ion size per ion certifies the same paste, as
[the CEM IV page](@ref ex-cem4-pozzolanic) measures. CEMDATA18 carries five
zeolites, and the families a blended cement actually needs — clinoptilolite,
heulandite, mordenite, phillipsite, analcime, stilbite and the
gismondine/faujasite/LTA series, in both their Na and their K forms — are not
among them.

The data are transcribed, number by number, from two open-access papers by the
laboratory that produced CEMDATA18 itself: the Na series from
[MaLothenbach2020](@cite) and the K series from [MaLothenbach2021](@citet).

```@example datapath
using ChemistryLab

base = build_species(datapath("cemdata18-thermofun.json"); verbose = false)
ext = build_species(datapath("cemdata18-zeolites.json"); verbose = false)
length(base), length(ext)
```

The extension is a strict superset: every symbol of the base file is present with
the same values, so a script that loads it instead of the base gets the same
answer unless it names one of the new phases.

```@example datapath
added = setdiff(Set(symbol.(ext)), Set(symbol.(base)))
println(length(added), " phases added")
for s in sort(collect(added))
    print(s, "  ")
end
```

!!! note "Nothing is overwritten, so both values stay available"
    [MaLothenbach2020, MaLothenbach2021](@citet) revise `ΔfG⁰` for phases
    CEMDATA18 already carries — natrolite by about 20 kJ/mol, from new
    solubility measurements. The new entries therefore take **distinct
    symbols**: `NAT-Na` sits beside `natrolite`, and the caller chooses.
    Declaring both in one system would count the same substance twice.

!!! warning "What makes the merge defensible, and what would make it wrong"
    Merging two thermodynamic datasets is only legitimate if they share a
    reference state, and the usual failure is silent: an offset of a few kJ/mol
    on `Na+` moves every dissolution equilibrium by an order of magnitude without
    any solver complaining.

    [MaLothenbach2020, MaLothenbach2021](@citet) both publish `log Ksp` **and**
    `ΔfG⁰` for each phase, referred to the CEMDATA18 primary species. The build
    recomputes one from the other through CEMDATA18's own aqueous Gibbs energies
    and **refuses to produce the database** if any phase misses by more than
    0.05 log units (`data/zeolites/regenerate.jl` prints the comparison). All 28
    agree to within 0.026. It refuses on two further grounds: a symbol that
    would overwrite a CEMDATA18 entry, and a dissolution reaction that does not
    balance in elements and charge when re-derived from the formula string.

    `data/zeolites/README.md` records the whole provenance, including the three
    other candidate datasets that were examined and rejected on measurement.

## [The chloride extension of CEMDATA18](@id sec-chloride-extension)

`cemdata18-chloride.json` is CEMDATA18 unchanged, with two substances appended.

The first is an end member of C-S-H, `CSHQ-Cl` = (CaCl₂)₀.₅. The phase `CSHQ_Cl`
of `data/solid_solutions.toml` is CSHQ with it, and is built from this file only.
Its Gibbs energy is fitted, not measured: on the chloride bound by C-S-H in the
sorption tests of [Hirao2005](@cite), by `data/chloride/regenerate.jl`. The
fitted value and its provenance are kept in `data/chloride/cshq_cl.json`, from
which the end member is built, and carried on the end member itself.

!!! warning "An end member fitted by ChemistryLab, not a published one"
    `CSHQ-Cl` is not part of CEMDATA18, and no paper gives it. One number is
    fitted, its Gibbs energy of formation from ½ Ca²⁺ + Cl⁻, on three points of
    one gel: C-S-H formed by alite with its portlandite, in NaCl solutions at
    20 °C, up to 1 mol/L. The fit is within 0.05 mmol/g at 0.5 and 1 mol/L and
    four times the measurement at 0.1 mol/L. Its dependence on the Ca/Si of the
    gel is not fitted, the tests holding one gel at portlandite saturation: it
    is what ideal mixing gives, the formula having been chosen so that the gel
    binds more chloride at a higher Ca/Si, as [Zibara2008](@citet) measured.
    Fitted on NaCl solutions alone, it does not reproduce the effect of the
    cation: in the CEM II concrete of [Tran2018](@citet), a CaCl₂ solution binds
    two to three times as much chloride as a NaCl one at the same free
    chloride, through the calcium the C-S-H adsorbs.
    It is an effective description: chloride does not adsorb specifically on
    C-S-H [Plusquellec2016](@cite), and the end member lumps the chloride that
    accompanies the calcium the surface adsorbs, saying nothing of where it sits.
    The entropy, volume and heat capacity of the end member are estimates, which
    enter only between 20 °C and the temperature of a calculation.

The second is Fe-Friedel's salt, `C4FCl2H10` = Ca₄Fe₂Cl₂(OH)₁₂·4H₂O, which the
Cemdata18 paper tabulates [Lothenbach2019; Tables 1 and 2](@cite) and its
ThermoFun export lacks. Its record is the row of Table 1, transcribed in
`data/literature/Lothenbach2019.json`, and the build refuses to write it unless
its `log Ks0` recomputed through the aqueous species of Cemdata18 is the one of
Table 2 within 0.05 (it is within 0.011). With Friedel's salt it forms the ideal
solid solution `Friedel_AlFe` of `data/solid_solutions.toml`, note k of Table 1.

```@example datapath
ext_cl = build_species(datapath("cemdata18-chloride.json"); verbose = false)
setdiff(symbol.(ext_cl), symbol.(base))
```

[Chloride binding in CEM III/A and CEM III/B](@ref sec-example-chloride-blended) uses
it, and says what it describes and what it does not. `data/chloride/README.md`
records the data, the fit, and the candidate that was rejected.

## [The products of the alkali-silica reaction](@id sec-asr-extension)

`cemdata18-asr.json` is CEMDATA18 unchanged, with K-shlykovite,
KCaSi₄O₈(OH)₃·2H₂O, and Na-shlykovite, NaCaSi₄O₈(OH)₃·2.3H₂O, appended: the
crystalline products of the alkali-silica reaction that [Jin2023](@citet) give
standard properties for. Their solubility products were measured at 80 °C
[ShiLothenbach2019](@cite) and carried to 25 °C by the authors, with an entropy
and a heat capacity estimated from the volume of the formula unit; the records
are the rows of their Table 1, transcribed in `data/literature/Jin2023.json`.
The build refuses to write them unless their `log Ks0` at 25 °C recomputed
through the aqueous species of Cemdata18 is the published one within 0.05 (they
agree to 0.02), and at 80 °C the database gives back the values the authors
refined, −28.0 and −30.8 against −28.3 ± 1.2 and −30.8 ± 0.4.

The enthalpy of each record is the one its Gibbs energy and entropy give with
the element entropies of Cemdata18. [Jin2023; Table 1](@citet) agrees with it
to 0.5 kJ/mol; its corrigendum [Jin2024](@cite) gives enthalpies 572 kJ/mol
higher with unchanged Gibbs energies and entropies, and the record uses the
enthalpy obtained from ``\Delta_f G^\circ``, ``S^\circ`` and the element
entropies. The equilibrium does not depend on it, the Gibbs energy at
temperature being built from ``\Delta_f G^\circ``, ``S^\circ`` and ``C_p^\circ``; a
heat would. Both published values are kept on each record.

```@example datapath
ext_asr = build_species(datapath("cemdata18-asr.json"); verbose = false)
setdiff(symbol.(ext_asr), symbol.(base))
```

## Filtering species with `speciation`

In practice, only a small subset of the database is relevant to a given problem. `speciation` filters a species list to those whose atomic composition is a subset of the atoms found in a set of *seed* species:

```@example database
using ChemistryLab #hide
all_species = build_species(datapath("cemdata18-thermofun.json")) #hide
# Keep only species that can form from the calcite / water system
species_calcite = speciation(all_species, split("Cal H2O@");
                             aggregate_state=[AS_AQUEOUS],
                             exclude_species=split("H2@ O2@ CH4@"))
dict_species_calcite = Dict(symbol(s) => s for s in species_calcite)
```

The `aggregate_state` keyword restricts results to aqueous species (the seed species themselves — `Cal` and `H2O@` — are always kept through `include_species` internally). For example, the properties of Ca(HCO₃)⁺ can then be read as:

```@example database
dict_species_calcite["Ca(HCO3)+"]
```

### `speciation` signatures

`speciation` accepts seed arguments in three forms:

| Seed argument | Description |
|:--------------|:------------|
| `Vector{Symbol}` | Explicit list of atom symbols |
| `Vector{<:AbstractSpecies}` | Species objects — their union of atoms defines the space |
| `Vector{<:AbstractString}` | Species symbol strings — looked up in `species_list` |

Common keyword arguments:

| Keyword | Default | Description |
|:--------|:--------|:------------|
| `aggregate_state` | all states | restrict to `[AS_AQUEOUS]`, `[AS_CRYSTAL]`, etc. |
| `class` | all classes | restrict to `[SC_AQSOLUTE]`, etc. |
| `exclude_species` | `[]` | species (or symbols) to always exclude |
| `include_species` | `[]` | species to always include regardless of composition |

## Building solid solution phases from a TOML file

Solid solution phases (e.g. C-S-H gel, AFm) group several end-member species with a
mixing model. Because database species carry `class = SC_COMPONENT` rather than
`class = SC_SSENDMEMBER`, they cannot be passed directly to [`SolidSolutionPhase`](@ref).

[`build_solid_solutions`](@ref) automates the full pipeline:

1. Reads phase definitions from a TOML file.
2. Looks up each end-member in a pre-built species dictionary.
3. Requalifies end-members to `SC_SSENDMEMBER` via [`with_class`](@ref).
4. Constructs and returns a `Vector{SolidSolutionPhase}`.

### TOML format

Each `[[solid_solution]]` entry specifies a phase name, the list of end-member
symbols (as they appear in the database), and the mixing model:

```toml
# Ideal solid solution (any number of end-members)
[[solid_solution]]
name        = "CSHQ"
end_members = ["CSHQ-TobD", "CSHQ-TobH", "CSHQ-JenH", "CSHQ-JenD",
               "KSiOH", "NaSiOH"]
model       = "ideal"
source      = "Lothenbach2015"

# Binary Redlich-Kister (exactly 2 end-members): published dimensionless
# Guggenheim parameters, read from data/literature (a = α R T at 298.15 K), or
# a0, a1, a2 in J/mol. It unmixes, so it takes two coexisting compositions.
[[solid_solution]]
name        = "AFm_SO4_OH"
end_members = ["C4AH13", "monosulphate12"]
model       = "redlich_kister"
guggenheim  = "Lothenbach2019:AFm SO4/OH"
instances   = 2
source      = "Lothenbach2019"
```

The file `data/solid_solutions.toml` shipped with ChemistryLab.jl declares the main
cemdata18 solid solutions (CSHQ, AFm_SO4_OH, AFt_SO4_CO3, Hydrogarnet, Hydrotalcite,
and others). `AFm_SO4_OH` and `AFt_SO4_CO3` carry the non-ideal parameters of
[Lothenbach2019](@cite), with the end-member order that reproduces the gaps the
article prints; the others are ideal, and the comments of the file say which of
those ideal models is an assumption. Monosulfate and monocarbonate are pure phases
in Cemdata18, and the file does not mix them.

### Usage

```julia
using ChemistryLab, DynamicQuantities

substances = build_species(datapath("cemdata18-thermofun.json"))
dict       = Dict(symbol(s) => s for s in substances)

# Build all solid solution phases defined in the TOML
ss_phases = build_solid_solutions(datapath("solid_solutions.toml"), dict)

# Take the ones the system needs, by name: the file holds alternatives
chosen = [p for p in ss_phases if name(p) in ("CSHQ", "C3(AF)S0.84H", "Ettringite_ss")]
cs = ChemicalSystem(species, CEMDATA_PRIMARIES; solid_solutions = chosen)
```

The file is a catalog, not a phase model to take whole. `CSHQ` and `CNASH_ss`
are two models of one C-S-H gel, which `ChemicalSystem` refuses together, and
`Ettringite_ss` and `AFt_SO4_CO3` describe the same aluminate sulfate under two
normalizations: `ettringite03_ss` is ettringite divided by three, and their
Gibbs energies agree to that factor within a few J/mol. With ideal mixing on
both, declaring the two is harmless, and the CEM I example does so on purpose.
With a non-ideal model on either, `ChemicalSystem` warns and names the pair:
moving the substance from one phase to the other then changes the Gibbs energy
by next to nothing, a flat direction on which the certified search can stop
short of the solution. Keep the phase whose mixing the problem needs, and drop
the other description.

Phases whose end-members are not found in `dict` are skipped with a warning
(pass `skip_missing = false` to raise an error instead).

!!! note "Manual requalification"
    If you only need one or two end-members, you can requalify them individually
    with [`with_class`](@ref) instead of going through the TOML file:

```julia
em = with_class(dict["CSHQ-TobD"], SC_SSENDMEMBER)
```

### [The models of the C-S-H gel](@id sec-csh-models)

The calcium silicate hydrate of a cement paste is one gel, and the file ships six
models of it. They differ in what their end-members can hold, and that decides
what a calculation can say: a model with no aluminum member puts every atom of
aluminum in another phase, and one with no potassium member leaves the potassium
in solution. `data/gel_models.toml` lists the six as models of one gel, so that a
system declaring two of them is refused. The table is read from the two files,
and the formulas from the database that carries every end-member, Cemdata18 with
those of CASH+:

```@example gels
using ChemistryLab, Printf, TOML
subs = Dict(symbol(s) => s for s in build_species(datapath("cemdata18-cashplus.json"); verbose = false))
entries = TOML.parsefile(datapath("solid_solutions.toml"))["solid_solution"]
@printf("%-9s %-8s %-16s %-15s %s\n", "model", "members", "mixing", "source", "elements besides Ca, Si, O, H")
for gel in TOML.parsefile(datapath("gel_models.toml"))["gel_model"]
    entry = only(e for e in entries if e["name"] == gel["model"])
    members = entry["end_members"]
    others = sort(unique(String(el) for m in members for el in keys(atoms(subs[m]))
                         if !(el in (:Ca, :Si, :O, :H))))
    @printf("%-9s %-8d %-16s %-15s %s\n", gel["model"], length(members), entry["model"],
            entry["source"], isempty(others) ? "none" : join(others, ", "))
end
```

`CNASH_ss` mixes on the sites of its formula unit, as its authors define it
[Myers2014](@cite). `CASH+` does too, and adds the energies of the reciprocal
reactions between its end-members and interactions on each site
[Kulik2022](@cite); its twelve-member form `CASH+NK`, with sodium and potassium
[Miron2022a](@cite), is in the same file and is the one a cement paste needs (see
[the CASH+ page](@ref ex-cashplus-csh)). The other four mix their end-members
ideally, which is how Cemdata18 ships them, and `sublattice_model("Kulik2011:csh3t",
members)` gives the site form of `CSH3T` (see [Solid solutions](@ref
sec-theory-solid-solutions)). [The CEM IV page](@ref ex-cem4-pozzolanic) computes
one paste with `CSHQ` and with `CNASH_ss`, and shows how far apart the two answers
are.

## Where to go next

How a species read from a database carries its temperature dependence is the
subject of [Thermodynamic Functions](@ref sec-thermodynamics), and
[ChemicalSystem and ChemicalState](@ref sec-system-state) turns a list of such
species into the system an equilibrium is computed on. Selecting that list for
a cement is a decision of its own, treated in
[Choosing the species list](@ref man-choosing-species).
