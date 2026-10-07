# [Database Interoperability](@id sec-databases)

!!! info "Before this page"
    [Species](@ref sec-species).

So far, we have looked at the possibility of creating and manipulating any species, whether they exist or not. If we wanted to create a H₂O⁺⁴ molecule, it would not be a problem. However, you will admit that it is a little strange...

This is why ChemistryLab relies on existing databases, in particular [Cemdata18](https://www.empa.ch/web/s308/thermodynamic-data) and [PSI-Nagra-12-07](https://www.psi.ch/en/les/thermodynamic-databases). Cemdata18 is a chemical thermodynamic database for hydrated Portland cements and alkali-activated materials. PSI-Nagra is a Chemical Thermodynamic Database. The formalism adopted for these databases is that of [Thermofun](https://thermohub.org/thermofun/thermofun/) which is a universal open-source client that delivers thermodynamic properties of substances and reactions at the temperature and pressure of interest. The information is stored in json files.

## [Where the databases come from](@id sec-datapath)

The thermodynamic databases ChemistryLab reads are published by their authors,
under their own terms: Cemdata18, PSI/Nagra 12/07, aq17 and SUPCRT's slop98
[Johnson1992](@cite), in the ThermoFun format that the ThermoHub project
distributes, and Cemdata18 again in PHREEQC format [ParkhurstAppelo2013](@cite),
distributed by Empa. ChemistryLab obtains each one from its publisher the first
time it is needed, checks it, and keeps it. You name a database by its file name
with [`datapath`](@ref), and the rest is automatic:

```@example datapath
using ChemistryLab

path = datapath("cemdata18-thermofun.json")
isfile(path)
```

The value is an absolute path, so it is not shown here: it depends on the
machine. [`database_info`](@ref) says where each database currently resolves,
without downloading anything:

```@example datapath
database_info();
nothing # hide
```

The statuses read as follows:

| status | meaning |
|:--|:--|
| `local` | found in the directory named by `ENV["CHEMISTRYLAB_DATABASE_DIR"]` |
| `cached` | in ChemistryLab's cache, ready |
| `downloadable` | not yet obtained; the first call that needs it will download it |
| `manual` | to be downloaded by hand once (see [the PHREEQC file of Empa](@ref sec-empa-dat)) |
| `buildable` | a database ChemistryLab builds; the first call that needs it will build it |

### What happens on the first call

A call such as `build_species(datapath("cemdata18-thermofun.json"))` looks for
the file in three places, in this order:

1. the directory named by `ENV["CHEMISTRYLAB_DATABASE_DIR"]`, if you set one,
   and its immediate subdirectories;
2. ChemistryLab's cache, a directory of the Julia depot that
   `ChemistryLab.database_cache()` returns;
3. the publisher. The ThermoFun files come from release v1.1.1 of ThermoHub
   ([github.com/thermohub/thermohub](https://github.com/thermohub/thermohub),
   archived as [10.5281/zenodo.7385311](https://doi.org/10.5281/zenodo.7385311)),
   from GitHub or, failing that, from the jsDelivr mirror of the same commit.

Every file is identified by the SHA-256 of the version ChemistryLab was
validated with. A download whose content differs is discarded, and a cached
file that no longer matches is removed and obtained again, so a damaged or
different file is never used without your knowing it. The first download of
each file prints one message naming its source, the license its publisher
states and the reference to cite. Later calls read the cache and print nothing.

Four databases are **built** by ChemistryLab on first use rather than downloaded:
`cemdata18-zeolites.json` ([below](@ref sec-zeolites)),
`cemdata18-chloride.json` ([below](@ref sec-chloride-extension)),
`cemdata18-asr.json` ([below](@ref sec-asr-extension)) and
`cemdata18-cashplus.json` ([the models of the C-S-H gel](@ref sec-csh-models)).
Each is the downloaded Cemdata18 file, copied through unchanged, with
ChemistryLab's own additions appended; the CASH+ build also replaces the aqueous
complex `CaSiO3@` by the value refitted together with that model
[Kulik2022](@cite), since neither holds with the other's Cemdata18 value. The build takes a few
seconds once. It is cached, and it is redone automatically whenever Cemdata18 or
the added data change.

The cache can be deleted at any time: what it holds is obtained or built again
when next needed.

### [Working without network access](@id sec-database-offline)

Three ways, from the simplest:

- **In advance.** [`fetch_databases`](@ref) obtains every database that can be
  downloaded and builds the derived ones. Run it once while connected;
  everything afterwards reads the cache.
- **By hand.** Download a file yourself from an address that
  [`database_info`](@ref) prints, on any machine, then call
  [`install_database`](@ref) with its path. The checksum is verified before the
  file enters the cache.
- **Your own directory.** Point `ENV["CHEMISTRYLAB_DATABASE_DIR"]` at a
  directory holding the files, for instance one shared by a team, before the
  first call. Nothing is then downloaded for the files it holds.

### [The PHREEQC file of Empa](@id sec-empa-dat)

`CEMDATA18-31-03-2022-phaseVol.dat` is Cemdata18 in PHREEQC format. Empa
distributes it from its
[thermodynamic data page](https://www.empa.ch/web/s308/thermodynamic-data),
which a web browser can use and a program cannot. It is therefore the one
database you download by hand, once:

```julia
# after downloading the file with a web browser
install_database(expanduser("~/Downloads/CEMDATA18-31-03-2022-phaseVol.dat"))
```

Nothing ChemistryLab does by default needs it. It is the input of
[`merge_json`](@ref) and of the PHREEQC readers when you want to work with that
file itself. If Empa has published a newer version since, `install_database`
refuses it with the message explained below. Pass `force = true` to use it
anyway, knowing that results may then differ from the documented ones.

### [Reading the messages](@id sec-database-messages)

Every message names the file concerned. What each one means, and what to do:

- **`Info: ChemistryLab downloaded 'cemdata18-thermofun.json' (…) and verified
  its checksum.`**, followed by `source`, `license`, `cite` and `stored`. This
  is the first use of that file on this machine. Nothing to do; the lines say
  where it came from, under which terms, what to cite, and where it is kept.
- **`DatabaseUnavailable: ChemistryLab could not obtain '…'`** followed by
  **`Downloading it failed:`**. There is one line per address tried, with the
  reason:
  - a network error, such as `Could not resolve host` or a timeout, means the
    machine is offline or behind a proxy. Use one of the three ways of
    [working without network access](@ref sec-database-offline), which the
    message lists too;
  - **`the content differs from the validated version (SHA-256 …)`** means that
    something between you and the publisher altered the file, or that the
    address serves another version. Download the file from the publisher's page
    and install it with [`install_database`](@ref): if it is not the validated
    version either, that function will say so.
- **`DatabaseUnavailable: … a web page that a program cannot use`**. This is the
  Empa PHREEQC file, which has to be
  [downloaded by hand](@ref sec-empa-dat) once.
- **`Warning: '…' in CHEMISTRYLAB_DATABASE_DIR is not the version ChemistryLab
  was validated with`**, with the `found` and `expected` checksums. Your copy is
  used as you asked, but it is a different release, so results may differ from
  the documented ones. If a subdirectory also holds the validated version,
  that one is chosen and no warning appears.
- **`ArgumentError: '…' is not the version of '…' ChemistryLab was validated
  with`**, from `install_database`. The file you downloaded is another release.
  Check that it is the file you meant; to use it anyway, call
  `install_database(path; force = true)`.
- **`Warning: '…' installed in a version ChemistryLab was not validated with`**.
  This confirms that `force = true` installed another release.
- **`ArgumentError: '…' is not a database ChemistryLab obtains`**, or **`no
  such data file`**. The name is misspelled, or it is a file of your own that is
  not where you said. Both messages list the names that exist.

!!! note "A file of your own takes precedence"
    The readers try a path against the working directory **first**, so a
    database file of your own, given by its path, is always the one read. A
    database name alone, as in `datapath("cemdata18-thermofun.json")`, goes
    through the three places above.

## ThermoFun metadata

Classification fields contain exact enum names, such as `AS_AQUEOUS` and
`SC_AQSOLUTE`. Missing, malformed, or unknown labels retain the undefined-state
fallback (`AS_UNDEF` or `SC_UNDEF`).

Unit fields accept numbers, registered unit symbols, arithmetic (`+`, `-`, `*`,
`/`, `//`, `^`), square and cube roots, and `Constants` names. Examples include
`J/(mol*K)`, `1e-05/K`, and `K^(1//2)`. Other function calls and executable Julia
syntax are rejected before unit parsing. An invalid or unsupported unit retains
the caller's default unit; rejection does not establish that the supplied data
are scientifically valid. Database metadata cannot define custom Julia code.

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

## Loading species from a database

The simplest way to load species from a ThermoFun-compatible JSON file is `build_species`, which reads the file and directly returns a `Vector{Species}` with compiled thermodynamic functions:

```julia
using ChemistryLab
all_species = build_species(datapath("cemdata18-thermofun.json"))
```

Each species already carries its molar mass and temperature-dependent thermodynamic functions (Cp⁰, ΔₐH⁰, ΔₐS⁰, ΔₐG⁰, logK⁰) as `SymbolicFunc`s and `NumericFunc`s.

!!! note "Low-level access"
    If you need the raw DataFrames (e.g. to inspect metadata or filter on database columns), the lower-level function `read_thermofun_database` is still available and returns three DataFrames `(df_elements, df_substances, df_reactions)`:
    ```julia
    df_elements, df_substances, df_reactions = read_thermofun_database(datapath("cemdata18-thermofun.json"))
    ```
    `build_species(df_substances)` can then be called on the filtered DataFrame.

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

## Primary species extraction

It is also possible to retrieve primary species from the Cemdata18 database in
PHREEQC format. Primary species are a minimal subset such that every other
species can be expressed as their linear combination. The file is the one
[installed by hand from Empa](@ref sec-empa-dat), which is why this block is not
run when the documentation is built:

```julia
df_primaries = extract_primary_species(datapath("CEMDATA18-31-03-2022-phaseVol.dat"))
show(df_primaries, allcols=true, allrows=true)
```

---

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
