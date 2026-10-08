# [Importing thermodynamic databases](@id sec-importing-databases)

!!! info "Before this page"
    [Where the numbers come from](@ref sec-manual-numbers).

A chemical equilibrium is computed from the standard properties of the species
that may form, and those are not the package's: they are measured, fitted and
published by others, in thermodynamic databases. ChemistryLab reads the
databases of the main formats in use, builds their species, and leaves the
databases themselves where their publishers distribute them. None is stored in
the package. This page says which formats are read, how a file is obtained, what
is read from it and what is not, and how a database is extended.

!!! tip "Questions this page answers"
      - Which formats can be read, and with which function?
      - Where does a database file come from, and under what terms?
      - Why can two databases not be mixed in one system?
      - Which activity model goes with a database?
      - How is a database extended with a species of one's own?

## Two families of databases

A **database of formation properties** gives each species its standard Gibbs
energy, enthalpy and entropy of formation from the elements at a reference
state, with a model of their dependence on temperature and pressure: a
heat-capacity polynomial for a mineral, the revised Helgeson-Kirkham-Flowers
equations for an aqueous species [TangerHelgeson1988](@cite). Its energies are
counted from the elements, and any two such databases can be combined, as far
as their data are consistent.

A **database of reactions** gives instead, for each species, a reaction that
forms it from a few master species, one for each element, and the equilibrium
constant of that reaction as a function of temperature. It determines the
energies of its species only relative to its master species, which the package
sets to zero at every temperature: a choice of zero, a *gauge*, which leaves an
equilibrium unchanged within the database ([Energies counted from the
primaries](@ref sec-theory-gauge)).

| format | family | function |
|:--|:--|:--|
| PHREEQC (`.dat`) | reactions | [`read_phreeqc_database`](@ref) |
| The Geochemist's Workbench (`.tdat`) | reactions | [`read_gwb_database`](@ref) |
| EQ3/6 (`data0`) | reactions | [`read_eq36_database`](@ref) |
| ThermoFun (`.json`) | formation properties | [`read_thermofun_database`](@ref) |
| Reaktoro (`.yaml`) | formation properties | [`read_reaktoro_database`](@ref) |

Every reader returns the same three tables, `DataFrame`s: the elements, the
substances (one row per species, with its formula, aggregate state, charge and
the data of its standard properties) and the reactions the file writes out.
What concerns the database as a whole (its format, its source, what was read and
not used) is the metadata of the table of substances. [`import_database`](@ref)
reads any of the five formats, choosing the reader by the name of the file and
its first lines, and [`build_species`](@ref) builds the species of the table of
substances, whichever it is.

## [Where the databases come from](@id sec-datapath)

No database is stored in the package. Each is published by its authors, under
their own terms, and ChemistryLab obtains a copy from its publisher the first time
it is needed, when the publisher allows a program to do so, checks it, and keeps
it. You name a database by its file name with [`datapath`](@ref):

| database | publisher | terms | obtained |
|:--|:--|:--|:--|
| `phreeqc.dat`, `llnl.dat`, `minteq.v4.dat`, `wateq4f.dat`, `sit.dat`, `pitzer.dat` | U.S. Geological Survey, PHREEQC 3.7.3 [ParkhurstAppelo2013](@cite) | user rights notice of PHREEQC: copy, modification and distribution allowed with the notice | automatically |
| ThermoFun files of Cemdata18, PSI/Nagra 12/07, aq17 and SUPCRT's slop98 [Johnson1992](@cite) | ThermoHub, release v1.1.1 | GPL-3.0 | automatically |
| any other file: Cemdata18 in PHREEQC format, Thermoddem, ThermoChimie, PSI/Nagra 2020, the thermo datasets of The Geochemist's Workbench, data0 files of EQ3/6, the YAML files of Reaktoro | its publisher | the publisher's | by hand, once |

A file obtained by hand is downloaded with a web browser from its publisher's
page, whose terms the user accepts there, and given to [`install_database`](@ref)
or read from its own path: the package never downloads it, and never
redistributes it.

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
3. the publisher. The PHREEQC files come from the repository of PHREEQC at
   version 3.7.3 ([github.com/phreeqc-dev/phreeqc3](https://github.com/phreeqc-dev/phreeqc3)),
   the ThermoFun files from release v1.1.1 of ThermoHub
   ([github.com/thermohub/thermohub](https://github.com/thermohub/thermohub),
   archived as [10.5281/zenodo.7385311](https://doi.org/10.5281/zenodo.7385311)),
   from GitHub or, failing that, from the jsDelivr mirror of the same commit.

Every file is identified by the SHA-256 of the version ChemistryLab was
validated with. A download whose content differs is discarded, and a cached
file that no longer matches is removed and obtained again, so a damaged or
different file is never used without your knowing it. The first download of
each file prints one message naming its source, the license its publisher
states and the reference to cite. Later calls read the cache and print nothing.

Four databases are **built** by ChemistryLab on first use rather than downloaded
([Database extensions, filters and solid solutions](@ref sec-databases)):
`cemdata18-zeolites.json` ([the zeolites](@ref sec-zeolites)),
`cemdata18-chloride.json` ([the chloride phases](@ref sec-chloride-extension)),
`cemdata18-asr.json` ([the products of the alkali-silica reaction](@ref sec-asr-extension)) and
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

Nothing ChemistryLab does by default needs it; [`read_phreeqc_database`](@ref)
reads it like any PHREEQC database. If Empa has published a newer version since, `install_database`
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


## Reading a database of reactions

A PHREEQC database is read whole: its master species, its species and phases,
each resolved against the masters, with the log K of its formation as a function
of temperature, as the manual of PHREEQC defines it (`log_k`, `delta_h` in any of
its units, the six-term analytical expression, named expressions)
[ParkhurstAppelo2013](@cite). The composition of a species is that of the
balance of its reaction, not of its name.

```@example importing
using ChemistryLab, DynamicQuantities
_, db, _ = read_phreeqc_database(datapath("phreeqc.dat"))
first(db[:, [:symbol, :name, :aggregate_state, :charge, :gamma]], 6)
```

The column `formation` holds the log K of the formation of each species from the
master species, as a function of temperature. What is read and not used is
listed with its line in the metadata `notes`: transport data, species made of the
pseudo-elements a database defines to keep a gas out of redox equilibrium, input
blocks a database may carry.

```@example importing
using DataFrames: metadata
first(metadata(db, "notes"), 3)
```

[`build_species`](@ref) builds the species named, by their PHREEQC names or by
ChemistryLab's symbols (a neutral solute takes `@`, as `CO2@`):

```@example importing
species = build_species(db, ["H2O", "H+", "OH-", "Ca+2", "CO3-2", "HCO3-", "CO2", "CaCO3", "CaHCO3+", "CaOH+", "Calcite"])
cs = ChemicalSystem(species)
model = database_activity_model(db)
st = ChemicalState(cs; T = 298.15u"K", P = 1.0u"bar")
set_quantity!(st, "H2O@", 1.0u"kg")
set_quantity!(st, "Calcite", 1.0u"mol")
eq = equilibrate(st; model)
round(pH(eq, model); digits = 4)
```

## The activity model of a database

The constants of a database are fitted with an activity model, and a database is
used with that model. [`database_activity_model`](@ref) returns the one PHREEQC
applies, which the blocks of the database decide:

| the database carries | model | as PHREEQC computes it |
|:--|:--|:--|
| `PITZER` | [`PitzerActivityModel`](@ref) | a pair the block leaves out taken as zero |
| `SIT` | [`SITActivityModel`](@ref) | the activity of water from the osmotic coefficient |
| `LLNL_AQUEOUS_MODEL_PARAMETERS` | [`LLNLActivityModel`](@ref) | B-dot with A, B, Ḃ interpolated on the grid; carbon dioxide by its own formula |
| none of them | [`TruesdellJonesActivityModel`](@ref) | WATEQ for a species with `-gamma`, Davies otherwise |

In the last two, the activity of water is PHREEQC's,
``a_w = 1 - 0.017 \sum_i m_i`` [ParkhurstAppelo1999](@cite). These conventions
were measured against PHREEQC 3.7.3 before being written, and the package is
tested against it on the six databases it distributes, with a solution of seven
ions and with calcite and gypsum in water, at 25 and 60 °C
(`test/phreeqc_databases.jl`). Given the Debye–Hückel A and B PHREEQC computes
from its own model of water, the molalities agree within ``2\times10^{-5}``,
except at 60 °C, where PHREEQC's van 't Hoff equation takes a gas constant of
8.31470 J/(mol K) (``7\times10^{-5}``), and with `pitzer.dat`, whose `-ZETA` terms
the Pitzer model does not carry (``8\times10^{-5}``). With the A and B of the
package's own model of water, 0.27 % apart at 25 °C, they agree within a
percent. Under Pitzer, PHREEQC reports single-ion activities on the MacInnes
scale, which moves its pH, by 0.03 in that solution, and none of the molalities.

## Energies, and why two databases do not mix

A species read from a database of reactions has the standard Gibbs energy of its
formation from the master species, which are at zero:
``\mu_i^\circ(T) = -RT \ln 10\, \log_{10} K_i(T)``. A database of formation
properties counts its energies from the elements. Within one database the choice
does not change an equilibrium, beside the species of another it is
meaningless ([Energies counted from the primaries](@ref sec-theory-gauge)). Each
species records the zero of its database under `:gauge`, and
[`ChemicalSystem`](@ref) refuses a system that mixes two:

```@example importing
thermofun = build_species(datapath("cemdata18-thermofun.json"), ["Portlandite"])
try
    ChemicalSystem(vcat(species, thermofun))
catch err
    print(first(sprint(showerror, err), 160), "…")
end
```

## Reading a database of formation properties

A ThermoFun file gives each substance its standard properties at a reference
state and the methods that carry them to other temperatures and pressures:
heat-capacity polynomials over one or several intervals, the HKF equations of an
aqueous species, the equation of state of water for the solvent. A substance the
file defines by a reaction (440 in PSI/Nagra 12/07, 17 in Cemdata18) is computed
from that reaction at every temperature, from the log K of the reaction and the
properties of its other species, as ThermoFun computes it.

```@example importing
_, db18, _ = import_database(datapath("cemdata18-thermofun.json"))
size(db18)
```

The minerals of the aq17 database carry the volume and the order-disorder
transitions of [HollandPowell1998](@cite), its dissolved gases the equation of
[AkinfievDiamond2003](@cite); both are computed as these articles define them,
and tested against ThermoFun up to 150 °C and 100 bar
(`test/thermofun_methods.jl`). The minerals agree to 0.05 J/mol once one
convention is accounted for: ThermoFun integrates the volume from zero pressure,
which adds the molar volume times one bar away from the reference state, 10 J/mol
for albite. The dissolved gases agree to 8 J/mol, the order of the difference
between the two models of water they are computed with. A method the
package does not implement is not ignored: the substance is built for its
reference state, 25 °C and 1 bar, and anywhere else its properties raise an
error naming the method; [`build_species`](@ref) warns which substances are in
that case.

A YAML file of Reaktoro is read the same way, for the models the package
computes: `HKF`, `MaierKelley`, `HollandPowell` (with the modified Tait equation
of state of [HollandPowell2011](@cite)) and water
([`read_reaktoro_database`](@ref)), tested against Reaktoro
(`test/reaktoro_yaml.jl`).

### What a ThermoFun file may hold

Classification fields contain exact enum names, such as `AS_AQUEOUS` and
`SC_AQSOLUTE`. Missing, malformed, or unknown labels retain the undefined-state
fallback (`AS_UNDEF` or `SC_UNDEF`).

Unit fields accept numbers, registered unit symbols, arithmetic (`+`, `-`, `*`,
`/`, `//`, `^`), square and cube roots, and `Constants` names. Examples include
`J/(mol*K)`, `1e-05/K`, and `K^(1//2)`. Other function calls and executable Julia
syntax are rejected before unit parsing. An invalid or unsupported unit retains
the caller's default unit; rejection does not establish that the supplied data
are scientifically valid. Database metadata cannot define custom Julia code.

[`build_species`](@ref) also reads a file directly, without the tables:

```julia
all_species = build_species(datapath("cemdata18-thermofun.json"))
```

## Thermo datasets of The Geochemist's Workbench, data0 files of EQ3/6

Both are databases of reactions, read into the same representation as a PHREEQC
database. A thermo dataset gives its log K as a table at its principal
temperatures, which the package evaluates as the dataset's applications do
between those temperatures, by the polynomial of degree four fitted to the table,
or as the polynomial of format "jan19" [BethkeFarrell2026](@cite). A data0 file
gives its log K on a grid of eight temperatures, which the package evaluates as
EQPT does, by interpolating polynomials through the valid points of 0–100 °C and
of 100–300 °C [DavelerWolery1992](@cite). Their activity models are not built
from them: their conventions for neutral species and for water rest on formulas
their format documentation does not state, and a model is chosen explicitly.

## Primary species

The master species of a database of reactions, one per element with the proton,
the electron and water, are the primary species from which every other species
forms. [`extract_primary_species`](@ref) lists them, with the parameters of their
activity coefficients:

```@example importing
first(extract_primary_species(datapath("phreeqc.dat")), 6)
```

## Extending a database

A database is extended by a database built from it, never by editing its file:
the package's own extensions of Cemdata18 (the zeolites, the chloride phases, the
products of the alkali-silica reaction, the CASH+ model of the C-S-H gel) are
built on first use from the downloaded file and data published in articles,
transcribed into `data/literature/` with their source ([Database extensions,
filters and solid solutions](@ref sec-databases)). A species of one's own enters
a system as a [`Species`](@ref) with its standard properties, in the gauge of the
database it is combined with.

## Where to go next

The extensions ChemistryLab builds on Cemdata18, the filters that select the
species of a problem, and the solid solutions are the subject of [Database
extensions, filters and solid solutions](@ref sec-databases). How a species
carries its temperature dependence is that of [Thermodynamic Functions](@ref
sec-thermodynamics), and [ChemicalSystem and ChemicalState](@ref
sec-system-state) turns a list of species into the system an equilibrium is
computed on.
