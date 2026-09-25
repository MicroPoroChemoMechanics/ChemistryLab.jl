# The oracle bench

Three established codes, used as references rather than as authorities. Each
script here *generates* a JSON fixture that the tests read, with the version of
the code and the identity of the database it read written into it as fields.
Nothing in this directory runs during `Pkg.test()`.

## The rule that comes before any of them

**Derived closed forms first, codes second.** An agreement between two programs
sharing the same wrong formulation proves nothing, so a new model is checked
against an answer derived by hand — a Langmuir isotherm, a competitive
denominator, a reduced Gibbs energy and its Hessian — before any cross-code
comparison is attempted. The codes then answer a different question: does the
implementation reproduce what the field's tools produce on a case nobody can do
by hand?

And a comparison is only worth its name when the inputs match. The header of
`test/equilibrium_reference.jl` records what that costs when they do not:
leaving Reaktoro on its own HKF model instead of the ideal one shifted
∂Ca²⁺/∂(CO₂) by 35 %, which says nothing about either code.

## Two environments

They are separate because Reaktoro and GEMS3K have no reason to share a solve,
and because a broken dependency solve should not take the working oracle with
it.

```bash
# Reaktoro — already used by reaktoro_calcite_co2.py and reaktoro_coupling.py
conda create -n reaktoro-env -c conda-forge reaktoro thermofun

# PHREEQC and GEMS3K
conda create -n mpcm-oracles -c conda-forge python=3.12 xgems gems3k
"$(conda info --base)/envs/mpcm-oracles/bin/python" -m ensurepip --upgrade
"$(conda info --base)/envs/mpcm-oracles/bin/python" -m pip install phreeqpython
```

`conda create` does not install `pip` unless asked, which is why `ensurepip`
appears above rather than a bare `pip install`.

## What each code arbitrates

| code | env | what it settles | state |
|:--|:--|:--|:--|
| Reaktoro 2.13 | `reaktoro-env` | aqueous speciation; **ion exchange in both conventions** (`ActivityModelIonExchangeGainesThomas`, `…Vanselow`); reactive-area models, including a power law matching `ShrinkingCoreArea` | working |
| PHREEQC (IPhreeqc, via phreeqpython) | `mpcm-oracles` | **surface complexation**: `SURFACE` with and without a diffuse layer, `EXCHANGE`, and the Dzombak & Morel HFO parameters shipped in `phreeqc.dat` | working |
| GEMS3K | `mpcm-oracles` | Gibbs minimization with multi-site sorption phases, and CEMDATA18 — the cementitious case | **loads, cannot run**: see below |

## What is not measurable today, and why it is said rather than worked around

The plan wanted a tolerance floor from running the *same* ClaySor model in
PHREEQC and in GEMS, since the deposit ships both. Half of that works: the
deposit's `PHREEQC/claysor23_v0.7.dat` is a database IPhreeqc loads. The other
half does not, because `xgems.ChemicalEngine` needs a GEMS3K **system export**
(`*-dat.lst`) and the deposit ships a GEM-Selektor **database**. Producing the
export needs GEM-Selektor itself, a graphical application.

So that floor is open, and no substitute number is invented for it. See
`gems_bench.py`, which states the blocker and what unblocks it.

## The fixtures are JSON, and why

Each generator writes a `*.json` file beside itself, which the test reads with
`reference_oracle` (see `test/reference_species.jl`). They used to print a block
of Julia for a human to paste into a test file, and that was wrong in three
ways:

  - pasting is a **manual transcription**, which is the error class these
    generators exist to remove;
  - a data file emitted as Julia has to be written in the *formatter's* dialect.
    `phreeqc_hfo_surface.py` carried a `jl()` whose only job was to rewrite
    `6.0e-09` as `6.0e-9` and `5e-06` as `5.0e-6`, because Runic gates CI and
    Python renders neither the way Julia does;
  - provenance written as `#` comments is readable by nobody and assertable by
    nothing.

Regenerate a fixture by running its script; each names itself in the file's own
`generator` field. The JSON is committed, the generators are not run in CI —
running them needs PHREEQC, Reaktoro or GEMS, which CI does not have.

## What is vendored here, and what deliberately is not

`phreeqc.dat` **is** committed, unmodified, with the USGS User Rights Notice
beside it — see `PHREEQC-PROVENANCE.md`. That notice covers data as well as
code, and requires the notice to travel with any copy.

`sit.dat` is **not**, and the difference is not an oversight. It is the ANDRA/RWM
*ThermoChimie-TDB* compilation, redistributed with PHREEQC but not USGS-authored,
so the notice above says nothing about its terms. `phreeqc_sit.py` therefore
takes `--database` pointing at a copy the caller already has, and writes into its
fixture only the handful of `ε` the comparison uses, with the database named and
hashed. A test stays runnable in CI without a compilation being redistributed.

## Provenance

Every generator reads its own versions at run time and writes them into the
fixture as fields — never a version assumed from a package list. Databases are
identified by md5, because two files with one name are the usual way a
cross-code comparison quietly stops comparing.

The ClaySor 2023 deposit is `doi:10.5281/zenodo.15095062`, CC-BY-4.0, and its
archives are **not** committed here. Fetch them into a scratch directory — never
into the working tree — and check the sums before using anything from them:

```sh
curl -LO https://zenodo.org/records/15095062/files/PHREEQC-ClaySor2023.zip
curl -LO https://zenodo.org/records/15095062/files/GEMS-ClaySor2023.zip
md5sum PHREEQC-ClaySor2023.zip GEMS-ClaySor2023.zip
```

```
0ed43657662b71d03d34b5903e7f6bdb  PHREEQC-ClaySor2023.zip
22af072618380f5bc0086c2836ae8954  GEMS-ClaySor2023.zip
```

The model file the `phreeqc_claysor.py` generator reads out of the first archive
is recorded in its fixture by sha256,
`6df84ee9cdbcf49ba7a2310a7224e7bf5b448f49de17e4dd80b1c5c4c9d7b387`, so a
re-fetch that differs is caught rather than compared.

`sit.dat` is the same kind of thing and is likewise **not** committed. It is the
ANDRA/RWM ThermoChimie-TDB in PHREEQC format, and its own header identifies the
copy this fixture was made from:

```
# Thermodynamic database ANDRA/RWM - THERMOCHIMIE-TDB (www.thermochimie-tdb.com)
# Version 9b0
# BDD Date: 10/8/2015
```

It comes from the ThermoChimie project, `www.thermochimie-tdb.com`, whose terms
are the project's own. `phreeqc_sit.py` defaults to looking for it beside
`phreeqpython`'s other databases and **it is not shipped there**, so point the
generator at wherever the file was put:

```sh
python3 test/reference/phreeqc_sit.py --database /path/to/sit.dat
sha256sum /path/to/sit.dat   # 427d6114ed3f3135054683882319a0852b593d2685f0ccc80855f10ae1c4b840
```

The generator writes that sha256 into the fixture, and the three `ε` the
comparison actually uses travel in the fixture with it — so `test/sit.jl` runs
in CI without the compilation being present at all, which is the point.

Attribution goes to Marinich et al. (2025), `doi:10.1016/j.apgeochem.2025.106510`,
wherever its data is used.

## The scripts

| script | env | produces |
|:--|:--|:--|
| `reaktoro_calcite_co2.py` | `reaktoro-env` | `reaktoro_calcite_co2.json`, read by `test/equilibrium_reference.jl`: calcite, CO₂ and water, amounts and their sensitivities to the CO₂ added |
| `reaktoro_coupling.py` | `reaktoro-env` | `reaktoro_coupling.json`, read by `test/coupling_reference.jl`: the aqueous partition along a constant-rate calcite dissolution |
| `phreeqc_hfo_surface.py` | `mpcm-oracles` | a Zn sorption edge on hydrous ferric oxide, strong and weak sites, `-no_edl`; `--edl` switches to the diffuse layer, which is a **different model** |
| `phreeqc_evolving_surface.py` | `mpcm-oracles` | a sorbent that DISSOLVES: a `SURFACE` coupled to an `EQUILIBRIUM_PHASES` mineral, titrated to exhaustion. Captures the coupling law, not a shared surface model |
| `gems_bench.py` | `mpcm-oracles` | the GEMS3K status report above |

### A trap `phreeqc_hfo_surface.py` exists to not fall into twice

Imposing a pH in PHREEQC is not what it looks like. `pH 4.0 charge` asks the
program to *compute* the pH from the charge balance and treats the number as a
starting guess, so a sweep written that way returns one equilibrium six times,
silently. Writing `pH 4.0` alone sets the solution before the surface is brought
into contact with it, and protolysis then moves it — a requested 4.0 came back
as 7.08. Holding a pH through a surface equilibration needs `Fix_H+` with a
titrant reservoir, which is the experiment the parameters were fitted to.

Both mistakes produce a plausible table. The script therefore **asserts that the
pH it got is the pH it asked for**, and that assertion is what caught them.
