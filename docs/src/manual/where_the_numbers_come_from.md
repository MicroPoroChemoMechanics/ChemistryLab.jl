# [Where the numbers come from](@id sec-manual-numbers)

A geochemical calculation is made almost entirely of numbers somebody measured:
gas constants, molar masses, standard Gibbs energies, equilibrium constants. Most
of them are already inside this package, and the single most useful habit when
writing a script with it is this one:

!!! tip "The rule, in one line"
    **If the package can tell you a number, ask it. Do not type it.**

This page shows how to ask, for each kind of number, and is explicit about the
few cases where typing one is the right thing to do.

## Why this matters more than it looks

Typing a number that exists elsewhere is not a style preference. Three things go
wrong, and all three are quiet:

  - **It drifts.** A value recalled from memory or copied from an old script is
    usually *nearly* right. The water in some of this package's own test files
    read `-237181 J/mol` against the database's `-237183.0`; portlandite read
    `-897010` against CEMDATA18's `-897013.0`. Neither gap is large enough to
    fail anything, and both are wrong.
  - **It asserts a unit system silently.** `8.31446261815324` is the gas
    constant *in joules per mole per kelvin*. Nothing in the code says so, and
    nothing checks.
  - **It cannot be traced.** A reviewer reading `-552790.0` has to find out what
    it is and where it came from. A reviewer reading `db["Ca+2"]` does not.

## Physical constants

They live in one place, and each comes in two forms. Use the **dimensional**
one whenever you are working with quantities, and the stripped one inside a
tight loop where a `Quantity` would cost an allocation per evaluation.

```@example numbers
using ChemistryLab, DynamicQuantities

R_GAS_Q     # the molar gas constant, with its dimensions
```

```@example numbers
R_GAS       # the same constant, as a bare Float64 in J/(mol·K)
```

The set currently exported is `R_GAS`/`R_GAS_Q` (gas constant),
`FARADAY`/`FARADAY_Q` and `VACUUM_PERMITTIVITY`/`VACUUM_PERMITTIVITY_Q`.
`ChemistryLab.RT_over_F(T)` gives the Nernst scale `RT/F` in volts; it is not
exported, so it is reached through the module. They are taken from
`DynamicQuantities.Constants`, which carries the CODATA values, so they are not
retyped anywhere — including here:

```@example numbers
FARADAY     # note the digits: 96485.33212 is NOT this number
```

If you need a constant that is missing, add it to `src/utils/constants.jl` in
both forms rather than writing it where you need it.

## Molar masses, and anything else a species knows

A species carries its own properties. Ask it.

```@example numbers
db = Dict(symbol(s) => s for s in
          build_species(datapath("slop98-inorganic-thermofun.json"); verbose = false))
water = db["H2O@"]
water[:M]        # molar mass, with units
```

```@example numbers
charge(water), atoms(water), aggregate_state(water)
```

Never build a table of molar masses by hand. The package computes them from the
formula and the periodic table, and it does so for species you invent as well as
for species it ships.

### "How many moles of water is a kilogram?"

`55.5` is the number everyone writes, and it is not one kilogram — it is
0.99983 kg. Derive it instead:

```@example numbers
n_water = ustrip(us"mol", 1.0u"kg" / water[:M])
```

The difference is 1.7 × 10⁻⁴ relative. That sounds negligible until you compare
against another code: matching Reaktoro's water basis exactly is what took one
of this package's cross-code comparisons from 5.5 × 10⁻⁵ agreement to
4.1 × 10⁻⁹. Molality is *per kilogram of solvent*, so the whole calculation
inherits whatever you chose here.

## Standard Gibbs energies and equilibrium constants

These come from a thermodynamic database. The package ships several; `datapath`
locates them and `build_species` reads them.

```@example numbers
G(s) = ustrip(us"J/mol", s[:ΔₐG⁰](T = 298.15u"K", P = 1.0e5u"Pa"; unit = true))
G(db["Ca+2"]), G(db["Cl-"]), G(db["OH-"])
```

A standard energy is generally a **function of temperature and pressure**, which
is why it is called rather than read. Asking for it at one temperature and using
it at another is a separate mistake, and this form makes it visible.

### Which databases are there, and what is in them?

```@example numbers
using Printf
for f in sort(readdir(dirname(datapath("cemdata18-thermofun.json"))))
    endswith(f, ".json") && @printf("%-40s %6.1f kB\n", f, filesize(datapath(f)) / 1024)
end
```

Each is a different body of work with a different scope, and they do not always
agree — `OH-` is `-157297.2 J/mol` in slop98 and `-157270.2` in CEMDATA18,
because they were fitted in different contexts. **Do not mix two databases in
one calculation** unless you know why you are doing it; see
[Database Interoperability](@ref sec-databases) for what each one covers.

To find out whether a database has what you need:

```@example numbers
c18 = Dict(symbol(s) => s for s in
           build_species(datapath("cemdata18-thermofun.json"); verbose = false))
haskey(c18, "Portlandite"), haskey(c18, "Ca(OH)2")
```

Note the answer. The name is the **database's**, not a spelling you choose, and
looking a phase up under the formula you would have written finds nothing.
[`speciation`](@ref) is the tool for selecting a set of species without
guessing at names one at a time.

## What you *should* type

The rule is not "never write a number". It is "never write a number that already
exists". Three kinds legitimately do not:

  - **The subject of the calculation.** A surface complexation constant under
    test, a rate coefficient being fitted, a capacitance you are exploring: these
    are inputs to the question, and writing them is the point. Where this
    package's tests build a species like `XsOH2+`, its standard energy is
    *constructed* to be the log K under test — the oracle, derived rather than
    looked up — and every such place says so in a comment.
  - **Data with no database.** A ground granulated blast-furnace slag glass has
    no thermodynamic entry; its composition comes from a datasheet, and that
    datasheet is the source. Label it where it is used.
  - **A published value you are comparing against.** A number from a paper or
    from another code is evidence, not a parameter. This package keeps those in
    `test/reference/*.json`, generated by the scripts beside them, carrying the
    code version and the database checksum as fields — so that a comparison can
    say *which* version of *which* code it matched.

And in every one of these cases, write **where it came from** next to it. A
number with a provenance is data; a number without one is a guess that someone
will have to re-derive.

## A checklist

Before a number goes into a script:

 1. Is it a physical constant? → `R_GAS`, `FARADAY`, `VACUUM_PERMITTIVITY`.
 2. Is it a property of a substance — mass, charge, volume, standard energy? →
    get the species from a database and ask it.
 3. Is it derived from one of those, like moles per kilogram? → compute it from
    the property, not from the number you remember.
 4. Is it the thing you are actually studying, or evidence from elsewhere? →
    write it, and write where it came from beside it.
 5. None of the above? → it is probably case 2 and you have not found it yet.

## See also

  - [Database Interoperability](@ref sec-databases) — what each shipped database covers.
  - [Species](@ref sec-species) — building one from scratch, including
    when you genuinely have to.
  - [Thermodynamic Functions](@ref sec-thermodynamics) — how a standard
    property becomes a function of temperature and pressure.
