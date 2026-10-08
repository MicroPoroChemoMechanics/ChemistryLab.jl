# [Where the numbers come from](@id sec-manual-numbers)

!!! info "Before this page"
    [Getting started](@ref sec-quickstart). This is the first page of the
    Manual, and it assumes nothing else.

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
`FARADAY`/`FARADAY_Q` and `VACUUM_PERMITTIVITY`/`VACUUM_PERMITTIVITY_Q`, with
`P_STANDARD`/`P_STANDARD_Q`, the standard-state pressure of 1 bar, which is a
convention rather than a measured constant.
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
against another code: matching the water basis of Reaktoro [Leal2017](@cite)
exactly is what took one of this package's cross-code comparisons from
5.5 × 10⁻⁵ agreement to 4.1 × 10⁻⁹. Molality is *per kilogram of solvent*, so
the whole calculation inherits whatever you chose here.

## Standard Gibbs energies and equilibrium constants

These come from a thermodynamic database, which `datapath` obtains from its
publisher and `build_species` reads.

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
agree — `OH-` is `-157297.2 J/mol` in slop98 [Johnson1992](@cite) and
`-157270.2` in CEMDATA18 [Lothenbach2019](@cite), because they were fitted in
different contexts. **Do not mix two databases in one calculation** unless you
know why you are doing it; see [Importing thermodynamic databases](@ref
sec-importing-databases) for what each one covers.

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
    datasheet is the source. Label it where it is used — see below for the type
    that does the labeling.
  - **A published value you are comparing against.** A number from a paper or
    from another code is evidence, not a parameter. This package keeps those in
    `test/reference/*.json`, generated by the scripts beside them, carrying the
    code version and the database checksum as fields — so that a comparison can
    say *which* version of *which* code it matched.

And in every one of these cases, write **where it came from** next to it. A
number with a provenance is data; a number without one is a guess that someone
will have to re-derive.

## [Published values live in data files](@id man-literature-data)

A value the package takes from an article (a fitted parameter, a measured
coefficient, a composition, a table) is kept in `data/literature/<key>.json`, one
file per source, the key being the entry of the bibliography. The file names the
place in the source the value was read from, says how it was transcribed and
checked, and gives its unit and, when the source states one, its uncertainty. The
code reads it rather than repeating it:

```@example numbers
using ChemistryLab
record = literature("Powers1948")
```

```@example numbers
record["w_c_sealed"], literature_value("Powers1948", "w_c_sealed")
```

Each quantity comes back as a [`Traced`](@ref) value, so its source and its kind
travel with it, and [`literature_value`](@ref) drops them at the point where the
number enters a calculation. [`literature_table`](@ref) returns a whole table as
columns carrying their units, and [`literature_table_info`](@ref) what the file
says about it: its kind, its location in the source and, for values read off a
figure rather than printed as numbers, how they were read and the relative error
the reading adds. The files follow the format described by
[`literature`](@ref), which checks them as it reads them; the test suite also
checks that every key is an entry of the bibliography and repeats its DOI.

Coefficients that define an equation of state or a published model, such as the
constants of the water equation of state or of the HKF model
[Helgeson1981](@cite), are part of the model rather than data taken from it, and
they remain in the code with their source in a comment.

## When a number has to carry its own history

Everything above is about *finding* a number. This is about not losing what you
found out about it.

`data/pitzer-reardon1990.toml` already does this in data. It carries
`origin = "estimated:<analog>"` on every coefficient [Reardon1990](@citet) had
to borrow, and says why in its own header: *"a caller reading beta0 for Al(OH)4-
would otherwise have no way of knowing that no aluminate solution was ever
measured to obtain it"*. [`Traced`](@ref) is that idea in a type, so it survives
out of the file and into a table, a figure or a fitted result.

```@example numbers
using ChemistryLab: value, source     # too generic to export; reach for them by name

ε = Traced(0.05, PROV_ESTIMATED, "analog: SO4-2")
value(ε), provenance(ε), is_evidence(ε)
```

Six kinds, ordered from the weakest claim to the strongest:
`PROV_UNSTATED`, `PROV_PLACEHOLDER`, `PROV_ESTIMATED`, `PROV_FITTED`,
`PROV_PUBLISHED`, `PROV_MEASURED`. Two of them are worth a sentence each.

**`PROV_UNSTATED` is the weakest, not the middle.** A number that forgot to say
where it came from must never strengthen a result, so it ranks below a declared
placeholder.

**`PROV_FITTED` is not evidence.** A fitted value may be excellent, and whether
it is depends on the data, the model, and whether the parameter was
*identifiable* from that data at all — three questions a predicate cannot
answer. [`is_evidence`](@ref) is therefore true only for `PROV_MEASURED` and
`PROV_PUBLISHED`.

```@example numbers
weakest(Traced(1.0, PROV_MEASURED, "this work"), Traced(2.0, PROV_ESTIMATED, "an analog"))
```

[`weakest`](@ref) is what a derived quantity can honestly claim about itself.
Nothing propagates it automatically — `Traced` is deliberately **not** a `Real`,
so it cannot flow silently into arithmetic and arrive at the far end with its
history gone. Unwrapping with `ChemistryLab.value` is an act a reader can see.

### Reading a whole table at once

```@example numbers
coefficients = vcat(
    [Traced(0.1i, PROV_PUBLISHED, "a compilation") for i in 1:9],
    [Traced(0.0, PROV_PLACEHOLDER, "pending a measurement")],
)
r = provenance_report(coefficients)
(r.total, r.weakest, r.all_evidence)
```

A table where nine coefficients are published and one is a placeholder is a
different object from one where all ten are published, and **the difference
shows in none of the numbers**. This is what to print beside a result.

[`SITActivityModel`](@ref) is the first model built this way: every `ε` carries
its own standing, a borrowed one does not inherit the compilation's, and
`missing_epsilon_pairs` says which pairs are resting on the convention that an
unlisted coefficient is zero.

## When the number does not exist yet

Some numbers have not been measured. A surface site density for a phase that has
never been titrated, an interaction coefficient for a pair no compilation
carries: they are real quantities and there is no source to take them from.

The answer is not to invent one. It is to make the calculation able to
**identify** it from data — and then to be exact about what the identification
established, which is usually less than the number of parameters that came out
of it.

```@example numbers
# A model with a collinearity put in on purpose: `a` and `c` enter only as
# their product, so no measurement of y can separate them.
t = range(0, 5; length = 40)
decay(p) = @. p[1] * p[3] * exp(-p[2] * t)
θ = [2.0, 0.7, 1.0]

id = identifiability(decay, θ; names = ["a", "b", "c"], observed = decay(θ))
id
```

Three things are worth reading there.

**The spectrum falls off a cliff.** With a noise level `σ`,
[`identifiable_rank`](@ref) counts a direction as determined when its standard
error in `log θ`, `σ` over its singular value, is below one. `σ` is the
instrument's when it is given (`noise`), and otherwise the residual standard
deviation of the fit on its `n − p` degrees of freedom. Here the fit is exact,
the residual is zero, and `σ` takes its floor, `√eps` times the
root-mean-square of the curve: only the direction at the rounding of the
arithmetic is left out. Without a noise level, the rank is read off the largest
*ratio* between consecutive singular values: a threshold would have units, and
a gap has none.

**The empty direction names the trade-off.**

```@example numbers
round.(id.V[:, end]; digits = 3)     # equal and opposite in a and c, nothing in b
```

**And the correlation says it more directly**, which is why it is the instrument
to reach for when two parameters trade off rather than one being invisible.
`scripts/hydration_calibration.jl` found a rate constant and an Avrami exponent
correlated at −0.96 — the data see a product, not its factors — and therefore
fitted one of the two. *Which* one is a modeling judgement and not a statistical
one: it kept the rate constant, because that is the quantity a different clinker
plausibly changes.

### What comes out is a claim, and it says so

```@example numbers
as_traced(id, θ; source = "a synthetic fit")
```

One `fitted` and **two** `placeholder`, and which two is the interesting part.

`a` and `c` enter only as their product, so neither is determined on its own —
both are placeholders — while `b` is determined and is reported as fitted. That
answer comes from [`null_participation`](@ref), which asks how much of each
parameter lies in the directions the data do not constrain:

```@example numbers
null_participation(id)     # a and c, half each; b, none
```

**Not** from the parameter's position relative to the rank. A rank of 2 counts
*directions in parameter space*, and the parameter order is whatever the caller
packed; reading it as "the first two are fine" would have cleared `a` and
flagged `c` — the wrong answer, reached by a plausible route.

!!! warning "A fit is not a mechanism"
    Reproducing a measurement establishes that a model *can* reproduce it.
    Adjusting a site density can absorb a denticity or a lateral interaction and
    still fit, so a parameter constrained independently, and data held out of
    the fit, are what separate a mechanism from a curve that passes through the
    points.

### The same question on a real model

The collinearity above was put in by hand. Here is one nobody put in.

A rate law here scales with the binder's fineness, and passing a
[`ShrinkingCoreArea`](@ref) makes that factor follow the grains as they are
consumed, `(n/n₀)^p`. Is `p` a new parameter, or a rewriting of one the law
already had? Over the range where the shell-formation branch controls, the law
reads `k₃(1-ξ)^{n₃} · (1-ξ)^p = k₃(1-ξ)^{n₃+p}`, so on paper only the **sum** is
visible. Now ask the model rather than the algebra.

```@example numbers
# Remaining fraction n/n₀ over the range where that branch is the active one.
fracs = collect(range(0.8, 0.1; length = 30))
idx = Dict("C3S" => 1)
lna0 = StateView([0.0], idx)
nini = StateView([1.0], idx)

function pk_curve(q)
    pr = merge(PK84_PARAMS_C3S, (k₃ = q[1] * u"1/d", n₃ = q[2]))
    law = parrott_killoh_avrami(
        pr, "C3S";
        blaine = ShrinkingCoreArea(BlaineSurfaceArea(PK_BLAINE_REF); exponent = q[3]),
    )
    return [law(293.15, 1.0e5, 86400.0, StateView([f], idx), lna0, nini) for f in fracs]
end

θpk = [1.1, 3.3, 2 / 3]     # k₃ [1/d], n₃, and the shrinking-core exponent
identifiability(pk_curve, θpk; names = ["k₃", "n₃", "p"])
```

Two of three directions, and the pair named is `n₃` / `p` — the degeneracy the
algebra predicted, found from the model rather than asserted about it.

#### The condition number of an exact degeneracy

The sensitivities are exact, by forward-mode differentiation, so an exact
degeneracy shows as one: the smallest singular value is ten orders of magnitude
below the largest, and the correlation of the pair is −1 to the digits printed.

```@example numbers
idpk = identifiability(pk_curve, θpk; names = ["k₃", "n₃", "p"])
(condition = idpk.condition, r_n3_p = idpk.correlation[2, 3])
```

Until 0.28.2 these sensitivities were central differences. Their default 5 %
step moved `n₃ = 3.3` by 0.165 **in an exponent**, far enough for the
second-order error of the differences to differ between two parameters that are
exactly collinear, and the same model came out with a condition number of about
eighty where it is 1e10: a degeneracy partly hidden by the method that was meant
to find it.

#### What it says

```@example numbers
round.(idpk.V[:, end]; digits = 4)     # the direction the data cannot see
```

Nothing in `k₃` and a trade-off between `n₃` and `p`, as written on paper.

```@example numbers
null_participation(idpk)
```

**Here the two instruments part company, and that is worth understanding.** The
correlation says `n₃` and `p` trade off exactly. The participation does *not*
split the blame evenly, because it is computed in **logarithms of the
parameters**: only `n₃ + p` is visible, so a one-percent change in `p` moves the
curve five times less than a one-percent change in `n₃` — `n₃` is 3.3 and `p` is
2/3. In relative terms `p` really is the less determined of the two, and a
reported uncertainty is a relative statement.

So `as_traced` keeps `k₃` and `n₃` as fitted and marks `p` a placeholder:

```@example numbers
as_traced(idpk, θpk; source = "a synthetic hydration curve")
```

Which is the right advice. Fitting `p` and `n₃` from a single hydration curve is
fitting a sum; fitting `p` with `n₃` held at its published value is a different
and legitimate thing, and what comes out is a `ShrinkingCoreArea` exponent for
that binder, not a geometry.

## A checklist

Before a number goes into a script:

 1. Is it a physical constant? → `R_GAS`, `FARADAY`, `VACUUM_PERMITTIVITY`.
 2. Is it a property of a substance — mass, charge, volume, standard energy? →
    get the species from a database and ask it.
 3. Is it derived from one of those, like moles per kilogram? → compute it from
    the property, not from the number you remember.
 4. Is it the thing you are actually studying, or evidence from elsewhere? →
    write it, and write where it came from beside it.
 5. Is it a value whose standing someone downstream will need to know — a
    placeholder, an estimate, a fitted parameter? → wrap it in [`Traced`](@ref),
    so the claim travels with it.
 6. Does the number not exist yet? → make it identifiable from data, and let
    [`as_traced`](@ref) say which parameters the data actually determined.
 7. None of the above? → it is probably case 2 and you have not found it yet.

## See also

  - [Importing thermodynamic databases](@ref sec-importing-databases) — which databases are read, and where each comes from.
  - [Species](@ref sec-species) — building one from scratch, including
    when you genuinely have to.
  - [Thermodynamic Functions](@ref sec-thermodynamics) — how a standard
    property becomes a function of temperature and pressure.
