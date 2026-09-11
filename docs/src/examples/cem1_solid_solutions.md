# A CEM I at equilibrium, with every solid solution declared

A thermodynamic model of a cement is only as good as the list of phases it is
allowed to form. Leave one out and the calculation cannot report it; put one in
that does not belong and it can take mass that belongs elsewhere. The list is
therefore a **modeling decision**, and the honest way to make it is to declare
everything the database defines and let the minimization decide.

This page does exactly that on a CEM I at `w/c = 0.5`, from the oxide analysis
on the datasheet through the Bogue conversion to the anhydrous charge, and then
to the hydrated assemblage. It also measures what the alternatives cost, because
the usual practice — declaring the two or three solid solutions one expects to
find — turns out to be the risky one.

```@example ss
using ChemistryLab
using DynamicQuantities
using OptimaSolver
using OrderedCollections
using Printf
using Plots
default(fontfamily = "Computer Modern", framestyle = :box, grid = false)
nothing # hide
```

## 1. The datasheet, and the Bogue conversion

The input is what a cement works publishes: an oxide analysis, in grams per
hundred grams. Nothing else about the cement is assumed.

```@example ss
oxides = OrderedDict(
    "CaO" => 65.03, "SiO2" => 21.4, "Al2O3" => 3.84, "Fe2O3" => 4.49,
    "MgO" => 1.0, "K2O" => 0.46, "Na2O" => 0.13, "SO3" => 2.3, "CO2" => 0.0,
)
water_g = 50.0                      # w/c = 0.5
sum(values(oxides))                 # the analysis does not close at 100 g
```

It sums to 98.65 g, and the missing 1.35 g is loss on ignition and minor oxides
the sheet does not report. Scaling the analysis to 100 g is a choice, and it
moves every element by 1.37 %; it is made here so that "100 g of cement" means
100 g of the oxides that are modeled.

The database is loaded now rather than later, because the molar masses used
below are read from it. None is written down here: an oxide formula goes through
`Species`, which computes the mass from the formula and the element data, while a
cement phase is a database species read by name — `Species("C3S")` would read
that name as the *formula* C₃S, three carbons and a sulfur.

```@example ss
substances = build_species(datapath("cemdata18-thermofun.json"); verbose = false)
byname = Dict(symbol(s) => s for s in substances)

oxide_mass(ox) = ustrip(us"g/mol", Species(ox)[:M])
phase_mass(name) = ustrip(us"g/mol", byname[name][:M])
nothing # hide
```

The Bogue conversion turns those four main oxides into the four clinker phases.
It is not reimplemented here: [Bogue Calculation](@ref) obtains it by inverting
the **mass stoichiometric matrix** of the clinker phases over the oxide
components, which is built from the species data, so no molar mass and no
tabulated coefficient appears anywhere.

```@example ss
clinker = CemSpecies.(split("C3S C2S C3A C4AF"))
oxide4 = CemSpecies.(split("C S A F"))
bogue = inv(mass_matrix(CanonicalStoichMatrix(clinker)).A)

for (i, sp) in enumerate(clinker)
    @printf("  %-5s = %+8.4f CaO %+8.4f SiO2 %+8.4f Al2O3 %+8.4f Fe2O3\n",
            symbol(sp), bogue[i, 1], bogue[i, 2], bogue[i, 3], bogue[i, 4])
end
```

Those are the classical Bogue coefficients, recovered rather than quoted. The
lime combined with the sulfate is not available to the silicates, so it comes
off before the inversion:

```@example ss
f = 100 / sum(values(oxides))

CaO_free = oxides["CaO"] -
    oxides["SO3"] * oxide_mass("CaO") / oxide_mass("SO3") -
    oxides["CO2"] * oxide_mass("CaO") / oxide_mass("CO2")
grams = bogue * [CaO_free, oxides["SiO2"], oxides["Al2O3"], oxides["Fe2O3"]]

charge = OrderedDict{String, Float64}(
    symbol(sp) => grams[i] * f / phase_mass(symbol(sp))
        for (i, sp) in enumerate(clinker)
)
# CEMDATA18 has no periclase and no anhydrous SO3, so the magnesia enters as
# brucite and the sulfate as gypsum.  The water they bring is part of the cement
# and comes off the mixing water.
n_ox = Dict(k => v * f / oxide_mass(k) for (k, v) in oxides)
charge["Gp"] = n_ox["SO3"]
charge["Brc"] = n_ox["MgO"]
charge["K2O"] = n_ox["K2O"]
charge["Na2O"] = n_ox["Na2O"]
free_water = water_g / oxide_mass("H2O") - 2n_ox["SO3"] - n_ox["MgO"]

for (k, v) in charge
    @printf("  %-5s %10.6f mol   (%6.2f g)\n", k, v, v * phase_mass(k))
end
@printf("  %-5s %10.6f mol\n", "H2O@", free_water)
```

```@example ss
names_c = collect(keys(charge))
mass_c = [charge[k] * phase_mass(k) for k in names_c]
bar(names_c, mass_c;
    legend = false, ylabel = "g per 100 g of oxides",
    title = "The anhydrous charge, from the oxide analysis",
    color = :steelblue, size = (720, 380), bottom_margin = 6Plots.mm,
    left_margin = 6Plots.mm)
```

## 2. The phase list, and how many solid solutions there really are

The pure phases are the CEMDATA18 list a Portland cement is given. The solid
solutions are **every multi-end-member phase the database defines** — not a
selection.

```@example ss
pure = split(
    "AlOHmic Kln Gr C12A7 C2S C3A C3S C4AF CA CA2 C2AH7.5 C3AH6 CAH10 " *
        "monosulphate10.5 monosulphate12 monosulphate14 monosulphate16 " *
        "monosulphate9 chabazite zeoliteP_Ca straetlingite5.5 monocarbonate9 " *
        "hemicarbonat10.5 hemicarbonate hemicarbonate9 monocarbonate " *
        "ettringite13 ettringite9 Arg Cal C3FH6 C4FH13 C3FS1.34H3.32 " *
        "Fe-hemicarbonate Femonocarbonate Dis-Dol Ord-Dol Lim Portlandite " *
        "Anh Gp hemihydrate Fe Sd Mag FeOOHmic Py Tro Melanterite K2SO4 " *
        "syngenite K2O hydrotalcite Mgs Brc Na2SO4 natrolite zeoliteX " *
        "zeoliteY Na2O Sulfur Amor-Sl"
)

solutions = [
    "CSHQ" => ["CSHQ-JenD", "CSHQ-JenH", "CSHQ-TobD", "CSHQ-TobH",
        "KSiOH", "NaSiOH"],
    "C3(AF)S0.84H" => ["C3AFS0.84H4.32", "C3FS0.84H4.32"],
    "AFt_SO4" => ["ettringite", "ettringite30"],
    "AFt_SO4_CO3" => ["tricarboalu03", "ettringite03_ss"],
    "AFm_SO4_OH" => ["C4AH13", "monosulphate12"],
    "straetlingite" => ["straetlingite", "straetlingite7"],
    "hydrotalc-pyro" => ["Mg3AlC0.5OH", "Mg3FeC0.5OH"],
    "MSH" => ["M075SH", "M15SH"],
]
length(solutions), length(reduce(vcat, last.(solutions)))
```

Eight phases over twenty end-members. A reader coming from GEM-Selektor will
count **ten** phases in its CEMDATA18 setup, and the difference is worth
stating because it is a difference of formulation and not of chemistry.

!!! note "Why eight here and ten there"
    Two of the ten repeat their end-members with the roles exchanged:
    `SO4_OH_AFm` and `OH_SO4_AFm` are both `C4AH13` + `monosulphate12`, and
    `SO4_CO3_AFt` and `CO3_SO4_AFt` are both `tricarboalu03` +
    `ettringite03_ss`. Declaring the same binary twice is how that formulation
    opens a **miscibility gap**: two immiscible compositions of one binary can
    then coexist as two phases. It requires a species to belong to two phases at
    once.

    Here — and in several other codes — the composition vector has one entry per
    species, so a species belongs to exactly one phase and the gap cannot be
    expressed. The eight distinct phases are the whole of the chemistry; what is
    given up is the ability to represent two coexisting compositions of the same
    binary, which matters only where such a gap actually opens.

Note also that `AFt_SO4_CO3` describes the same aluminate sulfate as `AFt_SO4`
in a different normalization: `ettringite03_ss` is ettringite divided by three.
Both are declared, because deciding in advance which normalization the answer
will use is exactly the kind of prior knowledge this page is trying to do
without.

```@example ss
members = reduce(vcat, last.(solutions))

species = speciation(
    substances, vcat(pure, members, ["CO2@"]); aggregate_state = [AS_AQUEOUS]
)
cs = ChemicalSystem(
    species, CEMDATA_PRIMARIES;
    solid_solutions = [
        SolidSolutionPhase(nm, [byname[m] for m in mem]) for (nm, mem) in solutions
    ],
)
length(cs.species), length(cs.solid_solutions)
```

## 3. The solve

The composition is the anhydrous charge and the mixing water; the element vector
is fixed once from it and passed to the solve, so whatever starting point is
used can only decide whether the optimum is reached, never which optimum it is.

```@example ss
state = ChemicalState(cs)
for (sym, n) in charge
    n > 0 && set_quantity!(state, sym, n * u"mol")
end
set_quantity!(state, "H2O@", free_water * u"mol")
set_quantity!(state, "CO2@", 1.0e-9u"mol")
b = Float64.(cs.SM.A) * ustrip.(us"mol", state.n)

# The activity model of a cement pore solution: CEMDATA18 carries no ion-size
# parameter, so the Debye-Hückel limiting law with the non-ideality in the B-dot
# term and none of it on the neutral species.
model = HKFActivityModel(å = 0.0, Ḃ = 0.097637, Kₙ = 0.0)
eq, cert = equilibrate_certified(state; model = model, b = b)

@printf("certificate: optimal=%s  worst supersaturation=%.3e  balance=%.1e\n",
        cert.optimal, cert.worst_supersaturation, cert.balance)
@printf("pH = %.4f     total volume = %.3f cm3\n",
        pH(eq, model), ustrip(uconvert(us"cm^3", volume(eq).total)))
```

`optimal = true` is a proof, not a report of successful iteration: the Gibbs
minimization is convex when the mixing terms are, so the KKT conditions are
sufficient and [`optimality_certificate`](@ref) decides. A negative worst
supersaturation says every absent phase is undersaturated at the answer — which
is the statement that the phase list was searched, not merely used.

```@example ss
n = ustrip.(us"mol", eq.n)
present = [
    (symbol(cs.species[i]), n[i]) for i in cs.idx_crystal if n[i] > 1.0e-8
]
sort!(present; by = last, rev = true)
for (s, x) in present
    @printf("  %-18s %10.6f mol\n", s, x)
end
```

```@example ss
bar([s for (s, _) in present], [x for (_, x) in present];
    legend = false, ylabel = "mol per 100 g of oxides", yscale = :log10,
    title = "The hydrated assemblage", color = :seagreen,
    xrotation = 45, size = (760, 420),
    bottom_margin = 14Plots.mm, left_margin = 8Plots.mm)
```

## 4. Which solid solutions the calculation keeps, and why

Three of the eight carry mass. The other five are absent — and the point of this
page is that **the calculation says so**, rather than being told.

The criterion is exact for an ideal solid solution and worth deriving, because
it is what makes "declare everything" a sound procedure rather than a hopeful
one. Let the phase have end-members ``i`` with standard chemical potentials
``\mu_i^0``, and let ``\nu_i = \sum_c \lambda_c A_{ci}`` be the potential the
components make available to end-member ``i`` at the answer. If the phase is
present at mole fractions ``x_i``, then

```math
\mu_i^0 + RT \ln x_i = \nu_i
\qquad\Longleftrightarrow\qquad
x_i = 10^{\mathrm{SI}_i},
\qquad
\mathrm{SI}_i = \frac{\nu_i - \mu_i^0}{RT\ln 10},
```

where ``\mathrm{SI}_i`` is the saturation index of the **pure** end-member. A
composition must satisfy ``\sum_i x_i = 1``, so the phase is exactly saturated
when

```math
\Omega = \sum_i 10^{\mathrm{SI}_i} = 1 ,
```

supersaturated above and undersaturated below. Two things follow. The first is
that ``\Omega`` is computable from a single equilibrium, with no iteration over
phase lists. The second is a bound: since ``\Omega \le k \max_i
10^{\mathrm{SI}_i}`` for ``k`` end-members, a phase whose best end-member is
more than ``\log_{10} k`` below saturation cannot form — ``0.301`` for a binary,
whatever its mixing.

[`saturation_indices`](@ref) reports the index of an end-member **in** its
phase, which is ``\mathrm{SI}_i - \log_{10} x_i``; adding ``\ln a_i / \ln 10``
back removes the mole fraction, and with it the ``0/0`` an absent phase would
otherwise produce.

```@example ss
si = saturation_indices(eq, model)
lna = log_activities(eq, model)
inv_ln10 = inv(log(10))

Ω = OrderedDict(
    name(ss) => sum(
        10^(si[symbol(cs.species[i])] + lna[symbol(cs.species[i])] * inv_ln10)
            for i in grp
    )
        for (grp, ss) in zip(cs.ss_groups, cs.solid_solutions)
)
for (nm, w) in Ω
    @printf("  %-16s log10 Omega = %+8.3f   %s\n", nm, log10(w),
            abs(log10(w)) < 1.0e-6 ? "saturated, present" :
            log10(w) > 0 ? "supersaturated" : "undersaturated")
end
```

The three present phases come out at exactly zero, which is the check on the
criterion itself: a phase that the minimization put in must be saturated, and
this formula, derived independently of the solve, says it is.

```@example ss
nms = collect(keys(Ω))
vals = [log10(w) for w in values(Ω)]
bar(nms, vals;
    legend = false, ylabel = "log₁₀ Ω",
    title = "How far each solid solution is from forming",
    color = [v > -1.0e-6 ? :seagreen : :steelblue for v in vals],
    xrotation = 30, size = (760, 420),
    bottom_margin = 12Plots.mm, left_margin = 8Plots.mm)
hline!([0.0]; color = :black, linewidth = 1.5, label = "")
```

One number on that chart deserves attention. `AFt_SO4_CO3` sits at
**−0.03**: thirty thousandths of a log unit from forming. A phase list chosen by
expectation would have dropped it as "obviously absent", and on a cement only
slightly richer in carbonate it would not be. That is the argument for declaring
everything, in one number.

## 5. What declaring only some of them costs

The common practice is to declare the solid solutions one expects to find and
leave the rest out. There are two ways to do that, and they behave very
differently.

**Leaving the other end-members in the system as pure phases.** This is the
variant that looks safest — "their saturation is still tested" — and it is the
one that fails.

```@example ss
three = solutions[1:3]
cs_pure = ChemicalSystem(
    species, CEMDATA_PRIMARIES;                # the same species list
    solid_solutions = [
        SolidSolutionPhase(nm, [byname[m] for m in mem]) for (nm, mem) in three
    ],
)
st = ChemicalState(cs_pure)
for (sym, x) in charge
    x > 0 && set_quantity!(st, sym, x * u"mol")
end
set_quantity!(st, "H2O@", free_water * u"mol")
set_quantity!(st, "CO2@", 1.0e-9u"mol")
b2 = Float64.(cs_pure.SM.A) * ustrip.(us"mol", st.n)
eq2, cert2 = equilibrate_certified(st; model = model, b = b2)

@printf("optimal=%s  worst supersaturation=%+.3f  pH=%.4f  V=%.3f cm3\n",
        cert2.optimal, cert2.worst_supersaturation, pH(eq2, model),
        ustrip(uconvert(us"cm^3", volume(eq2).total)))
```

The certificate refuses, and the answer is wrong by a quarter of the volume. The
reason is visible in the assemblage:

```@example ss
n2 = ustrip.(us"mol", eq2.n)
bad = [(symbol(cs_pure.species[i]), n2[i])
       for i in cs_pure.idx_crystal if n2[i] > 1.0e-8]
sort!(bad; by = last, rev = true)
for (s, x) in bad
    @printf("  %-18s %10.6f mol\n", s, x)
end
```

`ettringite03_ss` and `C4AH13` have taken mass. Offered as a pure phase,
`ettringite03_ss` — ettringite divided by three — competes with ettringite and
wins, because as a pure solid its activity is one whatever its amount, and the
mixing term that would have held it at a small mole fraction inside
`AFt_SO4_CO3` is not there. **An end-member of a solid solution is not a
candidate pure phase**, and offering it as one is not a conservative choice.

**Keeping the other end-members out of the system entirely.** This is what a
hand-written species list usually does, and it does reproduce the answer:

```@example ss
species3 = speciation(
    substances, vcat(pure, reduce(vcat, last.(three)), ["CO2@"]);
    aggregate_state = [AS_AQUEOUS],
)
cs3 = ChemicalSystem(
    species3, CEMDATA_PRIMARIES;
    solid_solutions = [
        SolidSolutionPhase(nm, [byname[m] for m in mem]) for (nm, mem) in three
    ],
)
st3 = ChemicalState(cs3)
for (sym, x) in charge
    x > 0 && set_quantity!(st3, sym, x * u"mol")
end
set_quantity!(st3, "H2O@", free_water * u"mol")
set_quantity!(st3, "CO2@", 1.0e-9u"mol")
b3 = Float64.(cs3.SM.A) * ustrip.(us"mol", st3.n)
eq3, cert3 = equilibrate_certified(st3; model = model, b = b3)

@printf("optimal=%s  worst supersaturation=%+.3f  pH=%.4f  V=%.3f cm3\n",
        cert3.optimal, cert3.worst_supersaturation, pH(eq3, model),
        ustrip(uconvert(us"cm^3", volume(eq3).total)))
```

Same pH, same volume, same assemblage as the eight-phase run — but the
certificate is a **weaker statement**. Its worst supersaturation is far from
zero because the five phases it never declared are not among those it tested;
the eight-phase certificate, whose margin is only `0.03`, is the one that
actually proves the five are out. A comfortable certificate over a small phase
list is not stronger evidence than a tight one over a complete list.

```@example ss
labels = ["all eight\ndeclared", "three declared,\nothers pure",
          "three declared,\nothers absent"]
vols = [ustrip(uconvert(us"cm^3", volume(e).total)) for e in (eq, eq2, eq3)]
phs = [pH(e, model) for e in (eq, eq2, eq3)]
ok = [cert.optimal, cert2.optimal, cert3.optimal]

p1 = bar(labels, vols; legend = false, ylabel = "total volume (cm³)",
    color = [o ? :seagreen : :firebrick for o in ok], title = "Volume")
p2 = bar(labels, phs; legend = false, ylabel = "pH",
    color = [o ? :seagreen : :firebrick for o in ok], title = "pH",
    ylims = (12.0, 13.4))
plot(p1, p2; layout = (1, 2), size = (860, 400),
    bottom_margin = 12Plots.mm, left_margin = 8Plots.mm)
```

Green is a certified answer, red is not.

## 6. A second opinion

[Reaktoro](https://reaktoro.org) solves the same problem from the same
literature ([Leal2017](@cite)) and reads the same ThermoFun file this package
ships, which makes it a useful independent check — see
[Validation against Reaktoro](@ref) for the small case where every knob is
matched and the agreement is quantified species by species.

On this cement the cross-check was run with
`scripts/crosscheck/cem1_solid_solutions_reaktoro.py`, over the same database,
the same species, the same ideal mixing and the same Debye-Hückel convention,
and — the point of the exercise — the **same element vector**, transferred as
moles rather than as grams of oxide so that the two codes do not each convert
the analysis with their own atomic masses. The element vectors agree to eight
decimal places, so what follows is a comparison of chemistry.

| | ChemistryLab | Reaktoro |
|:--|--:|--:|
| solid solutions declared | 8 | 3 |
| pH | 13.0994 | 13.1425 |
| total volume | 74.190 cm³ | 74.250 cm³ |

The two agree to 0.043 units of pH and 0.08 % of volume, on the same twelve
solid phases, with the majors within a few percent.

Two observations about the run itself, both of which are properties of the
problem rather than of either code:

- **The iteration count matters.** A full cement assemblage needs several
  hundred iterations, above the default cap, and a run stopped at the cap
  returns an intermediate iterate rather than an answer. Raising
  `EquilibriumOptions.optima.maxiters` is what makes the comparison possible at
  all; the converged run takes 345 iterations.
- **The eight-phase configuration is genuinely hard.** A solid solution whose
  end-members are *all* at zero has no mole fractions, so its ideal-mixing term
  is undefined there and its gradient depends on the direction of approach — a
  cone singularity at the origin of that phase's subspace. Every code has to
  regularize it somehow. Here the regularization is the `ϵ` floor of the
  activity closure combined with the multi-start of
  [`equilibrate_certified`](@ref), which reaches and certifies the answer; in
  the Reaktoro run the single Newton path did not converge within 4000
  iterations for any of the floors tried, so its numbers above are from the
  three-phase configuration, which is its documented setup for this system.

## 7. What this does not settle

**The mixing is ideal, and that is an assumption.** All eight phases are
declared with ideal mixing here. CEMDATA18 documents non-ideal parameters for
some of them; they are not used, because a mixing parameter one cannot source is
worse than an ideal model honestly labeled. A non-ideal phase changes ``\Omega``
— the criterion of §4 becomes a fixed point, ``x_i = 10^{\mathrm{SI}_i} /
\gamma_i(x)`` — and can open a miscibility gap that this formulation cannot
represent at all.

**The five absent phases are absent on *this* cement.** `AFt_SO4_CO3` at −0.03
is the warning: a different sulfate or carbonate content moves it across. The
procedure transfers; the result does not.

**No kinetics.** This is the assemblage the cement would reach given unlimited
time and complete reaction. What a paste actually reaches, and why it stops
short, is [The hydrating paste, end to end](@ref sec-coupled-hydration) and
[Self-desiccation](@ref sec-self-desiccation).
