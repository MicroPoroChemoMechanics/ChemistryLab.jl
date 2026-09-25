# [A composite binder: two glasses at once](@id ex-cem5-composite)

!!! info "Before this page"
    [A blastfurnace cement](@ref ex-cem3-slag) and [A pozzolanic binder](@ref
    ex-cem4-pozzolanic), whose difficulties this page combines.

A CEM V carries **both** a blastfurnace slag and a pozzolana, each between 18 %
and 30 % for a CEM V/A, leaving 40 % to 64 % clinker. It is the binder in which
every difficulty of the preceding pages arrives together:

- two constituents with **no formula**, so the budget comes from two oxide
  analyses rather than from phases ([the CEM III page](@ref ex-cem3-slag));
- **magnesium**, from the slag, which makes hydrotalcite and nothing else;
- **aluminum in excess of what the aluminate phases can hold**, from the fly
  ash, which is why the C-S-H has to be the one that carries it
  ([the CEM IV page](@ref ex-cem4-pozzolanic));
- **sulfur at two oxidation states**, S(-II) from the slag against S(+VI) from
  the clinker's calcium sulfate, which is what makes charge a conserved quantity
  of its own ([the redox chapter](@ref theory-redox)).

None of these is new here. What is new is that they must hold simultaneously,
and the point of the page is that the same element budget answers all four.

## 0. Why two constituents are not one constituent twice

A natural first reaction to a composite binder is that 24 % slag plus 24 % fly
ash should behave like 48 % of "something in between". It does not, and the
reason is worth stating because it is the reason the family exists.

The two glasses differ where it matters most — in **calcium**:

| | CaO | SiO₂ | Al₂O₃ | MgO |
|:--|--:|--:|--:|--:|
| blastfurnace slag | ~41 % | ~36 % | ~11 % | ~8 % |
| siliceous fly ash | ~4 % | ~53 % | ~26 % | ~2 % |

A slag is *latently hydraulic*: it carries nearly as much calcium as it needs and
mostly wants an alkaline trigger. A fly ash is *pozzolanic*: it carries almost
none, so it must **take** calcium from the paste — which means consuming the
portlandite the clinker makes.

Put them together and the slag's calcium partly feeds the ash's appetite, so the
portlandite lasts further than it would with the same mass of ash alone. That is
why a CEM V reaches a replacement level a CEM IV cannot, and it is visible in the
budget table of section 3: read the `Ca+2` row across the three columns.

The aluminum goes the other way. Both glasses bring it, the ash twice as much per
gram, and there are only so many aluminate hydrates to hold it — which is what
makes the C-A-S-H model the load-bearing choice here rather than a refinement.

!!! note "Nothing here reacts completely, and the page says how much does"
    A Gibbs minimization reacts whatever it is given. What it is given is
    therefore a modeling decision, not a detail: hand it the whole binder and it
    answers *what this paste becomes if everything reacts*, which no paste does —
    not the glasses, and not the clinker either. Section 1 fixes a reacted
    fraction for each constituent, from the water available and from published
    28-day measurements, and [section 6](@ref cem5-dor) sweeps it.

```@example cem5
using ChemistryLab
using DynamicQuantities
using OptimaSolver
using OrderedCollections
using Printf
using Plots
default(framestyle = :box, grid = false)

# The ZEOLITE-EXTENDED database, and that is not a detail of convenience.
# [The CEM IV page](@ref ex-cem4-pozzolanic) establishes why: past roughly a
# third replacement the aluminum and the alkalis the pozzolana brings exceed what
# the C-A-S-H and the aluminate hydrates can hold, and with no phase left to
# receive them the minimization has no admissible assemblage at all. A CEM V/A at
# the midpoint of its range is 48 % replaced, well inside that regime.
substances = build_species(datapath("cemdata18-zeolites.json"); verbose = false)
byname = Dict(symbol(s) => s for s in substances)
molar_mass(n) = ustrip(us"g/mol", byname[n][:M])

# Only the zeolites this binder could form. Two of the twenty-eight are a
# chloride and a nitrate sodalite, and this paste carries neither element:
# declaring them would pull every aqueous chloride and nitrate species into the
# system on a budget of exactly zero.
base = build_species(datapath("cemdata18-thermofun.json"); verbose = false)
carries(sp, el) = haskey(atoms(sp), Symbol(el))
ZEOLITES = [z for z in sort(collect(setdiff(Set(symbol.(substances)),
                                            Set(symbol.(base)))))
            if !carries(byname[z], "Cl") && !carries(byname[z], "N")]
nothing # hide
```

## 1. The record, and the three numbers in it

```@example cem5
rec = joinpath(datapath("experimental"),
               "smilauer2025-200-cemV-A-S-V-32.5R-prachovice.csv")
for line in eachline(rec)
    startswith(line, "#") || break
    occursin(r"cement:|blaine:|wb:|temperature:", line) && println(line)
end
```

`(S-V)` in the designation is the composition: **S** for blastfurnace slag, **V**
for siliceous fly ash. The deposit gives the fineness, the water/binder ratio and
the calorimetry. It gives neither the clinker composition nor the two replacement
levels, so those are assumed at the midpoint of the EN 197-1 range and labeled:

```@example cem5
# ASSUMED: a Bogue composition representative of a CEM I clinker.
CLINKER = OrderedDict("C3S" => 0.65, "C2S" => 0.11, "C3A" => 0.11, "C4AF" => 0.08)

# ASSUMED: midpoint of the EN 197-1 CEM V/A range, which is 40-64 % clinker
# with 18-30 % slag and 18-30 % pozzolana.
SLAG_FRACTION = 0.24
ASH_FRACTION = 0.24
GYPSUM = 0.046

# ASSUMED: the same two analyses used on the CEM III and CEM IV pages, so that
# the three calculations differ in their proportions and not in their inputs.
SLAG = OrderedDict("CaO" => 0.41, "SiO2" => 0.36, "Al2O3" => 0.11,
                   "MgO" => 0.08, "SO3" => 0.02)
FLYASH = OrderedDict("SiO2" => 0.53, "Al2O3" => 0.26, "Fe2O3" => 0.07,
                     "CaO" => 0.04, "MgO" => 0.02, "K2O" => 0.025,
                     "Na2O" => 0.008, "SO3" => 0.005)

# THE CLINKER'S OWN ALKALIS, which Bogue does not account for -- minor oxides
# outside the four-phase decomposition, and in a paste what fixes the pH: they
# dissolve almost completely and stay in solution, while the calcium is held down
# by portlandite at 12.5. ASSUMED at a usual industrial level, as a fraction of
# the CLINKER mass. The fly ash brings alkalis of its own, and section 3's table
# shows which constituent brought what. [The CEM III page](@ref cem3-alkali)
# sweeps this input over the industrial range: it, and almost nothing else, is
# what moves the pH.
ALKALIS = OrderedDict("K2O" => 0.008, "Na2O" => 0.002)

WB = 0.40            # MEASURED, from the record above
BINDER_G = 100.0

# HOW MUCH OF EACH CONSTITUENT HAS REACTED. Two ceilings, and the reacted
# fraction is the lower of them.
#
# THE WATER CEILING is the same for every constituent, because it is a property
# of the pore space and not of the grain. Powers (1948): a gram of cement needs
# about 0.42 g of water to hydrate completely -- 0.23 g written into the hydrate
# formulae, and 0.19 g held in the gel pores those hydrates create. Below that
# the paste desiccates itself and stops WITH WATER STILL IN IT, the remaining
# water being in pores too fine to reach an unhydrated grain at a scale orders
# of magnitude larger. Nothing about that argument mentions clinker, so it
# applies to a slag particle and an ash sphere exactly as it does to an alite
# grain. Under water curing the ceiling moves to 0.36, the volume emptied by
# chemical shrinkage being refilled from the bath -- `curing = :saturated`.
CURING = :sealed                                   # ASSUMED: the record is silent
ALPHA_WATER = powers_alpha_max(WB; curing = CURING)

# THE KINETIC CEILING is each constituent's own dissolution rate at the age
# considered, and for the two glasses it is far the lower of the two. The RILEM
# TC 238-SCM round robin [Durdzinski2017](@cite) measured exactly this geometry
# -- Portland cement blended with slag and with a siliceous fly ash, w/b 0.40 --
# in seven laboratories. Its Table 4 at 28 days, by SEM image analysis, the one
# technique the study found consistent: 38 % and 48 % for one slag, 45 % and
# 49 % for the other, 20 % for the siliceous fly ash. The study's own verdict on
# the precision is worth carrying: "at best +/- 5 %".
#
# ASSUMED from that table, and [section 6](@ref cem5-dor) sweeps both.
ALPHA_SLAG = min(0.45, ALPHA_WATER)
ALPHA_ASH = min(0.20, ALPHA_WATER)

# The clinker is taken AT its water ceiling. At 28 days it has not quite got
# there, so the assemblage below is an upper bound on what the clinker
# contributes; section 6 moves it too.
ALPHA_CLINKER = ALPHA_WATER

@printf("water ceiling (%s, w/b = %.2f) : %.3f\n", CURING, WB, ALPHA_WATER)
@printf("reacted: clinker %.0f %%, slag %.0f %%, fly ash %.0f %%\n",
        100ALPHA_CLINKER, 100ALPHA_SLAG, 100ALPHA_ASH)
```

## 2. The system: C-A-S-H, hydrotalcite and the sulfur ladder

```@example cem5
pure = split(
    "C3S C2S C3A C4AF Gp Anh Cal Portlandite ettringite monosulphate12 " *
    "monocarbonate hemicarbonate C4AH13 C3AH6 C3FH6 straetlingite " *
    "hydrotalcite Mg2AlC0.5OH Brc FeOOHmic AlOHmic Amor-Sl Mgs " *
    # Aluminum sinks that CEMDATA18 documents and an earlier version of this
    # list simply did not declare. The siliceous hydrogarnet `C3AS0.84H4.32`
    # is the one that matters most here: it is the ALUMINUM end-member of the
    # family whose iron end-member was already present, and a blended binder
    # puts a great deal of aluminum into it. Leaving it out does not make the
    # calculation conservative -- it makes it insoluble, because the element
    # has to go somewhere.
    "C3AS0.84H4.32 C3AS0.41H5.18 straetlingite7 Gbs AlOHam " *
    "M4A-OH-LDH M6A-OH-LDH M8A-OH-LDH C2AH7.5 C4AH11 C4AH19 " *
    "Tro Py Sulfur K2SO4 syngenite Na2SO4"
)
pure = vcat(String.(pure), ZEOLITES)
gel = ["T2C-CNASHss", "T5C-CNASHss", "TobH-CNASHss",
       "5CA", "5CNA", "INFCA", "INFCN", "INFCNA"]
    # THE DECLARED SOLID SOLUTION CANNOT REACH AN ALUMINUM-RICH COMPOSITION,
    # so the aluminum end-member is declared beside it. The siliceous
    # hydrogarnet is a substitution of Al and Fe(III) on TWO sites:
    #
    #   C3AS0.84H4.32   (AlAlO3)[...]       x(Al) = 1.0
    #   C3AFS0.84H4.32  (AlFe|3|O3)[...]    x(Al) = 0.5
    #   C3FS0.84H4.32   (Fe|3|Fe|3|O3)[...] x(Al) = 0.0
    #
    # CEMDATA18 declares the binary between the middle and the iron end
    # (`data/solid_solutions.toml`, source Lothenbach2019), which spans
    # x(Al) from 0.5 down to 0. A CEM I is iron-rich through its ferrite
    # phase and never needs more. A binder whose pozzolana brings twice as
    # much aluminum as iron does, and the declared phase cannot go there.
    #
    # Declaring the aluminum end-member as a separate pure phase is how that
    # half of the series is reachable at all. It is an approximation, and the
    # approximation is named: as a pure phase it carries no mixing entropy,
    # where a site-fraction model over x(Al) in [0,1] would. Extending the
    # solid solution to three end-members would be WORSE, not better --
    # three compositions of a two-site substitution are not three independent
    # end-members, and an ideal ternary over them gets the configurational
    # entropy wrong.
feal = ["C3AFS0.84H4.32", "C3FS0.84H4.32"]
redox_species = ["HS-", "H2S@", "SO4-2", "SO3-2", "S2O3-2", "O2@", "H2@", "CO2@"]

species = speciation(substances, vcat(pure, gel, feal, redox_species);
                     aggregate_state = [AS_AQUEOUS])
ss = [SolidSolutionPhase("CNASH_ss", [byname[m] for m in gel]),
      SolidSolutionPhase("C3(AF)S0.84H", [byname[m] for m in feal])]
cs = ChemicalSystem(species, CEMDATA_PRIMARIES; solid_solutions = ss)
model = HKFActivityModel(å = 0.0, Ḃ = 0.097637, Kₙ = 0.0)

components = String.(symbol.(cs.SM.primaries))
@printf("%d species, %d components: %s\n",
        length(cs.species), length(components), join(components, " "))
@printf("charge (`Zz`) kept as a conservation component: %s\n",
        "Zz" in components ? "yes" : "no")
```

## 3. The budget is additive

Three contributions, and the only reason they can be added is that they are
expressed on the **same components**: the clinker through its phases, each glass
through [`oxide_budget`](@ref).

```@example cem5
clinker_frac = 1 - SLAG_FRACTION - ASH_FRACTION - GYPSUM

"""
    paste(α_slag, α_ash; α_clinker) -> (; state, clinker, slag, ash, total)

The fresh state and the three element contributions, at given reacted fractions.

The state carries ONLY what reacts. What does not is still in the specimen --
unhydrated clinker cores and undissolved glass, both of them measurable -- but it
is no part of the minimization, and putting it in would be a different and false
statement: that an unreacted grain is at equilibrium with the pore solution it is
sitting in.

Two things carry no ceiling. The calcium sulfate is soluble and is gone within
hours, ceiling or no ceiling. And ALL of the mixing water enters: the ceiling
limits how far the reaction can go, not how much water was poured in -- the water
that cannot reach a grain is still in the specimen, and still in the balance.
"""
function paste(α_slag, α_ash; α_clinker = ALPHA_CLINKER)
    st = ChemicalState(cs)
    for (phase, frac) in CLINKER
        set_quantity!(st, phase,
            α_clinker * BINDER_G * clinker_frac * frac / molar_mass(phase) * u"mol")
    end
    set_quantity!(st, "Gp", BINDER_G * GYPSUM / molar_mass("Gp") * u"mol")
    set_quantity!(st, "H2O@", BINDER_G * WB / molar_mass("H2O@") * u"mol")

    clinker = Float64.(cs.SM.A) * ustrip.(us"mol", st.n)
    # The alkalis leave the grain as it dissolves: same fraction as the clinker.
    clinker .+= oxide_budget(ALKALIS, cs.SM.primaries;
                             mass = BINDER_G * clinker_frac * α_clinker * u"g")
    slag = oxide_budget(SLAG, cs.SM.primaries;
                        mass = BINDER_G * SLAG_FRACTION * α_slag * u"g")
    ash = oxide_budget(FLYASH, cs.SM.primaries;
                       mass = BINDER_G * ASH_FRACTION * α_ash * u"g")
    return (; state = st, clinker, slag, ash, total = clinker .+ slag .+ ash)
end

"""
    solve_paste(α_slag, α_ash; start) -> (state, certificate)

The certified equilibrium of the paste at given reacted fractions.

`start` is where the search begins. A cement equilibrium is **hard to start
cold** -- 135 species, an assemblage that is not known in advance, and a pore
solution four orders of magnitude more dilute than the solids -- and the standard
remedy is to walk to it from a state that is easier. Here the walk is in the
reacted fraction itself: a younger paste has released less of everything, so it
is a smaller perturbation of pure water, and its answer is a good start for an
older one.

That is safe here for a reason that is **checked rather than assumed**. Both
solid solutions of section 2 carry the default ideal mixing model, and
`SolidSolutionPhase` refuses a model whose mixing energy has a spinodal — so the
Gibbs function is convex, its minimum is unique, and a continuation **cannot
change what is found**, only whether the search finds it.

Waive that refusal with `check_convexity = false` and none of it holds: inside a
spinodal the minimum is two coexisting compositions rather than one, the
certificate loses the sufficiency that rests on convexity, and the starting point
would then decide which branch the answer lands on — which is exactly the
path-dependence [the miscibility gap page](@ref ex-miscibility-gap) is about. The
certificate decides at every step here, and a start is reused only after it has
been certified.
"""
function solve_paste(α_slag, α_ash; start = nothing)
    pa = paste(α_slag, α_ash)
    return equilibrate_certified(something(start, pa.state); model = model, b = pa.total)
end

p28 = paste(ALPHA_SLAG, ALPHA_ASH)
state, b = p28.state, p28.total

@printf("of 100 g of binder: clinker %.1f g, slag %.1f g, fly ash %.1f g, gypsum %.1f g\n",
        100clinker_frac, 100SLAG_FRACTION, 100ASH_FRACTION, 100GYPSUM)
@printf("of which reacted   : clinker %.1f g, slag %.1f g, fly ash %.1f g, gypsum %.1f g\n\n",
        100clinker_frac * ALPHA_CLINKER, 100SLAG_FRACTION * ALPHA_SLAG,
        100ASH_FRACTION * ALPHA_ASH, 100GYPSUM)
@printf("%-8s %10s %10s %10s %10s\n", "", "clinker", "slag", "fly ash", "total")
for (i, c) in enumerate(components)
    abs(b[i]) > 1.0e-6 &&
        @printf("%-8s %10.5f %10.5f %10.5f %10.5f\n",
                c, p28.clinker[i], p28.slag[i], p28.ash[i], b[i])
end
```

That table is the page in one object. The aluminum comes mostly from the ash,
the magnesium only from the slag, the sulfur from all three at two different
oxidation states, the alkalis from the clinker and the ash together — and the
calcium overwhelmingly from the clinker even at 48 % replacement.

Read the `K+` and `Na+` rows in particular. They are a few hundredths of a mole
against nearly a mole of calcium, and they are what sets the pH: the calcium is
buffered by portlandite at 12.5 and cannot rise, while the alkalis dissolve
almost entirely and stay in solution. A budget that leaves them out — and a Bogue
calculation leaves them out, since they sit outside the four-phase decomposition
— returns a portlandite floor and calls it a pore solution.

## 4. The equilibrium

```@example cem5
# Continued from a seven-day paste rather than started cold -- see `solve_paste`
# above for why, and section 6 for the ages this walk passes through.
eq_early, cert_early = solve_paste(0.35, 0.10)
eq, cert = solve_paste(ALPHA_SLAG, ALPHA_ASH;
                       start = cert_early.optimal ? eq_early : nothing)

@printf("started from a 7-day paste: certified %s\n", cert_early.optimal)
@printf("certificate: optimal=%s  worst SI=%.2e  element balance=%.1e\n",
        cert.optimal, cert.worst_supersaturation, cert.balance)
# The volume is that of the REACTED system and its pore solution. The
# unhydrated clinker and the undissolved glass occupy volume too; the clinker's
# is computable from its phases, the glass's would need a density this page has
# not been given, so neither is added rather than one of them being.
@printf("pH = %.3f   volume of the reacted system = %.2f cm3\n",
        pH(eq, model), ustrip(uconvert(us"cm^3", volume(eq).total)))
```

```@example cem5
n = ustrip.(us"mol", eq.n)
present = sort([(symbol(cs.species[i]), n[i])
                for i in cs.idx_crystal if n[i] > 1.0e-4]; by = last, rev = true)
for (name, amount) in present
    @printf("  %-18s %9.5f mol\n", name, amount)
end
```

Where each element ends up is the more useful reading, since the C-A-S-H holds
several of them at once:

```@example cem5
"""Moles of one conservation component held by each solid phase, largest first."""
function component_in_solids(component; tol = 1.0e-4)
    row = findfirst(==(component), String.(symbol.(cs.SM.primaries)))
    row === nothing && error("$component is not a component of this system")
    A = Float64.(cs.SM.A)
    out = [(symbol(cs.species[i]), A[row, i] * n[i]) for i in cs.idx_crystal]
    return sort(filter(p -> last(p) > tol, out); by = last, rev = true)
end

# The primary species carry one atom of their element each, so these rows count
# aluminum, magnesium and sulfur however the species that hold them are written.
for (component, element) in ("AlO2-" => "Al", "Mg+2" => "Mg", "SO4-2" => "S")
    println(element, ":")
    for (name, amount) in component_in_solids(component)
        @printf("  %-18s %9.5f mol %s\n", name, amount, element)
    end
    println()
end
```

## 5. The oxidation state

```@example cem5
r = half_reaction(eq, "SO4-2", "HS-")
println("half-reaction : ", r.equation)
@printf("log K at 25 C : %.2f\n", r.logK⁰(T = 298.15))
@printf("pe            : %+.2f\n", pe(eq, model))
@printf("Eh            : %+.3f V\n", Eh(eq, model))
@printf("aqueous S(VI)  : %.3e mol\naqueous S(-II) : %.3e mol\n",
        ustrip(us"mol", moles(eq, "SO4-2")), ustrip(us"mol", moles(eq, "HS-")))
```

The caution of [the CEM III page](@ref ex-cem3-slag) applies unchanged, and for
the same reason: a couple buffers a potential only while both members are in
solution, and a cement paste puts nearly all of its sulfur into solids. The
number is a bound on the reducing power the slag brings, not a measurement of the
paste's redox state — and even a buffered value would be the *mutual* equilibrium
of every couple, which sulfate reduction is far too slow to reach during
hydration.

## [6. What the reacted fractions are worth](@id cem5-dor)

The three fractions of section 1 are the page's only real assumption, so the
page is obliged to say what turns on them. The round robin gives the answer its
own way: it reports the same two glasses at 7, 28 and 90 days, so sweeping the
*age* rather than an abstract parameter keeps every point on measured ground.
The pairs below are read off its Table 4, SEM image analysis, rounded to the two
laboratories' agreement.

```@example cem5
AGES = [(" 7 days", 0.35, 0.10), ("28 days", ALPHA_SLAG, ALPHA_ASH),
        ("90 days", 0.52, 0.25)]

@printf("%-9s %7s %7s   %-9s %9s %7s %9s\n",
        "age", "slag", "ash", "certified", "balance", "pH", "portlandite")
i_portlandite = findfirst(sp -> symbol(sp) == "Portlandite", cs.species)

# `let` rather than a bare loop: a top-level `for` that assigns to a name of the
# enclosing scope makes a NEW LOCAL, so `prev` would be read before it is ever
# written. Wrapping the sweep gives it a scope of its own.
let prev = nothing
    for (label, a_sl, a_as) in AGES
        # Ordered by age and continued, for the reason `solve_paste` gives: each
        # answer is the next one's start. The certificate decides every point.
        e, c = solve_paste(a_sl, a_as; start = prev)
        c.optimal && (prev = e)
        nn = ustrip.(us"mol", e.n)
        @printf("%-9s %6.0f %% %6.0f %%   %-9s %9.1e %7.3f %9.5f\n",
                label, 100a_sl, 100a_as, c.optimal, c.balance, pH(e, model),
                nn[i_portlandite])
    end
end
```

The package's own kinetics would answer this differently, and the disagreement is
worth knowing about: the Waller sigmoid shipped for a slag and a fly ash gives
0.29 and 0.32 at 28 days where the round robin measures 0.38–0.49 and 0.20. Both
are fits to particular materials, and "a slag" is not a substance — [the rate law
chapter](@ref sec-theory-kinetics) sets the two side by side.

Read the columns against each other. The **portlandite** falls by nearly 40 %
across the three ages — that is the pozzolanic reaction, and the reason the
family exists. The **pH** moves the other way and by almost nothing, six
hundredths of a unit, because it is the alkalis that set it and the glasses
release them only slowly; the calcium hydroxide is a floor beneath, not a lever.
And every one of the three **certifies**, with the element balances the table
prints, all below 10⁻¹⁰.

That is the useful conclusion, and it is worth stating as a limit on what the
page claims: between 7 and 90 days the slag's reacted fraction changes by half
again, and the assemblage's *character* does not. **The answer is sensitive to
the reacted fraction in its amounts and robust in its identity** — so a reader
who disagrees with the fractions assumed here can move them and keep the
qualitative reading, while a reader who wants the amounts must supply a measured
degree of reaction for their own materials.

!!! warning "Where this stops being true"
    Push the fractions to 1 — every grain of slag and every ash sphere fully
    dissolved — and the calculation stops having an answer at all: the
    minimization reports supersaturated hydrotalcite and layered double
    hydroxides it has no room to precipitate, an element balance off by 3·10⁻¹,
    and a pH of 14.4 that no cement paste has ever had. That is not a solver
    failure and not a gap in CEMDATA18. It is the formulation being asked an
    unphysical question: a 48 %-replaced binder whose glasses have entirely
    dissolved would have to place alkalis and aluminum that a real paste never
    releases, and no assemblage the database can form will hold them. The
    remedy is not a better minimizer. It is the reacted fraction.

## 7. Where a CEM V sits among the others

```@example cem5
function final_heat(file)
    t, Q = 0.0, 0.0
    for line in eachline(joinpath(datapath("experimental"), file))
        (startswith(line, "#") || startswith(line, "time")) && continue
        f = split(line, ',')
        t, Q = parse(Float64, f[1]), parse(Float64, f[3])
    end
    return t, Q
end

records = [
    ("smilauer2025-122-cemI-52.5R-cizkovice.csv", "CEM I 52.5 R", 0.0),
    ("smilauer2025-165-cemII-A-LL-42.5R-hranice.csv", "CEM II/A-LL", 0.13),
    ("smilauer2025-149-cemII-B-S-32.5R-mokra.csv", "CEM II/B-S", 0.28),
    ("smilauer2025-184-cemIII-A-42.5N-hranice.csv", "CEM III/A", 0.50),
    ("smilauer2025-200-cemV-A-S-V-32.5R-prachovice.csv", "CEM V/A (S-V)", 0.48),
    ("smilauer2025-121-cemIII-B-32.5N-mokra.csv", "CEM III/B", 0.73),
]
labels = String[]
heats = Float64[]
repl = Float64[]
for (file, label, r) in records
    _, Q = final_heat(file)
    push!(labels, label); push!(heats, Q); push!(repl, 100r)
    @printf("  %-16s replacement %4.0f %% (assumed)   Q = %6.1f J/g\n", label, 100r, Q)
end
```

```@example cem5
fig = scatter(repl, heats; legend = false, markersize = 7, color = :seagreen,
              xlabel = "clinker replaced (%, assumed at the EN 197-1 midpoint)",
              ylabel = "measured heat at ~28 days (J/g of binder)",
              title = "Every joule comes from the clinker",
              size = (780, 420),
              bottom_margin = 10Plots.mm, left_margin = 10Plots.mm)
for (x, y, l) in zip(repl, heats, labels)
    annotate!(fig, x, y + 8, text(l, 8, :left, :bottom))
end
plot!(fig; xlims = (-8, 85), ylims = (200, 410))
savefig(fig, "cem5-heat.svg"); nothing # hide
```

![](cem5-heat.svg)

!!! warning "The abscissa is assumed, the ordinate is measured"
    The heats are measured on one instrument at 20 °C. The replacement levels are
    **not in the deposit** — each is the midpoint of the EN 197-1 range for its
    designation, the same assumption every page here makes. So the figure shows a
    real trend read against an estimated axis, and the scatter about it is as much
    the spread of the ranges as it is chemistry. What it does establish is the
    ordering and the magnitude: half the clinker removed costs roughly a third of
    the heat.

The CEM V/A and the CEM III/A sit almost on top of each other, at comparable
replacement and comparable heat, which is what makes the family interesting: the
same reduction in heat and clinker, reached with two *different* constituents, and
therefore with a different assemblage, a different aluminum balance and a
different long-term reactivity. The equilibrium calculation above distinguishes
them where the calorimeter does not.

## Where to go next

The pages of this group compute where a binder ends; how it gets there is the
subject of the applications in time, beginning with
[Cement clinker hydration kinetics](@ref), and of
[The full Portland cement, through its pore solution](@ref ex-ionic-opc) for the
coupled route.
