# [CEM IV/A (V) and CEM IV/B (V): a pozzolanic binder, and the C-S-H that has to carry the aluminum](@id ex-cem4-pozzolanic)

!!! info "Before this page"
    [CEM III/A, a blastfurnace cement](@ref ex-cem3-slag), whose skeleton this page
    follows, and [Solid solutions](@ref sec-theory-solid-solutions).

A CEM IV replaces 11 % to 55 % of the clinker with a pozzolana — siliceous fly
ash, natural or calcined pozzolana, silica fume, or a mixture of them. The
replacement is not a filler: a pozzolana **consumes portlandite** and makes more
C-S-H out of it, which is why a pozzolanic binder is specified where the paste
must be denser and less alkaline.

Two things follow for the calculation, and this page is about both.

1. **The aluminum has somewhere to go, and the model must let it.** Fly ash
   brings almost as much Al as Si. In a real paste most of that aluminum ends up
   *inside* the calcium silicate hydrate, as C-A-S-H. `CSHQ`, the C-S-H model
   the CEM I pages use, has **no aluminum end-member at all**, so a calculation
   that keeps it forces every atom of aluminum into the other aluminum-bearing
   phases: the hydrogarnets, AFt and AFm.
   `CNASH_ss` is the model that can take it.
2. **Portlandite becomes the limiting reagent.** Past a certain replacement it
   runs out, and what the paste can do afterwards changes.

## 0. What a pozzolana actually does

A pozzolana is not a binder on its own. Mix siliceous fly ash with water and
nothing happens. What makes it work is that a Portland clinker, while hydrating,
produces a large amount of **portlandite** — calcium hydroxide, ``\mathrm{Ca(OH)_2}``
— which is a by-product of the silicates:

```math
2\,\mathrm{C_3S} + 6\,\mathrm{H} \;\longrightarrow\; \mathrm{C_3S_2H_3} + 3\,\mathrm{CH}
```

Roughly a fifth of a hydrated CEM I paste is portlandite, and it contributes
little strength while being the phase most easily leached. The **pozzolanic
reaction** puts it to work: the ash's amorphous silica consumes it and makes more
of the phase that does carry the strength, the calcium silicate hydrate.

```math
\mathrm{S} + 1.7\,\mathrm{CH} + \text{water} \;\longrightarrow\; \mathrm{C_{1.7}SH_x}
```

So a pozzolanic binder trades a weak, soluble phase for a strong one, at the cost
of a slower reaction. Section 6 measures exactly that trade, by sweeping the
replacement level and watching the portlandite go.

But the ash brings something else that the equation above ignores, and it is what
the rest of the page is about: **aluminum**, almost as much as silicon. Where it
goes decides the answer, and it depends entirely on which C-S-H model the
calculation is given.

```@example cem4
using ChemistryLab
using DynamicQuantities
using Logging
using OptimaSolver
using OrderedCollections
using Printf
using Plots
default(framestyle = :box, grid = false)

substances = build_species(datapath("cemdata18-thermofun.json"); verbose = false)
byname = Dict(symbol(s) => s for s in substances)
molar_mass(n) = ustrip(us"g/mol", byname[n][:M])
nothing # hide
```

## 1. What is assumed here, and why all of it is

!!! warning "This page has no measured specimen behind it"
    The calorimetry deposit this package ships [Smilauer2025data](@cite) carries
    CEM I, CEM II, CEM III and CEM V records, and **no CEM IV**. So unlike
    [the CEM III page](@ref ex-cem3-slag) and [the CEM II page](@ref ex-cem2-blended),
    nothing below is anchored to a specimen. It is a **model study of a family**,
    run on a composition chosen inside the EN 197-1 range, and the numbers
    describe what a pozzolanic binder of that composition does — not what any
    particular cement contains.

    That distinction is worth keeping, because the aluminum question this page is
    about does not depend on the exact analysis, while the assemblage does.

```@example cem4
# ASSUMED: a Bogue composition representative of a CEM I clinker, here
# the Bogue composition of the CEM I 52.5 N of [Lavergne2018](@citet), Table 9.
bogue = literature_table("Lavergne2018", "cement_bogue")
CLINKER = OrderedDict(zip(bogue.phase, bogue.percent ./ 100))

# ASSUMED: a siliceous (class V) fly ash analysis of the kind European standards
# admit — low calcium, high silica and alumina, and the alkalis that make the
# aluminum question interesting.
FLYASH = OrderedDict("SiO2" => 0.53, "Al2O3" => 0.26, "Fe2O3" => 0.07,
                     "CaO" => 0.04, "MgO" => 0.02, "K2O" => 0.025,
                     "Na2O" => 0.008, "SO3" => 0.005)

# THE CLINKER'S OWN ALKALIS, which Bogue does not account for. They are minor
# oxides outside the four-phase decomposition, and in a paste they are what fixes
# the pH: they dissolve almost completely and stay in solution while the calcium
# is held down by portlandite at 12.5. On this page they matter twice over, since
# the whole question of section 6 is what the pozzolana does to the alkalis.
# ASSUMED at a usual industrial level, as a fraction of the CLINKER mass.
# [The CEM III page](@ref cem3-alkali) sweeps this over the industrial range and
# shows that it, and almost nothing else, is what moves the pH.
ALKALIS = OrderedDict("K2O" => 0.008, "Na2O" => 0.002)

# ASSUMED: the midpoint of the EN 197-1 range for a CEM IV/A, which is 65-89 %
# clinker and 11-35 % pozzolana. Section 7 goes to a CEM IV/B, and asks whether
# its full-reaction limit needs phases this list does not declare.
ASH_FRACTION = 0.23
ASH_FRACTION_B = 0.45
# The gypsum of the same cement, 4.6 % of the binder.
GYPSUM = literature_value("Lavergne2018", "gypsum_percent") / 100
WB = 0.50
BINDER_G = 100.0

# HOW MUCH REACTS. Two ceilings, and the reacted fraction is the lower.
#
# The WATER ceiling is Powers (1948): about 0.42 g of water per gram of cement
# is needed for complete hydration, 0.23 g written into the hydrates and 0.19 g
# held in gel pores too fine to feed a grain. It is a property of the pore space
# rather than of the grain, so it caps every constituent alike. At w/b = 0.50
# it does not bind at all -- water is abundant here, and the function correctly
# returns 1. It binds on [the CEM V page](@ref ex-cem5-composite), mixed at 0.40.
ALPHA_WATER = powers_alpha_max(WB)              # 1.0 at w/b = 0.50

# The KINETIC ceiling is what actually limits a fly ash, and by a long way. Only
# the glassy fraction reacts at all -- the crystalline mullite and quartz do not
# dissolve on any relevant time scale -- and the glass itself is slow. The RILEM
# TC 238-SCM round robin [Durdzinski2017](@cite) measured a siliceous fly ash at
# 30 % replacement and w/b 0.40 in seven laboratories; its Table 5 at 28 days
# gives 20 % by SEM image analysis, the technique the study found most
# consistent, with XRD-PONKCS between 19 % and 23 %. The study's verdict on the
# precision of any of these: "at best +/- 5 %".
#
# ASSUMED from that table. Section 7 drives it to 1 and shows what changes.
# Degrees of reaction measured by SEM image analysis on sealed pastes, in
# percent: [Durdzinski2017](@citet), Table 5, from data/literature/Durdzinski2017.json.
sem(material, age) = literature_table("Durdzinski2017", "degree_of_reaction";
    technique = "SEM-IA", material, curing = "sealed", age_days = age).degree_percent
mean_percent(x) = sum(x) / length(x)
ALPHA_ASH = min(only(sem("SFA", 28)) / 100, ALPHA_WATER)
ALPHA_CLINKER = ALPHA_WATER

@printf("water ceiling at w/b = %.2f : %.3f\n", WB, ALPHA_WATER)
@printf("reacted: clinker %.0f %%, fly ash %.0f %%\n",
        100ALPHA_CLINKER, 100ALPHA_ASH)
```

An equilibrium calculation cannot tell glass from mullite, and it cannot tell a
dissolved grain from an intact one: it reacts whatever budget it is handed. So
the reacted fraction is not a refinement to be added later — it is part of
posing the problem, and a page that leaves it at 1 has quietly asked what the
paste becomes after every ash sphere has dissolved, which is a question about
geological time, not about a specimen at 28 days.

## 2. Two species lists, differing in one phase

Everything is held fixed between the two systems except the C-S-H model. That is
what makes the comparison mean something:

```@example cem4
pure = split(
    "C3S C2S C3A C4AF Gp Anh Cal Portlandite ettringite monosulphate12 " *
    "monocarbonate hemicarbonate C4AH13 C3AH6 C3FH6 straetlingite " *
    "hydrotalcite Brc FeOOHmic AlOHmic Amor-Sl Mgs " *
    # Aluminum sinks that CEMDATA18 documents and an earlier version of this
    # list simply did not declare. The siliceous hydrogarnet `C3AS0.84H4.32`
    # is the one that matters most here: it is the ALUMINUM end-member of the
    # family whose iron end-member was already present, and a blended binder
    # puts a great deal of aluminum into it. Leaving it out does not make the
    # calculation conservative -- it makes it insoluble, because the element
    # has to go somewhere.
    "C3AS0.84H4.32 C3AS0.41H5.18 straetlingite7 Gbs AlOHam " *
    "M4A-OH-LDH M6A-OH-LDH M8A-OH-LDH C2AH7.5 C4AH11 C4AH19 " *
    "K2SO4 syngenite Na2SO4"
)
aqueous = ["SO4-2", "CO2@", "O2@"]

CSHQ = ["CSHQ-JenD", "CSHQ-JenH", "CSHQ-TobD", "CSHQ-TobH", "KSiOH", "NaSiOH"]
CNASH = ["T2C-CNASHss", "T5C-CNASHss", "TobH-CNASHss",
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
FEAL = ["C3AFS0.84H4.32", "C3FS0.84H4.32"]

function system(gel_name, gel_members)
    sp = speciation(substances, vcat(pure, gel_members, FEAL, aqueous);
                    aggregate_state = [AS_AQUEOUS])
    # CNASH_ss mixes on the sites of Myers et al. (2014), as it ships in
    # data/solid_solutions.toml; CSHQ and the hydrogarnet mix their end-members.
    members = [byname[m] for m in gel_members]
    mixing = gel_name == "CNASH_ss" ? sublattice_model("Myers2014:cnash", members) : IdealSolidSolutionModel()
    ss = [SolidSolutionPhase(gel_name, members; model = mixing),
          SolidSolutionPhase("C3(AF)S0.84H", [byname[m] for m in FEAL])]
    return ChemicalSystem(sp, CEMDATA_PRIMARIES; solid_solutions = ss)
end

cs_q = system("CSHQ", CSHQ)
cs_n = system("CNASH_ss", CNASH)
# The activity model Cemdata18 prescribes (its Eq. C.1): extended Debye-Hückel,
# with the common ion size and B-dot the paper gives for KOH solutions (it also
# gives them for NaOH). The alkalis of the clinker, and those of the fly ash in a similar ratio, are
# mostly potassium.
model = cemdata18_activity_model(:KOH)
# The molar K/Na ratio of its alkalis, which is why the KOH set applies.
Mox(ox) = ustrip(us"g/mol", Species(ox)[:M])
println("molar K/Na of the alkalis: ", round((2 * ALKALIS["K2O"] / Mox("K2O")) / (2 * ALKALIS["Na2O"] / Mox("Na2O")); digits = 1))

@printf("CSHQ system     : %d species\n", length(cs_q.species))
@printf("CNASH_ss system : %d species\n", length(cs_n.species))
```

!!! danger "Never both at once"
    `CSHQ`, `CNASH_ss` and the `ECSH` family are three *models of one gel*, not
    three phases. Declaring two of them counts the same calcium silicate hydrate
    twice, and `ChemicalSystem` refuses the pair: by composition when two
    end-members are one substance (`KSiOH` of `CSHQ` is `ECSH1-KSH`), and by the
    models listed in `data/gel_models.toml` otherwise — `CSHQ` and `CNASH_ss`
    share no composition, and until 0.25.1 that pair was accepted without a word.
    The two systems above are alternatives, built separately and compared, which
    is the only correct way to use them.

## 3. The budget

```@example cem4
# `α_ash` is a keyword rather than a constant so that section 7 can drive it,
# and the state carries only the reacted clinker: what has not reacted is still
# in the specimen but is not at equilibrium with the pore solution.
function budget(cs; ash, wb = WB, α_ash = ALPHA_ASH, α_clinker = ALPHA_CLINKER)
    clinker_frac = 1 - ash - GYPSUM
    state = ChemicalState(cs)
    for (phase, frac) in CLINKER
        set_quantity!(state, phase,
            α_clinker * BINDER_G * clinker_frac * frac / molar_mass(phase) * u"mol")
    end
    # The calcium sulfate is soluble and carries no ceiling, and all of the
    # mixing water enters: the ceiling limits how far the reaction goes, not how
    # much water was poured in.
    set_quantity!(state, "Gp", BINDER_G * GYPSUM / molar_mass("Gp") * u"mol")
    set_quantity!(state, "H2O@", BINDER_G * wb / molar_mass("H2O@") * u"mol")
    b = Float64.(cs.SM.A) * ustrip.(us"mol", state.n)
    # The alkalis leave the grain as it dissolves: same fraction as the clinker.
    b .+= oxide_budget(ALKALIS, cs.SM.primaries;
                       mass = BINDER_G * clinker_frac * α_clinker * u"g")
    ash > 0 && (b .+= oxide_budget(FLYASH, cs.SM.primaries;
                                   mass = BINDER_G * ash * α_ash * u"g"))
    return state, b
end

st_q, b_q = budget(cs_q; ash = ASH_FRACTION)
st_n, b_n = budget(cs_n; ash = ASH_FRACTION)

comps = String.(symbol.(cs_n.SM.primaries))
for (c, v) in zip(comps, b_n)
    abs(v) > 1.0e-6 && @printf("  %-8s %10.5f mol\n", c, v)
end
```

The silicon-to-aluminum ratio of that budget is the whole issue. A CEM I paste
has roughly seven times more Si than Al; this one has under three.

## 4. The same paste, two C-S-H models

```@example cem4
eq_q, c_q = equilibrate_certified(st_q; model = model, b = b_q)
eq_n, c_n = equilibrate_certified(st_n; model = model, b = b_n)

for (label, cs, eq, c) in (("CSHQ", cs_q, eq_q, c_q),
                           ("CNASH_ss", cs_n, eq_n, c_n))
    @printf("%-10s optimal=%-5s  worst SI=%+.2e  balance=%.1e  pH=%.3f  V=%.2f cm3\n",
            label, c.optimal, c.worst_supersaturation, c.balance,
            pH(eq, model), ustrip(uconvert(us"cm^3", volume(eq).total)))
end
```

Where does the aluminum end up? The conservation matrix already answers it. The
row of the primary species `AlO2-` counts one unit per aluminum atom, so
multiplying that row by the amounts distributes the element over the phases that
hold it:

```@example cem4
"""Moles of one conservation component held by each solid phase, largest first."""
function component_in_solids(cs, eq, component; tol = 1.0e-4)
    n = ustrip.(us"mol", eq.n)
    row = findfirst(==(component), String.(symbol.(cs.SM.primaries)))
    row === nothing && error("$component is not a component of this system")
    A = Float64.(cs.SM.A)
    out = [(symbol(cs.species[i]), A[row, i] * n[i]) for i in cs.idx_crystal]
    return sort(filter(p -> last(p) > tol, out); by = last, rev = true)
end

for (label, cs, eq) in (("CSHQ", cs_q, eq_q), ("CNASH_ss", cs_n, eq_n))
    println(label, " — where the aluminum is:")
    for (name, al) in component_in_solids(cs, eq, "AlO2-")
        @printf("  %-18s %9.5f mol Al\n", name, al)
    end
    println()
end
```

## 5. The assemblages, side by side

```@example cem4
function assemblage(cs, eq; tol = 1.0e-4)
    n = ustrip.(us"mol", eq.n)
    sort([(symbol(cs.species[i]), n[i]) for i in cs.idx_crystal if n[i] > tol];
         by = last, rev = true)
end

for (label, cs, eq) in (("CSHQ", cs_q, eq_q), ("CNASH_ss", cs_n, eq_n))
    println(label, ":")
    for (name, amount) in assemblage(cs, eq)
        @printf("  %-18s %9.5f mol\n", name, amount)
    end
    # The gel as one phase: its atomic ratios, from the amounts of its members.
    el = solid_solution_totals(eq, label).elements
    @printf("  gel Ca/Si = %.2f, Al/Si = %.2f\n\n", el[:Ca] / el[:Si], get(el, :Al, 0.0) / el[:Si])
end
```

!!! note "Read the C-A-S-H as a total, and its ratios"
    The eight `CNASH_ss` end-members span a space of rank 5 in their elements, so
    the element balance fixes only five combinations of their amounts. The mixing
    on the sites of [Myers2014](@citet) fixes the rest: its energy is
    strictly convex (the occupancy has rank 8), so the eight amounts are
    determined. They are still a description of one gel, and the total, the
    Ca/Si and the Al/Si are what compare with a measurement.

    Here the two gels differ most in their Ca/Si, 1.60 against 1.17: `CNASH_ss`
    takes less calcium, and its paste keeps more portlandite, 0.380 mol against
    0.237.

## 6. Portlandite is the limiting reagent

The pozzolanic reaction consumes calcium hydroxide, so sweeping the replacement
level says how much of it survives. The sweep is run **twice** — at the reacted
fraction of section 1, and in the limit where the ash has entirely dissolved —
because the two answer different questions and are routinely confused:

```@example cem4
fractions = 0.0:0.10:0.30
curves = Dict{Float64, Tuple{Vector{Float64}, Vector{Float64}, Vector{Bool}}}()
i_ch = findfirst(s -> symbol(s) == "Portlandite", cs_n.species)
for α in (ALPHA_ASH, 1.0)
    ch, phs, ok = Float64[], Float64[], Bool[]
    prev = nothing
    for f in fractions
        st, b = budget(cs_n; ash = f, α_ash = α)
        # CONTINUATION along the sweep: each point starts from its neighbor's
        # answer rather than from a fresh paste.
        #
        # That is safe here for a reason that is CHECKED rather than assumed.
        # CNASH_ss mixes on its sites and the hydrogarnet ideally, both convex
        # (`mixing_convexity`), and `SolidSolutionPhase` refuses a model whose
        # mixing energy is concave -- so the Gibbs function is convex, its minimum is unique,
        # and a continuation cannot change WHAT is found, only whether the search
        # finds it, which on a 109-species cement is the whole difficulty. Waive
        # that refusal with `check_convexity = false` and none of it holds: inside
        # a spinodal the minimum is two compositions, the certificate loses the
        # sufficiency that rests on convexity, and the start would then decide
        # which branch you land on. The certificate still decides every point
        # here, and a start is reused only once it has been certified.
        #
        # The certificate is printed below; the warning of a refusal would
        # only repeat it.
        eq, c = with_logger(NullLogger()) do
            equilibrate_certified(something(prev, st); model = model, b = b)
        end
        c.optimal && (prev = eq)
        n = ustrip.(us"mol", eq.n)
        push!(ch, n[i_ch])
        push!(phs, pH(eq, model))
        push!(ok, c.optimal)
        # The certificate is carried through to the figure, not just printed. A
        # point that does not certify is a point whose portlandite and pH mean
        # nothing, and a sweep that hides one draws a curve through a number the
        # solver never stood behind.
        @printf("  ash reacted %3.0f %% of %3.0f %%   certified %-5s   balance %8.1e   portlandite %8.5f mol   pH %.3f\n",
                100α, 100f, c.optimal, c.balance, n[i_ch], phs[end])
    end
    curves[α] = (ch, phs, ok)
end
```

```@example cem4
x = 100 .* collect(fractions)
p1 = plot(; xlabel = "fly ash (% of binder)", ylabel = "portlandite (mol / 100 g)",
          title = "Calcium hydroxide consumed", legend = :bottomleft)
p2 = plot(; xlabel = "fly ash (% of binder)", ylabel = "pH",
          title = "Pore solution pH", legend = false)
for (α, color, lab) in ((ALPHA_ASH, :seagreen, "ash reacted 20 % (28 days)"),
                        (1.0, :indianred, "ash reacted 100 % (the limit)"))
    ch, phs, ok = curves[α]
    plot!(p1, x, ch; marker = :circle, color, label = lab)
    plot!(p2, x, phs; marker = :circle, color, label = lab)
    # A point the certificate refused is drawn hollow and black, so that a
    # reader sees the gap in the evidence rather than a smooth curve through it.
    bad = .!ok
    if any(bad)
        scatter!(p1, x[bad], ch[bad]; marker = :circle, markersize = 8,
                 markercolor = :white, markerstrokecolor = :black,
                 label = "not certified")
        scatter!(p2, x[bad], phs[bad]; marker = :circle, markersize = 8,
                 markercolor = :white, markerstrokecolor = :black, label = "")
    end
end
fig = plot(p1, p2; layout = (1, 2), size = (900, 380),
           bottom_margin = 10Plots.mm, left_margin = 10Plots.mm)
savefig(fig, "cem4-sweep.svg"); nothing # hide
```

![](cem4-sweep.svg)

Every point of both branches certifies. A point the certificate refused would be
drawn hollow and black, and would not be a result: its portlandite and its pH
would be whatever the iteration stopped at. It would be drawn rather than
dropped, because dropping it would put a smooth curve where the evidence has a
hole.

The **portlandite falls** because the ash's silica turns it into more C-S-H —
that is the pozzolanic reaction, and it is the property the family is specified
for. **How far it falls is entirely a matter of how much ash has reacted**: in
the limit a CEM IV/A exhausts its portlandite inside the EN 197-1 range, at 30 %
ash, while at 28 days the same binder still has about half of it. Reporting the
first as though it described a specimen is the single easiest mistake to make
with an equilibrium code, and it is not a small one — portlandite is what buffers
the pH and what protects the reinforcement.

The **pH moves much less** than the portlandite, because in a cement paste it is
the alkalis that set it, not the calcium hydroxide; portlandite only fixes a
floor around 12.5 at 25 °C. At 28 days it falls by 0.17 unit across the sweep. In
the limit it *rises*, by 0.14, because the dissolved ash brings its own alkalis,
mostly potassium, and `CNASH_ss`, the gel of this sweep, has no potassium member
to take them. Each model expresses a different part of what a pozzolana does to
the alkalis: `CSHQ` holds potassium and sodium through its `KSiOH` and `NaSiOH`
end-members but no aluminum, `CNASH_ss` holds sodium and aluminum but no
potassium. Neither holds all three, which is one more reason to compute both
rather than trust one.

## 7. A CEM IV/B, and its full-reaction limit

Everything so far was a CEM IV/**A**, 23 % ash. Take it to a CEM IV/**B** — the
midpoint of 36–55 % — and ask the same question twice: once at the reacted
fraction a specimen has at 28 days, and once in the limit where every ash sphere
has dissolved.

!!! tip "How to read `optimal = false`, because it is not one thing"
    The certificate reports **three** residuals, and which one is too large says
    what went wrong. A novice reading only the `optimal` flag will mistake a
    solver difficulty for a statement about chemistry — this page did, for two
    releases.

    | the residual that is large | what it means | what to do |
    |:--|:--|:--|
    | **worst supersaturation**, clearly positive | a phase the system could form is not in your species list | declare it |
    | **element balance** | the answer does not conserve matter — it is not an answer at all | start somewhere else |
    | **stationarity** alone, the other two small | a genuinely hard point | nothing simple; report it as unresolved |

    A supersaturation of `+3e+02` names a missing phase. One of `+5e-02` beside a
    broken balance names a failed solve. Read them, not the flag.

```@example cem4
eq_b, c_b = nothing, nothing          # the full-reaction case, kept below

# The 28-day fraction, each model started from its own fresh paste.
for (label, cs) in ("CSHQ" => cs_q, "CNASH_ss" => cs_n)
    st, b = budget(cs; ash = ASH_FRACTION_B, α_ash = ALPHA_ASH)
    eq, c = equilibrate_certified(st; model = model, b = b)
    @printf("%2.0f %% ash reacted %3.0f %%  %-10s optimal=%-5s balance=%.1e  pH=%.3f\n",
            100ASH_FRACTION_B, 100ALPHA_ASH, label, c.optimal, c.balance, pH(eq, model))
end
```

Both models certify, and they disagree — which is the whole point of the section
and is now a comparison of **two proved answers** rather than of a success
against a failure.

```@example cem4
# The full-reaction limit, every ash sphere dissolved. The certificate carries
# the ionic strength it was reached at, and whether that lies within the range
# the activity model is stated for.
for (label, cs) in ("CSHQ" => cs_q, "CNASH_ss" => cs_n)
    st, b = budget(cs; ash = ASH_FRACTION_B, α_ash = 1.0)
    eq, c = equilibrate_certified(st; model = model, b = b)
    (label == "CNASH_ss") && (global eq_b, c_b = eq, c)
    @printf("%2.0f %% ash reacted 100 %%  %-10s optimal=%-5s balance=%.1e  pH=%.3f  I=%.2f mol/kg  within range: %s\n",
            100ASH_FRACTION_B, label, c.optimal, c.balance, pH(eq, model),
            c.ionic_strength, c.within_activity_range)
end
@printf("range Cemdata18 states for its model: I up to about %.1f mol/kg\n",
        activity_model_range(model))
```

At 28 days a CEM IV/B is an ordinary calculation for **both** models: what the
ash has released by then fits in the phases the paste can form, and both answers
certify. They differ by about 0.2 unit of pH, and that difference — not a
refusal — is what says the aluminum matters.

**In the limit both certify as well**, a whole unit of pH apart. With `CSHQ` the
pH is 12.10, below the floor of about 12.5 that portlandite holds, so the
portlandite is gone and the alkalis have the `KSiOH` and `NaSiOH` members to go
to; with `CNASH_ss` it is 13.20, since that gel has no potassium member.

Until this version the page ran the limiting law, `å = 0`, the configuration of
one GEM-Selektor run that the tests reproduce, and in the limit neither model
certified: the solve stopped with element balances of 0.26 and 0.35 mol. The
activity model Cemdata18 prescribes certifies both from the first solve, and for
`CSHQ` nothing else changed. What did not close at the limit was the limiting
law, not the minimization and not the phase list.

Read the `CNASH_ss` limit for what it is. Its ionic strength, 1.09 mol/kg, is
past the range Cemdata18 states for its model, about 1 mol/kg, and the
certificate says so (`within_activity_range` is `false`): it is a composition
consistent with an extrapolated activity model, certified as such and no more.
How much the extrapolation weighs is measured by changing the activity model
alone. Give every ion its own size, which is what `HKFActivityModel()` does with
its default parameters, and leave the budget, the phase list and the solver as
they are:

```@example cem4
perion = HKFActivityModel()
for (label, cs) in ("CSHQ" => cs_q, "CNASH_ss" => cs_n)
    st, b = budget(cs; ash = ASH_FRACTION_B, α_ash = 1.0)
    eq, c = equilibrate_certified(st; model = perion, b = b)
    @printf("100 %% ash reacted  %-10s optimal=%-5s balance=%.1e  pH=%.3f  I=%.2f mol/kg\n",
            label, c.optimal, c.balance, pH(eq, perion), ionic_strength(eq))
end
@printf("range the manual states for this model: I up to %.1f mol/kg\n",
        activity_model_range(perion))
```

For either gel the pH moves by 0.02 unit (12.10 to 12.12, 13.20 to 13.18): at
these ionic strengths the answer hardly depends on the ion sizes, and what sets
the two limits apart is the gel, not the activity model.

Whether the limit needs phases this species list does not contain is a claim
about chemistry, and the way to test it is to add the phases and let the
certificate decide. The limit is not an idle question: it is where a pozzolanic
binder is heading over years, and the regime an alkali-activated system is in
from the start. An alkaline aluminosilicate can precipitate **zeolites**, and
CEMDATA18 carries five, none of the families this binder could form. That is
what [the zeolite extension](@ref sec-zeolites) was built for.

```@example cem4
zeo_db = build_species(datapath("cemdata18-zeolites.json"); verbose = false)
zeo_byname = Dict(symbol(s) => s for s in zeo_db)
added = sort(collect(setdiff(Set(keys(zeo_byname)), Set(keys(byname)))))

# Only the ones this paste could form. Two of the twenty-eight are a chloride
# and a nitrate sodalite, and this binder carries neither element: declaring them
# would widen the species list to every aqueous chloride and nitrate species in
# the database, all of them on a budget of exactly zero, for no phase that can
# form. Dropping them is not a modeling choice, it is arithmetic.
carries(sp, el) = haskey(atoms(sp), Symbol(el))
zeolites = [z for z in added
            if !carries(zeo_byname[z], "Cl") && !carries(zeo_byname[z], "N")]

@printf("%d phases added by the extension, %d of them usable here\n",
        length(added), length(zeolites))
println("left out (no Cl and no N in this binder): ",
        join(setdiff(added, zeolites), ", "))
for z in zeolites
    print(z, "  ")
end
```

```@example cem4
pure_z = vcat(String.(pure), zeolites)
sp_z = speciation(zeo_db, vcat(pure_z, CNASH, FEAL, aqueous);
                  aggregate_state = [AS_AQUEOUS])
ss_z = [SolidSolutionPhase("CNASH_ss", [zeo_byname[m] for m in CNASH];
                           model = sublattice_model("Myers2014:cnash", [zeo_byname[m] for m in CNASH])),
        SolidSolutionPhase("C3(AF)S0.84H", [zeo_byname[m] for m in FEAL])]
cs_z = ChemicalSystem(sp_z, CEMDATA_PRIMARIES; solid_solutions = ss_z)

# The full-reaction limit, through the same `budget` as `c_b`, so that the two
# answers are to one question and differ only in the phases declared.
st_z, b_z = budget(cs_z; ash = ASH_FRACTION_B, α_ash = 1.0)

eq_z, c_z = equilibrate_certified(st_z; model = model, b = b_z)
@printf("with zeolites: optimal=%-5s  worst SI=%+.2e  pH=%.3f\n",
        c_z.optimal, c_z.worst_supersaturation, pH(eq_z, model))
@printf("without      : optimal=%-5s  worst SI=%+.2e  pH=%.3f\n",
        c_b.optimal, c_b.worst_supersaturation, pH(eq_b, model))

nz = ustrip.(us"mol", eq_z.n)
formed = sort([(symbol(cs_z.species[i]), nz[i])
               for i in cs_z.idx_crystal
               if nz[i] > 1.0e-6 && symbol(cs_z.species[i]) in zeolites];
              by = last, rev = true)
if isempty(formed)
    println("\nno zeolite is stable in this paste")
else
    println("\nzeolites formed:")
    for (name, amount) in formed
        @printf("  %-16s %9.5f mol\n", name, amount)
    end
end
```

With the zeolites declared the limit certifies again, at the same pH, and every
zeolite is undersaturated: **no zeolite is stable in this paste**. That is a
result because the certificate tests every declared phase, so it proves the
absence. Without the extension the calculation could say nothing about a
zeolite: a phase absent from the species list is not reported as undersaturated,
it is not reported at all.

An earlier version of this page reported FAU-Y-K here. It answered a different
question: its budget for this block left out the alkalis of the clinker, and its
gel and its activity model were not those Cemdata18 publishes. The block now
takes its budget from the same `budget` as `c_b`.

!!! danger "An equilibrium at 45 % replacement is not a 28-day paste"
    Everything above is the state the paste *tends to*, with the whole ash taken
    as reactive. A real fly ash binder reaches a fraction of it in a month. The
    figure is a map of the family's limit, and the way to a date on the calendar
    is the coupled route of [the kinetic pages](@ref ex-ionic-opc), which needs a
    rate law for the ash — `waller` ships published parameters for fly ash, and
    that is the subject of the coupled runs rather than of this page.

## Where to go next

The composite binder, which carries a slag and a pozzolana at once, is
[CEM V/A (S-V): a composite binder, two glasses at once](@ref ex-cem5-composite).
