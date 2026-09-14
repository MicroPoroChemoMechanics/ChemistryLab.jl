# =============================================================================
#  GENERATED FILE — do not edit.
#
#  Extracted from docs/src/examples/cem5_composite.md, whose executed blocks it reproduces in
#  order. Edit the page, then regenerate:
#
#      julia --project=. scripts/_page_to_script.jl
#
#  `test/scripts.jl` fails if this file and its page disagree.
#
#  Usage:
#      julia --project=docs scripts/cem5_composite.jl
#
#  No `Pkg.activate` here, deliberately: the active project is global process
#  state and this file is meant to be `include`d as well as run.
#
#  Formatted with Runic like the rest of the repository. `test/scripts.jl`
#  compares SYNTAX TREES rather than text, so formatting is free to differ
#  from the page while a real divergence still fails the suite.
# =============================================================================

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
ZEOLITES = [
    z for z in sort(
            collect(
                setdiff(
                    Set(symbol.(substances)),
                    Set(symbol.(base))
                )
            )
        )
        if !carries(byname[z], "Cl") && !carries(byname[z], "N")
]
nothing # hide

rec = joinpath(
    datapath("experimental"),
    "smilauer2025-200-cemV-A-S-V-32.5R-prachovice.csv"
)
for line in eachline(rec)
    startswith(line, "#") || break
    occursin(r"cement:|blaine:|wb:|temperature:", line) && println(line)
end

# ASSUMED: a Bogue composition representative of a CEM I clinker.
CLINKER = OrderedDict("C3S" => 0.65, "C2S" => 0.11, "C3A" => 0.11, "C4AF" => 0.08)

# ASSUMED: midpoint of the EN 197-1 CEM V/A range, which is 40-64 % clinker
# with 18-30 % slag and 18-30 % pozzolana.
SLAG_FRACTION = 0.24
ASH_FRACTION = 0.24
GYPSUM = 0.046

# ASSUMED: the same two analyses used on the CEM III and CEM IV pages, so that
# the three calculations differ in their proportions and not in their inputs.
SLAG = OrderedDict(
    "CaO" => 0.41, "SiO2" => 0.36, "Al2O3" => 0.11,
    "MgO" => 0.08, "SO3" => 0.02
)
FLYASH = OrderedDict(
    "SiO2" => 0.53, "Al2O3" => 0.26, "Fe2O3" => 0.07,
    "CaO" => 0.04, "MgO" => 0.02, "K2O" => 0.025,
    "Na2O" => 0.008, "SO3" => 0.005
)

# THE CLINKER'S OWN ALKALIS, which Bogue does not account for -- minor oxides
# outside the four-phase decomposition, and in a paste what fixes the pH: they
# dissolve almost completely and stay in solution, while the calcium is held down
# by portlandite at 12.5. ASSUMED at a usual industrial level, as a fraction of
# the CLINKER mass. The fly ash brings alkalis of its own, and section 3's table
# shows which constituent brought what. [The CEM III page](@ref cem3-alkali)
# sweeps this input over the industrial range: it, and almost nothing else, is
# what moves the pH.
ALKALIS = OrderedDict("K2O" => 0.008, "Na2O" => 0.002)

WB = 0.4            # MEASURED, from the record above
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
ALPHA_ASH = min(0.2, ALPHA_WATER)

# The clinker is taken AT its water ceiling. At 28 days it has not quite got
# there, so the assemblage below is an upper bound on what the clinker
# contributes; section 6 moves it too.
ALPHA_CLINKER = ALPHA_WATER

@printf("water ceiling (%s, w/b = %.2f) : %.3f\n", CURING, WB, ALPHA_WATER)
@printf(
    "reacted: clinker %.0f %%, slag %.0f %%, fly ash %.0f %%\n",
    100ALPHA_CLINKER, 100ALPHA_SLAG, 100ALPHA_ASH
)

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
gel = [
    "T2C-CNASHss", "T5C-CNASHss", "TobH-CNASHss",
    "5CA", "5CNA", "INFCA", "INFCN", "INFCNA",
]
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

species = speciation(
    substances, vcat(pure, gel, feal, redox_species);
    aggregate_state = [AS_AQUEOUS]
)
ss = [
    SolidSolutionPhase("CNASH_ss", [byname[m] for m in gel]),
    SolidSolutionPhase("C3(AF)S0.84H", [byname[m] for m in feal]),
]
cs = ChemicalSystem(species, CEMDATA_PRIMARIES; solid_solutions = ss)
model = HKFActivityModel(å = 0.0, Ḃ = 0.097637, Kₙ = 0.0)

components = String.(symbol.(cs.SM.primaries))
@printf(
    "%d species, %d components: %s\n",
    length(cs.species), length(components), join(components, " ")
)
@printf(
    "charge (`Zz`) kept as a conservation component: %s\n",
    "Zz" in components ? "yes" : "no"
)

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
        set_quantity!(
            st, phase,
            α_clinker * BINDER_G * clinker_frac * frac / molar_mass(phase) * u"mol"
        )
    end
    set_quantity!(st, "Gp", BINDER_G * GYPSUM / molar_mass("Gp") * u"mol")
    set_quantity!(st, "H2O@", BINDER_G * WB / molar_mass("H2O@") * u"mol")

    clinker = Float64.(cs.SM.A) * ustrip.(us"mol", st.n)
    # The alkalis leave the grain as it dissolves: same fraction as the clinker.
    clinker .+= oxide_budget(
        ALKALIS, cs.SM.primaries;
        mass = BINDER_G * clinker_frac * α_clinker * u"g"
    )
    slag = oxide_budget(
        SLAG, cs.SM.primaries;
        mass = BINDER_G * SLAG_FRACTION * α_slag * u"g"
    )
    ash = oxide_budget(
        FLYASH, cs.SM.primaries;
        mass = BINDER_G * ASH_FRACTION * α_ash * u"g"
    )
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

@printf(
    "of 100 g of binder: clinker %.1f g, slag %.1f g, fly ash %.1f g, gypsum %.1f g\n",
    100clinker_frac, 100SLAG_FRACTION, 100ASH_FRACTION, 100GYPSUM
)
@printf(
    "of which reacted   : clinker %.1f g, slag %.1f g, fly ash %.1f g, gypsum %.1f g\n\n",
    100clinker_frac * ALPHA_CLINKER, 100SLAG_FRACTION * ALPHA_SLAG,
    100ASH_FRACTION * ALPHA_ASH, 100GYPSUM
)
@printf("%-8s %10s %10s %10s %10s\n", "", "clinker", "slag", "fly ash", "total")
for (i, c) in enumerate(components)
    abs(b[i]) > 1.0e-6 &&
        @printf(
        "%-8s %10.5f %10.5f %10.5f %10.5f\n",
        c, p28.clinker[i], p28.slag[i], p28.ash[i], b[i]
    )
end

# Continued from a seven-day paste rather than started cold -- see `solve_paste`
# above for why, and section 6 for the ages this walk passes through.
eq_early, cert_early = solve_paste(0.35, 0.1)
eq, cert = solve_paste(
    ALPHA_SLAG, ALPHA_ASH;
    start = cert_early.optimal ? eq_early : nothing
)

@printf("started from a 7-day paste: certified %s\n", cert_early.optimal)
@printf(
    "certificate: optimal=%s  worst SI=%.2e  element balance=%.1e\n",
    cert.optimal, cert.worst_supersaturation, cert.balance
)
# The volume is that of the REACTED system and its pore solution. The
# unhydrated clinker and the undissolved glass occupy volume too; the clinker's
# is computable from its phases, the glass's would need a density this page has
# not been given, so neither is added rather than one of them being.
@printf(
    "pH = %.3f   volume of the reacted system = %.2f cm3\n",
    pH(eq, model), ustrip(uconvert(us"cm^3", volume(eq).total))
)

n = ustrip.(us"mol", eq.n)
present = sort(
    [
        (symbol(cs.species[i]), n[i])
            for i in cs.idx_crystal if n[i] > 1.0e-4
    ]; by = last, rev = true
)
for (name, amount) in present
    @printf("  %-18s %9.5f mol\n", name, amount)
end

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

r = half_reaction(eq, "SO4-2", "HS-")
println("half-reaction : ", r.equation)
@printf("log K at 25 C : %.2f\n", r.logK⁰(T = 298.15))
@printf("pe            : %+.2f\n", pe(eq, model))
@printf("Eh            : %+.3f V\n", Eh(eq, model))
@printf(
    "aqueous S(VI)  : %.3e mol\naqueous S(-II) : %.3e mol\n",
    ustrip(us"mol", moles(eq, "SO4-2")), ustrip(us"mol", moles(eq, "HS-"))
)

AGES = [
    (" 7 days", 0.35, 0.1), ("28 days", ALPHA_SLAG, ALPHA_ASH),
    ("90 days", 0.52, 0.25),
]

@printf(
    "%-9s %7s %7s   %-9s %9s %7s %9s\n",
    "age", "slag", "ash", "certified", "balance", "pH", "portlandite"
)
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
        @printf(
            "%-9s %6.0f %% %6.0f %%   %-9s %9.1e %7.3f %9.5f\n",
            label, 100a_sl, 100a_as, c.optimal, c.balance, pH(e, model),
            nn[i_portlandite]
        )
    end
end

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
    ("smilauer2025-184-cemIII-A-42.5N-hranice.csv", "CEM III/A", 0.5),
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

fig = scatter(
    repl, heats; legend = false, markersize = 7, color = :seagreen,
    xlabel = "clinker replaced (%, assumed at the EN 197-1 midpoint)",
    ylabel = "measured heat at ~28 days (J/g of binder)",
    title = "Every joule comes from the clinker",
    size = (780, 420),
    bottom_margin = 10Plots.mm, left_margin = 10Plots.mm
)
for (x, y, l) in zip(repl, heats, labels)
    annotate!(fig, x, y + 8, text(l, 8, :left, :bottom))
end
plot!(fig; xlims = (-8, 85), ylims = (200, 410))
savefig(fig, "cem5-heat.svg"); nothing # hide
