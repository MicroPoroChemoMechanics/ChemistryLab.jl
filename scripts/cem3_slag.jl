# =============================================================================
#  GENERATED FILE — do not edit.
#
#  Extracted from docs/src/examples/cem3_slag.md, whose executed blocks it reproduces in
#  order. Edit the page, then regenerate:
#
#      julia --project=. scripts/_page_to_script.jl
#
#  `test/scripts.jl` fails if this file and its page disagree.
#
#  Usage:
#      julia --project=docs scripts/cem3_slag.jl
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
using Logging
using OptimaSolver
using OrderedCollections
using Printf
using Plots
default(framestyle = :box, grid = false)

substances = build_species(datapath("cemdata18-thermofun.json"); verbose = false)
byname = Dict(symbol(s) => s for s in substances)
nothing # hide

rec = joinpath(datapath("experimental"), "smilauer2025-184-cemIII-A-42.5N-hranice.csv")
for line in eachline(rec)
    startswith(line, "#") || break
    occursin(r"cement:|blaine:|wb:|temperature:", line) && println(line)
end

# ASSUMED, not measured: the midpoint of the EN 197-1 range for CEM III/A,
# which is 35-64 % clinker and 36-65 % slag.
CLINKER_FRACTION = 0.5
SLAG_FRACTION = 0.5

# ASSUMED: a Bogue composition representative of a CEM I clinker. The deposit
# does not report one for this cement.
CLINKER = OrderedDict("C3S" => 0.65, "C2S" => 0.11, "C3A" => 0.11, "C4AF" => 0.08)

# ASSUMED: a European ground granulated blastfurnace slag analysis. The sulfur
# is the part that matters here, and it is the part a datasheet reports least
# consistently.
SLAG = Dict(
    "CaO" => 0.41, "SiO2" => 0.36, "Al2O3" => 0.11,
    "MgO" => 0.08, "SO3" => 0.02
)

# THE CLINKER'S ALKALIS, which Bogue does not account for and which set the pH.
#
# A Bogue calculation returns four phases and no sodium or potassium: they are
# minor oxides, a fraction of a percent, and they sit outside the four-phase
# decomposition. They are also, in a cement paste, **what fixes the pH** -- they
# dissolve almost completely into the pore solution and stay there, where the
# calcium is held down by portlandite at 12.5. Leaving them out does not make the
# calculation conservative: it makes it report a portlandite floor as though it
# were a pore solution.
#
# ASSUMED at a usual industrial level, as a fraction of the CLINKER mass.
ALKALIS = Dict("K2O" => 0.008, "Na2O" => 0.002)

WB = 0.4            # MEASURED, from the record above
BINDER_G = 100.0

# HOW MUCH REACTS. Two ceilings, and the reacted fraction is the lower.
#
# The WATER ceiling is Powers (1948): about 0.42 g of water per gram of cement
# is needed for complete hydration -- 0.23 g written into the hydrate formulae
# and 0.19 g held in the gel pores those hydrates create. Below that the paste
# stops WITH WATER STILL IN IT, what remains being in pores far too fine to
# reach an unhydrated grain. That argument is about the pore space and not about
# the grain, so it caps the slag exactly as it caps the alite. Water curing
# moves it to 0.36, the volume emptied by chemical shrinkage being refilled from
# outside -- `powers_alpha_max(WB; curing = :saturated)`.
ALPHA_WATER = powers_alpha_max(WB)                 # sealed, ASSUMED

# The KINETIC ceiling is the slag's own dissolution rate, and at 28 days it is
# the lower of the two by a wide margin. The RILEM TC 238-SCM round robin
# [Durdzinski2017](@cite) measured two ground granulated slags at 40 %
# replacement and w/b 0.40 in seven laboratories -- this page's geometry. Its
# Table 4 at 28 days, by SEM image analysis, the technique the study found most
# consistent: 38 % and 48 % for the first slag, 45 % and 49 % for the second.
# The study's verdict on the precision of any technique: "at best +/- 5 %".
#
# ASSUMED from that table. [The CEM V page](@ref cem5-dor) sweeps the same
# quantity across the round robin's 7-, 28- and 90-day columns.
ALPHA_SLAG = min(0.45, ALPHA_WATER)
ALPHA_CLINKER = ALPHA_WATER

@printf("water ceiling at w/b = %.2f : %.3f\n", WB, ALPHA_WATER)
@printf("reacted: clinker %.0f %%, slag %.0f %%\n", 100ALPHA_CLINKER, 100ALPHA_SLAG)

pure = split(
    "C3S C2S C3A C4AF Gp Anh Cal Portlandite ettringite monosulphate12 " *
        "monocarbonate hemicarbonate C4AH13 C3AH6 C3FH6 straetlingite " *
        "hydrotalcite Mg2AlC0.5OH Brc FeOOHmic AlOHmic Amor-Sl Tro Py Sulfur Mgs"
)
solutions = [
    "CSHQ" => ["CSHQ-JenD", "CSHQ-JenH", "CSHQ-TobD", "CSHQ-TobH", "KSiOH", "NaSiOH"],
    "C3(AF)S0.84H" => ["C3AFS0.84H4.32", "C3FS0.84H4.32"],
]
members = reduce(vcat, last.(solutions))

# The sulfur ladder must be in the species list, or "holding both oxidation
# states" is a claim about species that are not there.
redox_species = ["HS-", "H2S@", "SO4-2", "SO3-2", "S2O3-2", "O2@", "H2@"]

species = speciation(
    substances, vcat(pure, members, redox_species);
    aggregate_state = [AS_AQUEOUS]
)
ss = [SolidSolutionPhase(n, [byname[m] for m in ms]) for (n, ms) in solutions]
cs = ChemicalSystem(species, CEMDATA_PRIMARIES; solid_solutions = ss)

@printf(
    "%d species, %d conservation components\n",
    length(cs.species), size(cs.SM.A, 1)
)

components = String.(symbol.(cs.SM.primaries))
@printf(
    "charge (`Zz`) kept as a conservation component: %s\n",
    "Zz" in components ? "yes — the oxidation state is conserved" : "no"
)

molar_mass(n) = ustrip(us"g/mol", byname[n][:M])

"""
    paste(; alkali = 1.0) -> (; state, clinker, slag, total)

The fresh state and the element budget, with `alkali` scaling the clinker's
sodium and potassium so that section 6 can sweep the one input that sets the pH.

The state carries only what has reacted. The rest -- unhydrated clinker cores and
undissolved glass -- is still in the specimen, and is no part of the
minimization: an intact grain is not at equilibrium with the solution around it.
All of the mixing water enters, though: the ceiling limits how far the reaction
can go, not how much water was poured in, and the water that cannot reach a grain
is still in the balance.
"""
function paste(; alkali = 1.0)
    st = ChemicalState(cs)
    for (phase, frac) in CLINKER
        set_quantity!(
            st, phase,
            ALPHA_CLINKER * BINDER_G * CLINKER_FRACTION * frac / molar_mass(phase) * u"mol"
        )
    end
    set_quantity!(st, "H2O@", BINDER_G * WB / molar_mass("H2O@") * u"mol")

    clinker = Float64.(cs.SM.A) * ustrip.(us"mol", st.n)
    # The alkalis follow the clinker, and its reacted fraction: they leave the
    # grain as it dissolves.
    clinker .+= alkali * oxide_budget(
        ALKALIS, cs.SM.primaries;
        mass = BINDER_G * CLINKER_FRACTION * ALPHA_CLINKER * u"g"
    )
    slag = oxide_budget(
        SLAG, cs.SM.primaries;
        mass = BINDER_G * SLAG_FRACTION * ALPHA_SLAG * u"g"
    )
    return (; state = st, clinker, slag, total = clinker .+ slag)
end

p0 = paste()
state, b = p0.state, p0.total

for (comp, v) in zip(components, b)
    abs(v) > 1.0e-6 && @printf("  %-8s %10.5f mol\n", comp, v)
end

model = HKFActivityModel(å = 0.0, Ḃ = 0.097637, Kₙ = 0.0)
eq, cert = equilibrate_certified(state; model = model, b = b)

@printf(
    "certificate: optimal=%s  worst SI=%.2e  element balance=%.1e\n",
    cert.optimal, cert.worst_supersaturation, cert.balance
)
@printf("pH = %.3f\n", pH(eq, model))

n = ustrip.(us"mol", eq.n)
present = sort(
    [(symbol(cs.species[i]), n[i]) for i in cs.idx_crystal if n[i] > 1.0e-4];
    by = last, rev = true
)
for (name, amount) in present
    @printf("  %-18s %9.5f mol\n", name, amount)
end

r = half_reaction(eq, "SO4-2", "HS-")
println("half-reaction : ", r.equation)
@printf("log K at 25 C : %.2f\n", r.logK⁰(T = 298.15))
# `pe` warns when a member of the couple sits at the solver's floor; the two
# sulfur amounts printed below say the same, so the warning is not repeated.
pe_eq, eh_eq = with_logger(NullLogger()) do
    pe(eq, model), Eh(eq, model)
end
@printf("pe            : %+.2f\n", pe_eq)
@printf("Eh            : %+.3f V\n", eh_eq)

s6 = ustrip(us"mol", moles(eq, "SO4-2"))
s2 = ustrip(us"mol", moles(eq, "HS-"))
@printf("\naqueous S(VI)  : %.3e mol\naqueous S(-II) : %.3e mol\n", s6, s2)

i_ch = findfirst(sp -> symbol(sp) == "Portlandite", cs.species)
na2o_eq(scale) = 100 * scale * (ALKALIS["Na2O"] + 0.658 * ALKALIS["K2O"])

@printf(
    "%-14s %10s %11s %8s %13s\n",
    "Na2O eq (%)", "certified", "balance", "pH", "portlandite"
)

# `let` rather than a bare loop: a top-level `for` that assigns to a name of the
# enclosing scope makes a NEW LOCAL, so `prev` would be read before it is ever
# written. Wrapping the sweep gives it a scope of its own, which is cleaner than
# reaching for `global`.
let prev = nothing
    for scale in (0.5, 1.0, 1.5)
        pa = paste(; alkali = scale)
        e, c = equilibrate_certified(
            something(prev, pa.state);
            model = model, b = pa.total
        )
        c.optimal && (prev = e)
        nn = ustrip.(us"mol", e.n)
        @printf(
            "%-14.2f %10s %11.1e %8.3f %13.5f\n",
            na2o_eq(scale), c.optimal, c.balance, pH(e, model), nn[i_ch]
        )
    end
end

function final_heat(file)
    t, Q = 0.0, 0.0
    for line in eachline(joinpath(datapath("experimental"), file))
        startswith(line, "#") && continue
        startswith(line, "time") && continue
        f = split(line, ',')
        t, Q = parse(Float64, f[1]), parse(Float64, f[3])
    end
    return t, Q
end

for (file, label) in (
        ("smilauer2025-122-cemI-52.5R-cizkovice.csv", "CEM I 52.5 R"),
        ("smilauer2025-184-cemIII-A-42.5N-hranice.csv", "CEM III/A 42.5 N"),
        ("smilauer2025-121-cemIII-B-32.5N-mokra.csv", "CEM III/B 32.5 N"),
    )
    t, Q = final_heat(file)
    @printf("  %-18s %6.0f h   %6.1f J/g\n", label, t, Q)
end
