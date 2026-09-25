# =============================================================================
#  GENERATED FILE — do not edit.
#
#  Extracted from docs/src/examples/cem2_blended.md, whose executed blocks it reproduces in
#  order. Edit the page, then regenerate:
#
#      julia --project=. scripts/_page_to_script.jl
#
#  `test/scripts.jl` fails if this file and its page disagree.
#
#  Usage:
#      julia --project=docs scripts/cem2_blended.jl
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
molar_mass(n) = ustrip(us"g/mol", byname[n][:M])
nothing # hide

function header(file)
    for line in eachline(joinpath(datapath("experimental"), file))
        startswith(line, "#") || break
        occursin(r"cement:|blaine:|wb:|temperature:", line) && println(line)
    end
    return println()
end

header("smilauer2025-165-cemII-A-LL-42.5R-hranice.csv")
header("smilauer2025-149-cemII-B-S-32.5R-mokra.csv")

# ASSUMED: a Bogue composition representative of a CEM I clinker. The deposit
# reports none.
CLINKER = OrderedDict("C3S" => 0.65, "C2S" => 0.11, "C3A" => 0.11, "C4AF" => 0.08)

# ASSUMED: midpoints of the EN 197-1 ranges. CEM II/A-LL is 80-94 % clinker
# with 6-20 % limestone; CEM II/B-S is 65-79 % clinker with 21-35 % slag.
LL_LIMESTONE = 0.13
BS_SLAG = 0.28

# ASSUMED: a European ground granulated blastfurnace slag analysis.
SLAG = OrderedDict(
    "CaO" => 0.41, "SiO2" => 0.36, "Al2O3" => 0.11,
    "MgO" => 0.08, "SO3" => 0.02
)

# THE CLINKER'S ALKALIS, which Bogue does not account for and which set the pH.
# They are minor oxides, a fraction of a percent, sitting outside the four-phase
# decomposition -- and in a cement paste they are what fixes the pH, dissolving
# almost completely into the pore solution and staying there while the calcium is
# held down by portlandite at 12.5. Omitting them does not make the calculation
# conservative: it makes it report a portlandite floor as a pore solution.
# ASSUMED at a usual industrial level, as a fraction of the CLINKER mass.
# [The CEM III page](@ref cem3-alkali) sweeps this over the industrial range and
# shows that it, and almost nothing else, is what moves the pH.
ALKALIS = OrderedDict("K2O" => 0.008, "Na2O" => 0.002)

# The calcium sulfate ground in with every Portland clinker, as a mass fraction
# of the binder. ASSUMED at a usual industrial level.
GYPSUM = 0.046

BINDER_G = 100.0

# HOW MUCH REACTS, and it is not everything. Two ceilings; the reacted fraction
# is the lower of them.
#
# The WATER ceiling, Powers (1948): complete hydration needs about 0.42 g of
# water per gram of cement, 0.23 g written into the hydrates and 0.19 g held in
# the gel pores those hydrates create. Below that the paste stops with water
# still in it, in pores too fine to reach an unhydrated grain. It is a property
# of the pore space, so it caps the slag as it caps the clinker; and it depends
# on the mix, which is why the two pastes here do not get the same number --
# the limestone paste was mixed at w/b 0.45 and the slag paste at 0.40.
#
# The KINETIC ceiling, for the slag only: [Durdzinski2017](@cite), Table 4, a
# ground granulated slag at 28 days by SEM image analysis, 38-49 % across two
# slags and two laboratories, with a stated precision of "at best +/- 5 %".
# ASSUMED at 45 %. The limestone needs none: calcite is a declared phase, and
# the minimization dissolves exactly as much of it as is stable.
ALPHA_SLAG_KINETIC = 0.45
reacted(wb) = (
    clinker = powers_alpha_max(wb),
    slag = min(ALPHA_SLAG_KINETIC, powers_alpha_max(wb)),
)

for wb in (0.45, 0.4)
    r = reacted(wb)
    @printf(
        "w/b = %.2f : water ceiling %.3f -> clinker %.0f %%, slag %.0f %%\n",
        wb, powers_alpha_max(wb), 100r.clinker, 100r.slag
    )
end

pure = split(
    "C3S C2S C3A C4AF Gp Anh Cal Portlandite ettringite monosulphate12 " *
        "monocarbonate hemicarbonate monocarbonate9 hemicarbonat10.5 " *
        "C4AH13 C3AH6 C3FH6 straetlingite Femonocarbonate Fe-hemicarbonate " *
        "hydrotalcite Mg2AlC0.5OH Brc FeOOHmic AlOHmic Amor-Sl Tro Py Sulfur Mgs"
)
solutions = [
    "CSHQ" => ["CSHQ-JenD", "CSHQ-JenH", "CSHQ-TobD", "CSHQ-TobH", "KSiOH", "NaSiOH"],
    "C3(AF)S0.84H" => ["C3AFS0.84H4.32", "C3FS0.84H4.32"],
]
members = reduce(vcat, last.(solutions))
redox_species = ["HS-", "H2S@", "SO4-2", "SO3-2", "S2O3-2", "O2@", "H2@", "CO2@"]

species = speciation(
    substances, vcat(pure, members, redox_species);
    aggregate_state = [AS_AQUEOUS]
)
ss = [SolidSolutionPhase(n, [byname[m] for m in ms]) for (n, ms) in solutions]
cs = ChemicalSystem(species, CEMDATA_PRIMARIES; solid_solutions = ss)
model = HKFActivityModel(å = 0.0, Ḃ = 0.097637, Kₙ = 0.0)

components = String.(symbol.(cs.SM.primaries))
@printf(
    "%d species, %d conservation components: %s\n",
    length(cs.species), length(components), join(components, " ")
)

"""Element budget of one paste, in moles per 100 g of binder."""
function budget(; clinker_frac, limestone = 0.0, slag = 0.0, wb)
    α = reacted(wb)
    state = ChemicalState(cs)
    # Only the reacted clinker is posed to the minimization; the unhydrated
    # cores are still in the specimen but are not at equilibrium with it.
    for (phase, frac) in CLINKER
        set_quantity!(
            state, phase,
            α.clinker * BINDER_G * clinker_frac * frac / molar_mass(phase) * u"mol"
        )
    end
    # The calcium sulfate is soluble and the limestone is a declared phase, so
    # neither carries a ceiling. All of the mixing water enters: the ceiling
    # says how far the reaction goes, not how much water was poured in.
    set_quantity!(state, "Gp", BINDER_G * GYPSUM / molar_mass("Gp") * u"mol")
    limestone > 0 && set_quantity!(
        state, "Cal",
        BINDER_G * limestone / molar_mass("Cal") * u"mol"
    )
    set_quantity!(state, "H2O@", BINDER_G * wb / molar_mass("H2O@") * u"mol")

    b = Float64.(cs.SM.A) * ustrip.(us"mol", state.n)
    # The alkalis leave the grain as it dissolves, so they follow the clinker and
    # its reacted fraction.
    b .+= oxide_budget(
        ALKALIS, cs.SM.primaries;
        mass = BINDER_G * clinker_frac * α.clinker * u"g"
    )
    slag > 0 && (
        b .+= oxide_budget(
            SLAG, cs.SM.primaries;
            mass = BINDER_G * slag * α.slag * u"g"
        )
    )
    return state, b
end

# The clinker fraction is what is left once the replacement and the gypsum are
# taken out, so the three pastes really are 100 g of binder each.
st_ll, b_ll = budget(
    clinker_frac = 1 - LL_LIMESTONE - GYPSUM,
    limestone = LL_LIMESTONE, wb = 0.45
)
st_ref, b_ref = budget(clinker_frac = 1 - LL_LIMESTONE - GYPSUM, wb = 0.45)
st_bs, b_bs = budget(
    clinker_frac = 1 - BS_SLAG - GYPSUM,
    slag = BS_SLAG, wb = 0.4
)

@printf("%-10s %s\n", "", join((@sprintf("%8s", c) for c in components), ""))
for (label, b) in (
        "CEM II/A-LL" => b_ll, "no limestone" => b_ref,
        "CEM II/B-S" => b_bs,
    )
    @printf("%-12s%s\n", label, join((@sprintf("%8.3f", v) for v in b), ""))
end

function solve(state, b)
    eq, cert = equilibrate_certified(state; model = model, b = b)
    return eq, cert
end

eq_ll, c_ll = solve(st_ll, b_ll)
eq_ref, c_ref = solve(st_ref, b_ref)
eq_bs, c_bs = solve(st_bs, b_bs)

for (label, eq, c) in (
        ("CEM II/A-LL", eq_ll, c_ll),
        ("no limestone", eq_ref, c_ref),
        ("CEM II/B-S", eq_bs, c_bs),
    )
    @printf(
        "%-13s optimal=%-5s  worst SI=%+.2e  balance=%.1e  pH=%.3f  V=%.2f cm3\n",
        label, c.optimal, c.worst_supersaturation, c.balance,
        pH(eq, model), ustrip(uconvert(us"cm^3", volume(eq).total))
    )
end

function assemblage(eq; tol = 1.0e-4)
    n = ustrip.(us"mol", eq.n)
    return sort(
        [(symbol(cs.species[i]), n[i]) for i in cs.idx_crystal if n[i] > tol];
        by = last, rev = true
    )
end

for (label, eq) in ("CEM II/A-LL" => eq_ll, "no limestone" => eq_ref)
    println(label, ":")
    for (name, amount) in assemblage(eq)
        @printf("  %-18s %9.5f mol\n", name, amount)
    end
    println()
end

aluminates = [
    "ettringite", "monosulphate12", "monocarbonate", "hemicarbonate",
    "C4AH13", "straetlingite",
]
n_ll = ustrip.(us"mol", eq_ll.n)
n_ref = ustrip.(us"mol", eq_ref.n)
idx = Dict(symbol(s) => i for (i, s) in enumerate(cs.species))
get_n(n, s) = haskey(idx, s) ? n[idx[s]] : 0.0

v_ref = [get_n(n_ref, s) for s in aluminates]
v_ll = [get_n(n_ll, s) for s in aluminates]
top = 1.1 * maximum(vcat(v_ref, v_ll))

p_ref = bar(
    aluminates, v_ref; legend = false, color = :steelblue, ylims = (0, top),
    ylabel = "mol per 100 g of binder", xrotation = 30,
    title = "no limestone"
)
p_ll = bar(
    aluminates, v_ll; legend = false, color = :seagreen, ylims = (0, top),
    xrotation = 30, title = "CEM II/A-LL, 13 % calcite"
)
fig = plot(
    p_ref, p_ll; layout = (1, 2), size = (900, 420),
    bottom_margin = 16Plots.mm, left_margin = 8Plots.mm,
    plot_title = "The carbonate moves the sulfate out of the AFm and into the AFt",
    plot_titlefontsize = 11
)
savefig(fig, "cem2-aluminates.svg"); nothing # hide

println("CEM II/B-S:")
for (name, amount) in assemblage(eq_bs)
    @printf("  %-18s %9.5f mol\n", name, amount)
end

# `pe` warns when a member of the couple sits at the solver's floor; the two
# sulfur amounts printed below say the same, so the warning is not repeated.
pe_bs, eh_bs = with_logger(NullLogger()) do
    pe(eq_bs, model), Eh(eq_bs, model)
end
@printf("\npe  = %+.2f\nEh  = %+.3f V\n", pe_bs, eh_bs)
@printf(
    "aqueous S(VI)  : %.3e mol\naqueous S(-II) : %.3e mol\n",
    ustrip(us"mol", moles(eq_bs, "SO4-2")), ustrip(us"mol", moles(eq_bs, "HS-"))
)

labels = ["CEM II/A-LL", "no limestone", "CEM II/B-S"]
eqs = [eq_ll, eq_ref, eq_bs]
vols = [ustrip(uconvert(us"cm^3", volume(e).total)) for e in eqs]
phs = [pH(e, model) for e in eqs]

p1 = bar(
    labels, vols; legend = false, ylabel = "total volume (cm³)",
    color = :seagreen, title = "Volume of the hydrated paste", xrotation = 15
)
p2 = bar(
    labels, phs; legend = false, ylabel = "pH", ylims = (12.0, 13.6),
    color = :steelblue, title = "Pore solution pH", xrotation = 15
)
fig2 = plot(
    p1, p2; layout = (1, 2), size = (900, 400),
    bottom_margin = 14Plots.mm, left_margin = 8Plots.mm
)
savefig(fig2, "cem2-summary.svg"); nothing # hide

function final_heat(file)
    t, Q = 0.0, 0.0
    for line in eachline(joinpath(datapath("experimental"), file))
        (startswith(line, "#") || startswith(line, "time")) && continue
        f = split(line, ',')
        t, Q = parse(Float64, f[1]), parse(Float64, f[3])
    end
    return t, Q
end

for (file, label) in (
        ("smilauer2025-122-cemI-52.5R-cizkovice.csv", "CEM I 52.5 R"),
        ("smilauer2025-165-cemII-A-LL-42.5R-hranice.csv", "CEM II/A-LL 42.5 R"),
        ("smilauer2025-149-cemII-B-S-32.5R-mokra.csv", "CEM II/B-S 32.5 R"),
    )
    t, Q = final_heat(file)
    @printf("  %-20s %6.0f h   %6.1f J/g\n", label, t, Q)
end
