# =============================================================================
#  GENERATED FILE — do not edit.
#
#  Extracted from docs/src/examples/miscibility_gap.md, whose executed blocks it reproduces in
#  order. Edit the page, then regenerate:
#
#      julia --project=. scripts/_page_to_script.jl
#
#  `test/scripts.jl` fails if this file and its page disagree.
#
#  Usage:
#      julia --project=docs scripts/miscibility_gap.jl
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

RT = R_GAS * 298.15

# Cemdata18's dimensionless Guggenheim parameters, from
# data/literature/Lothenbach2019.json: a₀ = α₀RT and a₁ = α₁RT.
function published_rk(pair)
    p = literature_row("Lothenbach2019", "guggenheim_parameters", pair)
    return RedlichKisterModel(a0 = p.alpha0 * RT, a1 = p.alpha1 * RT)
end

models = OrderedDict(
    "ideal" => IdealSolidSolutionModel(),
    "AFm SO4/OH (published)" => published_rk("AFm SO4/OH"),
    "AFt SO4/CO3 (published)" => published_rk("AFt SO4/CO3"),
    "regular, W = 1.9 RT" => RegularSolutionModel([0.0 1.9RT; 1.9RT 0.0]),
    "regular, W = 2.1 RT" => RegularSolutionModel([0.0 2.1RT; 2.1RT 0.0]),
)

for (label, m) in models
    gap = spinodal_interval(m, 2)
    if gap === nothing
        @printf("  %-26s convex everywhere\n", label)
    else
        @printf("  %-26s concave for x in [%.3f, %.3f]\n", label, gap[1], gap[2])
    end
end

xs = range(0.001, 0.999; length = 400)
g(m, x) = begin
    # molar Gibbs energy of mixing, in units of RT
    ideal = x * log(x) + (1 - x) * log(1 - x)
    excess = if m isa RedlichKisterModel
        x * (1 - x) * (m.a0 + m.a1 * (2x - 1) + m.a2 * (2x - 1)^2) / RT
    else
        0.0
    end
    ideal + excess
end

afm = models["AFm SO4/OH (published)"]
gap = spinodal_interval(afm, 2)

fig = plot(
    xs, [g(afm, x) for x in xs]; label = "published Redlich-Kister",
    color = :firebrick, linewidth = 2,
    xlabel = "x (sulfate end-member)", ylabel = "g / RT",
    title = "The mixing energy of the AFm sulfate/hydroxide binary",
    size = (760, 420), bottom_margin = 8Plots.mm, left_margin = 8Plots.mm
)
plot!(
    fig, xs, [g(IdealSolidSolutionModel(), x) for x in xs];
    label = "ideal mixing", color = :steelblue, linewidth = 2, linestyle = :dash
)
vspan!(fig, [gap[1], gap[2]]; color = :firebrick, alpha = 0.12, label = "spinodal")
savefig(fig, "gap-energy.svg"); nothing # hide

# The Bogue composition of the CEM I 52.5 N of [Lavergne2018](@cite), Table 9,
# and its gypsum.
bogue = literature_table("Lavergne2018", "cement_bogue")
CLINKER = OrderedDict(zip(bogue.phase, bogue.percent ./ 100))
GYPSUM = literature_value("Lavergne2018", "gypsum_percent") / 100
WB = 0.5
BINDER_G = 100.0

pure = split(
    "C3S C2S C3A C4AF Gp Anh Cal Portlandite ettringite monocarbonate " *
        "hemicarbonate C3AH6 C3FH6 straetlingite FeOOHmic AlOHmic Amor-Sl Brc"
)
CSHQ = ["CSHQ-JenD", "CSHQ-JenH", "CSHQ-TobD", "CSHQ-TobH", "KSiOH", "NaSiOH"]
AFM = ["C4AH13", "monosulphate12"]
aqueous = ["SO4-2", "CO2@"]

# Debye-Hückel limiting law with a B-dot term, as GEM-Selektor runs CEMDATA18.
using JSON
# The B-dot is GEM-Selektor's, identified from the activity coefficients it
# printed on a CEMDATA18 Portland paste (test/reference/gems_cemdata18_portland.json):
# the two lowest charge classes fix the limiting-law slope and the B-dot, about 0.0976.
gems = JSON.parsefile(joinpath(pkgdir(ChemistryLab), "test", "reference", "gems_cemdata18_portland.json"))
lg1, lg2 = log10(gems["gamma"]["z1"]), log10(gems["gamma"]["z2"])
Ḃ_gems = (lg1 + (lg1 - lg2) / 3) / gems["ionic_strength_mol_per_kg"]
model = HKFActivityModel(å = 0.0, Ḃ = Ḃ_gems, Kₙ = 0.0)

"""
One system, differing only in how the AFm binary is declared.

`autostart` is a keyword because two of the three cases below are *meant* not to
certify. For those the multi-start cascade -- every backend, then the ideal
pre-solve, then the homotopy continuation -- cannot help: there is no admissible
single-composition minimum for it to find, so it spends minutes arriving at the
answer the first route already gave. The certificate reported is the same one
either way; declining the cascade declines only the search for a better start.
"""
function run_case(afm_phase; autostart = true)
    sp = speciation(
        substances, vcat(pure, CSHQ, AFM, aqueous);
        aggregate_state = [AS_AQUEOUS]
    )
    ss = [SolidSolutionPhase("CSHQ", [byname[m] for m in CSHQ]), afm_phase]
    cs = ChemicalSystem(sp, CEMDATA_PRIMARIES; solid_solutions = ss)

    st = ChemicalState(cs)
    for (phase, frac) in CLINKER
        set_quantity!(
            st, phase,
            BINDER_G * (1 - GYPSUM) * frac / molar_mass(phase) * u"mol"
        )
    end
    set_quantity!(st, "Gp", BINDER_G * GYPSUM / molar_mass("Gp") * u"mol")
    set_quantity!(st, "H2O@", BINDER_G * WB / molar_mass("H2O@") * u"mol")
    b = Float64.(cs.SM.A) * ustrip.(us"mol", st.n)

    # Every case prints its certificate; the warning of a refusal would only repeat it.
    eq, cert = with_logger(NullLogger()) do
        equilibrate_certified(st; model = model, b = b, autostart = autostart)
    end
    return cs, eq, cert
end

afm_em = [byname[m] for m in AFM]
published = published_rk("AFm SO4/OH")
nothing # hide

cs1, eq1, c1 = run_case(SolidSolutionPhase("AFm_SO4_OH", afm_em))
@printf(
    "optimal=%-5s  worst SI=%+.2e  balance=%.1e  pH=%.4f\n",
    c1.optimal, c1.worst_supersaturation, c1.balance, pH(eq1, model)
)

err = try
    SolidSolutionPhase("AFm_SO4_OH", afm_em; model = published)
    nothing
catch e
    e
end
println(err.msg)

cs2, eq2, c2 = run_case(
    SolidSolutionPhase(
        "AFm_SO4_OH", afm_em;
        model = published, check_convexity = false
    );
    autostart = false
)
@printf(
    "optimal=%-5s  worst SI=%+.2e  balance=%.1e  pH=%.4f\n",
    c2.optimal, c2.worst_supersaturation, c2.balance, pH(eq2, model)
)
if hasproperty(c2, :worst_violation_split)
    @printf(
        "phases reported as wanting to split: %d   worst split measure %+.2e\n",
        length(c2.split_phases), c2.worst_violation_split
    )
end

cs3, eq3, c3 = run_case(
    SolidSolutionPhase(
        "AFm_SO4_OH", afm_em;
        model = published, instances = 2
    );
    autostart = false
)
@printf(
    "optimal=%-5s  worst SI=%+.2e  balance=%.1e  pH=%.4f\n",
    c3.optimal, c3.worst_supersaturation, c3.balance, pH(eq3, model)
)

n3 = ustrip.(us"mol", eq3.n)
for (grp, ph) in zip(cs3.ss_groups, cs3.solid_solutions)
    startswith(name(ph), "AFm") || continue
    tot = sum(n3[i] for i in grp)
    x = tot > 1.0e-12 ? n3[grp[2]] / tot : NaN
    @printf(
        "  %-14s total %8.5f mol   x(sulfate) = %s\n",
        name(ph), tot, isnan(x) ? "—" : @sprintf("%.4f", x)
    )
end

ct = common_tangent(published, 2)
sp = spinodal_interval(published, 2)
@printf("common tangent : x = %.4f and %.4f\n", ct[1], ct[2])
@printf(
    "spinodal       : x = %.4f to %.4f  (contained in it, as it must be)\n",
    sp[1], sp[2]
)

A = 3.0                                  # W/RT, well past the threshold of 2
sym = RegularSolutionModel([0.0 A * RT; A * RT 0.0])
a, b = common_tangent(sym, 2)
@printf(
    "symmetric model: x = %.6f and %.6f, and b = 1 - a to %.1e\n",
    a, b, abs(b - (1 - a))
)
@printf(
    "residual of  ln(x/(1-x)) + A(1-2x) = 0  :  %+.2e and %+.2e\n",
    log(a / (1 - a)) + A * (1 - 2a), log(b / (1 - b)) + A * (1 - 2b)
)

# The composition of the AFm of this paste, as the single-composition answer of
# case 2 returns it: x is the fraction of C4AH13, the first end-member.
n2 = ustrip.(us"mol", eq2.n)
grp2 = only(g for (g, ph) in zip(cs2.ss_groups, cs2.solid_solutions) if startswith(name(ph), "AFm"))
x̄_paste = n2[grp2[1]] / sum(n2[i] for i in grp2)

for x̄ in (0.2, 0.4, x̄_paste, 0.8, 0.96)
    r = miscibility_split(published, x̄, 2)
    if r.f_beta == 0
        @printf("x̄ = %.4f : homogeneous (outside the pair)\n", x̄)
    else
        # Two calls rather than one: `@printf` takes a LITERAL format string,
        # and `"a" * "b"` is an expression -- `ArgumentError: First argument to
        # @printf after io must be a format string`, raised while the block is
        # lowered, so it kills the build rather than the line.
        @printf(
            "x̄ = %.4f : %.1f %% at x=%.4f and %.1f %% at x=%.4f",
            x̄, 100r.f_alpha, r.x_alpha, 100r.f_beta, r.x_beta
        )
        @printf(
            "   —  releases %6.1f J/mol   (mass balance %+.0e)\n",
            r.Δg, r.f_alpha * r.x_alpha + r.f_beta * r.x_beta - x̄
        )
    end
end

labels = ["ideal\nmixing", "published,\none composition", "published,\ntwo compositions"]
certs = [c1, c2, c3]
eqs = [eq1, eq2, eq3]
vols = [ustrip(uconvert(us"cm^3", volume(e).total)) for e in eqs]
phs = [pH(e, model) for e in eqs]
ok = [c.optimal for c in certs]

p1 = bar(
    labels, vols; legend = false, ylabel = "total volume (cm³)",
    color = [o ? :seagreen : :firebrick for o in ok], title = "Volume"
)
p2 = bar(
    labels, phs; legend = false, ylabel = "pH", title = "Pore solution pH",
    color = [o ? :seagreen : :firebrick for o in ok],
    ylims = (12.5, 13.5)
)
fig3 = plot(
    p1, p2; layout = (1, 2), size = (900, 420),
    bottom_margin = 16Plots.mm, left_margin = 8Plots.mm
)
savefig(fig3, "gap-summary.svg"); nothing # hide
