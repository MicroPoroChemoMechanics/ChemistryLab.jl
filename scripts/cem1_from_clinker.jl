# =============================================================================
#  GENERATED FILE — do not edit.
#
#  Extracted from docs/src/examples/cem1_from_clinker.md, whose executed blocks it reproduces in
#  order. Edit the page, then regenerate:
#
#      julia --project=. scripts/_page_to_script.jl
#
#  `test/scripts.jl` fails if this file and its page disagree.
#
#  Usage:
#      julia --project=docs scripts/cem1_from_clinker.jl
#
#  No `Pkg.activate` here, deliberately: the active project is global process
#  state and this file is meant to be `include`d as well as run.
#
#  Formatted with Runic like the rest of the repository. `test/scripts.jl`
#  compares SYNTAX TREES rather than text, so formatting is free to differ
#  from the page while a real divergence still fails the suite.
# =============================================================================

using ChemistryLab
using OrdinaryDiffEq
using DynamicQuantities
using OrderedCollections
using Plots
using Printf

# mass fractions of the anhydrous cement. C3S + C2S is held at 0.8126.
CLINKERS = OrderedDict(
    "measured (Baroghel-Bouny CO)" =>
        (C3S = 0.5728, C2S = 0.2398, C3A = 0.0303, C4AF = 0.0759, Gp = 0.0439),
    "alite-rich (constructed)" =>
        (C3S = 0.7126, C2S = 0.1, C3A = 0.0303, C4AF = 0.0759, Gp = 0.0439),
    "belite-rich (constructed)" =>
        (C3S = 0.4128, C2S = 0.3998, C3A = 0.0303, C4AF = 0.0759, Gp = 0.0439),
)
WC = 0.45          # enough water that the arrest is not the subject here

for (name, c) in CLINKERS
    @printf(
        "  %-30s C3S %.4f  C2S %.4f  C3A %.4f  C4AF %.4f  Gp %.4f\n",
        name, c.C3S, c.C2S, c.C3A, c.C4AF, c.Gp
    )
end

substances = build_species(datapath("cemdata18-thermofun.json"); verbose = false)
selection = split(
    "C3S C2S C3A C4AF Gp Portlandite Jennite ettringite monosulphate12 " *
        "C3AH6 C3FH6 H2O@"
)
cs = ChemicalSystem(
    speciation(substances, selection; aggregate_state = [AS_AQUEOUS]), CEMDATA_PRIMARIES
)
sp(n) = cs[n]
println("$(length(cs.species)) species in the system")

law = VanGenuchten(; a = 37.5479e6, m = 1 / 2.1684)   # Baroghel-Bouny, mix CO

function build(compo; humidity = true, tend = 90 * 86400.0)
    st = ChemicalState(cs)
    for (name, frac) in pairs(compo)
        set_quantity!(st, string(name), frac * u"kg")
    end
    set_quantity!(st, "H2O@", WC * sum(values(compo)) * u"kg")
    h = humidity ? PoreHumidity(law, cs; reference = st) : nothing

    specs = (
        (
            "C3S", PK84_PARAMS_C3S,
            OrderedDict(sp("C3S") => 1.0, sp("H2O@") => 103 / 30),
            OrderedDict(sp("Jennite") => 1.0, sp("Portlandite") => 4 / 3),
        ),
        (
            "C2S", PK84_PARAMS_C2S,
            OrderedDict(sp("C2S") => 1.0, sp("H2O@") => 73 / 30),
            OrderedDict(sp("Jennite") => 1.0, sp("Portlandite") => 1 / 3),
        ),
        # No gypsum in this one, deliberately — see §9.
        (
            "C3A", PK84_PARAMS_C3A,
            OrderedDict(sp("C3A") => 1.0, sp("H2O@") => 6.0),
            OrderedDict(sp("C3AH6") => 1.0),
        ),
        (
            "C4AF", PK84_PARAMS_C4AF,
            OrderedDict(sp("C4AF") => 1.0, sp("Portlandite") => 2.0, sp("H2O@") => 10.0),
            OrderedDict(sp("C3AH6") => 1.0, sp("C3FH6") => 1.0),
        ),
    )
    rxns = AbstractReaction[]
    for (nm, pk, reac, prod) in specs
        rx = Reaction(reac, prod; symbol = nm)
        rx[:rate] = parrot_killoh_avrami(
            pk, nm; α_max = powers_alpha_max(WC), blaine = 380.0u"m^2/kg",
            humidity = h
        )
        push!(rxns, rx)
    end
    kp = KineticsProblem(cs, rxns, st, (0.0, tend))
    ks = KineticsSolver(; ode_solver = Rodas5P(), reltol = 1.0e-6, abstol = 1.0e-10)
    return kp, integrate(kp, ks), st, h
end

runs = OrderedDict(name => build(c) for (name, c) in CLINKERS)
println("integrated: ", join(keys(runs), ", "))

kp, sol, st0, h = runs["measured (Baroghel-Bouny CO)"]
ts = exp10.(range(log10(3600.0), log10(90 * 86400.0); length = 120))
dh = degrees_of_hydration(sol, kp; times = ts)
days = ts ./ 86400

p1 = plot(;
    xscale = :log10, xlabel = "time (days)", ylabel = "degree of hydration α",
    title = "Clinker phases, measured composition", legend = :topleft, ylims = (0, 1)
)
for (nm, col) in (
        ("C3S", :steelblue), ("C2S", :seagreen),
        ("C3A", :firebrick), ("C4AF", :darkorange),
    )
    plot!(p1, days, dh[nm]; label = nm, linewidth = 2, color = col)
end
plot(p1; size = (720, 430), left_margin = 8Plots.mm, bottom_margin = 8Plots.mm)

p2 = plot(;
    xscale = :log10, xlabel = "time (days)",
    ylabel = "silicate degree of hydration",
    title = "Alite/belite balance, at constant total silicate",
    legend = :topleft, ylims = (0, 1)
)
for ((name, (kpi, soli, _, _)), col) in zip(runs, (:black, :firebrick, :seagreen))
    d = degrees_of_hydration(soli, kpi; times = ts)
    c = CLINKERS[name]
    blend = (c.C3S .* d["C3S"] .+ c.C2S .* d["C2S"]) ./ (c.C3S + c.C2S)
    plot!(p2, days, blend; label = name, linewidth = 2, color = col)
end
plot(p2; size = (760, 430), left_margin = 8Plots.mm, bottom_margin = 8Plots.mm)

states = [state_at(sol, kp, t) for t in ts]
amount(name) = [
    ustrip(us"mol", s.n[findfirst(x -> symbol(x) == name, cs.species)])
        for s in states
]

p3 = plot(;
    xscale = :log10, xlabel = "time (days)", ylabel = "amount (mol)",
    title = "Hydrate assemblage, measured composition", legend = :topleft
)
for (nm, lbl, col) in (
        ("Portlandite", "portlandite CH", :steelblue),
        ("Jennite", "C-S-H (Jennite)", :seagreen),
        ("C3AH6", "hydrogarnet C3AH6", :firebrick),
        ("C3FH6", "ferrite hydrogarnet", :darkorange),
    )
    plot!(p3, days, amount(nm); label = lbl, linewidth = 2, color = col)
end
plot(p3; size = (720, 430), left_margin = 8Plots.mm, bottom_margin = 8Plots.mm)

V0 = volume(st0)
poro = [porosity(s, st0).total for s in states]
# m³ → cm³ by hand: `us"cm^3"` is a symbolic unit and does not convert a
# quantity carried in SI dimensions.
shrink = [1.0e6 * ustrip(us"m^3", V0.total - volume(s).total) for s in states]
rh = [h(ustrip.(us"mol", s.n)) for s in states]

p4 = plot(
    days, 100 .* poro; xscale = :log10, xlabel = "time (days)",
    ylabel = "total porosity (%)", label = "porosity", linewidth = 2,
    color = :steelblue, title = "Pore space and self-desiccation", legend = :left
)
p4b = twinx(p4)
plot!(
    p4b, days, rh; xscale = :log10, ylabel = "internal relative humidity",
    label = "internal RH", linewidth = 2, color = :firebrick, legend = :right,
    ylims = (0.75, 1.01)
)
hline!(p4b, [0.8]; label = "humidity_factor cut", linestyle = :dash, color = :gray)
plot(
    p4; size = (760, 430), left_margin = 10Plots.mm, right_margin = 12Plots.mm,
    bottom_margin = 8Plots.mm
)

p5 = plot(
    days, shrink; xscale = :log10, xlabel = "time (days)",
    ylabel = "volume lost (cm³)", label = "chemical shrinkage", linewidth = 2,
    color = :seagreen, title = "Chemical shrinkage", legend = :topleft
)
plot(p5; size = (720, 400), left_margin = 10Plots.mm, bottom_margin = 8Plots.mm)
