# =============================================================================
#  GENERATED FILE — do not edit.
#
#  Extracted from docs/src/examples/cement_wc_ratio.md, whose executed blocks it reproduces in
#  order. Edit the page, then regenerate:
#
#      julia --project=. scripts/_page_to_script.jl
#
#  `test/scripts.jl` fails if this file and its page disagree.
#
#  Usage:
#      julia --project=docs scripts/cement_wc_ratio.jl
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

substances = build_species(datapath("cemdata18-thermofun.json"))

input_species = split("C3S C2S C3A C4AF Gp Anh Portlandite Jennite H2O@ ettringite monosulphate12 C3AH6 C3FH6 C4FH13")
species = speciation(substances, input_species; aggregate_state = [AS_AQUEOUS])

cs = ChemicalSystem(species, CEMDATA_PRIMARIES)
nothing # hide

cs

compo = ["C3S" => 0.678, "C2S" => 0.166, "C3A" => 0.04, "C4AF" => 0.072, "Gp" => 0.028]
c = sum(last.(compo))   # cement mass fraction (= 0.984 here)

using OptimaSolver

sp_idx = Dict(symbol(s) => i for (i, s) in enumerate(cs.species))
wc_range = range(0.3, 0.6; length = 13)

function fresh_paste(wc)
    w, mtot = wc * c, c + wc * c
    st = ChemicalState(cs)
    for (sym, mfrac) in compo
        set_quantity!(st, sym, mfrac / mtot * u"kg")
    end
    set_quantity!(st, "H2O@", w / mtot * u"kg")
    V = volume(st)
    set_quantity!(st, "H+", 1.0e-7u"mol/L" * V.liquid)   # charge seed, pH-neutral
    set_quantity!(st, "OH-", 1.0e-7u"mol/L" * V.liquid)
    return st
end

pH_vals = Float64[]
ϕ_liquid = Float64[]
ϕ_void = Float64[]
ϕ_total = Float64[]
n_portl = Float64[]
n_mono = Float64[]
n_ett = Float64[]
n_jennite = Float64[]
n_clinker = Float64[]

certified = Bool[]

for wc in wc_range
    fresh = fresh_paste(wc)
    eq, cert = equilibrate_certified(deepcopy(fresh))   # the solve may mutate its argument
    push!(certified, cert.optimal)
    ϕ = porosity(eq, fresh)                  # sealed-curing convention
    amount(sym) = ustrip(eq.n[sp_idx[sym]])
    push!(pH_vals, pH(eq))
    push!(ϕ_liquid, ϕ.liquid)
    push!(ϕ_void, ϕ.void)
    push!(ϕ_total, ϕ.total)
    push!(n_portl, amount("Portlandite"))
    push!(n_mono, amount("monosulphate12"))
    push!(n_ett, amount("ettringite"))
    push!(n_jennite, amount("Jennite"))
    push!(n_clinker, sum(amount(s) for s in ("C3S", "C2S", "C3A", "C4AF")))
end
nothing # hide

using Plots

p1 = plot(
    collect(wc_range), pH_vals;
    xlabel = "w/c ratio", ylabel = "Pore solution pH", label = "pH",
    linewidth = 2, marker = :circle, markersize = 4, color = :steelblue,
    title = "Pore solution pH", ylims = (11.5, 13.5), legend = :bottomright
)

p2 = plot(
    collect(wc_range), ϕ_total .* 100;
    xlabel = "w/c ratio", ylabel = "Porosity (%)", label = "total",
    linewidth = 2, marker = :circle, markersize = 4, color = :firebrick,
    title = "Porosity referred to the fresh volume",
    ylims = (0, 50), legend = :topleft
)
plot!(
    p2, collect(wc_range), ϕ_liquid .* 100;
    label = "water-filled", linewidth = 2, marker = :square, markersize = 3,
    color = :steelblue
)
plot!(
    p2, collect(wc_range), ϕ_void .* 100;
    label = "empty (Le Chatelier)", linewidth = 2, marker = :diamond,
    markersize = 3, color = :seagreen
)

plot(
    p1, p2; layout = (1, 2), left_margin = 8Plots.mm,
    bottom_margin = 8Plots.mm, size = (950, 410)
)

p3 = plot(
    collect(wc_range), n_portl;
    xlabel = "w/c ratio", ylabel = "Amount (mol / kg of paste)",
    label = "Portlandite  Ca(OH)₂", linewidth = 2, marker = :circle,
    markersize = 4, color = :steelblue, title = "Hydrate assemblage at equilibrium",
    legend = :right
)
plot!(
    p3, collect(wc_range), n_jennite;
    label = "Jennite (C-S-H)", linewidth = 2, marker = :square,
    markersize = 3, color = :firebrick
)
plot!(
    p3, collect(wc_range), n_mono;
    label = "monosulphate12 (AFm)", linewidth = 2, marker = :diamond,
    markersize = 3, color = :seagreen
)
plot!(
    p3, collect(wc_range), n_ett;
    label = "ettringite (AFt)", linewidth = 2, marker = :utriangle,
    markersize = 3, color = :darkorange
)
plot(p3; left_margin = 8Plots.mm, bottom_margin = 8Plots.mm, size = (700, 420))

using Printf
@printf "every point certified                    : %s\n" all(certified)
@printf "ettringite, largest value over the scan  : %.3e mol\n" maximum(n_ett)
@printf "clinker left, largest value over the scan: %.3e mol\n" maximum(n_clinker)

# Below the scanned range the water runs out. What the minimum does then is the
# point of the last two columns.
low = [0.15, 0.2, 0.25, 0.28, 0.3]
clinker(st) = sum(
    ustrip(us"kg", st.n[sp_idx[s]] * cs.species[sp_idx[s]][:M])
        for s in ("C3S", "C2S", "C3A", "C4AF")
)
println(" w/c   clinker left (%)   certified   free water (mol)   x(solvent)   I (mol/kg)")
for wc in low
    fr = fresh_paste(wc)
    # `autostart = false`: the point of this table is a configuration that does
    # NOT certify, and the multi-start cascade is exactly what cannot help there
    # -- it would try every backend, then the ideal pre-solve, then the homotopy
    # continuation, and arrive at the same answer several minutes later. The
    # certificate reported below is the same one; only the search for a better
    # start is declined.
    #
    # The warnings are what the table reports; they are not the transcript.
    eq, cert = Base.CoreLogging.with_logger(Base.CoreLogging.NullLogger()) do
        equilibrate_certified(deepcopy(fr); autostart = false)
    end
    @printf(
        "%5.2f   %16.1f   %9s   %16.3e   %10.3f   %10.4g\n",
        wc, 100 * clinker(eq) / clinker(fr), cert.optimal,
        ustrip(us"mol", eq.n[sp_idx["H2O@"]]), solvent_fraction(eq),
        ionic_strength(eq),
    )
end

function arrested(wc, α)
    mtot = c + wc * c
    st = ChemicalState(cs)
    for (sym, mfrac) in compo
        set_quantity!(st, sym, α * mfrac / mtot * u"kg")   # only α reacts
    end
    set_quantity!(st, "H2O@", wc * c / mtot * u"kg")       # all the water
    V = volume(st)
    set_quantity!(st, "H+", 1.0e-7u"mol/L" * V.liquid)
    set_quantity!(st, "OH-", 1.0e-7u"mol/L" * V.liquid)
    return st
end

println(" w/c   alpha   certified   x(solvent)   I (mol/kg)     pH   clinker %   porosity")
for wc in (0.25, 0.3, 0.35, 0.42)
    α = powers_alpha_max(wc)
    fresh = fresh_paste(wc)                 # the volume reference, all of it
    eq, cert = equilibrate_certified(arrested(wc, α))
    # The unreacted clinker is put back for the volume and porosity accounting.
    n = collect(eq.n)
    for (sym, _) in compo
        n[sp_idx[sym]] += (1 - α) * fresh.n[sp_idx[sym]]
    end
    final = ChemicalState(cs, n)
    @printf(
        "%5.2f  %6.3f   %9s   %10.4f   %10.4f  %5.2f   %9.1f   %8.4f\n",
        wc, α, cert.optimal, solvent_fraction(eq), ionic_strength(eq), pH(eq),
        100 * clinker(final) / clinker(fresh), porosity(final, fresh).total,
    )
end

wc_cure = 0.4
fresh_c = fresh_paste(wc_cure)

# Sealed: the ceiling is 0.42-based, the shrinkage volume stays empty.
α_sealed = powers_alpha_max(wc_cure)
eq_s, c_s = equilibrate_certified(arrested(wc_cure, α_sealed))

# Cured: the ceiling is 0.36-based because the bath refills what the reaction
# empties, AND the specimen is genuinely open to water -- `SaturatedCuring`
# holds its total volume at the fresh paste's and reports how much it drank.
α_cured = powers_alpha_max(wc_cure; curing = :saturated)
q = Ref(Float64[])
eq_c, c_c = equilibrate_certified(
    arrested(wc_cure, α_cured);
    constraint = SaturatedCuring(; reference = fresh_c), parameters = q,
)

function with_unreacted(eq, α)
    n = collect(eq.n)
    for (sym, _) in compo
        n[sp_idx[sym]] += (1 - α) * fresh_c.n[sp_idx[sym]]
    end
    return ChemicalState(cs, n)
end

@printf(
    "%-10s %7s %10s %12s %14s %12s\n",
    "curing", "alpha", "certified", "pH", "water in (mol)", "void"
)
@printf(
    "%-10s %7.3f %10s %12.2f %14s %12.4f\n",
    "sealed", α_sealed, c_s.optimal, pH(eq_s), "—",
    porosity(with_unreacted(eq_s, α_sealed), fresh_c).void
)
@printf(
    "%-10s %7.3f %10s %12.2f %14.4f %12.4f\n",
    "under water", α_cured, c_c.optimal, pH(eq_c), only(q[]),
    porosity(with_unreacted(eq_c, α_cured), fresh_c).void
)

# The molar mass comes from the database, never from a table typed here.
M_H2O = ustrip(us"g/mol", cs.species[sp_idx["H2O@"]][:M])
m_cement = 1000 / (1 + wc_cure)         # g: `fresh_paste` normalizes to 1 kg of paste

# The thermodynamic shrinkage: a difference of molar volumes, per gram REACTED.
ΔV = ustrip(uconvert(us"cm^3", volume(fresh_c).total)) -
    ustrip(uconvert(us"cm^3", volume(with_unreacted(eq_s, α_sealed)).total))
shrink_vol = ΔV / (α_sealed * m_cement)

# The same quantity as the water a cured specimen drinks, per gram reacted.
shrink_mass = only(q[]) * M_H2O / (α_cured * m_cement)

@printf(
    "chemical shrinkage, from molar volumes : %.4f cm3 per g of reacted cement\n",
    shrink_vol
)
@printf(
    "water drawn in by the cured specimen   : %.4f g   per g of reacted cement\n",
    shrink_mass
)
# Powers' two coefficients, read back out of the function that applies them --
# below either one the cap is w/c divided by it, so the page cannot quote a
# number the code does not use.
w_probe = 0.25
k_sealed = w_probe / powers_alpha_max(w_probe)
k_saturated = w_probe / powers_alpha_max(w_probe; curing = :saturated)
@printf(
    "Powers, as the gap between his two coefficients (%.2f - %.2f): %.4f g/g\n",
    k_sealed, k_saturated, k_sealed - k_saturated
)

@printf("%8s %12s %14s\n", "alpha", "certified", "cm3 per g reacted")
for α in (0.6, α_sealed, 1.0)
    e, cert = equilibrate_certified(arrested(wc_cure, α))
    dv = ustrip(uconvert(us"cm^3", volume(fresh_c).total)) -
        ustrip(uconvert(us"cm^3", volume(with_unreacted(e, α)).total))
    @printf("%8.3f %12s %14.4f\n", α, cert.optimal, dv / (α * m_cement))
end
