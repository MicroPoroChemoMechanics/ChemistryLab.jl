# [The cements of Lavergne et al. (2018): degrees of hydration, phase volumes and heat](@id ex-lavergne2018)

!!! info "Before this page"
    [A complete CEM I 52.5 N, through its pore solution](@ref ex-ionic-opc),
    whose model this page runs on the cement of the article.

[Lavergne2018](@citet) estimate the mechanical properties of hydrating cement
pastes from a hydration model: the kinetic law of [ParrottKilloh1984](@citet)
for each clinker phase, a set of stoichiometric reactions for the hydrates, and
an energy balance for the temperature. Before the mechanics, they check that
model against measurements from the literature and from their own laboratory:
degrees of hydration of the clinker phases by X-ray diffraction, heat in
isothermal and semi-adiabatic calorimeters. This page takes those comparisons
one by one and puts the package's calculation beside each: the same kinetic law
where the article uses it, and, for what the hydrates are, a Gibbs minimization
in place of the stoichiometric reactions.

Every number the page compares with is read from the article: the compositions of
its Table 7, the curves of its figures read exactly from their vector drawings,
measured points and model estimates alike (`data/literature/Lavergne2018.json`).

!!! tip "Questions this page answers"
      - How close does the package's Parrott–Killoh law come to the article's
        estimates, and to the measured degrees of hydration?
      - What phases does a Gibbs minimization give a CEM I paste as it hydrates,
        against the stoichiometric reactions of the article?
      - How does the same paste heat a semi-adiabatic calorimeter?

## The cements

The article gathers in its Table 7 the cements of every test it reproduces: their
Bogue composition, their calcium sulfates and their Blaine fineness. The first
seven serve the degrees of hydration below.

```@example lavergne
using ChemistryLab, DynamicQuantities, OptimaSolver, OrdinaryDiffEq, Printf, Plots, Logging
gr()

t7 = literature_table("Lavergne2018", "table7_cements")
row(c) = findfirst(==(c), t7.name)
println("cement   C3S   C2S   C3A  C4AF  sulfates        Blaine (m²/kg)")
for c in ("A", "B", "C", "OPCN", "OPCS", "OPC", "CRC")
    i = row(c)
    pct(x) = ismissing(x) ? "    –" : @sprintf("%5.1f", x)
    @printf("%-6s %s %s %s %s  %-14s  %.1f\n", c, pct(t7.C3S[i]), pct(t7.C2S[i]),
            pct(t7.C3A[i]), pct(t7.C4AF[i]), t7.calcium_sulfates[i], ustrip(us"m^2/kg", t7.blaine[i]))
end
```

## Degrees of hydration of the clinker phases

Each clinker phase hydrates by the law of [ParrottKilloh1984](@citet), the
slowest of three rates: nucleation and growth, diffusion through the hydrates,
and the formation of a shell around the grain, with the parameters of the
article's Table 3 ([`PK84_PARAMS_C3S`](@ref) and its siblings). The rate scales
with the Blaine fineness of the cement and, through the activation energy of each
phase (Table 4), with the temperature ([`parrott_killoh_avrami`](@ref)). The
degree of a phase at any time follows by integrating that rate alone, as the
article does: no other phase and no hydrate enters it.

```@example lavergne
const PK = Dict("C3S" => PK84_PARAMS_C3S, "C2S" => PK84_PARAMS_C2S,
                "C3A" => PK84_PARAMS_C3A, "C4AF" => PK84_PARAMS_C4AF)

# The degree of hydration of one clinker phase of a cement of Blaine fineness B,
# at the times `days`, under the temperature history T(t) (t in seconds, T in K).
function pk_degree(phase, B, T, days)
    law = parrott_killoh_avrami(PK[phase], phase; blaine = B)
    idx = Dict(phase => 1)
    n0, lna = StateView([1.0], idx), StateView([0.0], idx)
    f(n, _, t) = [-law(T(t), 1.0e5, t, StateView(n, idx), lna, n0)]
    sol = solve(ODEProblem(f, [1.0], (0.0, 86400 * maximum(days))), Rodas5P();
                abstol = 1.0e-12, reltol = 1.0e-9)
    return [1 - sol(86400 * d)[1] for d in days]
end
nothing # hide
```

[Lavergne2018; Fig. 4](@citet) compares that law with the degrees measured by
XRD/Rietveld analysis on three cements, A, B and C, at 20 °C. Below, the
package's integration (solid lines), the article's own estimate (dashed) and the
measurements (markers), on the four clinker phases:

```@example lavergne
PHASES = ("C3S", "C3A", "C2S", "C4AF")
days = 10 .^ range(-1, 3; length = 120)
at20(t) = 293.15

function degree_panels(meas, est, groups; T = g -> at20, label = g -> g.name, title = "")
    panels = map(PHASES) do ph
        p = plot(; xscale = :log10, xlims = (0.1, 1000), ylims = (0, 1), legend = false,
                 title = ph, xlabel = "t [days]", ylabel = "α")
        for (k, g) in enumerate(groups)
            c = palette(:tab10)[k]
            sm(t) = [i for i in eachindex(t.phase) if t.phase[i] == ph && g.select(t, i)]
            ie, im = sm(est), sm(meas)
            isempty(ie) && isempty(im) && continue
            plot!(p, days, pk_degree(ph, g.blaine, T(g), days); color = c, lw = 2, label = label(g))
            plot!(p, ustrip.(u"d", est.time[ie]), est.alpha[ie]; color = c, ls = :dash, label = "")
            scatter!(p, ustrip.(u"d", meas.time[im]), meas.alpha[im]; color = c, ms = 3, label = "")
        end
        p
    end
    plot(panels...; layout = (2, 2), size = (860, 640), legend = :topleft,
         plot_title = title, left_margin = 4Plots.mm, bottom_margin = 4Plots.mm)
end

fig4m = literature_table("Lavergne2018", "fig4_hydration_degree_measured")
fig4e = literature_table("Lavergne2018", "fig4_hydration_degree_estimated")
cements4 = [(name = c, blaine = t7.blaine[row(c)], select = (t, i) -> t.cement[i] == c) for c in ("A", "B", "C")]
degree_panels(fig4m, fig4e, cements4; label = g -> "cement " * g.name,
              title = "Fig. 4: computed (solid), the article's estimate (dashed), measured")
```

The two integrations of the same law agree to a few hundredths. The
measurements are another matter, as the article says itself: the law has the
right characteristic times for C₃S, C₂S and C₃A, but its degree can be off by
0.2 or more, and its C₄AF is uncertain on both sides, model and measurement.

```@example lavergne
# The largest difference with the article's estimate, phase by phase.
function gap(meas, est, groups; T = g -> at20)
    out = Dict{String, Float64}()
    for g in groups, ph in PHASES
        ie = [i for i in eachindex(est.phase) if est.phase[i] == ph && g.select(est, i)]
        isempty(ie) && continue
        d = maximum(abs.(pk_degree(ph, g.blaine, T(g), ustrip.(u"d", est.time[ie])) .- est.alpha[ie]))
        out[ph] = max(get(out, ph, 0.0), d)
    end
    return out
end
for (ph, d) in sort(collect(gap(fig4m, fig4e, cements4)); by = first)
    @printf("%-5s largest |Δα| with the article's estimate: %.3f\n", ph, d)
end
```

### Away from 20 °C

[Lavergne2018; Fig. 5](@citet) take the measurements of two cements, OPCN and
OPCS, cured at 20 °C for a day and then stored at 10, 20, 30 or 40 °C. The same
law follows the temperature through the activation energy of each phase:

```@example lavergne
fig5m = literature_table("Lavergne2018", "fig5_hydration_degree_measured")
fig5e = literature_table("Lavergne2018", "fig5_hydration_degree_estimated")
history(Tc) = t -> t < 86400 ? 293.15 : Tc + 273.15
cured(c) = [(name = "$(Int(Tc)) °C", Tc = Tc, blaine = t7.blaine[row(c)],
             select = (t, i) -> t.cement[i] == c && t.temperature_C[i] == Tc) for Tc in (10.0, 20.0, 30.0, 40.0)]
degree_panels(fig5m, fig5e, cured("OPCN"); T = g -> history(g.Tc), label = g -> g.name,
              title = "Fig. 5, OPCN: computed (solid), the article's estimate (dashed), measured")
```

```@example lavergne
degree_panels(fig5m, fig5e, cured("OPCS"); T = g -> history(g.Tc), label = g -> g.name,
              title = "Fig. 5, OPCS: computed (solid), the article's estimate (dashed), measured")
```

The temperature moves the computed degrees as it moves the measured ones, which
is what the activation energies of Table 4 are for; the article's estimate and
this integration stay within a few hundredths of each other again.

### Two more cements

[Lavergne2018; Fig. 6](@citet) end with an ordinary cement and a cement of
another composition, OPC and CRC:

```@example lavergne
fig6m = literature_table("Lavergne2018", "fig6_hydration_degree_measured")
fig6e = literature_table("Lavergne2018", "fig6_hydration_degree_estimated")
cements6 = [(name = c, blaine = t7.blaine[row(c)], select = (t, i) -> t.cement[i] == c) for c in ("OPC", "CRC")]
degree_panels(fig6m, fig6e, cements6; title = "Fig. 6: computed (solid), the article's estimate (dashed), measured")
```

```@example lavergne
for (ph, d) in sort(collect(gap(fig6m, fig6e, cements6)); by = first)
    @printf("%-5s largest |Δα| with the article's estimate: %.3f\n", ph, d)
end
```

Here the two calculations part by up to a quarter, and the figure shows where:
the article's estimates level off and stop rising after about a hundred days,
which is what its model does when the internal relative humidity of the paste
falls below 80 % (its Eq. 10). That term, which depends on the water to cement
ratio of a test the article does not give, is left out of the integration above;
before the article's curves level off, the two agree to a few hundredths.

## The phases of a hydrating paste

The kinetic law says how much of each clinker phase has reacted, not what it has
become. The article answers with stoichiometric reactions written in advance
(its Table 2); the package answers with a Gibbs minimization of everything that
is not clinker, at every step of the integration
([A complete CEM I 52.5 N, through its pore solution](@ref ex-ionic-opc)). The
paste of [Lavergne2018; Fig. 1](@citet), left, is that of their own cement, the
CEM I 52.5 N of Table 9, at w/c = 0.5 and sealed; the package's model of that
cement runs here for a year.

```@example lavergne
include(joinpath(pkgdir(ChemistryLab), "scripts", "ionic_hydration.jl"))
quiet(f) = with_logger(f, NullLogger())
paste = quiet() do
    run_ionic_hydration(; wb = 0.5, tend = 365 * 86400.0)
end
tdays = [0.1, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 5, 7, 10, 14, 28, 56, 90, 180, 365]
states = quiet() do
    speciated_states(paste.sol, paste.kp; times = tdays .* 86400)
end
_, fractions, _, _ = ionic_phase_history(paste, tdays .* 86400; states)

# The degree of hydration of the cement: the mass of clinker consumed.
clinker = ("C3S", "C2S", "C3A", "C4AF")
M = Dict(c => ustrip(us"kg/mol", paste.cs[c][:M]) for c in clinker)
mass(st) = sum(ustrip(us"mol", moles(st, c)) * M[c] for c in clinker)
α = [1 - mass(st) / mass(paste.state0) for st in states]
@printf("%s; α = %.2f after 28 days, %.2f after a year\n", paste.sol.retcode, α[13], α[end])
```

The article's phases and the package's families are put side by side: its
"Aft, Afm" are the ettringite and the AFm phases, its C-S-H gel holds its gel
water, as the package's does, and its "air" is the empty porosity a sealed paste
gains as it shrinks. The package's hydrogarnet and iron hydroxide, which the
article does not draw, are left out of the comparison.

```@example lavergne
fig1 = literature_table("Lavergne2018", "fig1_volume_fractions_estimated")
families = [
    "unhydrated cement" => ["anhydrous"], "gypsum" => ["gypsum"], "calcite" => ["calcite"],
    "Aft, Afm" => ["AFt", "AFm"], "Portlandite" => ["CH"], "C-S-H gel" => ["C-S-H"],
    "capillary water" => ["water"],
]
push!(families, "air" => ["void"])
mine(f, keys) = sum(get(f, k, 0.0) for k in keys)
theirs(name) = (i = [j for j in eachindex(fig1.phase) if fig1.paste[j] == "plain paste, w/c = 0.5" && fig1.phase[j] == name];
                (fig1.alpha[i], fig1.volume_fraction[i]))
panels = map(families) do (name, keys)
    a, v = theirs(name)
    p = plot(a, v; color = :black, ls = :dash, label = "article",
             title = name, xlabel = "α", ylabel = "volume fraction", xlims = (0, 1), legend = false)
    plot!(p, α, [mine(f, keys) for f in fractions]; color = :steelblue, lw = 2, marker = :circle, ms = 3, label = "computed")
end
plot(panels...; layout = (2, 4), size = (1000, 520), plot_title = "Fig. 1, left: computed (blue), the article's model (dashed)",
     left_margin = 3Plots.mm, bottom_margin = 4Plots.mm)
```

At the degree the paste reaches in 28 days, the two compositions read, in
fractions of the volume of the fresh paste (the article's curve taken between its
two nearest vertices):

```@example lavergne
k28 = findfirst(==(28), tdays)
interp(x, xs, ys) = (j = searchsortedlast(xs, x); j == length(xs) ? ys[end] :
                     ys[j] + (ys[j + 1] - ys[j]) * (x - xs[j]) / (xs[j + 1] - xs[j]))
@printf("α = %.2f\n%-18s %9s %9s\n", α[k28], "", "computed", "article")
for (name, keys) in families
    a, v = theirs(name)
    @printf("%-18s %9.3f %9.3f\n", name, mine(fractions[k28], keys), interp(α[k28], a, v))
end
```

The anhydrous cement, the portlandite and the capillary water agree within a
hundredth of the volume. The minimization puts more of it in the C-S-H gel and
less in the AFt and AFm phases and in the empty porosity. The early sulfate and
carbonate phases tell the two approaches apart most: the minimization consumes
the calcite at once, into carboaluminates, and the gypsum from the first hours,
where the reactions of the article form ettringite while gypsum lasts and then
monocarboaluminate while calcite does, so that its calcite lasts to the end.

## [A CEM I 52.5 N mortar in a semi-adiabatic calorimeter, inside the kinetics](@id ex-semiadiabatic)

An isothermal calorimeter holds the sample at one temperature, and what it
records can be read off a calculated trajectory afterwards, as the enthalpy the
states lose ([CEM I 52.5 N and slag in an isothermal calorimeter, read off the states](@ref sec-example-isothermal)). A
semi-adiabatic calorimeter cannot be read that way. It lets the heat of
hydration raise the temperature of the sample against the losses of the vessel,
and the temperature raises the rates of the reactions in turn. The temperature
is then an unknown of the kinetic problem: the enthalpy of the cell, the paste
over its whole composition and the vessel, changes only by what leaves through
the walls, and the temperature is the root of that balance, solved with the
equilibrium of the paste at every evaluation
([The semi-adiabatic cell](@ref sec-theory-pe-calorimeters)). The losses of the
vessel are quadratic in the temperature difference,
``\mathcal{L}(\Delta T) = a\,\Delta T + b\,\Delta T^2``.

This section reproduces the test of [Lavergne2018](@citet) on the plain-cement
mortar `C100` at w/b = 0.5, in their calorimeter, with the model of the paste
above and nothing adjusted.

### The cell

The cell is NF EN 196-9's, as [Lavergne2018](@citet) calibrated it:
their Eq. (23) for the losses, their Table 11 for the mix. The sand keeps the
temperature moderate and takes no part in the chemistry; it enters with the
vessel and the water it absorbs as a fixed heat capacity, while the paste's own
heat capacity is summed over its composition at every step.

```@example lavergne
mix = CALORIMETRY_MIX_C100
T_env = 293.15u"K"
cell = semiadiabatic_cell(; mix, T0 = T_env, T_env)
@printf("mortar: %.0f g binder, %.0f g sand, %.0f g water, of which %.1f g in the sand (w/b %.3f)\n",
        ustrip(us"g", mix.binder), ustrip(us"g", mix.sand), ustrip(us"g", mix.water),
        ustrip(us"g", mix.absorbed), mix.wb)
@printf("fixed heat capacity: vessel %.0f + sand %.0f + absorbed water %.0f = %.0f J/K\n",
        CALORIMETRY_VESSEL_CP, sand_heat_capacity(mix.sand), water_heat_capacity(mix.absorbed),
        ustrip(us"J/K", cell.Cp))
@printf("losses: a = %.4f W/K, b = %.2e W/K²\n", CALORIMETRY_LOSS_A, CALORIMETRY_LOSS_B)
```

The heat capacity of the vessel is taken as 380 J/K; [Lavergne2018](@citet) give
"about 380 kJ/K", and the measured temperature rises below correspond to 380 J/K
(with 380 kJ/K they would stay below one kelvin). The note on
`CALORIMETRY_VESSEL_CP` in `scripts/ionic_hydration.jl` gives the arithmetic.
The cell starts at 20 °C, where the measured curve starts and where the kinetic
parameters are referred.

### The run

The mortar's binder, 371 g, is integrated with the cell in the state: the
dissolution of the four clinker phases, the heat lost through the walls, and,
inside every evaluation, the Gibbs minimization of everything else and the
temperature it is in balance with.

```@example lavergne
binder = mix.binder
# The solver's warnings are summarized by what they are about, printed below:
# the re-speciations that failed, and the worst element balance of the accepted
# steps, in moles.
semi = quiet() do
    run_ionic_hydration(; wb = mix.wb, binder_mass = binder, calorimeter = cell, tend = 5 * 86400.0)
end
balance(run) = (run.sol.prob.p.eq_failures[], run.sol.prob.p.eq_worst_abs_acc[])
t, T = temperature_profile(semi.sol, cell)
j = argmax(T)
@printf("%s, %d accepted steps, %d failed re-speciations, worst balance %.1e mol\n",
        semi.sol.retcode, length(t), balance(semi)...)
@printf("maximum %.1f °C at %.2f d\n", T[j] - 273.15, t[j] / 86400)
```

### Against the measurement

The measured temperature is read from [Lavergne2018; Fig. 15(a)](@citet), from
its maximum to 3.5 days; before, the figure draws it as crosses that overlap
those of the other mixes, and it is not transcribed. The figure is a raster image
in the article, read as the notes of `data/literature/Lavergne2018.json`
describe.

```@example lavergne
meas = literature_table("Lavergne2018", "semi_adiabatic_C100_wb050_temperature")
tm = ustrip.(u"d", meas.time)
Tm = meas.temperature_C
Tc = last(temperature_profile(semi.sol, cell; times = tm .* 86400)) .- 273.15
k = argmax(Tm)
@printf("maximum: measured %.1f °C at %.2f d, computed %.1f °C at %.2f d\n",
        Tm[k], tm[k], T[j] - 273.15, t[j] / 86400)
println("  t (d)   measured   computed")
for i in 1:6:length(tm)
    @printf("  %5.2f   %7.1f    %7.1f\n", tm[i], Tm[i], Tc[i])
end
```

### What the feedback does

The same paste at 20 °C in an isothermal calorimeter gives the heat rate a
calculation without feedback would feed the cell. Integrated afterwards, as
`langavant_temperature` does, it gives the temperature the cell would reach if
the reactions ignored it.

```@example lavergne
iso_cal = IsothermalCalorimeter(T_env)
iso = quiet() do
    run_ionic_hydration(; wb = mix.wb, binder_mass = binder, calorimeter = iso_cal, tend = 5 * 86400.0)
end
ti, q = heat_flow(iso.sol, iso_cal)                          # W, for 371 g of binder
# The paste's heat capacity at the start, held, per kilogram of binder as
# `langavant_temperature` takes it: the cell's fixed part is most of the total,
# and the paste's own changes by a few percent as it hydrates.
m_g = ustrip(us"g", binder)
fresh = ChemicalState(semi.cs, semi.state0.n ./ (m_g / 1000); T = T_env)
T_off = langavant_temperature(ti, q ./ m_g, fill(fresh, length(ti)); mix)
i = argmax(T_off)
@printf("without feedback: maximum %.1f °C at %.2f d (%d failed re-speciations, worst balance %.1e mol)\n",
        T_off[i] - 273.15, ti[i] / 86400, balance(iso)...)

plot(t ./ 86400, T .- 273.15; lw = 2, label = "computed, coupled",
     xlabel = "time [days]", ylabel = "T [°C]", size = (720, 400), legend = :topright)
plot!(ti ./ 86400, T_off .- 273.15; lw = 2, ls = :dash, label = "computed, without feedback")
scatter!(tm, Tm; ms = 3, label = "measured (Lavergne et al. 2018)")
```

The coupled calculation reaches 56.4 °C at 0.81 day, where the measurement
peaks at 52.1 °C at 0.75 day, and it stays 2.8 to 5.4 K above the measured
curve through the cooling that follows. Nothing has been adjusted on this curve: the
kinetic parameters are those of the model of the paste, the cell and its losses
those [Lavergne2018](@citet) calibrated. The run without feedback, the heat flow
of the paste held at 20 °C integrated afterwards through the same cell, peaks at
40.1 °C only, and later, at 1.02 day. The difference, 16 K out of a rise of
36 K, is the acceleration of the reactions by the temperature they raise,
through their activation energies; a semi-adiabatic test is therefore a test of
those energies as much as of the heat, and it cannot be read off an isothermal
calculation.

## What this page does not reproduce yet

The article compares its model with other tests, which this page does not take
up yet: the isothermal calorimetry of the 27 cements of [Lavergne2018; Fig. 7](@citet),
whose compositions its Table 7 gives (`table7_cements`), the bound water and the
portlandite of blended pastes (Fig. 8), the adiabatic calorimetry of concretes
with silica fume and fly ash (Fig. 9), the paste with silica fume of Fig. 1, and
the semi-adiabatic tests of the other mixes of Table 11 (Figs. 15 and 16).

## See also

  - [A complete CEM I 52.5 N, through its pore solution](@ref ex-ionic-opc), the
    model of the paste and its assumptions.
  - [Hydration kinetics of a CEM I 52.5 R clinker](@ref sec-clinker-kinetics), the
    same law on the stoichiometric formulation.
  - [CEM I 52.5 N and slag in an isothermal calorimeter, read off the states](@ref sec-example-isothermal)
    and [Bound water of CEM I 52.5 N and slag pastes, and the thermogram it integrates to](@ref sec-example-tga),
    the measurements that are outputs of a calculation rather than part of it.
