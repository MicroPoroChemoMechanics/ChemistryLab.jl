# [Glasses of slag, fly ash and silica fume dissolving at pH 13](@id ex-glass-dissolution)

!!! info "Before this page"
    [Rate laws, and every parameter in them](@ref sec-theory-kinetics), for the
    transition-state rate laws and their inhibitors ([`RateModelInhibitor`](@ref)).

What reacts in a slag or a fly ash is a glass, an amorphous calcium
aluminosilicate: the slag of the [ternary pastes at three temperatures](@ref ex-slag-temperature-pastes)
holds 35 % SiO₂, 42 % CaO and 11 % Al₂O₃ and is 95 % amorphous.
[Snellings2013](@citet) synthesized six such glasses, from the composition of a
slag (G1, G2) through those of fly ashes (G3, G4) and of a natural pozzolan (G5)
to silica (G6), dissolved them in NaOH solutions at pH 13 and 20 °C, so dilute
that nothing precipitated, and measured the initial rate of each per square
meter of its BET surface. Measured in the same way on every glass, these rates
compare the glasses with one another and show what the solution does to them.

The rate grows with the calcium of the glass, which breaks up its silicate
network. The paper relates the logarithm of the rate to the molar ratio of the
calcium to the network formers, a straight line drawn on its Fig. 8, which
[`snellings2013_glass`](@ref) returns as a rate constant:

```@example glass
using ChemistryLab, DynamicQuantities, Printf
include(joinpath(pkgdir(ChemistryLab), "scripts", "snellings2013_glass.jl"))

g = literature_table(SN13, "glasses")
t = literature_table(SN13, "initial_rates")
println("glass   CaO  Al2O3  SiO2     x   log10 k (law)   measured")
for i in eachindex(g.glass)
    ox = Dict("CaO" => ustrip(g.CaO_percent[i]) / 100, "Al2O3" => ustrip(g.Al2O3_percent[i]) / 100,
              "SiO2" => ustrip(g.SiO2_percent[i]) / 100)
    x = ChemistryLab._snellings2013_abscissa(ox, ("CaO",))
    logk = log10(snellings2013_glass(ox; Ea = 0.0).k(; T = 293.15))   # at 20 °C, Ea plays no part
    j = findfirst(k -> t.glass[k] == g.glass[i] && iszero(ustrip(t.Al_initial[k])) &&
                       iszero(ustrip(t.Ca_initial[k])) && iszero(ustrip(t.Si_initial[k])), eachindex(t.glass))
    @printf("%-5s %5.1f %6.1f %5.1f %6.3f %10.2f %13.2f\n", g.glass[i], ustrip(g.CaO_percent[i]),
            ustrip(g.Al2O3_percent[i]), ustrip(g.SiO2_percent[i]), x, logk, ustrip(t.log_rate[j]))
end
```

The rates are in mol of the glass's cations per m² and per second; the law
meets every measurement within 0.1, inside the 0.15 the paper gives as its
error. The ratio on the axis of the figure is labeled Ca/(Al + Si), but the
abscissas of its six points are ``n_\text{Ca}/(2\,n_\text{Al} + n_\text{Si})``
computed from the compositions, aluminum counted twice, and the printed line
holds with that ratio only (with the ratio of the label, a line through the same
rates would have a slope of 2.05, not 2.74): the law is used with the abscissa
of the figure. It is measured on glasses holding nothing but CaO, Al₂O₃ and
SiO₂, from silica to G1, at pH 13 and 20 °C; outside that range the function
refuses unless told to extrapolate, and it takes no activation energy for
granted.

## Calcium and aluminum in solution

Added to the solution, calcium slows every glass, and aluminum the glasses in
which calcium only balances aluminum (G3 to G6, the tectosilicate glasses);
silicon up to 11 mM moves the rates by 0.27 at most, which the paper finds not
significant. The paper gives the rates, not a law. A
factor ``(1 + K a)^{-1}``, ``a`` the activity of Ca²⁺ or of aluminate (AlO₂⁻ in
the database), is one where the ion is absent and falls as ``1/a`` where it is
plentiful; its ``K`` is fitted here, one for calcium on all glasses and one for
aluminum on the tectosilicate glasses, to the change of each rate against the
same glass in NaOH alone. The activities are those of the 51 solutions, each
computed at pH 13 with the sodium hydroxide that holds it (the paper does not
give its concentration):

```@example glass
rows = sn13_rows()
fit = sn13_fit(rows)
@printf("solutions certified: %d of %d; NaOH for pH 13: %.3f to %.3f mol/kg\n",
        count(r -> r.certified, rows), length(rows), extrema(r.NaOH for r in rows)...)
@printf("calcium,  all glasses:          K = %7.0f kg/mol, rms %.2f, worst %.2f (%d rates)\n",
        fit.ca.K, fit.ca.rms, fit.ca.worst, fit.ca.n)
@printf("aluminum, tectosilicate glasses: K = %7.0f kg/mol, rms %.2f, worst %.2f (%d rates)\n",
        fit.al.K, fit.al.rms, fit.al.worst, fit.al.n)
@printf("aluminum, G1 and G2 apart:       K = %7.0f kg/mol, rms %.2f, worst %.2f (%d rates)\n",
        fit.al_percalcic.K, fit.al_percalcic.rms, fit.al_percalcic.worst, fit.al_percalcic.n)
@printf("silicon, no factor:              rms %.2f, worst %.2f\n", fit.si_rms, fit.si_worst)
```

Calcium is described within the error of the measurements: a Ca²⁺ activity of
``1/K``, about ``10^{-4}``, halves the rate of any of these glasses. Aluminum is
not: one factor cannot follow both G3, slowed by 0.4 from the smallest addition
on and no further, and G6, slowed by 0.4, 0.6 and 1.1 as the aluminum grows.
On the slag-like glasses aluminum moves the rate by 0.26 to 0.40, which the paper
reads as largely within its error; the law leaves them without that factor.

```@example glass
using Plots
default(framestyle = :box, grid = false)
law(r) = r.log_base - log10(1 + fit.ca.K * r.a_Ca) -
         (r.glass in SN13_TECTOSILICATE ? log10(1 + fit.al.K * r.a_Al) : 0.0)
series = (("NaOH alone", r -> r.al == r.ca == r.si == 0, :circle), ("Al added", r -> r.al > 0, :utriangle),
          ("Ca added", r -> r.ca > 0, :square), ("Si added", r -> r.si > 0, :diamond))
fig = plot([-9.8, -6.4], [-9.8, -6.4]; color = :black, label = "", xlabel = "log₁₀ r measured (mol m⁻² s⁻¹)",
           ylabel = "log₁₀ r computed", size = (560, 480), legend = :topleft)
plot!(fig, [-9.8, -6.4], [-9.8, -6.4] .+ 0.15; color = :gray, ls = :dash, label = "± 0.15, the error of the paper")
plot!(fig, [-9.8, -6.4], [-9.8, -6.4] .- 0.15; color = :gray, ls = :dash, label = "")
for (lab, keep, m) in series
    s = filter(keep, rows)
    scatter!(fig, [r.log_rate for r in s], law.(s); marker = m, label = lab)
end
savefig(fig, "glass_dissolution.svg"); nothing # hide
```

![](glass_dissolution.svg)

These are initial rates far from equilibrium, at pH 13, 20 °C and calcium and
aluminum below 2.5 and 5 mM. A glass in a paste meets a pore solution at a pH
of 13 to 13.8, with calcium and aluminum set by the hydrates around it, for
months: whether this law, fitted on dilute solutions, follows the slag of a
paste is a question for the measurements of a paste, not one it answers.
