# [CEM I 52.5 N and slag in an isothermal calorimeter, read off the states](@id sec-example-isothermal)

!!! info "Before this page"
    [CEM III/A, a blastfurnace cement](@ref ex-cem3-slag), for a slag entered through its
    oxides.

An isothermal calorimeter holds the sample at one temperature and records the
heat it gives off. Enthalpy is a state function, so that heat is what the
enthalpy of the sample loses between two states at that temperature,

```math
Q(t) \;=\; H(t_0) - H(t), \qquad H = \sum_i n_i\,\Delta_f H_i(T) ,
```

with the reactants, the ions and the hydrates each counted once and no reaction
to write down (Eqs. 17–21 of [Lavergne2018](@citet)). It is an output: once the
states are computed, the heat is read off them. Along a kinetic trajectory that
is [`heat_release`](@ref), compared with measured curves on
[the calibration page](@ref ex-hydration-calibration); this page reads it at
measured degrees of hydration instead, on the pastes of [Gruyaert2010](@citet).

## The pastes

A CEM I 52.5 N and a blast-furnace slag, pastes at w/b = 0.5 with 0 to 85 %
slag. [Gruyaert2010](@citet) report the oxide composition of both, the degree of
hydration of the cement and of the slag by image analysis at 2 days and 28
months, and the heat by isothermal calorimetry at 20 °C, all in
`data/literature/Gruyaert2010.json`. `scripts/gruyaert2010.jl` builds a paste
from them; its assumptions are written there.

```@example isocal
using ChemistryLab, DynamicQuantities, Printf, Plots
include(joinpath(pkgdir(ChemistryLab), "scripts", "gruyaert2010.jl"))

c = gruyaert_bogue("OPC-CAL")
for (ph, f) in c.clinker
    @printf("  %-5s %5.1f %%\n", ph, 100f)
end
@printf("  gypsum %4.1f %%, calcite %4.1f %%\n", 100c.gypsum, 100c.calcite)
```

## The cement alone

The heat of the plain paste at the degree of hydration the image analysis gives,
against the heat the calorimeter measured at the same age. At 2 days that is the
reaction degree of [Gruyaert2010; Table 6](@cite) times the total heat of its
Table 2; at 28 months the cement is at 74 %, the ultimate degree the article
computes for w/c = 0.5, and the total heat, extrapolated to infinite time, is
the comparison.

```@example isocal
cs = gruyaert_system()
J(st) = ustrip(us"J", enthalpy(st))
heat(p) = (J(p.initial) - J(p.eq)) / 100                 # J per g of binder
G10(table, column; where...) = only(getproperty(literature_table("Gruyaert2010", table; where...), column))
α(days, sb) = G10("hydration_degree", :alpha_cement; age_days = days * u"d", slag_to_binder = sb) / 100
Q_total(sb) = ustrip(us"J/g", G10("isothermal_total_heat", :heat_infinity; slag_to_binder = sb))
r2 = G10("reaction_degree", :r; age_days = 2u"d", slag_to_binder = 0.0) / 100

p2 = gruyaert_paste(cs; slag = 0.0, alpha_cement = α(2, 0.0), alpha_slag = 0.0)
p28 = gruyaert_paste(cs; slag = 0.0, alpha_cement = α(852, 0.0), alpha_slag = 0.0)
@printf("certified: %s, %s\n", p2.certificate.optimal, p28.certificate.optimal)
@printf("2 days:    α = %.2f, computed %5.0f J/g, measured %5.0f J/g\n",
        α(2, 0.0), heat(p2), r2 * Q_total(0.0))
@printf("28 months: α = %.2f, computed %5.0f J/g, measured total %5.0f J/g\n",
        α(852, 0.0), heat(p28), Q_total(0.0))
@printf("per unit degree of hydration: computed %.0f and %.0f J/g, measured %.0f and %.0f J/g\n",
        heat(p2) / α(2, 0.0), heat(p28) / α(852, 0.0), r2 * Q_total(0.0) / α(2, 0.0),
        Q_total(0.0) / α(852, 0.0))
```

```@example isocal
fig = bar([0.8, 1.8], [heat(p2), heat(p28)]; bar_width = 0.4, label = "computed from the states", color = :steelblue,
          xticks = ([1, 2], ["2 days", "28 months"]), ylabel = "heat (J per g of binder)", legend = :topleft,
          size = (700, 380), left_margin = 6Plots.mm, bottom_margin = 5Plots.mm,
          title = "CEM I 52.5 N paste, w/b 0.5, 20 °C")
bar!(fig, [1.2, 2.2], [r2 * Q_total(0.0), Q_total(0.0)]; bar_width = 0.4, color = :black,
     label = "measured (at 28 months: the total extrapolated)")
plot!(fig; ylims = (0, 1.15 * ylims(fig)[2]))
```

What the minimization made of the reacted cement and its water, at 28 months,
in grams per 100 g of binder:

```@example isocal
mass_of(st) = Dict(symbol(st.system.species[i]) => ustrip(us"g", st.n[i] * st.system.species[i][:M])
                   for i in eachindex(st.n) if !(i in st.system.idx_aqueous) || symbol(st.system.species[i]) == "H2O@")
before, after = mass_of(p28.initial), mass_of(p28.eq)
for d in (before, after)
    d["water"] = pop!(d, "H2O@")
end
keep = sort([k for k in union(keys(before), keys(after)) if max(get(before, k, 0.0), get(after, k, 0.0)) > 0.5];
            by = k -> -max(get(before, k, 0.0), get(after, k, 0.0)))
xs = collect(eachindex(keep))
fig = bar(xs .- 0.2, [get(before, k, 0.0) for k in keep]; bar_width = 0.4, label = "reacted cement and water",
          color = :grey60, xticks = (xs, keep), xrotation = 30, ylabel = "g per 100 g of binder", legend = :topright,
          size = (860, 400), left_margin = 6Plots.mm, bottom_margin = 12Plots.mm, title = "The plain paste, α = 0.74")
bar!(fig, xs .+ 0.2, [get(after, k, 0.0) for k in keep]; bar_width = 0.4, label = "at equilibrium", color = :steelblue)
plot!(fig; ylims = (0, 1.1 * ylims(fig)[2]))
```

At 2 days the computed heat is 7 % above the measured one; at 28 months it is
18 % below the total extrapolated from the measurements. The calculation carries
one heat per unit degree of hydration, fixed by the enthalpies of the phases,
and it comes out at 488 and 480 J/g at the two ages. The measurements carry 457
at 2 days and 585 at infinite time, and the two cannot hold together with the
degrees of hydration of the image analysis: [Gruyaert2010](@citet) note it
themselves, their reaction degree at 2 days being low against the degree of
hydration once the ultimate degree of 74 % is accounted for. The calculation,
which has no parameter to adjust, falls between them.

## The slag, whose enthalpy no database holds

A glass has no formula, and none of the databases carries its enthalpy of
formation: the slag enters the budget through its oxides, and its heat cannot be
computed. What the measurements can give is the enthalpy the glass must have for
the computed heat to be the measured one. At 28 months, with the cement and the
slag at their measured degrees and the total heat of
[Gruyaert2010; Table 2](@cite),

```math
h_{\rm glass} \;=\; \frac{Q\,m_b + H_{\rm eq} - H_0}{m_{\rm s}\,\alpha_{\rm s}} ,
```

per gram of reacted slag, with `H₀` the enthalpy of the reacted cement and the
water and `H_eq` that of the computed paste. Two blends give it twice, and one
glass should give one number.

```@example isocal
αs(sb) = G10("hydration_degree_slag", :alpha_slag; age_days = 852u"d", slag_to_binder = sb) / 100
h_glass = map((0.5, 0.85)) do sb
    p = gruyaert_paste(cs; slag = sb, alpha_cement = α(852, sb), alpha_slag = αs(sb))
    m_c, m_s = 100 * (1 - sb) * α(852, sb), 100 * sb * αs(sb)
    h = (Q_total(sb) * 100 + J(p.eq) - J(p.initial)) / m_s
    @printf("slag %.2f: reacted cement %4.1f g, slag %4.1f g, certified %s, h_glass = %6.0f J/g\n",
            sb, m_c, m_s, p.certificate.optimal, h)
    h
end
```

The two blends give the glass to within 0.5 %, although the reacted cement
weighs 1.3 times the reacted slag in the first and 0.41 times in the second: an
error in the heat computed for the cement would move the two apart, in that
proportion. The enthalpy is derived, and it carries the uncertainty of the
extrapolated totals and of the degrees of hydration; fed back into a calculation
of these pastes, it would reproduce their heat by construction, and says
nothing then.

The number means more against the crystalline oxides the glass is made of,
taken from SLOP98 at the same temperature, the sulfate as anhydrite less its
lime, and against two crystalline calcium silicates of the database, gehlenite
and åkermanite, whose half-and-half mixture has nearly the slag's composition.

```@example isocal
slop = Dict(symbol(s) => s for s in
            build_species(datapath("slop98-inorganic-thermofun.json"); verbose = false))
T0, P0 = temperature(p28.eq), pressure(p28.eq)
h(s) = ustrip(us"J/mol", slop[s][:ΔₐH⁰](T = T0, P = P0; unit = true))
M(s) = ustrip(us"g/mol", slop[s][:M])
oxide = Dict("CaO" => "Lim", "SiO2" => "Qtz", "Al2O3" => "Crn", "Fe2O3" => "Hem", "MgO" => "Per")
per_gram = Dict(o => h(s) / M(s) for (o, s) in oxide)
per_gram["SO3"] = (h("Anh") - h("Lim")) / (M("Anh") - M("Lim"))
bfs = gruyaert_oxides("BFS-CAL")
h_oxides = sum(ustrip(getproperty(bfs, Symbol(o))) / 100 * v for (o, v) in per_gram)
@printf("glass less its oxides: %4.0f and %4.0f J/g\n", (h_glass .- h_oxides)...)

# A crystal less its oxides, per gram, and its oxide mass fractions, from its formula.
function from_oxides(s)
    a = atoms(slop[s])
    n = Dict("CaO" => a[:Ca], "SiO2" => a[:Si], "Al2O3" => get(a, :Al, 0) / 2, "MgO" => get(a, :Mg, 0))
    dh = (h(s) - sum(Float64(k) * h(oxide[o]) for (o, k) in n)) / M(s)
    return dh, Dict(o => Float64(k) * M(oxide[o]) / M(s) for (o, k) in n)
end
(dh_gh, w_gh), (dh_ak, w_ak) = from_oxides("Gh"), from_oxides("Ak")
@printf("gehlenite less its oxides %4.0f J/g, åkermanite %4.0f J/g\n", dh_gh, dh_ak)
for o in ("CaO", "SiO2", "Al2O3", "MgO")
    @printf("  %-6s half and half %5.1f %%, slag %5.1f %%\n", o,
            50 * (w_gh[o] + w_ak[o]), ustrip(getproperty(bfs, Symbol(o))))
end
```

```@example isocal
p_h = bar(["glass, 50 % slag", "glass, 85 % slag", "gehlenite", "åkermanite"],
    [(h_glass .- h_oxides)..., dh_gh, dh_ak]; legend = false, color = [:steelblue, :steelblue, :grey60, :grey60],
    ylabel = "enthalpy less the oxides (J/g)", size = (720, 380), left_margin = 6Plots.mm, bottom_margin = 5Plots.mm,
    title = "Each less the crystalline oxides of its composition, SLOP98")
plot!(p_h; ylims = (1.1 * ylims(p_h)[1], 0))
```

The glass lies 330 to 400 J/g below the crystalline oxides of its composition,
and above the two crystals, which lie 500 and 670 J/g below theirs. A glass
holds more enthalpy than the crystals of its composition, and that is the order
found. The half-and-half mixture is richer in alumina than the slag, so the
difference between them, some 200 J/g, gives an order of magnitude and not an
enthalpy of vitrification.

## See also

  - [A CEM I 52.5 N mortar in a semi-adiabatic calorimeter, inside the kinetics](@ref ex-semiadiabatic),
    where the heat raises the temperature and the temperature the rates.
  - [Calibrating the hydration kinetics of a CEM I on measured calorimetry](@ref ex-hydration-calibration),
    the heat along a trajectory against measured curves.
  - [Bound water of CEM I 52.5 N and slag pastes, and the thermogram it integrates to](@ref sec-example-tga),
    the other measurement read off the same pastes.
