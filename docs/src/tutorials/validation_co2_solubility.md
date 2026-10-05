# [Carbon dioxide in water under pressure (Wiebe and Gaddy 1940)](@id sec-validation-co2-solubility)

!!! info "Before this page"
    [Real gases and pressure](@ref sec-theory-real-gases), for the equation of
    state of the gas and the pressure terms of the solution.

[WiebeGaddy1940](@citet) measured how much carbon dioxide dissolves in pure
water from 12 to 40 °C and from 25 to 500 atm. At these pressures the gas is far
from ideal: below 31 °C it condenses into a liquid above some 50 to 70 atm, and
above, a dense supercritical fluid. This page computes their 42 solubilities with
the gas ideal and with the equation of state of [PengRobinson1976](@citet), and
separates what each term of the calculation contributes. Nothing is fitted: the
critical constants of the gas are those of `phreeqc.dat`
[ParkhurstAppelo2013](@cite), the standard energies those of Cemdata18.

## 1. The measurements

Table I of the paper gives the volume of gas, reduced to 0 °C and 1 atm,
dissolved per gram of water. Divided by the molar volume of the ideal gas at
0 °C and 1 atm, it is the molality of the dissolved carbon dioxide
(`wg40_measured`); converted to mole fractions, it is the table
[Spycher2003](@citet) give in their Appendix A, to the last printed digit.

```@example co2
using ChemistryLab, DynamicQuantities, OptimaSolver, Printf
include(joinpath(pkgdir(ChemistryLab), "scripts", "wiebe_gaddy1940_co2.jl"))

d = wg40_measured()
temps = unique(d.T_C)
println("P (atm)  ", join((@sprintf("%8.2f °C", t) for t in temps)))
for P in sort(unique(round.(Int, d.P ./ 101325)))
    cells = map(temps) do t
        i = findfirst(k -> d.T_C[k] == t && round(Int, d.P[k] / 101325) == P, eachindex(d.T))
        i === nothing ? " "^11 : @sprintf("%11.3f", d.m[i])
    end
    println(rpad(P, 9), join(cells))
end
```

The molality, in mol/kg, rises steeply up to 50 to 75 atm and then slowly: past
that point the carbon dioxide is dense, and compressing it further raises its
fugacity much less than its pressure.

## 2. The calculation

A kilogram of water is put in contact with twenty moles of carbon dioxide, which
leaves a gas phase at every condition of the table, and brought to equilibrium
(`wg40_dissolved`). The aqueous species are the water, its ions, the dissolved
carbon dioxide and the bicarbonate and carbonate ions, with the extended
Debye–Hückel model; the bicarbonate is a twentieth of a percent of the dissolved
carbon. The gas is either the ideal gas of the database or the same species
given the critical constants of `phreeqc.dat` by [`peng_robinson`](@ref):

```@example co2
co2 = peng_robinson(WG40_DB["CO2"])
for (T_C, P_atm) in ((25.0, 50), (25.0, 100), (40.0, 100))
    st = ChemicalState(ChemicalSystem([co2], ["CO2"]); T = (T_C + 273.15) * u"K", P = P_atm * 101325.0u"Pa", n = [1.0u"mol"])
    @printf("%4.0f °C, %3d atm:  φ = %.3f,  Z = %.3f\n", T_C, P_atm, fugacity_coefficients(st)["CO2"], compressibility_factor(st))
end
```

At 25 °C the carbon dioxide is a vapor at 50 atm and a liquid at 100 atm, under
a quarter of the volume of the ideal gas; at 40 °C the fluid has no such jump. In
either case its fugacity is a fraction of its pressure, and that fraction is
what the dissolved amount follows.

```@example co2
r = wg40_compute()
dev_real = 100 .* (r.real ./ r.m .- 1)
dev_ideal = 100 .* (r.ideal ./ r.m .- 1)
@printf("Peng–Robinson: mean %+.1f %%, mean absolute %.1f %%, largest %.1f %%\n",
    sum(dev_real) / 42, sum(abs, dev_real) / 42, maximum(abs, dev_real))
@printf("ideal gas:     from %+.0f %% to %+.0f %%\n", minimum(dev_ideal), maximum(dev_ideal))
```

## 3. Measured and computed

```@example co2
using Plots
default(framestyle = :box, grid = false)
colors = palette(:viridis, length(temps) + 1)
cs_real, cs_ideal = wg40_system(:real), wg40_system(:ideal)
fig = plot(; xlabel = "P (atm)", ylabel = "dissolved CO₂ (mol/kg)", ylims = (0, 2.2), legend = :bottomright, size = (700, 430))
for (k, t) in enumerate(temps)
    T = t + 273.15
    Ps = 101325.0 .* (t < 15 ? (10:10:300) : (10:10:500))
    plot!(fig, Ps ./ 101325, [wg40_dissolved(cs_real, T, P) for P in Ps]; color = colors[k], label = @sprintf("%.4g °C", t))
    plot!(fig, Ps ./ 101325, [wg40_dissolved(cs_ideal, T, P) for P in Ps]; color = colors[k], linestyle = :dash, label = false)
    i = findall(==(t), r.T_C)
    scatter!(fig, r.P[i] ./ 101325, r.m[i]; color = colors[k], markersize = 4, label = false)
end
savefig(fig, "validation-co2-solubility.svg"); nothing # hide
```

![](validation-co2-solubility.svg)

The markers are the measurements, the solid lines the real gas, the dashed lines
the ideal gas, which leaves the frame above some 75 atm. The real gas follows
the shape of every isotherm: the sharp rise to the saturation of the gas at 12,
18 and 25 °C and the plateau of the liquid above it, the smoother rise of the
supercritical isotherms. The ideal gas, whose activity is its pressure, has
neither: it dissolves 12 % too much at 25 atm and 40 °C, and 5.6 times the
measurement at 300 atm and 12 °C.

## 4. What each term contributes

Under the pure gas and with unit activity coefficient, the molality of the
dissolved carbon dioxide is the product of four factors, and its logarithm the
sum of four terms (`wg40_terms`):

```math
\ln m = \underbrace{\frac{\mu^\circ_{\mathrm{gas}} - \mu^\circ_{\mathrm{aq}}}{RT}\Big|_{P^\circ}}_{\ln K_H}
      + \ln\frac{P}{P^\circ} + \ln\varphi(T, P)
      - \underbrace{\frac{\mu^\circ_{\mathrm{aq}}(T, P) - \mu^\circ_{\mathrm{aq}}(T, P^\circ)}{RT}}_{\text{volume}} ,
```

the Henry constant at the standard pressure, the pressure, the fugacity
coefficient of the gas, and the work of the partial molar volume of the
dissolved carbon dioxide, which its HKF equation of state gives as about
33 cm³/mol. The four add up to the molality computed at equilibrium to 0.1 %,
the share of the bicarbonate.

```@example co2
i40 = findall(==(40.0), r.T_C)
println("P (atm)   ln KH    ln P/P°   ln φ    volume   ln m: ideal  +φ      +volume  measured")
for i in i40
    t = r.terms[i]
    @printf("%5.0f   %7.3f  %7.3f  %7.3f  %7.3f      %7.3f  %7.3f  %7.3f  %7.3f\n", r.P[i] / 101325,
        t.henry, t.pressure, t.fugacity, t.volume,
        t.henry + t.pressure, t.henry + t.pressure + t.fugacity,
        t.henry + t.pressure + t.fugacity + t.volume, log(r.m[i]))
end
```

At 40 °C the fugacity coefficient takes 0.25 from the logarithm at 50 atm and
1.41 at 500 atm, a factor of four; the volume of the solute takes 0.06 and 0.65,
a factor of two at the highest pressure. Neither is a correction: without the
first the ideal gas misses the measurement by a factor of four, without the
second the real gas overshoots it by 80 % at 500 atm.

```@example co2
fig2 = plot(; xlabel = "P (atm)", ylabel = "ln m", legend = :bottomright, size = (700, 400))
t40 = r.terms[i40]
P40 = r.P[i40] ./ 101325
plot!(fig2, P40, [t.henry + t.pressure for t in t40]; label = "ideal gas", linestyle = :dashdot, marker = :circle, markersize = 3)
plot!(fig2, P40, [t.henry + t.pressure + t.fugacity for t in t40]; label = "+ ln φ", linestyle = :dash, marker = :circle, markersize = 3)
plot!(fig2, P40, [t.henry + t.pressure + t.fugacity + t.volume for t in t40]; label = "+ volume", marker = :circle, markersize = 3)
scatter!(fig2, P40, log.(r.m[i40]); label = "measured, 40 °C", color = :black, markersize = 5)
savefig(fig2, "validation-co2-terms.svg"); nothing # hide
```

![](validation-co2-terms.svg)

## 5. What remains

The real gas falls 2.1 % short of the measurements on average, 3.0 % in absolute
value. Two systematic patterns remain, and one point:

  - **Below the critical temperature, above saturation**, at 12 °C, the
    calculation is 2 to 4 % high: the solubility follows the fugacity of liquid
    carbon dioxide, which the equation of [PengRobinson1976](@citet), built for the
    vapor pressures of hydrocarbons, describes less closely than the gas.
  - **Above 31 °C, from 50 atm up**, it is 3 to 7 % low. The gas phase is
    taken as pure carbon dioxide, whereas it holds a few parts per thousand of
    water, 2.3 to 4.8 ‰ at 31 °C between 25 and 500 bar according to the
    measurements [Spycher2003](@citet) compile; counted, that water would lower
    the fugacity of the carbon dioxide by as much, and widen the gap rather than
    close it. The fugacity coefficient of the dense fluid, the partial molar
    volume of the solute and its activity coefficient, unit here, are the terms
    left to account for it.
  - **The point at 40 °C and 200 atm**, where the calculation falls 7.3 % short,
    is also off the trend of its own isotherm: it is within 0.3 % of the
    measurement at 35 °C, whereas at every other pressure the two isotherms are
    2.5 to 11 % apart.

[Spycher2003](@citet) fit these data, among others, with equilibrium constants
and partial molar volumes of their own and a separate equilibrium constant for
liquid carbon dioxide below 31 °C. The calculation here takes every constant from
elsewhere, keeps one standard state for the gas whatever its density, and lets
the equation of state carry the liquid.
