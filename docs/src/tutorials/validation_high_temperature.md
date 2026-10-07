# [Water and quartz from 0 to 1000 °C](@id sec-validation-high-temperature)

!!! info "Before this page"
    [Thermochemistry](@ref sec-theory-water-eos), §3, for the standard state of
    the solvent and [the domain of the HKF equations](@ref sec-theory-hkf-domain).

The standard states of the ThermoFun databases reach 1000 °C and 5000 bar: the
solvent by the equation of state of water of [Haar1984](@citet), the aqueous
species by the HKF equations [Helgeson1981, Shock1992](@cite), the minerals by
their heat-capacity functions. This page takes two equilibria that involve
nothing but water and, for the second, one mineral, and compares them over that
range with formulations fitted to the measurements: the ionization constant of
water of [BanduraLvov2006](@citet) and the solubility of quartz of
[Manning1994](@citet). The constants are those of slop98 [Johnson1992](@cite),
the database written for this range; nothing is fitted.

## 1. The ionization constant of water

```@example hightemp
using ChemistryLab, DynamicQuantities, Printf
include(joinpath(pkgdir(ChemistryLab), "scripts", "water_high_temperature.jl"))

r = hw_bandura_lvov()
Δ = r.ours .- r.ref
@printf("%d states, water denser than 350 kg/m³ and at most 500 MPa\n", length(Δ))
println("  T (°C)  states   largest |Δ|   mean Δ")
for T_C in sort(unique(r.T_C))
    i = findall(==(T_C), r.T_C)
    @printf("%7.0f %7d %13.2f %+9.2f\n", T_C, length(i), maximum(abs, Δ[i]), sum(Δ[i]) / length(i))
end
```

``Δ`` is this package's pKw minus theirs. From 25 to 400 °C the two differ by
0.20 at most at every state of the table, close to the 0.16 standard deviation of
the fit of [BanduraLvov2006](@citet) to the measurements, and within 0.3 at
0 °C. Above 400 °C the difference grows with the temperature, this package
giving the higher pKw, that is the less dissociated water: 0.4 on average at
600 °C, 0.8 at 800 °C, 1.4 at 1000 °C, where the table reaches densities of 0.4
to 0.6 g/cm³.

```@example hightemp
using Plots
default(framestyle = :box, grid = false)
fig = plot(; xlabel = "T (°C)", ylabel = "pKw", legend = :topright, size = (700, 430))
for (k, P_MPa) in enumerate((25.0, 100.0, 300.0, 500.0))
    i = findall(P -> P ≈ P_MPa * 1e6, r.P)
    keep = [T for T in range(0, 1000; length = 101) if hw_density(T + 273.15, P_MPa * 1e6) >= 0.35]
    plot!(fig, keep, [hw_pkw(T + 273.15, P_MPa * 1e6) for T in keep]; color = k, label = "$(Int(P_MPa)) MPa")
    scatter!(fig, r.T_C[i], r.ref[i]; color = k, markersize = 3, label = false)
end
savefig(fig, "validation-high-temperature-pkw.svg"); nothing # hide
```

![](validation-high-temperature-pkw.svg)

The lines are this package, the markers the formulation of
[BanduraLvov2006](@citet). Along each isobar pKw falls with the temperature,
then rises again as the expanding solvent screens the ions less.

## 2. The solubility of quartz

```@example hightemp
println("  T (°C)   P (bar)   log m (here)   log m (Manning)")
for (T_C, P_bar) in ((25, 1), (100, 1), (200, 500), (300, 1000), (400, 2000), (500, 2000), (600, 5000), (800, 5000))
    T, P = T_C + 273.15, P_bar * 1e5
    @printf("%7d %9d %13.3f %15.3f\n", T_C, P_bar, hw_quartz_log_m(T, P), hw_manning(T, P))
end
```

```@example hightemp
fig2 = plot(; xlabel = "T (°C)", ylabel = "log m SiO₂(aq)", legend = :bottomright, size = (700, 430))
for (k, P_bar) in enumerate((500, 1000, 2000, 5000))
    P = P_bar * 1e5
    Ts = [T for T in range(25, 900; length = 90) if hw_density(T + 273.15, P) >= 0.35]
    plot!(fig2, Ts, [hw_quartz_log_m(T + 273.15, P) for T in Ts]; color = k, label = "$(P_bar) bar")
    plot!(fig2, Ts, [hw_manning(T + 273.15, P) for T in Ts]; color = k, linestyle = :dash, label = false)
end
savefig(fig2, "validation-high-temperature-quartz.svg"); nothing # hide
```

![](validation-high-temperature-quartz.svg)

Solid lines, the database; dashed lines, the equation of [Manning1994](@citet).
The two differ by 0.004 to 0.022 in ``\log m`` up to 400 °C, by 0.06 at 500 °C
and 2 kbar, and part by 0.12 at 800 °C and 5 kbar, the database dissolving less.
The change of slope at 575 °C is the transition of α- to β-quartz, which the
record of quartz places at the top of its first heat-capacity interval and which
the package carries into the second.
