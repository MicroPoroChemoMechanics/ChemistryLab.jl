# [A semi-adiabatic calorimeter, inside the kinetics](@id ex-semiadiabatic)

!!! info "Before this page"
    [The full Portland cement, through its pore solution](@ref ex-ionic-opc),
    whose model this page puts in a calorimeter.

An isothermal calorimeter holds the sample at one temperature, and what it
records can be read off a calculated trajectory afterwards, as the enthalpy the
states lose ([An isothermal calorimeter, read off the states](@ref sec-example-isothermal)). A
semi-adiabatic calorimeter cannot be read that way. It lets the heat of
hydration raise the temperature of the sample against the losses of the vessel,
and the temperature raises the rates of the reactions in turn. The temperature
is then an unknown of the kinetic problem, integrated with the amounts:

```math
C_{\rm tot}\,\frac{\mathrm{d}T}{\mathrm{d}t} \;=\; \dot q \;-\; \varphi(T-T_{\rm env}),
\qquad
\varphi(\Delta T) \;=\; a\,\Delta T + b\,\Delta T^2 ,
```

with `q̇` the heat the reactions release at the current temperature and `C_tot`
the heat capacity of the vessel and of everything in it. Under partial
equilibrium the heat is the rate at which the enthalpy of the whole composition
falls, the hydrates precipitated by the minimization included (see
[`cumulative_heat`](@ref)).

This page reproduces the test of [Lavergne2018](@cite) on the plain-cement mortar
`C100` at w/b = 0.5, in their calorimeter, with the model of
[the preceding page](@ref ex-ionic-opc) and nothing adjusted.

## The cell

The cell is NF EN 196-9's, as [Lavergne2018](@cite) calibrated it:
their Eq. (23) for the losses, their Table 11 for the mix. The sand keeps the
temperature moderate and takes no part in the chemistry; it enters with the
vessel and the water it absorbs as a fixed heat capacity, while the paste's own
heat capacity is summed over its composition at every step.

```@example semiad
using ChemistryLab, DynamicQuantities, OptimaSolver, OrdinaryDiffEq, Printf, Plots, Logging
gr()
include(joinpath(pkgdir(ChemistryLab), "scripts", "ionic_hydration.jl"))

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

The vessel is read as 380 J/K where the article prints "about 380 kJ/K": the
note on `CALORIMETRY_VESSEL_CP` in `scripts/ionic_hydration.jl` gives the
arithmetic, and the measured temperatures below leave no doubt either way, since
380 kJ/K would let them rise by less than a kelvin. The cell starts at 20 °C,
where the measured curve starts and where the kinetic parameters are referred.

## The run

The mortar's binder, 371 g, is integrated with the cell in the state: the
dissolution of the four clinker phases, the Gibbs minimization of everything
else at every accepted step, and the temperature.

```@example semiad
binder = mix.binder
# The solver's warnings are summarized by what they are about, printed below:
# the re-speciations that failed, and the worst element balance of the accepted
# steps, in moles.
quiet(f) = with_logger(f, NullLogger())
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

## Against the measurement

The measured temperature is read from Fig. 15(a) of [Lavergne2018](@cite), from
its maximum to 3.5 days; before, the figure draws it as crosses that overlap
those of the other mixes, and it is not transcribed
(`data/literature/Lavergne2018.json`).

```@example semiad
meas = literature_table("Lavergne2018", "semi_adiabatic_C100_wb050_temperature")
tm = ustrip.(u"d", meas.time)
Tm = meas.temperature_C
Tc = [semi.sol(x * 86400)[end] - 273.15 for x in tm]
k = argmax(Tm)
@printf("maximum: measured %.1f °C at %.2f d, computed %.1f °C at %.2f d\n",
        Tm[k], tm[k], T[j] - 273.15, t[j] / 86400)
println("  t (d)   measured   computed")
for i in 1:6:length(tm)
    @printf("  %5.2f   %7.1f    %7.1f\n", tm[i], Tm[i], Tc[i])
end
```

## What the feedback does

The same paste at 20 °C in an isothermal calorimeter gives the heat rate a
calculation without feedback would feed the cell. Integrated afterwards, as
`langavant_temperature` does, it gives the temperature the cell would reach if
the reactions ignored it.

```@example semiad
iso_cal = IsothermalCalorimeter(T_env)
iso = quiet() do
    run_ionic_hydration(; wb = mix.wb, binder_mass = binder, calorimeter = iso_cal, tend = 5 * 86400.0)
end
ti, q = heat_flow(iso.sol, iso_cal)                          # W, for 371 g of binder
# The paste's heat capacity at the start, held, per kilogram of binder as
# `langavant_temperature` takes it: the cell's fixed part is most of the total,
# and the paste's own changes by a few percent as it hydrates.
m_g = ustrip(us"g", binder)
paste = ChemicalState(semi.cs, semi.state0.n ./ (m_g / 1000); T = T_env)
T_off = langavant_temperature(ti, q ./ m_g, fill(paste, length(ti)); mix)
i = argmax(T_off)
@printf("without feedback: maximum %.1f °C at %.2f d (%d failed re-speciations, worst balance %.1e mol)\n",
        T_off[i] - 273.15, ti[i] / 86400, balance(iso)...)

plot(t ./ 86400, T .- 273.15; lw = 2, label = "computed, coupled",
     xlabel = "time [days]", ylabel = "T [°C]", size = (720, 400), legend = :topright)
plot!(ti ./ 86400, T_off .- 273.15; lw = 2, ls = :dash, label = "computed, without feedback")
scatter!(tm, Tm; ms = 3, label = "measured (Lavergne et al. 2018)")
```

The coupled calculation reaches 56.4 °C at 0.82 day, where the measurement
peaks at 52.1 °C at 0.75 day, and it stays 3 to 5 K above the measured curve
through the cooling that follows. Nothing has been adjusted on this curve: the
kinetic parameters are those of the preceding page, the cell and its losses those
[Lavergne2018](@cite) calibrated. The run without feedback, the heat flow of the
paste held at 20 °C integrated afterwards through the same cell, peaks at 40.3 °C
only, and later, at 1.03 day. The difference, 16 K out of a rise of 36 K, is the
acceleration of the reactions by the temperature they raise, through their
activation energies; a semi-adiabatic test is therefore a test of those energies
as much as of the heat, and it cannot be read off an isothermal calculation.

## See also

  - [The full Portland cement, through its pore solution](@ref ex-ionic-opc), the
    model and its assumptions.
  - [Cement clinker hydration kinetics](@ref), the same kind of cell on the
    stoichiometric formulation.
  - [An isothermal calorimeter, read off the states](@ref sec-example-isothermal)
    and [Bound water, and the thermogram it integrates to](@ref sec-example-tga),
    the measurements that are outputs of a calculation rather than part of it.
