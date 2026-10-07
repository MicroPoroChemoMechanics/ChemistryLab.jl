# [The heat of CEM I 52.5 R with slag and limestone, from its degrees of reaction](@id ex-slag-heat)

!!! info "Before this page"
    [The enthalpy of a glass](@ref sec-theory-glass), which builds the enthalpy
    of the slag glass this page needs, and
    [CEM I 52.5 R with slag and limestone, integrated in time at 5, 20 and 40 °C](@ref ex-slag-temperature-pastes),
    which computes the same pastes.

[Snellings2022](@citet) measured three things on a cement of 50 % CEM I 52.5 R,
40 % slag and 10 % limestone: how far each of its clinker phases and its slag
had reacted, by X-ray diffraction; the water its hydrates bound and the
portlandite they held, by thermogravimetry; and the heat it released, by
isothermal calorimetry over 28 days. The heat of a paste is the fall of its
enthalpy, so the measurements can be set against each other once every part of
the paste has an enthalpy, the unreacted slag glass included. This page
computes the heat from the mixing to each age at the degrees of reaction the
paper measured, with nothing fitted, and compares it with the calorimetry. At
four weeks the two agree. In the first days they do not, and the page sets the
three measurements side by side to see what the difference belongs to.

```@example slagheat
using ChemistryLab, DynamicQuantities, Printf, Plots
default(framestyle = :box, grid = false)
include(joinpath(pkgdir(ChemistryLab), "scripts", "snellings2022_heat.jl"))
cs = sn22h_system()
model = cemdata18_activity_model(:KOH)
nothing # hide
```

## 1. What is computed, and what is assumed

The paste is 100 g of the binder, its materials those of
[the previous pages](@ref ex-slag-temperature-pastes): the cement by its
Rietveld phases, the slag by its calcite, its quartz and its glass, the
limestone by its calcite, quartz and dolomite. At each age each of the four
clinker phases has reacted to its own degree, one minus its content over its
content at the mixing, the contents read off [Snellings2022; Fig. 5](@cite),
which is how the authors form the degree of the clinker; the glass of the slag
has reacted to the degree given for the slag (Fig. 6b); everything else is at
equilibrium from the mixing, the sulfates and the limestone included, so that
their heat falls before the first state and is not counted; the quartz is inert.
The heat is that of [`heat_release`](@ref) from the state at mixing, with
nothing reacted, to the state at each age, per gram of Portland cement, as the
paper reports its calorimetry. The measured heat is read off Fig. 3 of the paper
(`data/literature/Snellings2022.json`). At w/b 0.5 the degrees are:

```@example slagheat
ages = (1, 2, 7, 28)
println(" T     age    alite  belite   C3A   C4AF    slag")
for T in (5, 20, 40), a in ages
    α = sn22h_phase_degrees(T, 0.5, a)
    @printf("%2d °C %3d d %7.2f %7.2f %6.2f %6.2f %7.2f\n", T, a,
            α["Alite"], α["Belite"], α["C3A"], α["C4AF"], sn22h_degree("slag", T, 0.5, a))
end
```

The hydrates and the solution have their Cemdata18 [Lothenbach2019](@cite)
enthalpies, the clinker phases those of their database records, and the glass
the enthalpy of [`glass_enthalpy`](@ref), the oxides on the scale of Cemdata18
where it holds them. Its titanium and manganese, which the system has no element
for, are set aside as the glass reacts, at the enthalpy of their crystals, which
is what the glass enthalpy counts them at; its phosphorus, left out of the
glass, is left out of the residue too.

```@example slagheat
for T in (5, 20, 40)
    g = sn22h_glass(T)
    @printf("glass at %2d °C: %9.1f J/g, of which vitrification %5.1f ± %3.1f J/g\n", T,
            ustrip(us"J/g", g.enthalpy), ustrip(us"J/g", g.vitrification), ustrip(us"J/g", g.uncertainty))
end
```

## 2. The heat at the degrees of the diffraction

At w/b 0.5, at the three temperatures, the heat computed at each age against
the heat measured:

```@example slagheat
runs = Dict{Int, Any}()
for T in (5, 20, 40)
    o = sn22h_origin(cs, T, 0.5; model)
    states = [o.state(sn22h_phase_degrees(T, 0.5, a), sn22h_degree("slag", T, 0.5, a)) for a in ages]
    runs[T] = (; o, states, q = [o.heat(rs) for rs in states],
               meas = [sn22h_measured(T, 0.5, a) for a in ages])
end
println(" T     age   computed   measured (age read)   difference")
for T in (5, 20, 40), (k, a) in enumerate(ages)
    r = runs[T]
    @printf("%2d °C %3d d %8.0f %10.0f  (%4.1f d) %12.0f %%\n", T, a, r.q[k], r.meas[k].heat,
            r.meas[k].age, 100 * (r.q[k] / r.meas[k].heat - 1))
end
```

At four weeks the computed heat is 1 to 7 % above the measured one. At one day
it is above it by a quarter at 40 °C, by half at 20 °C and by more than the
measurement itself at 5 °C, and the difference closes with age.

## 3. The bound water at the same degrees

The thermogravimetry of [Snellings2022](@citet) was run on the powders the
diffraction was measured on, the calorimetry on a paste of its own, sealed in
its ampoule from the mixing. The states of section 2 give the bound water and
the portlandite the thermogravimetry reports (Fig. 8, at w/b 0.5), over the mass
at 550 °C, every hydrate assumed to have lost all its water by then and nothing
else (`sn22h_tga` in `scripts/snellings2022_heat.jl`). The ratio of each
measurement to its computed value is:

```@example slagheat
println(" T     age   heat: measured/computed   bound water: measured/computed")
for T in (5, 20, 40), (k, a) in enumerate(ages)
    r = runs[T]
    bw = sn22h_tga(r.states[k]).bound_water
    @printf("%2d °C %3d d %18.2f %30.2f\n", T, a, r.meas[k].heat / r.q[k],
            sn22h_measured_tga("bound water", T, a) / bw)
end
```

The two measurements, on two specimens and by two techniques, find the same
fraction of the reaction the degrees of the diffraction imply, a fraction of
0.4 to 0.8 at one day and above 0.9 at four weeks. The calculation converts
the degrees into a heat and a bound water by the same hydrates; what the two
ratios share is therefore the degrees, and in the first days the heat and the
water correspond to less reaction than the diffraction reports.

## 4. One factor on the progress, fitted on the bound water

Here the diffraction is kept for how the reaction is shared among the five
constituents, and the thermogravimetry for how far it has gone. At each age
and temperature, the amounts reacted of the four clinker phases and of the slag
are those of the diffraction multiplied by one factor ``\lambda``, fitted so
that the computed bound water is the measured one (`sn22h_progress` in the
script). The heat at that progress is then compared with the
calorimetry, which nothing was fitted on, and the portlandite with the
thermogravimetry. The last column is the heat with the glass at the enthalpy
of the crystals of its composition, at the same progress.

```@example slagheat
glass = only(c for c in material_template("slag (Snellings 2022)", SN22_DB).constituents if c.name == "glass")
m_glass = sn22_value("slag_percent") * glass.mass_fraction
pc_g = sn22_value("pc_percent")
fits = Dict(T => [sn22h_progress(cs, T, a; model, origin = runs[T].o) for a in ages] for T in (5, 20, 40))
println(" T     age     λ    heat   measured  difference   portlandite  measured   as crystals")
for T in (5, 20, 40), (k, a) in enumerate(ages)
    f, r = fits[T][k], runs[T]
    q = r.o.heat(f.state)
    vit = f.λ * sn22h_degree("slag", T, 0.5, a) * m_glass * ustrip(us"J/g", sn22h_glass(T).vitrification) / pc_g
    @printf("%2d °C %3d d %6.3f %6.0f %8.0f %9.1f %% %11.2f %10.2f %11.0f\n", T, a, f.λ, q, r.meas[k].heat,
            100 * (q / r.meas[k].heat - 1), f.tga.portlandite, sn22h_measured_tga("portlandite", T, a), q - vit)
end
```

```@example slagheat
colors = Dict(5 => :blue, 20 => :green, 40 => :red)
fig = plot(; xlabel = "age (days)", ylabel = "heat (J per g of Portland cement)", legend = :bottomright,
           size = (720, 430), xscale = :log10, xticks = ([1, 2, 7, 28], ["1", "2", "7", "28"]))
for T in (5, 20, 40)
    r = runs[T]
    plot!(fig, collect(ages), r.q; color = colors[T], ls = :dash, lw = 1, label = "$T °C, degrees of the diffraction")
    plot!(fig, collect(ages), [r.o.heat(f.state) for f in fits[T]]; color = colors[T], marker = :circle,
          lw = 2, label = "$T °C, progress of the bound water")
    scatter!(fig, [m.age for m in r.meas], [m.heat for m in r.meas]; color = colors[T], marker = :square,
             markersize = 5, label = "$T °C, measured")
end
savefig(fig, "slag-heat.svg"); nothing # hide
```

![The heat released from the mixing against the age at 5, 20 and 40 °C: computed at the degrees of reaction of the diffraction, computed at the progress fitted on the bound water, and measured by calorimetry.](slag-heat.svg)

## 5. At 28 days, every paste

The thermogravimetry is reported at w/b 0.5 only. At four weeks, where the
degrees of the diffraction and the heat come closest, the nine pastes are
computed at those degrees:

```@example slagheat
println(" T    w/b   computed   measured   difference")
for T in (5, 20, 40), wb in (0.4, 0.5, 0.6)
    q = only(sn22h_heat(cs, T, wb; ages = (28,), model))
    m = sn22h_measured(T, wb, 28)
    @printf("%2d °C %.1f %9.0f %10.0f %+9.1f %%\n", T, wb, q, m.heat, 100 * (q / m.heat - 1))
end
```

The nine pastes are from 0.2 % below the calorimetry to 7.1 % above it, all
but one above, as the paste at w/b 0.5 is at each temperature, where the bound water implies 5 to
10 % less reaction than the diffraction at four weeks (section 4).

## 6. What the comparison says

At the progress the bound water gives, the computed heat follows the
calorimetry from the first day to the fourth week at the three temperatures,
1 to 12 % below it, where at the degrees of the diffraction it was 26 to 136 %
above it at one day; no parameter is fitted on the heat. The enthalpy of
vitrification of the glass is part of that agreement: with the glass as the
crystals of its composition, the heat at the same progress is lower by 6 to
62 J per gram of Portland cement, the share of the slag that has reacted, and
at four weeks at 20 °C it would be 10 % below the measurement instead of 1 %.

What the factor fitted on the bound water measures is the difference between
the progress the diffraction reports and the one the heat and the water
correspond to: more than half of the reaction at 5 °C and a third at 20 and
40 °C at one day, a tenth at most at four weeks. It does not say which
constituent the difference belongs to: one factor is applied to all five, and
the bound water alone cannot share it among them. The authors themselves put
the precision of the degree of the slag at 10 points, and note that the degrees
of the phases present in small amounts, the aluminate, the ferrite and the
belite, scatter most, the uncertainty of their contents being magnified when it
is converted into a degree [Snellings2022](@cite). The comparison here does not
decide where the difference comes from; it measures it, and shows that the
heat and the bound water agree on it.

Two limits of the comparison remain. The heat at the fitted progress is below
the measurement at every age, the most at 40 °C, so that the bound water
corresponds to a little less reaction than the heat does: either the hydrates
of the calculation hold more water per mole than those of the paste, or the
drying that stops the hydration removes some of it. And the portlandite at the
fitted progress is that of the thermogravimetry at one day at 20 and 40 °C,
below it at 5 °C in the first two days, and above it from a week on by 0.6 to
1.8 points, as on [the previous page](@ref ex-slag-temperature-pastes), where
it is set against the same measurement over six months.
