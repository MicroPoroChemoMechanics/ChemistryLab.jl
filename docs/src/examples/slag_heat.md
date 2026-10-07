# [The heat of CEM I 52.5 R with slag and limestone, from its degrees of reaction](@id ex-slag-heat)

!!! info "Before this page"
    [The enthalpy of a glass](@ref sec-theory-glass), which builds the enthalpy
    of the slag glass this page needs, and
    [CEM I 52.5 R with slag and limestone, integrated in time at 5, 20 and 40 °C](@ref ex-slag-temperature-pastes),
    which computes the same pastes.

[Snellings2022](@citet) measured, for a cement of 50 % CEM I 52.5 R, 40 % slag
and 10 % limestone, both how far its clinker and its slag had reacted, by X-ray
diffraction and image analysis, and the heat it released, by isothermal
calorimetry over 28 days. The heat of a paste is the fall of its enthalpy, so
the two measurements can be set against each other once every part of the
paste has an enthalpy, the unreacted slag glass included. This page computes
the heat from the mixing to each age, at the degrees of reaction the paper
measured, with nothing fitted, and compares it with the calorimetry.

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
limestone by its calcite, quartz and dolomite. At each age the four clinker
phases have reacted to the degree the paper gives for the clinker, all four
alike, since the paper reports one degree for the clinker (assumed); the glass
of the slag has reacted to the degree given for the slag; everything else is at
equilibrium from the mixing, the sulfates and the limestone included, so that
their heat falls before the first state and is not counted; the quartz is inert.
The heat is that of [`heat_release`](@ref) from the state at mixing, with
nothing reacted, to the state at each age, per gram of Portland cement, as the
paper reports its calorimetry. The measured heat is read off Fig. 3 of the
paper (`data/literature/Snellings2022.json`).

The hydrates and the solution have their Cemdata18 enthalpies, the clinker
phases those of their database records, and the glass the enthalpy of
[`glass_enthalpy`](@ref), the oxides on the scale of Cemdata18 where it holds
them. Its titanium and manganese, which the system has no element for, are set
aside as the glass reacts, at the enthalpy of their crystals, which is what the
glass enthalpy counts them at; its phosphorus, left out of the glass, is left
out of the residue too.

```@example slagheat
for T in (5, 20, 40)
    g = sn22h_glass(T)
    @printf("glass at %2d °C: %9.1f J/g, of which vitrification %5.1f ± %3.1f J/g\n", T,
            ustrip(us"J/g", g.enthalpy), ustrip(us"J/g", g.vitrification), ustrip(us"J/g", g.uncertainty))
end
```

## 2. The heat against the calorimetry

At w/b 0.5, at the three temperatures, with the heat the glass would release if
it were the crystals of its composition, the part of the enthalpy of
vitrification:

```@example slagheat
ages = (1, 2, 7, 28)
glass = only(c for c in material_template("slag (Snellings 2022)", SN22_DB).constituents if c.name == "glass")
m_glass = sn22_value("slag_percent") * glass.mass_fraction
pc_g = sn22_value("pc_percent")
runs = Dict{Int, Any}()
for T in (5, 20, 40)
    q = sn22h_heat(cs, T, 0.5; ages, model)
    v = ustrip(us"J/g", sn22h_glass(T).vitrification)
    vit = [sn22h_degree("slag", T, 0.5, a) * m_glass * v / pc_g for a in ages]
    meas = [sn22h_measured(T, 0.5, a) for a in ages]
    runs[T] = (; q, vit, meas)
end
println(" T   age   computed   as crystals   measured (age read)")
for T in (5, 20, 40), (k, a) in enumerate(ages)
    r = runs[T]
    @printf("%2d °C %3d d %8.0f %12.0f %10.0f  (%.1f d)\n", T, a, r.q[k], r.q[k] - r.vit[k], r.meas[k].heat, r.meas[k].age)
end
```

```@example slagheat
colors = Dict(5 => :blue, 20 => :green, 40 => :red)
fig = plot(; xlabel = "age (days)", ylabel = "heat (J per g of Portland cement)", legend = :bottomright,
           size = (720, 430), xscale = :log10, xticks = ([1, 2, 7, 28], ["1", "2", "7", "28"]))
for T in (5, 20, 40)
    r = runs[T]
    plot!(fig, collect(ages), r.q; color = colors[T], marker = :circle, lw = 2, label = "$T °C, computed")
    plot!(fig, collect(ages), r.q .- r.vit; color = colors[T], ls = :dash, lw = 1, label = "$T °C, glass as crystals")
    scatter!(fig, [m.age for m in r.meas], [m.heat for m in r.meas]; color = colors[T], marker = :square,
             markersize = 5, label = "$T °C, measured")
end
savefig(fig, "slag-heat.svg"); nothing # hide
```

![The heat released from the mixing against the age at 5, 20 and 40 °C, computed at the measured degrees of reaction, with the glass at its enthalpy and as crystals, and measured.](slag-heat.svg)

## 3. At 28 days, every paste

```@example slagheat
println(" T    w/b   computed   measured   difference")
for T in (5, 20, 40), wb in (0.4, 0.5, 0.6)
    q = only(sn22h_heat(cs, T, wb; ages = (28,), model))
    m = sn22h_measured(T, wb, 28)
    @printf("%2d °C %.1f %9.0f %10.0f %+9.1f %%\n", T, wb, q, m.heat, 100 * (q / m.heat - 1))
end
```

## 4. What the comparison says

At 28 days the computed heat falls within 4 % of the measured one for the nine
pastes, from 5 to 40 °C and from w/b 0.4 to 0.6, and it does so with the glass
at its enthalpy: as the crystals of its composition, the slag would release
between 35 and 70 J per gram of Portland cement less, and the heat would be 6 to
13 % short. The enthalpy of vitrification of the glass is not a correction at
the margin here, it is the part of the heat of a slag that a crystalline
reference misses, and the measured glasses give it without a parameter fitted
on a cement.

At one and two days the computed heat runs ahead of the measurement, by 80 to
130 J per gram at one day, and the gap closes with age. It is not resolved
here. It is the same at the three water-to-binder ratios, so the water is not
its cause, and it is largest in proportion at 5 °C, where the computed heat at
one day is two to three times the measured one, for a clinker the diffraction
finds reacted to 24 to 29 %. The two measurements of the paper do not count
the same thing at early age: the calorimeter integrates the heat from the
moment the paste, just mixed, is sealed in its ampoule, the diffraction reports the clinker phases missing
from a sample stopped at that age, and the four phases are put here at the one
degree reported for the clinker, which weighs the aluminate phase, the one that
releases the most heat per gram, as much as the alite.
