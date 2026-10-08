# [CEM I 52.5 N HTS and CEM II/A-L 42.5 R from 0 to 60 °C](@id ex-hydrates-temperature)

!!! info "Before this page"
    [Validation on pore solutions from 7 to 80 °C](@ref ex-validation-temperature),
    the temperature dependence of the solubility products tested on measured
    pore solutions, and
    [CEM I 52.5 R with slag and limestone at 5, 20 and 40 °C](@ref ex-slag-temperature),
    the rate laws at three temperatures.

[Lothenbach2008](@citet) calculated the hydrates of a sulfate-resisting
CEM I 52.5 N HTS (SRPC) and of a Portland-limestone CEM II/A-L 42.5 R (PLC)
from 0 to 60 °C, with the thermodynamic data of cemdata2007, at the degree of
hydration the clinker reaches in 150 days, kept the same at every temperature
(their Figs. 5 and 6). Above about 48 °C they find monosulfate more stable than
ettringite and monocarbonate; taking the solubility products uncertain by 0.1
log units spreads that temperature from 42 to 54 °C. This page computes the
same pastes with Cemdata18 [Lothenbach2019](@cite) and the phase list the
package uses for a Portland paste, and compares the two calculations.

## 1. The pastes

The cements are those of [Lothenbach2008; Table 1](@citet): the four clinker
phases of its normative composition, the free lime, the calcite, the calcium
sulfate (as anhydrite, the formula the table prints) and the alkali sulfates.
Each clinker phase is at the degree the law of [ParrottKilloh1984](@citet) gives
it after 150 days at 20 °C, with the constants of
[Lothenbach2008; Table 2](@cite), its w/c factor and the fineness of Table 1;
the SRPC is a paste at w/c 0.4, the PLC a mortar at 0.58. ASSUMED: the minor
oxides held in the clinker phases are released at the overall degree of the
clinker, the article giving their amounts and not how the phases share them.

```@example hydrates-temperature
using ChemistryLab, DynamicQuantities, OptimaSolver, Printf, Plots
include(joinpath(pkgdir(ChemistryLab), "scripts", "lothenbach2008_temperature.jl"))

cements = ("SRPC", "PLC")
for c in cements
    @printf("%-4s at w/c %.2f: %s; overall %.0f %%, the article about %.0f %%\n", c, l08t_water_cement(c),
            join((@sprintf("%s %.0f %%", p, 100a) for (p, a) in l08t_degrees(c)), ", "),
            100 * l08t_overall_degree(c), l08t_value(c == "SRPC" ? "overall_degree_srpc" : "overall_degree_plc"))
end
```

The same law with the same constants leaves more of the SRPC unreacted than
[Lothenbach2008](@citet) find, 71 % against about 80 %; the PLC is within two
points. The article does not detail how it applies the fineness, "on the initial
hydration", and its degrees could not be traced further: the SRPC below holds
about two cubic centimeters more clinker than its Fig. 5.

## 2. The hydrates from 0 to 60 °C

The equilibrium of each paste at every five degrees, the volumes in cm³ per
100 g of cement, the species gathered as the figures of [Lothenbach2008](@citet)
label their layers:

```@example hydrates-temperature
cs = l08t_system()
materials = Dict(c => l08t_material(c) for c in cements)
temperatures = 0.0:5.0:60.0
shown = ("ettringite", "monocarbonate", "monosulfate", "calcite", "portlandite", "hydrotalcite", "C-S-H", "iron hydroxide", "unhydrated clinker")
volumes = Dict(c => [l08t_grouped(l08t_volumes(l08t_state(c, T; cs, material = materials[c]))) for T in temperatures]
               for c in cements)
for c in cements
    println(c, ", cm³/100 g")
    println(" °C  ", join((lpad(first(g, 11), 12) for g in shown)))
    for (T, v) in zip(temperatures, volumes[c])
        @printf("%3.0f  %s\n", T, join((@sprintf("%12.2f", v[g]) for g in shown)))
    end
end
```

The same volumes stacked, from the unhydrated clinker up, as the figures of
[Lothenbach2008](@citet) draw them:

```@example hydrates-temperature
layers = reverse([g for g in shown if any(v -> v[g] > 0.01, vcat(volumes["SRPC"], volumes["PLC"]))])
panel(c) = areaplot(collect(temperatures), reduce(hcat, [[v[g] for v in volumes[c]] for g in layers]);
                    label = permutedims(layers), xlabel = "temperature (°C)", ylabel = "cm³ per 100 g of cement",
                    title = "$c, w/c $(l08t_water_cement(c))", titlefontsize = 10, legend = false,
                    palette = :tab10, lw = 0.5)
# The legend as a panel of its own, so that the two plots keep one width.
key = plot(fill(NaN, 1, length(layers)); label = permutedims(layers), palette = :tab10, lw = 8,
           framestyle = :none, legend = :left)
plot(panel("SRPC"), panel("PLC"), key; layout = @layout([a b c{0.2w}]),
     size = (1000, 420), left_margin = 5Plots.mm, bottom_margin = 5Plots.mm)
```

Between the transitions nothing changes but by fractions of a cubic
centimeter: the portlandite, the calcite, the C-S-H and the hydrotalcite are
there at every temperature, as in [Lothenbach2008](@citet).

## 3. The transition

The temperature above which each paste holds monosulfate, by bisection to
0.05 °C (`l08t_transition`):

```@example hydrates-temperature
transition = Dict(c => l08t_transition(c; cs, material = materials[c]) for c in cements)
for c in cements
    @printf("%-4s monosulfate from %.1f °C\n", c, transition[c])
end
@printf("the article: about %.0f °C, from %.0f to %.0f °C for solubility products uncertain by 0.1 log units\n",
        l08t_value("transition_temperature"), l08t_value("transition_temperature_low"), l08t_value("transition_temperature_high"))
```

This calculation puts the change five degrees above that of
[Lothenbach2008](@citet), at the upper end of the range the article gives for
the uncertainty of its own data, and at the same temperature in both cements
within a degree, as the article finds it.

## 4. Against the article on either side

The volumes [Lothenbach2008](@citet) calculate, read from their Figs. 5 and 6 at
5 and 58 °C to 0.2 cm³ (`data/literature/Lothenbach2008.json`), against those
above:

```@example hydrates-temperature
fig = literature_table(L08T, "hydrate_volumes")
for c in cements, T in (5, 58)
    v = l08t_grouped(l08t_volumes(l08t_state(c, float(T); cs, material = materials[c])))
    println(c, " at ", T, " °C     article  computed")
    for g in shown
        k = findfirst(i -> fig.cement[i] == c && ustrip(fig.temperature_C[i]) == T && fig.phase[i] == g, eachindex(fig.phase))
        a = k === nothing ? "      –" : @sprintf("%7.2f", ustrip(us"cm^3", fig.volume[k]))
        (k === nothing && v[g] < 0.01) && continue
        @printf("  %-17s %s  %8.2f\n", g, a, v[g])
    end
end
```

## 5. What the comparison says

**The sequence is that of [Lothenbach2008](@citet), the temperature five degrees
higher.** Below the transition both pastes hold ettringite and monocarbonate,
above it monosulfate and more calcite, as in the article; this calculation moves
the change from about 48 to 53 °C, within the 42 to 54 °C that the article's
0.1 log units of uncertainty allow. On either side the portlandite and the
calcite agree within 1.1 cm³ per 100 g of cement, the hydrotalcite within 0.3,
and the ettringite of both pastes below the transition within 0.6; the C-S-H is
larger here, by 2.1 to 2.7 cm³.

**The monocarbonate is where the two calculations part.** Below the transition
the package's list forms a third of the monocarbonate of
[Lothenbach2008](@citet) in the SRPC and a little over half in the PLC; above
it, the PLC of the article keeps monocarbonate and no ettringite, and the
package's keeps a little ettringite and no monocarbonate, and the SRPC more
ettringite than the article's. Two differences between the calculations act on
the aluminum and the iron the AFm and AFt phases share with the sulfate. The
article's ettringite, monocarbonate, monosulfate and hydrotalcite are solid
solutions of their aluminum and iron analogues (its footnote 1), where the
package's Portland list puts the iron in an iron hydroxide (0.6 to 0.9 cm³
here); and its clinker is less hydrated in the SRPC, while the sulfate is all
dissolved in both.

## Where to go next

[Validation on pore solutions from 7 to 80 °C](@ref ex-validation-temperature)
tests the same temperature dependence on measured solutions, and
[CEM I 52.5 R with slag and limestone at 5, 20 and 40 °C](@ref ex-slag-temperature)
the rate laws that set the degrees of hydration.
