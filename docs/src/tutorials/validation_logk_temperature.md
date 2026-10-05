# [Equilibrium constants from 5 to 90 °C](@id sec-validation-logk-temperature)

!!! info "Before this page"
    [Validation against published data](@ref), the constants of the database
    at 25 °C against their sources.

Cemdata18 tabulates its constants at 25 °C. At another temperature they follow
from the standard properties of its species: the equation of state of
[Helgeson1981](@citet) for the aqueous species, the heat capacities of the
solids, and the water's own equation of state. The PHREEQC database
[ParkhurstAppelo2013](@cite) gives many of the same reactions an analytical
expression fitted to measured constants over a temperature range, with its
sources in comments where it gives them. The two are
independent descriptions of the same equilibria, and this page compares them
from 5 to 90 °C, reaction by reaction (`scripts/logk_temperature_check.jl`, on
the copy of `phreeqc.dat` 3.7.3 the test suite carries).

## 1. The comparison

For each reaction, its constant at 25 °C in both databases, and the difference
of their temperature dependences at 50 °C and at its largest between 5 and
90 °C: each constant taken from its own value at 25 °C, the two differences
compared, so that what the two databases already disagree on at 25 °C is not
counted twice. Where `phreeqc.dat` gives a reaction two expressions (gypsum),
the second is used, which its comment calls the better fit.

```@example logk
using ChemistryLab, DynamicQuantities, Printf
include(joinpath(pkgdir(ChemistryLab), "scripts", "logk_temperature_check.jl"))

rows = lkt_compare()
println(rpad("reaction", 22), "  log K at 25 °C        difference of the temperature dependences")
println(rpad("", 22), "  PHREEQC  Cemdata18    at 50 °C   largest, 5-90 °C")
for r in rows
    @printf("%-22s %8.3f %9.3f    %+7.3f   %+6.3f at %2d °C\n", r.name, r.phreeqc25, r.cemdata25, r.at50, r.worst, r.worst_T)
end
```

The constants of two reactions along the range, the most and the least soluble
of the calcium sulfates and carbonates:

```@example logk
println(" °C   calcite: PHREEQC Cemdata18    anhydrite: PHREEQC Cemdata18")
for T_C in (5, 25, 50, 75, 90)
    T = T_C + 273.15
    c = lkt_analytic("Calcite"); a = lkt_analytic("Anhydrite")
    @printf("%3d          %7.3f %9.3f               %7.3f %9.3f\n", T_C,
            lkt_phreeqc(c, T), lkt_cemdata(Dict("Ca+2" => 1, "CO3-2" => 1, "Cal" => -1), T),
            lkt_phreeqc(a, T), lkt_cemdata(Dict("Ca+2" => 1, "SO4-2" => 1, "Anh" => -1), T))
end
```

## 2. What the comparison says

**At 25 °C the two databases share their constants.** Ten of the twelve agree
within 0.01; anhydrite differs by 0.08 and quartz by 0.23, Cemdata18 making
quartz the more soluble.

**Up to 50 °C their temperature dependences agree within 0.07,** the largest
difference being that of the bisulfate ion, and by 0.05 at most for the
minerals.
Up to 90 °C they part by at most 0.24, for calcite, whose PHREEQC expression is
a fit over 0 to 250 °C (its comment names Ellis 1959 and
[PlummerBusenberg1982](@citet)); aragonite, which PHREEQC gives another
expression, agrees with Cemdata18 within 0.02 over the whole range. Gypsum and
anhydrite stay within 0.11, the carbonate equilibria within 0.09, the water
within 0.02. For a cement cured up to 50 °C, the equations of state of
Cemdata18 carry its constants as the fits to measurements do, to a few
hundredths of a log unit; at 90 °C the calcium carbonates and sulfates are
known to a tenth or two, less well than at 25 °C.
