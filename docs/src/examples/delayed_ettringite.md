# [Delayed ettringite formation: a CEM I 52.5 N HTS from 5 to 80 °C and back](@id ex-delayed-ettringite)

!!! info "Before this page"
    [CEM I 52.5 N HTS and CEM II/A-L 42.5 R from 0 to 60 °C](@ref ex-hydrates-temperature),
    the paste and the temperature at which monosulfate replaces ettringite.

A paste heated above about 50 °C loses ettringite to monosulfate, the sulfate it
held going to the solution; cooled, it forms ettringite again, in a hardened
paste whose pores were not made for it: the delayed ettringite formation of
steam-cured concrete. [Lothenbach2007](@citet) hydrated the sulfate-resisting
CEM I 52.5 N HTS (SRPC) at w/c 0.4 at 5, 20 and 50 °C and measured its pore
solution and its phases. This page computes the paste of the companion paper
[Lothenbach2008](@citet), at the degree of hydration 150 days give it, at each
temperature, then heated to 80 °C and brought back to 20 °C.

```@example def
using ChemistryLab, DynamicQuantities, Printf
include(joinpath(pkgdir(ChemistryLab), "scripts", "delayed_ettringite.jl"))

cs = l08t_system()
paste = Dict(20.0 => de_state(20.0; cs))
for T in (5.0, 50.0, 80.0)
    paste[T] = de_state(T; start = paste[20.0].state, cs)
end
println(" T (°C)   SO4 (mmol/L)   OH (mol/L)   ettringite   monosulfate  (g per 100 g of cement)")
for T in (5.0, 20.0, 50.0, 80.0)
    o = de_observables(paste[T])
    @printf("%6.0f %8.3f / %5s %7.3f / %5s %10.2f %13.2f\n", T, 1000o.SO4,
            isnan(de_measured("SO4", T)) ? "-" : @sprintf("%.2f", 1000de_measured("SO4", T)),
            o.OH, isnan(de_measured("OH-", T)) ? "-" : @sprintf("%.2f", de_measured("OH-", T)),
            o.ettringite, o.monosulfate)
end
```

Here / measured at 150 days. The sulfate of the solution rises with the
temperature, as measured, but from far lower: 13 to 30 times below the
measurement, with a rise of 35 from 5 to 50 °C where it was 15. The hydroxide
is half the measured one, the sodium a quarter of it and the potassium two
fifths. At 50 °C the paste holds no monosulfate yet, the transition of this
calculation being at 53 °C
([the paste against temperature](@ref ex-hydrates-temperature)), where the
measured paste held some.

At 80 °C a third of the ettringite has given way to monosulfate. Brought back
to 20 °C, the same budget gives the same equilibrium as before the heating:

```@example def
back = de_state(20.0; start = paste[80.0].state, cs)
o, o20 = de_observables(back), de_observables(paste[20.0])
@printf("ettringite %.2f g, monosulfate %.2f g, as before the heating: %.2f g and %.2f g\n",
        o.ettringite, o.monosulfate, o20.ettringite, o20.monosulfate)
println("ettringite declared from ", temperature_range(L08T_DB["ettringite"])[1], " to ",
        temperature_range(L08T_DB["ettringite"])[2], " K")
```

The equilibrium is reversible: the ettringite comes back, which is what makes
the delayed formation possible, and the calculation cannot say where or when.
At 80 °C the heat capacity of ettringite is extrapolated, Cemdata18
[Lothenbach2019](@cite) declaring it to 60 °C only
([`temperature_range`](@ref)).
