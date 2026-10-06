# [Leaching: the C-S-H and a CEM I paste in renewed water](@id ex-leaching)

!!! info "Before this page"
    [A CEM I paste replaced by fly ash, carbonated, salted and leached](@ref ex-cement-processes),
    §5, for [`leach`](@ref), and the CEM I 42.5 N of
    [Validation against a measured paste](@ref ex-validation).

Water that renews itself around a cement paste takes its calcium. Portlandite
dissolves first, at the calcium concentration that saturates it; the C-S-H then
loses calcium incongruently, its Ca/Si falling with the calcium of the solution;
the AFm phases and ettringite dissolve on the way. This page computes both at
equilibrium, in zero dimensions: the C-S-H alone, against the solubility data
[Berner1992](@citet) compiled, and a paste leached by successive renewals of its
solution, against the zones [AdenotBuil1992](@citet) observed in a paste leached
by deionized water. Nothing is fitted.

## 1. The C-S-H alone

Lime and amorphous silica in water, at a bulk Ca/Si from 0.5 to 3, with the
C-S-H of Cemdata18 (CSHQ) and portlandite and amorphous silica beside it, at
25 °C:

```@example leaching
using ChemistryLab, DynamicQuantities, Printf
include(joinpath(pkgdir(ChemistryLab), "scripts", "berner1992_leaching.jl"))

sweep = be92_sweep(vcat(0.5:0.05:2.0, 2.1:0.1:3.0))
println("bulk Ca/Si   gel Ca/Si   Ca (mmol/L)   Si (mmol/L)     pH   beside the gel")
for p in sweep[1:5:end]
    beside = p.portlandite ? "portlandite" : p.silica ? "amorphous silica" : ""
    @printf("%10.2f %11.3f %13.2f %13.4f %6.2f   %s\n", p.bulk, p.Ca_Si, p.Ca, p.Si, p.pH, beside)
end
println("all certified: ", all(p -> p.certified, sweep))
```

The gel spans a Ca/Si from 0.68, where amorphous silica appears beside it, to
1.63, where portlandite does, at 20.3 mmol/L of calcium. Berner's Appendix A
gathers six sets of measurements on synthetic C-S-H, the Ca/Si of the solid
against the calcium and the silicon of the solution, from 17 to 30 °C. At the
Ca/Si of each measured solid within the range of the gel, the computed
concentration is interpolated and compared:

```@example leaching
m = be92_measured()
ratio(field) = [
    log10(getproperty(m, field)[i] / be92_interpolate(sweep, field, m.Ca_Si[i]))
        for i in eachindex(m.Ca_Si)
        if !ismissing(m.Ca_Si[i]) && !ismissing(getproperty(m, field)[i]) &&
            !isnan(be92_interpolate(sweep, field, m.Ca_Si[i]))
]
med(x) = sort(x)[cld(length(x), 2)]
for field in (:Ca, :Si)
    r = ratio(field)
    @printf("%-2s: %3d points, measured / computed at the median %.2f, within a factor 2 at %d\n",
            field, length(r), 10^med(r), count(abs.(r) .< log10(2)))
end
```

The calcium is reproduced: 164 of the 176 points lie within a factor of 2 of
the computed curve, the measurements 1.23 times higher at the median, within the
scatter of the six sets, which differ from one another by a factor of 2 to 3 at
a given Ca/Si. The silicon is not: CSHQ dissolves 2.3 times more silicon than
was measured at the median, and only 34 of the 101 points lie within a factor
of 2. The same model overestimates the silicon of the C-S-H syntheses in alkali
hydroxide ([Alkali uptake by C-A-S-H](@ref sec-validation-alkali-uptake)).

```@example leaching
using Plots
default(framestyle = :box, grid = false)
gel = [p for p in sweep if !p.portlandite && !p.silica]
measured(field) = [(m.Ca_Si[i], getproperty(m, field)[i]) for i in eachindex(m.Ca_Si)
                   if !ismissing(m.Ca_Si[i]) && !ismissing(getproperty(m, field)[i])]
f1 = plot(; xlabel = "Ca/Si of the solid", ylabel = "Ca (mmol/L)", yscale = :log10, legend = :bottomright)
scatter!(f1, first.(measured(:Ca)), last.(measured(:Ca)); ms = 2, label = "Berner (1992)")
plot!(f1, [p.Ca_Si for p in gel], [p.Ca for p in gel]; lw = 2, label = "CSHQ")
f2 = plot(; xlabel = "Ca/Si of the solid", ylabel = "Si (mmol/L)", yscale = :log10, legend = :topright)
scatter!(f2, first.(measured(:Si)), last.(measured(:Si)); ms = 2, label = "Berner (1992)")
plot!(f2, [p.Ca_Si for p in gel], [p.Si for p in gel]; lw = 2, label = "CSHQ")
savefig(plot(f1, f2; layout = (1, 2), size = (900, 380), left_margin = 5Plots.mm, bottom_margin = 5Plots.mm),
        "leaching-csh.svg"); nothing # hide
```

![](leaching-csh.svg)

## 2. A paste in renewed water

The CEM I 42.5 N of [LothenbachWinnefeld2006](@citet) at 28 days, its solution
removed and replaced by pure water forty times by 500 g per 100 g of cement,
then forty times by 5000 g:

```@example leaching
rows = be92_leaching()
println("certified: ", count(r -> r.certified, rows), " of ", length(rows))
println("   water (g/100 g)   Ca (mmol/L)    pH   portlandite   AFm   ettringite   C-S-H   Ca/Si")
for r in rows[[1, 3, 31, 33, 41, 44, 47, 54, 61, 81]]
    @printf("%17.0f %13.2f %6.2f %12.2f %6.2f %11.2f %8.1f %7.3f\n",
            r.water, r.Ca, r.pH, r.portlandite, r.afm, r.ettringite, r.csh, r.Ca_Si)
end
```

The masses are in g per 100 g of cement. The alkalis leave with the first
renewals, and the pH falls from 13.6 to 12.7. Portlandite then holds the calcium
at 21 mmol/L and the gel at its Ca/Si of 1.63 until it is gone; only then does
the C-S-H lose calcium, and the AFm phases and ettringite dissolve in turn:

```@example leaching
println("gone           water (g/100 g)   Ca (mmol/L)    pH   gel Ca/Si")
for e in be92_events(rows)
    @printf("%-12s %17.0f %13.2f %6.2f %11.3f\n", e.phase, e.water, e.Ca, e.pH, e.Ca_Si)
end
```

The order is the zoning of [AdenotBuil1992](@citet), read from the core of their
leached paste to its surface: the portlandite front, then the front of the
monosulfoaluminate, here the AFm of a cement that also holds calcite, then the
ettringite front. The water each step takes is not theirs: a renewed solution
removes every solute at once, where their paste lost its ions by diffusion,
faster the more mobile the ion, and the comparison is the sequence and the
Ca/Si of the gel against the calcium of the solution, not the position of the
fronts.

```@example leaching
f3 = plot(; xlabel = "water passed (g per 100 g of cement)", ylabel = "g per 100 g of cement",
          xscale = :log10, legend = :topright, size = (700, 400))
w = [max(r.water, 100.0) for r in rows]
for (field, lab) in ((:portlandite, "portlandite"), (:afm, "AFm"), (:ettringite, "ettringite"), (:csh, "C-S-H"))
    plot!(f3, w, [getproperty(r, field) for r in rows]; lw = 2, label = lab)
end
savefig(f3, "leaching-paste.svg"); nothing # hide
```

![](leaching-paste.svg)
