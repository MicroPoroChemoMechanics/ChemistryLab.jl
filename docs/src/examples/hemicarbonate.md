# [Hemicarbonate and monocarbonate: the carbonate of the AFm phases](@id ex-hemicarbonate)

!!! info "Before this page"
    [Carbonation of a CEM I 52.5 N and of its blends](@ref sec-cement-carbonation), for the carbonation of
    a whole paste, of which this is the part that happens in the AFm.

A paste that takes up carbonate turns its AFm phases from hydroxide to
hemicarbonate, then to monocarbonate, before calcite forms. [Georget2022](@citet)
hydrated 10 g of C₃A with portlandite and calcite in the proportions of these
phases, 28 days at 20 °C, replacing portlandite by calcite from one sample to
the next, with ζ the molar fraction of calcium carbonate among the two; they
published with their data set [Georget2021data](@cite) the same series computed
with GEMS and Cemdata18 in 101 steps. This page computes those steps, the water
of the GEMS run being unstated and taken at the 2.7 times the solids of the
experiments.

```@example hemicarbonate
using ChemistryLab, DynamicQuantities, Printf
include(joinpath(pkgdir(ChemistryLab), "scripts", "georget2022_hemicarbonate.jl"))

rows = ge22_series()
t = ge22_steps()
g(f, k) = ustrip(us"g", getproperty(t, f)[k])
println("certified: ", count(r -> r.certified, rows), " of ", length(rows))
println("    ζ    katoite    hemicarbonate   monocarbonate      calcite     portlandite   (g, here / GEMS)")
for k in 1:10:101
    r = rows[k]
    @printf("%5.3f  %5.2f/%5.2f   %6.2f/%6.2f   %6.2f/%6.2f   %5.2f/%5.2f   %5.2f/%5.2f\n", r.zeta,
            r.katoite, g(:katoite, k), r.hemicarbonate, g(:hemicarbonate, k), r.monocarbonate, g(:monocarbonate, k),
            r.calcite, g(:calcite, k), r.portlandite, g(:portlandite, k))
end
gap = maximum(maximum(abs(getproperty(rows[k], f) - g(f, k)) for f in (:katoite, :hemicarbonate, :monocarbonate, :calcite, :portlandite)) for k in eachindex(rows))
@printf("largest difference over the 101 steps: %.3f g\n", gap)
```

The two calculations agree to 0.03 g of a solid of 20 g at every step. Below
ζ = 0.452 the carbonate is all in hemicarbonate, beside katoite; up to 0.857
hemicarbonate turns into monocarbonate; above, calcite remains:

```@example hemicarbonate
b = ge22_breakpoints(rows)
@printf("katoite gone at ζ = %.4f, hemicarbonate gone at ζ = %.4f\n", b.katoite_gone, b.hemicarbonate_gone)
meas = literature_table("Georget2022", "carbonate_series")
println("  ζ     measured (Table 7)      here")
for i in eachindex(meas.pH)
    z = ustrip(meas.zeta_CO3_nominal[i])
    @printf("%4.2f   %-22s %s\n", z, meas.identified_phases[i], join(ge22_phases(rows, z), ", "))
end
@printf("pH %.2f here, %.2f to %.2f by GEMS, %.2f to %.2f measured\n", rows[1].pH, minimum(ustrip.(t.pH)),
        maximum(ustrip.(t.pH)), minimum(ustrip.(meas.pH)), maximum(ustrip.(meas.pH)))
```

The phases are those the samples held, but for the calcite found beside
monocarbonate and hemicarbonate at ζ = 0.75, which the equilibrium has not yet
formed there. The pH is 0.05 to 0.1 below both: the sodium the calcium
carbonate brought, which no one measured, is not in this calculation.
