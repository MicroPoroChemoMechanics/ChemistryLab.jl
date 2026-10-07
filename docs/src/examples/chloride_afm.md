# [Chloride and the AFm phases: a model Portland system titrated by CaCl₂](@id ex-chloride-afm)

!!! info "Before this page"
    [Friedel's and Kuzel's salts from 0 to 99 °C](@ref ex-chloride-temperature),
    the same mixture heated, on which this page builds.

Chloride entering a paste is bound by its AFm phases: the sulfate AFm becomes
Kuzel's salt, which becomes Friedel's salt as chloride grows; with carbonate,
the monocarbonate becomes Friedel's salt directly. [Balonis2010](@citet)
equilibrated at 25 °C the model mixture of a Portland paste, 0.01 mol of C₃A,
0.015 mol of portlandite and 0.01 mol of calcium sulfate in 60 mL of water,
without and with 0.0075 mol of calcite, with CaCl₂ up to 2Cl/Al₂O₃ = 1, and
identified the phases by XRD. This page computes the titration with
Cemdata18 [Lothenbach2019](@cite), on the system of the companion page on
temperature.

```@example chloride
using ChemistryLab, DynamicQuantities, Printf
include(joinpath(pkgdir(ChemistryLab), "scripts", "balonis2010_chloride.jl"))

ratios = 0.0:0.02:1.1
free = ba10_titration(ratios)
carb = ba10_titration(ratios; calcite = true)
println("certified: ", count(r -> r.certified, free) + count(r -> r.certified, carb), " of ", 2 * length(ratios))
b = ba10_boundaries(free)
@printf("without calcite: sulfate AFm gone at %.2f, Kuzel's salt from %.2f to %.2f (2Cl/Al2O3)\n",
        b.monosulfate_gone, b.kuzel_first, b.kuzel_gone)
xrd = literature_table("Balonis2010", "xrd_phases")
row(rows, r) = rows[argmin([abs(x.ratio - r) for x in rows])]
println("system              2Cl/Al2O3   XRD              here")
for i in eachindex(xrd.phases)
    rows = xrd.system[i] == "carbonate-free" ? free : carb
    r = ustrip(xrd.cl2_al2o3_nominal[i])
    @printf("%-19s %8.1f   %-16s %s\n", xrd.system[i], r, xrd.phases[i], ba10_phases(row(rows, r)))
end
```

E is ettringite, Ms the sulfate AFm, Ks and Fs Kuzel's and Friedel's salts,
Mc monocarbonate, P portlandite and Cc calcite. Without calcite, the sulfate
AFm is gone at 0.44 here and between 0.41 and 0.45 in the calculation of
[Balonis2010](@citet), Kuzel's salt at 0.78 here and between 0.68 and 0.72
there, on a grid of 0.02; the phases of every sample are found, but for traces
of calcite, which the carbonate-free mixture held from the air, and the calcite
of the carbonate-bearing mixture before any chloride, which the equilibrium
keeps and the XRD did not report.

```@example chloride
meas = literature_table("Balonis2010", "fig8b_measured_solution")
println("2Cl/Al2O3    Ca (mmol/L)     Cl (mmol/L)        pH     (here / measured)")
for i in eachindex(meas.pH)
    r = ustrip(meas.cl2_al2o3[i]); x = row(free, r)
    @printf("%8.2f   %5.1f / %5.1f   %6.1f / %5.1f   %5.2f / %5.2f\n", r, x.Ca, ustrip(meas.Ca[i]),
            x.Cl, ustrip(meas.Cl[i]), x.pH, ustrip(meas.pH[i]))
end
@printf("solid volume, 2Cl/Al2O3 0 to 1: %+.0f %% without calcite, %+.0f %% with\n",
        100 * (row(free, 1.0).solids / row(free, 0.0).solids - 1), 100 * (row(carb, 1.0).solids / row(carb, 0.0).solids - 1))
```

The solution of the samples, read on the figure of [Balonis2010](@citet), is
reproduced while the sulfate AFm remains; once Kuzel's salt alone holds the
chloride the calculation leaves twice the measured chloride in solution, and at
2Cl/Al₂O₃ = 1 a fifth less. Binding chloride swells the solids, by about a
quarter without calcite and by a twentieth with it, against 29 % and 7 % in the
calculation of the article.
