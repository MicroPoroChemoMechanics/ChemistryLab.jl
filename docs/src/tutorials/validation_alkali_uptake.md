# [Alkali uptake by C-A-S-H (L'Hôpital et al. 2016)](@id sec-validation-alkali-uptake)

!!! info "Before this page"
    [Validation against published data](@ref), whose section *Alkali uptake
    by C-S-H* compares the isotherms of [HongGlasser1999](@citet), on which the
    alkali end members of CSHQ were fitted.

[LHopital2016b](@citet) synthesized C-A-S-H, with an Al/Si of 0.05, at six Ca/Si
from 0.6 to 1.6, and three C-S-H without aluminum: 2 g of lime, silica fume and
monocalcium aluminate in 90 mL of water or of a KOH or NaOH solution, from 0.01
to 0.5 mol/L, each sample equilibrated at 20 °C for 91, 182 or 364 days. They
report the solutions, every element and the pH (their Appendix B), and the
C-S-H, its Ca/Si and the alkali it took up, from the fall of the dissolved
alkali (Appendix C). This page computes the same batches at equilibrium with the
two models of the gel that carry the alkalis, and compares.

!!! warning "Which model was fitted to these data"
    The alkali end members of CSHQ [Kulik2011](@cite) were fitted on the
    isotherms of [HongGlasser1999](@citet) [Lothenbach2019; §2.7](@cite), not on
    these: for CSHQ the comparison is a prediction. CASH+NK was fitted on these
    very data, the NaOH series at 0.05 and 0.1 mol/L and the KOH series from
    0.01 to 0.5 mol/L, among others [Miron2022a; Table 4](@cite): for CASH+NK it
    checks that the model is applied as its authors fitted it.

## 1. The batches

The solids of each batch, in g per 90 mL of solution (their Appendix A):

```@example alkali
using ChemistryLab, DynamicQuantities, OptimaSolver, Printf
include(joinpath(pkgdir(ChemistryLab), "scripts", "lhopital2016_alkali.jl"))

mix = lh16_table("mixing_proportions")
println("Al/Si  Ca/Si   CaO, g   SiO2, g   CaO·Al2O3, g")
for i in eachindex(mix.CaO)
    @printf("%5.2f  %5.1f  %7.3f  %8.3f  %10.3f\n", mix.Al_Si_target[i], mix.Ca_Si_target[i],
            ustrip(us"g", mix.CaO[i]), ustrip(us"g", mix.SiO2[i]), ustrip(us"g", mix.CaO_Al2O3[i]))
end
```

Each batch is computed as it was prepared (`lh16_state`): the lime, the silica
and the aluminate as the Cemdata18 records of the three solids, which dissolve,
in 90 g of water with the hydroxide dissolved in it, at 20 °C. Beside the gel the
system holds portlandite, amorphous silica and the calcium aluminate hydrates of
Cemdata18, hydrogarnet and its siliceous member, strätlingite and
microcrystalline gibbsite among them, for the aluminum neither gel takes. The
activity model is the one the authors used to speciate their solutions, the
extended Debye–Hückel equation of Cemdata18 with the ion size and the `b_γ` of
the hydroxide of the batch [Lothenbach2019](@cite).

ASSUMED: the 90 mL of solution are 90 g of water; the hydroxide solutions are
pure, where the potassium ones carry 0.1 to 3.5 mmol/L of sodium.

## 2. The two gels

CSHQ with the alkali end member of the hydroxide of the batch, `KSiOH` or
`NaSiOH`, and CASH+NK [Miron2022a](@cite) on `cemdata18-cashplus.json`, whose
twelve end members carry both alkalis, without the neutral aqueous complexes of
the alkalis, as its authors fitted it. The three times of one solution are three
samples of one composition, which an equilibrium does not tell apart: each of
the 49 compositions is computed once, each series of one Ca/Si walked down from
its most concentrated solution, each solve starting from the last.

```@example alkali
results = Dict{String, Any}()
for gel in ("CSHQ", "CASH+NK")
    results[gel] = lh16_compute(gel)
end
batches = sort(collect(keys(results["CSHQ"])); by = b -> (b[3], b[2], b[1], b[4]))
for gel in ("CSHQ", "CASH+NK")
    println(gel, ": ", length(results[gel]), " compositions, ",
            count(b -> results[gel][b].certified, batches), " certified")
end
```

## 3. The alkali the gel takes up

The alkali over silicon of the gel, computed over measured, in the batches
whose every measurement is positive. Above Ca/Si 1.2 the gel takes up so little
that the fall of the dissolved alkali, the measured quantity, is within its own
error, and several measurements are negative.

```@example alkali
println("              CSHQ             CASH+NK")
for ca in (0.6, 0.8, 1.0, 1.2, 1.4, 1.6)
    q, n = lh16_uptake_ratios(results["CSHQ"], ca), lh16_uptake_ratios(results["CASH+NK"], ca)
    @printf("Ca/Si %.1f  %5.2f to %5.2f    %5.2f to %5.2f   (%d batches)\n", ca,
            minimum(q), maximum(q), minimum(n), maximum(n), length(q))
end
```

```@raw html
<details><summary>Every batch: the alkali over silicon of the gel, and the dissolved alkali and silicon</summary>
```

```@example alkali
mean(x) = sum(x) / length(x)
println("                         alkali/Si of the gel                    alkali, mmol/L             Si, mmol/L")
println("Al/Si Ca/Si       mol/L  measured             CSHQ  CASH+NK    measured  CSHQ  CASH+NK    measured   CSHQ  CASH+NK")
for b in batches
    b[3] == "none" && continue
    m = lh16_measured(b, :alkali_Si)
    a, s = lh16_measured_solution(b, b[3] == "NaOH" ? "Na" : "K"), lh16_measured_solution(b, "Si")
    @printf("%4.2f  %3.1f %-4s %5.2f  %-19s %5.3f  %5.3f   %8.1f %6.1f %6.1f   %8.3f %6.3f %6.3f\n",
            b[1], b[2], b[3], b[4], join((@sprintf("%.3g", x) for x in m), " "),
            results["CSHQ"][b].alkali_Si, results["CASH+NK"][b].alkali_Si,
            isempty(a) ? NaN : mean(a), results["CSHQ"][b].alkali, results["CASH+NK"][b].alkali,
            isempty(s) ? NaN : mean(s), results["CSHQ"][b].Si, results["CASH+NK"][b].Si)
end
```

```@raw html
</details>
```

## 4. The solution

The dissolved alkali, silicon and calcium, computed over the mean of the
measured ones, as the smallest, the median and the largest of the batches, and
the difference of pH; below and above Ca/Si 1.1, where the gel takes up the
alkali and where it hardly does.

```@example alkali
med(x) = sort(x)[cld(length(x), 2)]
span(x) = @sprintf("%5.2f %5.2f %6.2f", minimum(x), med(x), maximum(x))
println("                     alkali               Si                    Ca               |ΔpH|")
println("                     min   median max     min   median max      min   median max   median  max")
for gel in ("CSHQ", "CASH+NK"), low in (true, false)
    r = lh16_solution_ratios(results[gel]; low)
    @printf("%-8s Ca/Si %s  %s    %s   %s    %4.2f  %4.2f\n", gel, low ? "< 1.1" : "> 1.1",
            span(r.alkali), span(r.Si), span(r.Ca), med(r.ΔpH), maximum(r.ΔpH))
end
```

## 5. What the comparison says

**CSHQ takes up too little alkali where the gel is poor in calcium.** At a Ca/Si
of 0.8 it takes up a quarter to a half of what was measured, at 0.6 and 1.0 as
little as an eighth and a third; only in 0.5 mol/L KOH does it reach the
measurement. This is the opposite of what CSHQ does on the isotherms it was
fitted on, where it binds more than the gel did at 46 of 48 points
([Validation against published data](@ref)): its two alkali end members, fixed
on one dataset, cannot follow both. The alkali it leaves in solution, at the
median a third more than was measured below Ca/Si 1.1, raises the pH, by up to
0.6 in 0.05 mol/L KOH at Ca/Si 0.6. Above Ca/Si 1.1 the gel takes up little and
CSHQ follows the dissolved alkali within 13 %. The largest difference of pH,
0.97, is the batch without alkali at Ca/Si 0.8, measured at 10.27, which CASH+NK
misses by 0.55 as well.

**CSHQ also dissolves too much silicon at a high pH**, at the median 2.5 times
the measurement below Ca/Si 1.1 and 3.6 times above, up to thirty times in the
most concentrated hydroxides. The authors found the same with the same model
(their Section 3.2, Fig. 10) and concluded that the model, built on alkali-free
C-S-H, has to be extended to high pH.

**CASH+NK reproduces the data it was fitted on.** Its alkali over silicon is
0.62 to 1.55 times the measured one at every Ca/Si up to 1.2, within the error
the authors give the indirect method at high concentration, up to 100 %; the
dissolved alkali is at the median the measured one, and the silicon within 4 %
of it below Ca/Si 1.1 and 38 % above. The calcium is where both models differ
most from the measurement, high below Ca/Si 1.1, by a factor of 2.4 at the
median with CSHQ and 1.7 with CASH+NK, and low above.

For a user, the page draws the line this way: in a gel of a Portland cement,
above Ca/Si 1.4, the two models take up the same small amount of alkali and give
the same solution; in a gel poor in calcium, below a Ca/Si of about 1.1, the gel of
a blend rich in silica fume, CSHQ underestimates the alkali the gel holds, and
CASH+NK is the model whose alkali uptake was fitted there.

## Where to go next

[The models of the C-S-H gel](@ref sec-csh-models) lists the gels of the
package and the databases they come from, and
[CEM II/B-V and CEM II/B-M (V-LL), with their CEM I, integrated in time](@ref ex-ternary-kinetics)
runs CASH+NK in the pastes of a blended cement.
