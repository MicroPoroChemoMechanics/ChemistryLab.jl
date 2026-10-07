# [Aluminum uptake by C-S-H (L'Hôpital et al. 2016)](@id sec-validation-aluminum-uptake)

!!! info "Before this page"
    [Alkali uptake by C-A-S-H](@ref sec-validation-alkali-uptake), the same
    syntheses in alkali hydroxide solutions.

[LHopital2016a](@citet) synthesized C-S-H and C-A-S-H at six Ca/Si from 0.6 to
1.6 and an Al/Si from 0 to 0.33: 2 g of lime, silica fume and monocalcium
aluminate in 90 mL of water, equilibrated at 20 °C for six months. They report
the C-S-H by mass balance, its Ca/Si and Al/Si, the other solids weighed by
thermogravimetry and diffraction (their Appendix B), and the solutions
(Appendix D). Their finding is that the gel takes up all the aluminum up to an
Al/Si of about 0.05, and that beyond it strätlingite and katoite precipitate and
hold the Al/Si of the gel near 0.15 ± 0.05, whatever its Ca/Si.

This page computes the 34 syntheses at equilibrium with the two models of the
gel the package ships for Cemdata18: `CNASH_ss` [Myers2014](@cite), the one that
takes aluminum, and `CSHQ` [Kulik2011](@cite), which takes none. Neither was
fitted on these data.

## 1. The syntheses

```@example aluminum
using ChemistryLab, DynamicQuantities, OptimaSolver, Printf
include(joinpath(pkgdir(ChemistryLab), "scripts", "lhopital2016_aluminum.jl"))

mix = lh16a_table("mixing_proportions")
for ca in unique(mix.Ca_Si_target)
    println("Ca/Si ", ca, ": Al/Si ", join([al for (c, al) in zip(mix.Ca_Si_target, mix.Al_Si_target) if c == ca], ", "))
end
```

Each synthesis is computed as it was prepared (`lh16a_state`): the lime, the
silica fume and the monocalcium aluminate of [LHopital2016a; Appendix A](@cite),
as the Cemdata18 records of the three solids, which dissolve, in 90 g of water
at 20 °C, beside portlandite, amorphous silica and the calcium aluminate
hydrates of Cemdata18, among them katoite with its siliceous member,
strätlingite and microcrystalline gibbsite. The activity model is the extended
Debye–Hückel equation of Cemdata18 [Lothenbach2019](@cite).

ASSUMED: the 90 mL are 90 g of water; the comparison is with the samples
hydrated for 182 days, or for 364 where a synthesis has no other.

```@example aluminum
results = Dict{String, Any}()
for gel in ("CNASH_ss", "CSHQ")
    results[gel] = lh16a_compute(gel)
end
batches = sort(collect(keys(results["CSHQ"])))
for gel in ("CNASH_ss", "CSHQ")
    println(gel, ": ", length(results[gel]), " syntheses, ", count(k -> results[gel][k].certified, batches), " certified")
end
```

## 2. The aluminum the gel takes up

The Al/Si of the gel, measured and computed with `CNASH_ss`, at each Ca/Si (the
rows) and target Al/Si (the columns); `CSHQ` takes none.

```@example aluminum
targets = sort(unique(mix.Al_Si_target))[2:end]
cell(x) = ismissing(x) ? "    –" : @sprintf("%5.3f", x)
println("Ca/Si    ", join((@sprintf("   Al/Si %-5s     ", t) for t in targets)))
println("         ", join(("measured  CNASH_ss  " for _ in targets)))
for ca in unique(mix.Ca_Si_target)
    cells = map(targets) do al
        haskey(results["CNASH_ss"], (ca, al)) || return " "^20
        @sprintf("%s     %5.3f    ", cell(lh16a_measured(ca, al, :Al_Si)), results["CNASH_ss"][(ca, al)].Al_Si)
    end
    println(rpad(ca, 9), join(cells))
end
```

## 3. Where the rest of the aluminum goes

At an Al/Si of 0.2, the other solids, in wt.%: measured (by thermogravimetry
and diffraction, on the freeze-dried powder), and computed with each model
(over every solid of the equilibrium, its water included, `lh16a_observables`).
Then, with each model, where the aluminum is.

```@example aluminum
println("            katoite            strätlingite       portlandite")
println("Ca/Si   meas  CNASH  CSHQ    meas  CNASH  CSHQ    meas  CNASH  CSHQ")
fmt(x) = ismissing(x) ? "   –" : @sprintf("%4.1f", x)
for ca in unique(mix.Ca_Si_target)
    k = (ca, 0.2)
    a, b = results["CNASH_ss"][k], results["CSHQ"][k]
    @printf("%-5s  %s   %4.1f  %4.1f    %s   %4.1f  %4.1f    %s   %4.1f  %4.1f\n", ca,
            fmt(lh16a_measured(k..., :katoite)), a.katoite, b.katoite,
            fmt(lh16a_measured(k..., :straetlingite)), a.straetlingite, b.straetlingite,
            fmt(lh16a_measured(k..., :portlandite)), a.portlandite, b.portlandite)
end
println()
println("CSHQ at Al/Si 0.03, wt.%   strätlingite   katoite")
for ca in unique(mix.Ca_Si_target)
    v = results["CSHQ"][(ca, 0.03)]
    @printf("Ca/Si %.1f %26.1f %9.1f\n", ca, v.straetlingite, v.katoite)
end
println()
println("aluminum, % of the whole   in the gel   katoite   strätlingite   gibbsite   other solids   solution")
for gel in ("CNASH_ss", "CSHQ"), ca in (0.6, 1.0, 1.6)
    s = lh16a_aluminum_split(results[gel][(ca, 0.2)].state, gel)
    @printf("%-8s Ca/Si %.1f        %8.1f %9.1f %14.1f %10.1f %14.1f %10.2f\n", gel, ca,
            s.gel, s.katoite, s.straetlingite, s.gibbsite, s.other, s.solution)
end
```

## 4. The calcium of the gel

The Ca/Si of the gel and the portlandite of the syntheses at Al/Si 0.05, where
the gel holds the aluminum alone in every model that takes it:

```@example aluminum
println("Ca/Si target    gel Ca/Si: measured  CNASH_ss  CSHQ      portlandite, wt.%: measured  CNASH_ss  CSHQ")
for ca in unique(mix.Ca_Si_target)
    k = (ca, 0.05)
    a, b = results["CNASH_ss"][k], results["CSHQ"][k]
    q = lh16a_measured(k..., :Ca_Si_qualifier)
    @printf("%-14s %16s %9.2f %6.2f %24s %9.1f %6.1f\n", ca,
            (q == "lower_bound" ? "≥" : "") * @sprintf("%.2f", lh16a_measured(k..., :Ca_Si)), a.Ca_Si, b.Ca_Si,
            fmt(lh16a_measured(k..., :portlandite)), a.portlandite, b.portlandite)
end
```

## 5. What the comparison says

**CNASH_ss reproduces the uptake.** Up to an Al/Si of 0.05 its gel takes all the
aluminum, as the syntheses do: 0.030 and 0.050 from Ca/Si 1.0 up, 0.032 and
0.053 at 0.8. At Ca/Si 0.6
it gives 0.042 and 0.071 where 0.031 and 0.052 were measured: its gel holds less
silicon there, at a Ca/Si of 0.85 where the measured one is at least 0.67, the
amorphous silica taking the rest, and the same aluminum is more per silicon.
Beyond, its Al/Si stops between 0.10 and 0.12 whatever the Ca/Si, which is the
finding of [LHopital2016a](@citet), 0.15 ± 0.05; the measured values run from
0.05 to 0.23, with an error the authors put at ±0.1 when other phases are
present.

**It does not reproduce where the rest goes, nor the calcium of the gel above
Ca/Si 1.2.** At an Al/Si of 0.2 the syntheses hold strätlingite at every Ca/Si,
3.0 to 13 wt.%, and katoite, 2.3 to 5.9. `CNASH_ss` forms no strätlingite at
all: below Ca/Si 1.2 the aluminum its gel leaves goes to gibbsite, 61 % of it at
Ca/Si 0.6 and 39 % at 1.0, and above to katoite, half of it at 1.6. Its gel
stops at a Ca/Si of 1.19, where the measured one reaches 1.32 and 1.43, and the
calcium it does not take precipitates as portlandite, 3.4 and 10.5 wt.% at Ca/Si
1.4 and 1.6 where the syntheses hold none and 1.8. This is the limit of its end
members that [the blended pastes](@ref sec-validation-blended-gels) of
[DeWeerdt2011](@citet) also meet.

**CSHQ holds the calcium and none of the aluminum.** Its gel differs from the
measured Ca/Si by 0.052 at most, at Ca/Si 1.2, and no portlandite forms, where
1.8 wt.% was measured at Ca/Si 1.6. All the aluminum goes out of the gel, from
the smallest Al/Si: to gibbsite at Ca/Si 0.6, to strätlingite at 0.8 and 1.0 and
to katoite from 1.2, where the syntheses show neither below an Al/Si of 0.1.

For a user, the page draws the line this way: where the aluminum held by the
C-S-H matters, in a blend with slag, metakaolin or fly ash, `CNASH_ss` carries
it in the right amount, at a Ca/Si that cannot exceed about 1.2; `CSHQ` carries
the calcium of a Portland cement's gel and puts the aluminum in the hydrates.
Neither forms the strätlingite these syntheses hold beside the gel.

## Where to go next

[The models of the C-S-H gel](@ref sec-csh-models) lists the gels of the
package, and [Alkali uptake by C-A-S-H](@ref sec-validation-alkali-uptake)
compares the two that carry the alkalis on the companion paper
[LHopital2016b](@cite).
