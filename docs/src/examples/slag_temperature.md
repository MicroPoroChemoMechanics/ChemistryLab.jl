# [CEM I 52.5 R with slag and limestone at 5, 20 and 40 °C](@id ex-slag-temperature)

!!! info "Before this page"
    [CEM II/B-V and CEM II/B-M (V-LL), with their CEM I, integrated in time](@ref ex-ternary-kinetics),
    the law of the clinker tested on another clinker at 20 °C.

[Snellings2022](@citet) blended a CEM I 52.5 R with 40 % ground granulated
blast-furnace slag and 10 % limestone, each ground on its own, and cured the
pastes under water at 5, 20 and 40 °C, at w/b 0.4, 0.5 and 0.6. They followed
the degrees of reaction of the clinker and of the slag by X-ray diffraction for
six months (their Fig. 6), the slag's to within 10 points by their estimate.
The laws the package gives these two constituents carry the temperature by an
Arrhenius factor, and this page tests them at the three temperatures. Nothing
in the law of the clinker is fitted. The slag has no published constants of its
own: they are fitted here, and their activation energy is compared with the one
the authors fit with another law.

## 1. The measurements

Fig. 6 is a raster image. Each of its 108 markers was found by its color and
its shape, and placed where the lines drawn through it meet its age, to 0.3
point (`data/literature/Snellings2022.json`, under `digitization`).

```@example slag-temperature
using ChemistryLab, DynamicQuantities, Printf
include(joinpath(pkgdir(ChemistryLab), "scripts", "snellings2022_kinetics.jl"))

slag, clinker = sn22_degrees("slag"), sn22_degrees("clinker")
temperatures, ratios, days = sort(unique(slag.temperature_C)), sort(unique(slag.w_b)), sort(unique(slag.age))
measured(d, T, w_b, t) = only(d.degree_percent[(d.temperature_C .== T) .& (d.w_b .== w_b) .& (d.age .== t)])
@printf("%d degrees of the slag and %d of the clinker, at %s °C and w/b %s\n",
        length(slag.age), length(clinker.age), join(Int.(temperatures), ", "), join(ratios, ", "))
```

## 2. The clinker under the published law

The four clinker phases of the cement (Table 1, X-ray diffraction) dissolve
under the law of [ParrottKilloh1984](@citet) with the constants and the
activation energies [Lavergne2018](@citet) tabulate, at the Blaine fineness of
the cement and at the temperature of the paste; the degree of the clinker is
their sum weighted by their masses, as the authors weight it. The water/cement
ratios of these pastes, 0.8 to 1.2, are above the range where the law's water
factor acts, so the law gives the same clinker at the three w/b:

```@example slag-temperature
println(" °C   days    law   measured at w/b ", join(ratios, ", "))
for T in temperatures
    law = sn22_clinker_degree(T, first(ratios), days)
    for (k, t) in enumerate(days)
        @printf("%3d %6d  %5.1f   %s\n", T, t, law[k], join((@sprintf("%5.1f", measured(clinker, T, w, t)) for w in ratios), " "))
    end
end
```

**The law has the clinker right only late.** At 5 °C it is within three
points of the measurement at one day and then falls behind it, by 21 to 22
points at seven days; at 20 and 40 °C it is too slow from the first day, by 19 to 21
points at one day. At 20 °C the measured clinker reacts in one day as far as
the law takes it in four. At 90 and 180 days the law is within four points of
the measurement at the three temperatures. The clinker of
[the De Weerdt pastes](@ref ex-ternary-kinetics) was also too slow under this
law at one day, by 7 to 10 points; this one is by 20. The authors find the w/b
of no effect on the clinker before seven days, as the law has it.

## 3. The slag under the Waller law

[Waller1999](@citet) wrote the degree of reaction of a pozzolan as
`α = 1/(1 + (τ/t)ⁿ)` at the reference temperature; the package's
[`waller`](@ref) writes its rate as a function of `α` instead, multiplies it by
the Arrhenius factor `A = exp(−(Ea/R)(1/T − 1/T_ref))`, and may stop it at a
ceiling `α_max`. At a constant temperature the law then integrates to

```math
α(t) = \frac{α_\mathrm{max}}{1 + \left(α_\mathrm{max}\,τ/(A\,t)\right)^{n}},
```

which `sn22_waller_degree` evaluates and the test suite checks against the
rate of `waller` itself. The only constants the package ships for the law are
those of a fly ash, `n = 0.7` and `Ea = 83.14 kJ/mol` ([`WALLER_PARAMS_FLY_ASH`](@ref)),
with no ceiling. The fits below start from them and free the constants named
in the first column (`sn22_slag_fit`), on the eighteen degrees of each w/b or
on the 54 at once:

```@example slag-temperature
fits = [
    ("τ, each w/b", [sn22_slag_fit((w,); shared = (:τ,)) for w in ratios]),
    ("τ, n, Ea, each w/b", [sn22_slag_fit((w,)) for w in ratios]),
    ("τ, n, Ea, α_max, each w/b", [sn22_slag_fit((w,); shared = (:τ, :n, :Ea, :α_max), τ₀ = 10.0) for w in ratios]),
    ("τ; one α_max per w/b", [sn22_slag_fit(; shared = (:τ,), per_wb = (:α_max,), τ₀ = 10.0)]),
    ("τ, n; one α_max per w/b", [sn22_slag_fit(; shared = (:τ, :n), per_wb = (:α_max,), τ₀ = 10.0)]),
    ("τ, n, Ea; one α_max per w/b", [sn22_slag_fit(; per_wb = (:α_max,), τ₀ = 10.0)]),
]
println(rpad("constants fitted", 30), "  misfit, points of degree")
for (label, f) in fits
    println(rpad(label, 30), join((@sprintf("%6.1f", x.rms) for x in f), ""))
end
```

The first three rows fit each w/b alone, one column each, the last three the
54 degrees together with one ceiling per w/b. The fly-ash shape, a sigmoid that goes on to
complete reaction, misses the measurement by 15 to 18 points of degree, more
than the 10 points the measurement is good to: the slag of these pastes slows
down where the sigmoid does not, the more so the hotter the paste. Freeing its
exponent and its activation energy brings the misfit to 4 or 5 points with an
exponent near 0.3, a sigmoid spread over decades. A ceiling does better with
fewer constants. With one ceiling per w/b and the fly-ash exponent and
activation energy kept, the misfit is 3.5 points on the 54 degrees; freeing
the exponent and then the activation energy brings it to 3.1 and 2.6. The
four constants fitted on each w/b alone, twelve in all, do barely better, 2.0
to 2.7 points, than these six. The water acts on the slag through its ceiling
alone, as the authors observe it: the w/b changes the degree of the slag only
at the later ages.

The last fit, its constants and how well the 54 degrees determine them, for a
measurement good to the 10 points the authors give ([`identifiability`](@ref)):

```@example slag-temperature
fit = only(last(fits)[2])
for (name, value, se) in zip(fit.names, fit.values, fit.identifiability.stderr)
    @printf("%-16s %8.4g  ± %2.0f %%\n", name, value, 100 * se)
end
println("singular values: ", round.(fit.identifiability.S; sigdigits = 3))
println()
println(" °C   days   computed/measured at w/b ", join(ratios, ", "))
for T in temperatures, t in days
    i = [findfirst((fit.temperature_C .== T) .& (fit.w_b .== w) .& (fit.age .== t)) for w in ratios]
    @printf("%3d %6d   %s\n", T, t, join((@sprintf("%5.1f/%5.1f", fit.model[j], fit.measured[j]) for j in i), "  "))
end
```

The six constants are all determined, the least well to 33 % (the time), the
activation energy to 17 %. The ceiling rises with the water, from 0.54 at w/b 0.4 to 0.70 at 0.6. The
exponent, 0.59, is near the fly ash's 0.7. The activation energy, 67 kJ/mol, is
the one the authors find, from the same degrees with their own law and with a
time of each temperature (their Table 2):

```@example slag-temperature
ea = literature_table(SN22, "activation_energies"; constituent = "slag")
println("Snellings et al. (2022), Table 2: ", join((@sprintf("%.0f at w/b %.1f", ustrip(us"kJ/mol", e), ustrip(w)) for (e, w) in zip(ea.Ea, ea.w_b)), ", "), " kJ/mol")
@printf("fitted here: %.1f kJ/mol; fly-ash set of the package: %.2f kJ/mol\n", fit.θ[first(ratios)].Ea, ustrip(us"kJ/mol", WALLER_PARAMS_FLY_ASH.Ea))
```

The fly-ash value raises the misfit from 2.6 to 3.1 points: the data prefer
the lower activation energy, though the 10 points of the measurement do not
exclude the higher one.

## 4. The same law on two other slags

The round robin of [Durdzinski2017](@citet) measured the degree of reaction of
two other slags, at 40 % in another cement at w/b 0.4, by image analysis of
electron micrographs in two laboratories (their Table 5). The law of w/b 0.4,
fitted above, against them, nothing refitted (`sn22_durdzinski`). ASSUMED:
the pastes cured at 20 °C, which the round robin does not state.

```@example slag-temperature
dz = sn22_durdzinski(fit.θ[0.4])
println("slag  lab  days  measured  computed")
for i in eachindex(dz.age)
    @printf("%-4s  %-3s %5.0f  %8.0f  %8.1f\n", dz.material[i], dz.lab[i], dz.age[i], dz.measured[i], dz.model[i])
end
for m in ("S1", "S2"), lab in ("B", "E")
    k = (dz.material .== m) .& (dz.lab .== lab)
    @printf("%s by %s: largest difference %.1f points\n", m, lab, maximum(abs, dz.model[k] .- dz.measured[k]))
end
```

As laboratory B measured them, both slags are within 6.2 points of the law at
every age, and the first one as laboratory E measured it within 7.2. The two
laboratories differ by up to 13 points on the same paste, and the round robin
puts the precision of any technique at ±5 points at best. Laboratory E finds
the second slag faster, by 11 points at seven days and 16 at ninety, beyond the
ceiling of 0.54 the law takes from the pastes of [Snellings2022](@citet) at that w/b.

## 5. What the comparison says

**The clinker law has the temperature right only late.** With its published
constants, the law of [ParrottKilloh1984](@citet) brings the clinker of this cement to
within four points of the measured degree at three and six months, at the
three temperatures, but it is too slow at 20 and 40 °C from the first day, and
at 5 °C after it.

**The slag needs a ceiling, and then the Waller law holds.** The sigmoid of
the fly ash cannot follow a slag that slows down by the third month; with a
ceiling that rises with the water, one time, one exponent and one activation
energy describe the three temperatures and the three w/b to 2.6 points of
degree, a quarter of the measurement's uncertainty. The activation energy,
67 kJ/mol, is the authors', found by another route, below the 83 kJ/mol of the
fly-ash set. These are fitted constants, of one slag in one cement; on the two
slags of the round robin they hold within 7.2 points in three series of four.
They are not shipped as the package's constants.

## Where to go next

[CEM II/B-V and CEM II/B-M (V-LL), with their CEM I, integrated in time](@ref ex-ternary-kinetics)
calibrates the alite of another clinker and integrates the pastes with their
hydrates, and [Kinetics under partial equilibrium](@ref sec-theory-pe-kinetics)
writes the formulation those runs integrate.
