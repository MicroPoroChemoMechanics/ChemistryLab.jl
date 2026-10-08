# [CEM I 52.5 R with slag, fly ash and limestone, integrated in time](@id ex-quaternary-kinetics)

!!! info "Before this page"
    [CEM I 52.5 N and slag pastes, integrated in time](@ref ex-blended-slag-kinetics),
    the first blended cement on the kinetic path, and
    [Recipes](@ref man-recipes) for the materials and their templates.

[Scholer2015](@citet) replaced half of a CEM I 52.5 R by blast-furnace slag,
siliceous fly ash and limestone powder in ten proportions (their Table 4), the
SO₃ of every mix brought to 3 % with anhydrite, and followed the pastes at 20 °C
by thermogravimetry from one day to six months: the bound water and the
portlandite of their Table 8. Their own calculations take the slag and the fly
ash at the degrees of reaction reported after a year in the literature they
cite, as long-term states. Here the ten pastes are integrated from the mixing,
the four clinker phases under the Parrott–Killoh law [ParrottKilloh1984](@cite),
the glass of each addition under the Waller law [Waller1999](@cite), the
limestone and everything else at equilibrium, and compared with Table 8, which
nothing below was fitted to.

## 1. The materials

The four materials come from their templates (`data/recipe_templates.toml`),
built from the analyses and the Rietveld phases of
[Scholer2015; Tables 1 and 2](@cite): the cement by its phases, the polymorphs
of C₂S and of C₃A each summed into one constituent; the slag and the fly ash as
their crystals and their glass, found by difference; the limestone by its
analysis. As in the authors' calculations, only the glass of the two additions
reacts, and [`glass_species`](@ref) gives each a formula, which
[`with_species`](@ref) puts in its place.

```@example quaternary
using ChemistryLab, DynamicQuantities, OptimaSolver, Printf, Plots
using Logging # hide
include(joinpath(pkgdir(ChemistryLab), "scripts", "scholer2015_kinetics.jl"))

setup = s15_setup()
cs, mats = setup.cs, setup.mats
for m in (mats.opc, mats.bfs, mats.fa)
    println(m.name, ": ", join((@sprintf("%s %.3f", c.name, c.mass_fraction) for c in m.constituents), ", "))
end
```

## 2. The time of each glass

No parameter set of the Waller law is published for these materials. Its
exponent is the one [Lavergne2018](@citet) fitted for a fly ash, assumed for
both glasses, and the characteristic time of each is calibrated on the degree
the article assumes for it after a year, 71.1 % of the slag glass and 43.6 % of
the fly-ash glass; Table 8 stops at six months.

```@example quaternary
for (which, key) in ((:bfs, "assumed_reaction_slag_glass"), (:fa, "assumed_reaction_fly_ash_glass"))
    @printf("%s glass: τ = %.0f days, %.1f %% after one year\n", which,
            ustrip(us"d", s15_glass_time(which)), literature_value("Scholer2015", key))
end
```

## 3. The ten pastes

Each paste is integrated over the six months of [Scholer2015; Table 8](@cite),
at 20 °C, in the activity model the authors used. Neither law reads the
equilibrium partition, so the trajectory does not depend on it: the run solves
it at each accepted step only to report it, with the interior point, and
everything below is computed on the certified replay of each run
([`speciated_states`](@ref)).

```@example quaternary
day = 86400.0
mixes = literature_table("Scholer2015", "mixes").mix
diagnostics = IOBuffer() # hide
runs = with_logger(ConsoleLogger(diagnostics)) do # hide
# `Any`: the type of a run is long enough that a dictionary specialized on it
# takes minutes to compile.
runs = Dict{String, Any}(name => s15_run(cs, mats, name) for name in mixes)
end # hide
occursin("re-speciation failed", String(take!(diagnostics))) && error("a re-speciation failed") # hide
for name in mixes
    r = runs[name]
    @printf("%-9s %s, %3d steps, anhydrite %.2f g\n", name, r.sol.retcode, length(r.sol.t), s15_anhydrite(r.mix))
end
```

## 4. Bound water and portlandite, against Table 8

What the thermobalance weighs is computed on the replayed states: the water of
every hydrate but portlandite, and the portlandite, in percent of the sample
dried at 500 °C, every hydrate taken to have lost all its water there
(`s15_tga`).

```@example quaternary
tga = with_logger(ConsoleLogger(diagnostics)) do # hide
tga = Dict(name => s15_tga(runs[name], s15_measured(name).days) for name in mixes)
end # hide
occursin("could not be certified", String(take!(diagnostics))) && error("an instant was not certified") # hide
ages = s15_measured(first(mixes)).days
mean(x) = sum(x) / length(x)
println(" days   bound water: computed  measured   portlandite: computed  measured   (mean over the ten pastes)")
for (k, d) in enumerate(ages)
    bw = [tga[n][k].bound_water for n in mixes]
    bwm = [s15_measured(n).bound_water[k] for n in mixes]
    ch = [tga[n][k].portlandite for n in mixes]
    chm = [s15_measured(n).portlandite[k] for n in mixes]
    @printf("%5.0f   %20.1f %9.1f   %22.1f %9.1f\n", d, mean(bw), mean(bwm), mean(ch), mean(chm))
end
```

Paste by paste, the two extremes of the substitution of fly ash by limestone at
each level of slag:

```@example quaternary
for name in ("20-30-0", "20-10-20", "30-20-0", "30-0-20")
    m = s15_measured(name)
    println(name, "   days  BW computed / measured   CH computed / measured")
    for (k, d) in enumerate(m.days)
        @printf("          %4.0f   %7.1f / %4.1f          %7.1f / %4.1f\n", d,
                tga[name][k].bound_water, m.bound_water[k], tga[name][k].portlandite, m.portlandite[k])
    end
end
```

The same four pastes in time, the lines computed at the ages of the
measurements:

```@example quaternary
four = ("20-30-0", "20-10-20", "30-20-0", "30-0-20")
function versus_age(field, ylabel, title)
    p = plot(; xscale = :log10, xlabel = "age (days)", ylabel, title, titlefontsize = 10, legend = false)
    for (k, name) in enumerate(four)
        plot!(p, ages, [getproperty(tga[name][j], field) for j in eachindex(ages)]; lw = 2, color = k,
              marker = :circle, ms = 3, label = name * ", computed")
        scatter!(p, ages, getproperty(s15_measured(name), field); color = k, marker = :diamond, label = name * ", measured")
    end
    return p
end
plot(versus_age(:bound_water, "% of the dry sample", "Bound water, 20 °C"),
     versus_age(:portlandite, "% of the dry sample", "Portlandite, 20 °C"),
     plot(fill(NaN, 1, 2length(four)); framestyle = :none, legend = :left,
          label = permutedims(vcat([[n * ", computed", n * ", measured"] for n in four]...)),
          color = permutedims(repeat(1:length(four); inner = 2)), seriestype = permutedims(repeat([:path, :scatter], length(four))),
          marker = permutedims(repeat([:circle, :diamond], length(four))), lw = 2);
     layout = @layout([a b c{0.2w}]), size = (1000, 400), left_margin = 5Plots.mm, bottom_margin = 6Plots.mm)
```

[Scholer2015](@citet) computed the two pastes without limestone at the degrees
of reaction they assume for the long term (their Table 7, in percent of the dry
hydrates rather than of the dry sample):

```@example quaternary
t7 = literature_table("Scholer2015", "modeled_long_term")
println("paste     portlandite: this page, 182 days   the authors' calculation   measured, 182 days")
for (k, name) in enumerate(t7.mix)
    @printf("%-9s %30.1f %26.1f %20.1f\n", name, tga[name][end].portlandite, t7.portlandite[k],
            s15_measured(name).portlandite[end])
end
```

And where the aluminum goes, in the paste richest in limestone after six months
(mol per 100 g of binder):

```@example quaternary
r = runs["20-10-20"]
st = with_logger(ConsoleLogger(diagnostics)) do # hide
st = only(speciated_states(r.sol, r.kp; times = [182day]))
end # hide
shown_phases = ("Cal", "ettringite", "monocarbonate", "hemicarbonate", "monosulphate12", "C3AFS0.84H4.32", "C3FS0.84H4.32")
amounts = [ustrip(us"mol", st.n[findfirst(x -> symbol(x) == s, r.kp.system.species)]) for s in shown_phases]
for (s, a) in zip(shown_phases, amounts)
    @printf("%-16s %.4f\n", s, a)
end
# The aluminum carriers, in mmol; the calcite, which holds none, is left out.
p_al = bar(collect(shown_phases[2:end]), 1000amounts[2:end]; legend = false, color = :steelblue, xrotation = 30,
    ylabel = "mmol per 100 g of binder", title = "20-10-20 after 182 days, 20 °C",
    size = (760, 400), left_margin = 6Plots.mm, bottom_margin = 12Plots.mm)
plot!(p_al; ylims = (0, 1.1 * ylims(p_al)[2]))
```

## 5. What the comparison says

**The bound water follows the measurements.** Its mean over the ten pastes is
within one point of the measured one at one day, 28 days and 91 days, about two
points low at two and seven days, and three points high at six months, where the
computed water keeps rising with the glasses and the measured one stops. Nothing
was fitted to it: the clinker law has its published parameters, and the time of
each glass is set on the degree [Scholer2015](@citet) assume after a year. The
direction of the gap at six months is the one the preparation of the samples
gives: the thermobalance weighs a sample dried at 40 °C after a solvent
exchange, which has already lost part of the water of the C–S–H and of the AFm
phases, and the computation counts all of it.

**The portlandite is where the computation and the measurement part.** Computed,
it is below the measurement from the first day, and its mean over the pastes
reaches its maximum at seven days and then declines, consumed by the pozzolanic reaction of the two glasses,
the more so the more fly ash the paste holds: the paste 20-30-0 has none left
at six months, and 30-0-20, without fly ash, keeps the most. Measured, it stays
near 12 % from two days on. The authors' own calculation, at their long-term
degrees, leaves as little in the two pastes without limestone. The gap is
therefore that of the equilibrium model, with these degrees of reaction, rather
than that of the kinetics: an equilibrium in which the C–S–H the glasses form
takes its calcium from portlandite cannot keep the portlandite, and neither
calculation describes what keeps it in these pastes. Section 6 runs the two
other gels the package ships.

**The limestone stays calcite.** In Cemdata18 [Lothenbach2019](@cite), used
here, the aluminum the clinker and the glasses release goes to ettringite and to
the siliceous hydrogarnet C₃(A,F)S₀.₈₄H₄.₃₂, and little of it to monocarbonate:
in the paste richest in limestone, a twentieth of the amount of hydrogarnet
after six months, and no hemicarbonate. [Scholer2015](@citet) calculated
hemicarbonate and monocarbonate in the presence of limestone, with the database
of the time (their Fig. 2), and found both by X-ray diffraction after six months
(their Fig. 4). The difference is in which phases the equilibrium may form,
which this page leaves as Cemdata18 does.

## 6. The two other gels

The paste with the most fly ash and the one without, with `CNASH_ss`
[Myers2014](@cite) and then `CASH+NK` [Miron2022a](@cite) in place of `CSHQ`,
everything else as in Section 3 (`s15_setup(; gel)`):

```@example quaternary
gel_runs = Dict{Tuple{String, String}, Any}()
cshq_last = Dict{String, Any}()
with_logger(ConsoleLogger(diagnostics)) do # hide
for name in ("20-30-0", "30-0-20")
    cshq_last[name] = only(speciated_states(runs[name].sol, runs[name].kp; times = [ages[end] * 86400.0]))
end
for gel in ("CNASH_ss", "CASH+NK")
    s = s15_setup(; gel)
    for name in ("20-30-0", "30-0-20")
        r = s15_run(s.cs, s.mats, name)
        gel_runs[(gel, name)] = (; tga = s15_tga(r, ages), last = only(speciated_states(r.sol, r.kp; times = [ages[end] * 86400.0])))
    end
end
end # hide
gel_log = String(take!(diagnostics)) # hide
occursin("re-speciation failed", gel_log) && error("a re-speciation failed") # hide
occursin("could not be certified", gel_log) && error("an instant was not certified") # hide
for name in ("20-30-0", "30-0-20")
    println(name, ", portlandite, % of the dry sample")
    println("  days      CSHQ  CNASH_ss  CASH+NK  measured")
    for (k, d) in enumerate(ages)
        @printf("%6.0f  %8.1f %8.1f %8.1f %9.1f\n", d, tga[name][k].portlandite,
                gel_runs[("CNASH_ss", name)].tga[k].portlandite, gel_runs[("CASH+NK", name)].tga[k].portlandite,
                s15_measured(name).portlandite[k])
    end
    ratio(st, gel) = (e = solid_solution_totals(st, gel).elements; @sprintf("Ca/Si %.2f, Al/Si %.3f", e[:Ca] / e[:Si], get(e, :Al, 0.0) / e[:Si]))
    println("  the gel at ", Int(ages[end]), " days: CSHQ ", ratio(cshq_last[name], "CSHQ"), "; CNASH_ss ",
            ratio(gel_runs[("CNASH_ss", name)].last, "CNASH_ss"), "; CASH+NK ", ratio(gel_runs[("CASH+NK", name)].last, "CASH+NK"))
end
```

```@example quaternary
panels = map(("20-30-0", "30-0-20")) do name
    p = plot(; xscale = :log10, xlabel = "age (days)", ylabel = "portlandite, % of the dry sample",
             title = name, titlefontsize = 10, legend = name == "20-30-0" ? :bottomleft : false)
    plot!(p, ages, [x.portlandite for x in tga[name]]; lw = 2, marker = :circle, ms = 3, label = "CSHQ")
    for gel in ("CNASH_ss", "CASH+NK")
        plot!(p, ages, [x.portlandite for x in gel_runs[(gel, name)].tga]; lw = 2, marker = :circle, ms = 3, label = gel)
    end
    scatter!(p, ages, s15_measured(name).portlandite; color = :black, marker = :diamond, label = "measured")
end
plot(panels...; layout = (1, 2), size = (900, 380), left_margin = 6Plots.mm, bottom_margin = 6Plots.mm)
```

`CNASH_ss` keeps the portlandite near the measurement in the paste with the
most fly ash: within 0.3 points at one, two and 182 days, 1.2 to 3.0 points high
between. Without fly ash it is 1.4 and 1.6 points low at one and two days, then
0.6 to 3.0 points high from seven. It does so with a gel at a Ca/Si of 1.17
holding aluminum, an Al/Si of 0.085 and 0.099, poorer in calcium than a gel
beside portlandite: the calcium it does not take is the portlandite that stays.
`CASH+NK` gives the portlandite of `CSHQ` within 0.3 points. The measurement
lies between the two models, as in the fly-ash pastes of [DeWeerdt2011](@citet)
([integrated in time here](@ref ex-ternary-kinetics)), whose gel was measured at
a Ca/Si of 1.4: neither gel reproduces the gel and the portlandite together,
`CNASH_ss` is the closer on the portlandite of these blends for a gel too poor
in calcium, and `CSHQ` stays the gel of this page.

## Where to go next

[Kinetics under partial equilibrium](@ref sec-theory-pe-kinetics) writes the
formulation these runs integrate, and
[CEM I 52.5 N and slag pastes, integrated in time](@ref ex-blended-slag-kinetics)
the same path on a slag cement, with its heat.
