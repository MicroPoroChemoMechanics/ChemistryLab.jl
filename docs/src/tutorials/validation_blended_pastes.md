# [Validation against measured blended pastes: CEM II/B-V and CEM II/B-M (V-LL), with their CEM I](@id ex-validation-blended)

!!! info "Before this page"
    [Validation against a measured paste: a CEM I 42.5 N through its first year](@ref ex-validation),
    whose method this page follows.

[DeWeerdt2011](@citet) blended one clinker, interground with
gypsum, with a siliceous fly ash and a limestone powder, in four pastes at a
water/binder ratio of 0.5 and 20 °C:

| paste | binder | designation |
|:--|:--|:--|
| OPC | clinker and gypsum | CEM I |
| OPC-L | 95 % OPC, 5 % limestone | CEM I, the limestone a minor additional constituent |
| OPC-FA | 65 % OPC, 35 % fly ash | CEM II/B-V |
| OPC-FA-L | 65 % OPC, 30 % fly ash, 5 % limestone | CEM II/B-M (V-LL) |

They measured how much of each clinker phase had reacted (XRD, their Table 7),
how much of the fly ash (image analysis, the fit printed on their Fig. 7), the
portlandite (Table 7) and the pore solution (Table 8), up to 180 days. With the
degrees of reaction taken from the measurement, what is left to compare is the
chemistry: what the reacted material becomes, and what stays in solution.

```@example blended
using ChemistryLab
using DynamicQuantities
using JSON
using Printf
using Plots
default(framestyle = :box, grid = false)

# The materials, recipes and system of the four pastes, shared with the test of
# this page and with the generator of the GEMS3K replay.
include(joinpath(pkgdir(ChemistryLab), "scripts", "de_weerdt_2011.jl"))
model = cemdata18_activity_model(:KOH)
nothing # hide
```

## 1. The materials

The clinker is described by its Rietveld phases (Table 2) and by what its oxide
analysis (Table 1) holds beyond them, its free lime, alkalis, magnesia and
sulfate; the fly ash by its crystals (Table 3), which do not react, and its glass,
found by difference from its analysis; the limestone by its calcite. All three
are templates of `data/recipe_templates.toml`:

```@example blended
for m in (dw11_clinker("OPC"), dw11_fly_ash(), dw11_limestone())
    println(m.name)
    for c in m.constituents
        @printf("  %-14s %5.1f %%\n", c.name, 100c.mass_fraction)
    end
end
```

## 2. How much has reacted

The degree of reaction of each clinker phase is one minus its content over its
content before hydration, both from Table 7; the glass of the fly ash reacts at
the rate of the fit on Fig. 7, the fly ash as a whole reaching
``-15 + 10\ln(t + 4.5)`` percent at ``t`` days. For the CEM II/B-V:

```@example blended
clinker, fly_ash = dw11_clinker("OPC-FA"), dw11_fly_ash()
println("            1 d     28 d    140 d")
for c in vcat(clinker.constituents, [c for c in fly_ash.constituents if c.name == "glass"])
    c.name in ("C3S", "C2S", "C3A", "C4AF", "glass") || continue
    @printf("  %-8s %6.2f  %6.2f  %6.2f\n", c.name, (extent(c.extent, t) for t in (1, 28, 140))...)
end
```

## 3. The pastes at each age

```@example blended
cs = dw11_system()
days = dw11_days()
# `Any`: the type of a state is long enough that a dictionary specialized on it
# takes minutes to compile.
pastes = Dict{String, Any}(mix => hydrate(dw11_recipe(mix), cs, days; model) for mix in DW11_MIXES)
for mix in DW11_MIXES
    @printf("%-9s %d ages, %d certified\n", mix, length(pastes[mix]), count(rs -> rs.certificate.optimal, pastes[mix]))
end
```

## 4. Portlandite

Table 7 gives the portlandite by XRD in weight percent of the dry paste, which
the paper states to ±1 wt.%; the model's is its mass over the mass of every solid
of the paste, the unreacted binder included.

```@example blended
println("            ", join((rpad(m, 16) for m in DW11_MIXES)))
println("            ", join((rpad("model/measured", 16) for _ in DW11_MIXES)))
for (k, d) in enumerate(days)
    measured = [dw11_measured_portlandite(mix, d) for mix in DW11_MIXES]
    all(isnan, measured) && continue
    @printf("%6.0f d    ", d)
    for (mix, m) in zip(DW11_MIXES, measured)
        print(rpad(@sprintf("%.1f/%.1f", dw11_portlandite(pastes[mix][k]), m), 16))
    end
    println()
end
```

## 5. The pore solution

```@example blended
elements = ["pH", "Na", "K", "S", "Si", "Al", "OH-"]
for mix in DW11_MIXES
    println(mix)
    println("   age    ", join((rpad(e, 14) for e in elements)))
    for (k, d) in enumerate(days)
        d in (1.0, 28.0, 140.0) || continue
        rs = pastes[mix][k]
        ps = pore_solution_mmol(rs.state)
        @printf("%6.0f d  ", d)
        for e in elements
            cell = e == "pH" ? @sprintf("%.2f/%.1f", pH(rs.state, model), dw11_measured(mix, d, e)) :
                @sprintf("%.3g/%.3g", ps[e], dw11_measured(mix, d, e))
            print(rpad(cell, 14))
        end
        println()
    end
end
```

```@example blended
colors = Dict("OPC" => :black, "OPC-L" => :gray, "OPC-FA" => :firebrick, "OPC-FA-L" => :darkorange)
p1 = plot(; xlabel = "time (d)", ylabel = "pH", title = "pH", legend = :bottomleft)
p2 = plot(; xlabel = "time (d)", ylabel = "portlandite (wt.%)", title = "Portlandite", legend = false)
for mix in DW11_MIXES
    plot!(p1, days, [pH(rs.state, model) for rs in pastes[mix]]; color = colors[mix], label = mix, linewidth = 2)
    scatter!(p1, days, [dw11_measured(mix, d, "pH") for d in days]; color = colors[mix], label = "")
    plot!(p2, days, [dw11_portlandite(rs) for rs in pastes[mix]]; color = colors[mix], linewidth = 2)
    scatter!(p2, days, [dw11_measured_portlandite(mix, d) for d in days]; color = colors[mix])
end
fig = plot(p1, p2; layout = (1, 2), size = (900, 380), left_margin = 6Plots.mm, bottom_margin = 6Plots.mm)
savefig(fig, "validation-deweerdt.svg"); nothing # hide
```

Lines are the model, markers the measurement.

![](validation-deweerdt.svg)

The C-S-H itself, against the SEM-EDX analyses of the paper (Ca/Si 1.8 ± 0.1
without fly ash; with it, 1.7 at 1 day falling to 1.4 at 140 days, the Al/Si
rising from 0.06 to 0.13):

```@example blended
println("            Ca/Si of the gel         Al/Si")
for mix in DW11_MIXES
    gel = solid_solution_totals(pastes[mix][end].state, "CSHQ").elements
    @printf("  %-9s %.2f at 140 d           %.2f\n", mix, gel[:Ca] / gel[:Si], get(gel, :Al, 0.0) / gel[:Si])
end
```

## 6. The same budgets through GEMS3K

`test/reference/xgems_deweerdt2011.json` holds the answer of GEMS3K on the twenty
budgets, with the phases of this system; the AFm sulfate and hydroxide are, in
both codes, the Guggenheim binary of Cemdata18, one composition unless the
certificate asks for two.

```@example blended
fixture = JSON.parsefile(joinpath(pkgdir(ChemistryLab), "test", "reference", "xgems_deweerdt2011.json"))
worst = Dict{String, Float64}()
worst_pH = 0.0
for row in fixture["rows"]
    rs = pastes[row["mix"]][findfirst(==(row["time_d"]), days)]
    ps = pore_solution_mmol(rs.state)
    for (e, x) in row["mmol_per_kg_water"]
        worst[e] = max(get(worst, e, 0.0), abs(ps[e] / x - 1))
    end
    global worst_pH = max(worst_pH, abs(pH(rs.state, model) - row["pH"]))
end
@printf("GEMS3K converged on all %d budgets: %s\n", length(fixture["rows"]), all(r -> r["converged"], fixture["rows"]))
@printf("largest difference in pH: %.4f\n", worst_pH)
for (e, w) in sort(collect(worst))
    @printf("  %-4s largest relative difference: %.2f %%\n", e, 100w)
end
```

## [7. A C-S-H that takes aluminum](@id sec-validation-blended-gels)

The package ships two other models of the gel: `CNASH_ss`
[Myers2014](@cite), mixed on the sublattices of its authors, which takes
aluminum and the alkalis, and `CASH+NK` [Miron2022a](@cite), whose end members
carry the alkalis but no aluminum. The same budgets, with each in place of
`CSHQ`, the rest of the phase list unchanged, and for `CASH+NK` without the four
aqueous ion pairs its authors left out when fitting it, NaOH⁰, KOH⁰, NaHSiO₃⁰
and KHSiO₃⁰ [Miron2022a; Section 3.2](@cite):

```@example blended
cashplus = build_species(datapath("cemdata18-cashplus.json"); verbose = false)
gel_systems = (
    "CNASH_ss" => phase_list_system(DW11_PHASES, DW11_SUBSTANCES; replace = Dict("CSHQ" => "CNASH_ss")),
    "CASH+NK" => phase_list_system(DW11_PHASES, cashplus; replace = Dict("CSHQ" => "CASH+NK"), exclude_aqueous = DW11_CASHPLUS_EXCLUDED),
)
other = Dict{String, Any}("CSHQ" => pastes)
for (gel, cs_g) in gel_systems
    other[gel] = Dict{String, Any}(mix => hydrate(dw11_recipe(mix), cs_g, days; model) for mix in DW11_MIXES)
end
k90 = findfirst(==(90.0), days)
ettringite(rs) = (m = phase_masses(rs); 100 * get(m, "ettringite", 0.0) / sum(values(m)))
println("at 90 days      portlandite, wt.%         ettringite, wt.%        gel at 140 days, Ca/Si and Al/Si")
for mix in DW11_MIXES
    println(mix, ", measured: portlandite ", dw11_measured_portlandite(mix, 90), ", ettringite ", something(dw11_phase_content(mix, 90, "ettringite"), NaN))
    for gel in ("CSHQ", "CNASH_ss", "CASH+NK")
        rs = other[gel][mix]
        e = solid_solution_totals(rs[end].state, gel).elements
        @printf("  %-9s %16.1f %24.1f %22.2f %6.3f\n", gel, dw11_portlandite(rs[k90]), ettringite(rs[k90]),
                e[:Ca] / e[:Si], get(e, :Al, 0.0) / e[:Si])
    end
end
```

Neither reproduces the four pastes. `CASH+NK` changes nothing that matters
here: it takes no aluminum, its gel is within 0.04 of the Ca/Si of that of
`CSHQ`, and its portlandite within 1.6 points. `CNASH_ss` takes aluminum, an Al/Si of 0.10 to 0.11 against
the 0.13 measured, but its gel sits at a Ca/Si of 1.16 in every paste, the plain
cement included; its end members reach 1.5 at most (`T2C-CNASHss`), where the
paper measures 1.8 without fly ash. The calcium the gel does not take goes to
portlandite: 31.7 against 21.8 wt.% at 90 days in the CEM I. With fly ash the
same excess brings the portlandite near the measurement, 14.7 against 12.5 and
15.7 against 12.2, and the ettringite of the CEM II/B-V is still lost, 0.7
against 6.6 wt.%.


## 8. What the comparison says

**The two codes agree.** On the twenty budgets the pore solutions differ by less
than 0.001 in pH and by at most 1.6 % on an element. What separates the model from
the pastes is the model's, and it would be found with either code.

**Without fly ash, the portlandite is reproduced** from the seventh day, within
1.5 wt.%, and within the ±1 wt.% of the measurement in the CEM I (19.0 against
18.9 at 7 days, 21.4 against 21.8 at 90). At one day the model holds 3.4 and
4.3 wt.% more than the two pastes: the clinker is dissolved as measured, and an equilibrium puts
its calcium hydroxide in portlandite at once.

**With fly ash, the model consumes the portlandite far faster than the pastes**:
4.5 against 12.5 wt.% at 90 days in the CEM II/B-V, 6.7 against 12.2 in the
CEM II/B-M (V-LL). The fly ash dissolves at the measured rate in both, so the
difference is what its silica becomes. The model's C-S-H stays at a Ca/Si of
1.58, the composition in equilibrium with portlandite, in every paste; the paper
measures 1.4 in the fly ash pastes at 140 days, with an Al/Si of 0.13, a gel
poorer in calcium beside a portlandite it has not reached equilibrium with, and
one that takes less calcium per silicon. `CSHQ` has no aluminum end-member, and
the model's gel holds none.

**The pore solution follows from both.** The pH is within 0.3 unit at every age.
In the CEM I the alkalis are low, as on the page of Lothenbach and Winnefeld
(sodium 114 against 302 mmol at 140 days, potassium 318 against 565), `CSHQ`
holding them. In the fly ash pastes the measured pH falls after 28 days, 13.6 to
13.4 and 13.5 to 13.3, while the model's rises, 13.55 to 13.58 and 13.57 to 13.60:
the pastes lose potassium after 28 days (271 to 227 mmol in the CEM II/B-V) and
the model gains it (242 to 268) from the clinker and the glass that go on
releasing it. The aluminum the gel does not take stays in solution, 1.07 against
0.27 mmol in the CEM II/B-V at 140 days, where the sulfate is a sixtieth of the
measurement (0.042 against 2.6); silicon is 4 to 10 times low throughout.

For a user, the line falls here. On the pastes without fly ash the model gives
the portlandite and the pH, with `CSHQ` or `CASH+NK`. With a siliceous fly ash,
`CSHQ` overstates the pozzolanic consumption of portlandite and cannot place the
aluminum, and `CNASH_ss`, which takes the aluminum, does so at a Ca/Si a
Portland cement's gel does not have ([Section 7](@ref sec-validation-blended-gels)):
the two bracket the measured portlandite, from below and from above, and
neither keeps the ettringite of the paste without limestone. What the pastes
need is a gel that takes aluminum at the Ca/Si of a Portland cement, which none
of the three is.

## Where to go next

The method, and a Portland cement followed from its first minute, are on
[Validation against a measured paste](@ref ex-validation).
