# [Carbonation of a CEM I 52.5 N and of its blends with limestone and metakaolin](@id sec-cement-carbonation)

!!! info "Before this page"
    [CO₂ Dissolution and the Carbonate System](@ref sec-co2-carbonate) and
    [Recipes: materials, extents and what has not reacted](@ref man-recipes).

**Carbonation** is the reaction of the carbon dioxide of the air with a hardened
paste. It dissolves in the pore water and takes the calcium of the hydrates to
form calcium carbonate: first that of portlandite,

```math
\mathrm{Ca(OH)_2 + CO_2 \longrightarrow CaCO_3 + H_2O},
```

then that of the C-S-H and of the aluminate hydrates. The pore solution loses
its alkalinity with them, and below a pH of about 9 the steel of a reinforced
concrete is no longer protected.

[Shi2016](@citet) carbonated four mortars of one white Portland cement
(CEM I 52.5 N) after 91 days of hydration, in air with 1 % CO₂:

| mortar | binder |
|:--|:--|
| P | the cement alone |
| L | 68.1 % cement, 31.9 % limestone |
| M | 68.1 % cement, 31.9 % metakaolin |
| ML | 68.1 % cement, 25.5 % metakaolin, 6.4 % limestone |

They measured the degree of hydration of the clinker and of the metakaolin by
NMR, and the portlandite and calcium carbonate by thermogravimetry. They also
computed each paste with GEMS and the Cemdata07 database, adding carbon dioxide
step by step. This page computes the same pastes from their published data with
the Cemdata18 database. It compares the pastes before carbonation with the
measurements, and the carbonation sequence with the authors' calculation.

```@example carbonation
using ChemistryLab
using DynamicQuantities
using JSON
using Printf
using Plots
default(framestyle = :box, grid = false)

# The materials, recipes and system of the four pastes, shared with the test of
# this page and with the generator of the GEMS3K replay.
include(joinpath(pkgdir(ChemistryLab), "scripts", "shi_2016.jl"))
model = cemdata18_activity_model(:KOH)
cs = shi16_system()
nothing # hide
```

## 1. The materials

The cement is described by its phases: alite and belite measured by NMR, the
aluminate by mass balance, gypsum, calcite and free lime. The paper gives no
ferrite phase: the little iron of a white cement is taken as part of the other
phases. Counted with their pure formulas, these phases leave 1.3 % of the
cement for the magnesia, alkalis and sulfate that the oxide analysis gives at
3.9 %. The template keeps those minor oxides at their analyzed amounts and
scales the phases down to make room (`remainder = "analysis"`). The limestone is
a chalk, 93.8 % calcite, and the metakaolin is described by its oxide analysis
alone.

```@example carbonation
for m in (material_template("white Portland cement (Shi 2016)", SHI16_DB),
          material_template("limestone (Shi 2016)", SHI16_DB),
          material_template("metakaolin (Shi 2016)", SHI16_DB))
    println(m.name)
    for c in m.constituents
        @printf("  %-14s %5.1f %%\n", c.name, 100c.mass_fraction)
    end
end
```

## 2. How much has reacted at 91 days

The paper measured the degrees of hydration at 28 and 180 days, not at 91, and
states that little hydration takes place after 91 days: the values at 180 days
are taken. The aluminate, the gypsum and the free lime are taken as reacted.

```@example carbonation
t = SHI16_AGE   # 91 days
println("        alite  belite  metakaolin")
for mix in SHI16_MIXES
    mk = mix in ("M", "ML") ? @sprintf("%6.2f", extent(shi16_extent(mix, "MK"), t)) : "     –"
    @printf("%-4s  %6.2f  %6.2f  %s\n", mix, extent(shi16_extent(mix, "alite"), t), extent(shi16_extent(mix, "belite"), t), mk)
end
```

## 3. The pastes before carbonation

The mortars are 1 part of binder to 3 of sand, at a water/binder ratio of 0.5.
The sand is inert and left out: the pastes are computed per 100 g of binder. The
paper counts in moles per 100 g of **ignited mortar**, the mortar heated to
800 °C, which is the sand and the binder less its loss on ignition.
`shi16_ignited_mortar` gives that mass per 100 g of binder, about 400 g,
and the helper below converts the model's amounts to the paper's unit.

```@example carbonation
pastes = Dict(mix => first(equilibrate_certified(shi16_recipe(mix), cs; t, model)) for mix in SHI16_MIXES)
per_mortar(x, mix) = 100 * x / shi16_ignited_mortar(mix)
portlandite(rs) = ustrip(us"mol", rs.state.n[findfirst(s -> symbol(s) == "Portlandite", cs.species)])
function csh_ca_si(rs)
    gel = solid_solution_totals(rs.state, "CSHQ").elements
    si = get(gel, :Si, 0.0)
    return si > 0 ? gel[:Ca] / si : NaN
end
println("       certified   pH      C-S-H Ca/Si    portlandite (mol/100 g ignited mortar)")
println("                   model   model          model    measured (TGA, core)")
for mix in SHI16_MIXES
    rs = pastes[mix]
    @printf("%-4s   %-9s  %5.2f   %5.2f          %6.3f   %6.3f\n", mix, rs.certificate.optimal, pH(rs.state, model),
            csh_ca_si(rs), per_mortar(portlandite(rs), mix), shi16_measured(mix, "core portlandite"))
end
```

The authors' calculation gave a pH of 13.4 and a Ca/Si of 1.63 for P and L, and
12.8 and 1.29 for the metakaolin pastes.

## 4. The calcium that carbonation can take

Every mole of calcium held by a hydrate can become a mole of calcium carbonate.
Table 5 of the paper counts it hydrate by hydrate, as a **CO₂ binding
capacity**. The calcite already present and the unreacted cement do not count.

```@example carbonation
println("       CO2 binding capacity (mol/100 g ignited mortar)")
println("       this model   the authors' model")
for mix in SHI16_MIXES
    @printf("%-4s   %8.3f   %10.3f\n", mix, per_mortar(shi16_hydrate_cao(pastes[mix])["total"], mix), shi16_computed(mix, "total"))
end
```

## 5. Carbonation step by step

[`carbonate`](@ref) adds carbon dioxide to each paste in steps of 5 g per 100 g
of binder, up to 50 g, the range of the authors' calculation. Each equilibrium
starts from the one before.

```@example carbonation
M_CO2 = ustrip(us"g/mol", Species("CO2")[:M])
grams = shi16_co2_grams()
sweeps = Dict(mix => carbonate(pastes[mix], grams ./ M_CO2) for mix in SHI16_MIXES)
for mix in SHI16_MIXES
    @printf("%-4s %d steps, %d certified\n", mix, length(sweeps[mix]), count(rs -> rs.certificate.optimal, sweeps[mix]))
end
```

```@example carbonation
colors = Dict("P" => :black, "L" => :gray, "M" => :firebrick, "ML" => :darkorange)
p1 = plot(; xlabel = "CO₂ added (g per 100 g of binder)", ylabel = "pH", legend = :bottomleft)
p2 = plot(; xlabel = "CO₂ added (g per 100 g of binder)", ylabel = "Ca/Si of the C-S-H", legend = false)
for mix in SHI16_MIXES
    plot!(p1, grams, [pH(rs.state, model) for rs in sweeps[mix]]; color = colors[mix], label = mix, linewidth = 2)
    plot!(p2, grams, [csh_ca_si(rs) for rs in sweeps[mix]]; color = colors[mix], linewidth = 2)
end
hline!(p1, [9.7]; color = :steelblue, linestyle = :dash, label = "9.7")
fig = plot(p1, p2; layout = (1, 2), size = (900, 380), left_margin = 6Plots.mm, bottom_margin = 6Plots.mm)
savefig(fig, "carbonation-shi2016.svg")
println("pH at    0 g    10 g    20 g    30 g    40 g    50 g of CO2")
for mix in SHI16_MIXES
    @printf("%-4s", mix)
    for g in 0:10:50
        @printf("  %6.2f", pH(sweeps[mix][findfirst(==(g), grams)].state, model))
    end
    println()
end
```

![](carbonation-shi2016.svg)

The authors found the pH to fall in two steps: to about 9.7 once the high-calcium
C-S-H is carbonated, and to 7.4 once the gel is gone, in contact with 1 % CO₂.
The carbon dioxide taken up when the pH first falls below 9.7 is their
**effective** binding capacity. It is read here between the two steps that
bracket that pH, interpolated linearly, so to within a step of 5 g, or
0.03 mol per 100 g of ignited mortar:

```@example carbonation
function below(rs_list, level)
    p = [pH(rs.state, model) for rs in rs_list]
    k = findfirst(<(level), p)
    (k === nothing || k == 1) && return NaN
    return grams[k - 1] + (grams[k] - grams[k - 1]) * (p[k - 1] - level) / (p[k - 1] - p[k])
end
println("       effective capacity (mol CO2/100 g ignited mortar)")
println("       this model (pH < 9.7)   the authors' model")
for mix in SHI16_MIXES
    @printf("%-4s   %10.3f   %16.3f\n", mix, per_mortar(below(sweeps[mix], 9.7) / M_CO2, mix), shi16_computed(mix, "effective"))
end
```

## 6. The same budgets through GEMS3K

`test/reference/xgems_shi2016.json` holds the answer of GEMS3K on the
carbonated budgets of the four pastes, with the phases of this system and the
database Cemdata18. An element at a trace, below 0.01 mmol per kg of water, is
left out of the relative differences, which mean nothing there.

```@example carbonation
fixture = JSON.parsefile(joinpath(pkgdir(ChemistryLab), "test", "reference", "xgems_shi2016.json"))
worst_pH = 0.0
worst = Dict{String, Float64}()
for row in fixture["rows"]
    rs = sweeps[row["mix"]][findfirst(==(row["co2_g"]), grams)]
    ps = pore_solution_mmol(rs.state)
    global worst_pH = max(worst_pH, abs(pH(rs.state, model) - row["pH"]))
    for (e, x) in row["mmol_per_kg_water"]
        x > 0.01 && (worst[e] = max(get(worst, e, 0.0), abs(ps[e] / x - 1)))
    end
end
@printf("GEMS3K converged on %d of %d budgets\n", count(r -> r["converged"], fixture["rows"]), length(fixture["rows"]))
@printf("largest difference in pH: %.4f\n", worst_pH)
for (e, w) in sort(collect(worst))
    @printf("  %-4s largest relative difference: %.2f %%\n", e, 100w)
end
```

## 7. What the comparison says

**Before carbonation, the pastes are those the paper measured.** The portlandite
of the cement alone is within 1 % of the thermogravimetric measurement (0.101
against 0.100 mol per 100 g of ignited mortar), and that of the limestone blend
within 5 % (0.073 against 0.070). In the two metakaolin blends the model has
consumed it all, where the measurement finds a little (0.008 and 0.004): the
paste keeps some portlandite beside a gel it is no longer in equilibrium with.
The C-S-H has the Ca/Si of the authors' calculation in P, L and ML (1.62 against
1.63, and 1.27 against 1.29). In M it holds more calcium, 1.46 against 1.29.
The pH is 0.2 and 0.3 lower than theirs in P and L (13.17 and 13.07 against
13.4), and within 0.3 in the metakaolin blends.

**The calcium there is to carbonate is the same.** The total CO₂ binding
capacity is within 1 % of the authors' in P and L, 2 % in ML, and 7 % below it
in M, the paste whose gel differs.

**The sequence is the same.** While portlandite carbonates, the pH stays above
13: it is set by the sodium and potassium of the pore solution. It then falls
as the C-S-H gives up its calcium, and holds at about 9.7 once the gel is at
its most calcium-poor, as the authors computed. The carbon dioxide taken up
when the pH falls below 9.7 orders the four pastes as theirs does, P, then L,
then the two metakaolin blends together, and exceeds their effective capacity
by 5 % (P) to 17 % (M). The last step does not compare: here every gram of
carbon dioxide stays in the paste, and at 50 g the pore water is acidified by
the excess to a pH near 5, where the authors kept the paste in contact with air
at 1 % CO₂ and reached 7.4.

**A mortar in air is far from that equilibrium.** Near the exposed surface the
thermogravimetry finds 0.12 to 0.13 mol of carbonate per 100 g of ignited mortar
after 91 days in all four mortars, below every capacity. The paper concludes the
same (its Section 4.3). An equilibrium calculation says how much calcium can be
carbonated and at what pH the paste then stands. How fast the carbon dioxide
gets there is a matter of transport, which the paper relates to the porosity
and to the water held in the finest pores.

**The two codes agree.** On the 44 budgets, GEMS3K and ChemistryLab give the
same pH to 0.011 and the same dissolved elements to 7 %, with the same phases
and the same database.

## Where to go next

The four processes of the recipe layer on one paste, carbonation among them, are
on [A CEM I paste replaced by fly ash, carbonated, salted and leached](@ref ex-cement-processes).
