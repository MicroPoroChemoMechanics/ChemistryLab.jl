# [C-S-H in water and in alkaline solutions: the CASH+ model](@id ex-cashplus-csh)

!!! info "Before this page"
    [Solid solutions](@ref sec-theory-solid-solutions), sections 7 and 8, where
    the site model and the compound energy formalism are derived.

CASH+ is a model of the C-S-H gel by
[Kulik2022](@citet), which [Miron2022a, Miron2022b](@citet) extended to sodium and
potassium. Like the CNASH gel of
[the site-mixing page](@ref ex-sublattice-csh), it mixes on the structural sites
of the silicate chain and of the interlayer. It also adds two terms: the energy
of the reciprocal reactions between its end-members, and interactions between
the species of one site. ChemistryLab writes it as a [`CompoundEnergyModel`](@ref).

This page does five things. It computes the gel in water, from the Ca/Si at
which amorphous silica stops forming to the one at which portlandite starts, and
compares both ends with [Kulik2022](@citet). It then adds sodium and potassium.
It checks the model against the calculation of 110 gel compositions by
[Miron2022a](@citet). It computes the pore solution of a hydrating Portland
cement with CASH+NK and with CSHQ, against its analysis. Finally it extends the
interlayer to the other alkali and alkaline-earth metals.

The model has two mixing sites. The **bridging tetrahedron** (BT) of the silicate
chain holds a silicate `S`, a vacancy `v` or a calcium `C`. The **interlayer
cation** (IC) holds a vacancy `v`, a calcium `C` and, with the alkalis, a sodium
`N` or a potassium `K`. An end-member is named by what it puts on each: TSCh has
a silicate on the bridging tetrahedron and a calcium in the interlayer (`T` and
`h` stand for the fixed dimeric unit and interlayer water).

```@example cashplus
using ChemistryLab
using DynamicQuantities
using OptimaSolver
using Printf
using Plots
default(framestyle = :box, grid = false)

# Cemdata18 with the thirty-three end-members of CASH+ and its extensions, the
# CaSiO3@ complex the model was fitted with, and the aqueous species the
# extensions add (built on first use).
substances = build_species(datapath("cemdata18-cashplus.json"); verbose = false)
byname = Dict(symbol(s) => s for s in substances)
phases = Dict(ChemistryLab.name(p) => p for p in build_solid_solutions(datapath("solid_solutions.toml"), byname))

# The activity model the paper fits with: extended Debye-Hückel with the ion size
# and B-dot Cemdata18 gives for KOH solutions.
model = cemdata18_activity_model(:KOH)

# The gel and the aqueous species of its elements, with portlandite and
# amorphous silica, which may form beside it. The neutral complexes of the
# alkalis (NaOH@, KOH@ ...) are left out, as Miron et al. left them out when they
# fitted the alkali end-members; the Ca(OH)2@ complex they derived for that fit
# is kept with the alkali model and left out of the core model, which Kulik et
# al. fitted without it.
function gel_system(gel)
    members = symbol.(phases[gel].end_members)
    sp = speciation(substances, vcat(["Portlandite", "Amor-Sl"], members); aggregate_state = [AS_AQUEOUS],
                    exclude_species = gel == "CASH+" ? ["Ca(OH)2@"] : String[])
    alkali_complex(s) = charge(s) == 0 && any(el -> haskey(atoms(s), el), (:Na, :K)) &&
        aggregate_state(s) == AS_AQUEOUS
    return ChemicalSystem(filter(!alkali_complex, sp), CEMDATA_PRIMARIES; solid_solutions = [phases[gel]])
end

# A kilogram of water with 0.05 mol of silica and lime at the Ca/Si given.
function paste(cs, ca_si)
    st = ChemicalState(cs)
    set_quantity!(st, "H2O@", 1.0u"kg")
    set_quantity!(st, "Amor-Sl", 0.05u"mol")
    set_quantity!(st, "Portlandite", ca_si * 0.05u"mol")
    b = budget(st)
    return equilibrate_certified(st; model, b)
end

# Moles of an element in solution per kilogram of water, and the amount of a phase.
function in_solution(eq, element)
    cs = eq.system
    n = ustrip.(us"mol", eq.n)
    iw = only(cs.idx_solvent)
    kg_water = n[iw] * ustrip(us"kg/mol", cs.species[iw][:M])
    return sum(get(atoms(cs.species[i]), element, 0) * n[i] for i in cs.idx_solutes) / kg_water
end
amount(eq, sp) = ustrip(us"mol", eq.n[findfirst(s -> symbol(s) == sp, eq.system.species)])

# The site fractions of the gel, from the amounts of its end-members.
function sites(eq, gel)
    p = phases[gel]
    x = [amount(eq, symbol(m)) for m in p.end_members]
    return site_fractions(ChemistryLab.model(p), x ./ sum(x))
end
nothing # hide
```

## 1. The C-S-H in water

Each paste is 0.05 mol of silica and some lime in a kilogram of water, from a
Ca/Si of 0.6 to 2.4. This is how the solubility of a C-S-H is measured: the gel is
equilibrated in water and the solution analyzed. The core model is enough here,
since there is no alkali.

```@example cashplus
cs = gel_system("CASH+")
lattice = ChemistryLab.model(phases["CASH+"]).lattice
ratios = 0.6:0.1:2.4
rows = []
println("Ca/Si   certified  silica  portl.  gel Ca/Si  Ca (mmol/kg)  Si (mmol/kg)   pH")
for r in ratios
    eq, cert = paste(cs, r)
    gel = solid_solution_totals(eq, "CASH+").elements
    y = sites(eq, "CASH+")
    push!(rows, (; r, gel = gel[:Ca] / gel[:Si], ca = 1000in_solution(eq, :Ca), si = 1000in_solution(eq, :Si),
                 pH = pH(eq, model), y, silica = amount(eq, "Amor-Sl"), portlandite = amount(eq, "Portlandite")))
    row = rows[end]
    @printf("%4.2f    %-9s  %-6s  %-6s  %7.3f    %10.3f    %10.4f    %6.3f\n", r, cert.optimal,
            row.silica > 0 ? "yes" : "", row.portlandite > 0 ? "yes" : "", row.gel, row.ca, row.si, row.pH)
end
```

Every paste certifies. Amorphous silica forms beside the gel at the lowest Ca/Si
and portlandite at the highest; between the two, the gel is the only solid. At
each end, three phases share three components (CaO, SiO₂ and H₂O) at a fixed
temperature and pressure, so the composition of every phase is fixed: adding
lime beyond the portlandite point, or silica below the silica point, only makes
more of the solid that forms. [Kulik2022](@citet) give both ends:

```@example cashplus
q(name) = literature_value("Kulik2022", name)
at_silica = rows[findlast(r -> r.silica > 0, rows)]
at_portlandite = rows[findfirst(r -> r.portlandite > 0, rows)]
@printf("gel Ca/Si beside amorphous silica: %.3f   (paper %.2f)\n", at_silica.gel, q("silica_saturation_Ca_Si"))
@printf("gel Ca/Si beside portlandite:      %.3f   (paper %.2f)\n", at_portlandite.gel, q("portlandite_saturation_Ca_Si"))
printed = literature_table("Kulik2022", "portlandite_saturated_site_fractions")
println("\nsite fractions beside portlandite:")
for (site, sp, f) in zip(printed.site, printed.species, printed.fraction)
    s = findfirst(==(site), lattice.sites)
    @printf("  %s %s  %.3f   (paper %.2f)\n", site, sp, at_portlandite.y[s][findfirst(==(sp), lattice.species[s])], f)
end
```

Both ends agree with [Kulik2022](@citet) to the digits they print. The figure
shows the solution and the sites between them:

```@example cashplus
gel_ratio = [r.gel for r in rows]
p1 = plot(gel_ratio, [r.si for r in rows]; label = "Si", yscale = :log10, marker = :circle,
          xlabel = "Ca/Si of the gel", ylabel = "in solution (mmol/kg)", legend = :right)
plot!(p1, gel_ratio, [r.ca for r in rows]; label = "Ca", marker = :square)
p2 = plot(; xlabel = "Ca/Si of the gel", ylabel = "site fraction", legend = :right, ylims = (0, 1))
for (s, style) in ((1, :solid), (2, :dash)), (i, sp) in enumerate(lattice.species[s])
    plot!(p2, gel_ratio, [r.y[s][i] for r in rows]; label = "$(lattice.sites[s]) $sp", linestyle = style, linewidth = 2)
end
fig = plot(p1, p2; layout = (1, 2), size = (900, 380), left_margin = 6Plots.mm, bottom_margin = 7Plots.mm)
savefig(fig, "cashplus-water.svg"); nothing # hide
```

![](cashplus-water.svg)

As the gel takes up calcium, the silicate leaves the bridging tetrahedra, first
for vacancies and then, from a Ca/Si of about 1, for calcium. The interlayer
fills with calcium up to a Ca/Si of about 1.1, then gives a little of it back to
vacancies. The silicon in solution falls by more than two orders of magnitude,
from 4.18 to 0.0097 mmol/kg, and the calcium rises to 20.4 mmol/kg, that of the
solution in equilibrium with portlandite.

### The same gel at 50 and 90 °C

[Kulik2022](@citet) fitted the model at 25 °C and gave its end-members and the CaSiO₃⁰
complex heat capacities for use up to 100 °C (their Section 3.5). They describe
what the model then predicts: the calcium and silicon in solution change
little with temperature, the pH falls by 2 to 2.5 units over a hundred degrees,
and beside portlandite the silicon rises a little. The pastes of Ca/Si 1.2 and
2.4 above, at three temperatures:

```@example cashplus
function paste_at(cs, ca_si, T)
    st = ChemicalState(cs; T = T)
    set_quantity!(st, "H2O@", 1.0u"kg")
    set_quantity!(st, "Amor-Sl", 0.05u"mol")
    set_quantity!(st, "Portlandite", ca_si * 0.05u"mol")
    b = budget(st)
    return equilibrate_certified(st; model, b)
end
println(" Ca/Si    T (°C)   certified   Ca (mmol/kg)   Si (mmol/kg)     pH")
for r in (1.2, 2.4), Tc in (25, 50, 90)
    eq, cert = paste_at(cs, r, (Tc + 273.15)u"K")
    @printf("%5.1f    %5d      %-9s  %10.3f     %10.4f     %6.3f\n", r, Tc, cert.optimal,
            1000in_solution(eq, :Ca), 1000in_solution(eq, :Si), pH(eq, model))
end
```

Every paste certifies. From 25 to 90 °C the pH falls by 1.8 units at a Ca/Si of
1.2 and by 1.8 beside portlandite, about 2.8 units per hundred degrees, somewhat
more than the 2 to 2.5 [Kulik2022](@citet) state. Beside portlandite the silicon
rises, from 0.0097 to 0.0137 mmol/kg, as the paper states, while the calcium
falls from 20.4 to 12.7 mmol/kg, with the solubility of portlandite.

## 2. Sodium and potassium: CASH+NK

`CASH+NK` puts sodium and potassium in the interlayer, as NaHOH⁺ and KHOH⁺. With
twelve end-members, it is the one to use for a cement paste. The pastes below
have a Ca/Si of 1 and of 1.6, in a solution of sodium and potassium hydroxides
with three times as much potassium as sodium, as in the pore solution of a
Portland cement. Potassium dominates, so the Debye-Hückel parameters are those
Cemdata18 [Lothenbach2019](@cite) gives for KOH:

```@example cashplus
cs_nk = gel_system("CASH+NK")
println("Ca/Si  Na, K added (mol/kg)  certified   Na, K in solution (mmol/kg)   gel Na/Si, K/Si   pH")
uptake = Dict()
for r in (1.0, 1.6), na in (0.01, 0.02, 0.05, 0.1)
    st = ChemicalState(cs_nk)
    set_quantity!(st, "H2O@", 1.0u"kg")
    set_quantity!(st, "Amor-Sl", 0.05u"mol")
    set_quantity!(st, "Portlandite", r * 0.05u"mol")
    set_quantity!(st, "Na+", na * u"mol")
    set_quantity!(st, "K+", 3na * u"mol")
    set_quantity!(st, "OH-", 4na * u"mol")
    b = budget(st)
    eq, cert = equilibrate_certified(st; model, b)
    gel = solid_solution_totals(eq, "CASH+NK").elements
    point = (na = 1000in_solution(eq, :Na), k = 1000in_solution(eq, :K), gna = gel[:Na] / gel[:Si], gk = gel[:K] / gel[:Si])
    push!(get!(uptake, r, []), point)
    @printf("%4.1f   %5.2f, %5.2f          %-9s  %8.2f, %8.2f              %6.4f, %6.4f   %6.3f\n",
            r, na, 3na, cert.optimal, point.na, point.k, point.gna, point.gk, pH(eq, model))
end
```

```@example cashplus
fig = plot(; xlabel = "alkali in solution (mmol/kg)", ylabel = "alkali/Si of the gel", legend = :topleft,
           xscale = :log10, size = (720, 400), left_margin = 6Plots.mm, bottom_margin = 7Plots.mm)
for (r, style) in ((1.0, :solid), (1.6, :dash))
    pts = uptake[r]
    plot!(fig, [p.na for p in pts], [p.gna for p in pts]; label = "Na, gel Ca/Si $r", color = :steelblue, linestyle = style, marker = :circle)
    plot!(fig, [p.k for p in pts], [p.gk for p in pts]; label = "K, gel Ca/Si $r", color = :firebrick, linestyle = style, marker = :square)
end
savefig(fig, "cashplus-alkali.svg"); nothing # hide
```

![](cashplus-alkali.svg)

Every paste certifies. The gel takes up both alkalis, a little more potassium
than their ratio in solution: at a Ca/Si of 1, with 96.8 mmol/kg of sodium and
288.8 of potassium in solution, it holds 0.065 Na and 0.228 K per Si. At a Ca/Si
of 1.6 and nearly the same solution it holds four times less sodium and eight
times less potassium, 0.016 and 0.027 per Si: the calcium that fills the
interlayer leaves the alkalis little room, the suppression of alkali uptake at
high Ca/Si that [Miron2022a](@citet) describe.

## 3. Against the authors' own calculation

[Miron2022a](@citet) also publish the model *discretized*: 110 gels of fixed composition,
each with the Gibbs energy their implementation computes for it, given as the
equilibrium constant of its dissolution. The formula of each fixes the site
fractions, so the Gibbs energy of our model at that composition can be set
against theirs, with no solver in between. The first gel is

```@example cashplus
dsp = literature_table("Miron2022a", "dsp_pseudocompounds")
println(dsp.name[1], ":  ", dsp.reaction[1])
```

The block below reads the site fractions from each formula, as the moieties of
the sites give them (silicon, alkali, hydrogen and calcium in turn), evaluates
the Gibbs energy of the model there, and subtracts the one the reaction gives:

```@example cashplus
nk = ChemistryLab.model(phases["CASH+NK"])
lat = nk.lattice
T = 298.15
RT = ChemistryLab.R_GAS * T
G25(sp) = ustrip(us"J/mol", byname[sp][:ΔₐG⁰](T = T * u"K", P = 1.0e5u"Pa"; unit = true))
# The energies of the end-members as Miron et al. (2022a) fit them: the database
# carries TCNh and TCKh as their second paper fine-tunes them.
alkali = literature_table("Miron2022a", "alkali_standard_properties")
Gm(m) = m in alkali.end_member ? ustrip(us"J/mol", alkali.G[findfirst(==(m), alkali.end_member)]) : G25(m)
g = [Gm(symbol(m)) / RT for m in phases["CASH+NK"].end_members]
# The Gibbs energy of the gel per formula unit at site fractions y, written out
# as in the theory: the reference surface, the configurational term of each site
# and the site interactions. A fraction the rounding of a formula makes -1e-4
# counts as an empty site.
xlogx(v) = v > 0 ? v * log(v) : 0.0
function gibbs(y)
    o = lat.occupancy
    reference = RT * sum(g[k] * y[1][o[1, k]] * y[2][o[2, k]] for k in axes(o, 2))
    configurational = RT * sum(sum(xlogx, ys) for ys in y)
    interactions = sum(W * y[s][i] * y[s][l] for (s, i, l, W) in nk.interactions)
    return reference + configurational + interactions
end
misfit = map(eachindex(dsp.name)) do r
    lhs, rhs = split(dsp.reaction[r], " = ")
    c = Dict(Symbol(m[1]) => parse(Float64, m[2]) for m in eachmatch(r"([A-Z][a-z]?)([0-9.]+)", first(split(lhs, " + "))))
    el(e) = get(c, e, 0.0)
    yN, yK = el(:Na), el(:K)
    yCIC = (el(:H) - 6 - yN - yK) / 2        # the interlayer calcium carries three H, an alkali two
    yS = el(:Si) - 2                         # the dimeric unit holds two Si
    yCBT = el(:Ca) - 2 - yCIC                # and two Ca
    y = [[yS, 1 - yS - yCBT, yCBT], [1 - yCIC - yN - yK, yCIC, yN, yK]]   # BT (S, v, C), IC (v, C, N, K)
    products = sum(split(rhs, " + ")) do term
        m = match(r"^\s*([0-9.]*)(\S+)$", term)
        (isempty(m[1]) ? 1.0 : parse(Float64, m[1])) * G25(String(m[2]))
    end
    gibbs(y) - (products - ustrip(us"J/mol", dsp.dG[r]))
end
@printf("%d gels: largest difference %.3f kJ/mol, median %.3f kJ/mol\n",
        length(misfit), maximum(abs, misfit) / 1000, sort(misfit)[end ÷ 2] / 1000)
```

The two implementations agree on every gel to within 0.064 kJ/mol, with a median
difference of 0.003 kJ/mol. That is the rounding of the formulas, which are
printed to four decimals; the energies of mixing and of the reciprocal reactions
are some ten kilojoules per mole.

## 4. A Portland cement with limestone (PC4)

[Miron2022b](@citet) tuned the calcium alkali end-members of CASH+NK
on the pore solutions of hydrated cements, among them the Portland cement with
4 % limestone (PC4) of [LothenbachLeSaout2008](@citet), whose
pore solution was analyzed from one day to 400 days. The same paste is computed
here twice, with its C-S-H as `CSHQ` and as `CASH+NK`, everything else equal:

- **the cement** is the normative composition of
  [LothenbachLeSaout2008; Table 1](@citet): the four clinker phases, periclase,
  free lime, calcite, gypsum and the readily soluble alkali sulfates;
- **the clinker phases hydrate** by the law of [ParrottKilloh1984](@citet) with the
  constants of [LothenbachLeSaout2008; Table 3](@citet), including the two it
  adapts for belite and the critical degree of hydration of each phase;
- **the minor oxides of the clinker** (0.052 g of K₂O, 0.31 g of Na₂O, 0.87 g of
  MgO and 0.11 g of SO₃ per 100 g) are released with the phases that hold them.
  [LothenbachLeSaout2008](@citet) give their totals, not how they are shared
  among the phases. This page assumes the sharing of
  [LothenbachWinnefeld2006](@citet), after Taylor, as a content per gram of each
  phase;
- **the phases that may form** are those of the Portland paste of
  [the validation page](@ref ex-validation), a paste of the same
  laboratory and the same modeling.

The recipe and both systems are written once, in `scripts/lothenbach_2008.jl`,
which the test of this page includes too.

```@example cashplus
include(joinpath(pkgdir(ChemistryLab), "scripts", "lothenbach_2008.jl"))
days = l08_days()
runs = Dict(gel => hydrate(l08_recipe("PC4"), l08_system(gel), days; model) for gel in (:CSHQ, :CASHNK))
for gel in (:CSHQ, :CASHNK)
    @printf("%-7s certified at %s of %d ages\n", gel, count(rs -> rs.certificate.optimal, runs[gel].states), length(days))
end
```

The table sets the two calculations against the analysis, at each age. The
model gives millimoles per kilogram of water and the analysis millimoles per liter
of solution; they differ by a few percent at these concentrations.

```@example cashplus
elements = ["Na", "K", "Ca", "Si", "S"]
println("  age   element   measured    CSHQ   CASH+NK")
for (k, d) in enumerate(days)
    with = Dict(gel => pore_solution_mmol(runs[gel].states[k].state) for gel in (:CSHQ, :CASHNK))
    for e in elements
        @printf("%5.0f d  %-7s  %9.3g  %8.3g  %8.3g\n", d, e, l08_measured("PC4", d, e), with[:CSHQ][e], with[:CASHNK][e])
    end
    @printf("%5.0f d  %-7s  %9.3g  %8.3g  %8.3g\n", d, "pH", l08_measured("PC4", d, "pH"),
            (pH(runs[gel].states[k].state, model) for gel in (:CSHQ, :CASHNK))...)
end
```

```@example cashplus
panels = map(("Na", "K")) do e
    p = plot(; xscale = :log10, title = e, xlabel = "time (days)", ylabel = "mmol/kg or mmol/L",
             legend = e == "Na" ? :topleft : false)
    for (gel, color) in ((:CSHQ, :darkorange), (:CASHNK, :steelblue))
        plot!(p, days, [pore_solution_mmol(rs.state)[e] for rs in runs[gel].states];
              label = gel === :CSHQ ? "CSHQ" : "CASH+NK", color, linewidth = 2, marker = :circle)
    end
    scatter!(p, days, [l08_measured("PC4", d, e) for d in days]; label = "measured", color = :black)
    p
end
fig = plot(panels...; layout = (1, 2), size = (900, 380), left_margin = 6Plots.mm, bottom_margin = 7Plots.mm)
savefig(fig, "cashplus-pc4.svg"); nothing # hide
```

![](cashplus-pc4.svg)

Both models certify the paste at every age. From 7 days on, CASH+NK leaves more
sodium and potassium in solution than CSHQ, and closer to the analysis: at 28
days 160 and 338 mmol/kg against 172 and 532 measured, where CSHQ gives 105 and
290. That is the improvement [Miron2022b](@citet) report. Neither model follows the rise
of both alkalis after 28 days, to 331 and 563 mmol/L at 400 days: CASH+NK stays
near 164 and 334, CSHQ near 108 and 281. The pH follows the same order, 13.7 with
CASH+NK, 13.6 with CSHQ and 13.7 to 13.8 measured. After the first day the
calcium is 1.1 mmol/kg with CASH+NK and 1.0 with CSHQ, against 1.6 to 1.0
measured; with CASH+NK it is carried in part by the Ca(OH)₂⁰ complex of the
alkali fit. The sulfate is the largest
difference of the two models alike: measured, it falls to 1.9 mmol/L at one day
and climbs to 34 to 41 after six months; computed, it stays in solution at one
day (53 and 85 mmol/kg) and between 2.7 and 6.2 afterwards.

## 5. The other cations: lithium to radium

[Miron2022a](@citet) extended the interlayer site further, to lithium, rubidium and
cesium and to magnesium, strontium, barium and radium, and fitted each on uptake
experiments [Miron2022a](@cite). `CASH+ext` is that model: 33 end-members, with
the interaction parameters of their Table A3. Their Table A2 gives the equilibrium
constant of the dissolution of each end-member into Ca²⁺, SiO₂⁰, water and its
interlayer cation, computed from their own energies; the database gives them back:

```@example cashplus
ext = phases["CASH+ext"]
printed = literature_table("Miron2022a", "dissolution_log_K")
cation = Dict("N" => "Na+", "K" => "K+", "Li" => "Li+", "Rb" => "Rb+", "Cs" => "Cs+",
              "Mg" => "Mg+2", "Sr" => "Sr+2", "Ba" => "Ba+2", "Ra" => "Ra+2")
# The energies `Gm` of section 3: the first paper's for the alkali members, whose
# TCNh and TCKh the database carries as the second paper fine-tunes them.
worst = maximum(ext.end_members) do m
    n = symbol(m)
    c = atoms(m)
    x = get(cation, n[3:(end - 1)], nothing)
    hplus = 2c[:Ca] + (x === nothing ? 0 : charge(byname[x]))
    products = c[:Ca] * G25("Ca+2") + c[:Si] * G25("SiO2@") + (c[:H] + hplus) / 2 * G25("H2O@") + (x === nothing ? 0.0 : G25(x))
    abs(-(products - Gm(n)) / (RT * log(10)) - printed.log_K[findfirst(==(n), printed.species)])
end
@printf("%d end-members: largest difference from Table A2, %.3f in log K\n", length(ext.end_members), worst)
```

With strontium and cesium at a millimole per kilogram beside the sodium and
potassium of section 2, at a Ca/Si of 1.2, the gel takes up both:

```@example cashplus
sp = speciation(substances, vcat(["Portlandite", "Amor-Sl"], symbol.(ext.end_members)); aggregate_state = [AS_AQUEOUS])
complexes(s) = charge(s) == 0 && any(el -> haskey(atoms(s), el), (:Na, :K, :Li, :Rb, :Cs)) && aggregate_state(s) == AS_AQUEOUS
# Cemdata18 has no component for the five cations the extension adds.
primaries = vcat(CEMDATA_PRIMARIES, ["Li+", "Rb+", "Cs+", "Ba+2", "Ra+2"])
cs_ext = ChemicalSystem(filter(!complexes, sp), primaries; solid_solutions = [ext])
st = ChemicalState(cs_ext)
set_quantity!(st, "H2O@", 1.0u"kg")
set_quantity!(st, "Amor-Sl", 0.05u"mol")
set_quantity!(st, "Portlandite", 0.06u"mol")
for (s, x) in ("Na+" => 0.02, "K+" => 0.06, "Sr+2" => 0.001, "Cs+" => 0.001)
    set_quantity!(st, s, x * u"mol")
end
set_quantity!(st, "OH-", 0.083u"mol")
eq, cert = equilibrate_certified(st; model)
gel = solid_solution_totals(eq, "CASH+ext").elements
n = ustrip.(us"mol", eq.n)
held(el) = gel[el] / sum(get(atoms(s), el, 0) * x for (s, x) in zip(cs_ext.species, n))
@printf("certified: %s; the gel holds %.0f %% of the Sr, %.0f %% of the Cs, %.0f %% of the Na and %.0f %% of the K\n",
        cert.optimal, 100held(:Sr), 100held(:Cs), 100held(:Na), 100held(:K))
```

The database gives back Table A2 for the 33 end-members, to 0.006 in log K. In
that solution the gel certifies, and holds 78 % of the strontium but 4 % of the
cesium, 2 % of the sodium and 2 % of the potassium. Strontium takes the place of
the interlayer calcium one for one (TSSrh is TSCh with Sr for one Ca), where a
monovalent cation takes it with one hydroxide less (TSNh). The
uptake experiments on which [Miron2022a](@citet) fitted each cation are not reproduced
here.

## Where to go next

The derivation of the model is in [Solid solutions](@ref sec-theory-solid-solutions),
section 8, and [The models of the C-S-H gel](@ref sec-csh-models) lists the models
of the gel the package ships, CASH+ among them.
