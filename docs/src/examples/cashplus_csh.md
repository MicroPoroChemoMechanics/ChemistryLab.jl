# [C-S-H in water and in alkaline solutions: the CASH+ model](@id ex-cashplus-csh)

!!! info "Before this page"
    [Solid solutions](@ref sec-theory-solid-solutions), sections 7 and 8, where
    the site model and the compound energy formalism are derived.

CASH+ is a model of the C-S-H gel by Kulik, Miron and Lothenbach
[Kulik2022](@cite), which Miron et al. extended to sodium and potassium
[Miron2022a, Miron2022b](@cite). Like the CNASH gel of
[the site-mixing page](@ref ex-sublattice-csh), it mixes on the structural sites
of the silicate chain and of the interlayer. It also adds two terms: the energy
of the reciprocal reactions between its end-members, and interactions between
the species of one site. ChemistryLab writes it as a [`CompoundEnergyModel`](@ref).

This page does three things. It computes the gel in water, from the Ca/Si at
which amorphous silica stops forming to the one at which portlandite starts, and
compares both ends with the paper. It then adds sodium and potassium. Finally it
checks the model against the authors' own calculation of 110 gel compositions.

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

# Cemdata18 with the twelve end-members of CASH+NK, and the CaSiO3@ complex the
# model was fitted with (built on first use).
substances = build_species(datapath("cemdata18-cashplus.json"); verbose = false)
byname = Dict(symbol(s) => s for s in substances)
phases = Dict(ChemistryLab.name(p) => p for p in build_solid_solutions(datapath("solid_solutions.toml"), byname))

# The activity model the paper fits with: extended Debye-Hückel with the ion size
# and B-dot Cemdata18 gives for KOH solutions.
model = cemdata18_activity_model(:KOH)

# The gel and the aqueous species of its elements, with portlandite and
# amorphous silica, which may form beside it. The neutral complexes of the
# alkalis (NaOH@, KOH@ ...) are left out, as Miron et al. left them out when they
# fitted the alkali end-members.
function gel_system(gel)
    members = symbol.(phases[gel].end_members)
    sp = speciation(substances, vcat(["Portlandite", "Amor-Sl"], members); aggregate_state = [AS_AQUEOUS])
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
    b = Float64.(cs.SM.A) * ustrip.(us"mol", st.n)
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
more of the solid that forms. The paper gives both ends:

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

Both ends agree with the paper to the digits it prints. The figure shows the
solution and the sites between them:

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

## 2. Sodium and potassium: CASH+NK

`CASH+NK` puts sodium and potassium in the interlayer, as NaHOH⁺ and KHOH⁺. With
twelve end-members, it is the one to use for a cement paste. The pastes below
have a Ca/Si of 1 and of 1.6, in a solution of sodium and potassium hydroxides
with three times as much potassium as sodium, as in the pore solution of a
Portland cement. Potassium dominates, so the Debye-Hückel parameters are those
Cemdata18 gives for KOH:

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
    b = Float64.(cs_nk.SM.A) * ustrip.(us"mol", st.n)
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
288.8 of potassium in solution, it holds 0.065 Na and 0.227 K per Si. At a Ca/Si
of 1.6 and nearly the same solution it holds four times less sodium and nine
times less potassium, 0.015 and 0.026 per Si: the calcium that fills the
interlayer leaves the alkalis little room, the suppression of alkali uptake at
high Ca/Si that Miron et al. describe.

!!! note "Both alkalis, or the core model"
    Declared in a system whose budget holds no potassium, or no sodium, the
    twelve-member phase keeps members of an element that cannot be there. Measured
    on these pastes with sodium alone, the search then stops at an element balance
    of about 2e-8 and does not certify at the lowest sodium. A system without
    alkalis takes `CASH+`, as in section 1; a cement paste holds both.

## 3. Against the authors' own calculation

Miron et al. also publish the model *discretized*: 110 gels of fixed composition,
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

## Where to go next

The derivation of the model is in [Solid solutions](@ref sec-theory-solid-solutions),
section 8, and [The models of the C-S-H gel](@ref sec-csh-models) lists the models
of the gel the package ships, CASH+ among them.
