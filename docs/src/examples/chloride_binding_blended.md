# [Chloride binding in blended cements](@id sec-example-chloride-blended)

!!! info "Before this page"
    [Chloride binding by C-S-H and Friedel's salt](@ref sec-example-csh-chloride)
    and [A blastfurnace cement](@ref ex-cem3-slag), whose paste this page salts.

A blended cement binds the chloride that enters it in its AFm phases, as
Kuzel's and Friedel's salts, and in its C-S-H. The salts are phases of
CEMDATA18. The C-S-H is its CSHQ solid solution, whose end members hold calcium
and alkalis but no chloride. This page salts the CEM III/A paste of the
[blastfurnace cement page](@ref ex-cem3-slag) along two routes that give the
C-S-H its share, then a CEM III/B paste without portlandite, along the one route
that applies there.

  - **Route A** keeps the surface model of the C-S-H that [Guo2018](@cite) use:
    silanol sites that bind calcium, alkalis and chloride. It puts them on a
    C-S-H whose amount and composition a first equilibrium with CSHQ has
    decided, and which [`freeze_solid_solution`](@ref) then sets aside.
  - **Route B** adds an end member to CSHQ, `CSHQ-Cl` = (CaCl₂)₀.₅, whose Gibbs
    energy was fitted on the sorption tests of [Hirao2005](@cite). It ships in
    `cemdata18-chloride.json`, and the phase is `CSHQ_Cl`.

!!! warning "What these calculations are, and what they are not"
      - Neither route says where the chloride sits in the gel.
        [Plusquellec2016](@cite) found that chloride does not adsorb
        specifically on C-S-H. What a depletion measurement counts as bound accompanies the
        calcium the surface adsorbs, in the diffuse layer that screens it. The
        surface complex of route A and the end member of route B are both
        effective descriptions of that uptake.
      - Route A freezes the gel at the composition of the chloride-free paste.
        That is a sound approximation only while portlandite fixes the
        activities of calcium and hydroxide, and
        [`freeze_solid_solution`](@ref) refuses a paste without it.
      - The end member of route B rests on three points of one gel at
        portlandite saturation. The fit is within 0.05 mmol/g at 0.5 and
        1 mol/L, and four times too high at 0.1 mol/L. Its dependence on the
        Ca/Si is the model's: the end member was chosen to follow the trend
        measured by [Zibara2008](@cite), not fitted to it
        (`data/chloride/README.md`).
      - The B-dot activity model of the CEMDATA18 pages holds up to an ionic
        strength of about 1 mol/kg. No solution below exceeds 0.7.

```@example clblend
using ChemistryLab, DynamicQuantities, JSON, OptimaSolver, OrderedCollections, Printf

substances = build_species(datapath("cemdata18-chloride.json"); verbose = false)
byname = Dict(symbol(s) => s for s in substances)

# The recipe of the CEM III page, with its assumptions: half clinker, half slag;
# the clinker of Lavergne2018, Table 9; a European slag analysis; the alkalis at
# a usual industrial level; w/b 0.40; the clinker reacted to the water ceiling
# and the slag to the mean of the round robin of Durdzinski2017 at 28 days.
bogue = literature_table("Lavergne2018", "cement_bogue")
CLINKER = OrderedDict(zip(bogue.phase, bogue.percent ./ 100))
SLAG = Dict("CaO" => 0.41, "SiO2" => 0.36, "Al2O3" => 0.11, "MgO" => 0.08, "SO3" => 0.02)
ALKALIS = Dict("K2O" => 0.008, "Na2O" => 0.002)
WB, BINDER_G = 0.40, 100.0
sem(material, age) = literature_table("Durdzinski2017", "degree_of_reaction";
    technique = "SEM-IA", material, curing = "sealed", age_days = age).degree_percent
ALPHA_CLINKER = powers_alpha_max(WB)
ALPHA_SLAG = min(sum([sem("S1", 28); sem("S2", 28)]) / 400, ALPHA_CLINKER)

gems = JSON.parsefile(joinpath(pkgdir(ChemistryLab), "test", "reference", "gems_cemdata18_portland.json"))
lg1, lg2 = log10(gems["gamma"]["z1"]), log10(gems["gamma"]["z2"])
model = HKFActivityModel(å = 0.0, Ḃ = (lg1 + (lg1 - lg2) / 3) / gems["ionic_strength_mol_per_kg"], Kₙ = 0.0)

pure = split(
    "C3S C2S C3A C4AF Gp Anh Cal Portlandite ettringite monosulphate12 " *
    "monocarbonate hemicarbonate C4AH13 C3AH6 C3FH6 straetlingite " *
    "hydrotalcite Mg2AlC0.5OH Brc FeOOHmic AlOHmic Amor-Sl Tro Py Sulfur Mgs"
)
redox = ["HS-", "H2S@", "SO4-2", "SO3-2", "S2O3-2", "O2@", "H2@"]
salts = ["C4AClH10", "C4AsClH12"]
phases = build_solid_solutions(datapath("solid_solutions.toml"), byname)
phase(n) = only(p for p in phases if name(p) == n)
hydrogarnet = phase("C3(AF)S0.84H")
members(p) = symbol.(end_members(p))
molar(s) = ustrip(us"g/mol", byname[s][:M])
M_Cl = ustrip(us"g/mol", calculate_molar_mass(Dict(:Cl => 1)))

"""
The reacted paste per 100 g of binder, with `dose` percent of chloride by mass of
binder added as NaCl: the state that carries the clinker and the water, and the
budget that adds the alkalis and the slag, on the components of `cs`.
"""
function paste(cs; dose = 0.0, clinker_fraction = 0.5)
    st = ChemicalState(cs)
    for (ph, frac) in CLINKER
        set_quantity!(st, ph, ALPHA_CLINKER * BINDER_G * clinker_fraction * frac / molar(ph) * u"mol")
    end
    set_quantity!(st, "H2O@", BINDER_G * WB / molar("H2O@") * u"mol")
    n_cl = dose / 100 * BINDER_G / M_Cl
    if n_cl > 0
        set_quantity!(st, "Na+", n_cl * u"mol")
        set_quantity!(st, "Cl-", n_cl * u"mol")
    end
    b = Float64.(cs.SM.A) * ustrip.(us"mol", st.n)
    b .+= oxide_budget(ALKALIS, cs.SM.primaries; mass = BINDER_G * clinker_fraction * ALPHA_CLINKER * u"g")
    b .+= oxide_budget(SLAG, cs.SM.primaries; mass = BINDER_G * (1 - clinker_fraction) * ALPHA_SLAG * u"g")
    return (; state = st, b)
end
nothing # hide
```

## The paste, before any chloride

The first equilibrium is that of the CEM III page, recomputed. The paste holds
portlandite, and its C-S-H holds part of the alkalis, which route A returns to
the solution when it freezes the gel.

```@example clblend
cshq = phase("CSHQ")
sp1 = speciation(substances, vcat(pure, members(cshq), members(hydrogarnet), redox);
                 aggregate_state = [AS_AQUEOUS])
cs1 = ChemicalSystem(sp1, CEMDATA_PRIMARIES; solid_solutions = [cshq, hydrogarnet])
p1 = paste(cs1)
eq1, cert1 = equilibrate_certified(p1.state; model, b = p1.b)
gel = solid_solution_totals(eq1, "CSHQ")
mol(eq, s) = ustrip(us"mol", moles(eq, byname[s]))
water(eq) = mol(eq, "H2O@") * molar("H2O@") / 1000
@printf("certified %s, pH %.3f, portlandite %.4f mol, pore water %.1f g\n",
        cert1.optimal, pH(eq1, model), mol(eq1, "Portlandite"), 1000water(eq1))
@printf("C-S-H: %.4f mol Si, Ca/Si %.3f, holding %.4f mol Na and %.4f mol K\n",
        gel.elements[:Si], gel.elements[:Ca] / gel.elements[:Si], gel.elements[:Na], gel.elements[:K])
```

## Route A: the published surface on a frozen gel

The second system carries the hydrates, the chloride salts and the aqueous
species of the first, and the sites. The C-S-H is no longer one of its phases.
The sites and their area are referred to the silicon of the frozen gel, through
the molar mass of Guo's C-S-H per silicon: 0.004 mol/g and 500 m²/g, times
191.4 g/mol. The gel of the first stage then carries as many sites per silicon
as Guo's does. The alkalis of NaSiOH and KSiOH go back to the solution as
hydroxides.

```@example clblend
sp2 = speciation(substances, vcat(pure, salts, members(hydrogarnet), redox, ["Na+", "K+", "Cl-"]);
                 aggregate_state = [AS_AQUEOUS])
guo(q) = literature_value("Guo2018", q)
t = literature_table("Guo2018", "csh_formula")
per_si = Dict(zip(t.oxide, t.coefficient ./ t.coefficient[findfirst(==("SiO2"), t.oxide)]))
M_guo = sum(k * ustrip(us"g/mol", Species(ox)[:M]) for (ox, k) in per_si)
n_sites = ustrip(u"mol/g", guo("csh_site_density")) * M_guo * gel.elements[:Si]
area = ustrip(u"m^2/g", guo("csh_specific_surface_area")) * M_guo * gel.elements[:Si]
rows = literature_table("Guo2018", "surface_reactions_phreeqc")
aqueous = [s for s in sp2 if aggregate_state(s) == AS_AQUEOUS]
family = site_family(
    "Csh_w", [eq => lk for (eq, lk) in zip(rows.reaction, rows.log_K)], aqueous;
    master = "Csh_w", site = "Xw", capacity = TotalSiteAmount(n_sites * u"mol"),
    support = SurfaceSupport("C-S-H", nothing, FixedSurfaceArea(area)), model = DiffuseLayer(; area),
)
sites = vcat([family.free_site], family.complexes)
cs2 = ChemicalSystem(vcat(sp2, sites),
                     vcat([s for s in sp2 if symbol(s) in CEMDATA_PRIMARIES], [family.free_site]);
                     solid_solutions = [hydrogarnet], site_families = [family])
frozen = freeze_solid_solution(eq1, "CSHQ", cs2; buffer = "Portlandite")
@printf("%.3f mol of sites on %.3g m², for %.1f g of pore water: a film %.2f nm thick\n",
        n_sites, area, 1000water(eq1), 1.0e9 * water(eq1) / 1000 / area)
```

At that area the pore solution is a film 0.57 nm thick on the gel, thinner than
the diffuse layer that screens the surface. The ions of the layer cannot be
counted by a [`DonnanLayer`](@ref), which assumes a free solution beyond it, and
the Gouy-Chapman potential is itself an extrapolation here. The layer's ions
are left implicit, as PHREEQC's default leaves them.

The paste is then salted with up to 0.4 % of chloride by mass of binder, added
as NaCl, and solved from the chloride-free point up, each solve starting from
the last. The table gives the free chloride (mol/kg of pore water), the ionic
strength, the pH, the portlandite, and the chloride in each place (mol per
100 g of binder): two per Friedel's salt, one per Kuzel's salt, one per surface
complex.

```@example clblend
DOSES = [0.0, 0.05, 0.1, 0.2, 0.4]          # percent of chloride by mass of binder
idx2 = Dict(symbol(s) => i for (i, s) in enumerate(cs2.species))
function salted(dose)
    n = ustrip.(us"mol", frozen.state.n)
    n[idx2["XwOH"]] = n_sites
    n_cl = dose / 100 * BINDER_G / M_Cl
    n[idx2["Na+"]] += n_cl
    n[idx2["Cl-"]] += n_cl
    return ChemicalState(cs2, n .* u"mol")
end
A2 = Float64.(cs2.SM.A)
route_a, certs_a = equilibrate_path(salted(0.0), [A2 * ustrip.(us"mol", salted(d).n) for d in DOSES]; model)
println("certified: ", count(c -> c.optimal, certs_a), " of ", length(certs_a))

chlorine(cs) = [Float64(get(atoms(s), :Cl, 0)) for s in cs.species]
function partition(eq, cs, dose; gel_species)
    n = ustrip.(us"mol", eq.n)
    cl = chlorine(cs)
    free = sum(cl[i] * n[i] for i in cs.idx_aqueous)
    added = dose / 100 * BINDER_G / M_Cl
    at(s) = n[findfirst(x -> symbol(x) == s, cs.species)]
    return (; free_molal = free / water(eq), friedel = 2at("C4AClH10"), kuzel = at("C4AsClH12"),
            gel = at(gel_species), bound = added > 0 ? 1 - free / added : 0.0,
            CH = at("Portlandite"), pH = pH(eq, model), I = ionic_strength(eq))
end
function report(states, cs, gel_species)
    println(" Cl %   free Cl   I      pH     CH      Friedel  Kuzel    C-S-H    bound")
    for (d, eq) in zip(DOSES, states)
        r = partition(eq, cs, d; gel_species)
        @printf("%5.2f  %7.4f  %5.3f  %6.3f  %.4f  %.5f  %.5f  %.5f  %5.1f %%\n",
                d, r.free_molal, r.I, r.pH, r.CH, r.friedel, r.kuzel, r.gel, 100r.bound)
    end
end
report(route_a, cs2, "XwOHCl-")
```

Two numbers belong to the model rather than to the paste. Before any chloride
is added, the surface binds 0.087 mol of calcium, taken from the portlandite:
81 % of it. The pH rises from 13.04 to 13.18 for that reason alone. The C-S-H of
the first stage already held that calcium in its Ca/Si. The published model
counts calcium on the surface on top of a bulk composition that, measured,
includes it, and freezing the gel keeps that assumption. The portlandite
nonetheless remains at every dose, so the condition the frozen gel rests on
holds.

The chloride goes to Kuzel's salt, which holds half a sulfate and one chloride
per formula unit. At 0.4 % it holds 0.0104 mol, and the 0.0056 mol of
monosulfate of the first stage carry sulfate for 0.0112. The surface holds
0.00018 mol from the lowest dose on: monosulfate and Kuzel's salt coexist and
fix the chloride activity, as in [the one-stage paste](@ref sec-example-csh-chloride).
No Friedel's salt forms at these doses.

## Route B: chloride in the solid solution

Route B needs no second stage: the paste and the salt are equilibrated
together, with `CSHQ_Cl` in place of CSHQ. The phase is read from
`data/solid_solutions.toml` like the others; with CEMDATA18 alone its entry is
skipped, since its end member is not there.

```@example clblend
cshq_cl = phase("CSHQ_Cl")
spB = speciation(substances, vcat(pure, salts, members(cshq_cl), members(hydrogarnet), redox, ["Cl-"]);
                 aggregate_state = [AS_AQUEOUS])
csB = ChemicalSystem(spB, CEMDATA_PRIMARIES; solid_solutions = [cshq_cl, hydrogarnet])
AB = Float64.(csB.SM.A)
route_b, certs_b = equilibrate_path(paste(csB).state, [paste(csB; dose = d).b for d in DOSES]; model)
println("certified: ", count(c -> c.optimal, certs_b), " of ", length(certs_b))
report(route_b, csB, "CSHQ-Cl")
```

The two routes agree on the salts where most of the chloride is: Kuzel's salt
differs by 2 % at 0.4 %, and by a third at 0.05 %, where little of it forms. They
differ on the C-S-H. The end member holds 0.00063 to 0.00066 mol of chloride
where the surface holds 0.00018, and the paste binds 83 to 96 % of the chloride
against 65 to 94 %. At the lowest dose the C-S-H holds half of the bound
chloride in route B and a fifth in route A; at 0.4 % the AFm phases hold 94 to
98 % of it in both. Route B leaves the portlandite where the first stage put it.

## Without portlandite

A CEM III/B at the lower bound of its clinker content, 20 %, with the same
clinker, slag and reacted fractions, has no portlandite left, and its C-S-H has
a Ca/Si of 1.52. Route A has nothing to rest on, and
[`freeze_solid_solution`](@ref) says so. Route B applies unchanged.

```@example clblend
CEM_III_B = 0.20        # clinker, the lower bound of CEM III/B in EN 197-1
p3 = paste(cs1; clinker_fraction = CEM_III_B)
eq3, cert3 = equilibrate_certified(p3.state; model, b = p3.b)
gel3 = solid_solution_totals(eq3, "CSHQ")
@printf("certified %s, pH %.3f, portlandite %.4f mol, Ca/Si of the gel %.3f\n", cert3.optimal,
        pH(eq3, model), mol(eq3, "Portlandite"), gel3.elements[:Ca] / gel3.elements[:Si])
try
    freeze_solid_solution(eq3, "CSHQ", cs2; buffer = "Portlandite")
catch e
    println(first(split(sprint(showerror, e), ", so ")), ".")
end
route_c, certs_c = equilibrate_path(paste(csB; clinker_fraction = CEM_III_B).state,
    [paste(csB; dose = d, clinker_fraction = CEM_III_B).b for d in DOSES]; model)
println("certified: ", count(c -> c.optimal, certs_c), " of ", length(certs_c))
report(route_c, csB, "CSHQ-Cl")
```

Without portlandite the salt moves the pH further, from 12.67 to 12.97. The gel
holds 0.0009 to 0.0014 mol of chloride, more than in the CEM III/A paste:
four fifths of the bound chloride at the lowest dose, a seventh at the highest.
No measurement on a slag cement checks these numbers; they are what the fitted
end member predicts for a gel it was not fitted on.

## See also

  - [Chloride binding by C-S-H and Friedel's salt](@ref sec-example-csh-chloride),
    the surface model on Guo's paste, against PHREEQC.
  - `data/chloride/README.md`, the data, the fit and the rejected candidate.
  - `test/chloride_blended_reference.jl` and `test/cshq_chloride.jl`.
