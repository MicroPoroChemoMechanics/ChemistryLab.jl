# [Mixing on sites: the CSH3T and CNASH gels](@id ex-sublattice-csh)

!!! info "Before this page"
    [Solid solutions](@ref sec-theory-solid-solutions), section 7, where the
    site model is derived.

A C-S-H gel is modeled as a solid solution of a few end-members, and there are
two ways of mixing them. The simpler one mixes the **end-members** themselves,
as if each were a molecule placed at random. The other mixes what actually
changes from one end-member to the next: the few **structural sites** of the
silicate chain and of the interlayer, each held by one species or another. The
second is how [Kulik2011](@citet) defines the CSH3T gel and
[Myers2014](@citet) the C-(N-)A-S-H gel, and it is what
[`SublatticeModel`](@ref) implements.

This page does two things. It solves the same CSH3T gel in water both ways and
shows how much the pore solution depends on the choice. Then it computes a
C-(N-)A-S-H gel with the model of [Myers2014](@citet) and reads from it the minimum chain
length of their Eq. (11).

```@example sublattice
using ChemistryLab
using DynamicQuantities
using OptimaSolver
using Printf
using Plots
default(framestyle = :box, grid = false)

substances = build_species(datapath("cemdata18-thermofun.json"); verbose = false)
byname = Dict(symbol(s) => s for s in substances)

# No alkali in the first part and little in the second: the ion size and B-dot of
# Cemdata18 for KOH solutions are used, as on the cement pages.
model = cemdata18_activity_model(:KOH)

# A gel phase and the aqueous species of its elements, plus the pure phases that
# may form beside it.
function gel_system(names, mixing; pure)
    sp = speciation(substances, vcat(pure, names); aggregate_state = [AS_AQUEOUS])
    gel = SolidSolutionPhase("C-S-H", [byname[m] for m in names]; model = mixing)
    return ChemicalSystem(sp, CEMDATA_PRIMARIES; solid_solutions = [gel])
end

# Moles of an element in solution, per kilogram of water.
function in_solution(eq, element)
    cs = eq.system
    n = ustrip.(us"mol", eq.n)
    iw = only(cs.idx_solvent)
    kg_water = n[iw] * ustrip(us"kg/mol", cs.species[iw][:M])
    return sum(get(atoms(cs.species[i]), element, 0) * n[i] for i in cs.idx_solutes) / kg_water
end
nothing # hide
```

## 1. CSH3T, mixed by end-members and by sites

CSH3T has three end-members: TobH, the silica-rich one, T2C, the calcium-rich
one, and T5C between them, with a Ca/Si of 1. Cemdata18 [Lothenbach2019](@cite)
ships them as an ideal mixture of end-members. Their site form puts Si or Ca on
two bridging positions of the silicate chain: TobH has Si on both, T2C Ca on
both, T5C Ca on one and Si on the other. The same three records serve both
declarations:

```@example sublattice
CSH3T = ["CSH3T-TobH", "CSH3T-T5C", "CSH3T-T2C"]
members = [byname[m] for m in CSH3T]
pure = ["Portlandite", "Amor-Sl"]
systems = ["end-members" => gel_system(CSH3T, IdealSolidSolutionModel(); pure),
           "sites" => gel_system(CSH3T, sublattice_model("Kulik2011:csh3t", members); pure)]
nothing # hide
```

Each paste below is 0.05 mol of silica and some lime in a kilogram of water, from
a Ca/Si of 0.7 to 1.6. That is how the solubility of a C-S-H is measured: the gel
is equilibrated in water and the solution analyzed.

```@example sublattice
ratios = 0.7:0.1:1.6
curves = Dict(label => (gel = Float64[], si = Float64[], ca = Float64[]) for (label, _) in systems)
println("Ca/Si   mixing        certified   Ca (mmol/kg)   Si (mmol/kg)    pH     gel Ca/Si")
for r in ratios, (label, cs) in systems
    st = ChemicalState(cs)
    set_quantity!(st, "Amor-Sl", 0.05u"mol")
    set_quantity!(st, "Portlandite", r * 0.05u"mol")
    set_quantity!(st, "H2O@", 1.0u"kg")
    b = budget(st)
    eq, cert = equilibrate_certified(st; model, b)
    gel = solid_solution_totals(eq, "C-S-H").elements
    c = curves[label]
    push!(c.gel, gel[:Ca] / gel[:Si]); push!(c.si, 1000in_solution(eq, :Si)); push!(c.ca, 1000in_solution(eq, :Ca))
    @printf("%4.2f    %-12s  %-9s  %10.3f     %10.4f     %6.3f   %6.3f\n",
            r, label, cert.optimal, c.ca[end], c.si[end], pH(eq, model), c.gel[end])
end
```

```@example sublattice
fig = plot(; xlabel = "Ca/Si of the gel", ylabel = "silicon in solution (mmol/kg)",
           yscale = :log10, title = "CSH3T in water", legend = :topright, size = (720, 420),
           left_margin = 8Plots.mm, bottom_margin = 8Plots.mm)
for (label, color) in (("end-members", :steelblue), ("sites", :firebrick))
    plot!(fig, curves[label].gel, curves[label].si; label = "mixed by $label",
          marker = :circle, color, linewidth = 2)
end
savefig(fig, "sublattice-csh3t.svg"); nothing # hide
```

![](sublattice-csh3t.svg)

Both declarations certify every paste, and they agree on the gel: its Ca/Si
follows the lime added, and the two differ by 0.05 at most. They disagree on the
solution. Mixed by sites, the gel leaves more silicon in solution at every Ca/Si
from 0.8 up, half as much again at a Ca/Si of 1.5 (0.063 against 0.041 mmol/kg),
and from 1.4 up less calcium. The two models give the three end-members different
activities at the same composition, so the solution in equilibrium with a given
gel is another one. Which of the two is closer to measured solubilities is not
settled here; what the table settles is how much the choice weighs.

## 2. The C-(N-)A-S-H gel and its chain length

The CNASH gel of [Myers2014](@citet) has eight end-members on six sites. Its formula
unit counts the silicate chain in *dreierketten* units, three tetrahedra each:
two paired ones, and a bridging one that can be vacant. Their Table 1 gives, for
each end-member, the fraction ``\nu`` of bridging sites that are vacant, and
their Eq. (11) turns it into the length of the chains, in tetrahedra, taken as not
crosslinked:

```math
\mathrm{CL} = \frac{3}{\sum_k \chi_k \nu_k} - 1 ,
```

with ``\chi_k`` the mole fractions of the end-members in the gel. A gel with
every bridging site vacant is made of dimers, ``\mathrm{CL} = 2``; one with none
vacant has infinite chains. The table reads each end-member's composition from
its database record and its ``\nu`` from the transcription of
[Myers2014; Table 1](@cite):

```@example sublattice
CNASH = ["T2C-CNASHss", "T5C-CNASHss", "TobH-CNASHss", "5CA", "5CNA", "INFCA", "INFCN", "INFCNA"]
cnash = sublattice_model("Myers2014:cnash", [byname[m] for m in CNASH])
table1 = literature_table("Myers2014", "cnash_bridging_vacancies")
ν = Dict(zip(table1.end_member, ustrip.(table1.bridging_vacancies)))

println("end-member     Ca/Si   Al/Si   Na/Si    ν     chain length")
for m in CNASH
    a = atoms(byname[m])
    ratio(el) = get(a, el, 0) / a[:Si]
    @printf("%-13s  %5.3f   %5.3f   %5.3f   %3.1f   %s\n", m, ratio(:Ca), ratio(:Al), ratio(:Na),
            ν[m], ν[m] == 0 ? "infinite" : @sprintf("%.0f", 3 / ν[m] - 1))
end
```

Four end-members have no vacant bridging site, and so infinite chains: TobH and
the three whose name starts with INF. T2C, with every bridging site vacant, is a
gel of dimers.

A gel of Ca/Si 1 in a dilute sodium hydroxide solution, with more and more
aluminum, shows what the aluminum does to the chains. Gibbsite, amorphous
aluminum hydroxide and the calcium aluminate hydrates are declared, so that the
aluminum the gel cannot take has somewhere to go:

```@example sublattice
pure_n = ["Portlandite", "Amor-Sl", "Gbs", "AlOHam", "C3AH6", "C4AH13", "straetlingite"]
cs_n = gel_system(CNASH, cnash; pure = pure_n)
println("Al/Si added  certified   gel Ca/Si  Al/Si  Na/Si   chain length    pH    other solids")
for alsi in (0.0, 0.05, 0.10, 0.15)
    st = ChemicalState(cs_n)
    set_quantity!(st, "Amor-Sl", 0.05u"mol")
    set_quantity!(st, "Portlandite", 0.05u"mol")
    alsi > 0 && set_quantity!(st, "Gbs", alsi * 0.05u"mol")
    set_quantity!(st, "Na+", 0.02u"mol")
    set_quantity!(st, "OH-", 0.02u"mol")
    set_quantity!(st, "H2O@", 1.0u"kg")
    b = budget(st)
    eq, cert = equilibrate_certified(st; model, b)
    gel = solid_solution_totals(eq, "C-S-H")
    χ = [gel.members[m] for m in CNASH] ./ gel.amount
    CL = 3 / sum(χ[k] * ν[CNASH[k]] for k in eachindex(CNASH)) - 1
    el = gel.elements
    n = ustrip.(us"mol", eq.n)
    others = [symbol(cs_n.species[i]) for i in cs_n.idx_crystal
              if n[i] > 1.0e-6 && !(symbol(cs_n.species[i]) in CNASH)]
    @printf("   %4.2f       %-9s   %5.3f     %5.3f  %5.3f     %5.2f       %6.3f   %s\n", alsi, cert.optimal,
            el[:Ca] / el[:Si], get(el, :Al, 0.0) / el[:Si], get(el, :Na, 0.0) / el[:Si], CL,
            pH(eq, model), isempty(others) ? "none" : join(others, ", "))
end
```

Without aluminum, the gel of Ca/Si 1 has a chain length of 5.2 tetrahedra. Each
step of aluminum lengthens it, to 7.7 at an Al/Si of 0.125, since the members that
take the aluminum have fewer vacant bridging sites than T2C. The gel also takes
sodium, about 0.18 per silicon. At this Ca/Si it holds no more aluminum than
0.125 per silicon: of the 0.15 added, the rest precipitates as gibbsite (`Gbs`).
Every paste certifies. The chain length is a lower bound, since Eq. (11) counts
the chains as not crosslinked, as [Myers2014](@citet) state.

## Where to go next

The derivation of the site model, its convexity and what the solver does with
it are in [Solid solutions](@ref sec-theory-solid-solutions).
[The CEM IV page](@ref ex-cem4-pozzolanic) uses the CNASH gel in a cement paste,
beside the `CSHQ` model of the same gel, and [The models of the C-S-H gel](@ref
sec-csh-models) lists the six models the package ships. [The CASH+ page](@ref
ex-cashplus-csh) adds to site mixing the energies of the compounds, in the model of
[Kulik2022](@citet).
