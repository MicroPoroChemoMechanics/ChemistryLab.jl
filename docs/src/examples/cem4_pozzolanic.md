# [A pozzolanic binder, and the C-S-H that has to carry the aluminum](@id ex-cem4-pozzolanic)

A CEM IV replaces 11 % to 55 % of the clinker with a pozzolana — siliceous fly
ash, natural or calcined pozzolana, silica fume, or a mixture of them. The
replacement is not a filler: a pozzolana **consumes portlandite** and makes more
C-S-H out of it, which is why a pozzolanic binder is specified where the paste
must be denser and less alkaline.

Two things follow for the calculation, and this page is about both.

1. **The aluminum has somewhere to go, and the model must let it.** Fly ash
   brings almost as much Al as Si. In a real paste most of that aluminum ends up
   *inside* the calcium silicate hydrate, as C-A-S-H. `CSHQ`, the C-S-H model
   the CEM I pages use, has **no aluminum end-member at all**, so a calculation
   that keeps it forces every atom of aluminum into the AFm and AFt phases.
   `CNASH_ss` is the model that can take it.
2. **Portlandite becomes the limiting reagent.** Past a certain replacement it
   runs out, and what the paste can do afterwards changes.

```@example cem4
using ChemistryLab
using DynamicQuantities
using OptimaSolver
using OrderedCollections
using Printf
using Plots
default(framestyle = :box, grid = false)

substances = build_species(datapath("cemdata18-thermofun.json"); verbose = false)
byname = Dict(symbol(s) => s for s in substances)
molar_mass(n) = ustrip(us"g/mol", byname[n][:M])
nothing # hide
```

## 1. What is assumed here, and why all of it is

!!! warning "This page has no measured specimen behind it"
    The calorimetry deposit this package ships [Smilauer2025data](@cite) carries
    CEM I, CEM II, CEM III and CEM V records, and **no CEM IV**. So unlike
    [the CEM III page](@ref ex-cem3-slag) and [the CEM II page](@ref ex-cem2-blended),
    nothing below is anchored to a specimen. It is a **model study of a family**,
    run on a composition chosen inside the EN 197-1 range, and the numbers
    describe what a pozzolanic binder of that composition does — not what any
    particular cement contains.

    That distinction is worth keeping, because the aluminum question this page is
    about does not depend on the exact analysis, while the assemblage does.

```@example cem4
# ASSUMED: a Bogue composition representative of a CEM I clinker.
CLINKER = OrderedDict("C3S" => 0.65, "C2S" => 0.11, "C3A" => 0.11, "C4AF" => 0.08)

# ASSUMED: a siliceous (class V) fly ash analysis of the kind European standards
# admit — low calcium, high silica and alumina, and the alkalis that make the
# aluminum question interesting.
FLYASH = OrderedDict("SiO2" => 0.53, "Al2O3" => 0.26, "Fe2O3" => 0.07,
                     "CaO" => 0.04, "MgO" => 0.02, "K2O" => 0.025,
                     "Na2O" => 0.008, "SO3" => 0.005)

# ASSUMED: the midpoint of the EN 197-1 range for a CEM IV/B, which is 45-64 %
# clinker and 36-55 % pozzolana.
ASH_FRACTION = 0.45
GYPSUM = 0.046
WB = 0.50
BINDER_G = 100.0
nothing # hide
```

Only the reactive part of a fly ash takes part in the chemistry — the glassy
fraction — and the crystalline mullite and quartz in it do not dissolve on any
relevant time scale. An equilibrium calculation has no way to distinguish them,
so it treats the whole analysis as available. **That overstates what the ash
contributes**, and it is stated here rather than corrected, because correcting it
would need a degree-of-reaction that only a kinetic description supplies.

## 2. Two species lists, differing in one phase

Everything is held fixed between the two systems except the C-S-H model. That is
what makes the comparison mean something:

```@example cem4
pure = split(
    "C3S C2S C3A C4AF Gp Anh Cal Portlandite ettringite monosulphate12 " *
    "monocarbonate hemicarbonate C4AH13 C3AH6 C3FH6 straetlingite " *
    "hydrotalcite Brc FeOOHmic AlOHmic Amor-Sl Mgs " *
    "K2SO4 syngenite Na2SO4"
)
aqueous = ["SO4-2", "CO2@", "O2@"]

CSHQ = ["CSHQ-JenD", "CSHQ-JenH", "CSHQ-TobD", "CSHQ-TobH", "KSiOH", "NaSiOH"]
CNASH = ["T2C-CNASHss", "T5C-CNASHss", "TobH-CNASHss",
         "5CA", "5CNA", "INFCA", "INFCN", "INFCNA"]
FEAL = ["C3AFS0.84H4.32", "C3FS0.84H4.32"]

function system(gel_name, gel_members)
    sp = speciation(substances, vcat(pure, gel_members, FEAL, aqueous);
                    aggregate_state = [AS_AQUEOUS])
    ss = [SolidSolutionPhase(gel_name, [byname[m] for m in gel_members]),
          SolidSolutionPhase("C3(AF)S0.84H", [byname[m] for m in FEAL])]
    return ChemicalSystem(sp, CEMDATA_PRIMARIES; solid_solutions = ss)
end

cs_q = system("CSHQ", CSHQ)
cs_n = system("CNASH_ss", CNASH)
model = HKFActivityModel(å = 0.0, Ḃ = 0.097637, Kₙ = 0.0)

@printf("CSHQ system     : %d species\n", length(cs_q.species))
@printf("CNASH_ss system : %d species\n", length(cs_n.species))
```

!!! danger "Never both at once"
    `CSHQ`, `CNASH_ss` and the `ECSH` family are three *models of one gel*, not
    three phases. Declaring two of them counts the same calcium silicate hydrate
    twice, and `ChemicalSystem` refuses the pair by name. The two systems
    above are alternatives, built separately and compared, which is the only
    correct way to use them.

## 3. The budget

```@example cem4
function budget(cs; ash, wb = WB)
    clinker_frac = 1 - ash - GYPSUM
    state = ChemicalState(cs)
    for (phase, frac) in CLINKER
        set_quantity!(state, phase,
            BINDER_G * clinker_frac * frac / molar_mass(phase) * u"mol")
    end
    set_quantity!(state, "Gp", BINDER_G * GYPSUM / molar_mass("Gp") * u"mol")
    set_quantity!(state, "H2O@", BINDER_G * wb / molar_mass("H2O@") * u"mol")
    b = Float64.(cs.SM.A) * ustrip.(us"mol", state.n)
    ash > 0 && (b .+= oxide_budget(FLYASH, cs.SM.primaries;
                                   mass = BINDER_G * ash * u"g"))
    return state, b
end

st_q, b_q = budget(cs_q; ash = ASH_FRACTION)
st_n, b_n = budget(cs_n; ash = ASH_FRACTION)

comps = String.(symbol.(cs_n.SM.primaries))
for (c, v) in zip(comps, b_n)
    abs(v) > 1.0e-6 && @printf("  %-8s %10.5f mol\n", c, v)
end
```

The silicon-to-aluminum ratio of that budget is the whole issue. A CEM I paste
has roughly seven times more Si than Al; this one has under three.

## 4. The same paste, two C-S-H models

```@example cem4
eq_q, c_q = equilibrate_certified(st_q; model = model, b = b_q)
eq_n, c_n = equilibrate_certified(st_n; model = model, b = b_n)

for (label, cs, eq, c) in (("CSHQ", cs_q, eq_q, c_q),
                           ("CNASH_ss", cs_n, eq_n, c_n))
    @printf("%-10s optimal=%-5s  worst SI=%+.2e  balance=%.1e  pH=%.3f  V=%.2f cm3\n",
            label, c.optimal, c.worst_supersaturation, c.balance,
            pH(eq, model), ustrip(uconvert(us"cm^3", volume(eq).total)))
end
```

Where does the aluminum end up? The conservation matrix already answers it. The
row of the primary species `AlO2-` counts one unit per aluminum atom, so
multiplying that row by the amounts distributes the element over the phases that
hold it:

```@example cem4
"""Moles of one conservation component held by each solid phase, largest first."""
function component_in_solids(cs, eq, component; tol = 1.0e-4)
    n = ustrip.(us"mol", eq.n)
    row = findfirst(==(component), String.(symbol.(cs.SM.primaries)))
    row === nothing && error("$component is not a component of this system")
    A = Float64.(cs.SM.A)
    out = [(symbol(cs.species[i]), A[row, i] * n[i]) for i in cs.idx_crystal]
    return sort(filter(p -> last(p) > tol, out); by = last, rev = true)
end

for (label, cs, eq) in (("CSHQ", cs_q, eq_q), ("CNASH_ss", cs_n, eq_n))
    println(label, " — where the aluminum is:")
    for (name, al) in component_in_solids(cs, eq, "AlO2-")
        @printf("  %-18s %9.5f mol Al\n", name, al)
    end
    println()
end
```

## 5. The assemblages, side by side

```@example cem4
function assemblage(cs, eq; tol = 1.0e-4)
    n = ustrip.(us"mol", eq.n)
    sort([(symbol(cs.species[i]), n[i]) for i in cs.idx_crystal if n[i] > tol];
         by = last, rev = true)
end

for (label, cs, eq) in (("CSHQ", cs_q, eq_q), ("CNASH_ss", cs_n, eq_n))
    println(label, ":")
    for (name, amount) in assemblage(cs, eq)
        @printf("  %-18s %9.5f mol\n", name, amount)
    end
    println()
end
```

!!! note "Read the C-A-S-H as a total, not as eight numbers"
    The eight `CNASH_ss` end-members span a space of rank 5: Myers' model carries
    site-occupancy constraints that an ideal eight-component mixture does not
    reproduce [Myers2014](@cite). The feasible set stays bounded and the solve is
    well posed, but the individual end-member amounts are not determined by the
    element balance alone — only their combinations are. The total, the Ca/Si and
    the Al/Si are the quantities to read.

## 6. Portlandite is the limiting reagent

The pozzolanic reaction consumes calcium hydroxide. Sweep the replacement level
and the point where it runs out is visible directly:

```@example cem4
fractions = 0.0:0.1:0.5
ch = Float64[]
phs = Float64[]
for f in fractions
    st, b = budget(cs_n; ash = f)
    eq, _ = equilibrate_certified(st; model = model, b = b)
    n = ustrip.(us"mol", eq.n)
    i = findfirst(s -> symbol(s) == "Portlandite", cs_n.species)
    push!(ch, n[i])
    push!(phs, pH(eq, model))
    @printf("  ash %3.0f %%   portlandite %8.5f mol   pH %.3f\n", 100f, n[i], phs[end])
end
```

```@example cem4
p1 = plot(100 .* collect(fractions), ch; marker = :circle, legend = false,
          xlabel = "fly ash (% of binder)", ylabel = "portlandite (mol / 100 g)",
          color = :seagreen, title = "Calcium hydroxide consumed")
p2 = plot(100 .* collect(fractions), phs; marker = :circle, legend = false,
          xlabel = "fly ash (% of binder)", ylabel = "pH",
          color = :steelblue, title = "Pore solution pH")
fig = plot(p1, p2; layout = (1, 2), size = (900, 380),
           bottom_margin = 10Plots.mm, left_margin = 10Plots.mm)
savefig(fig, "cem4-sweep.svg"); nothing # hide
```

![](cem4-sweep.svg)

Two things happen at once along that sweep, and they are worth separating.
The **portlandite** falls because the ash's silica turns it into more C-S-H —
that is the pozzolanic reaction, and it is the property the family is specified
for. The **pH** moves much less, because in a cement paste it is the alkalis that
set it, not the calcium hydroxide; portlandite only fixes a floor around 12.5 at
25 °C. A pozzolanic binder lowers the pH mainly by **binding alkalis into the
C-A-S-H**, and that is a mechanism only the `CNASH_ss` model can express at all.

## 7. The alkalis, and the phases that could take them

Section 6 left a loose end: the pH barely moves, and the reason given was that
the alkalis stay in solution. That is a statement about the **phase list**, not
about the chemistry — nothing in the species list above can hold potassium or
sodium except the C-A-S-H and two sulfates. A pozzolanic paste at high alkalinity
has another option, and it is the one CEMDATA18 does not carry:
[the zeolite extension](@ref sec-zeolites).

```@example cem4
zeo_db = build_species(datapath("cemdata18-zeolites.json"); verbose = false)
zeo_byname = Dict(symbol(s) => s for s in zeo_db)
zeolites = sort(collect(setdiff(Set(keys(zeo_byname)), Set(keys(byname)))))

# Only the ones this paste could form: its alkalis are K and Na, and it has no
# chloride or nitrate.
println(length(zeolites), " phases added by the extension:")
for z in zeolites
    print(z, "  ")
end
```

```@example cem4
pure_z = vcat(String.(pure), zeolites)
sp_z = speciation(zeo_db, vcat(pure_z, CNASH, FEAL, aqueous);
                  aggregate_state = [AS_AQUEOUS])
ss_z = [SolidSolutionPhase("CNASH_ss", [zeo_byname[m] for m in CNASH]),
        SolidSolutionPhase("C3(AF)S0.84H", [zeo_byname[m] for m in FEAL])]
cs_z = ChemicalSystem(sp_z, CEMDATA_PRIMARIES; solid_solutions = ss_z)

st_z = ChemicalState(cs_z)
for (phase, frac) in CLINKER
    set_quantity!(st_z, phase,
        BINDER_G * (1 - ASH_FRACTION - GYPSUM) * frac / molar_mass(phase) * u"mol")
end
set_quantity!(st_z, "Gp", BINDER_G * GYPSUM / molar_mass("Gp") * u"mol")
set_quantity!(st_z, "H2O@", BINDER_G * WB / molar_mass("H2O@") * u"mol")
b_z = Float64.(cs_z.SM.A) * ustrip.(us"mol", st_z.n)
b_z .+= oxide_budget(FLYASH, cs_z.SM.primaries;
                     mass = BINDER_G * ASH_FRACTION * u"g")

eq_z, c_z = equilibrate_certified(st_z; model = model, b = b_z)
@printf("with zeolites: optimal=%-5s  worst SI=%+.2e  pH=%.3f
",
        c_z.optimal, c_z.worst_supersaturation, pH(eq_z, model))
@printf("without      : optimal=%-5s  worst SI=%+.2e  pH=%.3f
",
        c_n.optimal, c_n.worst_supersaturation, pH(eq_n, model))

nz = ustrip.(us"mol", eq_z.n)
formed = sort([(symbol(cs_z.species[i]), nz[i])
               for i in cs_z.idx_crystal
               if nz[i] > 1.0e-6 && symbol(cs_z.species[i]) in zeolites];
              by = last, rev = true)
if isempty(formed)
    println("
no zeolite is stable in this paste")
else
    println("
zeolites formed:")
    for (name, amount) in formed
        @printf("  %-16s %9.5f mol
", name, amount)
    end
end
```

Whichever way that comes out, the calculation is now **able to answer the
question**, and before the extension it was not: a phase absent from the species
list is not reported as undersaturated, it is not reported at all. The
certificate's supersaturation figure is the one to read — it says whether
anything the system *could* form and did not is trying to.

!!! danger "An equilibrium at 55 % replacement is not a 28-day paste"
    Everything above is the state the paste *tends to*, with the whole ash taken
    as reactive. A real fly ash binder reaches a fraction of it in a month. The
    figure is a map of the family's limit, and the way to a date on the calendar
    is the coupled route of [the kinetic pages](@ref ex-ionic-opc), which needs a
    rate law for the ash — `waller` ships published parameters for fly ash, and
    that is the subject of the coupled runs rather than of this page.
