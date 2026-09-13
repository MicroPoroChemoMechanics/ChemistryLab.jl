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

## 0. What a pozzolana actually does

A pozzolana is not a binder on its own. Mix siliceous fly ash with water and
nothing happens. What makes it work is that a Portland clinker, while hydrating,
produces a large amount of **portlandite** — calcium hydroxide, ``\mathrm{Ca(OH)_2}``
— which is a by-product of the silicates:

```math
2\,\mathrm{C_3S} + 6\,\mathrm{H} \;\longrightarrow\; \mathrm{C_3S_2H_3} + 3\,\mathrm{CH}
```

Roughly a fifth of a hydrated CEM I paste is portlandite, and it contributes
little strength while being the phase most easily leached. The **pozzolanic
reaction** puts it to work: the ash's amorphous silica consumes it and makes more
of the phase that does carry the strength, the calcium silicate hydrate.

```math
\mathrm{S} + 1.7\,\mathrm{CH} + \text{water} \;\longrightarrow\; \mathrm{C_{1.7}SH_x}
```

So a pozzolanic binder trades a weak, soluble phase for a strong one, at the cost
of a slower reaction. Section 6 measures exactly that trade, by sweeping the
replacement level and watching the portlandite go.

But the ash brings something else that the equation above ignores, and it is what
the rest of the page is about: **aluminum**, almost as much as silicon. Where it
goes decides the answer, and it depends entirely on which C-S-H model the
calculation is given.

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

# ASSUMED: the midpoint of the EN 197-1 range for a CEM IV/A, which is 65-89 %
# clinker and 11-35 % pozzolana. Section 7 goes to a CEM IV/B, and shows what
# has to be added to the phase list before that is a question with an answer.
ASH_FRACTION = 0.23
ASH_FRACTION_B = 0.45
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
    # Aluminum sinks that CEMDATA18 documents and an earlier version of this
    # list simply did not declare. The siliceous hydrogarnet `C3AS0.84H4.32`
    # is the one that matters most here: it is the ALUMINUM end-member of the
    # family whose iron end-member was already present, and a blended binder
    # puts a great deal of aluminum into it. Leaving it out does not make the
    # calculation conservative -- it makes it insoluble, because the element
    # has to go somewhere.
    "C3AS0.84H4.32 C3AS0.41H5.18 straetlingite7 Gbs AlOHam " *
    "M4A-OH-LDH M6A-OH-LDH M8A-OH-LDH C2AH7.5 C4AH11 C4AH19 " *
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
fractions = 0.0:0.05:0.30
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

## 7. Push to a CEM IV/B, and the calculation stops having an answer

Everything so far was a CEM IV/**A**, 23 % ash. Take it to a CEM IV/**B** — the
midpoint of 36–55 % — and the equilibrium **fails to certify**, on either C-S-H
model. That is not a numerical accident and it is not a defect of `CNASH_ss`:

```@example cem4
eq_b, c_b = nothing, nothing          # the CNASH_ss case, kept for section 7
for (label, cs) in ("CSHQ" => cs_q, "CNASH_ss" => cs_n)
    st, b = budget(cs; ash = ASH_FRACTION_B)
    eq, c = equilibrate_certified(st; model = model, b = b)
    label == "CNASH_ss" && (global eq_b, c_b = eq, c)
    @printf("%-10s at %2.0f %% ash: optimal=%-5s  balance=%.1e  pH=%.3f\n",
            label, 100ASH_FRACTION_B, c.optimal, c.balance, pH(eq, model))
end
```

**The element budget has nowhere to put what the ash brings.** At 23 % the
aluminum fits in the C-A-S-H, the AFm/AFt phases and strätlingite, and the
alkalis fit in the C-A-S-H and the sulfates. At 45 % there is more aluminum and
more alkali than those phases can hold, portlandite is gone so the calcium
potential is no longer buffered by a pure phase, and the minimization is looking
for an assemblage that the declared phase list cannot form.

A real cement does not have this problem, because a real alkaline aluminosilicate
paste precipitates **zeolites** — and CEMDATA18 carries five of them, none of the
families this binder needs. That is exactly what
[the zeolite extension](@ref sec-zeolites) was built for.

```@example cem4
zeo_db = build_species(datapath("cemdata18-zeolites.json"); verbose = false)
zeo_byname = Dict(symbol(s) => s for s in zeo_db)
added = sort(collect(setdiff(Set(keys(zeo_byname)), Set(keys(byname)))))

# Only the ones this paste could form. Two of the twenty-eight are a chloride
# and a nitrate sodalite, and this binder carries neither element: declaring them
# would widen the species list to every aqueous chloride and nitrate species in
# the database, all of them on a budget of exactly zero, for no phase that can
# form. Dropping them is not a modeling choice, it is arithmetic.
carries(sp, el) = haskey(atoms(sp), Symbol(el))
zeolites = [z for z in added
            if !carries(zeo_byname[z], "Cl") && !carries(zeo_byname[z], "N")]

@printf("%d phases added by the extension, %d of them usable here\n",
        length(added), length(zeolites))
println("left out (no Cl and no N in this binder): ",
        join(setdiff(added, zeolites), ", "))
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
        BINDER_G * (1 - ASH_FRACTION_B - GYPSUM) * frac / molar_mass(phase) * u"mol")
end
set_quantity!(st_z, "Gp", BINDER_G * GYPSUM / molar_mass("Gp") * u"mol")
set_quantity!(st_z, "H2O@", BINDER_G * WB / molar_mass("H2O@") * u"mol")
b_z = Float64.(cs_z.SM.A) * ustrip.(us"mol", st_z.n)
b_z .+= oxide_budget(FLYASH, cs_z.SM.primaries;
                     mass = BINDER_G * ASH_FRACTION_B * u"g")

eq_z, c_z = equilibrate_certified(st_z; model = model, b = b_z)
@printf("with zeolites: optimal=%-5s  worst SI=%+.2e  pH=%.3f
",
        c_z.optimal, c_z.worst_supersaturation, pH(eq_z, model))
@printf("without      : optimal=%-5s  worst SI=%+.2e  pH=%.3f
",
        c_b.optimal, c_b.worst_supersaturation, pH(eq_b, model))

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
