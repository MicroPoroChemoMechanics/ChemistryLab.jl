# [Two CEM II, and the two different things a replacement can do](@id ex-cem2-blended)

!!! info "Before this page"
    [A CEM I from its clinker phases](@ref sec-cem1-from-clinker) and [The
    binders, and what distinguishes them](@ref man-binder-families).

A CEM II replaces between 6 % and 35 % of the clinker. What that does to the
chemistry depends entirely on **what the replacement is made of**, and the two
most common ones sit at opposite ends:

- **limestone (L, LL)** is a *crystalline* phase with a formula. It enters the
  calculation as calcite, and it is very nearly inert as a filler — but the
  carbonate it carries is not inert at all. It rewrites the aluminate sequence.
- **blastfurnace slag (S)** is a *glass*. It has no formula, so it enters as an
  oxide analysis, and it brings magnesium and reduced sulfur along with the
  calcium, silicon and aluminum.

This page runs both against measured records, and isolates the carbonate effect
by removing the limestone from an otherwise identical paste.

## 0. Two words you need first: AFm and AFt

Everything below turns on where the **aluminum** and the **sulfate** of a cement
end up, and the two families that hold them have names worth learning once.

They are both calcium aluminate hydrates built on the same idea: layers of
``\mathrm{Ca_2Al(OH)_6^+}`` with anions and water between them. What differs is
how many sulfates one aluminum carries.

| family | full name | anion per Al | archetype | formula |
|:--|:--|:--|:--|:--|
| **AFt** | alumino-ferrite-**tri** | 3 | ettringite | ``\mathrm{Ca_6Al_2(SO_4)_3(OH)_{12}\cdot 26\,H_2O}`` |
| **AFm** | alumino-ferrite-**mono** | 1 | monosulfate | ``\mathrm{Ca_4Al_2(SO_4)(OH)_{12}\cdot 6\,H_2O}`` |

The AFm layer is not fussy about *which* single anion it holds: sulfate gives
monosulfate, carbonate gives **monocarbonate**, hydroxide gives ``\mathrm{C_4AH_{13}}``.
That indifference is the whole mechanism of this page.

A Portland cement is ground with 3–5 % gypsum, so its early hydration makes
ettringite, which holds three sulfates per aluminum:

```math
\mathrm{C_3A} + 3\,\mathrm{C\bar{S}H_2} + 26\,\mathrm{H}
  \;\longrightarrow\; \mathrm{C_6A\bar{S}_3H_{32}} \quad (\text{ettringite, AFt})
```

Once the gypsum runs out, the remaining ``\mathrm{C_3A}`` attacks the ettringite
it just made, and three sulfates are spread over three aluminums instead of one:

```math
\mathrm{C_6A\bar{S}_3H_{32}} + 2\,\mathrm{C_3A} + 4\,\mathrm{H}
  \;\longrightarrow\; 3\,\mathrm{C_4A\bar{S}H_{12}} \quad (\text{monosulfate, AFm})
```

**Unless there is carbonate.** Calcite supplies an anion the AFm layer likes
better, so the aluminate goes to monocarbonate instead — and the sulfate it did
*not* consume has nowhere to go but back into ettringite:

```math
\mathrm{C_3A} + \mathrm{C\bar{C}} + 11\,\mathrm{H}
  \;\longrightarrow\; \mathrm{C_4A\bar{C}H_{11}} \quad (\text{monocarbonate, AFm})
```

That is why a few percent of ground limestone is not a filler. Ettringite carries
26 waters to monosulfate's 12, so keeping the sulfate in the AFt phase binds more
water into solids and the paste occupies more volume — with less clinker.

!!! note "Cement chemist notation"
    ``\mathrm{C = CaO}``, ``\mathrm{A = Al_2O_3}``, ``\mathrm{S = SiO_2}``,
    ``\mathrm{H = H_2O}``, ``\bar{\mathrm{S}} = \mathrm{SO_3}``,
    ``\bar{\mathrm{C}} = \mathrm{CO_2}``. The bar marks an acidic oxide, which
    is how ``\bar{\mathrm{S}}`` (sulfate) is told from ``\mathrm{S}`` (silica).
    [The notation page](@ref man-cement-notation) has the full alphabet.

**Nothing below assumes any of this.** The calculation is a Gibbs energy
minimization over an element budget: it is told which phases exist and how much
of each element the paste contains, and it finds the assemblage of lowest energy.
The reactions above are what the answer will turn out to mean, not what it was
told.

```@example cem2
using ChemistryLab
using DynamicQuantities
using Logging
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

## 1. What is measured

Two records from the CC-BY-4.0 deposit of [Smilauer2025data](@cite), one of each
family:

```@example cem2
function header(file)
    for line in eachline(joinpath(datapath("experimental"), file))
        startswith(line, "#") || break
        occursin(r"cement:|blaine:|wb:|temperature:", line) && println(line)
    end
    println()
end

header("smilauer2025-165-cemII-A-LL-42.5R-hranice.csv")
header("smilauer2025-149-cemII-B-S-32.5R-mokra.csv")
```

The fineness, the water/binder ratio and the calorimetry are measured. **The
clinker phase composition and the replacement level are not**, for either
record, and nothing in the deposit lets them be recovered. They are assumed
below, at the midpoint of the EN 197-1 range for each family, and the assumption
is stated where it is made rather than buried in a preamble.

```@example cem2
# ASSUMED: a Bogue composition representative of a CEM I clinker. The deposit
# reports none.
CLINKER = OrderedDict("C3S" => 0.65, "C2S" => 0.11, "C3A" => 0.11, "C4AF" => 0.08)

# ASSUMED: midpoints of the EN 197-1 ranges. CEM II/A-LL is 80-94 % clinker
# with 6-20 % limestone; CEM II/B-S is 65-79 % clinker with 21-35 % slag.
LL_LIMESTONE = 0.13
BS_SLAG = 0.28

# ASSUMED: a European ground granulated blastfurnace slag analysis.
SLAG = OrderedDict("CaO" => 0.41, "SiO2" => 0.36, "Al2O3" => 0.11,
                   "MgO" => 0.08, "SO3" => 0.02)

# THE CLINKER'S ALKALIS, which Bogue does not account for and which set the pH.
# They are minor oxides, a fraction of a percent, sitting outside the four-phase
# decomposition -- and in a cement paste they are what fixes the pH, dissolving
# almost completely into the pore solution and staying there while the calcium is
# held down by portlandite at 12.5. Omitting them does not make the calculation
# conservative: it makes it report a portlandite floor as a pore solution.
# ASSUMED at a usual industrial level, as a fraction of the CLINKER mass.
# [The CEM III page](@ref cem3-alkali) sweeps this over the industrial range and
# shows that it, and almost nothing else, is what moves the pH.
ALKALIS = OrderedDict("K2O" => 0.008, "Na2O" => 0.002)

# The calcium sulfate ground in with every Portland clinker, as a mass fraction
# of the binder. ASSUMED at a usual industrial level: the 4.6 % of the CEM I of
# [Lavergne2018](@cite), Table 9.
GYPSUM = literature_value("Lavergne2018", "gypsum_percent") / 100

BINDER_G = 100.0

# HOW MUCH REACTS, and it is not everything. Two ceilings; the reacted fraction
# is the lower of them.
#
# The WATER ceiling, Powers (1948): complete hydration needs about 0.42 g of
# water per gram of cement, 0.23 g written into the hydrates and 0.19 g held in
# the gel pores those hydrates create. Below that the paste stops with water
# still in it, in pores too fine to reach an unhydrated grain. It is a property
# of the pore space, so it caps the slag as it caps the clinker; and it depends
# on the mix, which is why the two pastes here do not get the same number --
# the limestone paste was mixed at w/b 0.45 and the slag paste at 0.40.
#
# The KINETIC ceiling, for the slag only: [Durdzinski2017](@cite), Table 5, a
# ground granulated slag at 28 days by SEM image analysis, 38-49 % across two
# slags and two laboratories, with a stated precision of "at best +/- 5 %".
# ASSUMED at their mean, 45 %. The limestone needs none: calcite is a declared
# phase, and the minimization dissolves exactly as much of it as is stable.
# Degrees of reaction measured by SEM image analysis on sealed pastes, in
# percent: [Durdzinski2017](@cite), Table 5, from data/literature/Durdzinski2017.json.
sem(material, age) = literature_table("Durdzinski2017", "degree_of_reaction";
    technique = "SEM-IA", material, curing = "sealed", age_days = age).degree_percent
mean_percent(x) = sum(x) / length(x)
ALPHA_SLAG_KINETIC = mean_percent([sem("S1", 28); sem("S2", 28)]) / 100
reacted(wb) = (clinker = powers_alpha_max(wb),
               slag = min(ALPHA_SLAG_KINETIC, powers_alpha_max(wb)))

for wb in (0.45, 0.40)
    r = reacted(wb)
    @printf("w/b = %.2f : water ceiling %.3f -> clinker %.0f %%, slag %.0f %%\n",
            wb, powers_alpha_max(wb), 100r.clinker, 100r.slag)
end
```

## 2. One species list for all three pastes

The comparison is only honest if the three calculations are allowed to form the
same phases. So the species list is built once, with the carbonate AFm phases,
the sulfate AFm and AFt phases, hydrotalcite for the magnesium, and the sulfur
ladder the slag needs:

```@example cem2
pure = split(
    "C3S C2S C3A C4AF Gp Anh Cal Portlandite ettringite monosulphate12 " *
    "monocarbonate hemicarbonate monocarbonate9 hemicarbonat10.5 " *
    "C4AH13 C3AH6 C3FH6 straetlingite Femonocarbonate Fe-hemicarbonate " *
    "hydrotalcite Mg2AlC0.5OH Brc FeOOHmic AlOHmic Amor-Sl Tro Py Sulfur Mgs"
)
solutions = [
    "CSHQ" => ["CSHQ-JenD", "CSHQ-JenH", "CSHQ-TobD", "CSHQ-TobH", "KSiOH", "NaSiOH"],
    "C3(AF)S0.84H" => ["C3AFS0.84H4.32", "C3FS0.84H4.32"],
]
members = reduce(vcat, last.(solutions))
redox_species = ["HS-", "H2S@", "SO4-2", "SO3-2", "S2O3-2", "O2@", "H2@", "CO2@"]

species = speciation(substances, vcat(pure, members, redox_species);
                     aggregate_state = [AS_AQUEOUS])
ss = [SolidSolutionPhase(n, [byname[m] for m in ms]) for (n, ms) in solutions]
cs = ChemicalSystem(species, CEMDATA_PRIMARIES; solid_solutions = ss)
# Debye-Hückel limiting law with a B-dot term, as GEM-Selektor runs CEMDATA18.
# The B-dot is identified from the activity coefficients GEMS printed on a
# Portland paste (test/reference/gems_cemdata18_portland.json), about 0.0976.
using JSON
gems = JSON.parsefile(joinpath(pkgdir(ChemistryLab), "test", "reference", "gems_cemdata18_portland.json"))
lg1, lg2 = log10(gems["gamma"]["z1"]), log10(gems["gamma"]["z2"])
Ḃ_gems = (lg1 + (lg1 - lg2) / 3) / gems["ionic_strength_mol_per_kg"]
model = HKFActivityModel(å = 0.0, Ḃ = Ḃ_gems, Kₙ = 0.0)

components = String.(symbol.(cs.SM.primaries))
@printf("%d species, %d conservation components: %s\n",
        length(cs.species), length(components), join(components, " "))
```

`Zz` is in that list because the sulfur ladder is: with sulfur at more than one
oxidation state, charge stops being a fixed combination of the element rows and
becomes a conserved quantity of its own. That is explained in
[the redox chapter](@ref theory-redox); here it simply means all three pastes
are posed on the same components, so their budgets are comparable term by term.

## 3. The three budgets

Each paste is 100 g of binder. The clinker enters through its phases, the
limestone as calcite — a species with a formula — and the slag through
[`oxide_budget`](@ref), which is what a glass needs:

```@example cem2
"""Element budget of one paste, in moles per 100 g of binder."""
function budget(; clinker_frac, limestone = 0.0, slag = 0.0, wb)
    α = reacted(wb)
    state = ChemicalState(cs)
    # Only the reacted clinker is posed to the minimization; the unhydrated
    # cores are still in the specimen but are not at equilibrium with it.
    for (phase, frac) in CLINKER
        set_quantity!(state, phase,
            α.clinker * BINDER_G * clinker_frac * frac / molar_mass(phase) * u"mol")
    end
    # The calcium sulfate is soluble and the limestone is a declared phase, so
    # neither carries a ceiling. All of the mixing water enters: the ceiling
    # says how far the reaction goes, not how much water was poured in.
    set_quantity!(state, "Gp", BINDER_G * GYPSUM / molar_mass("Gp") * u"mol")
    limestone > 0 && set_quantity!(state, "Cal",
        BINDER_G * limestone / molar_mass("Cal") * u"mol")
    set_quantity!(state, "H2O@", BINDER_G * wb / molar_mass("H2O@") * u"mol")

    b = Float64.(cs.SM.A) * ustrip.(us"mol", state.n)
    # The alkalis leave the grain as it dissolves, so they follow the clinker and
    # its reacted fraction.
    b .+= oxide_budget(ALKALIS, cs.SM.primaries;
                       mass = BINDER_G * clinker_frac * α.clinker * u"g")
    slag > 0 && (b .+= oxide_budget(SLAG, cs.SM.primaries;
                                    mass = BINDER_G * slag * α.slag * u"g"))
    return state, b
end

# The clinker fraction is what is left once the replacement and the gypsum are
# taken out, so the three pastes really are 100 g of binder each.
st_ll, b_ll = budget(clinker_frac = 1 - LL_LIMESTONE - GYPSUM,
                     limestone = LL_LIMESTONE, wb = 0.45)
st_ref, b_ref = budget(clinker_frac = 1 - LL_LIMESTONE - GYPSUM, wb = 0.45)
st_bs, b_bs = budget(clinker_frac = 1 - BS_SLAG - GYPSUM,
                     slag = BS_SLAG, wb = 0.40)

@printf("%-10s %s\n", "", join((@sprintf("%8s", c) for c in components), ""))
for (label, b) in ("CEM II/A-LL" => b_ll, "no limestone" => b_ref,
                   "CEM II/B-S" => b_bs)
    @printf("%-12s%s\n", label, join((@sprintf("%8.3f", v) for v in b), ""))
end
```

The middle row is not a cement anyone sells. It is the CEM II/A-LL with its
13 g of calcite removed and nothing put back — the *control*, which isolates
what the carbonate does from what the dilution does.

## 4. The three equilibria

```@example cem2
function solve(state, b)
    eq, cert = equilibrate_certified(state; model = model, b = b)
    return eq, cert
end

eq_ll, c_ll = solve(st_ll, b_ll)
eq_ref, c_ref = solve(st_ref, b_ref)
eq_bs, c_bs = solve(st_bs, b_bs)

for (label, eq, c) in (("CEM II/A-LL", eq_ll, c_ll),
                       ("no limestone", eq_ref, c_ref),
                       ("CEM II/B-S", eq_bs, c_bs))
    @printf("%-13s optimal=%-5s  worst SI=%+.2e  balance=%.1e  pH=%.3f  V=%.2f cm3\n",
            label, c.optimal, c.worst_supersaturation, c.balance,
            pH(eq, model), ustrip(uconvert(us"cm^3", volume(eq).total)))
end
```

## 5. What the carbonate does — the aluminate sequence

This is the reason a limestone addition is not a dilution. Put the two
aluminate-bearing assemblages side by side:

```@example cem2
function assemblage(eq; tol = 1.0e-4)
    n = ustrip.(us"mol", eq.n)
    sort([(symbol(cs.species[i]), n[i]) for i in cs.idx_crystal if n[i] > tol];
         by = last, rev = true)
end

for (label, eq) in ("CEM II/A-LL" => eq_ll, "no limestone" => eq_ref)
    println(label, ":")
    for (name, amount) in assemblage(eq)
        @printf("  %-18s %9.5f mol\n", name, amount)
    end
    println()
end
```

The sulfate has nowhere else to go when there is no carbonate: it ends up in
`monosulphate12`, the AFm phase. Add calcite and the carbonate takes the AFm
site instead — `monocarbonate`, and `hemicarbonate` at lower carbonate
availability — which leaves the sulfate to stay in **ettringite**. The paste
therefore holds more ettringite, and ettringite is the most voluminous and most
water-rich of the aluminate hydrates, so the solid volume goes up.

That mechanism is the one Lothenbach and co-workers established
[Matschei2007](@cite); this calculation reproduces it from the element budget
alone, with nothing fitted.

```@example cem2
aluminates = ["ettringite", "monosulphate12", "monocarbonate", "hemicarbonate",
              "C4AH13", "straetlingite"]
n_ll = ustrip.(us"mol", eq_ll.n)
n_ref = ustrip.(us"mol", eq_ref.n)
idx = Dict(symbol(s) => i for (i, s) in enumerate(cs.species))
get_n(n, s) = haskey(idx, s) ? n[idx[s]] : 0.0

v_ref = [get_n(n_ref, s) for s in aluminates]
v_ll = [get_n(n_ll, s) for s in aluminates]
top = 1.1 * maximum(vcat(v_ref, v_ll))

p_ref = bar(aluminates, v_ref; legend = false, color = :steelblue, ylims = (0, top),
            ylabel = "mol per 100 g of binder", xrotation = 30,
            title = "no limestone")
p_ll = bar(aluminates, v_ll; legend = false, color = :seagreen, ylims = (0, top),
           xrotation = 30, title = "CEM II/A-LL, 13 % calcite")
fig = plot(p_ref, p_ll; layout = (1, 2), size = (900, 420),
           bottom_margin = 16Plots.mm, left_margin = 8Plots.mm,
           plot_title = "The carbonate moves the sulfate out of the AFm and into the AFt",
           plot_titlefontsize = 11)
savefig(fig, "cem2-aluminates.svg"); nothing # hide
```

![](cem2-aluminates.svg)

## 6. What the slag does — magnesium, and a sulfur it cannot use

The slag paste is a different problem. The aluminum arrives with no sulfate of
its own, the magnesium has no clinker counterpart, and the sulfur arrives
**reduced**:

```@example cem2
println("CEM II/B-S:")
for (name, amount) in assemblage(eq_bs)
    @printf("  %-18s %9.5f mol\n", name, amount)
end

# `pe` warns when a member of the couple sits at the solver's floor; the two
# sulfur amounts printed below say the same, so the warning is not repeated.
pe_bs, eh_bs = with_logger(NullLogger()) do
    pe(eq_bs, model), Eh(eq_bs, model)
end
@printf("\npe  = %+.2f\nEh  = %+.3f V\n", pe_bs, eh_bs)
@printf("aqueous S(VI)  : %.3e mol\naqueous S(-II) : %.3e mol\n",
        ustrip(us"mol", moles(eq_bs, "SO4-2")), ustrip(us"mol", moles(eq_bs, "HS-")))
```

Read the potential as [the CEM III page](@ref ex-cem3-slag) reads it: a couple
buffers a potential only while both of its members are actually present in
solution, and `pe` warns rather than returning a number silently when one of them
is at the solver's floor. Compare the two aqueous sulfur figures printed above
before reading the potential as a property of the paste — the smaller they are,
the more the number is set by the floor and the less it is set by the chemistry.

The magnesium is the visible difference. It has no place in any calcium
aluminate hydrate, so it leaves the solution as **hydrotalcite** — a phase that
simply does not appear in a CEM I or a CEM II/A-LL calculation because the
element is not there to form it.

## 7. All three together, and the measured heat

```@example cem2
labels = ["CEM II/A-LL", "no limestone", "CEM II/B-S"]
eqs = [eq_ll, eq_ref, eq_bs]
vols = [ustrip(uconvert(us"cm^3", volume(e).total)) for e in eqs]
phs = [pH(e, model) for e in eqs]

p1 = bar(labels, vols; legend = false, ylabel = "total volume (cm³)",
         color = :seagreen, title = "Volume of the hydrated paste", xrotation = 15)
p2 = bar(labels, phs; legend = false, ylabel = "pH", ylims = (12.0, 13.6),
         color = :steelblue, title = "Pore solution pH", xrotation = 15)
fig2 = plot(p1, p2; layout = (1, 2), size = (900, 400),
            bottom_margin = 14Plots.mm, left_margin = 8Plots.mm)
savefig(fig2, "cem2-summary.svg"); nothing # hide
```

![](cem2-summary.svg)

The equilibrium says what the paste tends to. What it cannot say is how fast it
gets there, and that is exactly where the two replacements differ most:

```@example cem2
function final_heat(file)
    t, Q = 0.0, 0.0
    for line in eachline(joinpath(datapath("experimental"), file))
        (startswith(line, "#") || startswith(line, "time")) && continue
        f = split(line, ',')
        t, Q = parse(Float64, f[1]), parse(Float64, f[3])
    end
    return t, Q
end

for (file, label) in (
        ("smilauer2025-122-cemI-52.5R-cizkovice.csv", "CEM I 52.5 R"),
        ("smilauer2025-165-cemII-A-LL-42.5R-hranice.csv", "CEM II/A-LL 42.5 R"),
        ("smilauer2025-149-cemII-B-S-32.5R-mokra.csv", "CEM II/B-S 32.5 R"),
    )
    t, Q = final_heat(file)
    @printf("  %-20s %6.0f h   %6.1f J/g\n", label, t, Q)
end
```

Both replacements lower the heat, because every joule comes from the clinker and
neither calcite nor slag glass releases any of its own at this age. The slag
lowers it further than the limestone does, at a *larger* replacement — 28 %
against 13 % — so per unit of clinker removed the two are not far apart. What
separates them is that the limestone's contribution to the solid volume is
mostly the calcite itself plus the ettringite it stabilizes, while the slag's
contribution grows over months as the glass dissolves, which an equilibrium
calculation at a single instant cannot show at all.

!!! note "The 28-day heat is not the final heat"
    The records stop at about three weeks. A slag binder is still reacting
    there, so its curve is truncated much further from its asymptote than the
    CEM I's is. Comparing final values of truncated curves understates the slag
    binder, and the deposit's own fitted `DoHInf` of 0.85 is a reminder that the
    depositors extrapolated rather than measured the end.

## Where to go next

A binder in which the slag is the main constituent, and its sulfur can no longer
be taken as sulfate, is [A blastfurnace cement, and the oxidation state it needs](@ref ex-cem3-slag).
