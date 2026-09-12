# [Two CEM II, and the two different things a replacement can do](@id ex-cem2-blended)

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

```@example cem2
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

# The calcium sulfate ground in with every Portland clinker, as a mass fraction
# of the binder. ASSUMED at a usual industrial level.
GYPSUM = 0.046

BINDER_G = 100.0
nothing # hide
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
model = HKFActivityModel(å = 0.0, Ḃ = 0.097637, Kₙ = 0.0)

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
    state = ChemicalState(cs)
    for (phase, frac) in CLINKER
        set_quantity!(state, phase,
            BINDER_G * clinker_frac * frac / molar_mass(phase) * u"mol")
    end
    set_quantity!(state, "Gp", BINDER_G * GYPSUM / molar_mass("Gp") * u"mol")
    limestone > 0 && set_quantity!(state, "Cal",
        BINDER_G * limestone / molar_mass("Cal") * u"mol")
    set_quantity!(state, "H2O@", BINDER_G * wb / molar_mass("H2O@") * u"mol")

    b = Float64.(cs.SM.A) * ustrip.(us"mol", state.n)
    slag > 0 && (b .+= oxide_budget(SLAG, cs.SM.primaries;
                                    mass = BINDER_G * slag * u"g"))
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

@printf("\npe  = %+.2f\nEh  = %+.3f V\n", pe(eq_bs, model), Eh(eq_bs, model))
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
