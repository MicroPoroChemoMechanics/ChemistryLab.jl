# [A composite binder: two glasses at once](@id ex-cem5-composite)

A CEM V carries **both** a blastfurnace slag and a pozzolana, each between 18 %
and 30 % for a CEM V/A, leaving 40 % to 64 % clinker. It is the binder in which
every difficulty of the preceding pages arrives together:

- two constituents with **no formula**, so the budget comes from two oxide
  analyses rather than from phases ([the CEM III page](@ref ex-cem3-slag));
- **magnesium**, from the slag, which makes hydrotalcite and nothing else;
- **aluminum in excess of what the aluminate phases can hold**, from the fly
  ash, which is why the C-S-H has to be the one that carries it
  ([the CEM IV page](@ref ex-cem4-pozzolanic));
- **sulfur at two oxidation states**, S(-II) from the slag against S(+VI) from
  the clinker's calcium sulfate, which is what makes charge a conserved quantity
  of its own ([the redox chapter](@ref theory-redox)).

None of these is new here. What is new is that they must hold simultaneously,
and the point of the page is that the same element budget answers all four.

```@example cem5
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

## 1. The record, and the three numbers in it

```@example cem5
rec = joinpath(datapath("experimental"),
               "smilauer2025-200-cemV-A-S-V-32.5R-prachovice.csv")
for line in eachline(rec)
    startswith(line, "#") || break
    occursin(r"cement:|blaine:|wb:|temperature:", line) && println(line)
end
```

`(S-V)` in the designation is the composition: **S** for blastfurnace slag, **V**
for siliceous fly ash. The deposit gives the fineness, the water/binder ratio and
the calorimetry. It gives neither the clinker composition nor the two replacement
levels, so those are assumed at the midpoint of the EN 197-1 range and labeled:

```@example cem5
# ASSUMED: a Bogue composition representative of a CEM I clinker.
CLINKER = OrderedDict("C3S" => 0.65, "C2S" => 0.11, "C3A" => 0.11, "C4AF" => 0.08)

# ASSUMED: midpoint of the EN 197-1 CEM V/A range, which is 40-64 % clinker
# with 18-30 % slag and 18-30 % pozzolana.
SLAG_FRACTION = 0.24
ASH_FRACTION = 0.24
GYPSUM = 0.046

# ASSUMED: the same two analyses used on the CEM III and CEM IV pages, so that
# the three calculations differ in their proportions and not in their inputs.
SLAG = OrderedDict("CaO" => 0.41, "SiO2" => 0.36, "Al2O3" => 0.11,
                   "MgO" => 0.08, "SO3" => 0.02)
FLYASH = OrderedDict("SiO2" => 0.53, "Al2O3" => 0.26, "Fe2O3" => 0.07,
                     "CaO" => 0.04, "MgO" => 0.02, "K2O" => 0.025,
                     "Na2O" => 0.008, "SO3" => 0.005)

WB = 0.40            # MEASURED, from the record above
BINDER_G = 100.0
nothing # hide
```

## 2. The system: C-A-S-H, hydrotalcite and the sulfur ladder

```@example cem5
pure = split(
    "C3S C2S C3A C4AF Gp Anh Cal Portlandite ettringite monosulphate12 " *
    "monocarbonate hemicarbonate C4AH13 C3AH6 C3FH6 straetlingite " *
    "hydrotalcite Mg2AlC0.5OH Brc FeOOHmic AlOHmic Amor-Sl Mgs " *
    "Tro Py Sulfur K2SO4 syngenite Na2SO4"
)
gel = ["T2C-CNASHss", "T5C-CNASHss", "TobH-CNASHss",
       "5CA", "5CNA", "INFCA", "INFCN", "INFCNA"]
feal = ["C3AFS0.84H4.32", "C3FS0.84H4.32"]
redox_species = ["HS-", "H2S@", "SO4-2", "SO3-2", "S2O3-2", "O2@", "H2@", "CO2@"]

species = speciation(substances, vcat(pure, gel, feal, redox_species);
                     aggregate_state = [AS_AQUEOUS])
ss = [SolidSolutionPhase("CNASH_ss", [byname[m] for m in gel]),
      SolidSolutionPhase("C3(AF)S0.84H", [byname[m] for m in feal])]
cs = ChemicalSystem(species, CEMDATA_PRIMARIES; solid_solutions = ss)
model = HKFActivityModel(å = 0.0, Ḃ = 0.097637, Kₙ = 0.0)

components = String.(symbol.(cs.SM.primaries))
@printf("%d species, %d components: %s\n",
        length(cs.species), length(components), join(components, " "))
@printf("charge (`Zz`) kept as a conservation component: %s\n",
        "Zz" in components ? "yes" : "no")
```

## 3. The budget is additive

Three contributions, and the only reason they can be added is that they are
expressed on the **same components**: the clinker through its phases, each glass
through [`oxide_budget`](@ref).

```@example cem5
clinker_frac = 1 - SLAG_FRACTION - ASH_FRACTION - GYPSUM
state = ChemicalState(cs)
for (phase, frac) in CLINKER
    set_quantity!(state, phase,
        BINDER_G * clinker_frac * frac / molar_mass(phase) * u"mol")
end
set_quantity!(state, "Gp", BINDER_G * GYPSUM / molar_mass("Gp") * u"mol")
set_quantity!(state, "H2O@", BINDER_G * WB / molar_mass("H2O@") * u"mol")

b_clinker = Float64.(cs.SM.A) * ustrip.(us"mol", state.n)
b_slag = oxide_budget(SLAG, cs.SM.primaries;
                      mass = BINDER_G * SLAG_FRACTION * u"g")
b_ash = oxide_budget(FLYASH, cs.SM.primaries;
                     mass = BINDER_G * ASH_FRACTION * u"g")
b = b_clinker .+ b_slag .+ b_ash

@printf("clinker %.1f %%, slag %.1f %%, fly ash %.1f %%, gypsum %.1f %%\n\n",
        100clinker_frac, 100SLAG_FRACTION, 100ASH_FRACTION, 100GYPSUM)
@printf("%-8s %10s %10s %10s %10s\n", "", "clinker", "slag", "fly ash", "total")
for (i, c) in enumerate(components)
    abs(b[i]) > 1.0e-6 &&
        @printf("%-8s %10.5f %10.5f %10.5f %10.5f\n",
                c, b_clinker[i], b_slag[i], b_ash[i], b[i])
end
```

That table is the page in one object. The aluminum comes mostly from the ash,
the magnesium only from the slag, the sulfur from all three at two different
oxidation states, and the calcium overwhelmingly from the clinker even at 48 %
replacement.

## 4. The equilibrium

```@example cem5
eq, cert = equilibrate_certified(state; model = model, b = b)

@printf("certificate: optimal=%s  worst SI=%.2e  element balance=%.1e\n",
        cert.optimal, cert.worst_supersaturation, cert.balance)
@printf("pH = %.3f   total volume = %.2f cm3\n",
        pH(eq, model), ustrip(uconvert(us"cm^3", volume(eq).total)))
```

```@example cem5
n = ustrip.(us"mol", eq.n)
present = sort([(symbol(cs.species[i]), n[i])
                for i in cs.idx_crystal if n[i] > 1.0e-4]; by = last, rev = true)
for (name, amount) in present
    @printf("  %-18s %9.5f mol\n", name, amount)
end
```

Where each element ends up is the more useful reading, since the C-A-S-H holds
several of them at once:

```@example cem5
"""Moles of one conservation component held by each solid phase, largest first."""
function component_in_solids(component; tol = 1.0e-4)
    row = findfirst(==(component), String.(symbol.(cs.SM.primaries)))
    row === nothing && error("$component is not a component of this system")
    A = Float64.(cs.SM.A)
    out = [(symbol(cs.species[i]), A[row, i] * n[i]) for i in cs.idx_crystal]
    return sort(filter(p -> last(p) > tol, out); by = last, rev = true)
end

# The primary species carry one atom of their element each, so these rows count
# aluminum, magnesium and sulfur however the species that hold them are written.
for (component, element) in ("AlO2-" => "Al", "Mg+2" => "Mg", "SO4-2" => "S")
    println(element, ":")
    for (name, amount) in component_in_solids(component)
        @printf("  %-18s %9.5f mol %s\n", name, amount, element)
    end
    println()
end
```

## 5. The oxidation state

```@example cem5
r = half_reaction(eq, "SO4-2", "HS-")
println("half-reaction : ", r.equation)
@printf("log K at 25 C : %.2f\n", r.logK⁰(T = 298.15))
@printf("pe            : %+.2f\n", pe(eq, model))
@printf("Eh            : %+.3f V\n", Eh(eq, model))
@printf("aqueous S(VI)  : %.3e mol\naqueous S(-II) : %.3e mol\n",
        ustrip(us"mol", moles(eq, "SO4-2")), ustrip(us"mol", moles(eq, "HS-")))
```

The caution of [the CEM III page](@ref ex-cem3-slag) applies unchanged, and for
the same reason: a couple buffers a potential only while both members are in
solution, and a cement paste puts nearly all of its sulfur into solids. The
number is a bound on the reducing power the slag brings, not a measurement of the
paste's redox state — and even a buffered value would be the *mutual* equilibrium
of every couple, which sulfate reduction is far too slow to reach during
hydration.

## 6. Where a CEM V sits among the others

```@example cem5
function final_heat(file)
    t, Q = 0.0, 0.0
    for line in eachline(joinpath(datapath("experimental"), file))
        (startswith(line, "#") || startswith(line, "time")) && continue
        f = split(line, ',')
        t, Q = parse(Float64, f[1]), parse(Float64, f[3])
    end
    return t, Q
end

records = [
    ("smilauer2025-122-cemI-52.5R-cizkovice.csv", "CEM I 52.5 R", 0.0),
    ("smilauer2025-165-cemII-A-LL-42.5R-hranice.csv", "CEM II/A-LL", 0.13),
    ("smilauer2025-149-cemII-B-S-32.5R-mokra.csv", "CEM II/B-S", 0.28),
    ("smilauer2025-184-cemIII-A-42.5N-hranice.csv", "CEM III/A", 0.50),
    ("smilauer2025-200-cemV-A-S-V-32.5R-prachovice.csv", "CEM V/A (S-V)", 0.48),
    ("smilauer2025-121-cemIII-B-32.5N-mokra.csv", "CEM III/B", 0.73),
]
labels = String[]
heats = Float64[]
repl = Float64[]
for (file, label, r) in records
    _, Q = final_heat(file)
    push!(labels, label); push!(heats, Q); push!(repl, 100r)
    @printf("  %-16s replacement %4.0f %% (assumed)   Q = %6.1f J/g\n", label, 100r, Q)
end
```

```@example cem5
fig = scatter(repl, heats; legend = false, markersize = 7, color = :seagreen,
              xlabel = "clinker replaced (%, assumed at the EN 197-1 midpoint)",
              ylabel = "measured heat at ~28 days (J/g of binder)",
              title = "Every joule comes from the clinker",
              size = (780, 420),
              bottom_margin = 10Plots.mm, left_margin = 10Plots.mm)
for (x, y, l) in zip(repl, heats, labels)
    annotate!(fig, x, y + 8, text(l, 8, :left, :bottom))
end
plot!(fig; xlims = (-8, 85), ylims = (200, 410))
savefig(fig, "cem5-heat.svg"); nothing # hide
```

![](cem5-heat.svg)

!!! warning "The abscissa is assumed, the ordinate is measured"
    The heats are measured on one instrument at 20 °C. The replacement levels are
    **not in the deposit** — each is the midpoint of the EN 197-1 range for its
    designation, the same assumption every page here makes. So the figure shows a
    real trend read against an estimated axis, and the scatter about it is as much
    the spread of the ranges as it is chemistry. What it does establish is the
    ordering and the magnitude: half the clinker removed costs roughly a third of
    the heat.

The CEM V/A and the CEM III/A sit almost on top of each other, at comparable
replacement and comparable heat, which is what makes the family interesting: the
same reduction in heat and clinker, reached with two *different* constituents, and
therefore with a different assemblage, a different aluminum balance and a
different long-term reactivity. The equilibrium calculation above distinguishes
them where the calorimeter does not.
