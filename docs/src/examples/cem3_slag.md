# [A blastfurnace cement, and the oxidation state it needs](@id ex-cem3-slag)

A CEM III is between a third and four fifths blastfurnace slag. That single fact
changes the calculation in three ways, and this page follows each:

1. the slag is a **glass** — it has no phases to name, so it enters as an oxide
   analysis rather than through a Bogue conversion;
2. it brings **magnesium**, which the clinker barely has, and magnesium makes
   hydrotalcite;
3. it brings sulfur as **S(-II)** into a pore solution whose sulfur is
   **S(+VI)**, and a calculation that cannot hold both has decided the answer
   before it starts. That is what the [oxidation state](@ref theory-redox)
   machinery is for.

```@example cem3
using ChemistryLab
using DynamicQuantities
using OptimaSolver
using OrderedCollections
using Printf
using Plots
default(framestyle = :box, grid = false)

substances = build_species(datapath("cemdata18-thermofun.json"); verbose = false)
byname = Dict(symbol(s) => s for s in substances)
nothing # hide
```

## 1. What is measured, and what is assumed

This page is built around a **measured** record: CEM III/A 42.5 N from Hranice,
in the CC-BY-4.0 deposit of [Smilauer2025data](@cite). What that record gives,
and what it does not, decides how the rest must be written.

```@example cem3
rec = joinpath(datapath("experimental"), "smilauer2025-184-cemIII-A-42.5N-hranice.csv")
for line in eachline(rec)
    startswith(line, "#") || break
    occursin(r"cement:|blaine:|wb:|temperature:", line) && println(line)
end
```

Three numbers, all measured. **The clinker phase composition and the actual slag
content are not in the deposit**, and nothing in it lets them be recovered. They
have to be assumed, and the assumption is stated here rather than buried:

```@example cem3
# ASSUMED, not measured: the midpoint of the EN 197-1 range for CEM III/A,
# which is 35-64 % clinker and 36-65 % slag.
CLINKER_FRACTION = 0.50
SLAG_FRACTION = 0.50

# ASSUMED: a Bogue composition representative of a CEM I clinker. The deposit
# does not report one for this cement.
CLINKER = OrderedDict("C3S" => 0.65, "C2S" => 0.11, "C3A" => 0.11, "C4AF" => 0.08)

# ASSUMED: a European ground granulated blastfurnace slag analysis. The sulfur
# is the part that matters here, and it is the part a datasheet reports least
# consistently.
SLAG = Dict("CaO" => 0.41, "SiO2" => 0.36, "Al2O3" => 0.11,
            "MgO" => 0.08, "SO3" => 0.02)

WB = 0.40            # MEASURED, from the record above
BINDER_G = 100.0
nothing # hide
```

!!! warning "Read the results as a family, not as this specimen"
    Everything below follows from those three assumptions as much as from the
    thermodynamics. A different slag content inside the same EN 197-1 range
    moves the assemblage; a different clinker Bogue moves it again. What the
    calculation shows is how a **CEM III/A behaves**, not what this particular
    cement from Hranice contains.

## 2. The slag enters as an element budget

There is no `Species` for a slag, because a glass has no formula.
[`oxide_budget`](@ref) converts the analysis to the component totals the
constraint `A n = b` is written over, computing every molar mass from its
formula:

```@example cem3
pure = split(
    "C3S C2S C3A C4AF Gp Anh Cal Portlandite ettringite monosulphate12 " *
    "monocarbonate hemicarbonate C4AH13 C3AH6 C3FH6 straetlingite " *
    "hydrotalcite Mg2AlC0.5OH Brc FeOOHmic AlOHmic Amor-Sl Tro Py Sulfur Mgs"
)
solutions = [
    "CSHQ" => ["CSHQ-JenD", "CSHQ-JenH", "CSHQ-TobD", "CSHQ-TobH", "KSiOH", "NaSiOH"],
    "C3(AF)S0.84H" => ["C3AFS0.84H4.32", "C3FS0.84H4.32"],
]
members = reduce(vcat, last.(solutions))

# The sulfur ladder must be in the species list, or "holding both oxidation
# states" is a claim about species that are not there.
redox_species = ["HS-", "H2S@", "SO4-2", "SO3-2", "S2O3-2", "O2@", "H2@"]

species = speciation(substances, vcat(pure, members, redox_species);
                     aggregate_state = [AS_AQUEOUS])
ss = [SolidSolutionPhase(n, [byname[m] for m in ms]) for (n, ms) in solutions]
cs = ChemicalSystem(species, CEMDATA_PRIMARIES; solid_solutions = ss)

@printf("%d species, %d conservation components\n",
        length(cs.species), size(cs.SM.A, 1))
```

The component list is where the redox treatment announces itself. With sulfur at
two oxidation states in the species list, charge is no longer a fixed
combination of the element rows, so it survives as a component of its own:

```@example cem3
components = String.(symbol.(cs.SM.primaries))
@printf("charge (`Zz`) kept as a conservation component: %s\n",
        "Zz" in components ? "yes — the oxidation state is conserved" : "no")
```

## 3. The budget, clinker plus glass

```@example cem3
state = ChemicalState(cs)
molar_mass(n) = ustrip(us"g/mol", byname[n][:M])

for (phase, frac) in CLINKER
    set_quantity!(state, phase,
        BINDER_G * CLINKER_FRACTION * frac / molar_mass(phase) * u"mol")
end
set_quantity!(state, "H2O@", BINDER_G * WB / molar_mass("H2O@") * u"mol")

b_clinker = Float64.(cs.SM.A) * ustrip.(us"mol", state.n)
b_slag = oxide_budget(SLAG, cs.SM.primaries; mass = BINDER_G * SLAG_FRACTION * u"g")
b = b_clinker .+ b_slag

for (comp, v) in zip(components, b)
    abs(v) > 1.0e-6 && @printf("  %-8s %10.5f mol\n", comp, v)
end
```

## 4. The equilibrium

```@example cem3
model = HKFActivityModel(å = 0.0, Ḃ = 0.097637, Kₙ = 0.0)
eq, cert = equilibrate_certified(state; model = model, b = b)

@printf("certificate: optimal=%s  worst SI=%.2e  element balance=%.1e\n",
        cert.optimal, cert.worst_supersaturation, cert.balance)
@printf("pH = %.3f\n", pH(eq, model))
```

The assemblage, in moles per 100 g of binder:

```@example cem3
n = ustrip.(us"mol", eq.n)
present = sort([(symbol(cs.species[i]), n[i]) for i in cs.idx_crystal if n[i] > 1.0e-4];
               by = last, rev = true)
for (name, amount) in present
    @printf("  %-18s %9.5f mol\n", name, amount)
end
```

## 5. The oxidation state, which is the point

A CEM I has no answer to give here: its sulfur is all sulfate, so every couple
is degenerate. A CEM III does.

```@example cem3
r = half_reaction(eq, "SO4-2", "HS-")
println("half-reaction : ", r.equation)
@printf("log K at 25 C : %.2f\n", r.logK⁰(T = 298.15))
@printf("pe            : %+.2f\n", pe(eq, model))
@printf("Eh            : %+.3f V\n", Eh(eq, model))

s6 = ustrip(us"mol", moles(eq, "SO4-2"))
s2 = ustrip(us"mol", moles(eq, "HS-"))
@printf("\naqueous S(VI)  : %.3e mol\naqueous S(-II) : %.3e mol\n", s6, s2)
```

!!! danger "This is an equilibrium answer, and sulfate reduction is slow"
    The number above is what the thermodynamics gives if every redox couple is
    allowed to reach mutual equilibrium. In a real paste that is false on the
    time scale of hydration: sulfate reduction is kinetically frozen, so the
    sulfur stays closer to what the slag brought than to what equilibrium would
    make of it.

    A real slag paste is therefore **somewhere between** its initial oxidation
    state and this one, and nothing in a Gibbs minimization says where. The
    package can now pose the question; answering it needs a kinetic description
    of sulfate reduction, which it does not have.

## 6. What the measured calorimetry says, and what it does not

The record gives the heat, and the heat is the one quantity here that was
measured rather than assumed. It is worth putting beside the CEM I of the same
deposit:

```@example cem3
function final_heat(file)
    t, Q = 0.0, 0.0
    for line in eachline(joinpath(datapath("experimental"), file))
        startswith(line, "#") && continue
        startswith(line, "time") && continue
        f = split(line, ',')
        t, Q = parse(Float64, f[1]), parse(Float64, f[3])
    end
    return t, Q
end

for (file, label) in (
        ("smilauer2025-122-cemI-52.5R-cizkovice.csv", "CEM I 52.5 R"),
        ("smilauer2025-184-cemIII-A-42.5N-hranice.csv", "CEM III/A 42.5 N"),
        ("smilauer2025-121-cemIII-B-32.5N-mokra.csv", "CEM III/B 32.5 N"),
    )
    t, Q = final_heat(file)
    @printf("  %-18s %6.0f h   %6.1f J/g\n", label, t, Q)
end
```

Replacing clinker with slag lowers the heat, because **every joule comes from
the clinker** and slag reacts later and less exothermically. That is the
property a CEM III is specified for — a massive pour that would crack under the
thermal gradient of a CEM I.

What the calculation above does **not** do is predict that curve. It is an
equilibrium: it says what the paste tends to, not how fast. Coupling it to the
Waller kinetics that ship for slag is the subject of
[the coupled runs](@ref ex-ionic-opc), and doing it on a blend needs a degree of
reaction for the slag that this deposit does not report either.
