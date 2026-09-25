# [A blastfurnace cement, and the oxidation state it needs](@id ex-cem3-slag)

!!! info "Before this page"
    [Two CEM II](@ref ex-cem2-blended), whose skeleton this page follows, and
    [Oxidation state, and the potential conjugate to it](@ref theory-redox).

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

## 0. Why a glass is a different kind of input

A clinker phase has a formula. ``\mathrm{C_3S}`` is ``\mathrm{Ca_3SiO_5}``, and
the package can look up its molar mass, its Gibbs energy and its molar volume.
A **glass cannot be written that way at all.**

Ground granulated blastfurnace slag is quenched molten slag: its atoms are frozen
in a disordered network with no repeating unit, so there is no formula unit to
name and no tabulated thermodynamic data for "slag". What a datasheet reports is
an **oxide analysis** — how much CaO, SiO₂, Al₂O₃, MgO the material contains by
mass — and nothing else.

That is enough, because a Gibbs minimization does not need to know what the
starting material *was*. It needs the totals:

```math
\mathbf{A}\,\mathbf{n} = \mathbf{b}
```

where ``\mathbf{n}`` holds the amount of each species, ``\mathbf{A}`` says how
many units of each conserved component each species carries, and ``\mathbf{b}``
is the budget — how much of each component the paste contains in total. The
minimization finds the ``\mathbf{n}`` of lowest Gibbs energy subject to that.

A clinker enters through its phases, which have formulas, and ``\mathbf{b}``
follows as ``\mathbf{A}\mathbf{n}_0``. A glass enters by going straight to
``\mathbf{b}``: [`oxide_budget`](@ref) converts an oxide analysis into component
totals, and the two contributions are simply added.

!!! warning "A budget says what a material CONTAINS, not what it does"
    A slag and a quartz sand of the same analysis give the same ``\mathbf{b}``.
    Equilibrium has no notion of reactivity: it answers what the paste would
    become if everything reacted. How much of the glass actually dissolves, and
    how fast, is a **kinetic** question answered elsewhere — see
    [the coupled runs](@ref ex-ionic-opc).

    So read this page as the limit the paste tends to, not as a 28-day specimen.

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

# THE CLINKER'S ALKALIS, which Bogue does not account for and which set the pH.
#
# A Bogue calculation returns four phases and no sodium or potassium: they are
# minor oxides, a fraction of a percent, and they sit outside the four-phase
# decomposition. They are also, in a cement paste, **what fixes the pH** -- they
# dissolve almost completely into the pore solution and stay there, where the
# calcium is held down by portlandite at 12.5. Leaving them out does not make the
# calculation conservative: it makes it report a portlandite floor as though it
# were a pore solution.
#
# ASSUMED at a usual industrial level, as a fraction of the CLINKER mass.
ALKALIS = Dict("K2O" => 0.008, "Na2O" => 0.002)

WB = 0.40            # MEASURED, from the record above
BINDER_G = 100.0

# HOW MUCH REACTS. Two ceilings, and the reacted fraction is the lower.
#
# The WATER ceiling is Powers (1948): about 0.42 g of water per gram of cement
# is needed for complete hydration -- 0.23 g written into the hydrate formulae
# and 0.19 g held in the gel pores those hydrates create. Below that the paste
# stops WITH WATER STILL IN IT, what remains being in pores far too fine to
# reach an unhydrated grain. That argument is about the pore space and not about
# the grain, so it caps the slag exactly as it caps the alite. Water curing
# moves it to 0.36, the volume emptied by chemical shrinkage being refilled from
# outside -- `powers_alpha_max(WB; curing = :saturated)`.
ALPHA_WATER = powers_alpha_max(WB)                 # sealed, ASSUMED

# The KINETIC ceiling is the slag's own dissolution rate, and at 28 days it is
# the lower of the two by a wide margin. The RILEM TC 238-SCM round robin
# [Durdzinski2017](@cite) measured two ground granulated slags at 40 %
# replacement and w/b 0.40 in seven laboratories -- this page's geometry. Its
# Table 4 at 28 days, by SEM image analysis, the technique the study found most
# consistent: 38 % and 48 % for the first slag, 45 % and 49 % for the second.
# The study's verdict on the precision of any technique: "at best +/- 5 %".
#
# ASSUMED from that table. [The CEM V page](@ref cem5-dor) sweeps the same
# quantity across the round robin's 7-, 28- and 90-day columns.
ALPHA_SLAG = min(0.45, ALPHA_WATER)
ALPHA_CLINKER = ALPHA_WATER

@printf("water ceiling at w/b = %.2f : %.3f\n", WB, ALPHA_WATER)
@printf("reacted: clinker %.0f %%, slag %.0f %%\n", 100ALPHA_CLINKER, 100ALPHA_SLAG)
```

!!! warning "Read the results as a family, not as this specimen"
    Everything below follows from those assumptions as much as from the
    thermodynamics. A different slag content inside the same EN 197-1 range
    moves the assemblage; a different clinker Bogue moves it again; and the
    reacted fractions move the amounts more than either. What the calculation
    shows is how a **CEM III/A behaves**, not what this particular cement from
    Hranice contains.

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
molar_mass(n) = ustrip(us"g/mol", byname[n][:M])

"""
    paste(; alkali = 1.0) -> (; state, clinker, slag, total)

The fresh state and the element budget, with `alkali` scaling the clinker's
sodium and potassium so that section 6 can sweep the one input that sets the pH.

The state carries only what has reacted. The rest -- unhydrated clinker cores and
undissolved glass -- is still in the specimen, and is no part of the
minimization: an intact grain is not at equilibrium with the solution around it.
All of the mixing water enters, though: the ceiling limits how far the reaction
can go, not how much water was poured in, and the water that cannot reach a grain
is still in the balance.
"""
function paste(; alkali = 1.0)
    st = ChemicalState(cs)
    for (phase, frac) in CLINKER
        set_quantity!(st, phase,
            ALPHA_CLINKER * BINDER_G * CLINKER_FRACTION * frac / molar_mass(phase) * u"mol")
    end
    set_quantity!(st, "H2O@", BINDER_G * WB / molar_mass("H2O@") * u"mol")

    clinker = Float64.(cs.SM.A) * ustrip.(us"mol", st.n)
    # The alkalis follow the clinker, and its reacted fraction: they leave the
    # grain as it dissolves.
    clinker .+= alkali * oxide_budget(ALKALIS, cs.SM.primaries;
                                      mass = BINDER_G * CLINKER_FRACTION * ALPHA_CLINKER * u"g")
    slag = oxide_budget(SLAG, cs.SM.primaries;
                        mass = BINDER_G * SLAG_FRACTION * ALPHA_SLAG * u"g")
    return (; state = st, clinker, slag, total = clinker .+ slag)
end

p0 = paste()
state, b = p0.state, p0.total

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

Read that output carefully, because it contains a trap the calculation itself
warns about. **The aqueous sulfide is at the solver's floor** — the equilibrium
put essentially all of this paste's sulfur into the AFm phase
(`monosulphate12` above), leaving nothing in solution on either side of the
couple.

A couple buffers a potential only while both of its members are present. With
one at the floor, the `pe` printed above is set by the floor `ϵ` and not by the
chemistry, which is why [`pe`](@ref) emits a warning here rather than returning
the number silently. Read it as a **bound**, not as the redox state of the
paste.

That is not a defect of the calculation; it is the answer. A CEM III/A at this
sulfur content has no aqueous sulfur to speak of, so it has no sulfur redox
buffer either. A cement whose slag brings more sulfur, or a paste carbonated
enough to release the AFm sulfate, would.

!!! danger "And even a buffered answer would be an equilibrium answer"
    Suppose the couple were buffered. The potential would still be what the
    thermodynamics gives if every redox couple reaches **mutual** equilibrium,
    and on the time scale of hydration that is false: sulfate reduction is
    kinetically frozen, so the sulfur stays closer to what the slag brought than
    to what equilibrium would make of it.

    A real slag paste is therefore somewhere between its initial oxidation state
    and this one, and nothing in a Gibbs minimization says where. The package can
    now pose the question; answering it needs a kinetic description of sulfate
    reduction, which it does not have.

## [6. The one input that sets the pH](@id cem3-alkali)

Everything on this page is assumed except the water/binder ratio, the fineness
and the calorimetry — and of the assumptions, **one dominates the pore
solution**. The calcium is held at the portlandite floor and cannot rise; what
rises above it is the alkalis, which dissolve almost entirely and stay there. So
a page that assumes an alkali content owes the reader its sensitivity, and this
is it: the same paste with the clinker's Na₂O and K₂O scaled over the industrial
range, from a low-alkali cement to a high-alkali one.

```@example cem3
i_ch = findfirst(sp -> symbol(sp) == "Portlandite", cs.species)
na2o_eq(scale) = 100 * scale * (ALKALIS["Na2O"] + 0.658 * ALKALIS["K2O"])

@printf("%-14s %10s %11s %8s %13s\n",
        "Na2O eq (%)", "certified", "balance", "pH", "portlandite")

# `let` rather than a bare loop: a top-level `for` that assigns to a name of the
# enclosing scope makes a NEW LOCAL, so `prev` would be read before it is ever
# written. Wrapping the sweep gives it a scope of its own, which is cleaner than
# reaching for `global`.
let prev = nothing
    for scale in (0.5, 1.0, 1.5)
        pa = paste(; alkali = scale)
        e, c = equilibrate_certified(something(prev, pa.state);
                                     model = model, b = pa.total)
        c.optimal && (prev = e)
        nn = ustrip.(us"mol", e.n)
        @printf("%-14.2f %10s %11.1e %8.3f %13.5f\n",
                na2o_eq(scale), c.optimal, c.balance, pH(e, model), nn[i_ch])
    end
end
```

Read the two right-hand columns against each other. A **factor of three** on the
alkali content moves the pH by **0.40 unit** and the portlandite by under **5 %**.
That separation is the whole point: the calcium is held by portlandite and cannot
follow, so in a cement paste the alkalis *are* the pH and the calcium hydroxide
is only a floor beneath them.

The middle row is the page's own case, and it returns the 13.041 of section 4 —
which is worth checking rather than assuming, since a sweep that did not
reproduce its own nominal point would be measuring something else.

Two consequences for anyone using this page. A pH quoted from it is worth exactly
what the assumed alkali content is worth, so **substitute your own analysis** —
the constant is `ALKALIS` in section 1 and nothing else needs to change. And a
durability argument that turns on pore-solution pH — alkali-silica reaction,
steel passivation, leaching — cannot be settled by a calculation whose alkali
input was assumed. That is not a limitation of the minimization; it is what the
minimization is telling you about which measurement to go and make.

## 7. What the measured calorimetry says, and what it does not

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

## Where to go next

The pozzolanic binder, in which the C-S-H has to carry the aluminum, is
[A pozzolanic binder](@ref ex-cem4-pozzolanic).
