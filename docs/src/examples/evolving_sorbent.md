# [A sorbent that appears and disappears](@id sec-example-evolving-sorbent)

Everywhere else in this documentation, a surface is something that is simply
*there*: you say how many moles of binding sites the system holds, and that
number stays what you said. This page is about the case where it cannot —
where the solid carrying the sites is itself dissolving, or has not
precipitated yet.

It is the ordinary case in a cement paste, where the C-S-H that binds alkalis
grows as the paste hydrates. It is shown here on a simpler material whose
numbers are published, so that every step can be checked.

## What a site is, and why its budget is a separate quantity

A mineral grain in water is not a clean cut through the crystal. The atoms at
the surface have bonds that stop, and water completes them: an iron at the
surface of a ferric oxide ends in a hydroxyl, written `≡FeOH`. That hydroxyl is
a **site** — one place where something can happen. It can take a proton and
become positive, give one up and become negative, or hand its place to a metal
ion arriving from solution.

What makes a site different from an element is that it is **conserved without
being matter**. Turning `≡FeOH` into `≡FeOH₂⁺` changes no count of iron, oxygen
or hydrogen on the grain as a whole — the proton came from the water — but the
*number of sites* is unchanged too, and nothing in the element balance says so.
So the package carries a conservation row of its own for it, exactly as it does
for electric charge. [Surface complexation](@ref sec-theory-surface) explains
why that row is the right way to write it.

## Why a fixed number cannot describe a dissolving sorbent

Suppose you post the site budget as a constant: `2 × 10⁻⁴ mol`, say. Now add
acid. The oxide dissolves; at some point half of it is gone. Half the sites went
with it — but the budget still says `2 × 10⁻⁴`, so the calculation keeps half a
surface that no longer exists, and whatever was bound to it stays bound.

You could of course work out how much oxide survives, then re-post the budget.
But how much survives *depends on* the surface, because the surface consumes
acid too. You would need the answer before you could set up the question.

The way out is to stop posting a number and state the **relation** instead: one
mole of this solid carries `ν` moles of sites, whatever amount of it there
turns out to be. That relation is linear, it fits in the conservation matrix
next to the element balances, and the solver then finds the amount of solid and
the amount of sites together.

## The material, and where its numbers come from

Hydrous ferric oxide — amorphous `Fe(OH)₃`, the rust-colored precipitate that
forms wherever dissolved iron meets oxygen — is the reference sorbent of the
field. [DzombakMorel1990](@cite) measured its site densities and its binding
constants, and those are the numbers PHREEQC ships to this day.

Two numbers define the surface:

  - **0.2 mol of weak sites per mole of Fe(OH)₃.** That is `ν`. (There is a
    second, much rarer family of "strong" sites at `0.005 mol/mol`; leaving it
    out changes nothing about what this page shows, and the two-family case is
    in [Two families of sites, and a metal between them](@ref sec-example-hfo).)
  - **`log K = 7.29` and `−8.93`** for taking up and giving up a proton.

and one more for the metal: `log K = −3.50` for
`≡FeOH + Mn²⁺ ⇌ ≡FeOMn⁺ + H⁺`. Manganese is chosen because it binds firmly
enough that the release is visible, and because the database used here carries
its aqueous chemistry.

```@example sorbent
using ChemistryLab, DynamicQuantities, SciMLBase, OptimaSolver, LinearAlgebra

const RT = R_GAS * 298.15
lnK(logK) = -RT * log(10.0^logK)         # ΔrG° = −RT ln K
energy(v) = SymbolicFunc(v * u"J/mol")

ν = 0.2                                  # mol of sites per mole of Fe(OH)₃
logK_prot, logK_depr, logK_Mn = 7.29, -8.93, -3.50
nothing # hide
```

## The aqueous side

Everything that is not the surface comes from a database the package ships —
PSI/Nagra here, because it is the one that carries amorphous ferric hydroxide as
a solid. See [Where the numbers come from](@ref sec-manual-numbers).

```@example sorbent
subs = build_species(datapath("psinagra-12-07-thermofun.json"); verbose = false)
db = Dict(symbol(s) => s for s in subs)
host = db["Fe(OH)3(am)"]

# Iron in one oxidation state only: there is no reductant here, so the ferrous
# species would be an unused half of a redox couple. Carbonate, sulfate and
# fluoride are left out for the same reason — nothing in the recipe supplies
# them.
aqueous = speciation(
    subs, ["Fe(OH)3(am)", "Cl-", "Mn+2"];
    aggregate_state = [AS_AQUEOUS],
    exclude_species = split(
        "H2@ O2@ Fe+2 FeOH+ FeO+ FeCl+ Cls ClO4- " *
            "MnSeO4@ MnSeO3w2(cr) Mn(CO3)@ Mn(HCO3)+ MnF+ Mn(SO4)@"
    ),
)
symbol.(aqueous)
```

## The surface species

A surface species is an ordinary species. What marks it is its aggregate state
and the **site symbol in its formula** — `Xw` here — which is what puts it in
the site balance.

```@example sorbent
surface(sym, g) = (
    s = Species(sym; aggregate_state = AS_SURFACE, class = SC_SURFCOMPLEX);
    s[:ΔₐG⁰] = energy(g); s
)
nothing # hide
```

### The one number a novice will get wrong

The free site needs a standard energy, and the obvious choice — zero — is wrong
here, in a way worth understanding because it is invisible everywhere else.

With a **fixed** site budget, zero is genuinely free. Every surface reaction has
a site on both sides, so whatever you put there cancels; the package measures
that invariance as part of its test suite.

With a budget that follows its host it does **not** cancel, and the reason is
simple once said: `≡FeOH` **is made of something**. It carries a real oxygen and
a real hydrogen. Saying its formation costs nothing says that a hydroxyl appears
from the elements for free — and now that the host is the thing supplying the
sites, that free energy is charged to the host, which changes how soluble the
host is.

So the reference is not a convention to pick. It is the energy of the matter the
free site carries, and the package reads it off the same matrix the constraint
is built from:

```@example sorbent
G(s) = ustrip(us"J/mol", s[:ΔₐG⁰](T = 298.15u"K", P = 1.0e5u"Pa"; unit = true))
reference = G(db["H2O@"]) - G(db["H+"])        # μ°(H₂O) − μ°(H⁺)
reference / 1000                                # kJ/mol
```

```@example sorbent
free = surface("XwOH", reference)
prot = surface("XwOH2+", reference + lnK(logK_prot))
depr = surface("XwO-", reference + lnK(logK_depr))
# The manganese complex carries the aqueous manganese's own energy, because
# ≡FeOH + Mn²⁺ = ≡FeOMn⁺ + H⁺ puts it in the balance.
mn = surface("XwOMn+", reference + lnK(logK_Mn) + G(db["Mn+2"]))
atoms(free)
```

## Declaring the coupling

Three pieces: the support says *what* carries the sites and that they follow it;
the capacity says *how many* per unit of it; the family gathers the species.

```@example sorbent
M = ustrip(us"kg/mol", host[:M])               # kg per mole of Fe(OH)₃

support = SurfaceSupport(
    "hydrous ferric oxide", "Fe(OH)3(am)", FixedSurfaceArea(1.0);
    coupling = SITES_FOLLOW_HOST,
)
family = SiteFamily(
    "Xw", free, [prot, depr, mn];
    capacity = MassSiteDensity(ν / M),         # mol of sites per kg of host
    support,
)
ChemistryLab.sites_per_host(family, M)         # back to ν, by the same arithmetic
```

Naming a host does not couple anything by itself — a rate law has named one
since long before, to know what amount it scales with. The coupling is asked
for, with `coupling = SITES_FOLLOW_HOST`.

The capacity has to be one that is **proportional** to how much host there is,
which is why it is given per unit mass. The package does not take that on
trust: it evaluates the capacity at two host amounts and checks that the budget
scales with them. A capacity posted as a total number of moles is refused,
because a total is not a relation.

## The system, and one declaration that is not obvious

```@example sorbent
bare = Species("Xw+"; aggregate_state = AS_SURFACE, class = SC_SURFCOMPLEX)

cs = ChemicalSystem(
    AbstractSpecies[vcat(aqueous, [free, prot, depr, mn])...],
    AbstractSpecies[db["H2O@"], db["H+"], db["Fe+3"], db["Cl-"], db["Mn+2"], bare];
    site_families = [family],
)
symbol.(cs.SM.primaries)
```

`Xw+` is in that list of components and **not** in the list of species, which
looks like a mistake and is not. A component is a bookkeeping direction, not a
substance: it is the bare site, with no matter attached, and the host supplies
`ν` of them per mole. It carries the charge the free site carries with its site
symbol — `≡FeOH` is a site plus an `OH⁻`, so the site is positive.
[Surface areas](@ref sec-manual-surfaces) works through both ways of getting
this wrong, and the package refuses both by name.

Now the check that catches the reference-energy mistake before it costs you an
afternoon:

```@example sorbent
host_coupling_bias(cs)     # log units of solubility, on the host itself
```

Zero, because the reference was set. Had it been left at zero, the same call
would read `8.3` — and a family worth more than `0.05` log units is refused at
construction, with the value to use in the message.

## The starting state

One kilogram of water, one millimole of the oxide, a tenth of a millimole of
manganese, and enough chloride to balance it.

```@example sorbent
state = ChemicalState(cs)
set_quantity!(state, "H2O@", 1.0u"kg")        # `set_quantity!` converts
set_quantity!(state, "Fe(OH)3(am)", 1.0e-3u"mol")
set_quantity!(state, "Mn+2", 1.0e-4u"mol")
set_quantity!(state, "Cl-", 2.0e-4u"mol")

# Derive the free-site amount from the declaration rather than typing it: the
# state and the capacity then say the same thing by construction.
state = host_consistent_state(state)
site_budget_residual(state)
```

That residual is the difference between the sites the state carries and the
sites the capacity declares. Zero means they agree. If you set the amounts by
hand and they disagree, [`check_site_budget`](@ref) turns it into an error
rather than a quietly different calculation.

## Dissolving the sorbent

The experiment is a titration: add hydrochloric acid, in steps counted as moles
per mole of `Fe(OH)₃`. Three moles is the stoichiometric amount to dissolve it
all, and rather more than that is needed in practice, because the surface and
the dissolved iron both take a share.

```@example sorbent
A = conservation_matrix(cs)                    # the coupled matrix, not SM.A
b0 = A * Float64[ustrip(us"mol", x) for x in state.n]
prim = string.(symbol.(cs.SM.primaries))
iH, iCl = findfirst(==("H+"), prim), findfirst(==("Cl-"), prim)
idx = Dict(s => i for (i, s) in enumerate(symbol.(cs.species)))
model = DaviesActivityModel()

acid = [0.0, 0.05, 0.1, 0.2, 0.35, 0.5, 0.75, 1.0, 2.0, 4.0, 6.0, 6.5]

function titrate(x)
    b = copy(b0)
    b[iH] += x * 1.0e-3                        # HCl: a proton and a chloride
    b[iCl] += x * 1.0e-3
    eq, cert = equilibrate_certified(state; model = model, b = b)
    n = Float64[ustrip(us"mol", v) for v in eq.n]
    sites = sum(n[idx[s]] for s in ("XwOH", "XwOH2+", "XwO-", "XwOMn+"))
    return (
        eq = eq, optimal = cert.optimal,
        oxide = n[idx["Fe(OH)3(am)"]], sites = sites,
        bound = n[idx["XwOMn+"]] / 1.0e-4, pH = pH(eq, model),
    )
end

runs = titrate.(acid)
all(r -> r.optimal, runs)                      # every point is a certified minimum
```

```@example sorbent
using Printf
println(" HCl/Fe   Fe(OH)₃ mol   sites mol    sites/Fe(OH)₃   Mn bound %    pH")
for (x, r) in zip(acid, runs)
    @printf(
        "%6.2f   %10.3e   %10.3e   %12.6f   %9.2f   %6.3f\n",
        x, r.oxide, r.sites, r.sites / r.oxide, 100r.bound, r.pH
    )
end
```

Two different things happen, in that order, and telling them apart is the point
of the page.

**Up to about 0.2 mol of acid per mole of oxide, nothing dissolves.** The oxide
column does not move and neither does the site column. What the acid does is
protonate the surface, and in doing so it takes the manganese back off it: bound
manganese falls from half of the inventory to essentially none while the sorbent
is still entirely present. A fixed site budget describes this perfectly well.

**Past that, the oxide goes, and the sites go with it.** The two columns fall
together, and the ratio in the middle column stays at `ν` to six decimals all
the way down — that is the relation being enforced, not fitted.

```@example sorbent
using Plots

p = plot(
    layout = (2, 1), size = (700, 620), link = :x,
    left_margin = 5Plots.mm, bottom_margin = 4Plots.mm,
)
plot!(
    p[1], acid, [r.oxide for r in runs] .* 1e3;
    label = "Fe(OH)₃ present", color = :firebrick, linewidth = 2, marker = :circle,
    ylabel = "mmol", title = "The sorbent, and the sites it carries",
)
plot!(
    p[1], acid, [r.sites for r in runs] .* 1e3 ./ ν;
    label = "sites ÷ ν  (should lie on top)", color = :black,
    linewidth = 2, linestyle = :dash,
)
plot!(
    p[2], acid, [100r.bound for r in runs];
    label = "manganese bound", color = :steelblue, linewidth = 2, marker = :circle,
    xlabel = "mol HCl per mol Fe(OH)₃", ylabel = "% of the manganese",
    legend = :topright,
)
plot!(twinx(p[2]), acid, [r.pH for r in runs];
      label = "pH", color = :seagreen, linewidth = 2, linestyle = :dot,
      ylabel = "pH", legend = :right)
p
```

The dashed line is the site total divided by `ν`. It lies on the solid one
because the constraint says it must; if it ever did not, the calculation would
not be the one that was declared.

## The three things worth checking afterwards

A solve that converges, an answer that is physically admissible and a certified
minimum are three different statements. Here they are separately.

```@example sorbent
worst_ratio = maximum(abs(r.sites / r.oxide / ν - 1) for r in runs)
```

**The relation holds.** The site total is `ν` times the host amount, everywhere,
to that relative accuracy.

```@example sorbent
element(n, e) = sum(
    n[i] * Float64(get(atoms(cs.species[i]), e, 0)) for i in eachindex(n)
)
n0 = Float64[ustrip(us"mol", x) for x in state.n]
worst_element = maximum(
    abs(element(Float64[ustrip(us"mol", v) for v in r.eq.n], e) - element(n0, e)) /
        element(n0, e) for e in (:Fe, :Mn), r in runs
)
```

**The elements are conserved**, computed from the species' own declared
formulas rather than from the matrix, so this checks the matrix too. It is the
test that would have caught the tempting shortcut of subtracting the coupling
from the free site's row, which quietly invents `ν` moles of oxygen and `ν` of
hydrogen per mole of host.

```@example sorbent
worst_si = maximum(
    abs(saturation_indices(r.eq, model)["Fe(OH)3(am)"])
        for r in runs if r.oxide > 1.0e-9
)
```

**Every phase that is present sits at `log SI = 0`.** That is not a coincidence
to be admired: it is the definition of being at equilibrium with the solution,
and an index that disagreed with it would be a diagnostic that lies.

## What PHREEQC does, and where this differs

PHREEQC has been able to tie a surface to a mineral since version 2, written

```
Hfo_wOH   Fe(OH)3(a)   equilibrium_phase   0.2
```

which reads exactly as `SITES_FOLLOW_HOST` does here. The generator
`test/reference/phreeqc_evolving_surface.py` titrates a sorbent to exhaustion
with it, and the two codes agree on the coupling law to `2 × 10⁻¹⁰` over five
partially dissolved states.

They are **not** running one surface model, and the difference is worth knowing
before comparing anything else. PHREEQC scales the site totals from the phase
and stops there: the phase's own stability is untouched. This package carries
the coupling in the conservation matrix, so the host's saturation index picks up
the site potential — which is why the free site's reference energy matters here
and not there. The two therefore agree on the ratio of sites to phase and on the
amounts the acid budget fixes, and are expected to differ on dissolved iron.

## What this page does not cover

  - **The sorbent exhausted.** Past about seven moles of acid per mole of oxide
    the solid is gone, its whole family sits at the solver's floor, and the
    solve stops certifying. The sweep above stops before that on purpose; a
    surface with no host is a system to describe differently, not a limit to
    take.
  - **A sorbent that *appears*.** Everything here runs backwards as well — the
    same declaration describes a phase precipitating and bringing its sites with
    it — but the case that matters, a C-S-H forming as a paste hydrates, needs
    site densities for the C-S-H that nobody has published. That is a missing
    measurement, not a missing feature.
  - **Electrostatics.** The surface here has a charge and no potential; see
    [A charged surface, screened](@ref sec-example-diffuse-layer).
  - **The reference energy, taken seriously.** Setting it to the matter the free
    site carries is what this package does, and it is stated rather than
    inherited from a database. [Kulik2002](@cite) avoids needing it at all, by
    keeping the free site out of the balance as a *surface monolayer solvent* of
    fixed activity. That is a different formulation, and the one to read before
    pushing this to site densities where the approximation would show.

## See also

  - [Surface areas](@ref sec-manual-surfaces) — the syntax, and the two ways of
    declaring the component that are refused.
  - [Surface complexation](@ref sec-theory-surface) — why a site balance is a
    conservation row, and how Langmuir falls out of it.
  - [Two families of sites, and a metal between them](@ref sec-example-hfo) —
    the same oxide with a fixed budget, against PHREEQC.
