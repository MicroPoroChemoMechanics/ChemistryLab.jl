# [Effect of Water/Cement Ratio on Cement Hydration](@id sec-wc-ratio)

!!! info "Before this page"
    [Simplified Clinker Dissolution](@ref sec-clinker-dissolution).

The **water-to-cement ratio** (w/c) is the single most important mix-design parameter
of concrete. It controls workability, compressive strength, and durability simultaneously.
From a thermodynamic perspective, it determines how much water is available to hydrate the
clinker phases, which in turn governs the nature and amount of hydration products, the
porosity of the paste, and the pH of the pore solution.

This example scans w/c from 0.30 to 0.60 and tracks, at **full thermodynamic
equilibrium**, the pH of the pore solution, the hydrate assemblage, and the
porosity of the hardened paste in the sealed-curing convention.

!!! warning "Equilibrium answers a narrower question than mix design does"
    Read the scan below for what it is: the assemblage a paste would reach if
    every reaction ran to completion. It is **not** the state of a real paste at
    an age, and two of the effects usually attributed to w/c are absent from it
    by construction — over the range scanned no clinker survives, and there is no
    optimum. Those are **kinetic** limitations, carried by
    [`powers_alpha_max`](@ref) and the rate laws of the
    [kinetics tutorial](@ref sec-kinetics), not by the Gibbs minimum. A
    water-limited regime does exist in the Gibbs minimum, but it starts far below
    Powers' 0.42 — measured on this species list, between w/c = 0.28 and 0.30 —
    and the two limits are not the same statement, which the analysis below
    takes apart.
    This page is the reference the kinetic calculation converges toward; the
    [coupled hydration example](@ref sec-coupled-hydration) is the one to read
    for an age.

---

## System setup

The same clinker composition and species set as the "simplified clinker dissolution" example are used here.

```@example wc_setup
using ChemistryLab
using DynamicQuantities

substances = build_species(datapath("cemdata18-thermofun.json"))

input_species = split("C3S C2S C3A C4AF Gp Anh Portlandite Jennite H2O@ ettringite monosulphate12 C3AH6 C3FH6 C4FH13")
species = speciation(substances, input_species; aggregate_state = [AS_AQUEOUS])

cs = ChemicalSystem(species, CEMDATA_PRIMARIES)
nothing # hide
```

```@raw html
<details><summary>The chemical system in full — species, phases and the conservation matrix</summary>
```

```@example wc_setup
cs
```

```@raw html
</details>
```

The clinker composition (mass fractions of the anhydrous cement phases) is fixed throughout the scan. It is an assumed composition, typical of a CEM I; its source is not recorded.

| Phase | Symbol | Mass fraction |
|:------|:-------|:--------------|
| Alite | `C3S`  | 67.8 % |
| Belite | `C2S` | 16.6 % |
| Aluminate | `C3A` | 4.0 % |
| Ferrite | `C4AF` | 7.2 % |
| Gypsum | `Gp`  | 2.8 % |

```@example wc_setup
# An assumed CEM I composition, typical of the class; its source is not recorded.
compo = ["C3S" => 0.678, "C2S" => 0.166, "C3A" => 0.040, "C4AF" => 0.072, "Gp" => 0.028]
c     = sum(last.(compo))   # cement mass fraction (= 0.984 here)
```

---

## The solver

[`equilibrate_certified`](@ref) is used throughout, which needs nothing but
`OptimaSolver` loaded:

```@example wc_setup
using OptimaSolver
```

That route returns `(state, certificate)`, and the certificate is what makes the
numbers below quotable: for a convex problem the KKT conditions are *sufficient*,
so `cert.optimal == true` is a **proof** that the composition is the Gibbs
minimum and not merely the point an iteration stopped at. Every value on this
page is checked that way, and the check is not idle — on a cement the difference
between "the solver returned" and "the answer is proved" has been measured at
13 % of the total volume and four units of pH.

!!! note "Solving through Ipopt instead"
    An `EquilibriumSolver` around any nonlinear back end still works, and reaches
    the same assemblage here:

    ```julia
    using Optimization, OptimizationIpopt

    opt = IpoptOptimizer(
        acceptable_tol        = 1e-10,
        dual_inf_tol          = 1e-10,
        acceptable_iter       = 100,
        constr_viol_tol       = 1e-10,
        warm_start_init_point = "no",
    )
    solver = EquilibriumSolver(
        cs, DiluteSolutionModel(), opt;
        variable_space = Val(:linear), reltol = 1e-8,
    )
    eq = solve(solver, deepcopy(fresh))     # no certificate
    ```

    It is left out of the executed page for two reasons, neither of them about
    the quality of Ipopt as an optimizer. It is a **bare interior point**: it
    returns an iterate and no statement about it, so a caller has to audit the
    answer with [`optimality_certificate`](@ref) anyway. And it is an extra
    binary dependency for a calculation the package can already prove. The
    visible difference on this page is small but telling: the interior point
    leaves the absent phases at its lower bound, around `1e-8` mol, while the
    certified route puts them at exactly zero.

---

## Scanning the w/c ratio

For each value of w/c the fresh state is rebuilt from scratch, and the total mass is
normalized to 1 kg of paste (cement + water) so that all amounts are comparable
across the scan. The fresh state is kept: it is the volume reference the porosity
is referred to.

```@example wc_setup
sp_idx   = Dict(symbol(s) => i for (i, s) in enumerate(cs.species))
wc_range = range(0.30, 0.60; length = 13)

function fresh_paste(wc)
    w, mtot = wc * c, c + wc * c
    st = ChemicalState(cs)
    for (sym, mfrac) in compo
        set_quantity!(st, sym, mfrac / mtot * u"kg")
    end
    set_quantity!(st, "H2O@", w / mtot * u"kg")
    V = volume(st)
    set_quantity!(st, "H+",  1e-7u"mol/L" * V.liquid)   # charge seed, pH-neutral
    set_quantity!(st, "OH-", 1e-7u"mol/L" * V.liquid)
    return st
end

pH_vals    = Float64[]
ϕ_liquid   = Float64[]
ϕ_void     = Float64[]
ϕ_total    = Float64[]
n_portl    = Float64[]
n_mono     = Float64[]
n_ett      = Float64[]
n_jennite  = Float64[]
n_clinker  = Float64[]

certified = Bool[]

for wc in wc_range
    fresh    = fresh_paste(wc)
    eq, cert = equilibrate_certified(deepcopy(fresh))   # the solve may mutate its argument
    push!(certified, cert.optimal)
    ϕ = porosity(eq, fresh)                  # sealed-curing convention
    amount(sym) = ustrip(eq.n[sp_idx[sym]])
    push!(pH_vals,   pH(eq))
    push!(ϕ_liquid,  ϕ.liquid)
    push!(ϕ_void,    ϕ.void)
    push!(ϕ_total,   ϕ.total)
    push!(n_portl,   amount("Portlandite"))
    push!(n_mono,    amount("monosulphate12"))
    push!(n_ett,     amount("ettringite"))
    push!(n_jennite, amount("Jennite"))
    push!(n_clinker, sum(amount(s) for s in ("C3S", "C2S", "C3A", "C4AF")))
end
nothing # hide
```

!!! note "Why `porosity(eq, fresh)` and not `porosity(eq)`"
    The one-argument [`porosity`](@ref) divides the pore volume by the **current**
    total volume, which is right for a fixed-volume aqueous system and wrong for a
    setting binder: hydration products occupy less space than the reactants they
    consume, so the paste's own volume shrinks under the ratio. The two-argument
    form refers everything to the fresh volume and splits the result into the
    water-filled porosity and the empty porosity that the Le Chatelier
    contraction creates. On this mix at w/c = 0.50 the one-argument form returns
    0.285 against a total of 0.338, a gap of 5.3 points of porosity, the total
    volume having shrunk by 7.4 %.

---

## Results

### pH and porosity

```@example wc_setup
using Plots

p1 = plot(collect(wc_range), pH_vals;
    xlabel = "w/c ratio", ylabel = "Pore solution pH", label = "pH",
    linewidth = 2, marker = :circle, markersize = 4, color = :steelblue,
    title = "Pore solution pH", ylims = (11.5, 13.5), legend = :bottomright)

p2 = plot(collect(wc_range), ϕ_total .* 100;
    xlabel = "w/c ratio", ylabel = "Porosity (%)", label = "total",
    linewidth = 2, marker = :circle, markersize = 4, color = :firebrick,
    title = "Porosity referred to the fresh volume",
    ylims = (0, 50), legend = :topleft)
plot!(p2, collect(wc_range), ϕ_liquid .* 100;
    label = "water-filled", linewidth = 2, marker = :square, markersize = 3,
    color = :steelblue)
plot!(p2, collect(wc_range), ϕ_void .* 100;
    label = "empty (Le Chatelier)", linewidth = 2, marker = :diamond,
    markersize = 3, color = :seagreen)

plot(p1, p2; layout = (1, 2), left_margin = 8Plots.mm,
     bottom_margin = 8Plots.mm, size = (950, 410))
```

The pH is **flat to three decimals** over the whole range. That is the signature
of a buffered solution: portlandite is present at every w/c, so the calcium and
hydroxide activities are pinned by its saturation, and in a dilute-solution model
those activities do not know how large the pore volume is. Diluting a saturated
solution with more of its own solvent does not change its pH — it dissolves more
portlandite. The pH would start to move only once portlandite is exhausted, which
this composition never does.

### Phase assemblage

```@example wc_setup
p3 = plot(collect(wc_range), n_portl;
    xlabel = "w/c ratio", ylabel = "Amount (mol / kg of paste)",
    label = "Portlandite  Ca(OH)₂", linewidth = 2, marker = :circle,
    markersize = 4, color = :steelblue, title = "Hydrate assemblage at equilibrium",
    legend = :right)
plot!(p3, collect(wc_range), n_jennite;
    label = "Jennite (C-S-H)", linewidth = 2, marker = :square,
    markersize = 3, color = :firebrick)
plot!(p3, collect(wc_range), n_mono;
    label = "monosulphate12 (AFm)", linewidth = 2, marker = :diamond,
    markersize = 3, color = :seagreen)
plot!(p3, collect(wc_range), n_ett;
    label = "ettringite (AFt)", linewidth = 2, marker = :utriangle,
    markersize = 3, color = :darkorange)
plot(p3; left_margin = 8Plots.mm, bottom_margin = 8Plots.mm, size = (700, 420))
```

```@example wc_setup
using Printf
@printf "every point certified                    : %s\n" all(certified)
@printf "ettringite, largest value over the scan  : %.3e mol\n" maximum(n_ett)
@printf "clinker left, largest value over the scan: %.3e mol\n" maximum(n_clinker)
```

---

## Analysis

| Quantity | What the scan gives | Why |
|:--|:--|:--|
| **pH** | 12.39, constant to three decimals | buffered by portlandite saturation; independent of pore volume in a dilute model |
| **Portlandite, C-S-H, AFm** | decrease with w/c, in proportion to the cement fraction | amounts are per kg of *paste*; adding water dilutes the binder, it does not change what a gram of cement produces |
| **Ettringite** | exactly zero at every w/c | see below |
| **Clinker left** | exactly zero at every w/c of *this* scan | see below, and it is not "at any w/c" |
| **Total porosity** | 12.2 % to 41.0 %, monotone, no optimum | the excess water has nowhere to go but the pore space |
| **Empty porosity** | 9.8 % down to 6.6 % | the Le Chatelier contraction is roughly fixed per gram of cement, so it is a smaller fraction of a larger reference volume |

Two of those deserve stating plainly, because the intuition of mix design points
the other way.

!!! warning "Ettringite does not form here, and that is correct"
    AFt needs about three sulfates per aluminate; this clinker has 2.8 % gypsum
    against 4.0 % C₃A plus the ferrite, so at **equilibrium** the sulfate is all
    taken up by the AFm phase `monosulphate12`, and `ettringite` comes back at
    exactly zero — absence, not a small amount. Ettringite is the phase that forms **early**, while
    sulfate is still locally abundant, and then converts to AFm as it runs out.
    It is a kinetic intermediate for this mix, so an equilibrium scan cannot show
    it. Raise the gypsum content and AFt becomes stable — that is the
    sulfate-balance calculation, and it is worth doing before reading anything
    into an AFt amount.

!!! warning "No clinker survives *over this range*, and Powers' 0.42 is not the reason"
    The four anhydrous phases come back at exactly zero at every w/c of the scan,
    0.30 included, and the certificate proves it. That looks like it contradicts
    [Powers1948](@cite), whose ``\alpha_{\max} = w/c \,/\, 0.42`` leaves
    unreacted clinker in any paste below w/c = 0.42 — 0.36 with curing water. It
    does not, and the reason is worth stating because the two numbers measure
    different things.

    Powers' 0.42 g of water per gram of cement is **not** a stoichiometric
    demand. It is about 0.23 g of *non-evaporable* water, which is the water
    written into the hydrate formulae, plus about 0.19 g of **gel water** held in
    the C-S-H gel pores. Only the first is a mass balance that a Gibbs
    minimization must respect. The second is water that is physically there and
    chemically unavailable: in a sealed paste it is immobilized in pores too fine
    to feed further reaction, and hydration stops by self-desiccation with water
    still in the specimen. That is a statement about **transport and access**,
    which no equilibrium calculation contains, and it is what
    [`powers_alpha_max`](@ref) carries into the kinetic rate laws.

    So between roughly 0.23 and 0.42 the two disagree **and both are right**: the
    water suffices to write the hydrates, and a real sealed paste still cannot
    reach them. Below the stoichiometric demand they agree, because there the
    limit is mass balance and the Gibbs minimum obeys it like anything else.

    Where that crossover falls is a property of the hydrate assemblage, not a
    constant, so it is measured rather than quoted:

```@example wc_setup
# Below the scanned range the water runs out. What the minimum does then is the
# point of the last two columns.
low = [0.15, 0.20, 0.25, 0.28, 0.30]
clinker(st) = sum(
    ustrip(us"kg", st.n[sp_idx[s]] * cs.species[sp_idx[s]][:M])
        for s in ("C3S", "C2S", "C3A", "C4AF")
)
println(" w/c   clinker left (%)   certified   free water (mol)   x(solvent)   I (mol/kg)")
for wc in low
    fr = fresh_paste(wc)
    # `autostart = false`: the point of this table is a configuration that does
    # NOT certify, and the multi-start cascade is exactly what cannot help there
    # -- it would try every backend, then the ideal pre-solve, then the homotopy
    # continuation, and arrive at the same answer several minutes later. The
    # certificate reported below is the same one; only the search for a better
    # start is declined.
    #
    # The warnings are what the table reports; they are not the transcript.
    eq, cert = Base.CoreLogging.with_logger(Base.CoreLogging.NullLogger()) do
        equilibrate_certified(deepcopy(fr); autostart = false)
    end
    @printf(
        "%5.2f   %16.1f   %9s   %16.3e   %10.3f   %10.4g\n",
        wc, 100 * clinker(eq) / clinker(fr), cert.optimal,
        ustrip(us"mol", eq.n[sp_idx["H2O@"]]), solvent_fraction(eq),
        ionic_strength(eq),
    )
end
```

!!! warning "Read the last three columns: this is outside the model, not inside it"
    Below w/c = 0.30 the free water does not become small — it goes to **6e-9
    mol**, the solver's floor. The solids take all of it. The solvent then holds
    barely a **fifth** of its own aqueous phase, and the ionic strength is
    reported as **409 mol/kg** by a Debye-Huckel model valid to about one.

    That is not the solver failing to converge, and "hydrate a little and leave a
    large stock of anhydrous clinker" is not what a Gibbs minimum does here.
    Forming more hydrate always lowers the energy, and nothing in the model
    penalizes a solution concentrated past any physical meaning: the water
    activity of a real paste collapses as the pores empty and stops the reaction,
    while an activity model extrapolated to 409 mol/kg goes on returning finite
    numbers. The minimization runs off the end of its own domain, and the
    certificate cannot see it — a certificate proves the composition minimizes
    the problem *as posed*, not that the problem was posed inside the model.

    So this table is read for **one** thing: residual clinker exists in the Gibbs
    minimum below the stoichiometric water demand, which is a mass balance — 15 g
    of water cannot hydrate 100 g of cement whatever the algorithm, since the
    hydrates would need some 23 g. Its molar amounts, its pH, its ionic strength
    are **not** equilibrium values and must not be quoted as such. Since 0.15.2
    the package says so itself: [`equilibrate_certified`](@ref) checks
    [`solvent_fraction`](@ref) on the answer it returns and warns — or raises,
    under `STRICT_CONVERGENCE[]` — when the aqueous phase has effectively
    vanished.

    The previous version of this page said no clinker survives *at any* w/c and
    explained it by the minimum "always forming a less hydrous assemblage". The
    first half is true only over the range scanned, and the second is false: the
    least hydrous assemblage available still binds water, and when there is not
    enough, alite stays — and then, shortly after, the model stops applying.

    What remains true is the rest of the original claim, and it matters for mix
    design: over 0.30–0.60 this scan shows **no optimum w/c and no inflection**.
    The porosity rises monotonically and the minimum-porosity mix design does not
    appear, because it is set by the degree of hydration a paste actually
    reaches, not by the assemblage it would reach given time.

### A usable answer below the stoichiometric demand

A high-performance concrete is mixed at w/c between 0.25 and 0.35. That regime is
ordinary, not pathological, and it must be computable — with a certificate, and
with all three of the things one wants from it: how much clinker stays
unhydrated, which hydrates form, and what the pore solution ends up containing.

The way to get it is to stop the reaction where the physics stops it, instead of
asking the minimizer to discover an arrest point it has no term for. React a
fraction ``\alpha`` of the clinker with **all** the water; the rest stays
unhydrated, and the equilibrium is then computed on a system that still has a
solution in it. ``\alpha`` is exactly what [`powers_alpha_max`](@ref) supplies,
and imposing the reacted fraction is the standard construction of cement
thermodynamic modeling — it is how [LothenbachWinnefeld2006](@cite) computes a hydrating
paste.

```@example wc_setup
function arrested(wc, α)
    mtot = c + wc * c
    st = ChemicalState(cs)
    for (sym, mfrac) in compo
        set_quantity!(st, sym, α * mfrac / mtot * u"kg")   # only α reacts
    end
    set_quantity!(st, "H2O@", wc * c / mtot * u"kg")       # all the water
    V = volume(st)
    set_quantity!(st, "H+", 1e-7u"mol/L" * V.liquid)
    set_quantity!(st, "OH-", 1e-7u"mol/L" * V.liquid)
    return st
end

println(" w/c   alpha   certified   x(solvent)   I (mol/kg)     pH   clinker %   porosity")
for wc in (0.25, 0.30, 0.35, 0.42)
    α        = powers_alpha_max(wc)
    fresh    = fresh_paste(wc)                 # the volume reference, all of it
    eq, cert = equilibrate_certified(arrested(wc, α))
    # The unreacted clinker is put back for the volume and porosity accounting.
    n = collect(eq.n)
    for (sym, _) in compo
        n[sp_idx[sym]] += (1 - α) * fresh.n[sp_idx[sym]]
    end
    final = ChemicalState(cs, n)
    @printf(
        "%5.2f  %6.3f   %9s   %10.4f   %10.4f  %5.2f   %9.1f   %8.4f\n",
        wc, α, cert.optimal, solvent_fraction(eq), ionic_strength(eq), pH(eq),
        100 * clinker(final) / clinker(fresh), porosity(final, fresh).total,
    )
end
```

Every point certifies, the solvent holds 0.999 of its phase, and the ionic
strength is 0.035 mol/kg — a pore solution, not the 409 mol/kg of the table
above. The residual clinker is ``1 - \alpha`` by construction, which is the
point: the water limit enters as the closure it physically is, and everything
else is then a proved Gibbs minimum.

!!! note "What would remove the ``\alpha``"
    Predicting the arrest point instead of imposing it is a well-posed
    thermodynamic question, and this package cannot answer it yet. Two ingredients
    are missing, and both are about water that is present but unavailable. An
    **activity model valid at very high concentration** — Pitzer-class — because
    what physically stops hydration is the collapse of the water activity as the
    last of the pore solution is consumed, and an extended Debye-Huckel model
    extrapolated to 409 mol/kg goes on returning finite numbers instead of
    collapsing. And a **coupling between pore structure and water activity**, the
    Kelvin term, because in a fine pore water is held at a reduced activity
    whatever its composition; that is what self-desiccation is, and it is
    poromechanics, not solution chemistry.

    Until then, ``\alpha`` is not a fudge: it is where the missing physics is
    parameterized, measured on real pastes, and the calculation downstream of it
    is proved.

### Assumptions behind these numbers

  - **Complete reaction.** Full equilibrium, no time, no kinetic barrier.
  - **Sealed curing.** No water exchanged with the outside, so the empty
    porosity from the chemical shrinkage stays empty. An immersed specimen would
    draw water in and `ϕ.void` would fill.
  - **The fresh volume is the reference**, and it is held fixed — the specimen
    keeps its cast dimensions, the contraction showing up as internal void rather
    than as shrinkage of the outside. This is the usual convention for a set
    paste; it is wrong before setting, when the material still contracts
    externally.
  - **Ideal molar volumes.** Phase volumes are the sum of ``n_i V_i^0``, with no
    mixing term.
  - **Dilute solution model**, so no ionic-strength correction. The pore solution
    of a cement paste is around 0.1–0.3 mol/kg, where activity coefficients
    depart from unity by tens of percent — use [`HKFActivityModel`](@ref) or
    [`DaviesActivityModel`](@ref) if the ion concentrations themselves matter.
    The pH being buffered, it is the quantity least affected by this choice.
  - **The species list is closed.** Only the 14 species selected above may form;
    siliceous hydrogarnet, hydrotalcite and the alkali sulfates are absent, as
    are the alkalis themselves, which in a real paste raise the pore-solution pH
    to 13 or above.

## The other boundary condition: curing under water

The assumption list says an immersed specimen "would draw water in and `ϕ.void`
would fill". That is not a thought experiment — it is a constraint, and it is the
mirror image of the sealed convention everything above uses:

```@example wc_setup
wc_cure = 0.40
fresh_c = fresh_paste(wc_cure)

# Sealed: the ceiling is 0.42-based, the shrinkage volume stays empty.
α_sealed = powers_alpha_max(wc_cure)
eq_s, c_s = equilibrate_certified(arrested(wc_cure, α_sealed))

# Cured: the ceiling is 0.36-based because the bath refills what the reaction
# empties, AND the specimen is genuinely open to water -- `SaturatedCuring`
# holds its total volume at the fresh paste's and reports how much it drank.
α_cured = powers_alpha_max(wc_cure; curing = :saturated)
q = Ref(Float64[])
eq_c, c_c = equilibrate_certified(
    arrested(wc_cure, α_cured);
    constraint = SaturatedCuring(; reference = fresh_c), parameters = q,
)

function with_unreacted(eq, α)
    n = collect(eq.n)
    for (sym, _) in compo
        n[sp_idx[sym]] += (1 - α) * fresh_c.n[sp_idx[sym]]
    end
    return ChemicalState(cs, n)
end

@printf("%-10s %7s %10s %12s %14s %12s\n",
        "curing", "alpha", "certified", "pH", "water in (mol)", "void")
@printf("%-10s %7.3f %10s %12.2f %14s %12.4f\n",
        "sealed", α_sealed, c_s.optimal, pH(eq_s), "—",
        porosity(with_unreacted(eq_s, α_sealed), fresh_c).void)
@printf("%-10s %7.3f %10s %12.2f %14.4f %12.4f\n",
        "under water", α_cured, c_c.optimal, pH(eq_c), only(q[]),
        porosity(with_unreacted(eq_c, α_cured), fresh_c).void)
```

Two things to read there. The **ceiling moves**, so more of the clinker reacts —
that is [`powers_alpha_max`](@ref) and its `curing` keyword, and it would apply to
a kinetic run just as well. And the **void closes**, because the water drawn in
occupies the volume the chemical shrinkage emptied.

### The chemical shrinkage is computed, not supplied

That last column is worth dwelling on, because it is not an input anywhere in
this calculation.

**What the quantity is.** Hydration products are denser than the reagents that
made them, so a paste occupies less volume after reacting than before. Le
Chatelier measured it in 1900 by watching a sealed flask draw water in. Here it
is

```math
\Delta V = \sum_i \bar V_i\, n_i^{\text{fresh}} - \sum_i \bar V_i\, n_i^{\text{hydrated}} ,
```

a difference of **standard molar volumes** read from the database — the same
thermodynamic data that fixed the assemblage, and the same ``\bar V_i`` that
[`volume`](@ref) and [`porosity`](@ref) use everywhere else on this page. No
calorimetry, no shrinkage test, nothing fitted.

**What Powers says about the same quantity.** His two complete-hydration ratios
differ by ``0.42 - 0.36 = 0.06`` g of water per gram of cement, and that gap *is*
this quantity: the water an immersed specimen takes up and a sealed one must find
in itself. He measured it on pastes in 1948.

So the package can be asked a question it was never fitted to answer: **does it
reproduce the coefficient Powers measured?**

```@example wc_setup
# The molar mass comes from the database, never from a table typed here.
M_H2O = ustrip(us"g/mol", cs.species[sp_idx["H2O@"]][:M])
m_cement = 1000 / (1 + wc_cure)         # g: `fresh_paste` normalizes to 1 kg of paste

# The thermodynamic shrinkage: a difference of molar volumes, per gram REACTED.
ΔV = ustrip(uconvert(us"cm^3", volume(fresh_c).total)) -
     ustrip(uconvert(us"cm^3", volume(with_unreacted(eq_s, α_sealed)).total))
shrink_vol = ΔV / (α_sealed * m_cement)

# The same quantity as the water a cured specimen drinks, per gram reacted.
shrink_mass = only(q[]) * M_H2O / (α_cured * m_cement)

@printf("chemical shrinkage, from molar volumes : %.4f cm3 per g of reacted cement\n",
        shrink_vol)
@printf("water drawn in by the cured specimen   : %.4f g   per g of reacted cement\n",
        shrink_mass)
# Powers' two coefficients, read back out of the function that applies them --
# below either one the cap is w/c divided by it, so the page cannot quote a
# number the code does not use.
w_probe = 0.25
k_sealed = w_probe / powers_alpha_max(w_probe)
k_saturated = w_probe / powers_alpha_max(w_probe; curing = :saturated)
@printf("Powers, as the gap between his two coefficients (%.2f - %.2f): %.4f g/g\n",
        k_sealed, k_saturated, k_sealed - k_saturated)
```

#### Is the comparison circular? No, and that is testable rather than arguable

Powers' coefficients do enter the calculation, through `α_sealed = w/c / 0.42`.
An objection follows immediately: if his number sets how much reacts, is the
agreement above anything more than his number coming back out?

It is not, and the reason is the normalization. ``\Delta V`` scales with how much
clinker reacted, and it is divided by **that same reacted mass**, so `α` cancels
to first order. The way to settle it is not to argue but to move `α` and watch:

```@example wc_setup
@printf("%8s %12s %14s\n", "alpha", "certified", "cm3 per g reacted")
for α in (0.60, α_sealed, 1.00)
    e, cert = equilibrate_certified(arrested(wc_cure, α))
    dv = ustrip(uconvert(us"cm^3", volume(fresh_c).total)) -
         ustrip(uconvert(us"cm^3", volume(with_unreacted(e, α)).total))
    @printf("%8.3f %12s %14.4f\n", α, cert.optimal, dv / (α * m_cement))
end
```

The shrinkage per gram of reacted cement is **flat in `α`** — it moves by 0.3 %
while `α` moves by two thirds, and all three points certify. So Powers' 0.42,
which is what fixes `α = 0.952`, is demonstrably not what produces the answer.
What produces it is the table of molar volumes, and the comparison with his
0.06 g/g is therefore a check on that table rather than his own number coming
back out.

#### Every hypothesis behind the number

Stated in full, because a number that reproduces a measurement is worth exactly
what its assumptions are worth:

  - **The species list is closed** — the 14 species declared at the head of this
    page. No alkalis, no siliceous hydrogarnet, no hydrotalcite. A different
    assemblage has a different volume, and this is the assumption with the most
    room in it.
  - **Molar volumes are ideal**: the volume of a phase is ``\sum_i n_i \bar V_i``
    with no excess term, for the solids and for the solution alike. Real mixing
    volumes are small but not zero.
  - **The standard molar volumes are the database's**, at 25 °C and 1 bar, and
    they are pressure-independent in the shipped data — which is also why
    [`FixedVolume`](@ref) refuses a condensed system.
  - **The pore solution is treated with [`DiluteSolutionModel`](@ref)**, so no
    ionic-strength correction enters the volume either.
  - **The reference is the fresh paste's volume**, the specimen keeping its cast
    dimensions; the contraction appears as internal void rather than as external
    shrinkage. That is the usual convention for a set paste and wrong before
    setting.
  - **Powers' 0.06 is itself an average** over the cements he had, and the
    difference of two separately measured ratios, so it carries the uncertainty
    of both.

None of these is tuned. Two of them — the closed species list and the ideal
volumes — are the ones that would move the answer if the agreement were worse,
and they are the ones to revisit first if a reader's own mix disagrees.

Powers' two complete-hydration ratios differ by
`0.42 − 0.36 = 0.06` g of water per gram of cement, and that gap is his statement
of the same quantity: the water an immersed specimen takes up that a sealed one
must find in itself. It was measured on pastes in 1948. The numbers above come
from a table of standard molar volumes and a Gibbs minimization, with no
calorimetry, no shrinkage test and nothing fitted — so the comparison is a real
check on the volume data rather than a restatement.

The two are under no obligation to agree, and the reasons they need not are worth
naming: the assemblage here is the 14 declared species and not a real paste's,
the molar volumes are ideal with no mixing term, and Powers' coefficient is an
average over the cements he had in 1948. Whatever comes out is therefore a real
check on the volume data rather than a restatement of it — and it is what makes
an empirical coefficient *intelligible* rather than merely used.

Read the `void` column of the table above alongside it, and read it for what it
is. Sealed, that fraction of the specimen is gas-filled porosity, and it is a
**result** — where a sealed paste's self-desiccation comes from. Under water it
is zero **by construction**, since holding the total volume is exactly what the
constraint does. The independent quantity there is the column beside it: how much
water it took, which is what the comparison above is about.

!!! warning "Curing is not `FixedActivity(\"H2O@\", 1.0)`"
    A cement pore solution sits near ``a_w = 0.98`` from its dissolved salts
    alone, so prescribing unit water activity would draw water in until the
    solution was dilute enough to reach it — which never happens. A bath does not
    fix the activity inside the specimen; it fixes the **availability**, which is
    a statement about volume. That is why [`SaturatedCuring`](@ref) is written on
    the volume closure and not on an activity.

!!! note "Extending the scan"
    To study **supplementary cementitious materials**, substitute part of the
    clinker and add the corresponding species from `cemdata18` (e.g. `C2ASH8`).
    For **carbonation**, add `CO2@`, `HCO3-`, `CO3-2` and the carbonate phases —
    see the [cement carbonation example](@ref sec-cement-carbonation).

## Where to go next

The sweep above has no arrest criterion; where a sealed paste stops hydrating,
and why, is the subject of [Self-desiccation](@ref sec-self-desiccation). The
equilibrium view continues with
[Carbonation of a Cement Paste](@ref sec-cement-carbonation).
