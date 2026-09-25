# [Chemical Kinetics](@id sec-kinetics)

!!! info "Before this page"
    The tutorial [Chemical Equilibrium](@ref sec-equilibrium), since every step
    of a kinetic run re-equilibrates the solution, and the theory page
    [Rate laws, and every parameter in them](@ref sec-theory-kinetics) for the
    rate laws themselves.

This tutorial runs a complete kinetic calculation, the hydration of a Portland
clinker, in which the dissolution and precipitation of minerals are governed by
rate laws while the aqueous speciation, much faster, is kept at equilibrium,
following the methodology of [Leal2017](@cite).

## Background

The kinetics algorithm solves:

```math
\frac{d n_k}{dt} = r_k(t), \quad k \in \text{kinetic minerals}
```

where ``r_k`` [mol/s] is the net rate of reaction ``k`` (positive = dissolution,
negative stoichiometric coefficient for the mineral so `dn/dt < 0`).
The aqueous speciation can be re-equilibrated once per accepted step, which
provides the activity coefficients entering the saturation ratio
``\Omega = \text{IAP}/K``; how this coupling is arranged is described in
[The equilibrium–kinetics coupling](@ref) below.

## Two routes, and the one followed here

A kinetic calculation advances in time by one of two routes. The implicit step,
[`kinetic_step`](@ref), performs one Gibbs minimization per time step with the
reaction extents among the unknowns, so that the products are decided by
thermodynamics; the integration of an ordinary differential equation,
[`integrate`](@ref), follows the stoichiometry written in the kinetic reactions
and leaves the step length to an adaptive solver. The choice between them, their
respective limits and the implicit step in full are the subject of
[Which route: two ways to advance in time](@ref sec-kinetics-routes). The run
below follows the second route on a Portland clinker, and every call it makes is
described in [Writing a kinetic model](@ref sec-kinetics-syntax).

## Full workflow: OPC clinker hydration with KineticsProblem

Complete chain — [`ChemicalSystem`](@ref) → [`ChemicalState`](@ref) →
[`KineticReaction`](@ref) → [`KineticsProblem`](@ref) → [`integrate`](@ref):

```julia
using ChemistryLab, OrdinaryDiffEq, DynamicQuantities, Printf

# ── 1. ChemicalSystem from CEMDATA18 ────────────────────────────────────────
DATA_FILE = datapath("cemdata18-thermofun.json")
substances = build_species(DATA_FILE)

input_species = split(
    "C3S C2S C3A C4AF " *
    "Portlandite Jennite ettringite monosulphate12 C3AH6 C3FH6 H2O@",
)
species = speciation(substances, input_species; aggregate_state = [AS_AQUEOUS])
cs = ChemicalSystem(species, CEMDATA_PRIMARIES)

# ── 2. Initial state: 1 kg OPC (CEM I 52.5 R), w/c = 0.40 ─────────────────
WC          = 0.40
COMPOSITION = (C3S=0.619, C2S=0.165, C3A=0.080, C4AF=0.087)   # assumed, typical of the class

state0 = ChemicalState(cs)
for (name, frac) in pairs(COMPOSITION)
    set_quantity!(state0, string(name), frac * u"kg")
end
set_quantity!(state0, "H2O@", WC * u"kg")

# ── 3. [ParrotKilloh1984](@cite) rate functions with Powers α_max ───────────────────────
α_max   = powers_alpha_max(WC)
BLAINE  = 380.0u"m^2/kg"
pk_C3S  = parrot_killoh_avrami(PK84_PARAMS_C3S,  "C3S";  α_max, blaine = BLAINE)
pk_C2S  = parrot_killoh_avrami(PK84_PARAMS_C2S,  "C2S";  α_max, blaine = BLAINE)
pk_C3A  = parrot_killoh_avrami(PK84_PARAMS_C3A,  "C3A";  α_max, blaine = BLAINE)
pk_C4AF = parrot_killoh_avrami(PK84_PARAMS_C4AF, "C4AF"; α_max, blaine = BLAINE)

# ── 4. Kinetic reactions (reaction-centric) ─────────────────────────────────
# Reactions follow [LothenbachWinnefeld2006](@cite) — Jennite = Ca₉Si₆O₁₈(OH)₆·8H₂O
# Balanced hydration reactions — ΔᵣH⁰ computed from species ΔₐH⁰.
sp(name) = cs[name]

rxn_C3S = Reaction(
    OrderedDict(sp("C3S") => 1.0, sp("H2O@") => 103/30),
    OrderedDict(sp("Jennite") => 1.0, sp("Portlandite") => 4/3);
    symbol = "C₃S hydration",
)
rxn_C3S[:rate] = pk_C3S

rxn_C2S = Reaction(
    OrderedDict(sp("C2S") => 1.0, sp("H2O@") => 73/30),
    OrderedDict(sp("Jennite") => 1.0, sp("Portlandite") => 1/3);
    symbol = "C₂S hydration",
)
rxn_C2S[:rate] = pk_C2S

rxn_C3A = Reaction(
    OrderedDict(sp("C3A") => 1.0, sp("H2O@") => 6.0),
    OrderedDict(sp("C3AH6") => 1.0);
    symbol = "C₃A hydration",
)
rxn_C3A[:rate] = pk_C3A

rxn_C4AF = Reaction(
    OrderedDict(sp("C4AF") => 1.0, sp("Portlandite") => 2.0, sp("H2O@") => 10.0),
    OrderedDict(sp("C3AH6") => 1.0, sp("C3FH6") => 1.0);
    symbol = "C₄AF hydration",
)
rxn_C4AF[:rate] = pk_C4AF

# ── 5. Problem + semi-adiabatic calorimeter ─────────────────────────────────
cal = SemiAdiabaticCalorimeter(;
    Cp        = (1.0 * 800.0 + WC * 4186.0 + 1.0 * 900.0) * u"J/K",
    T_env     = 293.15u"K",
    heat_loss = ΔT -> 0.30 * ΔT + 0.003 * ΔT^2,
    T0        = 293.15u"K",
)

kp = KineticsProblem(
    cs, [rxn_C3S, rxn_C2S, rxn_C3A, rxn_C4AF], state0, (0.0, 7.0 * 86400.0);
    calorimeter = cal,
    equilibrium_solver = nothing,
)

# ── 6. Integrate and post-process ───────────────────────────────────────────
ks  = KineticsSolver(; ode_solver = Rodas5P(), reltol = 1e-6, abstol = 1e-9)
sol = integrate(kp, ks)

_, T_vec = temperature_profile(sol, cal)
_, Q_vec = cumulative_heat(sol, cal)

n0_kin = [sol.prob.p.n_initial_full[i] for i in kp.idx_kinetic]
n_kin  = [[u[i] for u in sol.u] for i in eachindex(n0_kin)]

function phase_alpha(cs, kp, n0_kin, n_kin, name)
    sp_idx = findfirst(s -> ChemistryLab.symbol(s) == name, cs.species)
    pos    = findfirst(==(sp_idx), kp.idx_kinetic)
    isnothing(pos) && return fill(NaN, length(sol.t))
    return 1.0 .- n_kin[pos] ./ n0_kin[pos]
end

α_C3S = phase_alpha(cs, kp, n0_kin, n_kin, "C3S")

@printf "ΔT_max = %.2f °C   Q_7d = %.1f kJ/kg   α(C3S) = %.4f\n" \
    maximum(T_vec .- 273.15) - 20.0  Q_vec[end]/1000  α_C3S[end]
```

!!! tip "Choosing `equilibrium_solver`"
    Setting `equilibrium_solver = nothing` skips the Gibbs minimization, which is
    appropriate for the [ParrotKilloh1984](@cite) model: its rate closure ignores
    the `lna` argument, so re-speciation cannot change `α(t)` or the calorimetry.
    It does change what the *products* are — with `nothing`, the hydrate
    assemblage is whatever the hand-written reaction stoichiometry says, not what
    thermodynamics gives. For `transition_state` models, whose rate depends on the
    saturation ratio `Ω = IAP/K`, an `EquilibriumSolver` is required for the rates
    themselves to mean anything.

    The solver may be passed on the [`KineticsProblem`](@ref) or on the
    [`KineticsSolver`](@ref); both work. If it is set on both, the one on the
    problem is used and the conflict is reported.

## The equilibrium–kinetics coupling

The coupling is the partitioned formulation of [Leal2017](@cite), the one
Reaktoro implements. Species are split into a **kinetic partition** — the
minerals carrying a rate law — and an **equilibrium partition**, everything
else: the aqueous phase and any mineral free to precipitate or dissolve
instantaneously.

The ODE state is `(bₑ, nₖ)`: the element amounts held by the equilibrium
partition, and the moles of the kinetic minerals. It advances as

```math
\frac{\mathrm{d} n_k}{\mathrm{d} t} = \nu_k^{\mathsf T} r,
\qquad
\frac{\mathrm{d} b_e}{\mathrm{d} t} = A_e \, \nu_e^{\mathsf T} r ,
```

and the composition of the equilibrium partition is recovered at each step by

```math
n_e = \varphi(b_e) \;=\; \arg\min_{n} \; G(n)
\quad \text{s.t.} \quad A_e\, n = b_e , \; n \ge 0 .
```

Three points are worth stating, because each is a way the coupling can be got
wrong and look plausible:

**The minimization runs over the equilibrium partition only.** Posing it on the
whole system would equilibrate the kinetic minerals instantaneously, which is
what a kinetic description exists to prevent. `KineticsProblem` therefore builds
a sub-system restricted to the partition, sharing the parent's primary species
so that `bₑ`, `dbₑ/dt` and the solve all live in one conservation basis.

**`bₑ` is integrated, not `nₑ`.** Along the way an individual species may want
to go negative — the generated dissolution reactions are written in `H⁺`, and a
cement paste contains no acid — and it is the minimizer, not the caller, that
redistributes the elements over a feasible set. Element amounts are what is
conserved; species amounts are what is solved for.

**The rate sign is fixed by the stoichiometry.** Each kinetic reaction is
normalized so its controlling mineral carries `ν = −1`: a positive rate is a
dissolution. A reaction generated from the nullspace comes out with an arbitrary
orientation, and taken as-is the ODE grows the clinker instead of consuming it.

!!! note "How the two are coupled in time — operator splitting"
    The equilibrium is **not** solved inside the ODE right-hand side. The ODE
    advances the kinetic minerals with the speciation held frozen, and
    `respeciate!` re-equilibrates the equilibrium partition once per accepted
    step, wired as a `DiscreteCallback`.

    This is deliberate. A stiff solver evaluates the right-hand side many times
    per step and differentiates it to build its Jacobian; a Gibbs minimization in
    there makes the cost per step unpredictable and leaves the Jacobian
    describing a different model from the one being integrated. Splitting the two
    is first-order accurate in the step size — the standard arrangement in
    reactive transport — and keeps residual and Jacobian consistent.

    The element amounts `bₑ` carried by the ODE state are handed to the solver
    as the constraint of the sub-problem, not derived from a starting
    composition. The composition passed alongside is a starting guess only.

!!! warning "A failed re-speciation is reported, not hidden"
    If the equilibrium solve fails, that step keeps its frozen composition, the
    first failure is reported with its exception, and `integrate` states how many
    steps were affected. A run in which re-speciation never succeeded must not
    look like a healthy one.

!!! tip "Calorimetry and ΔᵣH⁰"
    The calorimeter computes the heat generation rate as `q̇ = Σ rᵢ × (−ΔᵣH⁰ᵢ)`,
    where `ΔᵣH⁰ᵢ(T)` is the reaction enthalpy (thermodynamic convention:
    negative = exothermic). It is built automatically from species `ΔₐH⁰`
    properties via `complete_thermo_functions!` — **reactions must be
    mass-balanced** for this computation to be correct. For **custom species**
    that lack a `ΔₐH⁰` entry (GGBS, MK, …), set `:ΔᵣH⁰` directly on the
    reaction: `rxn[:ΔᵣH⁰] = NumericFunc((T,) -> -36_100.0, (:T,), u"J/mol")`.

## [Post-processing a kinetics run](@id kinetics-postprocessing)

The ODE state vector carries only the **kinetic** species. Every other amount is
held in a single buffer that the integrator mutates in place, so after a run it
reflects the last accepted step and nothing else. Recovering the composition at an
arbitrary time therefore means replaying the stoichiometry, which is what
[`reaction_extents`](@ref) and [`state_at`](@ref) do.

```julia
sol = integrate(kp, ks)

α  = degrees_of_hydration(sol, kp)          # Dict(symbol => α(t))
ᾱ  = mean_degree_of_hydration(sol, kp)      # mass-weighted binder average

ξ  = reaction_extents(sol, kp)              # extent of each reaction, mol
st = state_at(sol, kp, 28 * 86400.0)        # full ChemicalState at 28 days
```

[`state_at`](@ref) rebuilds **every** species from `n = n₀ + νᵀξ`, so the result
satisfies `A(n - n₀) = (Aνᵀ)ξ` to machine precision — whatever the quadrature
error. How well the elements themselves balance is a property of the reactions
you wrote, not of the reconstruction; check it with `maximum(abs, A * ν')`.

!!! warning "Re-speciation is not replayed"
    When the run used an equilibrium solver, the aqueous partition was
    re-speciated at every accepted step, and that redistribution cannot be
    recovered from the stoichiometry. `state_at` then returns the purely kinetic
    reconstruction; call [`equilibrate`](@ref) on it to speciate.

### From moles to volume fractions

[`volume_fractions`](@ref) turns a state into the input a mean-field
homogenization scheme consumes, using the standard molar volumes `V⁰` that
CEMDATA18 supplies for every cement phase.

```julia
groups = [
    "anhydrous" => ["C3S", "C2S", "C3A", "C4AF"],
    "C-S-H"     => "Jennite",
    "CH"        => "Portlandite",
    "AFt"       => "ettringite",
    "water"     => "H2O@",
]
f = volume_fractions(st, groups; reference = state0)
```

Passing `reference` selects the **sealed-curing** convention of
[Lavergne2018](@cite): fractions are referred to the initial volume, held fixed,
and the deficit left by the reactions appears as a `"void"` phase. That void is
the chemical shrinkage — hydration products occupy less space than the reactants
they consume — and it is a genuine phase of the microstructure, not a rounding
error. Without `reference`, fractions are referred to the current volume and the
void does not exist.

```julia
sum(values(f))                              # 1.0, void included
f["void"]                                   # Le Chatelier contraction
chemical_shrinkage(st, state0)              # the same thing, as a volume
```

!!! note "Species without a molar volume are invisible"
    A species carrying no `V⁰` contributes nothing to [`volume`](@ref),
    [`porosity`](@ref) or [`volume_fractions`](@ref) — a custom species added by
    hand, for instance. Call [`missing_molar_volumes`](@ref) before trusting a
    volume balance:

    ```julia
    isempty(missing_molar_volumes(st)) || @warn "incomplete volume balance"
    ```

Individual fractions may be slightly **negative** for aqueous solutes, whose
standard partial molar volumes are negative (electrostriction contracts the
solvent around an ion). They are kept so that `volume_fractions` and
[`volume`](@ref) always agree; grouping folds them back into the liquid.

## Where to go next

The syntax used above — rate functions, rate constants, kinetic reactions,
calorimeters — is described call by call in
[Writing a kinetic model](@ref sec-kinetics-syntax). The coupling between
kinetics and equilibrium is derived in
[Coupling kinetics and equilibrium](@ref sec-coupling), and the application
pages carry it further: [Cement clinker hydration kinetics](@ref) on the same
clinker, [The hydrating paste, end to end](@ref sec-coupled-hydration) with a
hydrate assemblage computed rather than imposed, and
[Calibrating hydration kinetics on measured calorimetry](@ref ex-hydration-calibration)
when the parameters themselves are the question.
