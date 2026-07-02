# Changelog

## v0.3.0 — Chemical kinetics module

### New features

**`src/kinetics/` — mineral dissolution/precipitation kinetics**

- `KineticReaction` — couples a reaction to a rate law; supports
  `transition_state`, `first_order_rate`, and the empirical
  Parrot–Killoh (1984) model for cement clinker hydration
- `KineticsProblem` / `KineticsSolver` — ODE problem formulation
  following the SciML `(u, p, t)` convention; integrated via
  `KineticsOrdinaryDiffEqExt` (weakdep, activated by `using OrdinaryDiffEq`)
- `SemiAdiabaticCalorimeter` / `IsothermalCalorimeter` — coupled
  heat-balance ODE; `temperature_profile`, `heat_flow`,
  `cumulative_heat` post-processing helpers
- `StateView{T,I}` — O(1) named access to species data vectors in the
  ODE hot path (no per-step allocation)
- `RateMechanism`, `RateModelCatalyst`, `BETSurfaceArea`,
  `FixedSurfaceArea` — building blocks for custom rate closures
- `parrot_killoh(params, mineral_name)` factory with built-in
  Schindler & Follliard (2005) Arrhenius correction and default
  parameters for C₃S, C₂S, C₃A, C₄AF
- ForwardDiff-compatible throughout; `KineticFunc` and `transition_state`
  closures accept `Dual` numbers

### Infrastructure
- Relicensed to LGPL-2.1-or-later
- Migrated to Codeberg (`MicroPoroChemoMechanics/ChemistryLab.jl`)
- Forgejo workflows: CI, Documentation, Release, Runic, Zenodo
- Registered in MPCM-Registry (OptimaSolver resolved via registry)
- Multi-version documentation deployment (`docs/deploy_docs.jl`)
- Runic.yml: switch to `workflow_dispatch` only

## v0.2.3 — Activity models & solid solutions

- Extended Debye–Hückel and Davies activity models
- Redlich–Kister solid solution phases
- `SolidSolutionPhase`, `build_solid_solutions` from TOML
- `with_class` for end-member requalification

## v0.2.0 — Equilibrium solver integration

- `EquilibriumProblem` / `EquilibriumSolver` / `equilibrate` API
- `OptimizationIpoptExt` and `OptimaSolverExt` extensions
- Implicit-differentiation sensitivity via OptimaSolver
- HKF aqueous solute thermodynamic model (`NumericFunc`)
- ThermoFun JSON and PHREEQC `.dat` database parsers
