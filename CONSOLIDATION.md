# Consolidation Roadmap

This roadmap organizes follow-up to the audit of version 0.18.0
(`1fcbb4a5`). Each item belongs on a dedicated branch and is accepted through
a PR: one maintainer approves routine changes; @jfbarthelemy approves scientific
or major changes. An unchecked item is planned work,
not a claim that a fix or scientific validation has been completed.

Track implementation and review in [issue #57](https://github.com/MicroPoroChemoMechanics/ChemistryLab.jl/issues/57).

## Immediate corrections

- [ ] **Safe database imports** — `fix/safe-database-import`.
  Interpret classification labels and unit expressions as data. Preserve the
  bundled formats and test that imported text cannot execute instructions.
- [ ] **Dimensional checks** — `fix/dimensional-checks`.
  Preserve errors for dimensioned arguments to dimensionless mathematical
  functions. Replace global extensions with explicit internal conversions.
- [ ] **Species identity** — `fix/species-identity`.
  Agree on identity versus composition before implementation. Align equality
  and hashing; cover polymorphs, aliases, phase instances, sets, and dictionaries.
- [ ] **Solver concurrency** — `fix/solver-concurrency`.
  Move temporary policies and diagnostics into each solve's context. Test
  overlapping tasks and verify that one solve cannot alter another's policy.
- [ ] **Water autoprotolysis** — `fix/water-autoprotolysis`.
  Reproduce the audit failure in `test/equilibrium_reference.jl` with
  `nullspace_step=false` on Julia 1.12.7/macOS ARM64 and OptimaSolver 0.6.0.
  Identify the dependency or numerical cause before changing tolerances.

## Scientific contracts

- [ ] **Energy, potentials, and certificates** —
  `refactor/thermodynamic-consistency`. Define the formulation and guarantee for
  each activity model and constraint. Check Gibbs-Duhem and gradient consistency;
  distinguish approximate equilibrium residuals from global optimality claims.
- [ ] **Results and provenance** — `refactor/solve-results`.
  Design a result carrying the reference budget, model, solver, residuals,
  warnings, and certificate scope, with an explicit compatibility plan.
- [ ] **Kinetic and thermal coupling** — `fix/coupled-kinetics`.
  Review rate-law dependencies, splitting error, and the complete energy balance.
  Validate saturation-dependent rates, exhausted reactants, and phase changes.

## Validation and maintenance

- [ ] **Extension and platform coverage** — `test/extension-matrix`.
  Exercise optional extensions separately and cover supported platforms. Test
  physical invariants rather than requiring an external solver to stay wrong.
- [ ] **Data validity and reference cases** — `test/model-validation`.
  Preserve parameter provenance and validity ranges. Extend independent checks
  to mixtures, trace elements, phase transitions, and kinetic trajectories.
- [ ] **Measured performance** — `perf/representative-benchmarks`.
  Record loading, allocations, and solve times for representative problems before
  optimizing types, buffers, or subsystem construction.

## Acceptance evidence

Each PR identifies its baseline, reproduction, final behavior, tests, and
limitations. Scientific changes include assumptions and before/after numerical
results. Record unresolved failures explicitly; passing syntax or documentation
reference checks does not establish numerical correctness.

The initial audit observed an equilibrium-test failure. The final outcome of a
separate utility/kinetics run could not be recovered after session interruption;
do not cite that run as passing. Re-run relevant validation for each PR.
