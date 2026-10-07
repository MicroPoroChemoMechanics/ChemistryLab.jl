# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

"""
    KineticsOrdinaryDiffEqExt

Extension activated when `OrdinaryDiffEq` is loaded. Provides the concrete
`ChemistryLab.integrate(::KineticsProblem, ::KineticsSolver; ...)` implementation
using `OrdinaryDiffEq.ODEProblem` and registers `Rodas5P()` as the default solver.

Follows the formulation of [Leal2017](@citet) with partial equilibrium and
optional thermal coupling (isothermal or semi-adiabatic calorimetry).

# Usage

```julia
using ChemistryLab, OrdinaryDiffEq
sol = integrate(kp, KineticsSolver(; ode_solver=Rodas5P()))
# or shortcut (uses default Rodas5P):
sol = integrate(kp)
```
"""
module KineticsOrdinaryDiffEqExt

using OrdinaryDiffEq
using SciMLBase
import ChemistryLab
import ChemistryLab:
    integrate,
    KineticsProblem,
    KineticsSolver,
    build_kinetics_ode,
    build_u0,
    build_kinetics_params,
    _DEFAULT_KINETICS_SOLVER_FACTORY,
    SemiAdiabaticCalorimeter,
    IsothermalCalorimeter,
    respeciate!,
    _with_equilibrium_solver,
    symbol

# ── Concrete integrate implementation ────────────────────────────────────────

"""
    integrate(kp::KineticsProblem, ks::KineticsSolver) -> ODESolution

Integrate the kinetics ODE using `OrdinaryDiffEq`, in the formulation of [Leal2017](@citet).

The ODE function, initial state, and parameters are built from `kp`.
Calorimetry (isothermal or semi-adiabatic) is integrated in the same ODE: under
partial equilibrium the state carries the change of the enthalpy of the cell
and the temperature is the root of its energy balance, solved with the
partition at every evaluation of the right-hand side (see the theory page
*Kinetics under partial equilibrium*).

Default tolerances: `reltol = 1e-8`, `abstol = 1e-10`.

`speciation` says where the right-hand side reads the equilibrium partition:
`:frozen`, as the last accepted step left it, re-speciated once per step;
`:rhs`, solved at every evaluation, with its derivative in the Jacobian; or
`:auto`, the default, which takes `:rhs` when a rate law reads the partition (an
activity or an amount of an equilibrium species, a saturation ratio) or a
semi-adiabatic cell takes its temperature from it, and `:frozen` otherwise,
where splitting is exact. `:rhs` needs the certified solver
of the partition (OptimaSolver, an aqueous phase with `H2O@`); in that mode a
step that leaves the kinetic amounts outside what the system holds is rejected.

A trajectory that reaches such amounts in any mode is returned with the retcode
`Unstable` and a warning, or an error under `STRICT_CONVERGENCE`, never with
`Success`.

# Examples

```julia
using ChemistryLab, OrdinaryDiffEq

ks  = KineticsSolver(; ode_solver=Rodas5P(), reltol=1e-8, abstol=1e-10)
sol = integrate(kp, ks)
```
"""
# `u_modified!` was renamed `derivative_discontinuity!` in SciMLBase; the new
# name where it exists, the old one on the versions the compat bound still
# admits. Chosen once, when the extension loads.
const _mark_modified! = isdefined(SciMLBase, :derivative_discontinuity!) ?
    SciMLBase.derivative_discontinuity! : SciMLBase.u_modified!

function integrate(kp::KineticsProblem, ks::KineticsSolver; speciation::Symbol = :auto, kwargs...)
    speciation in (:auto, :rhs, :frozen) || throw(
        ArgumentError("speciation must be :auto, :rhs or :frozen; got :$speciation"),
    )
    # `KineticsSolver` also carries an equilibrium solver, and its docstring
    # advertises passing one there. Honor it: without this the field is dead
    # and re-speciation silently never happens.
    kp = _with_equilibrium_solver(kp, ks.equilibrium_solver)

    f! = build_kinetics_ode(kp)
    u0 = build_u0(kp)
    p = build_kinetics_params(kp)

    # Warn for missing Cp° when semi-adiabatic
    if kp.calorimeter isa SemiAdiabaticCalorimeter
        missing_cp = String[]
        for (sp, cp_fn) in zip(kp.system.species, p.cp_fns)
            isnothing(cp_fn) && push!(missing_cp, string(symbol(sp)))
        end
        if !isempty(missing_cp)
            shown = join(missing_cp[1:min(5, length(missing_cp))], ", ")
            suffix = length(missing_cp) > 5 ? "…" : ""
            @warn "SemiAdiabaticCalorimeter: variable Cp_total requires Cp° data per " *
                "species. Missing for $(length(missing_cp)) species " *
                "($shown$suffix). Their contribution to Cp_total is treated as zero."
        end
    end

    # Precedence: call-site kwargs beat the solver's, which beat the defaults.
    # `integrate(kp; reltol = …)` forwards through the one-argument method and
    # used to hit a MethodError, this one taking no kwargs at all.
    defaults = (reltol = 1.0e-8, abstol = 1.0e-10)
    merged = merge(merge(defaults, ks.kwargs), values(kwargs))
    # `ode_solver = :auto` hands the choice to `OrdinaryDiffEq`, whose default
    # polyalgorithm detects stiffness at run time and switches. It is offered
    # rather than made the default, and the reason is measured: on the problems
    # here it lands on the same stiff method as `Rodas5P` and returns the same
    # answers, so switching the default would change nothing while removing a
    # reproducible one. `nothing` keeps `Rodas5P`.
    auto_solver = ks.ode_solver === :auto
    solver = (isnothing(ks.ode_solver) || auto_solver) ? Rodas5P() : ks.ode_solver

    prob = ODEProblem(f!, u0, kp.tspan, p)

    # The equilibrium partition is re-speciated once per accepted step by this
    # callback, and, when a rate law or a semi-adiabatic cell reads it, solved
    # at every evaluation of the right-hand side as well (`:rhs`, see
    # `build_kinetics_ode`). The callback never touches `u`, so
    # `save_positions = (false, false)`.
    if p.n_be > 0
        # A non-converged solve is a warning, not an exception, and its result is
        # used anyway — so it never reached `eq_failures`. Count it over the run.
        nonconv0 = ChemistryLab.NONCONVERGED[]
        respeciate!(p, u0)          # start from an equilibrated state
        # The reference of a calorimeter's energy balance, at the partition the
        # right-hand side will solve: the first equilibrium releases nothing.
        ChemistryLab._initialize_cell!(p, u0)
        reads = ChemistryLab._rates_read_speciation(p)
        p.rates_read_speciation[] = reads
        # A semi-adiabatic cell under partial equilibrium takes its temperature
        # from the partition at every evaluation, whatever the rate laws read.
        cell = p.has_T && p.heat_eq
        cell && speciation === :frozen && throw(
            ArgumentError(
                "a semi-adiabatic cell under partial equilibrium solves its temperature with " *
                    "the partition at every evaluation, which `speciation = :frozen` excludes; " *
                    "use `:auto` or `:rhs`.",
            ),
        )
        mode = speciation === :auto ? (reads || cell ? :rhs : :frozen) : speciation
        mode === :rhs && p.eq_dual === nothing && throw(
            ArgumentError(
                "speciation = :rhs solves the partition at every evaluation, with the " *
                    "certified solver: it needs OptimaSolver loaded and an aqueous phase " *
                    "with `H2O@` in the partition. Otherwise integrate by implicit steps, " *
                    "`kinetic_step_adaptive`, or force `speciation = :frozen`.",
            ),
        )
        p.rhs_mode[] = mode
        if mode === :rhs
            # A step leaving the kinetic amounts outside what the system holds is
            # rejected, and the user's own test kept.
            own = get(merged, :isoutofdomain, nothing)
            domain = own === nothing ? ((u, q, t) -> ChemistryLab._kinetic_state_infeasible(q, u)) :
                ((u, q, t) -> ChemistryLab._kinetic_state_infeasible(q, u) || own(u, q, t))
            merged = merge(merged, (isoutofdomain = domain,))
        end
        cb = DiscreteCallback(
            (u, t, integrator) -> true,
            integrator -> begin
                q = integrator.p
                # Flag the trajectory: only speciations computed here belong to
                # the solution. Everything else is a probe.
                q.on_accepted[] = true
                changed = respeciate!(q, integrator.u)
                q.on_accepted[] = false
                # The right-hand side reads what the re-speciation wrote into `p`
                # only when a law that reads the partition is run frozen; the
                # integrator then has to drop what it derived from the old
                # values. Otherwise it is a function of the state alone, and the
                # re-speciation is its warm start.
                reads_frozen = q.rates_read_speciation[] && q.rhs_mode[] === :frozen
                _mark_modified!(integrator, changed && reads_frozen)
            end;
            save_positions = (false, false),
        )
        sol = auto_solver ? solve(prob; callback = cb, merged...) :
            solve(prob, solver; callback = cb, merged...)

        if p.eq_failures[] > 0
            @warn "re-speciation failed on $(p.eq_failures[]) step(s); those steps kept a frozen composition."
        end
        nonconv = ChemistryLab.NONCONVERGED[] - nonconv0
        if nonconv > 0
            res = p.eq_worst_residual[]
            abs_res = p.eq_worst_abs_acc[]
            abs_all = p.eq_worst_abs[]
            @warn """$nonconv equilibrium solve(s) stopped short of the optimizer's \
            tolerance and were used anyway. Judge them on the element balance, not \
            on the retcode: worst |Aₑn − bₑ|∞ over the run was \
            $(round(abs_res; sigdigits = 3)) mol ON THE ACCEPTED STEPS, i.e. on the \
            trajectory itself; $(round(abs_all; sigdigits = 3)) mol counting also the \
            Jacobian probes and rejected steps, which never enter the solution. \
            Read the first figure: 1e-10 mol is machine precision whatever the \
            system, while 1e-2 mol against a 0.3 mol sulfate budget is not. \
            How much it matters depends on the RATE LAWS: `bₑ` is integrated from \
            the rates alone, so a law that reads only its own degree of reaction \
            (Parrott-Killoh, Waller) gives a trajectory independent of the \
            speciation, and this figure then bears on the reported composition \
            only — recover that with `speciated_states`, which certifies each instant against the KKT conditions. A \
            law reading log-activities (a saturation ratio) does feed the \
            speciation back into the trajectory, and there this figure is a \
            direct measure of the error. \
            Do NOT simply loosen the optimizer tolerance — on the calcite reference \
            case that degrades the speciation from 4 % to 250 % against Reaktoro."""
        end
    else
        sol = auto_solver ? solve(prob; merged...) : solve(prob, solver; merged...)
    end
    # A co-reactant is judged to the precision asked of the integration.
    return _flag_infeasible(sol, p, kp; rtol = max(ChemistryLab._FEASIBILITY_RTOL, float(maximum(merged.reltol))))
end

"""
    _flag_infeasible(sol, p, kp; rtol) -> ODESolution

`sol`, with the retcode `Unstable` when a saved state holds kinetic amounts
outside what the system can produce (`ChemistryLab._kinetic_state_infeasible`):
an amount negative beyond rounding, the kinetic species holding more of an
element than the system was given, or, without an equilibrium partition, a
co-reactant consumed beyond what the system holds, by more than `rtol` of the
most of it the element totals allow. A warning names the first such state, or,
under `STRICT_CONVERGENCE`, an error is raised.

A trajectory of that kind used to come back with `Success` and a warning on its
final state alone, judged against twice the total amount of matter, which the
water dominates. Measured on calcite dissolving under `r = k(1 − Ω)` over `10⁵ s`
with the partition frozen within a step, `Rodas5P` ended on 487 mol of calcite
from 0.05: frozen, the rate is constant within a step, and a step longer than the
relaxation of `Ω` overshoots the equilibrium and reverses it. Such a law now has
its partition solved in the right-hand side (`speciation = :rhs`), where this
state is not reached; the check stands for every run.
"""
function _flag_infeasible(sol, p, kp; rtol = ChemistryLab._FEASIBILITY_RTOL)
    k = findfirst(u -> ChemistryLab._kinetic_state_infeasible(p, u; rtol), sol.u)
    k === nothing && return sol
    u = sol.u[k]
    nb = p.n_be
    t_k = round(ChemistryLab._plain(sol.t[k]); sigdigits = 4)
    out = ChemistryLab._exhausted_coreactant(p, u; rtol)
    msg = if out !== nothing && !ChemistryLab._kinetic_amounts_infeasible(p, u)
        # A co-reactant no rate law reads: named, with what the extents take of
        # it beyond what there is by the end of the run.
        name = symbol(kp.system.species[first(out)])
        last_out = ChemistryLab._exhausted_coreactant(p, sol.u[end]; rtol)
        short = last_out === nothing ? -last(out) : -last(last_out)
        "`$name` runs out at t = $t_k s while the reactions consuming it go on: by the end " *
            "of the run their extents take $(round(short; sigdigits = 4)) mol of it beyond what " *
            "the system holds. A rate law that does not read a co-reactant cannot stop when it " *
            "runs out, and the amounts it reads are floored at zero, so the shortfall is created. " *
            "Gate the rate on that co-reactant, write the reaction without it, or attach an " *
            "equilibrium solver so that the assemblage follows from the element budgets. The run " *
            "is returned with the retcode `Unstable`."
    else
        amounts = join(
            (
                "$(symbol(kp.system.species[kp.idx_kinetic[j]])) = " *
                    "$(round(ChemistryLab._plain(u[nb + j]); sigdigits = 4)) mol"
                    for j in 1:(p.n_nk)
            ), ", ",
        )
        "the trajectory reaches kinetic amounts no chemistry can produce from what " *
            "the system holds, at t = $t_k s: $amounts. The run is returned with the " *
            "retcode `Unstable`."
    end
    ChemistryLab._strict_convergence() && throw(ErrorException(msg))
    @warn msg maxlog = 1
    return SciMLBase.successful_retcode(sol) ?
        SciMLBase.solution_new_retcode(sol, SciMLBase.ReturnCode.Unstable) : sol
end

# ── __init__: register default solver ────────────────────────────────────────

function __init__()
    return _DEFAULT_KINETICS_SOLVER_FACTORY[] =
        () -> KineticsSolver(; ode_solver = Rodas5P())
end

end  # module KineticsOrdinaryDiffEqExt
