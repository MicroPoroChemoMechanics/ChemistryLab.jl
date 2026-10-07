# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

# ── dual_solver.jl ────────────────────────────────────────────────────────────
#
# Chemical equilibrium by exact Newton on the KKT system, certified.
#
# THE ALGORITHM IS NOT HERE. It lives in `OptimaSolver.dual_newton_solve`, where
# it belongs: the active set on the bound-constrained variables, the degeneracy
# test on the rows of `A`, the two-level Newton and the certificate are all
# statements about a convex program, not about chemistry. What this file supplies
# is the four things that ARE chemistry:
#
#   1. the conservation matrix `A` and the reference potentials `g = Δ_a G⁰/RT`;
#   2. the activity model, as the callback `h` with `∇f = g + h(n)`;
#   3. which species are strictly positive at any solution — the aqueous ones,
#      parameterized by `ln n` — and which may vanish — the pure phases;
#   4. which species plays the reference role: the solvent, whose activity is a
#      mole fraction and therefore bounded above, so that its stationarity cannot
#      be inverted and must be carried by the outer system.
#
# The mathematics is documented in the manual under
# "Proving that an answer is the answer", and the general form in OptimaSolver's
# theory page.

using LinearAlgebra

"""
    DualEquilibriumSolver(system, model; tol, maxit, max_active_updates, si_tol,
                          inner_tol, inner_maxit, inner_fall_bound,
                          lenient_line_search, verbose)

Equilibrium by Newton on the KKT system, in element-potential space.

Where [`EquilibriumSolver`](@ref) minimizes `G` by an interior-point method over
all species and stops on `MaxIters` for a cement, this solves the stationarity
conditions directly and returns a result that
[`optimality_certificate`](@ref) can **prove** optimal — the Gibbs problem being
convex, the KKT conditions are sufficient.

Writing `u = −Aᵀy` for the element potentials, an aqueous species obeys the
mass-action law `aᵢ = exp(uᵢ − gᵢ)` and a pure phase is present exactly when
`uᵢ = gᵢ`, absent when undersaturated: the classical phase-stability criterion.

Starting from the answer of an [`EquilibriumSolver`](@ref) is the intended use —
the interior-point method reaches a neighborhood, this reaches the conditions.

The keywords are OptimaSolver's `DualNewtonOptions`: `maxit` and
`max_active_updates` bound the outer Newton and the active-set search,
`inner_tol` and `inner_maxit` the inner fixed point that recovers the phase
compositions. `tol` and `si_tol` are also the thresholds of the certificate this
solver issues, so loosening them loosens the proof, not only the search.
`inner_fall_bound` (30) and `lenient_line_search` (`false`) change the path of
the search and not the conditions of its answer; the steps of a kinetic run set
them to `Inf` and `true`, which is where they pay (see OptimaSolver's
`DualNewtonOptions`).

See also: [`optimality_certificate`](@ref), [`speciated_states`](@ref).
"""
struct DualEquilibriumSolver{L, M <: AbstractActivityModel, AT <: AbstractMatrix}
    system::ChemicalSystem
    lna::L
    model::M
    idx_aq::Vector{Int}
    idx_pure::Vector{Int}
    j_solvent::Int
    ss_groups::Vector{Vector{Int}}   # end-members of each declared solid solution
    site_groups::Vector{Vector{Int}} # members of each surface site family
    # `A` is the composition matrix with the site-coupling rows APPENDED, not
    # merged into it: `SM.A` is what encodes element conservation, and a family
    # whose sites follow its host adds an equation rather than changing a
    # composition. `n_element_rows` is where the original rows stop, which is
    # also the range the degeneracy criterion may judge — an added row has a
    # right-hand side of zero by construction and entries of both signs, and
    # asking a question meant for an element balance of it gets the wrong
    # answer.
    A::AT
    n_element_rows::Int
    opts::NamedTuple
    # The number type the activity and mixing models carry: `Float64`, or the
    # dual numbers of a differentiation with respect to their parameters. Found
    # once here rather than at every solve.
    captured::Type
end

function DualEquilibriumSolver(
        system::ChemicalSystem,
        model::AbstractActivityModel = DiluteSolutionModel();
        tol::Float64 = 1.0e-10,
        maxit::Int = 200,
        max_active_updates::Int = 200,
        si_tol::Float64 = 1.0e-8,
        inner_tol::Float64 = 1.0e-10,
        inner_maxit::Int = 200,
        inner_fall_bound::Float64 = 30.0,
        lenient_line_search::Bool = false,
        verbose::Bool = false,
    )
    idx_aq = [i for (i, s) in enumerate(system.species) if aggregate_state(s) == AS_AQUEOUS]
    idx_pure = [i for (i, s) in enumerate(system.species) if aggregate_state(s) != AS_AQUEOUS]
    isempty(idx_aq) && throw(
        ArgumentError(
            "DualEquilibriumSolver needs an aqueous phase: the element potentials are " *
                "determined by the mass-action laws of the aqueous species."
        )
    )
    jw = something(findfirst(i -> symbol(system.species[i]) == "H2O@", idx_aq), 0)
    jw == 0 && throw(
        ArgumentError(
            "DualEquilibriumSolver needs `H2O@` among the species: the solvent's " *
                "log activity is bounded above, by zero, so its stationarity " *
                "cannot be inverted and is carried by the outer system instead."
        )
    )
    # A solid solution is a MIXING phase, exactly like the aqueous one: its
    # end-members have composition-dependent activities, so none of them is ever
    # exactly absent while the phase exists, and the active set belongs at the
    # level of the phase. They are therefore removed from the bound-constrained
    # set, where they would have been treated as pure phases.
    ss = system.solid_solutions
    ss_groups = Vector{Int}[]
    if ss !== nothing
        byname = Dict(symbol(sp) => i for (i, sp) in enumerate(system.species))
        for phase in ss
            idx = [byname[symbol(em)] for em in phase.end_members if haskey(byname, symbol(em))]
            length(idx) == length(phase.end_members) && push!(ss_groups, idx)
        end
    end
    # A surface site family is a mixing phase too, and for a sharper reason: its
    # members' activities are site fractions, so none is ever exactly absent
    # while the family has a budget, and the budget is pinned by a conservation
    # row. Left in the bound-constrained set they would be pure phases of unit
    # activity — the site mixing never computed, the balance unsatisfiable.
    site_groups = copy(system.site_groups)
    in_mixing = Set(vcat(ss_groups..., site_groups...))
    idx_pure = [i for i in idx_pure if !(i in in_mixing)]

    # `conservation_matrix`, not `SM.A`: identical when nothing follows its host,
    # and carrying the `−ν` in each coupled family's (site row, host column)
    # entry when something does. The site row then states the coupling instead
    # of a fixed budget, which is what makes the budget follow the host at all.
    A_elem = conservation_matrix(system)
    # A site density being differentiated makes the matrix dual, as a model's
    # parameters make its activities.
    captured = promote_type(
        _captured_number_type(model), _captured_number_type(_MixingTerms(system)), _captured_number_type(A_elem),
    )
    return DualEquilibriumSolver(
        system, _scoped_lna(activity_model(system, model)), model,
        idx_aq, idx_pure, jw, ss_groups, site_groups, A_elem, size(A_elem, 1),
        (; tol, maxit, max_active_updates, si_tol, inner_tol, inner_maxit, inner_fall_bound, lenient_line_search, verbose),
        captured,
    )
end

activity_model(des::DualEquilibriumSolver) = des.model

"""
    _dual_problem(des, p, n0) -> DualNewtonProblem

Package the chemistry as the convex program `OptimaSolver` solves. Built per
solve because the reference potentials `Δ_a G⁰/RT` depend on temperature and
pressure.
"""
function _dual_phases(des::DualEquilibriumSolver, n0, p = nothing; invert = _aqueous_inverter(des))
    # One mixing phase for the aqueous solution — always present, the solvent as
    # its reference — and one more per declared solid solution, whose presence
    # the tangent-plane test decides.
    #
    # The reference is the member whose stationarity the OUTER system carries
    # instead of inverting, so it has to be the one the phase is mostly made of.
    # For the aqueous phase that is the solvent, and not by convention: water is
    # the only species there whose activity is a mole fraction, bounded above, so
    # its law cannot be inverted at all. Inside a solid solution every member is
    # a mole fraction and any of them could serve — but the choice is not
    # indifferent. The outer unknown is `ln x_ref`, and picking a member that
    # happens to be minor sends that unknown towards −∞ and the Jacobian with it.
    # So the reference is chosen by magnitude, from the composition the caller
    # supplied, which is what the solvent already is for the aqueous phase.
    # Every entry carries every field, `split_starts` empty here (the aqueous
    # phase cannot unmix: there is one solvent). The vector is `Any` because
    # the entries differ in the type of `local_h`, nothing or a closure.
    #
    # The aqueous phase carries `invert`, its solutes recovered at once through
    # the ionic strength (`_aqueous_inverter`), when the activity model allows it.
    phases = Any[
        (
            members = des.idx_aq, j_ref = des.j_solvent,
            always_present = true, mole_fraction = false,
            split_starts = Vector{Vector{Float64}}(),
            newton = false, bounded_members = Int[], local_h = nothing, invert = invert,
        ),
    ]
    models = _ss_models(des)
    for (k, grp) in pairs(des.ss_groups)
        j_ref = argmax(@view n0[grp])
        mdl = get(models, k, nothing)
        push!(
            phases,
            (
                members = grp, j_ref = j_ref,
                always_present = false, mole_fraction = true,
                split_starts = _split_starts(mdl, length(grp)),
                newton = _needs_newton_inversion(mdl),
                bounded_members = _bounded_members(mdl),
                local_h = _needs_newton_inversion(mdl) ? _local_log_activities(mdl, p, grp) : nothing,
            ),
        )
    end

    # One more per surface site family, and `always_present = true` is not a
    # convenience here. With `false`, the seeding test asks whether the members
    # already hold more than 1e-9 mol; from a cold start they all sit at the
    # floor, the phase is never activated, and the conservation row
    # `Σ n_member = N_t` with `N_t > 0` is infeasible from the first Newton
    # iteration. `true` also exempts it from the drop rule, which is right: a
    # phase whose total is pinned must not be dropped because its members are
    # small.
    #
    # The reference is the **free site**, the first member by construction. It is
    # the state a fresh surface is mostly in, which is the same criterion the
    # solvent meets for the aqueous phase, and it is stable: a family cannot run
    # out of free sites without the occupied ones taking their place.
    #
    # A site family does not unmix, so no split start is offered.
    for grp in des.site_groups
        push!(
            phases,
            (
                members = grp, j_ref = 1,
                always_present = true, mole_fraction = true,
                split_starts = Vector{Vector{Float64}}(),
                newton = false, bounded_members = Int[], local_h = nothing,
            ),
        )
    end
    return phases
end

# How the solver recovers the composition of a solid solution. The successive
# substitution is exact in one sweep for ideal mixing and contracts under a weak
# excess term; on a sublattice model it diverges (the gel of Myers et al. within
# six sweeps, at the potentials of a CEM I paste), so that phase is inverted by
# Newton's method instead. So is a phase with an excess term: the published AFm
# SO4/OH binary of Cemdata18, whose energy is concave on part of its range, held
# the search of a CEM I paste for minutes with the substitution and never
# certified, while Newton certifies it in one route, at the composition GEMS3K
# gives. Ideal mixing keeps the substitution, and its results bit for bit.
_needs_newton_inversion(::Any) = false
_needs_newton_inversion(::Union{SublatticeModel, CompoundEnergyModel}) = true
_needs_newton_inversion(
    ::Union{RedlichKisterModel, RegularSolutionModel, SubregularSolutionModel, MulticomponentRedlichKisterModel, VanLaarModel}
) = true

# The log activities of a solid solution's members from their own amounts, as
# `_solid_solution_lna!` computes them inside the activity closure (same ϵ, same
# temperature), for the Newton inversion to differentiate instead of the whole
# closure: eight variables instead of a hundred on a cement.
# The members' `ΔₐG⁰/RT` are those of the solve, read once, for a model whose
# activities depend on them (`CompoundEnergyModel`).
function _local_log_activities(mdl, p, grp)
    ϵ = (p !== nothing && hasproperty(p, :ϵ)) ? p.ϵ : _AMOUNT_FLOOR
    T = (p !== nothing && hasproperty(p, :T)) ? p.T : 298.15
    g = (p !== nothing && hasproperty(p, :ΔₐG⁰overRT)) ? collect(p.ΔₐG⁰overRT[grp]) : nothing
    return function (nm)
        tot = sum(nm) + ϵ
        x = nm ./ tot
        return _ss_log_activities!(similar(x), eachindex(x), x, mdl, T, ϵ, g)
    end
end

# The members of a sublattice phase that own no species on any site: their
# activity stays finite as they vanish, so the phase can be present without them
# and their condition is then an inequality (`SublatticeModel`).
_bounded_members(::Any) = Int[]
_bounded_members(m::SublatticeModel) = findall(iszero, m.exponents)
# A compound-energy phase selects the split of its members in which each
# member's activity goes as its mole fraction (`CompoundEnergyModel`): none is
# bounded.
_bounded_members(m::CompoundEnergyModel) = Int[]

"""
    _ss_models(des) -> Dict{Int, Any}

The mixing model of each entry of `des.ss_groups`, by position.

`ss_groups` is built by walking `system.solid_solutions` and **skipping** any
declaration whose end-members are not all in the species list, so the two lists
are the same length only when nothing was skipped. Rebuilding the correspondence
the same way is the only way to be sure a model is matched to its own group; a
positional zip would silently pair a phase with someone else's model the first
time a declaration is dropped.
"""
function _ss_models(des::DualEquilibriumSolver)
    out = Dict{Int, Any}()
    ss = des.system.solid_solutions
    ss === nothing && return out
    byname = Dict(symbol(sp) => i for (i, sp) in enumerate(des.system.species))
    k = 0
    for phase in ss
        idx = [byname[symbol(em)] for em in phase.end_members if haskey(byname, symbol(em))]
        length(idx) == length(phase.end_members) || continue
        k += 1
        out[k] = model(phase)
    end
    return out
end

"""
    _split_starts(model, nmembers) -> Vector{Vector{Float64}}

Where to look for the other lobe of a phase that may unmix, as mole fractions
over its members.

The tangent-plane search in `OptimaSolver` probes the **corners** of the
composition simplex and refines by successive substitution, which converges to
the stationary point nearest its start. A corner is usually on the right side of
the barrier and sometimes is not; where it is not, the iteration walks back to
the phase's own composition and reports nothing, though splitting would lower the
energy.

Measured on the AFm sulfate/hydroxide binary with the published Redlich-Kister
parameters (spinodal [0.631, 0.914], binodal [0.4999, 0.9700]):

| phase sits at | corners alone | with the binodal handed over |
|:--|--:|--:|
| x = 0.5268, where a CEM I settles | +3.99e-02, found | the same verdict |
| x = 0.95 | +2.4e-16, **missed** | +1.23e-01, trial x = 0.444 |
| x = 0.98, outside the binodal | stable | stable |

Both flagged compositions are metastable, so being metastable is not by itself
what defeats the corners — sitting in the lobe they lead back into is, and that
is not knowable in advance. Hence a start supplied unconditionally rather than a
rule for when to supply one.

[`common_tangent`](@ref) computes the binodal from the mixing model alone, in
microseconds and with no reference to the rest of the system, and the search then
refines it with the chemical potentials the system actually has — which is the
pair that matters. That division is the whole point: the model's binodal is a
good *place to look*, not the answer. Extra starts can only raise the maximum the
search returns, so they never take a verdict away and, as the last row shows, do
not invent one.

Empty for anything but a binary with a gap, which is the only case
`common_tangent` is defined for.
"""
function _split_starts(model, nmembers::Int)
    out = Vector{Vector{Float64}}()
    (model === nothing || nmembers != 2) && return out
    pair = try
        common_tangent(model)
    catch
        nothing
    end
    pair === nothing && return out
    for x in pair
        (isfinite(x) && 0 < x < 1) || continue
        push!(out, Float64[x, 1 - x])   # x is the mole fraction of the first end-member
    end
    return out
end

function _dual_problem(des::DualEquilibriumSolver, p, n0, blocks = nothing)
    bl = blocks === nothing ?
        _constraint_blocks(FixedTP(), des, nothing, p, n0) : blocks
    # The inversion of the aqueous phase has to see the parameters the activity
    # model sees. A constraint whose `hq` changes them (the temperature of an
    # adiabatic solve) says how with `pq`; one that does not say is left to the
    # sweeps, whatever it changes.
    pq = get(bl, :pq, nothing)
    invert = (bl.hq === nothing || pq !== nothing) ? _aqueous_inverter(des, pq) : nothing
    phases = _dual_phases(des, n0, p; invert)
    # `float`, not `Float64`: the potentials of a state carrying dual numbers
    # (a temperature being differentiated) keep them.
    return _optima_dual_problem(
        des.A, float.(p.ΔₐG⁰overRT), des.lna, phases, des.idx_pure, p,
        bl.gq, bl.hq, bl.cq, bl.q0, bl.qscale, bl.Aq, Int[],
        1:des.n_element_rows,
    )
end

"""
    SciMLBase.solve(des::DualEquilibriumSolver, state; b = nothing, ϵ = 1e-16)
        -> ChemicalState

Equilibrium composition. `state` supplies the temperature, the pressure and the
starting guess; `b` the component totals, defaulting to those of `state`.

The amounts are returned as computed, without a floor: an amount of `e⁻³⁰⁰` IS
the mass-action answer for a species that is not there, and raising it to `ϵ`
falsifies its activity by hundreds of `RT` units — which is enough to make
[`optimality_certificate`](@ref) report a residual of 74 on a composition solved
to `5e-12`.
"""
function SciMLBase.solve(
        des::DualEquilibriumSolver,
        state::ChemicalState;
        b = nothing,
        ϵ::Float64 = _AMOUNT_FLOOR,
        constraint::EquilibriumConstraint = FixedTP(),
        parameters::Union{Nothing, Base.RefValue} = nothing,
        surface_potential::Symbol = :auto,
    )
    surface_potential in (:auto, :unknown, :eliminated) || throw(
        ArgumentError(
            "surface_potential must be :auto, :unknown or :eliminated; " *
                "got :$surface_potential"
        ),
    )
    # A state or a budget carrying dual numbers: solved on their values, one
    # level of duals down, and the answer lifted to the caller's duals.
    D = _input_number_type(state, b; constraint = constraint, captured = des.captured)
    if D <: ForwardDiff.Dual
        Tg = ForwardDiff.tagtype(D)
        qv = Ref{Any}(Float64[])
        eq = with(_STRIP_TAGS => (_STRIP_TAGS[]..., Tg)) do
            SciMLBase.solve(
                _rebuilt(des), _strip_state(state, Tg); b = b === nothing ? nothing : _strip_tag(collect(b), Tg),
                ϵ = ϵ, constraint = constraint, parameters = qv, surface_potential = surface_potential,
            )
        end
        bd = b === nothing ? des.A * _build_n0(state) : collect(b)
        eq_d, q_d = _lift_equilibrium(
            des, state, eq, bd; ϵ = ϵ, constraint = constraint, q = qv[],
            surface_potential = surface_potential, strip_tag = Tg,
        )
        parameters === nothing || (parameters[] = q_d)
        return eq_d
    end
    p = _build_params(state; ϵ = ϵ)
    n0 = Float64[ustrip(us"mol", x) for x in state.n]
    bv = b === nothing ? des.A * n0 : Float64.(collect(b))

    blocks = _constraint_blocks(constraint, des, state, p, n0)
    if surface_potential !== :eliminated
        surf = _surface_potential_blocks(des, state, p, n0)
        surface_potential === :unknown && surf === nothing && throw(
            ArgumentError(
                "surface_potential = :unknown was asked for, and no site family of " *
                    "this system needs one. Only a `DiffuseLayer` or `ChargePlanes` does."
            ),
        )
        blocks = _compose_blocks(blocks, surf)
    end
    blocks = _scoped_blocks(blocks)
    prob = _dual_problem(des, p, n0, blocks)
    # The solve starts with no amount below `ϵ`, the amount of an absent species
    # in a cold state: an exact zero would start it at `exp(−700)`, from which a
    # solute climbs thirty log units a sweep. The budget stays that of the state.
    res = _optima_dual_solve(prob, bv, max.(n0, ϵ), des.opts)

    res.converged || begin
        Threads.atomic_add!(NONCONVERGED, 1)
        # Silent while a multi-start route is trying candidates: one of them not
        # converging is what the search is for, and the verdict belongs to the
        # certificate of the answer, not to a candidate.
        _EXPLORING_STARTS[] || @warn """the dual equilibrium solve did not certify \
        optimality; audit it with `optimality_certificate`.""" maxlog = 1
    end

    # The parameters the constrained solve FOUND — the temperature an adiabatic
    # solve reached, the titrant amount a prescribed pH required. They are part of
    # the answer, not a diagnostic, so the caller can ask for them.
    parameters === nothing || (parameters[] = copy(res.q))

    T_out, P_out = blocks.apply(temperature(state), pressure(state), res.q)
    x = res.converged ? _complete_floored_solutes(des, prob, res, _activity_floor(p)) : res.x
    return ChemicalState(
        des.system, [nᵢ * u"mol" for nᵢ in x];
        T = T_out, P = P_out,
    )
end

# ── differentiating through the certified solve ──────────────────────────────

_carries_duals(b) = b !== nothing && any(x -> x isa ForwardDiff.Dual, b)

# The same solver, built again in the current scope: inside `_STRIP_TAGS` its
# activity closure returns values.
_rebuilt(des::DualEquilibriumSolver) = DualEquilibriumSolver(des.system, des.model; des.opts...)

"""
    _lift_equilibrium(des, state, eq, b; ϵ, constraint, q, surface_potential, strip_tag,
                      floor = _CERTIFICATE_FLOOR) -> (ChemicalState, q)

The equilibrium `eq`, found on the values of the data, carrying the derivatives
the dual numbers of `state` (its temperature and pressure) and of the budget `b`
imply. The problem is the one `solve` poses, constraint and surface unknowns
included, built at the caller's duals; its answer is lifted by the
implicit-function theorem at `eq` with the active set frozen (OptimaSolver's
`dual_newton_tangent`), so the derivatives are those of the solution, not of the
iteration that found it, and each level of a nested differentiation is taken in
turn. `q` is what the solve found for the constraint's unknowns, returned with
its derivatives. A pure phase holding `floor` or less is absent: the floor of the
certificate for a certified answer, the solver's lower bound for an
interior-point one, which leaves an absent phase there.
"""
function _lift_equilibrium(
        des::DualEquilibriumSolver, state::ChemicalState, eq::ChemicalState, b;
        ϵ::Float64 = _AMOUNT_FLOOR, constraint::EquilibriumConstraint = FixedTP(),
        q = nothing, surface_potential::Symbol = :auto, strip_tag = nothing,
        floor::Real = _CERTIFICATE_FLOOR,
    )
    n = [ustrip(us"mol", x) for x in eq.n]
    # `local`: a closure shares the locals of the function it is defined in, and
    # the second call below, on the values, would otherwise replace the blocks
    # of the problem on duals that `apply` reads at the end.
    function problem(st)
        local blocks
        p = _build_params(st; ϵ = ϵ)
        d = _STRIP_TAGS[] == () ? des : _rebuilt(des)
        blocks = _constraint_blocks(constraint, d, st, p, n)
        if surface_potential !== :eliminated
            surf = _surface_potential_blocks(d, st, p, n)
            surf === nothing || (blocks = _compose_blocks(blocks, surf))
        end
        return _dual_problem(d, p, n, _scoped_blocks(blocks)), blocks
    end
    prob, blocks = problem(state)
    # The same problem on the values, for the Jacobian of the tangent: the duals
    # of the data, of the model or of the constraint live in closures the solver
    # cannot strip by itself.
    primal = strip_tag === nothing ? nothing :
        first(with(() -> problem(_strip_state(state, strip_tag)), _STRIP_TAGS => (_STRIP_TAGS[]..., strip_tag)))
    t = _optima_tangent(prob, b, n; q = (q === nothing || isempty(q)) ? nothing : q, primal = primal, floor = floor)
    T_out, P_out = blocks.apply(temperature(state), pressure(state), t.q)
    return ChemicalState(des.system, t.x .* u"mol"; T = T_out, P = P_out), t.q
end

"""
    _complete_floored_solutes(des, prob, res, floor) -> Vector{Float64}

The amounts of `res`, with each solute below the activity floor given the amount
the potentials of the answer give it.

The activity models read a solute through `max(n, floor)` (see
[`_ACTIVITY_FLOOR`](@ref)), so below the floor its log activity does not depend on
its amount. The dual solve recovers a solute from its own stationarity, `w ← w + r`
with `r = uᵢ − ∇fᵢ`, a step that assumes the log activity follows `ln n`: below
the floor it gains `r` per sweep, whatever the distance to its answer, and it can
be left there. The equation of the model is then solved exactly by
`n = floor·exp(r)` when `r > 0`; when `r ≤ 0` any amount below the floor
satisfies it, and the one computed is kept. The element balance moves by the
amounts added, of the order of the floor.

Measured with the floor at 1e-16, where it stood until 0.28.0, on a cement paste
at pH 14 in 7.8 g of water, whose equilibrium holds 1.2e-16 mol of H⁺: a solve
left it at 3e-100 mol, and `pH(eq, model)`, which reads the amount, came out 0.09
high. With the floor at 1e-30 such a solute is far above it, and this is the net
for one that is not.
"""
function _complete_floored_solutes(des::DualEquilibriumSolver, prob, res, floor::Real)
    des.j_solvent > 0 || return res.x
    jw = des.idx_aq[des.j_solvent]
    solutes = [i for i in des.idx_aq if i != jw && res.x[i] < floor]
    isempty(solutes) && return res.x
    return _optima_complete_floored(prob, res, solutes, floor)
end

"""
    optimality_certificate(des, state; b = nothing, ϵ = 1e-16, floor = 1e-25)
        -> (; stationarity, stationarity_floored, balance, balance_relative,
             worst_supersaturation, n_interior,
             n_absent_component, param_residual, worst_violation_split,
             split_phases, split_trials, optimal, scope, scope_reasons,
             reduced_curvature, stationarity_abs,
             ionic_strength, activity_range, within_activity_range)

Check the KKT conditions at a composition, independently of how it was obtained.

For a convex problem these conditions are sufficient, so `optimal = true` is a
proof of **global** optimality. Whether the problem is one is not always the case,
and `scope` says what the proof covers here: `:global_minimum`, `:local_minimum`,
`:kkt_point` or `:self_consistent`, with `scope_reasons` naming the property that
decided and `reduced_curvature` the smallest eigenvalue of the Hessian over the
directions that conserve matter (see `_certificate_scope`). The extended activity
models in general use (B-dot, Davies with neutral species) are not the gradient
of one Gibbs energy, and a certified equilibrium computed with them is a
composition consistent with its own activities rather than the minimum of an
energy. `stationarity_abs` is the stationarity residual in `RT` units, before
the scaling `stationarity` applies. Use it to audit any solver — including
[`EquilibriumSolver`](@ref), whose interior-point iteration reports `MaxIters` on
a cement equilibrium and cannot say whether the point it returns is the answer.

The three quantities are the stationarity of the interior species, the component
balance, and the worst saturation index among absent phases (negative when every
one of them is undersaturated, as optimality requires). The balance is reported
twice: `balance`, the worst row in moles, and `balance_relative`, the worst row
relative to what it holds. Each row is judged on the larger of the two below one
mole and in moles above it, so that a trace is held to its own amount, as
PHREEQC and GEMS hold a mass balance to its element total; judged in moles alone,
as until 0.29, a trace of 1e-9 mol could be 10 % wrong and certified. A component
whose budget is below 1e-12 of the largest is one nobody supplies, its carriers
held at the floor, and its row is judged in moles. `stationarity_floored` is
the first, one-sided, on the members of a present phase held below `floor`: such
a member may hold more than its exact amount, by truncation, but not less, and it
fails when the search left it far below what the multipliers give it.

`worst_violation_split` extends that last test to the phases that are **present**:
Michelsen's tangent-plane distance, which asks whether a mixing phase would lower
the Gibbs energy by separating into two compositions. It is `-Inf` when no phase
could be tested, negative when every one of them is stable, and positive when one
wants to unmix — `split_phases` then names them and `split_trials` carries, per
phase, the composition it wants to split into. That composition is what
[`equilibrate_split`](@ref) seeds a second instance with; it exists nowhere else,
being a property of the full system and not of the mixing model alone.

`ionic_strength` is the effective ionic strength of the composition (NaN without
an aqueous phase), `activity_range` the ionic strength up to which the activity
model is stated valid ([`activity_model_range`](@ref), `nothing` when no range is
stated), and `within_activity_range` whether the first lies within the second
(`nothing` when either is unknown). An answer can certify outside the range of
its model: the certificate says so, and the caller decides.
"""
function optimality_certificate(
        des::DualEquilibriumSolver, state::ChemicalState;
        b = nothing, ϵ::Float64 = _AMOUNT_FLOOR, floor::Float64 = _CERTIFICATE_FLOOR,
        constraint::EquilibriumConstraint = FixedTP(),
        q = nothing,
    )
    # A composition carrying dual numbers is audited at its values: the
    # conditions are about the point, and its derivatives are another question.
    state = _primal(state)
    p = _build_params(state; ϵ = ϵ)
    n = Float64[ustrip(us"mol", x) for x in state.n]
    bv = b === nothing ? des.A * n : Float64[_plain(x) for x in collect(b)]

    # The certificate has to audit the problem that was SOLVED, and a constraint
    # is part of that problem. Rebuilt with `FixedTP` — which is what this did
    # until 0.16 — the audit misses two things at once: the conservation rows
    # lose their `Aq q` term, so a prescribed activity or pH is measured against
    # a budget short by exactly the titrant amount; and `hq` is skipped, so a
    # constraint that shifts a chemical potential is measured against the
    # unshifted one and can never certify however right it is.
    blocks = _constraint_blocks(constraint, des, state, p, n)

    # A solve that carried the surface potential as an unknown returns it in `q`
    # after the constraint's own unknowns, and `equilibrate_certified` hands that
    # whole vector here. It is then the augmented problem that was solved, and the
    # one audited: the surface block is composed exactly as `solve` composes it,
    # so the closure of each potential is part of the parameter residual. Built
    # from the constraint alone, as this did until 0.23, a prescribed pH on a
    # diffuse layer met a `q` one entry too long and threw a `DimensionMismatch`
    # before measuring anything. Without the potentials in `q` — a solve that
    # eliminated them, or a caller auditing a bare composition — the eliminated
    # form is audited, which is exact: the activity model then evaluates each
    # potential from the composition itself.
    if q !== nothing && length(q) != blocks.nq
        surf = _surface_potential_blocks(des, state, p, n)
        nsurf = surf === nothing ? 0 : surf.nq
        length(q) == blocks.nq + nsurf || throw(
            ArgumentError(
                "q has $(length(q)) entries; the constraint has $(blocks.nq) " *
                    "unknowns and the system's surfaces $nsurf, so it should " *
                    "have $(blocks.nq) or $(blocks.nq + nsurf)."
            ),
        )
        blocks = _compose_blocks(blocks, surf)
    end
    qv = blocks.nq == 0 ? nothing :
        (q === nothing ? Float64[_plain(x) for x in blocks.q0] : Float64[_plain(x) for x in q])

    c = _optima_kkt_certificate(
        _dual_problem(des, p, n, _scoped_blocks(blocks)), n, bv, floor,
        des.opts.tol, des.opts.si_tol, qv,
    )
    scope, scope_reasons, reduced_curvature = _certificate_scope(
        des, p, n, constraint; Aq = blocks.Aq, floor = floor,
    )
    I, I_max, within = _activity_range_report(des.model, state)
    return (;
        # In moles, and relative to what each row holds; `optimal` judges the
        # larger of the two below one mole (OptimaSolver's `kkt_certificate`).
        stationarity = c.stationarity, balance = c.feasibility_abs,
        balance_relative = c.feasibility_rel,
        # The unscaled stationarity, in RT units. `stationarity` is divided by the
        # size of the potentials it is built from — they are of order 10²-10³, so
        # an absolute threshold on their residual would ask for thirteen digits of
        # cancellation — and this is the raw figure for reporting.
        stationarity_abs = c.stationarity_abs,
        # The same condition on the members of a present phase held below the
        # floor, one-sided: such a member may hold less than its multipliers
        # give it only by truncation, never by being left behind. Read through
        # `hasproperty` for a certificate from a back end that does not run it.
        stationarity_floored = hasproperty(c, :stationarity_floored) ? c.stationarity_floored : 0.0,
        worst_supersaturation = c.worst_violation, n_interior = c.n_interior,
        n_absent_component = c.n_forced_zero,
        # Zero on the unconstrained route, so it costs nothing there and is the
        # constraint's own residual when there is one.
        param_residual = hasproperty(c, :param_residual) ? c.param_residual : 0.0,
        # Michelsen's split verdict on the PRESENT mixing phases, forwarded
        # rather than dropped. `equilibrate_split` seeds a second instance from
        # `split_trials`, and it is the only place that composition exists: the
        # trial is a property of the full system, fixed jointly with the solution
        # the phase sits in, so a caller cannot recompute it from the mixing
        # model alone. Read through `hasproperty` because a certificate also
        # arrives from a back end that does not run the test.
        worst_violation_split = hasproperty(c, :worst_violation_split) ?
            c.worst_violation_split : -Inf,
        split_phases = hasproperty(c, :split_phases) ? c.split_phases : Int[],
        split_trials = hasproperty(c, :split_trials) ? c.split_trials :
            Dict{Int, NamedTuple{(:members, :x), Tuple{Vector{Int}, Vector{Float64}}}}(),
        optimal = c.optimal,
        # What `optimal = true` proves for this problem; see `_certificate_scope`.
        scope, scope_reasons, reduced_curvature,
        # Whether the answer lies where its activity model is stated valid. The
        # ionic strength is that of the free ions, the one the model itself
        # uses; `within_activity_range` is `nothing` when the model states no
        # range or there is no aqueous phase.
        ionic_strength = I, activity_range = I_max, within_activity_range = within,
    )
end

"""
    _activity_range_report(model, state) -> (I, I_max, within)

The effective ionic strength of `state` (NaN without an aqueous phase), the range
[`activity_model_range`](@ref) states for `model`, and whether the first lies
within the second (`nothing` when either is unknown).
"""
function _activity_range_report(model::AbstractActivityModel, state::ChemicalState)
    I_max = activity_model_range(model)
    I = isempty(state.system.idx_solvent) ? NaN : ionic_strength(state)
    within = (I_max === nothing || !isfinite(I)) ? nothing : I <= I_max
    return (I, I_max, within)
end

# Above this relative asymmetry of the Jacobian of the log activities, the
# activities are not the gradient of one Gibbs energy. The two populations sit
# far apart: an exact model measures 2e-16 (the ideal dilute model, Davies on
# ions, the Debye-Hückel form with a common ion size and no linear term), and
# the published extended forms measure 1.2e-2 for a solvent row taken as Raoult's
# mole fraction against molality solutes, 0.30 for the GEMS setting of the B-dot
# model and 1.0 for its default and for Davies with neutral species, which carry
# a salting-out term; the threshold sits eight decades above the first and six
# below the smallest of the others.
const _SCOPE_ASYMMETRY = 1.0e-8

"""
    _is_variational(constraint) -> Bool

Whether a constraint keeps the equilibrium a minimization whose optimality
conditions are sufficient. `CapillaryWater` shifts the activity of water by an
amount that is not the gradient of a convex energy for an arbitrary retention
law, and the constraints that make the temperature or the pressure an unknown
turn the solve into a KKT system whose convexity is not established here.
"""
_is_variational(::EquilibriumConstraint) = true
_is_variational(::CapillaryWater) = false
_is_variational(::FixedEnthalpy) = false
_is_variational(::Adiabatic) = false
_is_variational(::FixedVolume) = false
_is_variational(::SealedVolume) = false

# Below this smallest eigenvalue of the reduced Hessian, scaled so that an ideal
# solute contributes one, a point is not called a strict local minimum.
const _REDUCED_CURVATURE_TOL = 1.0e-8

"""
    _aqueous_convexity(model, system, p) -> (proved, reason)

Whether the Gibbs energy of the aqueous phase under `model` is proved convex over
its whole domain, every composition of the declared solutes at the temperature and
pressure of `p`, and if not, why. It is a sufficient condition, and it is what a
`:global_minimum` scope rests on; Gibbs–Duhem consistency, which the symmetry of
the Jacobian measures, is checked before and is not repeated here.

  - The ideal dilute model: its energy, `Σₛ nₛ (ln(nₛ/(n_w M_w)) − 1)` plus a
    linear term, has the Hessian `Σₛ nₛ (vₛ/nₛ − v_w/n_w)²`, positive.
  - A Debye–Hückel form with one function of `I` times `zᵢ²` on every ion, with
    the solvent row that integrates it: the HKF model with one ion size, `Ḃ = 0`
    and no salting-out, and Davies with neutral species carrying none. Its excess
    energy adds a term of rank one to the ideal Hessian, which the Cauchy–Schwarz
    inequality bounds by it as long as
    `ln(10)·A·z_max²·√I / (2(1 + κ√I)²) ≤ 1`, `κ = B å` (one for Davies, whose
    `b I` term only adds convexity), and the supremum over `I`, at `√I = 1/κ`,
    gives `ln(10)·A·z_max²/(8κ) ≤ 1`. See [What the certificate proves, and when](@ref sec-theory-certificate-scope)
    for the derivation.
  - Any other model: not proved, whether or not it is convex.
"""
_aqueous_convexity(model::AbstractActivityModel, system, p) =
    (false, "no convexity proof is available for $(nameof(typeof(model)))")
_aqueous_convexity(::DiluteSolutionModel, system, p) = (true, "")

function _ions_and_neutrals(system)
    ions = [i for i in system.idx_solutes if !iszero(charge(system.species[i]))]
    neutrals = [i for i in system.idx_solutes if iszero(charge(system.species[i]))]
    return ions, neutrals
end

# The bound of the rank-one Debye–Hückel term against the ideal Hessian.
function _debye_huckel_bound(A, κ, ions, system)
    isempty(ions) && return (true, "")
    κ > 0 || return (false, "without an ion size the Debye–Hückel term is not bounded by the ideal one")
    zmax = maximum(i -> abs(charge(system.species[i])), ions)
    bound = log(10.0) * A * zmax^2 / (8 * κ)
    bound <= 1 && return (true, "")
    return (
        false, "the convexity bound ln(10)·A·z²/(8κ) is $(round(bound; sigdigits = 3)) for " *
            "the charge $(zmax), above one",
    )
end

function _aqueous_convexity(model::HKFActivityModel, system, p)
    iszero(model.Ḃ) || return (false, "the linear term Ḃ I of the B-dot model has no convexity proof")
    ions, neutrals = _ions_and_neutrals(system)
    any(i -> !iszero(_setschenow(system.species[i], model)), neutrals) &&
        return (false, "the salting-out term of the neutral species has no convexity proof")
    å = unique(_hkf_lookup_å(system.species[i], model) for i in ions)
    length(å) <= 1 || return (false, "the ions do not share one ion size")
    AB = model.temperature_dependent ? hkf_debye_huckel_params(_plain(p.T), _plain(p.P)) : (A = model.A, B = model.B)
    return _debye_huckel_bound(_plain(AB.A), isempty(å) ? 1.0 : _plain(AB.B) * only(å), ions, system)
end

function _aqueous_convexity(model::DaviesActivityModel, system, p)
    model.b >= 0 || return (false, "a negative b of the Davies equation has no convexity proof")
    ions, neutrals = _ions_and_neutrals(system)
    !isempty(neutrals) && !iszero(model.bₙ) &&
        return (false, "the salting-out term of the neutral species has no convexity proof")
    A = model.temperature_dependent ? hkf_debye_huckel_params(_plain(p.T), _plain(p.P)).A : model.A
    return _debye_huckel_bound(_plain(A), 1.0, ions, system)
end

# Whether a site family's mixing energy is convex: ideal mixing on the sites,
# with a constant capacitance or without, whose charging energy is a convex
# quadratic of the surface charge.
_site_mixing_convex(::IdealSiteMixing) = true
_site_mixing_convex(m::ConstantCapacitance) = _site_mixing_convex(m.base)
_site_mixing_convex(::AbstractSiteMixingModel) = false

"""
    _reduced_curvature(J, A, n, Aq, floor) -> Float64

The smallest eigenvalue of the Hessian of the Gibbs energy over the directions
that conserve matter and move only the species present, `J` being the Jacobian of
the log activities at `n` (the Hessian of `G/RT`, symmetric where this is called).

The species present are those above `floor`; the directions are the null space of
their columns of the conservation matrix, with the columns `Aq` of a constraint's
titrant, along which the energy is linear. The coordinates are scaled by `√nᵢ`, so
that an ideal solute contributes exactly one whatever its amount, and the figure is
dimensionless. Positive, the point is a strict local minimum on its active set;
negative, a direction lowers the energy and the point is a saddle. `Inf` when no
direction is free.
"""
function _reduced_curvature(J, A, n, Aq, floor)
    F = findall(>(floor), n)
    isempty(F) && return Inf
    d = sqrt.(n[F])
    H = d .* ((J[F, F] .+ transpose(J[F, F])) ./ 2) .* transpose(d)
    C = A[:, F] .* transpose(d)
    nq = Aq === nothing ? 0 : size(Aq, 2)
    if nq > 0
        C = hcat(C, Aq)
        H = [H zeros(length(F), nq); zeros(nq, length(F) + nq)]
    end
    Z = nullspace(C)
    size(Z, 2) == 0 && return Inf
    return eigmin(Symmetric(transpose(Z) * H * Z))
end

"""
    _certificate_scope(des, p, n, constraint; Aq = nothing, floor = _CERTIFICATE_FLOOR)
        -> (scope, reasons, curvature)

What `optimal = true` proves for the problem at the composition `n`:

  - `:global_minimum` when the log activities are the gradient of one Gibbs
    energy, that energy is proved convex over the whole domain — the ideal
    phases, the aqueous models of [`_aqueous_convexity`](@ref), convex solid
    solutions and ideal site mixing — and the constraint is variational, so that
    the optimality conditions are sufficient;
  - `:local_minimum` when the activities are such a gradient and the constraint
    variational, convexity is not proved, and the Hessian of the energy over the
    directions that conserve matter is positive definite at the answer: a strict
    local minimum, which a non-convex model can hold beside others;
  - `:kkt_point` when that Hessian is not positive definite — a saddle if a
    direction lowers the energy, a minimum that is not strict if one leaves it
    flat — or when the constraint is not variational;
  - `:self_consistent` when the activities are not the gradient of one energy,
    which is the case of the extended activity models in general use and of a
    diffuse layer: the conditions then state a composition consistent with its
    own activities, not the minimum of an energy.

`reasons` says which property decided, in words, and `curvature` is the smallest
eigenvalue of [`_reduced_curvature`](@ref) (`NaN` where it is not computed: no
energy, or a constraint that is not variational).
"""
function _certificate_scope(des::DualEquilibriumSolver, p, n, constraint; Aq = nothing, floor::Real = _CERTIFICATE_FLOOR)
    reasons = String[]
    gradient = true
    cs = des.system
    J = ForwardDiff.jacobian(x -> des.lna(x, p), n)
    asym, (i, j) = _jacobian_asymmetry(J)
    if asym > _SCOPE_ASYMMETRY
        gradient = false
        push!(
            reasons,
            "the log activities are not the gradient of one Gibbs energy: their " *
                "Jacobian is asymmetric by $(round(asym; sigdigits = 2)) between " *
                "$(symbol(cs.species[i])) and $(symbol(cs.species[j]))",
        )
    end
    convex = true
    fams = cs.site_families
    if fams !== nothing
        for f in fams
            if !is_gradient_consistent(f.model)
                gradient = false
                push!(
                    reasons,
                    "the site family $(name(f)) carries a diffuse layer, whose potential " *
                        "depends on the ionic strength of a solution that does not depend " *
                        "on the surface in return",
                )
            elseif !_site_mixing_convex(f.model)
                convex = false
                push!(reasons, "the mixing energy of the site family $(name(f)) has no convexity proof")
            end
        end
    end
    gradient || return (:self_consistent, reasons, NaN)
    if !_is_variational(constraint)
        push!(
            reasons,
            "the constraint $(nameof(typeof(constraint))) leaves a system of " *
                "optimality conditions whose sufficiency is not established",
        )
        return (:kkt_point, reasons, NaN)
    end
    if !isempty(cs.idx_solvent)
        ok, why = _aqueous_convexity(des.model, cs, p)
        ok || (convex = false; push!(reasons, "the aqueous phase: " * why))
    end
    ss = cs.solid_solutions
    if ss !== nothing
        T = p.T
        for ph in ss
            # A compound-energy model is judged with the energies of its members
            # at the temperature of the solve.
            g = ph.model isa CompoundEnergyModel ?
                [Float64(ForwardDiff.value(p.ΔₐG⁰overRT[findfirst(s -> symbol(s) == symbol(m), cs.species)])) for m in ph.end_members] :
                nothing
            c = mixing_convexity(ph.model, length(ph.end_members); T = ForwardDiff.value(T), g)
            c.verdict === :convex && continue
            convex = false
            push!(
                reasons,
                c.verdict === :nonconvex ?
                    "the mixing energy of $(name(ph)) is concave on part of its range, so " *
                    "the present phases are tested against splitting but the minimum " *
                    "is not proved global" :
                    "the convexity of the mixing energy of $(name(ph)) could not be " *
                    "decided ($(c.how)), so the minimum is not proved global",
            )
        end
    end
    # A real gas alone is a pure phase, its energy linear in its amount whatever
    # the equation of state, the root of least energy having already chosen
    # between its vapor and its liquid. A mixture of real gases can split into
    # two fluids of different compositions, and no convexity is proved for it.
    if length(cs.idx_gas) > 1 && _gas_mixing(cs) !== nothing
        convex = false
        push!(
            reasons,
            "the gas phase is a mixture following the equation of state of Peng and Robinson, " *
                "which can split into two fluids, so the minimum is not proved global",
        )
    end
    λ = _reduced_curvature(J, des.A, n, Aq, floor)
    convex && return (:global_minimum, reasons, λ)
    if λ > _REDUCED_CURVATURE_TOL
        push!(
            reasons,
            "the Hessian of the energy over the directions that conserve matter is " *
                "positive definite (smallest eigenvalue $(round(λ; sigdigits = 3))): a strict " *
                "local minimum",
        )
        return (:local_minimum, reasons, λ)
    end
    push!(
        reasons,
        λ < -_REDUCED_CURVATURE_TOL ?
            "a direction that conserves matter lowers the energy (curvature " *
            "$(round(λ; sigdigits = 3))): the point is a saddle, not a minimum" :
            "the energy is flat along a direction that conserves matter: the minimum " *
            "is not strict",
    )
    return (:kkt_point, reasons, λ)
end

"""
    _judged_balance(cert) -> Float64

The balance as the certificate judges it: the larger of `balance` (moles) and
`balance_relative` (relative to what each row holds), which is the worst row on
`min(scale, 1)`. A certificate without the relative figure gives `balance`.
"""
_judged_balance(cert) = hasproperty(cert, :balance_relative) ?
    max(cert.balance, cert.balance_relative) : cert.balance

# The balance in words, for a message: both figures.
_balance_text(cert) = hasproperty(cert, :balance_relative) ?
    "element balance $(cert.balance) mol ($(cert.balance_relative) of the row that holds it)" :
    "element balance $(cert.balance)"

"""
    _kkt_error(cert) -> Float64

How far a composition is from satisfying the KKT conditions: the worst of the
three residuals the certificate reports, in one number.

All of them, and not the stationarity alone (with its one-sided form on the
members of a present phase below the floor). A composition can be stationary to
1e-3 while violating mass conservation by **moles** — measured, an answer with
stationarity 2.4e-3, an element balance off by 6.7 mol and a phase supersaturated
by 45 — and that is not a near-answer, it is not an answer to this problem at
all. Ranking on stationarity alone lets such a point beat a candidate that
conserves mass, which is how a multi-start search can discard the good answer it
just computed.

The constraint's own residual counts too, when there is one. Leaving it out would
rank a candidate that minimizes the Gibbs energy while violating the very
equation that makes it a *constrained* answer above one that satisfies both —
the same hole OptimaSolver records for the kinetic step, where "a march that
should have stopped at saturation dissolved everything and was proved optimal".
Read with `hasproperty` because a certificate also arrives from the kinetic
route, which builds its own.

`worst_supersaturation` is clamped at zero because a negative value is not an
error: it means every absent phase is undersaturated, as optimality requires.

The balance enters as the certificate judges it (`_judged_balance`): relative to
what each row holds below one mole, in moles above it, so that a candidate that
lost a trace ranks below one that kept it. The others are log-activities. It ranks
candidates of **one** problem; it is not a measure to compare two problems with.
"""
_kkt_error(cert) = max(
    cert.stationarity, _judged_balance(cert), max(cert.worst_supersaturation, 0.0),
    hasproperty(cert, :param_residual) ? cert.param_residual : 0.0,
    hasproperty(cert, :stationarity_floored) ? cert.stationarity_floored : 0.0,
)

"""
    solve_certified(des, starts; b = nothing, ϵ = 1e-16, floor = 1e-25, memo = nothing)
        -> (state, certificate)

Solve from each starting composition in `starts` and return the first answer
[`optimality_certificate`](@ref) **proves** optimal. If none is proved, return the
one with the smallest KKT error, together with its certificate, so the caller sees
what it is getting.

# Why trying more than one start is the rigorous thing to do, not a fudge

The certificate is an **oracle**: for a convex problem it does not rank answers, it
decides them. Given an oracle, running several routes and keeping a proved answer
is not guesswork — the proof is the same proof whichever route produced the point,
and nothing about it depends on having predicted the winner. What would be a fudge
is picking a route by taste and reporting its output unproved.

It is also necessary, because no single route dominates. Measured on an LC³
equilibrium at three degrees of reaction, with the dual Newton started from each
interior-point back end in turn:

| degree of reaction | from `OptimaOptimizer` | from `IpoptOptimizer` |
|---|---|---|
| 0.05 | certified, 1.8e-12 | certified, 9.1e-13 |
| 0.25 | **not certified**, 7.2 | certified, 4.2e-12 |
| 1.00 | certified, 1.2e-11 | **not certified**, 3.5e-3 |

Each back end solves what the other misses. With both offered, all three are
proved.

# Example

```julia
using Optimization, OptimizationIpopt      # for IpoptOptimizer

starts = [
    SciMLBase.solve(EquilibriumSolver(cs, model, OptimaOptimizer()), st; b = b),
    SciMLBase.solve(EquilibriumSolver(cs, model, IpoptOptimizer()), st; b = b),
]
eq, cert = solve_certified(des, starts; b = b)
cert.optimal || @warn "no start produced a certifiable answer" cert
```

The starts are supplied by the caller rather than built here, so this adds no
dependency: whichever back ends are loaded are the ones available.

# Offering the same start twice

`memo`, when given, is an `IdDict` that records the answer and the certificate
reached from each start, keyed on the start object itself, and a start met again
is answered from it instead of being solved a second time. The dual solve and
the certificate are deterministic functions of the start once `des`, `b`, `ϵ`,
`floor` and `constraint` are fixed, so the record is the result, bit for bit;
the caller owns the dictionary and must not share it between calls that differ
in any of those. [`equilibrate_certified`](@ref) keeps one for the duration of a
call, because its search offers the same cached starts again after the
continuation, after each restart and in each repair round: measured on a CEM I
paste that does not certify at once, twelve of the thirty-six dual solves of one
call were repeats of earlier ones.
"""
function solve_certified(
        des::DualEquilibriumSolver, starts;
        b = nothing, ϵ::Float64 = _AMOUNT_FLOOR, floor::Float64 = _CERTIFICATE_FLOOR,
        constraint::EquilibriumConstraint = FixedTP(),
        parameters::Union{Nothing, Base.RefValue} = nothing,
        memo::Union{Nothing, IdDict} = nothing,
        stop::Union{Nothing, Function} = nothing,
    )
    # A budget carrying dual numbers, or starts whose temperature, pressure or
    # amounts do: the search runs on the values, one level of duals down, and
    # only the answer it keeps is lifted to the caller's duals, at the
    # conditions the starts carry. A lazy sequence of starts is never inspected
    # here (the certified route lifts its own duals before calling this).
    explicit = starts isa Union{Tuple, AbstractVector} && !isempty(starts)
    D = explicit ?
        _input_number_type(first(starts), b; constraint = constraint, captured = des.captured) :
        promote_type(Float64, _carries_duals(b) ? mapreduce(typeof, promote_type, (x for x in b if x isa ForwardDiff.Dual)) : Float64)
    if D <: ForwardDiff.Dual
        Tg = ForwardDiff.tagtype(D)
        s1 = explicit ? first(starts) : nothing
        qv = Ref{Any}(Float64[])
        eq, cert = with(_STRIP_TAGS => (_STRIP_TAGS[]..., Tg)) do
            solve_certified(
                _rebuilt(des), explicit ? map(s -> _strip_state(s, Tg), starts) : starts;
                b = b === nothing ? nothing : _strip_tag(collect(b), Tg), ϵ = ϵ, floor = floor,
                constraint = constraint, parameters = qv, memo = memo, stop = stop,
            )
        end
        eq === nothing && return (eq, cert)
        bd = b !== nothing ? collect(b) : des.A * _build_n0(s1)
        at = explicit ?
            ChemicalState(des.system, eq.n; T = temperature(s1), P = pressure(s1)) : eq
        eq_d, q_d = _lift_equilibrium(des, at, eq, bd; ϵ = ϵ, constraint = constraint, q = qv[], strip_tag = Tg)
        parameters === nothing || (parameters[] = q_d)
        return (eq_d, cert)
    end
    best = nothing
    best_cert = nothing
    best_err = Inf
    best_q = Float64[]
    best_ok = false
    for s0 in starts
        # Each candidate's own parameters, captured here rather than written
        # straight into the caller's `Ref`. Written straight through, the `Ref`
        # would end up holding the LAST candidate's parameters while the state
        # returned is the best-by-error one — so a prescribed-pH scan would
        # report a titrant amount belonging to a different composition.
        hit = memo === nothing ? nothing : get(memo, s0, nothing)
        if hit === nothing
            qref = Ref(Float64[])
            eq = SciMLBase.solve(
                des, s0; b = b, ϵ = ϵ, constraint = constraint, parameters = qref,
            )
            # The certificate is evaluated at the T and P the constrained solve
            # FOUND, which `eq` carries — not at the ones the start had. `∇f`
            # depends on both, so certifying against the start's conditions would
            # measure the stationarity of a different problem. `q` is part of what
            # the solve found in exactly the same way, and is passed for the same
            # reason.
            cert = optimality_certificate(
                des, eq; b = b, ϵ = ϵ, floor = floor,
                constraint = constraint, q = qref[],
            )
            memo === nothing || (memo[s0] = (eq, cert, qref[]))
        else
            eq, cert, q = hit
            qref = Ref(q)
        end
        if cert.optimal
            parameters === nothing || (parameters[] = qref[])
            return (eq, cert)
        end
        # A caller may know that an uncertified answer is already what it needs:
        # a phase declared `instances = :auto` asking to split, which no other
        # start can make certify with one composition.
        if stop !== nothing && stop(cert)
            parameters === nothing || (parameters[] = qref[])
            return (eq, cert)
        end
        # Ranked on the KKT error, but only among answers the model can describe:
        # a composition whose solvent has been taken by the solids is outside the
        # formulation, and letting it win on a smaller residual hides every
        # candidate that conserved mass behind it. See `_within_domain`.
        err = _kkt_error(cert)
        ok = _within_domain(eq)
        if (ok && !best_ok) || ((ok == best_ok) && err < best_err)
            best_err = err
            best = eq
            best_cert = cert
            best_q = qref[]
            best_ok = ok
        end
    end
    parameters === nothing || (parameters[] = best_q)
    return (best, best_cert)
end

# ── hooks filled in by the OptimaSolver extension ────────────────────────────
#
# The algorithm is `OptimaSolver`'s, and that package is a weak dependency here,
# so the four entry points are indirected through functions the extension
# overrides. Called without it loaded, they say what is missing.

_optima_dual_problem(args...) = _need_optima()
_optima_dual_solve(args...) = _need_optima()
_optima_kkt_certificate(args...) = _need_optima()
_optima_lp(args...) = _need_optima()
_optima_complete_floored(args...) = _need_optima()
_optima_tangent(args...; kwargs...) = _need_optima()

_need_optima() = error(
    "DualEquilibriumSolver needs OptimaSolver ≥ 0.7: the KKT solver and its " *
        "certificate live there. Add `using OptimaSolver`."
)
