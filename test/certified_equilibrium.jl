# The certified equilibrium path.
#
# Every assertion here is on a measured property of the solvers, not on a target
# value: that the certificate decides, that it is not fooled by a start which is
# itself infeasible, and that the multi-start route certifies cases no single
# back end does.

include("reference_species.jl")

@testsection "Certified equilibrium" begin

    sp = Dict(
        symbol(s) => s for s in build_species(
                datapath("slop98-inorganic-thermofun.json")
            )
    )
    species = [sp[s] for s in split("H2O@ H+ OH- CO2@ HCO3- CO3-2 Ca+2 Cal")]
    cs = ChemicalSystem(species, ["H2O@", "H+", "Ca+2", "CO3-2", "Zz"])
    A = Matrix{Float64}(cs.SM.A)

    function calcite(; θ = 25.0, nco2 = 0.0)
        st = ChemicalState(cs)
        set_quantity!(st, "Cal", 1.0e-3u"mol")
        set_quantity!(st, "H2O@", 1.0u"kg")
        nco2 > 0 && set_quantity!(st, "CO2@", nco2 * u"mol")
        V = volume(st)
        set_quantity!(st, "H+", 1.0e-4u"mol/L" * V.liquid)
        set_quantity!(st, "OH-", 1.0e-10u"mol/L" * V.liquid)
        set_temperature!(st, (273.15 + θ) * u"K")
        return st
    end

    # Element balance, each row against its own budget. An absolute ∞-norm hides
    # the small rows: at ‖An − b‖∞ = 3e-6 the water row (111 mol) is satisfied to
    # 1.8e-8 while the charge row (budget 1e-4) is out by 3 %.
    function balance_rel(st0, eq)
        b = A * [ustrip(us"mol", x) for x in st0.n]
        r = A * [ustrip(us"mol", x) for x in eq.n] - b
        keep = abs.(b) .> 1.0e-8
        return any(keep) ? maximum(abs.(r[keep] ./ b[keep])) : maximum(abs, r)
    end

    @testsection "the certificate says whether the answer lies in its model's range" begin
        eq, cert = equilibrate_certified(calcite(); model = DaviesActivityModel())
        @test cert.optimal
        @test cert.ionic_strength ≈ ionic_strength(eq)
        @test cert.activity_range == 0.5
        @test cert.within_activity_range === true
        _, cert0 = equilibrate_certified(calcite())
        @test cert0.activity_range === nothing
        @test cert0.within_activity_range === nothing

        # A brine past the range of the Davies equation: certified, and flagged.
        brine_cs = ChemicalSystem([sp[s] for s in split("H2O@ H+ OH- Na+ Cl-")], ["H2O@", "H+", "Na+", "Cl-", "Zz"])
        brine = ChemicalState(brine_cs)
        set_quantity!(brine, "H2O@", 1.0u"kg")
        set_quantity!(brine, "Na+", 2.0u"mol")
        set_quantity!(brine, "Cl-", 2.0u"mol")
        set_quantity!(brine, "H+", 1.0e-7u"mol")
        set_quantity!(brine, "OH-", 1.0e-7u"mol")
        eqb, certb = equilibrate_certified(brine; model = DaviesActivityModel())
        @test certb.optimal
        @test certb.ionic_strength > 1.5
        @test certb.within_activity_range === false

        # With a fallback asked for out of range, the fallback's answer is
        # returned, said once, with the first certificate kept.
        eqf, certf = @test_logs (:warn, r"fallback activity model HKFActivityModel") equilibrate_certified(
            brine; model = DaviesActivityModel(), fallback_model = HKFActivityModel(), fallback_on = :out_of_range,
        )
        @test certf.optimal
        @test certf.fallback_used
        @test certf.fallback_reason === :out_of_range
        @test certf.primary_certificate.within_activity_range === false
        @test certf.activity_range == 1.0
        # Without `:out_of_range`, a certified answer is kept whatever its range.
        _, certk = equilibrate_certified(brine; model = DaviesActivityModel(), fallback_model = HKFActivityModel())
        @test !certk.fallback_used
        @test_throws ArgumentError equilibrate_certified(brine; fallback_model = HKFActivityModel(), fallback_on = :always)

        # When neither model certifies, the first answer is returned and said so,
        # and under the strict flag it raises.
        impossible = A * [ustrip(us"mol", x) for x in calcite().n]
        impossible[3] = -1.0                      # a negative calcium budget
        strict = ChemistryLab.STRICT_CONVERGENCE[]
        try
            ChemistryLab.STRICT_CONVERGENCE[] = false
            _, certn = @test_logs (:warn, r"did not certify either") match_mode = :any equilibrate_certified(
                calcite(); model = DaviesActivityModel(), fallback_model = HKFActivityModel(), b = impossible,
            )
            @test !certn.optimal
            @test !certn.fallback_used
            @test certn.fallback_reason === :refusal
            ChemistryLab.STRICT_CONVERGENCE[] = true
            @test_throws ErrorException equilibrate_certified(
                calcite(); model = DaviesActivityModel(), fallback_model = HKFActivityModel(), b = impossible,
            )
        finally
            ChemistryLab.STRICT_CONVERGENCE[] = strict
        end
    end

    @testsection "a temperature given to the solve is refused, not ignored" begin
        # It belongs to the state. Forwarded to the optimizer, as it used to be,
        # it was dropped there and the solve ran at the state's temperature.
        st = calcite()
        @test_throws ArgumentError equilibrate_certified(st; T = 293.15u"K")
        @test_throws ArgumentError equilibrate_certified(st; P = 2.0e5u"Pa")
        @test_throws ArgumentError equilibrate(st; temperature = 293.15)
        @test_throws ArgumentError equilibrate_path(st, [A * ustrip.(us"mol", st.n)]; T = 300.0)
        @test_throws ArgumentError EquilibriumSolver(cs, DiluteSolutionModel(), nothing; T = 300.0)
        err = try
            equilibrate_certified(st; T = 293.15u"K")
        catch e
            e
        end
        @test occursin("set_temperature!", err.msg)
        # What is meant is honored: the state's temperature moves the answer.
        @test pH(first(equilibrate_certified(calcite(; θ = 20.0)))) !=
            pH(first(equilibrate_certified(calcite(; θ = 25.0))))
    end

    @testsection "the certificate decides, and the answer is reproducible" begin
        st = calcite()
        eq1, c1 = equilibrate_certified(calcite())
        eq2, c2 = equilibrate_certified(calcite())
        @test c1.optimal
        @test pH(eq1) == pH(eq2)          # bit-for-bit, no path dependence
        @test balance_rel(st, eq1) < 1.0e-10
        @test c1.worst_supersaturation < 0 || c1.worst_supersaturation == -Inf
    end

    @testsection "the interior point alone is wrong on this case" begin
        # Not a target value: the point is that the certified route is two orders
        # of magnitude better on a balance the old default reported as fine.
        st = calcite()
        eq_bar = equilibrate(calcite(), OptimaOptimizer())
        eq_cert, cert = equilibrate_certified(calcite())
        @test balance_rel(st, eq_bar) > 1.0e-3      # measured at 3.0e-2
        @test balance_rel(st, eq_cert) < 1.0e-8
        @test cert.optimal
        # And they disagree on the answer, not merely on the residual.
        @test abs(pH(eq_bar) - pH(eq_cert)) > 1.0e-3
    end

    @testsection "b is fixed once, not taken from each start" begin
        # A start that violates the balance shifts the component totals by its own
        # infeasibility. Recomputing `b` from each start therefore poses a
        # different problem per start, and the dual solve certifies the answer to
        # the shifted one: measured, two "certified" compositions 0.2 % apart on
        # dissolved calcium, which a convex problem with one minimum cannot have.
        st = calcite()
        b0 = A * [ustrip(us"mol", x) for x in st.n]
        des = DualEquilibriumSolver(cs, DiluteSolutionModel())
        bar = equilibrate(calcite(), OptimaOptimizer())
        from_raw = SciMLBase.solve(des, calcite(); b = b0)
        from_bar = SciMLBase.solve(des, bar; b = b0)
        # With the same b both routes must land on the same minimum.
        @test isapprox(pH(from_raw), pH(from_bar); atol = 1.0e-6)
        # Without it, the second solves a different problem — the shift being
        # exactly the interior point's own infeasibility.
        shifted = SciMLBase.solve(des, bar)
        @test abs(pH(shifted) - pH(from_raw)) > 1.0e-3
    end

    @testsection "equilibrate certifies by default, and can be told not to" begin
        st = calcite()
        @test balance_rel(st, equilibrate(calcite())) < 1.0e-8
        @test balance_rel(st, equilibrate(calcite(); certify = false)) > 1.0e-3
    end

    @testsection "the dual route needs an aqueous phase and H2O@" begin
        # It parameterizes the interior variables by the solvent's potential, so
        # both conditions are structural. `equilibrate_certified` returns
        # `nothing` for the certificate on such a system rather than claiming a
        # proof it cannot give.
        @test ChemistryLab._dual_applicable(cs)

        solid_only = ChemicalSystem(
            [Species("NaCl"; aggregate_state = AS_CRYSTAL, class = SC_COMPONENT)],
            ["NaCl"],
        )
        @test !ChemistryLab._dual_applicable(solid_only)

        # Aqueous, but the solvent is not called `H2O@`.
        no_solvent = ChemicalSystem([sp["Ca+2"], sp["CO3-2"]], ["Ca+2", "CO3-2"])
        @test !ChemistryLab._dual_applicable(no_solvent)
    end


    @testsection "STRICT_CONVERGENCE is honored by the certified route" begin
        # The flag existed only on the interior-point retcode, so a caller who set
        # it — asking that a non-converged solve never pass as a result — still got
        # an uncertified answer back with a `@warn`. That answer can violate the
        # element balance by moles and still look like an ordinary `ChemicalState`:
        # measured on a cement paste, a balance off by 6.7 mol, every hydrate at
        # zero, and a table of amounts that reads as a result.
        #
        # `-b` is infeasible by construction: every component budget is negative and
        # every stoichiometric coefficient is non-negative, so no composition with
        # `n >= 0` can meet it and no route can certify.
        st = calcite()
        b_bad = -(A * [ustrip(us"mol", x) for x in st.n])

        strict = ChemistryLab.STRICT_CONVERGENCE[]
        try
            ChemistryLab.STRICT_CONVERGENCE[] = true
            @test_throws "no route produced a certifiable equilibrium" equilibrate_certified(
                st; b = b_bad, autostart = false
            )

            # The default stays a warning: the answer is still the best one found,
            # and `optimality_certificate` is there to audit it.
            ChemistryLab.STRICT_CONVERGENCE[] = false
            eq, cert = equilibrate_certified(st; b = b_bad, autostart = false)
            @test cert.optimal == false
            @test eq isa ChemicalState
        finally
            ChemistryLab.STRICT_CONVERGENCE[] = strict
        end

        # And the flag is left exactly as it was found.
        @test ChemistryLab.STRICT_CONVERGENCE[] == strict

    end
end

@testsection "ForwardDiff through the certified route" begin

    # The certified route solves in real arithmetic — an active set has a discrete
    # component, so the map `b ↦ n*(b)` is smooth only piecewise and pushing duals
    # through the search would be wrong. The derivative is attached at the answer,
    # from the optimality conditions.
    #
    # This has to be tested rather than assumed: the certified path converts the
    # component totals to `Float64` to fix them once, and before the dual branch
    # existed that silently dropped every partial.
    sp = Dict(
        symbol(s) => s for s in build_species(
                datapath("slop98-inorganic-thermofun.json")
            )
    )
    cs = ChemicalSystem(
        [sp[s] for s in split("H2O@ H+ OH- CO2@ HCO3- CO3-2 Ca+2 Cal")],
        ["H2O@", "H+", "Ca+2", "CO3-2", "Zz"],
    )
    idx = Dict(symbol(s) => i for (i, s) in enumerate(cs.species))

    function pH_of(x)
        n = Any[fill(0.0 * oneunit(x) * u"mol", length(cs.species))...]
        n[idx["H2O@"]] = 55.5 * oneunit(x) * u"mol"
        n[idx["Cal"]] = 0.05 * oneunit(x) * u"mol"
        n[idx["CO2@"]] = x * u"mol"
        return pH(equilibrate(ChemicalState(cs, n)))
    end

    x0 = 0.01
    ad = ForwardDiff.derivative(pH_of, x0)
    h = 1.0e-6
    fd = (pH_of(x0 + h) - pH_of(x0 - h)) / (2h)

    @test isfinite(ad)
    @test ad != 0                       # the partials are not dropped
    @test ad ≈ fd rtol = 1.0e-5         # measured at 4.3e-9
    @test ad < 0                        # more CO2, lower pH

    # And the primal value is the certified one, not something the dual path
    # computed separately.
    @test pH_of(x0) ≈ pH(
        first(
            equilibrate_certified(
                let
                    n = Any[fill(0.0u"mol", length(cs.species))...]
                    n[idx["H2O@"]] = 55.5u"mol"; n[idx["Cal"]] = 0.05u"mol"
                    n[idx["CO2@"]] = x0 * u"mol"
                    ChemicalState(cs, n)
                end
            )
        )
    ) atol = 1.0e-12

end

@testsection "the certified route restarts from its own answer" begin
    # `_keep_better` is the arbiter of every round of the multi-start search:
    # certified beats uncertified, and among uncertified the smaller KKT error
    # wins. It is unit-testable without a solve, and it is what keeps a restart
    # from ever making the answer worse.
    kb = ChemistryLab._keep_better
    a, b = :A, :B
    cert(opt, stat; bal = 0.0, si = -1.0) =
        (;
        optimal = opt, stationarity = stat, balance = bal,
        worst_supersaturation = si,
    )

    # Certified beats uncertified, in both directions and whatever the errors.
    @test first(kb(a, cert(false, 1.0e-16), b, cert(true, 1.0e-3))) === b
    @test first(kb(a, cert(true, 1.0e-3), b, cert(false, 1.0e-16))) === a

    # Among uncertified, the smaller KKT error wins.
    @test first(kb(a, cert(false, 1.0e-6), b, cert(false, 1.0e-9))) === b
    @test first(kb(a, cert(false, 1.0e-9), b, cert(false, 1.0e-6))) === a

    # And that error is the worst of ALL THREE residuals, not the stationarity
    # alone. This is the case that was wrong, and it is not academic: a Windows
    # run of a CEM I paste came back stationary to 2.4e-3 with an element
    # balance off by 6.7 mol and a phase supersaturated by 45, and it beat every
    # candidate the continuation produced — those being stationary to only 1e-2
    # while conserving mass. Ranked on stationarity, the answer that is not an
    # answer wins; ranked on the worst residual, the usable one does. It is also
    # the ranking `solve_certified` already used internally, so the two agree.
    bad = cert(false, 2.4e-3; bal = 6.7, si = 45.3)
    good = cert(false, 1.0e-2; bal = 1.0e-13, si = -0.5)
    @test ChemistryLab._kkt_error(bad) > ChemistryLab._kkt_error(good)
    @test first(kb(a, bad, b, good)) === b
    @test first(kb(a, good, b, bad)) === a

    # A negative worst supersaturation is not an error: every absent phase
    # undersaturated is what optimality requires, so it must not be counted.
    @test ChemistryLab._kkt_error(cert(false, 1.0e-9; si = -12.0)) == 1.0e-9

    # Between two certified answers the KKT error still decides, and a tie keeps
    # the incumbent: a round that buys nothing changes nothing, which is what
    # makes the restart loop safe to run and what stops it.
    @test first(kb(a, cert(true, 1.0e-12), b, cert(true, 1.0e-16))) === b
    @test first(kb(a, cert(false, 1.0e-9), b, cert(false, 1.0e-9))) === a
    @test first(kb(a, cert(true, 1.0e-16), b, cert(true, 1.0e-16))) === a

    # The bound exists so a case improving by a hair every round cannot loop.
    @test ChemistryLab._MAX_RESTARTS isa Integer
    @test ChemistryLab._MAX_RESTARTS >= 1

end

@testsection "a candidate's diagnostics are not the answer's" begin
    # `equilibrate_certified` runs every back end from several compositions
    # precisely because none of them works on every problem, and keeps whichever
    # answer the certificate proves. A candidate that does not converge is
    # therefore ordinary. Left unguarded it printed "returned `MaxIters`" and
    # "did not certify optimality" from candidates along the way, so a call that
    # ended `optimal = true` read as a failed solve.

    # The scope sets the flag and restores it, including when the body throws —
    # a start search that raises must not leave the whole session quiet.
    @test ChemistryLab._EXPLORING_STARTS[] == false
    inside = ChemistryLab._exploring_starts() do
        ChemistryLab._EXPLORING_STARTS[]
    end
    @test inside
    @test ChemistryLab._EXPLORING_STARTS[] == false
    @test_throws ErrorException ChemistryLab._exploring_starts() do
        error("a back end failed")
    end
    @test ChemistryLab._EXPLORING_STARTS[] == false

    # Nested scopes restore the previous value, not `false`: the walk runs inside
    # the search, and the inner scope ending must not un-quiet the outer one.
    ChemistryLab._exploring_starts() do
        ChemistryLab._exploring_starts() do
        end
        @test ChemistryLab._EXPLORING_STARTS[]
    end
    @test ChemistryLab._EXPLORING_STARTS[] == false

    # TWO SEARCHES THAT OVERLAP, which is what a save-and-restore around a module
    # global cannot do and a scoped value can. The ordering that breaks it is A
    # entering first and leaving first, with B still inside its own scope:
    # save-and-restore then took the flag away from B, and -- worse -- B's own
    # restore put back the value it had saved, which was A's `true`, leaving the
    # flag SET for the rest of the session. Every later solve then swallowed its
    # non-convergence warning. Measured before the fix: `false` inside B, `true`
    # afterwards; both assertions below failed.
    function interleaved(flag, wrap)
        ch1, ch2 = Channel{Nothing}(1), Channel{Nothing}(1)
        inside_b = Ref(false)
        @sync begin
            @async wrap() do                     # A: in first, out first
                put!(ch1, nothing)
                take!(ch2)
                nothing
            end
            @async begin                         # B: overlaps A
                take!(ch1)
                wrap() do
                    put!(ch2, nothing)
                    yield()
                    sleep(0.05)
                    inside_b[] = flag()
                    nothing
                end
            end
        end
        return inside_b[], flag()
    end

    inside_b, afterwards = interleaved(
        () -> ChemistryLab._EXPLORING_STARTS[], ChemistryLab._exploring_starts
    )
    @test inside_b                  # B keeps its own scope
    @test afterwards == false       # and nothing is left behind

    # Same for the strict-convergence suspension, where a leak is worse still: it
    # decides whether a solve RAISES, not merely whether it warns.
    inside_b, afterwards = interleaved(
        ChemistryLab._strict_convergence, ChemistryLab._relaxed_convergence
    )
    @test inside_b == false
    @test afterwards == false

    # The suspension is internal and must not touch the caller's setting, which
    # is the documented `Ref` and stays one.
    ChemistryLab.STRICT_CONVERGENCE[] = true
    try
        @test ChemistryLab._strict_convergence()
        ChemistryLab._relaxed_convergence() do
            @test ChemistryLab._strict_convergence() == false
            @test ChemistryLab.STRICT_CONVERGENCE[]          # untouched
        end
        @test ChemistryLab._strict_convergence()
        # It is restored on the way out of a body that throws, too.
        @test_throws ErrorException ChemistryLab._relaxed_convergence() do
            error("a back end failed")
        end
        @test ChemistryLab._strict_convergence()
    finally
        ChemistryLab.STRICT_CONVERGENCE[] = false
    end

    # The counter survives as an integer under `[]`, which is its documented use.
    @test ChemistryLab.NONCONVERGED[] isa Int
    let before = ChemistryLab.NONCONVERGED[]
        ChemistryLab.NONCONVERGED[] = 7
        @test ChemistryLab.NONCONVERGED[] == 7
        Threads.atomic_add!(ChemistryLab.NONCONVERGED, 1)
        @test ChemistryLab.NONCONVERGED[] == 8
        ChemistryLab.NONCONVERGED[] = before
    end

    # The contract as a caller sees it: a call that ends with a certificate emits
    # nothing at all. `@test_logs` installs a fresh logger, so the `maxlog = 1`
    # carried by those warnings does not make this depend on what ran before.
    #
    # This one assertion also guards the other silence in this release: SciMLBase
    # warns "arrays or dicts to store parameters of different types can hurt
    # performance" the moment a conservation matrix with an abstract element type
    # reaches the problem's parameters, so a regression in
    # `_concrete_conservation` shows up here as a stray record.
    # Built here rather than with `calcite()`, which is a local of another
    # testset in this file.
    sp3 = Dict(
        symbol(s) => s for s in build_species(
                datapath("slop98-inorganic-thermofun.json")
            )
    )
    cs3 = ChemicalSystem(
        [sp3[s] for s in split("H2O@ H+ OH- CO2@ HCO3- CO3-2 Ca+2 Cal")],
        ["H2O@", "H+", "Ca+2", "CO3-2", "Zz"],
    )
    st3 = ChemicalState(cs3)
    set_quantity!(st3, "H2O@", 1.0u"kg")
    set_quantity!(st3, "Cal", 1.0e-3u"mol")
    eq3, cert3 = @test_logs equilibrate_certified(st3)
    @test cert3.optimal

end

@testsection "the conservation matrix reaching the solver is concretely typed" begin
    # `ChemicalSystem` stores its stoichiometry as `Matrix{Real}` whenever
    # integer and rational coefficients coexist — a cement's does, through
    # `C3AFS0.84H4.32` and its kind — which is right for the chemistry and wrong
    # for the solver. An abstract element type boxes every entry and makes
    # `mul!(res, A, x)` the generic fallback with a dispatch per element, on a
    # product evaluated at every objective and constraint call. SciMLBase warns
    # about exactly this as soon as such an array reaches a problem's parameters,
    # and the warning was correct.
    abstract_A = Matrix{Real}([1 0 2 // 5; 0 1 3])
    @test !isconcretetype(eltype(abstract_A))
    narrowed = ChemistryLab._concrete_float(abstract_A)
    @test isconcretetype(eltype(narrowed))
    @test eltype(narrowed) <: AbstractFloat      # not Rational: the pipeline is float
    @test !ChemistryLab.SciMLBase.should_warn_paramtype(narrowed)

    # Lossy for a non-dyadic rational, and deliberately so: `2//5` becomes `0.4`,
    # which is not equal to it. It is the same rounding the rest of the solver
    # already applies, and the exact matrix stays in `system.SM.A`.
    @test narrowed ≈ abstract_A
    @test narrowed[1, 3] == 0.4
    @test narrowed[1, 3] != 2 // 5

    # A concrete array is returned untouched, identically — the narrowing is for
    # the abstract case and nothing else. An exact rational stoichiometry stays
    # exact, and a caller differentiating through `A` keeps their number type.
    for m in (Rational{Int}[1 0; 0 1], Int[1 2; 3 4], Float32[1 0; 0 1])
        @test ChemistryLab._concrete_float(m) === m
    end

    # And what the solver actually receives. The default `b = A * u0` inherits
    # the abstract element type from `A`, so it has to be narrowed too — this is
    # the assertion that caught it.
    ep = EquilibriumProblem(abstract_A, (n, q) -> n, [1.0, 1.0, 1.0])
    @test isconcretetype(eltype(ep.A))
    @test isconcretetype(eltype(ep.b))
    @test ep.A ≈ abstract_A
    @test !ChemistryLab.SciMLBase.should_warn_paramtype(
        (; A = ep.A, b = ep.b, T = 298.15)
    )

end

@testsection "the certificate names the missing phase, and the route puts it in" begin
    # A positive worst supersaturation means a phase sits at the lower bound
    # while the solution is supersaturated with respect to it. On a convex
    # problem that is a genuine KKT failure: the active set is wrong and the
    # answer is not the answer. The certificate reports it as one number;
    # `saturation_indices` says which phase, and `_repair_start` puts it in.
    #
    # The case this exists for is a phase SWAP, which an active-set loop that
    # admits one phase at a time cannot perform. Reproduced exactly on a CEM I
    # paste by solving it with `hydrotalcite` out of the phase list: all
    # 0.02515 mol of magnesium goes to brucite, the aqueous phase equilibrates
    # with that assemblage, and putting hydrotalcite back leaves it absent and
    # supersaturated by 5.58 log units while the solve is otherwise impeccable —
    # stationarity 5.9e-16, element balance 1.6e-14. That is the reported
    # failure, and `equilibrate_certified` from that state now certifies at
    # 74.1899 cm3 with the magnesium back where it belongs. Too heavy for this
    # suite; what is asserted here is the mechanism, on a system of eight
    # species.
    sp4 = Dict(
        symbol(s) => s for s in build_species(
                datapath("slop98-inorganic-thermofun.json")
            )
    )
    cs4 = ChemicalSystem(
        [sp4[s] for s in split("H2O@ H+ OH- CO2@ HCO3- CO3-2 Ca+2 Cal")],
        ["H2O@", "H+", "Ca+2", "CO3-2", "Zz"],
    )
    A4 = Float64.(cs4.SM.A)
    model4 = DiluteSolutionModel()

    # At a certified equilibrium every phase present is at LogSI = 0 and no
    # absent one is supersaturated. That is the built-in check on the whole
    # computation: if the present phases are not at zero, nothing else in the
    # result means anything.
    st4 = ChemicalState(cs4)
    set_quantity!(st4, "H2O@", 1.0u"kg")
    set_quantity!(st4, "Cal", 1.0e-3u"mol")
    eq4, cert4 = equilibrate_certified(st4; model = model4)
    @test cert4.optimal
    si4 = saturation_indices(eq4, model4)
    @test si4 isa AbstractDict
    @test length(si4) == length(cs4.species)
    @test abs(si4["Cal"]) < 1.0e-8            # present, hence saturated
    # Nothing to repair at an answer that certifies, and the route must say so
    # rather than invent a start.
    @test ChemistryLab._repair_start(
        eq4, model4, A4 * ustrip.(us"mol", eq4.n), 1.0e-16
    ) === nothing

    # Now a state where calcite is absent and the solution is grossly
    # supersaturated with respect to it — the shape of the reported failure,
    # without its size.
    st5 = ChemicalState(cs4)
    set_quantity!(st5, "H2O@", 1.0u"kg")
    set_quantity!(st5, "Ca+2", 0.1u"mol")
    set_quantity!(st5, "CO3-2", 0.1u"mol")
    si5 = saturation_indices(st5, model4)
    @test si5["Cal"] > 1.0
    b5 = A4 * ustrip.(us"mol", st5.n)

    fixed = ChemistryLab._repair_start(st5, model4, b5, 1.0e-16)
    @test fixed isa ChemicalState
    n5 = ustrip.(us"mol", fixed.n)
    i_cal = findfirst(s -> symbol(s) == "Cal", cs4.species)
    # The amount is what the recipe could make of it, scaled: a chemical bound,
    # not a guess at the answer. Calcite takes one Ca and one CO3, and there are
    # 0.1 mol of each.
    @test n5[i_cal] ≈ ChemistryLab._REPAIR_FRACTION * 0.1 rtol = 1.0e-8
    @test n5[i_cal] > 0
    # Everything else is left alone: the element balance of a start is not this
    # function's business, since `b` is fixed by the caller and every start is
    # projected onto it.
    @test all(
        n5[i] == ustrip(us"mol", st5.n[i]) for i in eachindex(n5) if i != i_cal
    )
    # And the temperature and pressure of the state it came from are carried.
    @test temperature(fixed) == temperature(st5)
    @test pressure(fixed) == pressure(st5)

    # A round of the repair, with the search injected. The situation it exists
    # for -- a back end converging onto the wrong active set -- needs a system of
    # some 135 species to arise, while the round's logic needs eight, so the
    # search is passed in rather than closed over.
    solve_from(f) = equilibrate_certified(f; model = model4, b = b5, autostart = false)
    # The certificate the round is handed stands for "the back ends failed",
    # which is what it is only ever called after.
    failed = (;
        optimal = false, stationarity = 1.0, balance = 1.0,
        worst_supersaturation = 10.0,
    )
    eq6, cert6, improved6 = ChemistryLab._repair_round(
        st5, failed, model4, b5, 1.0e-16, solve_from, false,
    )
    @test improved6                          # anything beats that certificate
    @test cert6.optimal                      # and here the repair certifies
    @test ustrip(us"mol", eq6.n[i_cal]) > 1.0e-6    # calcite precipitated
    @test abs(saturation_indices(eq6, model4)["Cal"]) < 1.0e-8   # and is saturated

    # Nothing to repair: the round says so and changes nothing, which is what
    # stops the loop.
    eq7, cert7, improved7 = ChemistryLab._repair_round(
        eq4, failed, model4, A4 * ustrip.(us"mol", eq4.n), 1.0e-16, solve_from, false,
    )
    @test !improved7
    @test eq7 === eq4
    @test cert7 === failed

end

@testsection "an infeasible budget is neither certified nor made to hang" begin
    # The automatic cascade — continuation, restart from the answer, repair of a
    # missing phase — runs only when the ordinary starts fail to certify, and on
    # a well-posed small system they never do. A budget no composition can meet
    # gets into it, and what must hold there is that the route terminates, says
    # plainly that it has no certificate, and does not fabricate one.
    #
    # `-b` is infeasible by construction: every component budget is negative and
    # every stoichiometric coefficient non-negative, so no `n >= 0` can meet it.
    sp8 = Dict(
        symbol(s) => s for s in build_species(
                datapath("slop98-inorganic-thermofun.json")
            )
    )
    cs8 = ChemicalSystem(
        [sp8[s] for s in split("H2O@ H+ OH- CO2@ HCO3- CO3-2 Ca+2 Cal")],
        ["H2O@", "H+", "Ca+2", "CO3-2", "Zz"],
    )
    A8 = Float64.(cs8.SM.A)
    st8 = ChemicalState(cs8)
    set_quantity!(st8, "H2O@", 1.0u"kg")
    set_quantity!(st8, "Cal", 1.0e-3u"mol")
    b8 = A8 * ustrip.(us"mol", st8.n)

    strict = ChemistryLab.STRICT_CONVERGENCE[]
    try
        ChemistryLab.STRICT_CONVERGENCE[] = false
        # `lp_start = false`: the cascade itself is under test here. With the
        # linear program, the budget is refused before it; see the next section.
        eq8, cert8 = equilibrate_certified(st8; b = -b8, lp_start = false)
        @test !cert8.optimal
        @test eq8 isa ChemicalState
        # A budget with a negative component offers nothing to lift a phase off
        # its bound with, so the repair declines rather than inventing an amount.
        @test ChemistryLab._repair_start(eq8, DiluteSolutionModel(), -b8, 1.0e-16) ===
            nothing

        # A back end that throws is reported under `verbose` and costs the search
        # only that candidate. Registered first, since a start is taken from the
        # first factory that answers; the logger swallows the dual solver's own
        # iteration trace, which `verbose` also turns on.
        factories = ChemistryLab._SOLVER_FACTORIES
        saved = copy(factories)
        try
            pushfirst!(factories, () -> error("this back end is unavailable"))
            # `Base.CoreLogging` rather than `using Logging`, which would have
            # to be declared in the test target for one call.
            eq9, cert9 = Base.CoreLogging.with_logger(Base.CoreLogging.NullLogger()) do
                equilibrate_certified(st8; b = -b8, verbose = true, lp_start = false)
            end
            @test !cert9.optimal
            @test eq9 isa ChemicalState
        finally
            empty!(factories)
            append!(factories, saved)
        end
        @test ChemistryLab._SOLVER_FACTORIES == saved
    finally
        ChemistryLab.STRICT_CONVERGENCE[] = strict
    end
end


@testsection "an infeasible budget is refused before the search, with its reason" begin
    # The linear program over the pure phases proves it: a Farkas vector z with
    # Aᵀz ≥ 0 and bᵀz < 0, a combination of the balances that every species
    # raises and the budget lowers. Until 0.25.2 such a budget went through the
    # whole cascade (continuation, restarts, repairs) and ended uncertified, with
    # the solver's non-convergence counter raised on the way.
    sp8 = Dict(symbol(s) => s for s in build_species(datapath("slop98-inorganic-thermofun.json")))
    cs8 = ChemicalSystem(
        [sp8[s] for s in split("H2O@ H+ OH- CO2@ HCO3- CO3-2 Ca+2 Cal")],
        ["H2O@", "H+", "Ca+2", "CO3-2", "Zz"],
    )
    A8 = Float64.(cs8.SM.A)
    st8 = ChemicalState(cs8)
    set_quantity!(st8, "H2O@", 1.0u"kg")
    set_quantity!(st8, "Cal", 1.0e-3u"mol")
    b8 = A8 * ustrip.(us"mol", st8.n)
    row(name) = findfirst(==(name), [symbol(p) for p in cs8.SM.primaries])
    quiet(f) = Base.CoreLogging.with_logger(f, Base.CoreLogging.NullLogger())
    strict = ChemistryLab.STRICT_CONVERGENCE[]
    try
        ChemistryLab.STRICT_CONVERGENCE[] = false

        # One balance alone: calcium, which every species holds with a
        # coefficient of zero or more, asked for at a negative amount.
        bCa = copy(b8)
        bCa[row("Ca+2")] = -1.0e-3
        before = ChemistryLab.NONCONVERGED[]
        eq, cert = quiet(() -> equilibrate_certified(st8; b = bCa))
        @test cert.budget_feasible === false
        @test cert.route === :infeasible
        @test !cert.optimal
        @test occursin("Ca+2", cert.unplaceable)
        @test cert.n_dual_solves == 0
        @test ChemistryLab.NONCONVERGED[] == before
        @test eq.n == st8.n && eq !== st8

        # Two balances together: no water, and less than no hydrogen ion. OH-
        # would lower the second only by raising the first, and CO2@, which
        # takes water, gives back twice the hydrogen; the reason names both rows.
        b2 = copy(b8)
        b2[row("H2O@")] = 0.0
        b2[row("H+")] = -1.0e-3
        eq, cert = quiet(() -> equilibrate_certified(st8; b = b2))
        @test cert.budget_feasible === false
        @test occursin("H2O@", cert.unplaceable) && occursin("H+", cert.unplaceable)

        # Under the strict flag it raises.
        ChemistryLab.STRICT_CONVERGENCE[] = true
        @test_throws "no non-negative amounts" equilibrate_certified(st8; b = bCa)
    finally
        ChemistryLab.STRICT_CONVERGENCE[] = strict
    end

    # A feasible budget says so, and which start the answer came from.
    eq, cert = equilibrate_certified(st8)
    @test cert.optimal && cert.budget_feasible === true
    @test cert.route in (:lp_start, :state)
    @test cert.n_dual_solves >= 1
    eq0, cert0 = equilibrate_certified(st8; lp_start = false)
    @test cert0.route === :state
    @test ustrip.(us"mol", eq.n) ≈ ustrip.(us"mol", eq0.n) rtol = 1.0e-8
end

@testsection "the start from the linear program puts each species in its phase" begin
    # A tenth of a kilogram of water and a millimole of calcite, beside
    # portlandite. The vertex holds the water and the calcite. A solute it leaves
    # out starts at the molality its multipliers give it, per kilogram of that
    # water, and portlandite, a pure phase outside the vertex, starts absent.
    # Until 0.28.0 every species outside the vertex started at its activity in
    # moles, portlandite included, which on a cement put the start tens of moles
    # off the budget and cost the search seconds from which nothing certified.
    sp = Dict(symbol(s) => s for s in build_species(datapath("cemdata18-thermofun.json"); verbose = false))
    cs = ChemicalSystem(
        [sp[s] for s in split("H2O@ H+ OH- CO2@ HCO3- CO3-2 Ca+2 CaOH+ Cal Portlandite")],
        ["H2O@", "H+", "Ca+2", "CO3-2", "Zz"],
    )
    st = ChemicalState(cs)
    set_quantity!(st, "H2O@", 0.1u"kg")
    set_quantity!(st, "Cal", 1.0e-3u"mol")
    b = Float64.(cs.SM.A) * ustrip.(us"mol", st.n)
    des = ChemistryLab.DualEquilibriumSolver(cs, DiluteSolutionModel())
    lp = ChemistryLab._linear_program(des, st, b)
    @test lp.start.status === :optimal
    ϵ = ChemistryLab._AMOUNT_FLOOR
    lifted = ustrip.(us"mol", ChemistryLab._lp_lifted_state(des, st, lp, ϵ).n)
    x = lp.start.x
    j(name) = findfirst(s -> symbol(s) == name, cs.species)
    @test lifted[j("H2O@")] == x[j("H2O@")] > 0 && lifted[j("Cal")] == x[j("Cal")] > 0
    @test x[j("Portlandite")] == 0 && lifted[j("Portlandite")] == 0
    kg = x[j("H2O@")] * ustrip(us"kg/mol", sp["H2O@"][:M])
    u = -(transpose(lp.prob.A) * lp.start.y)
    outside = [s for s in ("H+", "OH-", "CO2@", "HCO3-", "CO3-2", "Ca+2", "CaOH+") if x[j(s)] == 0]
    @test !isempty(outside)
    for s in outside
        @test lifted[j(s)] ≈ max(exp(min(u[j(s)] - lp.prob.g[j(s)], 0.0)) * kg, ϵ) rtol = 1.0e-12
    end
    # And no solute starts below the floor of the search, whatever the
    # multipliers give it: at the pH of a cement they give H+ about 1e-100, far
    # below anything the search resolves. A floor of a millimole makes the rule
    # visible on every solute outside the vertex; the vertex itself and the
    # absent pure phase are untouched.
    high = ustrip.(us"mol", ChemistryLab._lp_lifted_state(des, st, lp, 1.0e-3).n)
    @test all(high[j(s)] >= 1.0e-3 for s in outside)
    @test high[j("H2O@")] == x[j("H2O@")] && high[j("Cal")] == x[j("Cal")] && high[j("Portlandite")] == 0
    # The search from it gives the answer the search without it gives.
    eq, cert = equilibrate_certified(st)
    eq0, cert0 = equilibrate_certified(st; lp_start = false)
    @test cert.optimal && cert0.optimal
    @test ustrip.(us"mol", eq.n) ≈ ustrip.(us"mol", eq0.n) rtol = 1.0e-8

    # A solute held far below the amount the potentials give it is not certified.
    # The certificate leaves a member below its floor (1e-25 mol) out of the
    # equality, as the truncation of a smaller exact amount, and until
    # OptimaSolver 0.7.2 it held it to nothing: this answer with its H+ at
    # 1e-100 mol certified, and a cement answer did, its pH read 0.09 high.
    @test cert.stationarity_floored < 1.0e-10
    n = ustrip.(us"mol", eq.n)
    n[j("H+")] = 1.0e-100
    left = ChemistryLab.optimality_certificate(des, ChemicalState(cs; T = eq.T[1], P = eq.P[1], n = n .* u"mol"); b)
    @test !left.optimal
    @test left.stationarity_floored > 1.0e-6
    @test left.balance < 1.0e-10         # the 1e-11 mol of H+ taken out

    # A solute left below the activity floor at a converged answer is given the
    # amount the potentials of the answer give it, `floor·exp(r)`: the exact
    # solution of the floored model, here the H+ of the equilibrium itself.
    p = ChemistryLab._build_params(st)
    n0 = ustrip.(us"mol", st.n)
    blocks = ChemistryLab._constraint_blocks(FixedTP(), des, st, p, n0)
    prob = ChemistryLab._dual_problem(des, p, n0, blocks)
    res = ChemistryLab._optima_dual_solve(prob, b, max.(n0, p.ϵ), des.opts)
    @test res.converged
    floor = ChemistryLab._activity_floor(p)
    h = j("H+")
    stuck = merge(res, (; x = [i == h ? 1.0e-100 : v for (i, v) in enumerate(res.x)]))
    completed = ChemistryLab._complete_floored_solutes(des, prob, stuck, floor)
    @test completed[h] ≈ res.x[h] rtol = 1.0e-6
    @test completed[[i for i in eachindex(res.x) if i != h]] == res.x[[i for i in eachindex(res.x) if i != h]]
    # Nothing to complete on the answer itself.
    @test ChemistryLab._complete_floored_solutes(des, prob, res, floor) == res.x
end

# The three routes that make a complete phase list usable: refusing an answer
# outside the model's domain, the ideal model as a starting point, and offering a
# solid solution back by its own criterion. Their own fixtures, since each needs a
# system the earlier sections do not build.
@testsection "what makes a complete phase list usable" begin

    sp2 = Dict(
        symbol(s) => s for s in build_species(
                datapath("slop98-inorganic-thermofun.json"); verbose = false
            )
    )
    species2 = [sp2[s] for s in split("H2O@ H+ OH- CO2@ HCO3- CO3-2 Ca+2 Cal")]
    cs2 = ChemicalSystem(species2, ["H2O@", "H+", "Ca+2", "CO3-2", "Zz"])
    A2 = Matrix{Float64}(cs2.SM.A)

    function calcite2(; nco2 = 0.0)
        st = ChemicalState(cs2)
        set_quantity!(st, "Cal", 1.0e-3u"mol")
        set_quantity!(st, "H2O@", 1.0u"kg")
        nco2 > 0 && set_quantity!(st, "CO2@", nco2 * u"mol")
        V = volume(st)
        set_quantity!(st, "H+", 1.0e-4u"mol/L" * V.liquid)
        set_quantity!(st, "OH-", 1.0e-10u"mol/L" * V.liquid)
        return st
    end

    @testset "ranking refuses an answer outside the model's domain" begin
        # `_check_solvent` already reports a state whose solvent has been taken by
        # the solids; what it could not do is stop one from being CHOSEN. The
        # ranking now compares admissibility first, in both directions, as it
        # already did for the optimality flag.
        good = calcite2()
        @test ChemistryLab._within_domain(good)
        @test solvent_fraction(good) >= ChemistryLab.SOLVENT_FRACTION_FLOOR

        # The same system with the water taken out: `x_w` falls under the floor.
        starved = ChemicalState(cs2)
        set_quantity!(starved, "H2O@", 1.0e-3u"mol")
        set_quantity!(starved, "Ca+2", 1.0u"mol")
        set_quantity!(starved, "CO3-2", 1.0u"mol")
        @test solvent_fraction(starved) < ChemistryLab.SOLVENT_FRACTION_FLOOR
        @test !ChemistryLab._within_domain(starved)

        # A system with no solvent at all has no such question to answer.
        @test ChemistryLab._within_domain(
            ChemicalState(ChemicalSystem([sp2["Ca+2"], sp2["CO3-2"]], ["Ca+2", "CO3-2"]))
        )

        # And the ranking uses it: an inadmissible answer loses to an admissible
        # one even when its residuals are smaller.
        c(opt, err) = (;
            optimal = opt, stationarity = err, balance = err,
            worst_supersaturation = 0.0,
        )
        kb = ChemistryLab._keep_better
        @test first(kb(starved, c(false, 1.0e-12), good, c(false, 1.0))) === good
        @test first(kb(good, c(false, 1.0), starved, c(false, 1.0e-12))) === good
        # a certified answer still beats an admissible uncertified one
        @test first(kb(good, c(false, 1.0e-16), starved, c(true, 1.0))) === starved
    end

    @testset "the ideal model is a usable starting point" begin
        # `_ideal_start` answers the same question without activity coefficients,
        # which is better conditioned and certifies where the non-ideal model may
        # not. What it returns is an equilibrium of the ideal problem, so the
        # non-ideal solve that starts from it begins with the right active set.
        st = calcite2(; nco2 = 0.01)
        b = A2 * ustrip.(us"mol", st.n)

        ideal = ChemistryLab._ideal_start(
            st, HKFActivityModel(), b, 1.0e-16, FixedTP(), false
        )
        @test ideal isa ChemicalState
        _, cert = equilibrate_certified(ideal; model = DiluteSolutionModel(), b = b)
        @test cert.optimal

        # Already ideal: there is no easier question to fall back on.
        @test ChemistryLab._ideal_start(
            st, DiluteSolutionModel(), b, 1.0e-16, FixedTP(), false
        ) === nothing

        # A failure inside is swallowed: this builds a start, and the caller is
        # owed the outer verdict rather than an error raised in a heuristic.
        @test ChemistryLab._ideal_start(
            st, HKFActivityModel(), Float64[], 1.0e-16, FixedTP(), false
        ) === nothing
    end

    @testset "the ideal pre-solve is reached only when nothing else certifies" begin
        # The call site, not the helper: it runs only under `autostart`, only when
        # no route has certified, and only for a non-ideal model. An infeasible
        # budget — every component negative against non-negative stoichiometry —
        # never certifies, so the cascade is entered and the pre-solve with it.
        st = calcite2()
        b = A2 * ustrip.(us"mol", st.n)
        strict = ChemistryLab.STRICT_CONVERGENCE[]
        try
            ChemistryLab.STRICT_CONVERGENCE[] = false
            eq, cert = equilibrate_certified(st; model = HKFActivityModel(), b = -b, lp_start = false)
            @test !cert.optimal
            @test eq isa ChemicalState
        finally
            ChemistryLab.STRICT_CONVERGENCE[] = strict
        end
    end

    @testset "a cold cement is rescued by the ideal pre-solve" begin
        # The call site doing its job, on the case that showed it missing. The
        # stage had been unwired by accident when the linear-programming start
        # was removed, and nothing here noticed: the test above reaches it only
        # with a budget no stage can rescue. This is a 109-species paste of
        # clinker and gypsum at w/b = 0.5, with the C-N-A-S-H gel declared, solved
        # from the cast state. Neither back end hands the dual solve a start it
        # can certify from -- the search alone ends with an element balance of
        # 0.12 -- and the ideal answer does: the same paste then certifies, at
        # the pH the CEM IV page reports for it.
        substances = build_species(datapath("cemdata18-thermofun.json"); verbose = false)
        byname = Dict(symbol(s) => s for s in substances)
        M(n) = ustrip(us"g/mol", byname[n][:M])
        pure = split(
            "C3S C2S C3A C4AF Gp Anh Cal Portlandite ettringite monosulphate12 " *
                "monocarbonate hemicarbonate C4AH13 C3AH6 C3FH6 straetlingite " *
                "hydrotalcite Brc FeOOHmic AlOHmic Amor-Sl Mgs " *
                "C3AS0.84H4.32 C3AS0.41H5.18 straetlingite7 Gbs AlOHam " *
                "M4A-OH-LDH M6A-OH-LDH M8A-OH-LDH C2AH7.5 C4AH11 C4AH19 " *
                "K2SO4 syngenite Na2SO4"
        )
        gel = [
            "T2C-CNASHss", "T5C-CNASHss", "TobH-CNASHss",
            "5CA", "5CNA", "INFCA", "INFCN", "INFCNA",
        ]
        feal = ["C3AFS0.84H4.32", "C3FS0.84H4.32"]
        cs = ChemicalSystem(
            speciation(
                substances, vcat(pure, gel, feal, ["SO4-2", "CO2@", "O2@"]);
                aggregate_state = [AS_AQUEOUS],
            ),
            CEMDATA_PRIMARIES;
            solid_solutions = [
                SolidSolutionPhase("CNASH_ss", [byname[m] for m in gel]),
                SolidSolutionPhase("C3(AF)S0.84H", [byname[m] for m in feal]),
            ],
        )
        st = ChemicalState(cs)
        # The CEM I of Lavergne et al. (2018), Table 9: its gypsum and the
        # Bogue composition of its clinker.
        gypsum = literature_value("Lavergne2018", "gypsum_percent")   # g, of 100 g binder
        bogue = literature_table("Lavergne2018", "cement_bogue")
        clinker = 100.0 * (1 - gypsum / 100)                         # g, of 100 g binder
        for (p, f) in zip(bogue.phase, bogue.percent ./ 100)
            set_quantity!(st, p, clinker * f / M(p) * u"mol")
        end
        set_quantity!(st, "Gp", gypsum / M("Gp") * u"mol")
        set_quantity!(st, "H2O@", 50.0 / M("H2O@") * u"mol")
        b = Float64.(cs.SM.A) * ustrip.(us"mol", st.n)
        b .+= oxide_budget(
            OrderedDict("K2O" => 0.008, "Na2O" => 0.002), cs.SM.primaries;
            mass = clinker * u"g",
        )
        model = HKFActivityModel(å = 0.0, Ḃ = gems_bdot(), Kₙ = 0.0)

        strict = ChemistryLab.STRICT_CONVERGENCE[]
        try
            ChemistryLab.STRICT_CONVERGENCE[] = false
            quiet(f) = Base.CoreLogging.with_logger(f, Base.CoreLogging.NullLogger())
            _, declined = quiet(
                () -> equilibrate_certified(st; model = model, b = b, autostart = false)
            )
            # The premise, and it is a verdict of rounding rather than of
            # chemistry: the search alone is refused on Linux and on Windows alike,
            # but with element balances of 0.12 and 9.1e-5, so no threshold on it
            # is portable. A platform whose arithmetic lets it certify at once has
            # nothing to rescue, and that is reported rather than failed.
            declined.optimal && @info "the search certified without the " *
                "automatic stages on this platform; the rescue is not exercised" declined.balance

            @test ChemistryLab._ideal_start(
                st, model, b, 1.0e-16, FixedTP(), false
            ) isa ChemicalState
            # Without the linear-programming start, which reaches the answer
            # before the ideal pre-solve is needed.
            eq, cert = quiet(() -> equilibrate_certified(st; model = model, b = b, lp_start = false))
            @test cert.optimal
            @test cert.balance < 1.0e-10
            @test pH(eq, model) ≈ 13.444 atol = 1.0e-3

            # With it, the default: the same answer from the first start, in
            # fewer dual solves.
            eqlp, certlp = quiet(() -> equilibrate_certified(st; model = model, b = b))
            @test certlp.optimal
            @test certlp.route === :lp_start
            @test certlp.n_dual_solves < cert.n_dual_solves
            @test ustrip.(us"mol", eqlp.n) ≈ ustrip.(us"mol", eq.n) rtol = 1.0e-8 atol = 1.0e-14

            # The same paste with the CNASH_ss that ships, Myers' own model: ideal
            # on six sites rather than between the eight end-members. The solver
            # is told to invert that phase by Newton's method, the substitution
            # diverging on it, and which members may leave it; it certifies under
            # the activity model Cemdata18 prescribes. Under the limiting law
            # above (å = 0) it does not: measured, the search ends uncertified
            # at pH 14.8 after five minutes, where ideal mixing certifies at 13.44.
            shipped = build_solid_solutions(datapath("solid_solutions.toml"), byname)
            cnash = only(filter(p -> ChemistryLab.name(p) == "CNASH_ss", shipped))
            @test ChemistryLab.model(cnash) isa SublatticeModel
            cs_sl = ChemicalSystem(
                cs.species, cs.SM.primaries;
                solid_solutions = [cnash, SolidSolutionPhase("C3(AF)S0.84H", [byname[m] for m in feal])],
            )
            st_sl = ChemicalState(cs_sl; n = st.n)
            c18 = cemdata18_activity_model(:KOH)
            des = DualEquilibriumSolver(cs_sl, c18)
            handed = ChemistryLab._dual_phases(des, ustrip.(us"mol", st_sl.n))
            gel_phase = only(filter(p -> p.newton, handed))
            @test [symbol(cs_sl.species[gel_phase.members[j]]) for j in gel_phase.bounded_members] ==
                ["T5C-CNASHss", "5CA"]
            @test gel_phase.local_h isa Function
            eqs, certs = quiet(() -> equilibrate_certified(st_sl; model = c18, b = b))
            @test certs.optimal
            @test certs.balance < 1.0e-10
        finally
            ChemistryLab.STRICT_CONVERGENCE[] = strict
        end
    end

    @testset "a vanished aqueous phase is reported, and raised under the strict flag" begin
        # `_check_solvent` is the diagnosis `_within_domain` turned into a ranking.
        # Both halves of its verdict are exercised here.
        starved = ChemicalState(cs2)
        set_quantity!(starved, "H2O@", 1.0e-3u"mol")
        set_quantity!(starved, "Ca+2", 1.0u"mol")
        set_quantity!(starved, "CO3-2", 1.0u"mol")

        strict = ChemistryLab.STRICT_CONVERGENCE[]
        try
            ChemistryLab.STRICT_CONVERGENCE[] = false
            @test_logs (:warn, r"aqueous phase has effectively vanished") match_mode = :any (
                ChemistryLab._check_solvent(starved)
            )
            ChemistryLab.STRICT_CONVERGENCE[] = true
            @test_throws ErrorException ChemistryLab._check_solvent(starved)
        finally
            ChemistryLab.STRICT_CONVERGENCE[] = strict
        end

        # An answer inside the domain says nothing at all.
        @test ChemistryLab._check_solvent(calcite2()) === nothing
    end

    @testset "a missing solid solution is offered back by its own criterion" begin
        # `_repair_start` skips solid-solution end-members, and rightly: the
        # saturation index of a member at the bound reports a small mole fraction,
        # not a phase that should form. The PHASE has its own criterion —
        # `Omega = sum_i 10^SI_i(pure) > 1` — and without it a search that has lost
        # a solid solution cannot get it back, because no member is ever offered.
        cs_ss = ChemicalSystem(
            [sp2[s] for s in split("H2O@ H+ OH- CO2@ HCO3- CO3-2 Ca+2 Cal Arg")],
            ["H2O@", "H+", "Ca+2", "CO3-2", "Zz"];
            solid_solutions = [
                SolidSolutionPhase("carbonate", [sp2["Cal"], sp2["Arg"]]),
            ],
        )
        st = ChemicalState(cs_ss)
        set_quantity!(st, "H2O@", 1.0u"kg")
        set_quantity!(st, "Ca+2", 0.05u"mol")
        set_quantity!(st, "CO3-2", 0.05u"mol")
        bfix = Float64.(cs_ss.SM.A) * ustrip.(us"mol", st.n)

        # Both end-members sit at the floor while the solution is loaded with
        # calcium carbonate, so the phase is supersaturated and comes back in.
        repaired = ChemistryLab._repair_start(st, DiluteSolutionModel(), bfix, 1.0e-16)
        @test repaired isa ChemicalState
        n = ustrip.(us"mol", repaired.n)
        grp = only(cs_ss.ss_groups)
        @test any(n[i] > 1.0e-8 for i in grp)
        @test all(n[i] >= 0 for i in grp)
    end

end

@testsection "equilibrate_path — a sweep continued from its own answers" begin
    sp = Dict(
        symbol(s) => s for s in build_species(
                datapath("slop98-inorganic-thermofun.json"); verbose = false
            )
    )
    cs = ChemicalSystem(
        [sp[s] for s in split("H2O@ H+ OH- CO2@ HCO3- CO3-2 Ca+2 Cal")],
        ["H2O@", "H+", "Ca+2", "CO3-2", "Zz"],
    )
    st = ChemicalState(cs)
    set_quantity!(st, "H2O@", 1.0u"kg")
    set_quantity!(st, "Cal", 0.05u"mol")
    b0 = Float64.(cs.SM.A) * ustrip.(us"mol", st.n)
    budgets = [b0 .* f for f in (1.0, 1.1, 1.2)]

    states, certs = equilibrate_path(st, budgets)
    @test length(states) == 3
    @test length(certs) == 3
    @test all(c.optimal for c in certs)

    # THE ASSERTION THAT MATTERS. On a convex problem the minimum is unique, so
    # walking to it from a neighbor cannot change what is found -- only whether
    # it is found. If this drifts, the continuation is choosing answers rather
    # than reaching them.
    for (k, b) in enumerate(budgets)
        eq, c = equilibrate_certified(st; b = b)
        @test c.optimal
        @test ustrip.(us"mol", eq.n) ≈ ustrip.(us"mol", states[k].n) rtol = 1.0e-6
    end
end

@testsection "the split seed moves material without moving the budget" begin
    # THE ARITHMETIC OF THE SEED, on its own. `equilibrate_split` needs a system
    # that fails to certify AND reports an incipient composition before its loop
    # runs at all, which on this package's chemistry means a 91-species cement.
    # The transfer itself is pure arithmetic on a vector of moles, so it is
    # tested here as such — and it is the part that was wrong.
    #
    # WHAT WAS WRONG. The first version removed material at the DONOR's
    # composition and added it at the trial's. That moves the right number of
    # moles and the wrong mixture: the total of each end-member changes, so the
    # seed no longer satisfies the element budget the solver is about to be
    # measured against. Two end-members of one binary are different substances —
    # `C4AH13` and `monosulphate12` do not have the same sulfur.
    #
    # Species 1,2 are the base instance; 3,4 its `#2` twin.
    twin = Dict(1 => 3, 2 => 4)
    untwin = Dict(3 => 1, 4 => 2)
    trial(members, x) = Dict(1 => (members = members, x = x))

    # 1. THE CONSERVATION. Whatever moves, the total of each end-member across
    #    the pair is unchanged — to the last bit, since the same `t.x` leaves one
    #    instance and enters the other.
    let n = [0.3, 0.1, 0.0, 0.0]
        before = (n[1] + n[3], n[2] + n[4])
        @test ChemistryLab._seed_split!(n, trial([1, 2], [0.25, 0.75]), twin, untwin, 0.5)
        @test n[1] + n[3] ≈ before[1] atol = 0.0
        @test n[2] + n[4] ≈ before[2] atol = 0.0
        @test all(n .>= -1.0e-15)
        # And the receiving instance holds exactly the trial composition.
        tot = n[3] + n[4]
        @test tot > 0
        @test n[4] / tot ≈ 0.75 rtol = 1.0e-12
    end

    # 2. FROM THE FULLER INTO THE EMPTIER, whichever instance the trial names.
    #    The trial below names the TWIN, and the base still holds everything, so
    #    the base is the donor.
    let n = [0.3, 0.1, 0.0, 0.0]
        @test ChemistryLab._seed_split!(n, trial([3, 4], [0.25, 0.75]), twin, untwin, 0.5)
        @test n[3] + n[4] > 0                       # the twin received
        @test n[1] + n[2] < 0.4                    # the base gave
    end
    let n = [0.0, 0.0, 0.3, 0.1]                  # everything in the twin
        @test ChemistryLab._seed_split!(n, trial([1, 2], [0.25, 0.75]), twin, untwin, 0.5)
        @test n[1] + n[2] > 0                       # now the base receives
        @test n[3] + n[4] < 0.4
    end

    # 3. `move` IS BOUNDED so no end-member of the donor goes negative. Here the
    #    donor holds almost no member 2 and the trial asks for three quarters of
    #    it, so `share = 0.5` cannot be honored in full — and must not be.
    let n = [0.4, 0.001, 0.0, 0.0]
        @test ChemistryLab._seed_split!(n, trial([1, 2], [0.25, 0.75]), twin, untwin, 0.5)
        @test all(n .>= -1.0e-15)
        @test n[2] ≈ 0.0 atol = 1.0e-12             # drained, not overdrawn
        @test n[1] + n[3] ≈ 0.4 atol = 1.0e-15     # still conserved
    end

    # 4. A PAIR IS SEEDED ONCE even when both of its instances are flagged, and a
    #    trial on a phase with no twin is left alone.
    let n = [0.3, 0.1, 0.0, 0.0]
        both = Dict(
            1 => (members = [1, 2], x = [0.25, 0.75]),
            2 => (members = [3, 4], x = [0.25, 0.75]),
        )
        @test ChemistryLab._seed_split!(n, both, twin, untwin, 0.5)
        @test n[1] + n[3] ≈ 0.3 atol = 0.0
        # One transfer, not two — and of the BOUNDED amount. `share = 0.5` of a
        # 0.4 mol pair asks for 0.2, but the donor is 75/25 while the trial wants
        # 25/75, so the scarce end-member runs out first: 0.1 / 0.75 = 0.1333.
        # Seeded twice, this would read 0.2444.
        @test n[3] + n[4] ≈ 0.1 / 0.75 rtol = 1.0e-12
    end
    let n = [0.3, 0.1, 0.0, 0.0]
        @test !ChemistryLab._seed_split!(
            n, Dict(1 => (members = [7, 8], x = [0.5, 0.5])), twin, untwin, 0.5
        )
        @test n == [0.3, 0.1, 0.0, 0.0]           # untouched
    end

    # 5. Nothing to move: an empty pair is skipped rather than divided by zero.
    let n = [0.0, 0.0, 0.0, 0.0]
        @test !ChemistryLab._seed_split!(n, trial([1, 2], [0.25, 0.75]), twin, untwin, 0.5)
    end

    # 6. `_instance_pairs` finds the twins, and finds none when there are none.
    let
        sp = Dict(
            symbol(s) => s for s in build_species(
                    datapath("slop98-inorganic-thermofun.json"); verbose = false
                )
        )
        names = split("H2O@ H+ OH- CO2@ HCO3- CO3-2 Ca+2 Mg+2 Cal Mgs")
        comps = ["H2O@", "H+", "Ca+2", "Mg+2", "CO3-2", "Zz"]
        plain = ChemicalSystem([sp[s] for s in names], comps)
        @test isempty(first(ChemistryLab._instance_pairs(plain)))

        doubled = ChemicalSystem(
            [sp[s] for s in names], comps;
            solid_solutions = [
                SolidSolutionPhase(
                    "carbonate", [sp["Cal"], sp["Mgs"]];
                    model = RedlichKisterModel(a0 = 20_000.0), instances = 2,
                ),
            ],
        )
        tw, un = ChemistryLab._instance_pairs(doubled)
        @test length(tw) == 2                      # one pair per end-member
        @test length(un) == 2
        syms = String.(symbol.(doubled.species))
        for (i, j) in tw
            @test syms[j] == syms[i] * "#2"        # and they are the right ones
            @test un[j] == i
        end
    end
end

@testsection "equilibrate_split — the split loop, and the pair it is about" begin
    # THE REGRESSION THIS GUARDS. `equilibrate_split` reads the incipient
    # composition out of the certificate and seeds a second instance with it.
    # Two links in that chain were broken and neither could be seen by reading
    # the code: `optimality_certificate` rebuilt its NamedTuple and DROPPED
    # `split_trials`, so the loop always read `nothing` and returned on its first
    # pass; and the pass was accepted on `cert.worst_violation`, a field this
    # package's certificate does not have, so reaching that line at all raised a
    # `FieldError`. Both are the same defect in the end -- a function that had
    # never been run -- and a test that runs it is the only thing that catches
    # them.
    sp = Dict(
        symbol(s) => s for s in build_species(
                datapath("slop98-inorganic-thermofun.json"); verbose = false
            )
    )

    # 1. WITHOUT instances there is nothing to split into, and the function must
    #    be exactly `equilibrate_certified` -- same answer, not merely a good one.
    let
        names = split("H2O@ H+ OH- CO2@ HCO3- CO3-2 Ca+2 Cal Arg")
        cs = ChemicalSystem(
            [sp[s] for s in names], ["H2O@", "H+", "Ca+2", "CO3-2", "Zz"];
            solid_solutions = [SolidSolutionPhase("carbonate", [sp["Cal"], sp["Arg"]])],
        )
        st = ChemicalState(cs)
        set_quantity!(st, "H2O@", 1.0u"kg")
        set_quantity!(st, "Cal", 0.05u"mol")
        b = Float64.(cs.SM.A) * ustrip.(us"mol", st.n)
        eq_s, c_s = equilibrate_split(st; b = b)
        eq_c, c_c = equilibrate_certified(st; b = b)
        @test c_s.optimal == c_c.optimal
        @test ustrip.(us"mol", eq_s.n) ≈ ustrip.(us"mol", eq_c.n) rtol = 1.0e-8
    end

    # 2. THE PAIR ITSELF, against an analytic oracle. Calcite and magnesite are
    #    two different substances, so the element budget PINS the overall
    #    composition: 0.025 mol of each fixes x̄ = 1/2 whatever the energetics
    #    say. Put a Redlich-Kister gap on that binary and x̄ = 1/2 is inside it,
    #    so the Gibbs minimum is two coexisting compositions — and with two
    #    instances declared, the minimization finds them. They must be the
    #    common-tangent pair, which [`common_tangent`](@ref) computes from the
    #    mixing model alone and which nothing in the solve has been told.
    #
    #    This is the case the AFm binary of a real CEM I is NOT: there the
    #    sulfate has somewhere else to go (ettringite) and the hydroxide is
    #    abundant, so nothing pins the phase's composition and both instances sit
    #    at the same x. Pinned, the pair comes out by itself.
    names = split("H2O@ H+ OH- CO2@ HCO3- CO3-2 Ca+2 Mg+2 Cal Mgs")
    comps = ["H2O@", "H+", "Ca+2", "Mg+2", "CO3-2", "Zz"]
    for a0 in (8_000.0, 14_000.0, 20_000.0)
        gap = RedlichKisterModel(a0 = a0)
        @test spinodal_interval(gap, 2) !== nothing        # it IS a gap
        pair = common_tangent(gap)
        @test pair !== nothing

        cs = ChemicalSystem(
            [sp[s] for s in names], comps;
            solid_solutions = [
                SolidSolutionPhase(
                    "carbonate", [sp["Cal"], sp["Mgs"]]; model = gap, instances = 2
                ),
            ],
        )
        @test length(cs.species) == length(names) + 2      # the `#2` twins exist
        st = ChemicalState(cs)
        set_quantity!(st, "H2O@", 1.0u"kg")
        set_quantity!(st, "Cal", 0.025u"mol")
        set_quantity!(st, "Mgs", 0.025u"mol")
        b = Float64.(cs.SM.A) * ustrip.(us"mol", st.n)

        eq, cert = equilibrate_certified(st; b = b)
        @test cert.optimal
        n = ustrip.(us"mol", eq.n)
        xs = Float64[]
        totals = Float64[]
        for g in cs.ss_groups
            tot = sum(n[i] for i in g)
            tot > 1.0e-10 || continue
            push!(totals, tot)
            push!(xs, n[g[2]] / tot)
        end
        @test length(xs) == 2
        sort!(xs)
        # The two instances sit ON the binodal, computed independently.
        @test xs[1] ≈ pair[1] atol = 1.0e-3
        @test xs[2] ≈ pair[2] atol = 1.0e-3
        # And in the proportions the lever rule asks for at x̄ = 1/2.
        f = totals[1] / sum(totals)
        @test f ≈ (pair[2] - 0.5) / (pair[2] - pair[1]) atol = 5.0e-3

        # 3. THE WIRING. Whatever the verdict, the certificate must carry the
        #    split diagnostics through -- they are computed in `OptimaSolver` and
        #    were being discarded on the way out.
        @test hasproperty(cert, :split_trials)
        @test hasproperty(cert, :split_phases)
        @test hasproperty(cert, :worst_violation_split)
        for (_, t) in cert.split_trials
            @test length(t.x) == length(t.members)
            @test all(t.x .>= 0)
            @test sum(t.x) ≈ 1.0 rtol = 1.0e-6
            @test all(1 <= i <= length(cs.species) for i in t.members)
        end

        # 4. THE PROMISE. A pass is kept only when the KKT error improves, so the
        #    result is never worse than the answer without splitting. This is
        #    what the `FieldError` used to prevent from ever being evaluated.
        eq2, cert2 = equilibrate_split(st; b = b, maxpasses = 2)
        @test ChemistryLab._kkt_error(cert2) <=
            ChemistryLab._kkt_error(cert) * (1 + 1.0e-8)

        # And at the temperature it was given: the seeded state was built at the
        # default 25 °C until 0.24.0.
        st20 = ChemicalState(cs, st.n; T = 293.15u"K")
        eq20, _ = equilibrate_split(st20; b = b, maxpasses = 2)
        @test temperature(eq20) == temperature(st20)

        # 5. `share` is a fraction and is checked, not trusted.
        @test_throws ArgumentError equilibrate_split(st; b = b, share = 0.0)
        @test_throws ArgumentError equilibrate_split(st; b = b, share = 1.0)
    end

    # 6. THE SEED, reached. Between the binodal and the spinodal a composition is
    #    metastable: started with both instances at x̄ = 0.12, the minimization
    #    has no descent direction and stays there, the certificate refuses it and
    #    names the incipient phase, and only the seed leads to the pair. At 20 °C,
    #    and the pair is that of 20 °C: until 0.24.0 the seeded state was built at
    #    25 °C, where the binodal of this gap lies elsewhere (0.053 against 0.049).
    let gap = RedlichKisterModel(a0 = 8_000.0), xbar = 0.12, T = 293.15
        cs = ChemicalSystem(
            [sp[s] for s in names], comps;
            solid_solutions = [
                SolidSolutionPhase(
                    "carbonate", [sp["Cal"], sp["Mgs"]]; model = gap, instances = 2
                ),
            ],
        )
        pair = common_tangent(gap; T)
        @test pair[1] < xbar < spinodal_interval(gap, 2; T)[1]   # metastable, not unstable
        st = ChemicalState(cs; T = T * u"K")
        set_quantity!(st, "H2O@", 1.0u"kg")
        for (s, x) in (("Cal", 1 - xbar), ("Mgs", xbar), ("Cal#2", 1 - xbar), ("Mgs#2", xbar))
            set_quantity!(st, s, 0.025x * u"mol")
        end
        b = Float64.(cs.SM.A) * ustrip.(us"mol", st.n)
        x_of(eq) = sort([ustrip(us"mol", eq.n[g[2]]) / sum(ustrip(us"mol", eq.n[i]) for i in g) for g in cs.ss_groups])

        eq0, c0 = equilibrate_certified(st; b = b, autostart = false)
        @test !c0.optimal
        @test !isempty(c0.split_trials)
        @test x_of(eq0)[2] - x_of(eq0)[1] < 1.0e-6       # both instances where they started

        eq, c = equilibrate_split(st; b = b, autostart = false)
        @test c.optimal
        @test temperature(eq) == temperature(st)
        @test x_of(eq) ≈ collect(pair) atol = 1.0e-3

        # 7. UNDER THE STRICT FLAG. The first pass is, by construction, the one
        #    that does not certify (`c0` above). Until 0.25.1 the flag was
        #    honored there, and the function raised before it had seeded
        #    anything. The passes now run relaxed; only the final answer is
        #    judged strictly.
        strict = ChemistryLab.STRICT_CONVERGENCE[]
        try
            ChemistryLab.STRICT_CONVERGENCE[] = true
            eq_s, c_s = equilibrate_split(st; b = b, autostart = false)
            @test c_s.optimal
            @test x_of(eq_s) ≈ collect(pair) atol = 1.0e-3
            # And an answer that does not certify still raises: with no pass
            # allowed, the unseeded first answer is the final one.
            @test_throws ErrorException equilibrate_split(
                st; b = b, autostart = false, maxpasses = 0,
            )
        finally
            ChemistryLab.STRICT_CONVERGENCE[] = strict
        end
    end
end

@testsection "instances = :auto: a second composition only when the phase wants it" begin
    sp = Dict(
        symbol(s) => s for s in build_species(
                datapath("slop98-inorganic-thermofun.json"); verbose = false
            )
    )
    names = split("H2O@ H+ OH- CO2@ HCO3- CO3-2 Ca+2 Mg+2 Cal Mgs")
    comps = ["H2O@", "H+", "Ca+2", "Mg+2", "CO3-2", "Zz"]
    gap = RedlichKisterModel(a0 = 14_000.0)
    pair = common_tangent(gap)
    carbonate(instances) = SolidSolutionPhase("carbonate", [sp["Cal"], sp["Mgs"]]; model = gap, instances)
    system(instances) = ChemicalSystem([sp[s] for s in names], comps; solid_solutions = [carbonate(instances)])
    # `fresh`, not `st`: a closure assigning a name the enclosing testset also
    # assigns would rebind the testset's own `st`.
    function paste(cs, cal, mgs)
        fresh = ChemicalState(cs)
        set_quantity!(fresh, "H2O@", 1.0u"kg")
        set_quantity!(fresh, "Cal", cal * u"mol")
        set_quantity!(fresh, "Mgs", mgs * u"mol")
        return fresh, Float64.(cs.SM.A) * ustrip.(us"mol", fresh.n)
    end
    compositions(eq) = begin
        n = ustrip.(us"mol", eq.n)
        sort(
            [
                (n[g[2]] / sum(n[i] for i in g), sum(n[i] for i in g)) for g in eq.system.ss_groups
                    if sum(n[i] for i in g) > 1.0e-10
            ]
        )
    end

    # The declaration: one composition, room for a second, refused where there is no gap.
    auto = carbonate(:auto)
    @test auto.instances == 1 && auto.max_instances == 2
    @test occursin("(:auto)", sprint(show, auto))
    @test_throws ErrorException SolidSolutionPhase("ideal", [sp["Cal"], sp["Mgs"]]; instances = :auto)
    @test_throws ArgumentError carbonate(:many)

    # `with_instances` rebuilds the system: the copies exist, the budget reads the same.
    cs = system(:auto)
    @test length(cs.species) == length(names)
    cs2 = with_instances(cs, "carbonate" => 2)
    @test length(cs2.species) == length(names) + 2
    @test symbol.(cs2.SM.primaries) == symbol.(cs.SM.primaries)
    st, b = paste(cs, 0.025, 0.025)
    st2 = with_instances(st, cs2)
    @test Float64.(cs2.SM.A) * ustrip.(us"mol", st2.n) ≈ b rtol = 1.0e-14
    @test temperature(st2) == temperature(st)
    @test_throws ArgumentError with_instances(cs, "calcite" => 2)
    @test length(with_instances(cs2, "carbonate" => 2).species) == length(cs2.species)

    # Pinned inside the gap by the budget: the second instance is added, and the
    # pair is the common tangent, in the proportions of the lever rule -- the
    # answer of the phase declared with two instances from the start.
    eq, cert = equilibrate_certified(st; b = b)
    @test cert.optimal
    @test cert.instances_added == ["carbonate"]
    @test length(eq.system.species) == length(names) + 2
    (x1, t1), (x2, t2) = compositions(eq)
    @test x1 ≈ pair[1] atol = 1.0e-3
    @test x2 ≈ pair[2] atol = 1.0e-3
    @test t1 / (t1 + t2) ≈ (pair[2] - 0.5) / (pair[2] - pair[1]) atol = 5.0e-3
    st_two, b_two = paste(system(2), 0.025, 0.025)
    eq_two, _ = equilibrate_certified(st_two; b = b_two)
    @test first.(compositions(eq)) ≈ first.(compositions(eq_two)) atol = 1.0e-6

    # Never worse than one instance. A first answer that asks to split, but whose
    # error the second instance cannot lower (here made zero), sends the solve back
    # to the ordinary search on the one-instance system.
    first_eq, first_cert = ChemistryLab._relaxed_convergence() do
        ChemistryLab._certified_with_fallback(st; b = b)
    end
    @test ChemistryLab._wants_auto_split(cs, first_cert)
    unbeatable = merge(first_cert, (; stationarity = 0.0, balance = 0.0, worst_supersaturation = -1.0))
    back_eq, back_cert = ChemistryLab._relaxed_convergence() do
        ChemistryLab._after_first_search(st, first_eq, unbeatable; b = b)
    end
    @test back_eq.system === cs
    @test !hasproperty(back_cert, :instances_added)

    # A budget no amounts can meet is refused as usual, and under the strict flag
    # it raises once the second instance has been considered.
    bad = copy(b)
    bad[findfirst(==("Ca+2"), String.(symbol.(cs.SM.primaries)))] = -1.0e-3
    strict = ChemistryLab.STRICT_CONVERGENCE[]
    try
        ChemistryLab.STRICT_CONVERGENCE[] = true
        @test_throws ErrorException equilibrate_certified(st; b = bad)
    finally
        ChemistryLab.STRICT_CONVERGENCE[] = strict
    end

    # Outside the gap nothing is added: one instance, the system unchanged. The
    # pair of this model is close to the pure end-members, so outside is half way
    # between its calcite side and pure calcite.
    x_out = (1 + pair[2]) / 2
    st_out, b_out = paste(cs, 0.05x_out, 0.05(1 - x_out))
    eq_out, cert_out = equilibrate_certified(st_out; b = b_out)
    @test cert_out.optimal
    @test !hasproperty(cert_out, :instances_added)
    @test eq_out.system === cs

    # Under the strict flag the one-composition solve is a search, not a result.
    strict = ChemistryLab.STRICT_CONVERGENCE[]
    try
        ChemistryLab.STRICT_CONVERGENCE[] = true
        _, cert_strict = equilibrate_certified(st; b = b)
        @test cert_strict.optimal
    finally
        ChemistryLab.STRICT_CONVERGENCE[] = strict
    end

    # And from the solid-solution file.
    mktempdir() do dir
        path = joinpath(dir, "ss.toml")
        write(
            path, """
            [[solid_solution]]
            name        = "carbonate"
            end_members = ["Cal", "Mgs"]
            model       = "redlich_kister"
            a0          = 14000.0
            instances   = "auto"
            """
        )
        ss = only(build_solid_solutions(path, sp))
        @test ss.instances == 1 && ss.max_instances == 2
    end

    # A starting point stays in the system it is for. The ideal pre-solve of a
    # budget held inside the gap must not give the phase its second instance:
    # its answer is handed to a search in the one-instance system, and a state
    # with two more species made the dual solve index past the end of the
    # conservation matrix.
    cs_gap = system(:auto)
    st_gap, b_gap = paste(cs_gap, 0.025, 0.025)
    ideal = ChemistryLab._ideal_start(st_gap, HKFActivityModel(), b_gap, 1.0e-16, FixedTP(), false)
    @test ideal === nothing || length(ideal.system.species) == length(cs_gap.species)
    @test ChemistryLab._AUTO_SPLIT[]

    # The split test is read only on an answer that has reached the solution. The
    # first answer inside the gap is stationary and balanced, and asks to split;
    # the same certificate with an element balance of 0.17 is a search that has
    # not arrived, and asks nothing.
    _, first_cert = ChemistryLab._equilibrate_certified(st_gap; b = b_gap, split_early = true)
    @test ChemistryLab._wants_auto_split(cs_gap, first_cert)
    @test !ChemistryLab._wants_auto_split(cs_gap, merge(first_cert, (; balance = 0.17)))
    @test !ChemistryLab._wants_auto_split(cs_gap, merge(first_cert, (; stationarity = 1.0e-4)))
end

@testsection "the certificate says what it proves" begin
    # `optimal = true` is a proof of a global minimum only when the log
    # activities are the gradient of one Gibbs energy and the problem is convex.
    # The certificate measures the first at the composition it audits, from the
    # symmetry of the Jacobian of the log activities, and reads the second off
    # the declarations.
    sp = Dict(
        symbol(s) => s for s in build_species(
                datapath("slop98-inorganic-thermofun.json"); verbose = false
            )
    )
    names = split("H2O@ H+ OH- CO2@ HCO3- CO3-2 Ca+2 Mg+2 Cal Mgs")
    comps = ["H2O@", "H+", "Ca+2", "Mg+2", "CO3-2", "Zz"]
    exact = HKFActivityModel(å = 4.0, Ḃ = 0.0, Kₙ = 0.0)   # symmetric by construction
    function certified(model; gap = false)
        ss = gap ?
            [SolidSolutionPhase("carbonate", [sp["Cal"], sp["Mgs"]]; model = RedlichKisterModel(a0 = 8_000.0), instances = 2)] :
            nothing
        cs = ChemicalSystem([sp[s] for s in names], comps; solid_solutions = ss)
        st = ChemicalState(cs)
        set_quantity!(st, "H2O@", 1.0u"kg")
        set_quantity!(st, "Cal", 0.025u"mol")
        set_quantity!(st, "Mgs", 0.025u"mol")
        b = Float64.(cs.SM.A) * ustrip.(us"mol", st.n)
        return equilibrate_certified(st; model, b)
    end

    _, c = certified(exact)
    @test c.optimal
    @test c.scope === :global_minimum
    @test isempty(c.scope_reasons)

    # Davies carries a salting-out term on the neutral species and a
    # mole-fraction water activity: its activities are not one gradient.
    _, c = certified(DaviesActivityModel())
    @test c.optimal
    @test c.scope === :self_consistent
    @test occursin("not the gradient", only(c.scope_reasons))

    # A concave mixing energy, even with exact activities: a KKT point, stable
    # against splitting, not a proved global minimum.
    _, c = certified(exact; gap = true)
    @test c.optimal
    @test c.scope === :kkt_point
    # One reason per instance of the phase.
    @test !isempty(c.scope_reasons) && all(r -> occursin("concave", r), c.scope_reasons)
end
