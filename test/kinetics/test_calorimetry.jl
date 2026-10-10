using ForwardDiff
using JSON
using OrderedCollections

# Standard enthalpies of formation at 25 °C, read from the shipped CEMDATA18
# rather than recalled: these tests check the algebra of the heat terms, and a
# value typed from a table drifts from the database it came from.
const _CAL_H298 = let d = JSON.parsefile(datapath("cemdata18-thermofun.json"))
    Dict(
        s["symbol"] => float(only(s["sm_enthalpy"]["values"]))
            for s in d["substances"] if haskey(s, "sm_enthalpy")
    )
end
const H_WATER = _CAL_H298["H2O@"]
const H_CA2 = _CAL_H298["Ca+2"]
const H_CALCITE = _CAL_H298["Cal"]
const H_LIME = _CAL_H298["Lim"]
const H_PORTLANDITE = _CAL_H298["Portlandite"]

# ── IsothermalCalorimeter ─────────────────────────────────────────────────────

@testset "IsothermalCalorimeter" begin

    cal = IsothermalCalorimeter(298.15)
    @test cal isa IsothermalCalorimeter
    @test cal.T ≈ 298.15us"K"

    @test n_extra_states(cal) == 1

    u0 = [0.01, 0.05]
    u0_ext = extend_u0(u0, cal)
    @test length(u0_ext) == 3
    @test u0_ext[end] == 0.0
    @test u0_ext[1:2] == u0

end

# ── SemiAdiabaticCalorimeter ──────────────────────────────────────────────────

@testset "SemiAdiabaticCalorimeter" begin

    # Linear heat loss via L keyword
    cal_lin = SemiAdiabaticCalorimeter(;
        Cp = 4000.0u"J/K", T_env = 293.15u"K", L = 0.5u"W/K", T0 = 293.15u"K",
    )
    @test cal_lin isa SemiAdiabaticCalorimeter
    @test cal_lin.Cp ≈ 4000.0us"J/K"
    @test cal_lin.T_env ≈ 293.15us"K"
    @test cal_lin.T0 ≈ 293.15us"K"
    @test cal_lin.heat_loss(1.0) ≈ 0.5
    @test cal_lin.heat_loss(2.0) ≈ 1.0

    @test n_extra_states(cal_lin) == 1

    u0 = [0.01, 0.05]
    u0_ext = extend_u0(u0, cal_lin)
    @test length(u0_ext) == 3
    @test u0_ext[end] ≈ 293.15
    @test u0_ext[1:2] == u0

    # Quadratic heat loss (Lavergne et al. 2018)
    a, b = 0.48, 0.002
    cal_quad = SemiAdiabaticCalorimeter(;
        Cp = 4000.0u"J/K",
        T_env = 293.15u"K",
        heat_loss = ΔT -> a * ΔT + b * ΔT^2,
        T0 = 293.15u"K",
    )
    @test cal_quad isa SemiAdiabaticCalorimeter
    @test cal_quad.heat_loss(10.0) ≈ a * 10.0 + b * 100.0
    @test cal_quad.heat_loss(0.0) ≈ 0.0

    # Plain Real (SI) inputs
    cal_si = SemiAdiabaticCalorimeter(; Cp = 3500.0, T_env = 295.0, L = 0.3, T0 = 295.0)
    @test cal_si.Cp ≈ 3500.0us"J/K"
    @test cal_si.T_env ≈ 295.0us"K"

    # Constant heat loss
    cal_const = SemiAdiabaticCalorimeter(;
        Cp = 3000.0u"J/K",
        T_env = 293.15u"K",
        heat_loss = _ -> 2.5,
        T0 = 300.0u"K",
    )
    @test cal_const.T0 ≈ 300.0us"K"
    @test cal_const.heat_loss(0.0) ≈ 2.5
    @test cal_const.heat_loss(100.0) ≈ 2.5

    # Requires either heat_loss or L
    @test_throws ArgumentError SemiAdiabaticCalorimeter(;
        Cp = 4000.0, T_env = 293.15, T0 = 293.15,
    )

end

# ── heat_rate ─────────────────────────────────────────────────────────────────

@testset "heat_rate (constant-enthalpy reactions)" begin

    # Species with known ΔₐH⁰
    H2O = Species("H2O"; name = "Water", aggregate_state = AS_AQUEOUS, class = SC_AQSOLVENT)
    H2O.properties[:ΔₐH⁰] = NumericFunc((T) -> H_WATER, (:T,), u"J/mol")

    Ca2p = Species("Ca+2"; name = "Calcium ion", aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
    Ca2p.properties[:ΔₐH⁰] = NumericFunc((T) -> H_CA2, (:T,), u"J/mol")

    Calcite = Species("Calcite"; name = "Calcite", aggregate_state = AS_CRYSTAL, class = SC_COMPONENT)
    Calcite.properties[:ΔₐH⁰] = NumericFunc((T) -> H_CALCITE, (:T,), u"J/mol")

    # Explicit reactants/products so complete_thermo_functions! gives correct ΔᵣH⁰
    reaction = Reaction(
        OrderedDict(Calcite => 1.0),
        OrderedDict(Ca2p => 1.0);
        symbol = "calcite dissolution test",
    )
    dummy_fn = KineticFunc((T, P, t, n, lna, n0) -> 0.0, NamedTuple(), u"mol/s")
    kr = KineticReaction(reaction, dummy_fn, 1, [-1.0, 1.0])

    # Thermodynamic ΔᵣH⁰ = ΔₐH⁰(Ca²⁺) − ΔₐH⁰(Calcite) > 0 (endothermic)
    # heat_rate uses −ΔᵣH⁰: negative for endothermic (heat absorbed from calorimeter)
    ΔHr_thermo = H_CA2 - H_CALCITE
    rates = [1.0e-5]

    qdot = heat_rate([kr], rates, 298.15)
    @test isapprox(qdot, rates[1] * (-ΔHr_thermo); rtol = 1.0e-6)

    # With explicit heat_per_mol: should override stoichiometric path
    # heat_per_mol > 0 means heat generated (exothermic convention)
    kr_explicit = KineticReaction(reaction, dummy_fn, 1, [-1.0, 1.0]; heat_per_mol = 50_000.0)
    qdot_explicit = heat_rate([kr_explicit], rates, 298.15)
    @test isapprox(qdot_explicit, 1.0e-5 * 50_000.0; rtol = 1.0e-10)

    # AD smoke-test through heat_rate
    dqdot_dr = ForwardDiff.derivative(r -> heat_rate([kr], [r], 298.15), 1.0e-5)
    @test isapprox(dqdot_dr, -ΔHr_thermo; rtol = 1.0e-6)

    # A heat of reaction being fitted carries its derivative: q̇ = r ΔH, so
    # ∂q̇/∂ΔH = r. Only `Float64` had a method, and a dual raised.
    dq_dH = ForwardDiff.derivative(
        H -> heat_rate([KineticReaction(reaction, dummy_fn, 1, [-1.0, 1.0]; heat_per_mol = H)], rates, 298.15),
        50_000.0,
    )
    @test dq_dH == rates[1]

end

# ── extend_ode! for IsothermalCalorimeter ─────────────────────────────────────

@testset "extend_ode! IsothermalCalorimeter" begin

    cal = IsothermalCalorimeter(298.15)

    CaO = Species("CaO"; aggregate_state = AS_CRYSTAL, class = SC_COMPONENT)
    CaO.properties[:ΔₐH⁰] = NumericFunc((T) -> H_LIME, (:T,), u"J/mol")
    H2Osp = Species("H2O"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLVENT)
    H2Osp.properties[:ΔₐH⁰] = NumericFunc((T) -> H_WATER, (:T,), u"J/mol")
    Ca_OH_2 = Species("Ca(OH)2"; aggregate_state = AS_CRYSTAL, class = SC_COMPONENT)
    Ca_OH_2.properties[:ΔₐH⁰] = NumericFunc((T) -> H_PORTLANDITE, (:T,), u"J/mol")

    dummy_fn = KineticFunc((T, P, t, n, lna, n0) -> 0.0, NamedTuple(), u"mol/s")
    rxn = Reaction(
        OrderedDict(CaO => 1.0, H2Osp => 1.0),
        OrderedDict(Ca_OH_2 => 1.0);
        symbol = "portlandite formation",
    )
    kr = KineticReaction(rxn, dummy_fn, 1, [-1.0, -1.0, 1.0])

    du = [-0.001, 0.0]
    u = [0.01, 0.0]
    p = (kin_rxns = [kr], ϵ = 1.0e-30, rates_buf = [0.001])

    extend_ode!(du, u, p, 1, cal)
    # ΔᵣH⁰ = ΔₐH⁰(Ca(OH)₂) − ΔₐH⁰(CaO) − ΔₐH⁰(H₂O) < 0 (exothermic)
    # heat_rate uses −ΔᵣH⁰ > 0: positive qdot for exothermic (heat generated)
    ΔHr_thermo = H_PORTLANDITE - H_LIME - H_WATER
    @test isapprox(du[2], -0.001 * ΔHr_thermo; rtol = 1.0e-6)
    @test isfinite(du[2])

end

# ── SemiAdiabaticCalorimeter dT/dt ────────────────────────────────────────────

@testset "SemiAdiabaticCalorimeter dT/dt" begin

    cal = SemiAdiabaticCalorimeter(;
        Cp = 4000.0u"J/K", T_env = 293.15u"K", L = 0.5u"W/K", T0 = 293.15u"K",
    )

    CaO = Species("CaO"; aggregate_state = AS_CRYSTAL, class = SC_COMPONENT)
    CaO.properties[:ΔₐH⁰] = NumericFunc((T) -> H_LIME, (:T,), u"J/mol")
    H2Osp = Species("H2O"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLVENT)
    H2Osp.properties[:ΔₐH⁰] = NumericFunc((T) -> H_WATER, (:T,), u"J/mol")
    Ca_OH_2 = Species("Ca(OH)2"; aggregate_state = AS_CRYSTAL, class = SC_COMPONENT)
    Ca_OH_2.properties[:ΔₐH⁰] = NumericFunc((T) -> H_PORTLANDITE, (:T,), u"J/mol")

    dummy_fn = KineticFunc((T, P, t, n, lna, n0) -> 0.0, NamedTuple(), u"mol/s")
    rxn = Reaction(
        OrderedDict(CaO => 1.0, H2Osp => 1.0),
        OrderedDict(Ca_OH_2 => 1.0);
        symbol = "portlandite",
    )
    kr = KineticReaction(rxn, dummy_fn, 1, [-1.0, -1.0, 1.0])

    # ΔᵣH⁰ = ΔₐH⁰(Ca(OH)₂) − ΔₐH⁰(CaO) − ΔₐH⁰(H₂O) < 0 (exothermic)
    # heat_rate = r × (−ΔᵣH⁰) > 0: positive for exothermic → T rises
    ΔHr_thermo = H_PORTLANDITE - H_LIME - H_WATER
    r = 0.001
    qdot_expected = r * (-ΔHr_thermo)   # > 0, heat generated

    # At T = T_env: no heat loss (L * ΔT = 0); dT/dt = q̇ / Cp_total
    T_curr = 293.15
    du = [-0.001, T_curr]
    u = [0.01, T_curr]
    # p needs: kin_rxns, rates_buf, n_full, cp_fns
    p = (
        kin_rxns = [kr], ϵ = 1.0e-30, rates_buf = [0.001],
        n_full = [0.01], cp_fns = [nothing],
    )

    extend_ode!(du, u, p, 1, cal)
    # Cp_total = cal.Cp (all cp_fns are nothing)
    dTdt_expected = (qdot_expected - 0.0) / 4000.0   # positive → T rises
    @test isapprox(du[2], dTdt_expected; rtol = 1.0e-6)

    # With nonzero ΔT: heat loss L * ΔT reduces dT/dt
    T_hot = 303.15   # +10 °C above T_env
    du_hot = [-0.001, T_hot]
    u_hot = [0.01, T_hot]
    p_hot = (
        kin_rxns = [kr], ϵ = 1.0e-30, rates_buf = [0.001],
        n_full = [0.01], cp_fns = [nothing],
    )
    extend_ode!(du_hot, u_hot, p_hot, 1, cal)
    ΔT = T_hot - 293.15
    dTdt_hot = (qdot_expected - 0.5 * ΔT) / 4000.0
    @test isapprox(du_hot[2], dTdt_hot; rtol = 1.0e-6)

    # Variable Cp_total: add a Cp° function for the mineral
    cp_fn = NumericFunc((T) -> 100.0, (:T,), u"J/(mol*K)")   # 100 J/(mol·K)
    p_cp = (
        kin_rxns = [kr], ϵ = 1.0e-30, rates_buf = [0.001],
        n_full = [0.01], cp_fns = [cp_fn],
    )
    du_cp = [-0.001, T_curr]
    u_cp = [0.01, T_curr]
    extend_ode!(du_cp, u_cp, p_cp, 1, cal)
    Cp_total_expected = 4000.0 + 0.01 * 100.0   # = 4001 J/K
    dTdt_cp = (qdot_expected - 0.0) / Cp_total_expected
    @test isapprox(du_cp[2], dTdt_cp; rtol = 1.0e-6)

end

@testset "the per-species functions behave as the vector they wrap" begin

    f = NumericFunc((T) -> H_WATER, (:T,), u"J/mol")
    v = [f, nothing]
    fns = ChemistryLab._Heterogeneous(v)

    @test length(fns) == 2
    @test fns[1] === f && fns[2] === nothing
    @test eltype(fns) == eltype(v)
    @test collect(fns) == v
    # A sum over it reads exactly what a sum over the vector reads.
    term(g) = sum(n * h(; T = 298.15, unit = false) for (n, h) in zip([0.5, 1.0], g) if h !== nothing)
    @test term(fns) == term(v)

    # What it is for: the vector, heterogeneous by nature, is what SciMLBase
    # reads as badly typed parameters and warns about; the wrapper is not.
    @test ChemistryLab.SciMLBase.should_warn_paramtype((cp_fns = v,))
    @test !ChemistryLab.SciMLBase.should_warn_paramtype((cp_fns = fns,))

end

# ── The heat of a run, against the enthalpy of its states ─────────────────────
#
# A cell's first law, checked on the runs themselves. Adiabatic, the enthalpy of
# the cell, `H(t) = Σᵢ nᵢ ΔₐH⁰ᵢ(T(t)) + C_vessel T(t)`, does not change; with a
# heat loss it falls by exactly what left, `∫ φ(T − T_env) dt`. Isothermal, the
# heat the calorimeter integrates is the enthalpy drop of the certified states.
# Under partial equilibrium the state carries the change of the enthalpy of the
# cell, and the temperature is the root of its energy balance, solved with the
# partition at every evaluation: the right-hand side is a function of the state,
# and its Jacobian is exact.

const _CAL_SUBS = Dict(
    symbol(s) => s for s in build_species(datapath("cemdata18-thermofun.json"); verbose = false)
)

# Alite dissolving into ions at a rate proportional to what is left, with an
# Arrhenius factor so that the temperature of a semi-adiabatic cell feeds back.
function _alite_dissolution(cs; k = 2.0e-7, Ea = 40.0e3, T_ref = 293.15)
    rxn = Reaction(
        OrderedDict(cs["C3S"] => 1.0, cs["H2O@"] => 3.0),
        OrderedDict(cs["Ca+2"] => 3.0, cs["SiO2@"] => 1.0, cs["OH-"] => 6.0);
        symbol = "C3S dissolution",
    )
    rate(T, P, t, n, lna, n0) =
        k * n["C3S"] / n0["C3S"] * exp(-Ea / ChemistryLab.R_GAS * (1 / T - 1 / T_ref))
    rxn[:rate] = rate
    return rxn
end

function _partial_equilibrium_paste(calorimeter)
    sp = speciation(
        collect(values(_CAL_SUBS)), ["C3S", "Portlandite", "Jennite", "Amor-Sl"];
        aggregate_state = [AS_AQUEOUS], exclude_species = split("H2@ O2@ CH4@"),
    )
    cs = ChemicalSystem(sp, CEMDATA_PRIMARIES)
    st = ChemicalState(cs; T = 293.15u"K")
    set_quantity!(st, "H2O@", 0.5u"kg")
    set_quantity!(st, "C3S", 0.02u"mol")
    kp = KineticsProblem(
        cs, [_alite_dissolution(cs)], st, (0.0, 2 * 86400.0);
        calorimeter, activity_model = DaviesActivityModel(),
        equilibrium_solver = EquilibriumSolver(cs, DaviesActivityModel(), OptimaOptimizer()),
    )
    return kp, integrate(kp, KineticsSolver(; ode_solver = Rodas5P(), reltol = 1.0e-8, abstol = 1.0e-12))
end

@testset "under partial equilibrium the heat is the enthalpy the states lose" begin
    cal = IsothermalCalorimeter(293.15u"K")
    kp, sol = _partial_equilibrium_paste(cal)
    @test sol.retcode == ReturnCode.Success
    # Nothing leaves an isothermal cell but the heat the bath takes: the change
    # of its enthalpy, the state's last entry, stays zero.
    @test all(iszero(u[end]) for u in sol.u)
    ts, Q_all = cumulative_heat(sol, cal)
    ks = unique([findmin(abs.(ts .- x))[2] for x in (0.0, 3600.0, 6 * 3600.0, 86400.0, 2 * 86400.0)])
    _, Q_ref, q_ref = heat_release(sol, kp; times = ts[ks])
    # The portlandite and the C-S-H precipitate, and their heat is counted.
    @test Q_ref[end] > 1000
    # The heat is the enthalpy the paste has lost, `H₀ − H`, at the partition of
    # each instant; the certified replay solves the same partition. Measured:
    # equal to the third decimal of a joule over 2 kJ.
    @test Q_all[ks] ≈ Q_ref rtol = 1.0e-7 atol = 1.0e-6
    # And its rate, `−dH/dt` with the partition lifted along the run, is the one
    # the certified states give.
    _, qdot = heat_flow(sol, cal)
    @test qdot[ks] ≈ q_ref rtol = 1.0e-8
    # The enthalpy of a composition from the run's buffers, the kinetic amounts
    # read from `u` and the partition from `p.n_full`: at the start, where both
    # hold the initial state, its enthalpy.
    p0 = build_kinetics_params(kp)
    T0 = ustrip(us"K", temperature(kp.initial_state))
    @test system_enthalpy(p0, build_u0(kp), T0) ≈ ustrip(us"J", enthalpy(kp.initial_state)) rtol = 1.0e-12
end

@testset "an adiabatic cell under partial equilibrium conserves its enthalpy" begin
    C_vessel = 50.0
    cal = SemiAdiabaticCalorimeter(;
        Cp = C_vessel * u"J/K", T_env = 293.15u"K", heat_loss = ΔT -> zero(ΔT), T0 = 293.15u"K",
    )
    kp, sol = _partial_equilibrium_paste(cal)
    @test sol.retcode == ReturnCode.Success
    p = sol.prob.p
    t, T = temperature_profile(sol, cal)
    # The temperature is the root of the cell's balance, which Newton returns
    # within one step of `_CELL_T_TOL` (1e-9 K), and the root itself lies off
    # `T₀` by what the rounding of the partition leaves of `H − H₀`, over the
    # 2140 J/K of the cell: 1.07e-9 K in all on the Julia 1.12 job of the CI.
    @test T[1] ≈ 293.15 atol = 10 * ChemistryLab._CELL_T_TOL
    @test T[end] > T[1] + 0.2
    # The cell is closed to heat: its enthalpy, `H + C_vessel (T − T₀)`, is
    # what it was at the start, at every instant, with the certified partition
    # at the temperature the run solved. Measured: 1.5e-6 J (3e-7 J before the
    # solvent's enthalpy came from the equation of state of water) against the
    # 51 J the paste released.
    times = [0.0, 6 * 3600.0, 86400.0, 2 * 86400.0]
    states = speciated_states(sol, kp; times)
    _, Tt = temperature_profile(sol, cal; times)
    @test all(temperature(st) ≈ Ti * u"K" for (st, Ti) in zip(states, Tt))
    H = [ustrip(us"J", enthalpy(st)) + C_vessel * (Ti - 293.15) for (st, Ti) in zip(states, Tt)]
    released = C_vessel * (Tt[end] - Tt[1])
    @test maximum(abs, H .- p.H0[]) < 1.0e-7 * released
    # The heat the paste released is what warmed the vessel, exactly.
    _, Q = cumulative_heat(sol, cal)
    _, Q_ref, q_ref = heat_release(sol, kp; times)
    @test Q[[1, end]] ≈ Q_ref[[1, end]] rtol = 1.0e-7 atol = 1.0e-6
    # Its rate, `−dH/dt`, is `C_vessel dT/dt`, the paste's own heat capacity on
    # neither side: the certified states against the root of the balance.
    _, qdot = heat_flow(sol, cal)
    @test qdot[end] ≈ q_ref[end] rtol = 1.0e-6

    # The derivatives of the temperature the Jacobian holds, against a route that
    # shares nothing with them: the certified equilibrium differentiated with
    # respect to its temperature and its element amounts.
    u = Float64.(sol.u[end])
    T_end = ChemistryLab._cell_temperature(p, u)
    nb = p.n_be
    be = u[1:nb]
    h = [p.h_fns[i](; T = T_end, unit = false) for i in p.idx_equilibrium]
    function n_at(T, b)
        st = ChemicalState(p.eq_system[], p.n_full[p.idx_equilibrium] .* u"mol"; T = T * u"K", P = p.P_q[])
        eq, cert = solve_certified(p.eq_dual, (st,); b, ϵ = p.ϵ)
        @test cert.optimal
        return ustrip.(us"mol", eq.n)
    end
    n_e = n_at(T_end, be)
    nk = u[(nb + 1):(nb + p.n_nk)]
    # The heat capacity at fixed composition as the derivative of the enthalpy
    # the balance is written with: the database's own `Cp` agrees with it to
    # 2e-6 only (see above).
    cp(i, n) = n * ForwardDiff.derivative(x -> p.h_fns[i](; T = x, unit = false), T_end)
    C_eq = C_vessel + sum(cp(i, n_e[j]) for (j, i) in enumerate(p.idx_equilibrium)) +
        sum(cp(i, max(nk[j], p.ϵ)) for (j, i) in enumerate(p.idx_kinetic)) +
        h' * ForwardDiff.derivative(x -> n_at(x, be), T_end)
    # ∂T/∂ΔH = 1/C, the heat capacity of the cell at equilibrium.
    dT_dH = ForwardDiff.derivative(x -> ChemistryLab._cell_temperature(p, vcat(u[1:(end - 1)], x)), u[end])
    @info "heat capacity of the cell at equilibrium" C_eq inverse_dTdH = 1 / dT_dH
    @test 1 / dT_dH ≈ C_eq rtol = 1.0e-8
    # ∂T/∂bₑ = −(Σₑ hₑ ∂nₑ/∂bₑ)/C: the heat a change of the budget releases.
    dT_db = ForwardDiff.gradient(x -> ChemistryLab._cell_temperature(p, vcat(x, u[(nb + 1):end])), be)
    dn_db = ForwardDiff.jacobian(x -> n_at(T_end, x), be)
    @test dT_db ≈ -(h' * dn_db)' ./ C_eq rtol = 1.0e-5 atol = 1.0e-9 * maximum(abs, dT_db)

    # A state whose partition cannot be solved, a negative budget, has no
    # temperature: the right-hand side is NaN there, and an integrator rejects
    # the step rather than take a temperature the balance does not have.
    f! = build_kinetics_ode(kp)
    bad = copy(u)
    bad[1:nb] .*= -1
    du = similar(bad)
    quiet(f) = Base.CoreLogging.with_logger(f, Base.CoreLogging.NullLogger())
    quiet(() -> f!(du, bad, p, sol.t[end]))
    @test all(isnan, du)
    @test quiet(() -> ChemistryLab._cell_temperature_value(p, bad[1:nb], nk, bad[end])) === nothing

    # The temperature is solved with the partition at every evaluation, which a
    # frozen partition cannot do.
    @test_throws ArgumentError integrate(kp, KineticsSolver(; ode_solver = Rodas5P()); speciation = :frozen)
end

@testset "a cell with losses: any integrator, the same energy balance" begin
    cal = SemiAdiabaticCalorimeter(;
        Cp = 50.0u"J/K", T_env = 293.15u"K", L = 0.05u"W/K", T0 = 293.15u"K",
    )
    run_with(solver) = integrate(kp, KineticsSolver(; ode_solver = solver, reltol = 1.0e-8, abstol = 1.0e-12))
    kp, sol = _partial_equilibrium_paste(cal)
    # What left through the walls is the change of the enthalpy of the cell the
    # state carries, integrated on the dense output.
    ts = range(0.0, 2 * 86400.0; length = 4001)
    _, Tts = temperature_profile(sol, cal; times = ts)
    lost = sum((0.05 * (Tts[i] - 293.15) + 0.05 * (Tts[i + 1] - 293.15)) / 2 * (ts[i + 1] - ts[i]) for i in 1:(length(ts) - 1))
    @test -sol.u[end][end] ≈ lost rtol = 1.0e-5
    # The right-hand side is a function of the state and its Jacobian is exact:
    # a stiff Rosenbrock method, a BDF method and an explicit Runge–Kutta method
    # integrate the same trajectory. Measured: 1e-6 K apart at two days.
    tq = [0.25, 0.5, 1.0, 2.0] .* 86400.0
    T_ref = last(temperature_profile(sol, cal; times = tq))
    for solver in (FBDF(), Tsit5())
        s2 = run_with(solver)
        @test SciMLBase.successful_retcode(s2)
        @test last(temperature_profile(s2, cal; times = tq)) ≈ T_ref atol = 1.0e-5
    end
end

@testset "a stoichiometric cell loses exactly what leaves it" begin
    sp = speciation(
        collect(values(_CAL_SUBS)), ["C3S", "Portlandite", "Jennite"];
        aggregate_state = [AS_AQUEOUS], exclude_species = split("H2@ O2@ CH4@"),
    )
    cs = ChemicalSystem(sp, CEMDATA_PRIMARIES)
    rxn = Reaction(
        OrderedDict(cs["C3S"] => 1.0, cs["H2O@"] => 103 / 30),
        OrderedDict(cs["Jennite"] => 1.0, cs["Portlandite"] => 4 / 3);
        symbol = "C3S hydration",
    )
    hydration(T, P, t, n, lna, n0) =
        2.0e-7 * n["C3S"] / n0["C3S"] * exp(-40.0e3 / ChemistryLab.R_GAS * (1 / T - 1 / 293.15))
    rxn[:rate] = hydration
    st = ChemicalState(cs; T = 293.15u"K")
    set_quantity!(st, "H2O@", 0.5u"kg")
    set_quantity!(st, "C3S", 0.02u"mol")
    C_vessel, L = 50.0, 0.05
    cal = SemiAdiabaticCalorimeter(;
        Cp = C_vessel * u"J/K", T_env = 293.15u"K", L = L * u"W/K", T0 = 293.15u"K",
    )
    kp = KineticsProblem(cs, [rxn], st, (0.0, 2 * 86400.0); calorimeter = cal, equilibrium_solver = nothing)
    sol = integrate(kp, KineticsSolver(; ode_solver = Rodas5P(), reltol = 1.0e-10, abstol = 1.0e-12))
    t, T = temperature_profile(sol, cal)
    @test T == [u[end] for u in sol.u]
    # What left, integrated on the dense output.
    ts = range(0.0, 2 * 86400.0; length = 20001)
    lost = let φ = [L * (sol(x)[end] - 293.15) for x in ts]
        sum((φ[i] + φ[i + 1]) / 2 * (ts[i + 1] - ts[i]) for i in 1:(length(ts) - 1))
    end
    H(x) = ustrip(us"J", enthalpy(state_at(sol, kp, x))) + C_vessel * sol(x)[end]
    @info "stoichiometric cell" ΔTmax = maximum(T) - 293.15 lost drift = H(ts[end]) + lost - H(0.0)
    @test abs(H(ts[end]) + lost - H(0.0)) < 1.0e-3 * lost
    # The reconstruction from the vessel alone, as documented: C_vessel ΔT + ∫φ.
    tq, q = heat_flow(sol, cal)
    tQ, Q = cumulative_heat(sol, cal)
    @test length(q) == length(tq) == length(Q) == length(t)
    @test Q[end] ≈ C_vessel * (T[end] - T[1]) + lost rtol = 0.05
    # The heat rate of the paste, `−dH/dt` with the temperature term of a cell
    # that warms, is what the vessel takes and what leaves it:
    # `C_vessel dT/dt + L (T − T_env)`, on the same interpolant.
    idx = [findmin(abs.(sol.t .- x))[2] for x in (3600.0, 86400.0, 2 * 86400.0)]
    _, _, q_paste = heat_release(sol, kp; times = sol.t[idx])
    @test q_paste ≈ q[idx] rtol = 1.0e-4
end

@testset "a species without an enthalpy is refused under partial equilibrium" begin
    # Its term would drop out of `−dH/dt`, and its heat with it.
    sp = speciation(
        collect(values(_CAL_SUBS)), ["C3S", "Portlandite"];
        aggregate_state = [AS_AQUEOUS], exclude_species = split("H2@ O2@ CH4@"),
    )
    cs = ChemicalSystem(sp, CEMDATA_PRIMARIES)
    h = Any[sp[:ΔₐH⁰] for sp in cs.species]
    @test ChemistryLab._refuse_missing_enthalpy(cs, h) === nothing
    h[end] = nothing
    @test_throws ArgumentError ChemistryLab._refuse_missing_enthalpy(cs, h)
end

# ── the heat read back across an assemblage switch ───────────────────────────

isdefined(@__MODULE__, :run_ionic_hydration) ||
    include(joinpath(pkgdir(ChemistryLab), "scripts", "ionic_hydration.jl"))

@testset "the heat of a run is read back across an assemblage switch" begin
    # Cement c13 of Lerch and Ford (Lavergne et al. 2018, Table 7) at w/c = 0.4
    # and 23.9 °C, its anhydrite as the gypsum of the same sulfate. Near 4.2 h
    # hydrogarnet gives way to monosulfate, and the certified solve started from
    # the partition of the instant before stalls short of the certificate: the
    # accessor, which walks the run instant after instant, threw "the partition
    # cannot be solved at this state" where the run had passed, until the answer
    # of the interior point was offered as a start.
    t7 = literature_table("Lavergne2018", "table7_cements")
    i = findfirst(==("c13"), t7.name)
    s = parse(Float64, replace(t7.calcium_sulfates[i], "CS̅" => "")) / 100
    M(f) = ustrip(us"kg/mol", Species(f)[:M])
    g = s * (M("CaSO4") + 2 * M("H2O")) / M("CaSO4")
    binder = 1 - s + g
    f = 100 * (1 - s)
    clinker = (C3S = t7.C3S[i] / f, C2S = t7.C2S[i] / f, C3A = t7.C3A[i] / f, C4AF = t7.C4AF[i] / f)
    cal = IsothermalCalorimeter((23.9 + T_ZERO_CELSIUS)u"K")
    run = Base.CoreLogging.with_logger(Base.CoreLogging.NullLogger()) do
        run_ionic_hydration(;
            wb = 0.4 / binder, binder_mass = binder * u"kg", clinker, gypsum = g / binder,
            filler = 0.0, blaine = t7.blaine[i], calorimeter = cal, tend = 0.25 * 86400.0,
            tstops = [1 / 24, 1 / 6] .* 86400,
        )
    end
    @test SciMLBase.successful_retcode(run.sol)
    t, Q = cumulative_heat(run.sol, cal)
    @test length(t) == length(run.sol.t)
    # The same heat as the certified replay of the end of the run.
    _, Qr, _ = Base.CoreLogging.with_logger(Base.CoreLogging.NullLogger()) do
        heat_release(run.sol, run.kp; times = [0.0, t[end]])
    end
    @test Q[end] ≈ Qr[end] rtol = 1.0e-6
end

@testset "a heat-loss coefficient being fitted carries its derivative" begin
    # φ(ΔT) = L ΔT, so ∂φ/∂L = ΔT; `L` was converted to `Float64` and a dual
    # one raised. The cell's other data keep their number type as well.
    loss(L) = SemiAdiabaticCalorimeter(; Cp = 1000.0, T_env = 293.15, T0 = 293.15, L = L).heat_loss(2.0)
    @test ForwardDiff.derivative(loss, 0.5) == 2.0
    cal = SemiAdiabaticCalorimeter(; Cp = 1000.0, T_env = 293.15, T0 = ForwardDiff.Dual(300.0, 1.0), L = 0.5)
    @test ForwardDiff.partials(extend_u0([1.0, 2.0], cal)[end])[1] == 1.0
end
