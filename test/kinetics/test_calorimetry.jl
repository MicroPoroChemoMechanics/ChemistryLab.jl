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

# ── _total_enthalpy ────────────────────────────────────────────────────────────

@testset "_total_enthalpy" begin

    sp1 = Species("H2O"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLVENT)
    sp1.properties[:ΔₐH⁰] = NumericFunc((T) -> H_WATER, (:T,), u"J/mol")
    sp2 = Species("CaO"; aggregate_state = AS_CRYSTAL, class = SC_COMPONENT)

    h_fns = [sp1[:ΔₐH⁰], nothing]
    n_full = [0.5, 1.0]

    H = ChemistryLab._total_enthalpy(n_full, h_fns, 298.15)
    @test isapprox(H, 0.5 * H_WATER; rtol = 1.0e-10)
    @test isapprox(H - H, 0.0; atol = 1.0e-12)

    H_none = ChemistryLab._total_enthalpy(n_full, [nothing, nothing], 298.15)
    @test iszero(H_none)

    dHdn = ForwardDiff.derivative(n -> ChemistryLab._total_enthalpy([n, 1.0], h_fns, 298.15), 0.5)
    @test isfinite(dHdn)
    @test isapprox(dHdn, H_WATER; rtol = 1.0e-10)

    dHdT = ForwardDiff.derivative(T -> ChemistryLab._total_enthalpy(n_full, h_fns, T), 298.15)
    @test isfinite(dHdT)

end

@testset "the per-species functions behave as the vector they wrap" begin

    f = NumericFunc((T) -> H_WATER, (:T,), u"J/mol")
    v = [f, nothing]
    fns = ChemistryLab._Heterogeneous(v)

    @test length(fns) == 2
    @test fns[1] === f && fns[2] === nothing
    @test eltype(fns) == eltype(v)
    @test collect(fns) == v
    # The heat terms read it exactly as they read the vector.
    @test ChemistryLab._total_enthalpy([0.5, 1.0], fns, 298.15) ==
        ChemistryLab._total_enthalpy([0.5, 1.0], v, 298.15)

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
# Under partial equilibrium the heat is `−dH/dt` over the whole composition, the
# equilibrium partition followed through `∂nₑ/∂bₑ`; before 0.24.0 it was the heat
# of the kinetic dissolution alone, and these runs refused rather than warned.

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
    ts, Q_all = cumulative_heat(sol, cal)
    @test Q_all == [u[end] for u in sol.u]
    # At the accepted steps, where the partition has just been re-speciated and
    # the heat of what the linearized partition did not predict has been added.
    # Between them the heat is interpolated without it.
    ks = unique([findmin(abs.(ts .- x))[2] for x in (0.0, 3600.0, 6 * 3600.0, 86400.0, 2 * 86400.0)])
    _, Q_ref, _ = heat_release(sol, kp; times = ts[ks])
    @info "isothermal heat under partial equilibrium" Q = Q_all[ks] Q_ref
    # The portlandite and the C-S-H precipitate, and their heat is counted.
    @test Q_ref[end] > 1000
    # At an accepted step the heat IS the enthalpy the in-run partition lost; it
    # differs from the certified replay only as that partition does, which is
    # not certified. Measured: 10.5 J of 2165 near one hour, while the assemblage
    # forms, and 2.5e-3 J at the end.
    @test Q_all[ks[end]] ≈ Q_ref[end] rtol = 1.0e-5
    @test maximum(abs.(Q_all[ks] .- Q_ref)) < 1.0e-2 * Q_ref[end]
end

@testset "an adiabatic cell under partial equilibrium conserves its enthalpy" begin
    C_vessel = 50.0
    cal = SemiAdiabaticCalorimeter(;
        Cp = C_vessel * u"J/K", T_env = 293.15u"K", heat_loss = ΔT -> zero(ΔT), T0 = 293.15u"K",
    )
    kp, sol = _partial_equilibrium_paste(cal)
    @test sol.retcode == ReturnCode.Success
    t, T = temperature_profile(sol, cal)
    @test T[end] > T[1] + 0.2
    times = [0.0, 6 * 3600.0, 86400.0, 2 * 86400.0]
    states = speciated_states(sol, kp; times)
    Tt = [sol(x)[end] for x in times]
    @test all(temperature(st) ≈ Ti * u"K" for (st, Ti) in zip(states, Tt))
    H = [ustrip(us"J", enthalpy(st)) + C_vessel * Ti for (st, Ti) in zip(states, Tt)]
    released = ustrip(us"J/K", heat_capacity(states[end])) * (Tt[end] - Tt[1])
    @info "adiabatic cell under partial equilibrium" ΔT = Tt[end] - Tt[1] drift = H .- H[1] released
    @test maximum(abs, H .- H[1]) < 1.0e-3 * released
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
