# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

# The pastes of De Weerdt et al. (2011) on the coupled kinetic path
# (`scripts/deweerdt2011_kinetics.jl`): the clinker phases under Parrott–Killoh,
# the fly-ash glass at the degree of reaction the authors measured, everything
# else at equilibrium. The comparison with Table 7 is the page
# `examples/ternary_kinetics.md`; this file holds a week of two pastes to what
# the laws must do whatever the data, and pins the contents the page reports at
# seven days for the limestone cement.

isdefined(@__MODULE__, :dw11k_run) || include(joinpath(pkgdir(ChemistryLab), "scripts", "deweerdt2011_kinetics.jl"))

@testsection "a ternary cement with fly ash and limestone on the coupled kinetic path" begin
    day = 86400.0
    setup = dw11k_setup()

    # The glass is a species of the system of a fly-ash paste only.
    @test any(s -> symbol(s) == "FA", dw11k_system(setup, "OPC-FA").species)
    @test !any(s -> symbol(s) == "FA", dw11k_system(setup, "OPC-L").species)

    fa = dw11k_run(setup, "OPC-FA"; days = 7)
    @test SciMLBase.successful_retcode(fa.sol)
    # The glass follows the state law it was given. From an unreacted glass the
    # law integrates to y = a + b ln(t + exp(-a/b)) percent of the fly ash, the
    # fit of Fig. 7 with c = exp(1.5) = 4.48 days in place of the printed 4.5:
    # the fit leaves 0.04 % reacted at the mixing, which the glass is not.
    a, b = dw11_value("fly_ash_fit_a"), dw11_value("fly_ash_fit_b")
    share = only(c.mass_fraction for c in setup.fly_ash.constituents if c.name == "FA")
    i = findfirst(s -> symbol(s) == "FA", fa.kp.system.species)
    k = findfirst(==(i), fa.kp.idx_kinetic)
    n0 = ustrip(us"mol", fa.kp.initial_state.n[i])
    for d in (1.0, 7.0)
        y = 100 * share * (1 - fa.sol(d * day)[size(fa.kp.Ae, 1) + k] / n0)
        @test y ≈ a + b * log(d + exp(-a / b)) rtol = 1.0e-3
    end
    # The replayed states hold the elements of the mix.
    A = Float64.(fa.kp.system.SM.A)
    b0 = A * ustrip.(us"mol", fa.kp.initial_state.n)
    states = dw11k_replay(fa, [1, 7])
    for st in states
        @test A * ustrip.(us"mol", st.n) ≈ b0 rtol = 1.0e-9 atol = 1.0e-12
    end
    c = dw11k_contents(fa, states)
    @test c[2]["portlandite"] > c[1]["portlandite"] > 0
    @test 0 < c[2]["opc_reacted"] < 100

    # The limestone cement at seven days, as the page prints it (Table 7 in
    # the comments): portlandite 17.1 (18.0), C3S 10.7 (6.1), C2S 8.4 (15.5),
    # ettringite 12.0 (9.6), clinker reacted 71.6 (69.6).
    l = dw11k_run(setup, "OPC-L"; days = 7)
    @test SciMLBase.successful_retcode(l.sol)
    c7 = only(dw11k_contents(l, dw11k_replay(l, [7])))
    for (q, v) in (("portlandite", 17.1), ("C3S", 10.7), ("C2S", 8.4), ("ettringite", 12.0), ("opc_reacted", 71.6))
        @test c7[q] ≈ v atol = 0.1
    end
end
