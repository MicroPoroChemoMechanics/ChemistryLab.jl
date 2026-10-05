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

    # The alite calibrated on the plain cement, as Section 6 of the page reports
    # it: the misfit over the five ages of Table 7, two directions determined,
    # the interaction constant and the critical degree inactive at the fit.
    opc = dw11k_run(setup, "OPC")
    fit = dw11k_alite_fit(opc)
    @test fit.rms_published ≈ 6.4 atol = 0.01
    @test fit.rms ≈ 0.43 atol = 0.01
    @test fit.θ.k₃ ≈ 7.381 rtol = 1.0e-3
    @test fit.identifiability.S[1:2] ≈ [25.1, 1.97] rtol = 1.0e-2
    @test all(<(1.0e-8), fit.identifiability.S[3:end])
end

@testsection "the CEM II/B-V of De Weerdt et al. (2011) in time with the two other gels" begin
    # Section 7 of the page: the paste with fly ash and no limestone, its C-S-H
    # under CNASH_ss and under CASH+NK, whose compound-energy activities the
    # run evaluates. Table 7 and the SEM-EDX analyses in the comments.
    for (gel, ch, ett, early, late) in (
            # portlandite at 28 and 180 days (13.7, 11.3), ettringite (7.5, 6.6),
            # the gel's Ca/Si and Al/Si at 1 day (1.7, 0.06) and 140 days (1.4, 0.13)
            ("CNASH_ss", (15.8, 13.8), (5.4, 1.8), (1.12, 0.09), (1.16, 0.107)),
            ("CASH+NK", (8.0, 3.1), (0.0, 0.0), (1.57, 0.0), (1.55, 0.0)),
        )
        r = dw11k_run(dw11k_setup(; gel), "OPC-FA")
        @test SciMLBase.successful_retcode(r.sol)
        states = dw11k_replay(r, [1, 28, 140, 180])
        c = dw11k_contents(r, states[[2, 4]])
        # Each to the digit the page prints it to.
        for k in 1:2
            @test c[k]["portlandite"] ≈ ch[k] atol = 0.05
            @test c[k]["ettringite"] ≈ ett[k] atol = 0.05
        end
        for (k, (ca_si, al_si)) in ((1, early), (3, late))
            g = dw11k_gel(states[k], gel)
            @test g.Ca_Si ≈ ca_si atol = 0.005
            @test g.Al_Si ≈ al_si atol = 0.0005
        end
    end
end
