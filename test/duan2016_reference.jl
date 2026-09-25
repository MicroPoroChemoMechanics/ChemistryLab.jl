# The HKF model away from 25 °C and 1 bar.
#
# Duan, Feng, Zhong, Shang & Huang, "Thermodynamic Simulation of Carbonate
# Cements-Water-Carbon Dioxide Equilibrium in Sandstone", PLoS ONE 11 (2016)
# e0167035, https://doi.org/10.1371/journal.pone.0167035, Tables 4 and 6.
#
# `test/cemdata18_reference.jl` pins every HKF equation-of-state coefficient in
# the shipped file against the table it was published in. Nothing pins what the
# model DOES with them once the temperature and pressure leave the reference
# point, and that is the half that carries a burial or an autoclave calculation.
#
# Duan et al. is a poor source in most respects -- their carbonate data for
# ferrocalcite and ankerite are estimated rather than measured, their Table 1
# and Table 3 disagree with their own text in places, and their conclusion
# prints "-2.24 mmol/L mmol/L". But their Table 4 contains one column this
# package can be held to: two equilibrium constants computed by the HKF equation
# itself, taken from Yu, Dong & Ruan (2008), at two temperatures. That is an HKF
# oracle, and it is the only one here.

@testsection "HKF away from the reference point" begin

    sp = Dict(
        symbol(s) => s for s in
            build_species(datapath("cemdata18-thermofun.json"); verbose = false)
    )
    function logK(products, reactants, T, P)
        ΔrG = sum(sp[k].ΔₐG⁰(T = T, P = P) for k in products) -
            sum(sp[k].ΔₐG⁰(T = T, P = P) for k in reactants)
        return -ΔrG / (R_GAS * T * log(10))
    end

    # Duan's Tables 4 and 6, from data/literature/Duan2016.json. Table 4's first
    # row is at 298.15 K and its second at 373.15 K; its pressure column is
    # discussed below.
    t4 = literature_table("Duan2016", "equilibrium_constants")
    calcite_row(T, P_MPa) = only(
        literature_table("Duan2016", "calcite"; T = T * u"K", P = P_MPa * 1.0e6u"Pa").log_K
    )
    K(args...) = 10.0^logK(args...)

    # K3 is the second dissociation of carbonic acid, K4 water autoprotolysis,
    # in Duan's numbering.
    k3(T, P) = K(("H+", "CO3-2"), ("HCO3-",), T, P)
    k4(T, P) = K(("H+", "OH-"), ("H2O@",), T, P)

    @testset "Table 4 at the reference point" begin
        # 298.15 K, 1 bar. Both to 0.03 %, which is inside the two figures Duan
        # quote.
        @test k3(298.15, 1.0e5) ≈ t4.K3_HKF[1] rtol = 1.0e-3
        @test k4(298.15, 1.0e5) ≈ t4.K4_HKF[1] rtol = 1.0e-3
    end

    @testset "Table 4 at 373.15 K, and the unit its pressure column is in" begin
        # The column is headed "P(Pa)" and reads 1.0 and 1000.0. Neither can be
        # pascals: 1 Pa is a hard vacuum and no aqueous constant is tabulated
        # there. Bar fits, and K3 settles it. At 373.15 K the published value is
        # 2.42e-10, and this model gives
        #
        #     1 bar     8.240e-11   -66 %
        #     1000 bar  2.412e-10   -0.3 %
        #
        # so the row is at 1000 bar, and reproducing a factor of 2.9 in pressure
        # to three parts in a thousand is a real check of the volume terms.
        @test k3(373.15, 1.0e8) ≈ t4.K3_HKF[2] rtol = 5.0e-3
        @test k3(373.15, 1.0e5) / k3(373.15, 1.0e8) < 0.4      # the factor itself

        # K4 IN THE SAME ROW DOES NOT AGREE AT THAT PRESSURE. It matches at
        # 1 bar (5.50e-13 against 5.38e-13, 2.3 %) and is 14.7 % out at 1000
        # bar. Water autoprotolysis at 100 °C and 1 bar is pKw = 12.26, so the
        # 1 bar value is the physical one and the row mixes two pressures.
        # Pinned as found: the tolerance on the second is the disagreement.
        @test k4(373.15, 1.0e5) ≈ t4.K4_HKF[2] rtol = 0.03
        @test !isapprox(k4(373.15, 1.0e8), t4.K4_HKF[2]; rtol = 0.05)
    end

    @testset "calcite: the textbook value, and a method that is not HKF" begin
        # Table 6 is calcite dissolving to Ca2+ + CO3^2-, over 298-478 K and
        # 0.1-70 MPa. Duan compute it with their own Gibbs-energy integration
        # and SRK volumes, NOT with HKF, so this is a comparison of two methods
        # and a disagreement is not by itself a defect.
        #
        # At the reference point there is an independent arbiter, and it favors
        # this package: the accepted log Ksp of calcite at 25 °C and 1 bar is
        # -8.48, which is what comes out here. Duan print -8.53.
        @test logK(("Ca+2", "CO3-2"), ("Cal",), 298.15, 1.0e5) ≈
            literature_value("PlummerBusenberg1982", "log_K_calcite_25C") atol = 0.01

        # The two signs both papers agree on, and which any burial calculation
        # turns on: heating dissolves less, compressing dissolves more.
        hot = logK(("Ca+2", "CO3-2"), ("Cal",), 478.15, 1.5e7)
        cold = logK(("Ca+2", "CO3-2"), ("Cal",), 301.15, 1.5e7)
        @test hot < cold - 2.0                                  # Duan: -9.69 vs -8.53
        squeezed = logK(("Ca+2", "CO3-2"), ("Cal",), 301.15, 7.0e7)
        @test squeezed > cold + 0.1                             # Duan: -7.82 vs -8.53

        # How far apart the two methods run. At 301 K they are within 0.1 log
        # unit; by 478 K they are 1.4 apart, this package giving the lower
        # solubility. Worth recording because Duan's own Table 4 puts their
        # method within 3 % of HKF on K3 and K4 -- 0.01 log units -- so a
        # disagreement of 1.4 on calcite is far outside what they claim for it,
        # and the estimated carbonate data of their Table 3 are the likeliest
        # reason. Their ferrocalcite row carries a heat-capacity coefficient of
        # +2.09e6 where every other carbonate in the table has zero or a large
        # negative.
        @test abs(cold - calcite_row(301.15, 15.0)) < 0.15
        @test hot - calcite_row(478.15, 15.0) < -1.0
    end
end
