# Sulfate attack and thaumasite (docs/src/examples/sulfate_attack.md): the CEM I
# of Lothenbach et al. (2010) titrated by Na2SO4, and the batches of Schmidt et
# al. (2008) at 8 and 20 °C. The numbers the page states are pinned.

isdefined(@__MODULE__, :sa_paste) || include(joinpath(pkgdir(ChemistryLab), "scripts", "sulfate_attack.jl"))

@testsection "durability: sulfate attack and thaumasite" begin
    @testset "a paste titrated by the solution" begin
        paste = sa_paste(; w_b = 0.5, T = 293.15)
        s = _sa_solids(paste.state)
        @test s.SO3 / s.CaO ≈ 0.045 atol = 5.0e-4
        M = ustrip(us"g/mol", Species("Na2SO4")[:M])
        volumes = 10.0 .^ range(0, 5; length = 31)
        low = sa_titration(paste, 4 / M, volumes)
        high = sa_titration(paste, 44 / M, volumes)
        @test all(r -> r.certified, low) && all(r -> r.certified, high)
        e4, e44 = sa_events(low), sa_events(high)
        @test isnan(e4.gypsum) && e4.gypsum_max == 0
        @test e44.gypsum ≈ 2154.4 atol = 0.1
        @test e4.afm_gone ≈ 1467.8 atol = 0.1
        @test e44.afm_gone ≈ 146.8 atol = 0.1
        @test e4.portlandite_gone ≈ 21544.3 atol = 0.1
        @test e44.portlandite_gone ≈ 6812.9 atol = 0.1
        @test e4.ratio_max ≈ 0.274 atol = 5.0e-4
        @test e44.ratio_max ≈ 0.571 atol = 5.0e-4
        # The measured maxima the page quotes.
        q(name) = ustrip(literature_value("Lothenbach2010sulfate", name))
        @test (q("so3_max_4"), q("so3_max_44"), q("cao_4_270d")) == (10, 20, 45)
    end

    @testset "thaumasite at 8 and 20 °C, and the solutions without it" begin
        th(b) = 100 * b.thaumasite / b.solids
        for (ls, s, t8, t20) in ((0.05, "A", 20.3, 20.3), (0.05, "B", 20.6, 20.9), (0.25, "A", 67.1, 60.6), (0.25, "B", 32.8, 31.0))
            b8 = sa_batch(; limestone = ls, level = s, T = 281.15)
            b20 = sa_batch(; limestone = ls, level = s, T = 293.15)
            @test b8.certified && b20.certified
            @test th(b8) ≈ t8 atol = 0.05
            @test th(b20) ≈ t20 atol = 0.05
        end
        meas = literature_table("Schmidt2008", "solutions_9_months")
        mval(b, s, col) = only(
            ustrip(getproperty(meas, col)[i]) for i in eachindex(meas.binder)
                if meas.temperature_C[i] == 20 && meas.binder[i] == b && meas.subsystem[i] == s
        )
        S_ratio, Ca_ratio = Float64[], Float64[]
        for (b, ls) in (("P0", 0.0), ("P5", 0.05), ("P25", 0.25)), s in ("A", "B")
            r = sa_batch(; limestone = ls, level = s, T = 293.15, thaumasite = false)
            @test r.certified
            push!(S_ratio, r.S / mval(b, s, :S))
            push!(Ca_ratio, r.Ca / mval(b, s, :Ca))
            @test 0.2 <= r.pH - mval(b, s, :pH) <= 0.35
        end
        @test minimum(S_ratio) ≈ 1.6 atol = 0.05
        @test maximum(S_ratio) ≈ 3.9 atol = 0.05
        @test all(x -> 0.5 < x < 2, Ca_ratio)
    end
end
