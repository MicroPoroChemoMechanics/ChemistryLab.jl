# Leaching in zero dimensions (docs/src/examples/leaching.md): the C-S-H of
# Cemdata18 against the solubility data Berner (1992) compiled, and a CEM I
# paste leached by renewals of its solution, against the zoning of Adenot and
# Buil (1992). The numbers the page states are pinned.

isdefined(@__MODULE__, :be92_sweep) || include(joinpath(pkgdir(ChemistryLab), "scripts", "berner1992_leaching.jl"))

@testsection "durability: leaching" begin
    @testset "Berner's compilation, transcribed" begin
        m = be92_measured()
        @test length(m.Ca_Si) == 264
        @test count(i -> !ismissing(m.Ca_Si[i]) && !ismissing(m.Ca[i]), eachindex(m.Ca_Si)) == 244
        @test count(i -> !ismissing(m.Ca_Si[i]) && !ismissing(m.Si[i]), eachindex(m.Ca_Si)) == 145
        @test count(!ismissing, m.pH) == 42
    end

    sweep = be92_sweep(vcat(0.5:0.05:2.0, 2.1:0.1:3.0))
    @testset "the C-S-H alone, against the compilation" begin
        @test all(p -> p.certified, sweep)
        gel = [p for p in sweep if !p.portlandite && !p.silica]
        @test minimum(p.Ca_Si for p in gel) ≈ 0.68 atol = 0.005
        top = first(p for p in sweep if p.portlandite)
        @test top.Ca_Si ≈ 1.63 atol = 0.005
        @test top.Ca ≈ 20.3 atol = 0.05
        m = be92_measured()
        ratio(field) = [
            log10(getproperty(m, field)[i] / be92_interpolate(sweep, field, m.Ca_Si[i]))
                for i in eachindex(m.Ca_Si)
                if !ismissing(m.Ca_Si[i]) && !ismissing(getproperty(m, field)[i]) &&
                !isnan(be92_interpolate(sweep, field, m.Ca_Si[i]))
        ]
        med(x) = sort(x)[cld(length(x), 2)]
        ca, si = ratio(:Ca), ratio(:Si)
        @test length(ca) == 176
        @test count(abs.(ca) .< log10(2)) == 164
        @test 10^med(ca) ≈ 1.23 atol = 0.005
        @test length(si) == 101
        @test count(abs.(si) .< log10(2)) == 34
        @test 1 / 10^med(si) ≈ 2.3 atol = 0.05
    end

    @testset "a CEM I paste, leached" begin
        rows = be92_leaching()
        @test all(r -> r.certified, rows)
        @test length(rows) == 81
        @test rows[1].pH ≈ 13.6 atol = 0.05
        @test rows[3].pH ≈ 12.7 atol = 0.05
        @test rows[10].Ca ≈ 21.0 atol = 0.05
        @test rows[10].Ca_Si ≈ 1.63 atol = 0.005
        ev = be92_events(rows)
        @test [e.phase for e in ev] == ["portlandite", "AFm", "ettringite"]
        @test [e.water for e in ev] == [15000.0, 30000.0, 95000.0]
        @test issorted([e.Ca_Si for e in ev]; rev = true)
        @test rows[end].Ca_Si ≈ 0.705 atol = 0.005
    end
end
