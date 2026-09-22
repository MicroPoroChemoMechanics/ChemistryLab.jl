@testsection "Utils" begin
    @testsection "safe_ustrip" begin
        # Quantity → strip with target unit
        @test ChemistryLab.safe_ustrip(1u"m", 5u"m") ≈ 5.0
        @test ChemistryLab.safe_ustrip(1u"cm", 1u"m") ≈ 100.0

        # Plain number → returned unchanged
        @test ChemistryLab.safe_ustrip(1u"m", 3.0) ≈ 3.0
        @test ChemistryLab.safe_ustrip(1u"K", 42) == 42
    end

    @testsection "safe_uconvert" begin
        result = ChemistryLab.safe_uconvert(us"m", 5u"cm")
        @test isapprox(ustrip(result), 0.05; rtol = 1.0e-10)

        # Plain number → returned unchanged
        @test ChemistryLab.safe_uconvert(us"m", 3.0) == 3.0
    end

    @testsection "force_uconvert" begin
        result = ChemistryLab.force_uconvert(1u"m", 100u"cm")
        @test ustrip(result) ≈ 1.0
        @test dimension(result) == dimension(u"m")

        # Plain number scaled to unit
        result2 = ChemistryLab.force_uconvert(1u"K", 300.0)
        @test ustrip(result2) ≈ 300.0
    end

    @testsection "root_type" begin
        # Vector{Int} is an alias for Array{Int,1}, so root_type returns Array
        @test ChemistryLab.root_type(Vector{Int}) == Array
        @test ChemistryLab.root_type(Dict{Symbol, Int}) == Dict
    end

    @testsection "print_title does not error" begin
        @test_nowarn ChemistryLab.print_title("Test"; style = :none)
        @test_nowarn ChemistryLab.print_title("Test"; style = :underline)
        @test_nowarn ChemistryLab.print_title("Test"; style = :box)
    end
end

@testsection "a dimensionless argument is checked, not assumed" begin
    # `_adim` replaced some thirty `Base.f(::Quantity) = f(ustrip(x))` methods.
    # Those were type piracy — they applied to every package loaded beside this
    # one — and they did not merely drop a dimension check: `ustrip` returns the
    # value in SI BASE units, so the same quantity written two ways gave two
    # answers. Both halves are asserted here.

    # THE CHECK, restored. This raised nothing at all before.
    @test_throws DimensionError ChemistryLab._adim(2u"m")
    @test_throws DimensionError ChemistryLab._adim(1u"mol/L")
    @test_throws DimensionError log(2u"m")
    @test_throws DimensionError exp(1u"K")

    # The silent answer those methods used to give, recorded so that the reason
    # is not lost: one mole per liter and one millimole per liter are different
    # quantities, and `log` of either used to return the log of its SI value.
    @test ustrip(u"mol/L" |> x -> 1x) != ustrip(u"mmol/L" |> x -> 1x)

    # WHAT MUST STILL PASS. A ratio whose units cancel IS a number.
    @test ChemistryLab._adim(2u"m/m") == 2
    @test ChemistryLab._adim(1u"1") == 1
    @test ChemistryLab._adim(0.78) === 0.78
    @test ChemistryLab._adim(3) === 3

    # Dispatch, not a branch: `Float64`, `Dual` and anything else `<: Real`
    # take the identity method, so no dimension is looked for where none can be.
    @test ChemistryLab._adim(ForwardDiff.Dual(2.0, 1.0)) isa ForwardDiff.Dual
    @test ForwardDiff.derivative(x -> log(ChemistryLab._adim(x)), 2.0) ≈ 0.5

    # The one call site that needed it, with the activity given both ways.
    r_bare = kelvin_radius(0.78; γ = 0.0728, V_m = 1.807e-5, T = 298.15)
    r_quantity = kelvin_radius(0.78u"1"; γ = 0.0728, V_m = 1.807e-5, T = 298.15)
    @test r_bare ≈ r_quantity
    # Positive, and nanometric: `ln(a_w) < 0` makes the denominator negative, so
    # the radius comes out positive, and at 78 % relative humidity the Kelvin
    # radius of water is a few nanometers — which is the order the retention
    # literature works in, and a check on the formula rather than on its sign.
    @test 1.0e-9 < r_bare < 1.0e-8
end
