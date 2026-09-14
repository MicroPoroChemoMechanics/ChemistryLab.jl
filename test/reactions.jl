using ChemistryLab
using Test

@testsection "Reactions" begin
    @testsection "Constructor from equation" begin
        r = Reaction("H2O = H+ + OH-")
        @test reactants(r)[Species("H2O")] == 1
        @test products(r)[Species("H+")] == 1
        @test products(r)[Species("OH-")] == 1
        @test r.equal_sign == '='
    end

    @testsection "Arrow variants" begin
        # Forward arrow
        r_fwd = Reaction("H2O → H+ + OH-")
        @test r_fwd.equal_sign == '→'
        @test reactants(r_fwd)[Species("H2O")] == 1

        # Equilibrium arrow
        r_eq = Reaction("H2O ↔ H+ + OH-")
        @test r_eq.equal_sign == '↔'

        # Double arrow
        r_dbl = Reaction("H2O ⇌ H+ + OH-")
        @test r_dbl.equal_sign == '⇌'
    end

    @testsection "Getindex and properties" begin
        r = Reaction("CaCO3 = Ca+2 + CO3-2")
        @test r[Species("CaCO3")] == -1
        @test r[Species("Ca+2")] == 1
        r[:testprop] = "abc"
        @test r[:testprop] == "abc"
        @test r.testprop == "abc"
    end

    @testsection "Charge balance" begin
        # Balanced reaction: charge on left == charge on right
        r_balanced = Reaction("H2O = H+ + OH-")
        @test charge(r_balanced) == 0

        # Unbalanced: net charge ≠ 0
        r_unbalanced = Reaction("Fe+2 = Fe+3")
        @test charge(r_unbalanced) != 0
    end

    @testsection "Reaction iteration" begin
        r = Reaction("CaCO3 = Ca+2 + CO3-2")
        all_reac = collect(reactants(r))
        all_prod = collect(products(r))
        @test length(all_reac) == 1
        @test length(all_prod) == 2
        # Each item is a Pair{Species, Number}
        @test all(p -> p isa Pair, all_reac)
        @test all(p -> p isa Pair, all_prod)
    end

    @testsection "Simplify reaction" begin
        reac = Dict(Species("OH⁻") => -1, Species("H2O") => -1)
        prod = Dict(Species("H2O") => 1, Species("H⁺") => 1)
        r = Reaction(reac, prod)
        rs = simplify_reaction(r)
        @test length(reactants(rs)) == 1
        @test length(products(rs)) == 2
        @test haskey(reactants(rs), Species("OH⁻"))
        @test haskey(products(rs), Species("H⁺"))
        # H2O appears only once per side after simplification (not canceled here
        # because it has coeff -1 in reac and +1 in prod, creating net zero which
        # simplify_reaction keeps on the products side with coeff 1)
        @test !haskey(reactants(rs), Species("H2O"))
    end

    @testsection "Addition and subtraction" begin
        r1 = Reaction("H2O = H+ + OH-")
        r2 = Reaction("2H2O = H2 + 2OH⁻")
        rsum = r1 + r2
        @test reactants(rsum)[Species("H2O")] == 3

        rsub = r2 - r1
        # r2 - r1 keeps both H2O terms unsimplified: H2O appears in reactants with
        # coefficient 2 (from r2). Use simplify_reaction to obtain the net coefficient.
        rs = simplify_reaction(rsub)
        @test haskey(reactants(rs), Species("H2O")) || haskey(products(rs), Species("H2O"))
    end

    @testsection "split/merge species by stoich" begin
        s = Dict(Species("CaCl2") => -1, Species("Ca+2") => 1, Species("Cl-") => 2)
        reac, prod = ChemistryLab.split_species_by_stoich(s)
        @test haskey(reac, Species("CaCl2"))
        @test haskey(prod, Species("Ca+2"))
        merged = ChemistryLab.merge_species_by_stoich(reac, prod)
        @test haskey(merged, Species("CaCl2")) && haskey(merged, Species("Ca+2"))
    end

    @testsection "scale_stoich!" begin
        s = Dict(Species("O") => 2, Species("H") => 4)
        ChemistryLab.scale_stoich!(s)
        @test s[Species("O")] == 4 && s[Species("H")] == 8
    end
end

@testsection "a half-reaction keeps its electrons on every Julia" begin
    # THE REGRESSION. CI was green on Julia 1.13 and red on 1.12 from the same
    # commit, with `SO4-2/HS-` reported as "balancing with no electron" on 1.12
    # alone. The cause was not the chemistry and not the CPU: the `Reaction`
    # constructor stripped the `Zz` and `e` pseudo-species with `delete!`, which
    # looks a key up by `hash` and confirms with `isequal` -- and for
    # `AbstractSpecies` those two disagree.
    #
    # `isequal` compares formula, aggregate state and class; `hash` also mixes in
    # the SYMBOL. So `ELECTRON` and `Species("e")` are `==` while hashing
    # differently, and whether `delete!` reaches one through the other depends on
    # where the table puts them, hence on the hash function, hence on the Julia
    # version. Removing by symbol is version-independent.
    e = ChemistryLab.ELECTRON
    @test symbol(e) == "e-"
    # The disagreement itself, asserted so that a future fix to `isequal` or
    # `hash` shows up here as a deliberate change rather than as a surprise.
    @test e == Species("e")                       # equal...
    @test hash(e) != hash(Species("e"))           # ...and not hash-equal

    # And the property that matters: the electron survives the constructor.
    subs = build_species(datapath("cemdata18-thermofun.json"); verbose = false)
    by = Dict(symbol(s) => s for s in subs)
    r = Reaction([by["SO4-2"], by["H+"], e, by["HS-"], by["H2O@"]])
    n_e = 0
    for (sp, ν) in merge(r.reactants, r.products)
        symbol(sp) == symbol(e) && (n_e = ν)
    end
    @test n_e != 0                                # it is there at all
    @test occursin("e⁻", r.equation)
    # Eight of them, from the element and charge balance alone -- no coefficient
    # is transcribed here.
    @test occursin("8e⁻", r.equation)

    # A pseudo-species that SHOULD be stripped still is.
    r2 = Reaction([by["Cal"], by["Ca+2"], by["CO3-2"]])
    @test !any(
        String(symbol(sp)) in ("Zz", "e")
            for sp in keys(merge(r2.reactants, r2.products))
    )
end
