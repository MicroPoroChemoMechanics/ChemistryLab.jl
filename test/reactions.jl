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

    @testsection "a small coefficient stays in the reaction" begin
        # A coefficient below 1e-3 was rounded to zero and its species dropped.
        # (One within 1e-3 of a simple fraction is still read as that fraction,
        # as the formulas of the databases mean it: 4.0004 is 4.)
        r = Reaction("Ca2Si3Na0.0004O11H6 + 4.199H+ = 0.0004Na+ + 2Ca+2 + 3SiO2 + 5.0993H2O")
        @test products(r)[Species("Na+")] == 0.0004
        @test reactants(r)[Species("H+")] == 4.199
        @test occursin("0.0004", string(r))
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
        # The reaction itself yields its species with signed coefficients, as
        # pairs on both sides: the reactants came as tuples, so that `collect`
        # gave a `Vector{Any}`.
        terms = collect(r)
        @test terms isa Vector{Pair{Species{Int}, Int}}
        @test first.(terms) == collect(keys(r))
        @test last.(terms) == [-1, 1, 1] == collect(values(r))
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
        # The smallest integers in the same ratios. Until 0.37.0 the
        # coefficients were multiplied by their divisor: {2, 4} gave {4, 8}.
        s = Dict(Species("O") => 2, Species("H") => 4)
        ChemistryLab.scale_stoich!(s)
        @test s[Species("O")] == 1 && s[Species("H")] == 2
        r = Dict{Species, Rational{Int}}(Species("Ca") => -1, Species("O") => 1 // 2, Species("H") => 3 // 2)
        ChemistryLab.scale_stoich!(r)
        @test r[Species("Ca")] == -2 && r[Species("O")] == 1 && r[Species("H")] == 3
        f = Dict(Species("O") => 0.5, Species("H") => 1.0)
        @test ChemistryLab.scale_stoich!(f) == Dict(Species("O") => 0.5, Species("H") => 1.0)
    end
end

@testsection "a half-reaction keeps its electrons on every Julia" begin
    # THE REGRESSION, AND ITS ROOT CAUSE, NOW REMOVED. CI was green on Julia 1.13
    # and red on 1.12 from the same commit, with `SO4-2/HS-` reported as
    # "balancing with no electron" on 1.12 alone. The cause was not the chemistry
    # and not the CPU: the `Reaction` constructor stripped the `Zz` and `e`
    # pseudo-species with `delete!`, which looks a key up by `hash` and confirms
    # with `isequal` -- and for `AbstractSpecies` those two disagreed. `isequal`
    # compared formula, aggregate state and class; `hash` also mixed in the
    # SYMBOL. So `ELECTRON` and `Species("e")` were `==` while hashing
    # differently, and whether `delete!` reached one through the other depended
    # on where the table put them, hence on the hash function, hence on the Julia
    # version.
    #
    # This testset used to assert that disagreement, "so that a future fix shows
    # up here as a deliberate change rather than as a surprise". This is that
    # change: `hash` and `isequal` now rest on the same four things, so the two
    # agree and the version-dependence has nowhere to come from. Removing by
    # symbol remains what the constructor does, and stays correct either way.
    e = ChemistryLab.ELECTRON
    @test symbol(e) == "e-"
    # Agreement, in both directions, on the very pair that broke.
    @test e == Species("e")
    @test hash(e) == hash(Species("e"))
    @test Dict(e => 1)[Species("e")] == 1         # what `delete!` needs

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

@testsection "Reactions: the less common constructors and accessors" begin
    h2, o2, h2o = Species("H2"), Species("O2"), Species("H2O")

    @testsection "absent species, absent properties, and fixed fields" begin
        r = Reaction("2H2 + O2 = 2H2O")
        # A species that takes no part has a zero coefficient.
        @test r[Species("CO2")] == 0
        # A key that is neither a property nor a thermodynamic function is absent.
        @test !haskey(r, :not_a_property)
        # The structural fields are set by the constructors only.
        @test_throws "Cannot modify field 'equation'" setproperty!(r, :equation, "H2O = H2O")
        @test r.equation == "2H2 + O2 = 2H2O"
    end

    @testsection "a reaction is named after a participant when asked to" begin
        # `symbol = nothing` names it after its last reactant, after its first
        # product when it has no reactant, and leaves it unnamed when it has
        # neither: from an equation and from two dictionaries alike.
        @test symbol(Reaction("CaCO3 + H+ = Ca+2 + HCO3-"; symbol = nothing)) == "H+"
        @test symbol(Reaction("∅ = 2H2O + O2"; symbol = nothing)) == "H2O"
        @test symbol(Reaction("∅ = ∅"; symbol = nothing)) == ""
        none = OrderedDict{Species{Int}, Int}()
        @test symbol(Reaction(OrderedDict(h2 => 2, o2 => 1), OrderedDict(h2o => 2); symbol = nothing)) == "O2"
        @test symbol(Reaction(none, OrderedDict(h2o => 2, o2 => 1); symbol = nothing)) == "H2O"
        @test symbol(Reaction(none, none; symbol = nothing)) == ""
    end

    @testsection "one electron is written without a coefficient" begin
        # Fe²⁺ → Fe³⁺ gives up one electron. The balancing term is `e⁻`, not
        # `1e⁻`, in the plain and in the colored equation alike.
        r = Reaction(OrderedDict(Species("Fe+2") => 1), OrderedDict(Species("Fe+3") => 1))
        @test r.equation == "Fe²⁺ = Fe³⁺ + e⁻"
        @test endswith(colored(r), " + e⁻")
    end

    @testsection "a species on both sides with the same coefficient cancels" begin
        # The sodium ion is a spectator of the hydration of CO2: it leaves both
        # sides, and the rest is the reaction without it.
        r = simplify_reaction(Reaction("CO2 + H2O + Na+ = H2CO3 + Na+"))
        @test !haskey(reactants(r), Species("Na+")) && !haskey(products(r), Species("Na+"))
        @test r[Species("Na+")] == 0
        @test r[Species("CO2")] == -1 && r[Species("H2O")] == -1 && r[Species("H2CO3")] == 1
    end

    @testsection "balanced coefficients scaled to the smallest integers" begin
        # Balanced on water as the dependent species, water splitting reads
        # H2O = H2 + ½O2; scaled, 2H2O = 2H2 + O2.
        unscaled = Reaction([h2o, h2, o2])
        @test unscaled[h2o] == -1 && unscaled[h2] == 1 && unscaled[o2] == 1 // 2
        r = Reaction([h2o, h2, o2]; auto_scale = true)
        @test r[h2o] == -2 && r[h2] == 2 && r[o2] == 1
    end

    @testsection "every participant on one side keeps its signed coefficient" begin
        # Lime slaking, CaO + H2O = Ca(OH)2. With `side = :products` every
        # species is written on the right, the reactants with negative
        # coefficients, and the signed coefficient of each is unchanged.
        lime, water, port = Species("CaO"), Species("H2O"), Species("Ca(OH)2")
        ref = Reaction([lime, water], [port])
        @test ref[lime] == -1 && ref[water] == -1 && ref[port] == 1
        right = Reaction([lime, water], [port]; side = :products)
        @test isempty(reactants(right))
        @test all(right[s] == ref[s] for s in (lime, water, port))
        # The same from the equation.
        from_eq = Reaction("CaO + H2O = Ca(OH)2"; side = :products)
        @test isempty(reactants(from_eq))
        @test from_eq[lime] == -1 && from_eq[water] == -1 && from_eq[port] == 1
    end

    @testsection "parameters given on a reaction set its reference state" begin
        # A reaction carrying its own parameters, as a database reader stores
        # them, takes its reference temperature and pressure from that list when
        # it states none. The value of log K is made up: only where it lands
        # matters.
        r = Reaction("CaCO3 = Ca+2 + CO3-2")
        r[:thermo_params] = [:logKr => 1.25, :T => 298.15u"K", :P => 1.0e5u"Pa"]
        ChemistryLab.complete_thermo_functions!(r)
        @test r.Tref == 298.15u"K"
        @test r.Pref == 1.0e5u"Pa"
        @test r[:logKr_Tref] == 1.25
        @test r[:logKr](T = 350.0) ≈ 1.25
        @test !haskey(properties(r), :thermo_params)
    end
end

@testsection "the dictionary constructor leaves its arguments as they were" begin
    # The pseudo-species `e` and `Zz` are dropped from the reaction, not from the
    # caller's dictionaries, which a `filter!` on them used to empty of their
    # electrons.
    e = Species("e"; symbol = "e")
    reac = OrderedDict(Species("Fe+3") => 1, e => 1)
    prod = OrderedDict(Species("Fe+2") => 1)
    r = Reaction(reac, prod)
    @test length(reac) == 2 && haskey(reac, e)
    @test !any(String(symbol(sp)) in ("Zz", "e") for sp in keys(r.reactants))
end
