@testsection "Formulas" begin
    # Basic parsing and composition
    f = Formula("H2O")
    @test composition(f) == Dict(:H => 2, :O => 1)
    @test charge(f) == 0
    @test f[:H] == 2
    @test f[:He] == 0   # missing element returns zero
    @test length(f) == 2
    @test f == Formula("H2O")

    # String representations
    @test phreeqc(f) == "H2O"
    @test unicode(f) == "H₂O"

    # Charge handling
    na_plus = Formula("Na+")
    @test charge(na_plus) == 1
    @test composition(na_plus) == Dict(:Na => 1)
    @test phreeqc(Formula("Ca+2")) == "Ca+2"
    @test unicode(Formula("Ca+2")) == "Ca²⁺"

    # Build from composition Dict (neutral)
    comp = Dict(:Ca => 1, :C => 1, :O => 3)
    fc = Formula(comp)
    @test composition(fc) == comp
    @test fc[:Ca] == 1

    # Build from composition Dict with non-zero charge
    fc_charged = Formula(Dict(:Fe => 1), 3)
    @test charge(fc_charged) == 3
    @test fc_charged[:Fe] == 1

    # Arithmetic on formulas
    f2 = f * 2
    @test composition(f2) == Dict(:H => 4, :O => 2)
    f_half = f / 2
    @test composition(f_half) == Dict(:H => 1.0, :O => 0.5)

    # Rational division: H2O // 2 → H=2//2=1, O=1//2
    f_rat = f // 2
    @test f_rat[:H] == 1       # 2 // 2 = 1
    @test f_rat[:O] == 1 // 2  # 1 // 2
    @test f_rat[:O] isa Rational

    # Addition with AtomGroup
    f_plus_na = f + AtomGroup(:Na)
    @test composition(f_plus_na)[:Na] == 1
    f_plus_cl = f + AtomGroup(2, :Cl)
    @test composition(f_plus_cl)[:Cl] == 2

    # AtomGroup + AtomGroup → Formula
    agsum = AtomGroup(:H) + AtomGroup(2, :H)
    @test isa(agsum, Formula)
    @test composition(agsum)[:H] == 3

    # AtomGroup + Symbol
    ag_sym = AtomGroup(:O) + :H
    @test isa(ag_sym, Formula)
    @test ag_sym[:H] == 1

    # convert to Float64
    f_int = Formula(Dict(:H => 2, :O => 1))
    f_float = convert(Float64, f_int)
    vals = collect(values(composition(f_float)))
    @test all(x -> x isa Float64, vals)

    # apply: identity preserves composition
    f_applied = apply(x -> x, f_int)
    @test composition(f_applied) == composition(f_int)

    # apply: scaling by 3
    f_tripled = apply(x -> x * 3, f)
    @test f_tripled[:H] == 6
    @test f_tripled[:O] == 3

    # check_mendeleev
    @test check_mendeleev(Formula("H2O"))
    @test !check_mendeleev(Formula("Xx"))

    # hash consistency
    f_a = Formula("NaCl")
    f_b = Formula("NaCl")
    @test hash(f_a) == hash(f_b)
    @test f_a == f_b
    s = Set([f_a, f_b])
    @test length(s) == 1   # deduplicated

    # Formula with parentheses via parse_formula (used internally by Formula)
    f_ca_oh2 = Formula("Ca(OH)2")
    @test f_ca_oh2[:Ca] == 1
    @test f_ca_oh2[:O] == 2
    @test f_ca_oh2[:H] == 2
end

@testsection "Formula copies keep the charge" begin
    # Until 0.37.0 the copy was built from the composition alone, at charge 0.
    f = Formula("SO4-2")
    g = Formula(f)
    @test charge(g) == charge(f) == -2
    @test composition(g) == composition(f)
    @test expr(g) == expr(f)
end

@testsection "a negative coefficient is parenthesized, and apply keeps what it cannot map" begin
    # A negative coefficient, which removing atoms with an `AtomGroup` produces,
    # is written in parentheses, since a sign followed by a number is how a
    # charge is written. The composition and the charge are those of the
    # arithmetic: 2 − 3 hydrogens, and no charge.
    f = Formula("H2O") + AtomGroup(-3, :H)
    @test f[:H] == -1
    @test charge(f) == 0
    @test expr(f) == "H(-1)O"

    # `apply` maps every coefficient it can and keeps the others as they were:
    # `Int(1//2)` throws, so oxygen keeps its half while hydrogen is mapped.
    half = Formula("H2O") // 2
    g = apply(x -> Int(x), half)
    @test g[:H] === 1
    @test g[:O] === 1 // 2

    # Coefficients carrying a unit take their own branch; the identity leaves
    # them as they are.
    fq = Formula(OrderedDict(:H => 2.0u"mol", :O => 1.0u"mol"))
    @test composition(apply(identity, fq)) == composition(fq)
    # They are mapped: the function on the quantity, whose arithmetic keeps the
    # dimension, or on its magnitude when it takes plain numbers only. Both cases
    # used to leave them unchanged.
    @test composition(apply(x -> 2x, fq)) == OrderedDict(:H => 4.0u"mol", :O => 2.0u"mol")
    @test apply(sqrt, fq)[:H] ≈ sqrt(2.0u"mol")
    @test apply(x -> round(Int, x), fq)[:H] == 2u"mol"
end
