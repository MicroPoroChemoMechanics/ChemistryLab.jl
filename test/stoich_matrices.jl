using ChemistryLab
using Test
using LinearAlgebra

@testsection "Stoichiometric Matrix" begin
    @testsection "Basic matrix construction" begin
        h2o = Species("H2O")
        hplus = Species("H+")
        oh = Species("OH-")
        species = [h2o, hplus, oh]

        CSM = CanonicalStoichMatrix(species)
        A, atoms = CSM.A, CSM.primaries
        @test size(A) == (3, 3)   # H, O, Zz (charge) × 3 species
        @test atoms == [:H, :O, :Zz]
        @test A[1, 1] == 2   # H2O has 2 H
        @test A[2, 1] == 1   # H2O has 1 O
    end

    @testsection "Matrix with charged species" begin
        species = [Species("Fe2+"), Species("Fe3+"), Species("e-")]
        SM = StoichMatrix(species)
        A, indep_comp, dep_comp = SM.A, SM.primaries, SM.species

        @test length(dep_comp) == 3   # Fe2+, Fe3+, e-
        @test any(s -> charge(s) != 0, dep_comp)
    end

    @testsection "Reaction conversion" begin
        na2so4 = Species("Na2(SO4)")
        na = Species("Na+")
        so4 = Species("SO4-2")
        naso4 = Species("Na(SO4)-")
        species = [na2so4, na, so4, naso4]

        SM = StoichMatrix(species)
        list_reactions = reactions(SM)

        @test length(list_reactions) > 0
        @test all(r -> r isa Reaction, list_reactions)

        # Each reaction must have non-empty reactants and products
        @test all(r -> !isempty(reactants(r)) || !isempty(products(r)), list_reactions)
    end

    @testsection "Mass-based calculations" begin
        h2o = Species("H2O")
        h2 = Species("H2")
        o2 = Species("O2")
        species = [h2o, h2, o2]

        CSM = mass_matrix(CanonicalStoichMatrix(species))
        A, atoms = CSM.A, CSM.primaries
        @test size(A) == (2, 3)   # 2 elements (H, O) × 3 species

        # Mass conservation: each column (species) sums to 1.0 (normalized mass)
        h2o_idx = findfirst(s -> s == h2o, species)
        @test sum(A[:, h2o_idx]) ≈ 1.0 atol = 1.0e-10
    end

    @testsection "CemSpecies handling" begin
        species = [CemSpecies("C3S"), CemSpecies("CH"), CemSpecies("CSH")]
        CSM = CanonicalStoichMatrix(species)
        A, atoms = CSM.A, CSM.primaries
        @test size(A, 1) ≥ 2   # At least Ca(C) and Si(S) oxide components

        # Verify atoms list contains cement oxide symbols
        @test any(x -> x in [:C, :S, :H], atoms)
    end

    @testsection "Utility functions — union_atoms" begin
        d1 = Dict(:Ca => 1, :O => 1)
        d2 = Dict(:Si => 1, :O => 2)
        atoms = union_atoms([d1, d2])
        @test :Ca in atoms
        @test :Si in atoms
        @test :O in atoms
        # No duplicates
        @test length(atoms) == length(unique(atoms))
    end

    @testsection "Utility functions — same_components" begin
        # For regular Species → should return atoms_charge
        species = [Species("CaCO3")]
        f = ChemistryLab.same_components(species)
        @test f === atoms_charge

        # For CemSpecies → should return oxides_charge
        cem_species = [CemSpecies("C3S")]
        g = ChemistryLab.same_components(cem_species)
        @test g === oxides_charge

        # Verify functions produce the expected keys
        keys_atoms = keys(atoms_charge(Species("CaCO3")))
        @test :Ca in keys_atoms && :C in keys_atoms && :O in keys_atoms

        keys_oxides = keys(oxides_charge(CemSpecies("C3S")))
        @test :C in keys_oxides   # CaO present in C3S
        @test :S in keys_oxides   # SiO2 present in C3S
    end

    @testsection "the rank is exact, and so is the null space" begin
        # THE DEFECT THIS GUARDS, measured rather than argued.
        #
        # The rank of an element-count matrix decides two booleans: whether the
        # charge row survives as an independent conservation law, and which
        # species are independent enough to be components. Both were decided by
        # `rank(A; rtol = 1e-6)` -- a LAPACK SVD -- and the elimination behind
        # the null space was fraction-free with a `÷` that TRUNCATES.
        #
        # This 5x6 integer matrix is a product of a 5x4 and a 4x6, so its rank is
        # at most 4; its smallest singular value is exactly zero and every 5x5
        # minor is singular. The fraction-free pass left a stray entry in an
        # eliminated column, counted 5 pivots, and returned a null space ONE
        # VECTOR SHORT -- a missing conservation law, reported as a fact.
        A = [
            0 -3 -14 -6 9 -2
            -2 0 3 -4 -4 -1
            -2 -2 10 -1 -8 0
            5 3 6 4 4 -5
            -8 -15 8 -2 -15 8
        ]
        @test ChemistryLab._exact_rank(A) == 4
        N = ChemistryLab._rational_nullspace(A)
        @test size(N) == (6, 2)                       # 6 columns, rank 4
        @test all(iszero, Rational{BigInt}.(A) * N)   # and it IS the null space

        # The agreement is general, not anecdotal: over a spread of shapes,
        # ranks and rational entries, the exact rank matches the numerical one,
        # the null space has the complementary dimension, and it annihilates `A`
        # exactly. A deterministic generator, so a failure is reproducible from
        # the test alone.
        seed = Ref(UInt64(0x2026_0914_0000_0007))
        nextint(lo, hi) = begin
            x = seed[]
            x ⊻= x << 13
            x ⊻= x >> 7
            x ⊻= x << 17
            seed[] = x
            return lo + Int(x % UInt64(hi - lo + 1))
        end
        for _ in 1:200
            m, n = nextint(1, 7), nextint(1, 7)
            k = nextint(1, min(m, n))
            B = [nextint(-4, 4) for _ in 1:m, _ in 1:k]
            C = [nextint(-4, 4) for _ in 1:k, _ in 1:n]
            M = B * C
            r = ChemistryLab._exact_rank(M)
            @test r == rank(M)
            Nm = ChemistryLab._rational_nullspace(M)
            @test size(Nm, 2) == n - r
            size(Nm, 2) == 0 || @test all(iszero, Rational{BigInt}.(M) * Nm)
        end
    end

    @testsection "a species no component can carry is refused, not projected" begin
        # `pinv` is a projection and a projection never fails: a species holding
        # an element that no component has comes back written over the components
        # anyway. Measured: magnetite in a system with no iron decomposed as
        # `4 H2O@ - 8 H+`, the certificate confirmed the element balance to
        # 1e-11 -- because that balance is the one the matrix stated -- and the
        # equilibrium made 2.4 mol of magnetite from a budget holding no iron.
        #
        # The membership test is a RANK COMPARISON in exact rational arithmetic
        # and not a tolerance on the least-squares residual, which matters for a
        # reason found by measurement: jennite is `(SiO2)1(CaO)1.666667(H2O)2.1`
        # in CEMDATA18, the parser keeps `5//3` on the calcium row while the
        # oxygen row sums the decimal, and a perfectly expressible `SiO2` then
        # shows a residual of 1e-6. The pozzolanic reaction below is that case,
        # and it must build.
        by = Dict(
            symbol(s) => s for s in build_species(
                    datapath("slop98-inorganic-thermofun.json"); verbose = false
                )
        )
        carbonates = split("H2O@ H+ OH- CO2@ HCO3- CO3-2 Ca+2 Mg+2 Cal")
        comps = ["H2O@", "H+", "Ca+2", "Mg+2", "CO3-2", "Zz"]
        err = try
            ChemicalSystem([by[s] for s in vcat(carbonates, "Mag")], comps)
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("Mag", sprint(showerror, err))
        @test occursin("cannot be written over", sprint(showerror, err))

        # The same list with magnesite -- which the components DO carry -- is
        # built without complaint. The refusal is about the span, not about
        # being strict.
        cs = ChemicalSystem([by[s] for s in vcat(carbonates, "Mgs")], comps)
        @test "Mgs" in String.(symbol.(cs.species))

        # AND THE CASE A TOLERANCE WOULD HAVE REFUSED. Silica IS expressible
        # over portlandite, water and jennite -- jennite carries the silicon --
        # but the decimal stoichiometry of its CEMDATA18 formula leaves a 1e-6
        # residual in the numerical decomposition. The exact test sees through
        # it; a threshold would have to be placed by hand, above this and below
        # a missing element.
        cem = build_species(datapath("cemdata18-thermofun.json"); verbose = false)
        pozz = speciation(
            cem, split("Amor-Sl Portlandite Jennite H2O@"); aggregate_state = [AS_AQUEOUS]
        )
        cs_p = ChemicalSystem(pozz, CEMDATA_PRIMARIES)
        rxn = Reaction(
            [cs_p["Amor-Sl"], cs_p["Portlandite"], cs_p["H2O@"]], [cs_p["Jennite"]]
        )
        @test occursin("SiO", rxn.equation)
    end

    @testsection "pprint does not error" begin
        h2o = Species("H2O")
        hplus = Species("H+")
        species = [h2o, hplus]
        CSM = CanonicalStoichMatrix(species)
        @test_nowarn pprint(CSM.A, CSM.primaries, CSM.species)
    end

    @testsection "pprint labels the rows the way it was asked to" begin
        # `row_labels` read `col_label`. `row_label` was assigned from `label`
        # and then never used, so asking for one accessor on the rows and
        # another on the columns labeled BOTH by the column's — silently, since
        # a wrong label prints as happily as a right one.
        #
        # Shown on species passed straight to `pprint`, because that is where it
        # is visible at all: the primaries of a `CanonicalStoichMatrix` are plain
        # `Symbol`s (`[:H, :O, :Zz]`), `symbol` of a `Symbol` throws, and the
        # surrounding `try` falls back to the raw names whatever was asked for.
        water = Species("H2O"; name = "Water", symbol = "H2O@")
        proton = Species("H+"; name = "Proton", symbol = "H+")
        rows = [water, proton]
        A = [1 0; 0 1]

        captured(; kwargs...) = mktemp() do _, io
            redirect_stdout(io) do
                pprint(A, rows, rows; kwargs...)
            end
            flush(io)
            seekstart(io)
            read(io, String)
        end

        by_name = captured(; row_label = :name, col_label = :symbol)
        @test occursin("Water", by_name)        # the rows carry `row_label`
        @test occursin("H2O@", by_name)         # the columns still carry theirs

        # The regression stated without naming a string: asking for a different
        # accessor on the rows must change the output. It did not before.
        @test by_name != captured(; row_label = :symbol, col_label = :symbol)

        # An accessor a label does not support, and a name outside the set, both
        # fall back rather than raising — the behavior the `try` always gave.
        @test_nowarn captured(; row_label = :no_such_accessor)
        @test ChemistryLab._label_accessor(:no_such_accessor) === identity
        @test ChemistryLab._label_accessor(:identity) === identity
        @test ChemistryLab._label_accessor(:symbol) === symbol
        @test ChemistryLab._label_accessor(:name) === name
        @test ChemistryLab._label_accessor(:formula) === formula

        # `label` sets both accessors at once: the same table as asking for
        # each, with the names on both axes and no symbol left.
        both = captured(; label = :name)
        @test both == captured(; row_label = :name, col_label = :name)
        @test occursin("Water", both) && occursin("Proton", both)
        @test !occursin("H2O@", both)
    end

    @testsection "a null space over a matrix given by hand keeps simple fractions exact" begin
        # Primaries that are not species send the constructor to the rational
        # elimination. A species whose column is integer gets an integer
        # coefficient, the others what the elimination gives: a fraction of
        # denominator below ten exactly, anything else as a float rounded to
        # five digits, as `stoich_coef_round` documents. Ozone carries 3/2 of
        # the O2 unit, so 3 O2 = 2 O3, written with a unit coefficient on O2.
        o2, o3 = Species("O2"), Species("O3")
        SM = StoichMatrix([1.0 1.5], [:O2], [o2, o3])
        @test size(SM.N) == (2, 1)
        @test SM.N[1, 1] == 1
        @test SM.N[2, 1] ≈ -2 / 3 rtol = 1.0e-14
        @test abs(only(SM.A * SM.N)) < 1.0e-14
        # A species carrying 1.1 of the unit: the fraction is 10/11, whose
        # denominator is past the bound, so the coefficient is the rounded float.
        nacl, salt = Species("NaCl"), Species("Na1.1Cl1.1")
        SMx = StoichMatrix([1.0 1.1], [:NaCl], [nacl, salt])
        @test SMx.N[1, 1] == 1
        @test SMx.N[2, 1] == round(-10 / 11; digits = 5)
        @test abs(only(SMx.A * SMx.N)) < 1.0e-5
    end

    @testsection "no reaction among species that are all components" begin
        # Water and carbon dioxide have independent element vectors: both are
        # components, the null space is empty, and no reaction can be written
        # between them.
        SM = StoichMatrix([Species("H2O"), Species("CO2")])
        @test size(SM.N, 2) == 0
        @test isempty(reactions(SM))
    end

    @testsection "cement species decompose over the oxides they are written in" begin
        # C3S is three lime and one silica, CH one lime and one water: the
        # cement notation itself is the reference.
        SM = StoichMatrix(
            [CemSpecies("C3S"), CemSpecies("CH")],
            [CemSpecies("C"), CemSpecies("S"), CemSpecies("H")],
        )
        coef(sp, ox) = SM.A[
            findfirst(p -> symbol(p) == ox, SM.primaries),
            findfirst(s -> symbol(s) == sp, SM.species),
        ]
        @test coef("C3S", "C") == 3 && coef("C3S", "S") == 1 && coef("C3S", "H") == 0
        @test coef("CH", "C") == 1 && coef("CH", "H") == 1 && coef("CH", "S") == 0
        @test size(SM.N, 2) == 0
    end

    @testsection "primaries chosen by a pivoted QR still decompose every species" begin
        # The basis may differ from the greedy one, and from one machine to
        # another, so what is checked holds for any basis: as many primaries as
        # the rank, and each species' elements and charge equal to those of its
        # decomposition over them.
        species = Species.(split("H2O H+ OH- CO2 HCO3- CO3-2"))
        SM = StoichMatrix(species; optimize_primaries = true)
        @test length(SM.primaries) == length(StoichMatrix(species).primaries) == 3
        for (j, s) in enumerate(SM.species), el in (:H, :O, :C, :Zz)
            @test get(atoms_charge(s), el, 0) ≈
                sum(SM.A[i, j] * get(atoms_charge(p), el, 0) for (i, p) in enumerate(SM.primaries)) atol = 1.0e-4
        end
        @test maximum(abs, Float64.(SM.A) * Float64.(SM.N)) < 1.0e-12
    end

    @testsection "a matrix with symbols: the rank, the inverse and the printing fall back" begin
        # A C-S-H of symbolic composition, as in the manual: no rank or
        # pseudo-inverse can be computed numerically on it, so the rank falls
        # back to the dimensions and the decomposition to the inverse, and the
        # colored printing, which compares each entry with zero, falls back to
        # the plain table. The decomposition is checked at a composition
        # substituted afterwards, against the one solved by hand:
        # C3S = (a − g − 3) H + CSH + (3 − a) CH.
        Sym = ChemistryLab.Symbolics
        a, g = Sym.variable(:a), Sym.variable(:g)
        CSH = CemSpecies(Dict(:C => a, :S => one(a), :H => g))
        C3S, H, CH = CemSpecies("C3S"), CemSpecies("H"), CemSpecies("CH")
        SM = StoichMatrix([C3S], [CSH, H, CH])
        @test length(SM.primaries) == 3
        @test isempty(SM.N)
        at = Dict(a => 1.7, g => 4.0)
        num(x) = Float64(Sym.value(Sym.substitute(x, at)))
        iH = findfirst(p -> symbol(p) == "H", SM.primaries)
        iCH = findfirst(p -> symbol(p) == "CH", SM.primaries)
        iCSH = only(setdiff(1:3, [iH, iCH]))
        @test num(SM.A[iCSH, 1]) ≈ 1 atol = 1.0e-12
        @test num(SM.A[iCH, 1]) ≈ 3 - 1.7 atol = 1.0e-12
        @test num(SM.A[iH, 1]) ≈ 1.7 - 4.0 - 3 atol = 1.0e-12

        out = mktemp() do _, io
            redirect_stdout(io) do
                pprint(SM)
            end
            flush(io)
            seekstart(io)
            read(io, String)
        end
        @test occursin("C3S", out) && occursin("CH", out)
    end
end

@testsection "an abstract element type is rationalized on its entries" begin
    # `Matrix{Real}` is how a `ChemicalSystem` holds integer and rational
    # coefficients together, a cement's among them. A decimal coefficient beside
    # its rational form, as CEMDATA18 gives (jennite's 1.666667 calcium is 5//3 on
    # another row), is rationalized with the package's tolerance, so the two
    # columns below are one, as they are in a `Matrix{Float64}`. Decided on the
    # element type, a `Matrix{Real}` went to the floating-point rank instead,
    # which counted two.
    @test ChemistryLab._is_rationalizable(Real[1 2 // 3; 0.5 1])
    @test ChemistryLab._exact_rank(Real[5 // 3 1.6667; 1 1]) == 1
    @test ChemistryLab._exact_rank([5 / 3 1.6667; 1 1]) == 1
    @test ChemistryLab._exact_rank(Real[1 2 // 3; 3 2]) == 1
    # An entry that is not a number of fixed value keeps the matrix out.
    @test !ChemistryLab._is_rationalizable(Real[1 ForwardDiff.Dual(1.0, 1.0)])
end
