using JSON
using ChemistryLab.DataFrames: metadata

# The readers of thermo datasets of The Geochemist's Workbench and of data0 files
# of EQ3/6, on files written during the test from `llnl.dat`, which PHREEQC
# distributes: the same LLNL data in the two other formats. No dataset of either
# is stored or downloaded. The log K of `llnl.dat`, PHREEQC's analytic expression,
# is written exactly as the polynomial of GWB format "jan19", and as its values
# on the grid of a data0 file, which the reader interpolates as EQPT does: the
# equilibria at 25 and 60 °C, grid temperatures, are those of `llnl.dat`, and
# are compared with PHREEQC (`test/reference/phreeqc_databases.json`).

const _FORMATS_ELEMENTS = Set([:H, :O, :Na, :K, :Ca, :Mg, :Cl, :S, :C])

# The coefficients a-f of GWB's polynomial for a log K of formation, a
# combination of PHREEQC constants: a + b(T-Tr) + c(T²-Tr²) + d(1/T-1/Tr) +
# e(1/T²-1/Tr²) + f ln(T/Tr).
function _gwb_coefficients(formation)
    Tr = T_STANDARD
    acc = zeros(6)
    for (w, k) in formation.terms
        if k.analytic !== nothing
            A1, A2, A3, A4, A5, A6 = k.analytic
            acc .+= w .* [A1 + A2 * Tr + A3 / Tr + A4 * log10(Tr) + A5 / Tr^2 + A6 * Tr^2, A2, A6, A3, A5, A4 / log(10)]
        else
            c = 1000 * k.delta_h / (R_GAS * log(10))
            acc .+= w .* [k.log_k, 0, 0, -c, 0, 0]
        end
    end
    return acc
end

# A species as the sum of the master species: its dissociation reaction.
function _dissociation(sp, masters_of)
    rest = Dict{Symbol, Float64}(sp.atoms)
    terms = Tuple{Float64, String}[]
    z = sp.charge
    for (e, n) in sort!(collect(rest); by = first)
        e in (:H, :O) && continue
        m = masters_of[e]
        k = n / m.atoms[e]
        push!(terms, (k, m.name))
        z -= k * m.charge
        for (e2, n2) in m.atoms
            rest[e2] = get(rest, e2, 0.0) - k * n2
        end
        rest[e] = 0.0
    end
    nH = get(rest, :H, 0.0) - 2 * get(rest, :O, 0.0)
    nO = get(rest, :O, 0.0)
    iszero(nO) || push!(terms, (nO, "H2O"))
    iszero(nH) || push!(terms, (nH, "H+"))
    return terms
end

function _write_formats(db, dir)
    candidates = [r for r in eachrow(db) if r.aggregate_state == AS_AQUEOUS && all(in(_FORMATS_ELEMENTS), keys(r.atoms))]
    masters = [r for r in candidates if isempty(r.formation.terms)]
    masters_of = Dict{Symbol, Any}()
    for m in masters, e in keys(m.atoms)
        e in (:H, :O) || (masters_of[e] = m)
    end
    # The species of the oxidation states of the masters (sulfate, carbonate):
    # those whose dissociation into the masters balances the charge.
    balanced(sp) = abs(sum(c * (n == "H+" ? 1.0 : n == "H2O" ? 0.0 : only(m.charge for m in masters if m.name == n)) for (c, n) in _dissociation(sp, masters_of); init = 0.0) - sp.charge) < 1.0e-9
    master_names = Set(m.name for m in masters)
    others = [sp for sp in candidates if !(sp.name in master_names) && balanced(sp)]
    species = vcat(masters, others)
    phases = [r for r in eachrow(db) if r.name in ("Calcite", "Gypsum")]
    # The molar volumes of the two minerals, which llnl.dat does not give, from
    # phreeqc.dat; their molar masses from the formula.
    _, pdat, _ = read_phreeqc_database(datapath("phreeqc.dat"))
    vm(sp) = ustrip(us"cm^3/mol", only(pdat[pdat.name .== sp.name, :molar_volume]) * u"m^3/mol")
    mw(sp) = ustrip(us"g/mol", Species(String(sp.formula))[:M])
    T8 = [0.01, 25.0, 60.0, 100.0, 150.0, 200.0, 250.0, 300.0]
    logK(sp, tc) = -first(ChemistryLab._log10K(sp.formation, tc + T_ZERO_CELSIUS))
    num(x) = string(round(x; digits = 12) + 0.0)   # + 0.0: no -0.0

    gwb = joinpath(dir, "llnl.tdat")
    open(gwb, "w") do io
        println(io, "dataset of thermodynamic data for gwb programs\ndataset format: jan26\nactivity model: b-dot")
        println(io, "* temperatures")
        println(io, join(T8[1:4], " "), "\n", join(T8[5:8], " "))
        for name in ("pressures", "debye huckel a (adh)", "debye huckel b (bdh)")
            println(io, "* ", name, "\n1 1 1 1\n1 1 1 1")
        end
        println(io, "   $(length(_FORMATS_ELEMENTS)) elements")
        for e in sort!(collect(_FORMATS_ELEMENTS))
            println(io, "Element$(e) ($(e)) mole wt.= 1.0 g")
        end
        println(io, "-end-\n   $(length(masters)) basis species")
        for m in masters
            println(io, m.name, "\n  charge= ", Int(m.charge), "  ion size= 4.0 A  mole wt.= 1.0 g")
            println(io, "  $(length(m.atoms)) elements in species")
            println(io, "  ", join(("$(num(n)) $(e)" for (e, n) in m.atoms), "  "))
        end
        println(io, "-end-\n   0 redox couples\n-end-\n   $(length(others)) aqueous species")
        for sp in others
            terms = _dissociation(sp, masters_of)
            a = -_gwb_coefficients(sp.formation)
            println(io, sp.name, "\n  charge= ", Int(sp.charge), "  ion size= 4.0 A  mole wt.= 1.0 g")
            println(io, "  $(length(terms)) species in reaction")
            println(io, "  ", join(("$(num(c)) $(n)" for (c, n) in terms), "  "))
            println(io, "  a= $(num(a[1]))  b= $(num(a[2]))  c= $(num(a[3]))\n  d= $(num(a[4]))  e= $(num(a[5]))  f= $(num(a[6]))")
        end
        println(io, "-end-\n   0 free electron\n-end-\n   $(length(phases)) minerals")
        for sp in phases
            terms = _dissociation(sp, masters_of)
            a = -_gwb_coefficients(sp.formation)
            println(io, sp.name, "\n  formula= ", sp.formula, "\n  mole vol.= ", round(vm(sp); digits = 3), " cc  mole wt.= ", round(mw(sp); digits = 4), " g")
            println(io, "  $(length(terms)) species in reaction")
            println(io, "  ", join(("$(num(c)) $(n)" for (c, n) in terms), "  "))
            println(io, "  a= $(num(a[1]))  b= $(num(a[2]))  c= $(num(a[3]))\n  d= $(num(a[4]))  e= $(num(a[5]))  f= $(num(a[6]))")
        end
        println(io, "-end-\n   0 solid solutions\n-end-\n   0 gases\n-end-\n   0 oxides\n-end-")
    end

    eq36 = joinpath(dir, "data0.tst")
    open(eq36, "w") do io
        bar = "+" * "-"^68
        println(io, "data0.tst.R0\nwritten during a test from llnl.dat\n", bar)
        println(io, "data0 parameters\ntemperature limits\n   0.0100   300.0000\ntemperatures")
        println(io, join(T8[1:4], "  "), "\n", join(T8[5:8], "  "))
        println(io, "pressures\n1 1 1 1\n1 1 1 1\n", bar, "\nbdot parameters")
        for sp in species
            println(io, rpad(sp.name, 24), "   4.0000    0")
        end
        println(io, bar, "\nelements")
        for e in sort!(collect(_FORMATS_ELEMENTS))
            println(io, rpad(lowercase(string(e)), 8), "  1.00000")
        end
        println(io, bar)
        block(sp, kind, terms) = begin
            println(io, rpad(sp.name, 24), "\n    keys   = $kind")
            sp.aggregate_state == AS_AQUEOUS ? println(io, "     charge  =   $(sp.charge)") : println(io, "     V0PrTr =   $(round(vm(sp); digits = 3)) cm**3/mol")
            println(io, "     $(length(sp.atoms)) chemical elements =")
            println(io, "      ", join(("$(num(n)) $(lowercase(string(e)))" for (e, n) in sp.atoms), "    "))
            if !isempty(terms)
                println(io, "     $(length(terms) + 1) species in reaction =")
                all = vcat([(-1.0, sp.name)], terms)
                for k in 1:2:length(all)
                    println(io, join((" " * lpad(string(round(c; digits = 4)), 10) * "  " * rpad(n, 24) for (c, n) in all[k:min(k + 1, end)]), ""))
                end
                println(io, "*\n     log k grid (0-25-60-100/150-200-250-300 C) =")
                v = [logK(sp, t) for t in T8]
                println(io, "     ", join(num.(v[1:4]), "  "), "\n     ", join(num.(v[5:8]), "  "))
            end
            println(io, bar)
        end
        println(io, "basis species\n", bar)
        foreach(m -> block(m, "basis", Tuple{Float64, String}[]), masters)
        println(io, "auxiliary basis species\n", bar, "\naqueous species\n", bar)
        foreach(sp -> block(sp, "aqueous", _dissociation(sp, masters_of)), others)
        println(io, "solids\n", bar)
        foreach(sp -> block(sp, "solid", _dissociation(sp, masters_of)), phases)
        println(io, "liquids\n", bar, "\ngases\n", bar, "\nsolid solutions\n", bar, "\nreferences\n", bar, "\nstop.")
    end
    return gwb, eq36
end

@testsection "Thermo datasets of GWB and data0 files of EQ3/6" begin
    oracle = JSON.parsefile(joinpath(@__DIR__, "reference", "phreeqc_databases.json"))["databases"]["llnl.dat"]
    _, llnl, _ = read_phreeqc_database(datapath("llnl.dat"))
    model = database_activity_model(llnl)
    mktempdir() do dir
        gwb, eq36 = _write_formats(llnl, dir)
        # The excerpts the manual shows: lines of `llnl.dat`, and the calcite this
        # writer gives in the two other formats.
        page = read(joinpath(pkgdir(ChemistryLab), "docs", "src", "manual", "importing_databases.md"), String)
        excerpts = [m.captures[1] for m in eachmatch(r"```text\n(.*?)```"s, page)]
        @test length(excerpts) == 3
        @test all(in(Set(readlines(datapath("llnl.dat")))), filter(!isempty, split(chomp(excerpts[1]), "\n")))
        g = readlines(gwb)
        i = findfirst(==("Calcite"), g)
        lines(t) = rstrip.(split(chomp(t), "\n"))
        @test lines(excerpts[2]) == rstrip.(g[i:(i + 6)])
        e = readlines(eq36)
        i = findfirst(l -> startswith(l, "Calcite "), e)
        @test lines(excerpts[3]) == rstrip.(e[i:findnext(l -> startswith(l, "+---"), e, i)])
        # `import_database` chooses the reader by the name of the file, then by
        # its first lines.
        fmt = ChemistryLab._database_format
        @test fmt(gwb) == :gwb && fmt(eq36) == :eq36 && fmt(datapath("llnl.dat")) == :phreeqc
        @test fmt("any.json") == :thermofun && fmt("any.yml") == :reaktoro && fmt("any.yaml") == :reaktoro
        for (from, to, f) in ((gwb, "gwb.txt", :gwb), (eq36, "eq36.txt", :eq36), (datapath("llnl.dat"), "llnl.txt", :phreeqc))
            cp(from, joinpath(dir, to))
            @test fmt(joinpath(dir, to)) == f
        end
        write(joinpath(dir, "notes.txt"), "nothing a reader knows\n")
        @test_throws "is not recognized" fmt(joinpath(dir, "notes.txt"))
        @test metadata(import_database(joinpath(dir, "gwb.txt"))[2], "format") == "gwb"
        @test metadata(import_database(joinpath(dir, "eq36.txt"))[2], "format") == "eq36"
        @test metadata(import_database(joinpath(dir, "llnl.txt"))[2], "format") == "phreeqc"
        @test metadata(import_database(joinpath(dir, "gwb.txt"); format = :gwb)[2], "format") == "gwb"
        @test_throws "format is :auto" import_database(gwb; format = :csv)
        for (label, (_, db, _)) in (("GWB", read_gwb_database(gwb)), ("EQ3/6", read_eq36_database(eq36)))
            @testset "$label" begin
                @test isempty(filter(n -> !occursin("not read", n) || occursin("rests on", n) || occursin("balance", n), metadata(db, "notes")))
                # The log K of formation of every species is that of llnl.dat:
                # everywhere for the polynomial of GWB, at the grid temperatures
                # for the interpolating polynomials of a data0 file.
                # By symbol: llnl.dat names a solute and a phase alike (MgSO4).
                ref = Dict(r.symbol => r for r in eachrow(llnl))
                for sp in eachrow(db), tc in (25.0, 60.0, 150.0, 300.0)
                    T = tc + T_ZERO_CELSIUS
                    mine = first(ChemistryLab._log10K(sp.formation, T))
                    theirs = first(ChemistryLab._log10K(ref[sp.symbol].formation, T))
                    @test mine ≈ theirs atol = 1.0e-8
                end
                # Calcite and gypsum in water, against PHREEQC with llnl.dat.
                cs = ChemicalSystem(build_species(db))
                for T in (25.0, 60.0)
                    o = oracle["minerals"][string(T)]
                    st = ChemicalState(cs; T = (T_ZERO_CELSIUS + T) * u"K", P = 1.0u"Constants.atm")
                    set_quantity!(st, "H2O@", 1.0u"kg")
                    set_quantity!(st, "Calcite", 10.0u"mol")
                    set_quantity!(st, "Gypsum", 10.0u"mol")
                    eq = equilibrate(st; model)
                    for m in ("Calcite", "Gypsum")
                        @test 10 - ustrip(us"mol", moles(eq, m)) ≈ o["dissolved"][m] rtol = 2.0e-5
                    end
                    @test pH(eq, model) ≈ o["pH"] atol = 1.0e-4
                end
            end
        end
    end
end

@testsection "Thermo datasets of GWB and data0 files of EQ3/6: what is not read" begin
    mktempdir() do dir
        # A thermo dataset with a stray line, a section not read, a range in
        # Celsius, a species resting on another the dataset does not define, and
        # a basis species of a pseudo-element.
        gwb = joinpath(dir, "faults.tdat")
        write(
            gwb, """
            dataset of thermodynamic data for gwb programs
            dataset format: jan26
            * temperatures
            0.01 25.0 60.0 100.0
            150.0 200.0 250.0 300.0
               3 elements
            Calcium (Ca) mole wt.= 40.08 g
            Hydrogen (H) mole wt.= 1.008 g
            Oxygen (O) mole wt.= 16.0 g
            -end-
            a stray line
               4 basis species
            H+
              charge= 1  ion size= 9.0 A  mole wt.= 1.0 g
              1 elements in species
              1.000 H
            H2O
              charge= 0  ion size= 0.0 A  mole wt.= 18.0 g
              2 elements in species
              2.000 H  1.000 O
            Ca+2
              charge= 2  ion size= 6.0 A  mole wt.= 40.0 g
              1 elements in species
              1.000 Ca
            Qq+
              charge= 1  ion size= 4.0 A  mole wt.= 1.0 g
              1 elements in species
              1.000 Qq
            -end-
               2 aqueous species
            CaOH+
              charge= 1  ion size= 4.0 A  mole wt.= 57.0 g
              3 species in reaction
              1.000 Ca+2  1.000 H2O  -1.000 H+
              a= 12.8  b= 0  c= 0
              d= 0  e= 0  f= 0
              TminC= 0.0  TmaxC= 300.0
            Undef1+
              charge= 1  ion size= 4.0 A  mole wt.= 1.0 g
              1 species in reaction
              1.000 Xx+
              a= 1.0  b= 0  c= 0
              d= 0  e= 0  f= 0
            -end-
               1 oxides
            CaO
            -end-
            """,
        )
        _, g, _ = read_gwb_database(gwb)
        notes = join(metadata(g, "notes"), "\n")
        for what in ("`a stray line` outside any section", "section `oxides` is not read", "Undef1+ rests on a species", "Qq+ is made of Qq")
            @test occursin(what, notes)
        end
        caoh = only(eachrow(g[g.name .== "CaOH+", :]))
        k = only(caoh.formation.terms)[2]
        @test k.range == (T_ZERO_CELSIUS, 300.0 + T_ZERO_CELSIUS)
        @test first(ChemistryLab._log10K(caoh.formation, T_STANDARD)) ≈ -12.8

        # A data0 file with a species defined twice, one without its grid, one
        # resting on a species the file does not define, one whose reaction does
        # not balance its stated composition, and a basis species of a
        # pseudo-element.
        bar = "+" * "-"^68
        grid = "     log k grid (0-25-60-100/150-200-250-300 C) =\n     1.0  1.0  1.0  1.0\n     1.0  1.0  1.0  1.0"
        block(name, rest) = "$(rpad(name, 24))\n$rest\n$bar"
        caoh_block = block("CaOH+", "     charge  =   1.0\n     3 chemical elements =\n      1.0000 ca    1.0000 h    1.0000 o\n     4 species in reaction =\n  -1.0000 CaOH+    -1.0000 H+\n   1.0000 Ca+2      1.0000 H2O\n$grid")
        eq36 = joinpath(dir, "data0.flt")
        write(
            eq36, join(
                [
                    "data0.flt.R0\na data0 file with faults\n$bar\ndata0 parameters\ntemperatures\n0.01  25.0  60.0  100.0\n150.0  200.0  250.0  300.0\n$bar",
                    "basis species\n$bar",
                    block("H+", "     charge  =   1.0\n     1 chemical elements =\n      1.0000 h"),
                    block("H2O", "     charge  =   0.0\n     2 chemical elements =\n      2.0000 h    1.0000 o"),
                    block("Ca+2", "     charge  =   2.0\n     1 chemical elements =\n      1.0000 ca"),
                    block("Qq+", "     charge  =   1.0\n     1 chemical elements =\n      1.0000 qq"),
                    "aqueous species\n$bar",
                    caoh_block, caoh_block,
                    block("Nogrid+", "     charge  =   1.0\n     2 species in reaction =\n  -1.0000 Nogrid+   1.0000 Ca+2"),
                    block("Undef1+", "     charge  =   1.0\n     2 species in reaction =\n  -1.0000 Undef1+   1.0000 Xx+\n$grid"),
                    block("Unbal+2", "     charge  =   2.0\n     1 chemical elements =\n      2.0000 ca\n     2 species in reaction =\n  -1.0000 Unbal+2   1.0000 Ca+2\n$grid"),
                    "stop.\n",
                ], "\n",
            ),
        )
        _, e, _ = read_eq36_database(eq36)
        notes = join(metadata(e, "notes"), "\n")
        for what in (
                "CaOH+ is defined again", "Nogrid+ has no log K grid", "Undef1+ rests on a species",
                "reaction of Unbal+2 does not balance", "Qq+ is made of Qq",
            )
            @test occursin(what, notes)
        end
        @test "CaOH+" in e.name
    end
end
