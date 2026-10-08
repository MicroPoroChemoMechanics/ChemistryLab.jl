using JSON
using TOML

@testsection "Databases" begin
    @testset "progress bars only on a terminal" begin
        # A test runner's output is usually not a terminal, so the readers draw
        # no bar there; this checks the step they take when there is one.
        p = ChemistryLab.ProgressMeter.Progress(2; output = devnull)
        ChemistryLab._tick!(p)
        @test p.counter == 1
        @test ChemistryLab._tick!(nothing) === nothing
    end

    @testset "ThermoFun metadata is data" begin
        parse_unit = ChemistryLab.extract_unit
        classify = ChemistryLab.extract_classification
        for fallback in (AS_UNDEF, SC_UNDEF)
            # Every member of the enum round-trips through its own name. This
            # loop is why `AS_LIQUID` no longer appears in the list below: it
            # used to be an unsupported label that fell back, and is now a member
            # like any other, so `instances` covers it here.
            for value in instances(typeof(fallback))
                @test classify(Dict("0" => string(value)), fallback) == value
            end
            for value in (
                    missing, nothing, Dict(), Dict("0" => "unknown"),
                    Dict("0" => "AS_SOLID_SOLUTION"),   # a name no enum carries
                    Dict("0" => "AS_GAS", "1" => "AS_CRYSTAL"),
                )
                @test classify(value, fallback) == fallback
            end
        end
        for text in (
                "1", "1e-05/K", "J/(mol*K^0.5)", "J/(mol*bar)",
                "K^(-1)", "K^(1//2)", "1/√K", "sqrt(K)", "Constants.c^2 * Hz^2",
            )
            @test parse_unit(text) == uparse(text)
        end
        # A negative numeric literal such as -1 is folded by Julia's parser;
        # exercise actual unary and binary subtraction in unit expressions.
        @test parse_unit("-K") == -u"K"
        @test parse_unit("2*K - K") == u"K"
        @test parse_unit("-(K, K, K)", u"Pa") == u"Pa"
        for text in ("unknown_unit", "K[1]", "@time K", "K; mol", "K = mol", repr("K"), "(")
            @test parse_unit(text, u"Pa") == u"Pa"
        end
        @test parse_unit(missing, u"Pa") == u"Pa"

        # Scan every ThermoFun unit of the databases, including coefficient metadata.
        bundled_units = Set{String}()
        function collect_units(value)
            if value isa AbstractDict
                for (key, child) in value
                    if key == "units"
                        for unit in child
                            push!(bundled_units, unit)
                        end
                    else
                        collect_units(child)
                    end
                end
            elseif value isa AbstractVector
                foreach(collect_units, value)
            end
        end
        literature_dir = joinpath(datapath(), "literature")
        for (directory, _, files) in walkdir(datapath())
            # `data/literature` is not ThermoFun: its files are read and checked
            # by `read_literature` (test/literature.jl), and a text column of
            # one of their tables has no unit at all.
            startswith(directory, literature_dir) && continue
            for file in files
                endswith(file, ".json") || continue
                collect_units(JSON.parsefile(joinpath(directory, file); dicttype = Dict{String, Any}))
            end
        end

        for unit in bundled_units
            # Compare with the old reader, including its unknown-unit fallback.
            expected = try
                uparse(unit)
            catch
                u"1"
            end
            @test ChemistryLab.is_unit_expression(Meta.parse(unit))
            @test parse_unit(unit) == expected
        end

        mktempdir() do directory
            sentinel = joinpath(directory, "metadata-executed")
            payload = "touch($(repr(sentinel)))"
            for text in (
                    payload, "Base.$payload", "($payload; K)", "K * $payload",
                    "(x -> $payload)(K)", "getfield(Base, :touch)($(repr(sentinel)))",
                )
                @test parse_unit(text, u"Pa") == u"Pa"
                @test !ispath(sentinel)
            end
            database = JSON.parsefile(datapath("cemdata18-thermofun.json"); dicttype = Dict{String, Any})
            original = deepcopy(first(database["substances"]))
            database["substances"] = [deepcopy(original)]
            filename = joinpath(directory, "metadata.json")
            # Exercise the public file loader for each formerly evaluated field.
            for field in ("aggregate_state", "class_", "sm_gibbs_energy")
                substance = deepcopy(original)
                if field == "sm_gibbs_energy"
                    substance[field] = Dict("values" => [1.0], "units" => [payload])
                else
                    substance[field] = Dict("0" => "($payload; AS_GAS)")
                end
                database["substances"] = [substance]
                write(filename, JSON.json(database))
                species = only(build_species(filename))
                @test !ispath(sentinel)
                if field == "aggregate_state"
                    @test aggregate_state(species) == AS_UNDEF
                elseif field == "class_"
                    @test class(species) == SC_UNDEF
                end
            end
        end
    end

    @testset "a record without its reference state is read at 298.15 K and 1 bar" begin
        # REGRESSION. One substance of the slop98 organic database, `Eth@`,
        # carries no `Tst`; the reader took the missing value into `missing * u"K"`
        # and the whole database failed to load.
        db = JSON.parsefile(datapath("cemdata18-thermofun.json"); dicttype = Dict{String, Any})
        rec = only(s for s in db["substances"] if s["symbol"] == "Portlandite")
        full = only(build_species(ChemistryLab.DataFrame(ChemistryLab.Tables.dictrowtable(JSON.parse(JSON.json([rec]))))))
        bare = deepcopy(rec)
        delete!(bare, "Tst")
        delete!(bare, "Pst")
        sp = only(build_species(ChemistryLab.DataFrame(ChemistryLab.Tables.dictrowtable(JSON.parse(JSON.json([bare]))))))
        @test sp.Tref == 298.15u"K"
        @test sp.Pref == 1.0e5u"Pa"
        @test sp[:ΔₐG⁰](T = 298.15) ≈ full[:ΔₐG⁰](T = 298.15) rtol = 1.0e-12
    end

    @testset "a solid recorded with a zero molar volume has none" begin
        # Cemdata18 prints "not defined" for the volume of the amorphous and
        # microcrystalline Fe(OH)3, and its ThermoFun file writes 0.
        sp = Dict(symbol(s) => s for s in build_species(datapath("cemdata18-thermofun.json"); verbose = false))
        @test !haskey(sp["Fe(OH)3(am)"], :V⁰)
        @test !haskey(sp["Fe(OH)3(mic)"], :V⁰)
        @test haskey(sp["Portlandite"], :V⁰)
        cs = ChemicalSystem([sp[s] for s in ("H2O@", "Fe(OH)3(am)", "Portlandite")])
        st = ChemicalState(cs)
        set_quantity!(st, "Fe(OH)3(am)", 1.0e-3u"mol")
        set_quantity!(st, "Portlandite", 1.0e-3u"mol")
        @test missing_molar_volumes(st) == ["Fe(OH)3(am)"]
    end

    # No path in a script, a documentation block or a test may depend on the
    # working directory. `resolve_data_path` tries the working directory first,
    # so a call that already resolves keeps resolving to the very same file: the
    # fallbacks can only turn a failure into a success.
    @testset "data path resolution" begin
        @test isdir(datapath())
        @test isfile(datapath("solid_solutions.toml"))
        @test datapath("experimental", "README.md") ==
            joinpath(datapath(), "experimental", "README.md")

        resolve = ChemistryLab.resolve_data_path
        bundled = datapath("cemdata18-thermofun.json")
        calorimetry = "smilauer2025-116-cemI-52.5R-ladce.csv"

        mktempdir() do dir
            cd(dir) do
                # Every form a script or a documentation page has used, from a
                # working directory holding none of them.
                @test resolve("cemdata18-thermofun.json") == bundled
                @test resolve("data/cemdata18-thermofun.json") == bundled
                @test resolve("../../../data/cemdata18-thermofun.json") == bundled
                @test resolve(joinpath("experimental", calorimetry)) ==
                    datapath("experimental", calorimetry)

                # A name that is not a bundled data file fails loudly rather
                # than silently resolving to something else.
                @test_throws ArgumentError resolve("no-such-database.json")

                # The working directory wins over a bundled file of the same name.
                mkdir("data")
                local_copy = joinpath("data", "cemdata18-thermofun.json")
                write(local_copy, "{}")
                @test resolve(local_copy) == local_copy
            end
        end

        # The banner a reader sees stays machine-independent: a database by its
        # name, wherever the cache put it; a data file of the package relative
        # to the package root.
        @test ChemistryLab.display_data_path(bundled) == "cemdata18-thermofun.json"
        @test ChemistryLab.display_data_path(datapath("solid_solutions.toml")) ==
            joinpath("data", "solid_solutions.toml")
        # Outside the package, shown unchanged -- including a path built the way
        # the running platform builds one, since `tempdir()` is on another drive
        # from the checkout on a Windows runner and that is exactly the case a
        # `relpath`-based test would get wrong.
        outside = joinpath(tempdir(), "elsewhere", "foo.json")
        @test ChemistryLab.display_data_path(outside) == outside
        @test ChemistryLab.display_data_path("/elsewhere/foo.json") == "/elsewhere/foo.json"
        # A sibling of the package root whose name begins with it is not inside
        # it, which the prefix test has to get right.
        sibling = pkgdir(ChemistryLab) * "-elsewhere"
        @test ChemistryLab.display_data_path(joinpath(sibling, "f.json")) ==
            joinpath(sibling, "f.json")
    end

end

@testsection "build_solid_solutions" begin
    # ── Helpers: build a minimal species dict without loading real databases ──
    _make(sym, cls = SC_COMPONENT) =
        Species(sym; aggregate_state = AS_CRYSTAL, class = cls)

    # Build a fake dict that mimics what build_species returns
    dict = Dict(
        "TobD" => _make("TobD"),
        "TobH" => _make("TobH"),
        "JenH" => _make("JenH"),
        "JenD" => _make("JenD"),
        "Ms" => _make("Ms"),
        "Mc" => _make("Mc"),
    )

    # ── Two temp TOMLs: one clean, one with a missing phase ─────────────────
    toml_clean = """
    [[solid_solution]]
    name        = "CSHQ"
    end_members = ["TobD", "TobH", "JenH", "JenD"]
    model       = "ideal"

    [[solid_solution]]
    name        = "AFm"
    end_members = ["Ms", "Mc"]
    model       = "redlich_kister"
    a0          = 3000.0
    a1          = 500.0
    a2          = 0.0
    """
    toml_with_missing = toml_clean * """
        [[solid_solution]]
        name        = "Missing_phase"
        end_members = ["NonExistent1", "NonExistent2"]
        model       = "ideal"
        """
    tmp = tempname() * ".toml"
    tmp_missing = tempname() * ".toml"
    write(tmp, toml_clean)
    write(tmp_missing, toml_with_missing)

    # Pre-build phases from the clean TOML (no warnings, reused across sub-tests)
    phases = build_solid_solutions(tmp, dict)

    @testset "load two phases" begin
        @test length(phases) == 2

        cshq = first(filter(ss -> name(ss) == "CSHQ", phases))
        afm = first(filter(ss -> name(ss) == "AFm", phases))

        @test length(end_members(cshq)) == 4
        @test model(cshq) isa IdealSolidSolutionModel

        @test length(end_members(afm)) == 2
        @test model(afm) isa RedlichKisterModel
        @test model(afm).a0 ≈ 3000.0
        @test model(afm).a1 ≈ 500.0
    end

    @testset "end-members requalified to SC_SSENDMEMBER" begin
        for ss in phases
            for em in end_members(ss)
                @test class(em) == SC_SSENDMEMBER
            end
        end
    end

    @testset "original dict species class untouched" begin
        # SolidSolutionPhase requalifies internally; dict values must be unchanged
        @test all(class(s) == SC_COMPONENT for s in values(dict))
    end

    @testset "skip_missing = true warns and skips" begin
        phases_skip = @test_logs (:warn, r"skipping.*Missing_phase") match_mode = :any build_solid_solutions(
            tmp_missing, dict; skip_missing = true
        )
        @test length(phases_skip) == 2   # Missing_phase skipped
    end

    @testset "skip_missing = false raises on missing end-members" begin
        @test_throws ErrorException build_solid_solutions(
            tmp_missing, dict; skip_missing = false
        )
    end

    @testset "data/solid_solutions.toml parses without error" begin
        toml_path = datapath("solid_solutions.toml")
        # Just check the file is valid TOML and has solid_solution entries
        data = TOML.parsefile(toml_path)
        @test haskey(data, "solid_solution")
        @test length(data["solid_solution"]) >= 2
        @test any(e -> e["name"] == "CSHQ", data["solid_solution"])
        @test any(e -> e["name"] == "AFm_SO4_OH", data["solid_solution"])
        # Cemdata18 has no monosulfate-monocarbonate solid solution, and the
        # file no longer declares one.
        @test !any(e -> Set(e["end_members"]) == Set(["monosulphate12", "monocarbonate"]), data["solid_solution"])
    end

    @testset "published Guggenheim parameters, read rather than copied" begin
        # `guggenheim = "<key>:<pair>"` reads the dimensionless parameters from
        # data/literature and makes them a = α R T at 298.15 K; the gap they
        # open is the one the article prints, read from the same row.
        dcem = Dict(
            symbol(s) => s for s in build_species(datapath("cemdata18-thermofun.json"); verbose = false)
        )
        byss = Dict(name(p) => p for p in build_solid_solutions(datapath("solid_solutions.toml"), dcem))
        RT = ChemistryLab.R_GAS * 298.15
        # The AFm gap is printed in X(OH), the first end-member; the AFt gap in
        # X(SO4), the second; the two Al/Fe gaps in X(Al), the first.
        for (ss, pair, second) in (
                ("AFm_SO4_OH", "AFm SO4/OH", false), ("AFt_SO4_CO3", "AFt SO4/CO3", true),
                ("AFt_AlFe", "AFt Al/Fe", false), ("AFm_AlFe", "AFm Al/Fe", false),
            )
            p = literature_row("Lothenbach2019", "guggenheim_parameters", pair)
            m = model(byss[ss])
            @test m isa RedlichKisterModel
            @test m.a0 ≈ ustrip(p.alpha0) * RT
            @test m.a1 ≈ ustrip(p.alpha1) * RT
            @test byss[ss].instances == 2
            x = common_tangent(byss[ss])
            gap = second ? (1 - x[2], 1 - x[1]) : (x[1], x[2])
            @test gap[1] ≈ ustrip(p.gap_from) atol = 5.0e-3
            @test gap[2] ≈ ustrip(p.gap_to) atol = 5.0e-3
        end
        @test_throws ErrorException ChemistryLab._guggenheim_model("Lothenbach2019", "AFm")
        @test_throws KeyError ChemistryLab._guggenheim_model("Lothenbach2019:AFm SO4/CO3", "AFm")
    end

    rm(tmp; force = true)
    rm(tmp_missing; force = true)
end

@testsection "every unit the reader actually reads can be parsed" begin
    # A GUARD, not a discovery. `extract_unit` falls back on a default unit when
    # a string does not parse, and a fallback is silent: a value declared in one
    # unit would then be used as though it were in another. PR #59 recorded that
    # `cal` already takes that fallback with the installed DynamicQuantities, and
    # left the scientific question open. This settles it for the databases read.
    #
    # Measured over every ThermoFun database: seven distinct unit strings, of which two do
    # not parse -- `cal/(mol*bar)` on `eos_hkf_coeffs` and `kbar` on
    # `m_expansivity`. NEITHER REACHES THE PARSER.
    #
    #   * `eos_hkf_coeffs` is converted by `HKF_SI_CONVERSIONS`, an explicit
    #     SUPCRT-to-SI table, precisely because the JSON's own unit metadata is
    #     wrong for `a3` and `a4`. It never calls `extract_unit`.
    #   * `m_expansivity` is not read by this package at all.
    #
    # So no value the package uses is off by a factor of 4.184. This test is what
    # keeps that true when a database is added or a field starts being read.
    read_fields = Set(
        [
            "sm_gibbs_energy", "sm_enthalpy", "sm_entropy_abs", "sm_heat_capacity_p",
            "sm_volume", "drsm_gibbs_energy", "drsm_enthalpy", "drsm_entropy_abs",
            "drsm_volume", "logKr", "m_heat_capacity_ft_coeffs",
        ]
    )
    units = Set{String}()
    function collect_units!(o, key)
        if o isa AbstractDict
            if key in read_fields && haskey(o, "units") && o["units"] isa AbstractVector
                for u in o["units"]
                    u isa AbstractString && push!(units, u)
                end
            end
            for (k, v) in o
                collect_units!(v, k)
            end
        elseif o isa AbstractVector
            for v in o
                collect_units!(v, key)
            end
        end
        return nothing
    end
    for name in sort!(vcat(collect(keys(ChemistryLab.THIRD_PARTY_DATABASES)), collect(keys(ChemistryLab.DERIVED_DATABASES))))
        endswith(name, ".json") || continue
        collect_units!(JSON.parsefile(datapath(name); dicttype = Dict{String, Any}), "")
    end

    @test !isempty(units)                       # the walk found something
    parses(u) = ChemistryLab.is_unit_expression(Meta.parse(u)) && (
        try
            uparse(u)
            true
        catch
            false
        end
    )
    @test isempty(filter(!parses, collect(units)))

    # Fractional exponents are among them and must survive both the security
    # guard and `uparse`: `J/(mol*K^0.5)` is a heat-capacity coefficient.
    @test ChemistryLab.is_unit_expression(Meta.parse("J/(mol*K^0.5)"))
    @test dimension(uparse("J/(mol*K^0.5)")) == dimension(u"J/mol" / sqrt(u"K"))
end

@testsection "AS_LIQUID is read from the database instead of being lost" begin
    # One species in `slop98-inorganic-thermofun.json` declares `AS_LIQUID`:
    # metallic mercury. Until the enum carried that member the label matched
    # nothing and the import fell back on `AS_UNDEF`, silently discarding what
    # the file said.
    @test AS_LIQUID isa AggregateState
    @test Int(AS_LIQUID) == 4                   # appended, so nothing renumbered
    @test Int(AS_UNDEF) == 0 && Int(AS_GAS) == 3

    sp = build_species(datapath("slop98-inorganic-thermofun.json"); verbose = false)
    hg = only(s for s in sp if symbol(s) == "Hg")
    @test aggregate_state(hg) == AS_LIQUID

    # AND IT IS NOT FILTERED OUT ANYWHERE THE OTHERS ARE KEPT. Adding a member to
    # an enum turns every hand-written list of "all of them" into a filter,
    # silently. `idx_speciation` carried one: four states named when there were
    # four, so its intent was not to filter on state at all. It now says
    # `instances`, and this is the assertion that would have caught the day it
    # stopped meaning that.
    # `idx_speciation` returns a BOOLEAN MASK, not a list of indices, so the
    # assertion has to read it as one: `length` of the mask is the number of
    # candidates, not the number kept, and an empty mask never happens.
    @test only(ChemistryLab.idx_speciation([hg], collect(keys(atoms(hg)))))
    # Naming a state still filters, so the default is permissive rather than the
    # filter being inert.
    @test !any(
        ChemistryLab.idx_speciation(
            [hg], collect(keys(atoms(hg))); aggregate_state = [AS_GAS]
        )
    )
end

@testsection "every ThermoFun substance label has a member to land on" begin
    # THE GENERAL FORM OF THE `AS_LIQUID` DEFECT. `extract_classification` matches
    # a label by name against `instances`, and falls back when there is no match
    # -- silently, because a fallback is a valid value. So a label the enum does
    # not carry is data the import throws away without saying so, and `AS_LIQUID`
    # was one: metallic mercury read as `AS_UNDEF` until the member was added.
    #
    # Testing that one label would leave the next one to be found the same way.
    # This walks the `substances` of every database instead and requires
    # each label to resolve.
    #
    # `substances` and not the whole file, deliberately. The `elements` section
    # uses `ELEMENT` (221 occurrences) and `CHARGE` (7), which are ThermoFun's
    # classes for ELEMENTS; `Class` describes substances, so it is right not to
    # carry them, and widening this walk would turn that correctness into a
    # failure.
    agg_names = Set(string.(instances(AggregateState)))
    class_names = Set(string.(instances(Class)))
    unknown_agg, unknown_class = Set{String}(), Set{String}()
    for name in sort!(vcat(collect(keys(ChemistryLab.THIRD_PARTY_DATABASES)), collect(keys(ChemistryLab.DERIVED_DATABASES))))
        endswith(name, ".json") || continue
        d = JSON.parsefile(datapath(name); dicttype = Dict{String, Any})
        for item in get(d, "substances", Any[])
            item isa AbstractDict || continue
            a = get(item, "aggregate_state", nothing)
            a isa AbstractDict && for v in values(a)
                v isa AbstractString && !(v in agg_names) && push!(unknown_agg, v)
            end
            c = get(item, "class_", nothing)
            c isa AbstractDict && for v in values(c)
                v isa AbstractString && !(v in class_names) && push!(unknown_class, v)
            end
        end
    end
    @test isempty(unknown_agg)
    @test isempty(unknown_class)

    # The walk found something, so an empty result above means agreement and not
    # a broken loop.
    @test "AS_LIQUID" in agg_names
    @test "SC_AQSOLUTE" in class_names
end

@testsection "PHREEQC files: -gamma, and the options of PHASES" begin
    dat = datapath("phreeqc.dat")
    γp = phreeqc_gamma_parameters(dat)
    # Master and secondary species alike, by the symbol convention of the package.
    @test γp["H+"] == (9.0, 0.0)
    @test γp["Ca+2"] == (5.0, 0.165)
    # Na+ carries two -gamma lines in phreeqc.dat; PHREEQC keeps the second.
    @test γp["Na+"] == (4.08, 0.082)
    @test length(γp) > 80

    mktempdir() do dir
        f = joinpath(dir, "tiny.dat")
        write(
            f, """
            SOLUTION_SPECIES
            Ca+2 + H2O = CaOH+ + H+
                log_k -12.78
                -gamma 4 0.064.064
            CO3-2 + H+ = HCO3-
                -GAMMA 5.4 0.0
            CO3-2 + 2 H+ = CO2 + H2O
                -g 0 0.1
            PHASES
                -Vm 1.0
            Calcite
                CaCO3 = CO3-2 + Ca+2
                -LOG_K -8.48
                -Analytic 1 2 3 4.60517 5 6
                -Vm 36.9 cm3/mol
            zeoliteP_Ca
                CaAl2Si2O8(H2O)4.5 = Ca+2 + 2AlO2- + 2SiO2 + 4.5H2O
                -log_k -20.3
            Broken
                CaO + 2H+ = Ca+2 + H2O
                -log_k not-a-number
            END
            """,
        )
        # A -gamma that does not parse is reported and skipped; case and the
        # abbreviation `-g` are PHREEQC's.
        p = @test_logs (:warn, r"does not parse") phreeqc_gamma_parameters(f)
        @test !haskey(p, "CaOH+")
        @test p["HCO3-"] == (5.4, 0.0)
        @test p["CO2@"] == (0.0, 0.1)

    end
end

@testsection "the molar volumes of the Empa PHREEQC file are Cemdata18's" begin
    # Runs where the Empa file is installed; see `install_database`.
    dat = try
        database_path("CEMDATA18-31-03-2022-phaseVol.dat"; download = false)
    catch err
        err isa DatabaseUnavailable || rethrow()
        nothing
    end
    if dat === nothing
        @info "skipped: `CEMDATA18-31-03-2022-phaseVol.dat` is not installed (see `install_database`)"
        @test_skip false
    else
        _, table, _ = read_phreeqc_database(dat)
        phases = Dict(r.name => r for r in eachrow(table) if r.aggregate_state != AS_AQUEOUS)
        # The table gives volumes in m³/mol, ThermoFun in J/bar: both in cm³/mol.
        vm(p) = ismissing(p.molar_volume) ? nothing : ustrip(us"cm^3/mol", p.molar_volume * u"m^3/mol")
        withvm = [k for (k, p) in phases if vm(p) !== nothing]
        db = JSON.parsefile(datapath("cemdata18-thermofun.json"); dicttype = Dict{String, Any})
        V = Dict(
            String(s["symbol"]) => ustrip(us"cm^3/mol", Float64(s["sm_volume"]["values"][1]) * u"J/(bar*mol)")
                for s in db["substances"] if haskey(s, "sm_volume") && !isempty(s["sm_volume"]["values"])
        )
        # Some are absent from the ThermoFun file (gases, which PHREEQC names
        # `X(g)`, and hydrates it names differently or lacks) and three carry a
        # zero V° there (the two iron hydroxides and FeCO3(pr)).
        common = [k for k in withvm if haskey(V, k) && V[k] != 0]
        # To the rounding of the PHREEQC file, which prints one decimal or two,
        # except INFCNA: 64.51 cm³/mol in the PHREEQC file, whose header records
        # a correction of that end member in 2019, against 69.3 in the ThermoFun
        # file. The two iron hydroxides carry 34 cm³/mol here and "not defined"
        # in Cemdata18, hence none in the ThermoFun file.
        differ = sort([k for k in common if abs(vm(phases[k]) - V[k]) > 0.06])
        @test differ == ["INFCNA"]
        @test vm(phases["Fe(OH)3(am)"]) ≈ 34.0
        @test haskey(phases, "zeoliteP_Ca")
    end
end
