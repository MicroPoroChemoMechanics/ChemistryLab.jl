using JSON
using TOML

@testsection "Databases" begin
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

        # Scan every shipped ThermoFun unit, including coefficient metadata.
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
        for (directory, _, files) in walkdir(datapath())
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

    @testset "the merged database: what the .dat file actually adds" begin
        # `merge_json` exists because Cemdata18's ThermoFun file and PHREEQC's
        # `.dat` file carry DIFFERENT things about the same phases, and the
        # merged database this package ships is the result.
        #
        # What it adds is not species -- a reasonable reading of the word
        # "merged", and the wrong one. Both files describe the same 228
        # substances. What the `.dat` file brings is the REACTIONS, and with them
        # the phase-volume data that makes volumes and porosity available on one
        # consistent dataset.
        #
        # Asserted here so the claim in the manual is checked rather than
        # believed, and so that a future regeneration cannot quietly change it.
        base = JSON.parsefile(datapath("cemdata18-thermofun.json"); dicttype = Dict{String, Any})
        merged = JSON.parsefile(datapath("cemdata18-merged.json"); dicttype = Dict{String, Any})

        syms(db) = Set(String(s["symbol"]) for s in db["substances"])
        @test length(syms(base)) == 228
        @test syms(merged) == syms(base)          # identical in substances

        nrxn(db) = length(get(db, "reactions", []))
        @test nrxn(base) == 7                     # and the numbers the manual quotes
        @test nrxn(merged) == 148

        # Every reaction is usable as a reaction: it has a symbol, and it has
        # something on its left-hand side.
        for r in merged["reactions"]
            @test haskey(r, "symbol") && !isempty(String(r["symbol"]))
            @test !isempty(get(r, "reactants", []))
        end

        # A REACTANT IS NOT ALWAYS A DECLARED SUBSTANCE SYMBOL, and asserting
        # that it is was wrong here before: of the 751 reactant entries, 146
        # name their participant by FORMULA rather than by symbol
        # (`Mg6Al2(OH)18(H2O)3`, `Ca2Al(OH)7(H2O)3`, `(CaO)3Al2O3`), and one of
        # them is the electron, `e-`, which is no substance at all. The `.dat`
        # file names participants the way PHREEQC writes them, and the merge
        # carries that through rather than rewriting it.
        #
        # What is checked instead is that the two naming conventions are the
        # only ones: a reactant is either a declared symbol, or something the
        # formula parser accepts, or the electron.
        known = syms(merged)
        for r in merged["reactions"]
            for part in get(r, "reactants", [])
                sym = String(part["symbol"])
                sym in known && continue
                sym == "e-" && continue
                @test (
                    try
                        Species(sym)
                        true
                    catch
                        false
                    end
                )
            end
        end
    end

    # Test parse_reaction_stoich_cemdata
    @testset "parse_reaction_stoich_cemdata" begin
        # Test basic reaction parsing
        reaction = "CaCO3 = Ca+2 + CO3-2"
        reactants, equation, comment = ChemistryLab.parse_reaction_stoich_cemdata(reaction)
        @test length(reactants) == 3
        @test any(r -> r["symbol"] == "CaCO3" && r["coefficient"] == -1.0, reactants)
        @test any(r -> r["symbol"] == "Ca+2" && r["coefficient"] == 1.0, reactants)
        @test any(r -> r["symbol"] == "CO3-2" && r["coefficient"] == 1.0, reactants)

        # Test reaction with comment
        reaction_with_comment = "H2O = H+ + OH- # water dissociation"
        reactants, equation, comment = ChemistryLab.parse_reaction_stoich_cemdata(reaction_with_comment)
        @test comment == "water dissociation"
        @test length(reactants) == 3

        # Test reaction with coefficients
        reaction_with_coef = "2H2O = 2H+ + 2OH-"
        reactants, equation, comment = ChemistryLab.parse_reaction_stoich_cemdata(reaction_with_coef)
        @test length(reactants) == 3
        @test any(r -> r["symbol"] == "H2O" && r["coefficient"] == -2.0, reactants)
    end

    # Test parse_float_array
    @testset "parse_float_array" begin
        line = "-analytical_expression 1.23 -4.56 7.89 # some comment"
        values = ChemistryLab.parse_float_array(line)
        @test length(values) == 3
        @test values ≈ [1.23, -4.56, 7.89]

        # Test empty line
        @test isempty(ChemistryLab.parse_float_array(""))

        # Test line with only comments
        @test isempty(ChemistryLab.parse_float_array("# only comment"))
    end

    # Test parse_phases
    @testset "parse_phases" begin
        dat_content = """
        PHASES
        Calcite
        CaCO3 = Ca+2 + CO3-2
        -log_K -8.48
        -analytical_expression 1.23 -4.56 7.89

        Portlandite
        Ca(OH)2 = Ca+2 + 2OH-
        -log_K -5.2
        """

        phases = ChemistryLab.parse_phases(dat_content)
        @test haskey(phases, "Calcite")
        @test haskey(phases, "Portlandite")
        @test phases["Calcite"]["logKr"]["values"][1] ≈ -8.48
        @test length(phases["Calcite"]["analytical_expression"]) == 3
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

        # The banner a reader sees stays machine-independent.
        @test ChemistryLab.display_data_path(bundled) ==
            joinpath("data", "cemdata18-thermofun.json")
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
        @test any(e -> e["name"] == "AFm", data["solid_solution"])
    end

    rm(tmp; force = true)
    rm(tmp_missing; force = true)
end

@testsection "every unit the reader actually reads can be parsed" begin
    # A GUARD, not a discovery. `extract_unit` falls back on a default unit when
    # a string does not parse, and a fallback is silent: a value declared in one
    # unit would then be used as though it were in another. PR #59 recorded that
    # `cal` already takes that fallback with the installed DynamicQuantities, and
    # left the scientific question open. This settles it for the shipped data.
    #
    # Measured over `data/*.json`: seven distinct unit strings, of which two do
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
    for f in filter(endswith(".json"), readdir(datapath(); join = true))
        collect_units!(JSON.parsefile(f; dicttype = Dict{String, Any}), "")
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
end
