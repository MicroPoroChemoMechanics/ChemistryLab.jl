# The files of `data/literature/`, and the loader that reads them.

using JSON

@testset "published values are data files" begin
    root = dirname(@__DIR__)
    keys_ = available_literature()
    @test !isempty(keys_)

    # The bibliography, read as text: every key and, where there is one, its DOI.
    bib = read(joinpath(root, "docs", "src", "refs.bib"), String)
    entries = Dict{String, Union{Nothing, String}}()
    for m in eachmatch(r"@\w+\{([^,\s]+),(.*?)\n\}"s, bib)
        doi = match(r"\bdoi\s*=\s*\{([^}]*)\}"i, m.captures[2])
        entries[m.captures[1]] = doi === nothing ? nothing : lowercase(strip(doi.captures[1]))
    end

    @testset "$key" for key in keys_
        r = literature(key)
        @test r isa LiteratureRecord
        @test r.key == key
        # The source is a bibliography entry, and the DOI is the one it gives.
        @test haskey(entries, key)
        doi = get(r.source, "doi", nothing)
        @test (doi === nothing ? nothing : lowercase(doi)) == get(entries, key, missing)
        @test haskey(r.transcription, "from")
        @test haskey(r.transcription, "checked_against_source")
        for (_, q) in r.quantities
            @test q isa Traced
            @test startswith(ChemistryLab.source(q), key)
        end
    end

    @testset "the values moved out of the code are the values that were there" begin
        # Unchanged by the move, and still reaching the code that uses them.
        @test literature_value("Powers1948", "w_c_sealed") == 0.42
        @test literature_value("Powers1948", "w_c_saturated") == 0.36
        @test ChemistryLab.POWERS_W_SEALED === literature_value("Powers1948", "w_c_sealed")
        @test powers_alpha_max(0.21) ≈ 0.5
        @test provenance(literature("Powers1948")["w_c_sealed"]) == PROV_PUBLISHED
    end

    @testset "a malformed file is refused, with its name and the field" begin
        dir = mktempdir()
        good = JSON.parsefile(literature_path("Powers1948"); dicttype = Dict{String, Any})
        write_as(d, name) = (p = joinpath(dir, name * ".json"); write(p, JSON.json(d)); p)

        @test ChemistryLab.read_literature(write_as(good, "Powers1948")) isa LiteratureRecord

        bad = deepcopy(good); bad["schema"] = "something-else/0"
        @test_throws ArgumentError ChemistryLab.read_literature(write_as(bad, "Powers1948"))

        # The key must be the file name.
        @test_throws ArgumentError ChemistryLab.read_literature(write_as(good, "Other1948"))

        bad = deepcopy(good); bad["quantities"]["w_c_sealed"]["kind"] = "guessed"
        @test_throws ArgumentError ChemistryLab.read_literature(write_as(bad, "Powers1948"))

        # A unit is unit arithmetic, never code.
        bad = deepcopy(good); bad["quantities"]["w_c_sealed"]["unit"] = "run(`true`)"
        @test_throws ArgumentError ChemistryLab.read_literature(write_as(bad, "Powers1948"))

        bad = deepcopy(good); bad["quantities"]["w_c_sealed"]["value"] = "0.42"
        @test_throws ArgumentError ChemistryLab.read_literature(write_as(bad, "Powers1948"))

        bad = deepcopy(good)
        bad["tables"] = Dict(
            "t" => Dict(
                "columns" => ["a", "b"], "units" => ["1", "1"], "kind" => "published",
                "rows" => [[1.0, 2.0], [3.0]]
            )
        )
        @test_throws ArgumentError ChemistryLab.read_literature(write_as(bad, "Powers1948"))

        @test_throws ArgumentError literature("NoSuchSource2099")
    end

    @testset "tables carry their units, and text columns stay text" begin
        dir = mktempdir()
        d = JSON.parsefile(literature_path("Powers1948"); dicttype = Dict{String, Any})
        d["tables"] = Dict(
            "t" => Dict(
                "columns" => ["phase", "m", "ratio"], "units" => [nothing, "mol/kg", "1"],
                "kind" => "published", "location" => "Table 1",
                "rows" => [["A", 0.1, 1.5], ["B", 0.2, 2.5]]
            )
        )
        d["quantities"]["with_unit"] = Dict(
            "value" => 8.0, "unit" => "J/mol", "kind" => "fitted", "uncertainty" => 0.5
        )
        p = joinpath(dir, "Powers1948.json")
        write(p, JSON.json(d))
        r = ChemistryLab.read_literature(p)
        t = r.tables["t"]
        @test t.phase == ["A", "B"]
        @test t.m ≈ [0.1, 0.2] .* u"mol/kg"
        @test t.ratio == [1.5, 2.5]
        q = r["with_unit"]
        @test ChemistryLab.value(q) ≈ 8.0u"J/mol"
        @test uncertainty(q) ≈ 0.5u"J/mol"
        @test provenance(q) == PROV_FITTED
    end
end
