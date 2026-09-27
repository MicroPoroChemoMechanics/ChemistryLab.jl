# How a database published by others is obtained: the local directory, the
# cache, the download and its checksum, the manual installation, the derived
# databases, and the messages a user reads when something is missing.
#
# The sources are simulated with `file://` URLs in a temporary directory, so no
# test needs the network except the last one, which resolves the real Cemdata18
# file as any user's first call would.

using SHA

@testsection "Databases obtained from their publishers" begin
    CL = ChemistryLab
    digest(bytes) = bytes2hex(sha256(bytes))

    # A fake publisher: a file, its checksum, and the registry entry naming it.
    function publisher(dir; name = "fake-thermofun.json", content = "{\"substances\": []}", urls = nothing)
        src = joinpath(dir, "published-" * name)
        write(src, content)
        entry = CL.ThirdPartyDatabase(
            name, "a fake database", urls === nothing ? ["file://" * src] : urls,
            digest(content), "https://example.org/fake", "a fake license", "a fake citation",
        )
        return entry, src
    end

    @testset "download, verification and cache" begin
        mktempdir() do dir
            cache = joinpath(dir, "cache")
            bad = joinpath(dir, "damaged.json")
            write(bad, "not the right content")
            db, src = publisher(dir)
            db = CL.ThirdPartyDatabase(db.name, db.title, ["file://" * bad, "file://" * src], db.sha256, db.page, db.license, db.citation)

            # The first address serves the wrong content and is passed over; the
            # second is verified and cached, and the download is announced once.
            path = @test_logs (:info, r"downloaded `fake-thermofun.json`") CL._obtain(db, nothing, cache)
            @test path == joinpath(cache, db.name)
            @test read(path, String) == read(src, String)

            # Now cached: the publisher can disappear.
            rm(src)
            @test CL._obtain(db, nothing, cache) == path

            # A damaged cached file is removed rather than used, and nothing can
            # replace it here, so the failure says what was tried.
            write(path, "damaged in the cache")
            err = try
                CL._obtain(db, nothing, cache)
                nothing
            catch e
                e
            end
            @test err isa DatabaseUnavailable
            @test !isfile(path)
            msg = sprint(showerror, err)
            @test occursin("could not obtain `fake-thermofun.json`", msg)
            @test occursin("Downloading it failed", msg)
            @test occursin("differs from the validated version", msg)
            @test occursin("install_database", msg)
            @test occursin("CHEMISTRYLAB_DATABASE_DIR", msg)
            @test occursin("a fake citation", msg)

            # With `download = false` nothing is attempted.
            db2, _ = publisher(dir; name = "other-thermofun.json")
            @test_throws DatabaseUnavailable CL._obtain(db2, nothing, cache; download = false)
        end
    end

    @testset "a database no program can download" begin
        mktempdir() do dir
            db, _ = publisher(dir; name = "manual.dat", urls = String[])
            err = try
                CL._obtain(db, nothing, joinpath(dir, "cache"))
                nothing
            catch e
                e
            end
            @test err isa DatabaseUnavailable
            msg = sprint(showerror, err)
            @test occursin("web browser", msg)
            @test occursin("https://example.org/fake", msg)
            @test occursin("install_database", msg)
            @test !occursin("Downloading it failed", msg)
        end
    end

    @testset "the user's own directory comes first" begin
        mktempdir() do dir
            db, src = publisher(dir)
            mine = joinpath(dir, "mine")
            mkpath(joinpath(mine, "old"))
            mkpath(joinpath(mine, "new"))
            # Two copies under the same name: the validated one is chosen,
            # whatever the order of the subdirectories.
            write(joinpath(mine, "old", db.name), "an older release")
            cp(src, joinpath(mine, "new", db.name))
            @test CL._obtain(db, mine, joinpath(dir, "cache")) == joinpath(mine, "new", db.name)
            # Only a different version: it is used, and said so once.
            rm(joinpath(mine, "new"); recursive = true)
            @test_logs (:warn, r"not the version ChemistryLab") CL._obtain(db, mine, joinpath(dir, "cache"))
            @test CL._obtain(db, mine, joinpath(dir, "cache")) == joinpath(mine, "old", db.name)
        end
    end

    @testset "install_database" begin
        mktempdir() do dir
            saved = CL._CACHE_OVERRIDE[]
            try
                CL._CACHE_OVERRIDE[] = joinpath(dir, "cache")
                real = CL.THIRD_PARTY_DATABASES["CEMDATA18-31-03-2022-phaseVol.dat"]
                wrong = joinpath(dir, real.name)
                write(wrong, "# not the Empa file")
                err = try
                    install_database(wrong)
                    nothing
                catch e
                    e
                end
                @test err isa ArgumentError
                @test occursin("not the version", sprint(showerror, err))
                @test occursin("force = true", sprint(showerror, err))
                # Forced: installed, announced, and then recognized from the cache.
                path = @test_logs (:warn, r"not validated") install_database(wrong; force = true)
                @test isfile(path)
                @test CL._cached(real, CL.database_cache()) == path
                @test_throws ArgumentError install_database(wrong; name = "unknown.dat")
                @test_throws ArgumentError install_database(joinpath(dir, "absent.dat"); name = real.name)
            finally
                CL._CACHE_OVERRIDE[] = saved
            end
        end
    end

    @testset "names and derived databases" begin
        @test_throws ArgumentError database_path("no-such-database.json")
        @test CL.is_database_name("cemdata18-thermofun.json")
        @test CL.is_database_name("cemdata18-zeolites.json")
        @test !CL.is_database_name("solid_solutions.toml")
        # A data file of the package is not a database, and stays under data/.
        @test datapath("solid_solutions.toml") == joinpath(pkgdir(ChemistryLab), "data", "solid_solutions.toml")

        # The key of a derived database moves with its base.
        mktempdir() do dir
            d = CL.DERIVED_DATABASES["cemdata18-zeolites.json"]
            a, b = joinpath(dir, "a.json"), joinpath(dir, "b.json")
            write(a, "{}")
            write(b, "{ }")
            @test CL._derived_key(d, a) != CL._derived_key(d, b)
            @test CL._derived_key(d, a) == CL._derived_key(d, a)
        end
    end

    @testset "the real Cemdata18 file resolves to the validated version" begin
        path = datapath("cemdata18-thermofun.json")
        @test isfile(path)
        @test CL._digest(path) == CL.THIRD_PARTY_DATABASES["cemdata18-thermofun.json"].sha256
        @test CL.display_data_path(path) == "cemdata18-thermofun.json"
        rows = database_info(devnull)
        @test any(r -> r.name == "cemdata18-thermofun.json" && r.status in (:local, :cached), rows)
        @test any(r -> r.name == "CEMDATA18-31-03-2022-phaseVol.dat", rows)
    end
end
