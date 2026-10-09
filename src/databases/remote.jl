# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

using Downloads: Downloads
using SHA: sha256
using Scratch: get_scratch!

# ── Thermodynamic databases published by others ──────────────────────────────
#
# The thermodynamic databases ChemistryLab reads are published by their authors,
# under their own terms, and ChemistryLab obtains them from those publishers the
# first time they are needed. Each one is identified by the SHA-256 of the
# version the package was validated with, so a file that arrives damaged, or a
# different release, is recognized and never used silently.
#
# Where a file is looked for, in order:
#
#  1. `ENV["CHEMISTRYLAB_DATABASE_DIR"]`, a directory of the user's own copies
#     (searched together with its immediate subdirectories);
#  2. the package's cache, a Scratch space of the depot;
#  3. the publisher, downloaded into the cache after its checksum is verified.
#
# A file ChemistryLab does not download (the publisher's site demands a browser,
# or its terms leave the copy to the user) stops at step 2 with instructions,
# and `install_database` puts the copy the user downloaded into the cache. A
# program downloads only what its publisher states may be copied: the PHREEQC
# databases of the U.S. Geological Survey, under its User Rights Notice, and
# the ThermoHub files, under GPL-3.0. No database is stored in the package.

"""
    ThirdPartyDatabase

A thermodynamic database published outside ChemistryLab: its file name, what it
is, where it is obtained, the SHA-256 of the validated version, and the license
and citation its publisher states. `urls` is empty for a file that has to be
downloaded by hand from `page`.
"""
struct ThirdPartyDatabase
    name::String
    title::String
    urls::Vector{String}
    sha256::String
    page::String
    license::String
    citation::String
end

# ThermoHub publishes the ThermoFun database files in `thermofun/` of its
# repository. Release v1.1.1, archived on Zenodo as doi:10.5281/zenodo.7385311
# with the same five files, byte for byte.
const _THERMOHUB_COMMIT = "a50ef4177f5a7f447e1b1e92ee967b7357b9fa63"
const _THERMOHUB_PAGE = "https://github.com/thermohub/thermohub/tree/v1.1.1/thermofun " *
    "(archive: https://doi.org/10.5281/zenodo.7385311)"
const _THERMOHUB_LICENSE = "GPL-3.0, as stated by the ThermoHub repository"

_thermohub_urls(file) = [
    "https://raw.githubusercontent.com/thermohub/thermohub/$(_THERMOHUB_COMMIT)/thermofun/$file",
    "https://cdn.jsdelivr.net/gh/thermohub/thermohub@$(_THERMOHUB_COMMIT)/thermofun/$file",
]

function _thermohub(file, title, sha, citation)
    return ThirdPartyDatabase(
        file, title, _thermohub_urls(file), sha, _THERMOHUB_PAGE, _THERMOHUB_LICENSE, citation,
    )
end

# The databases distributed with PHREEQC, from `database/` of its repository at
# tag v3.7.3. That is the release of the IPhreeqc engine the PHREEQC oracles of
# the test suite run on (phreeqpython bundles 3.7.3), and a database of a later
# release uses keywords that engine refuses; v3.7.3's `phreeqc.dat` is byte for
# byte that of v3.7.1. The User Rights Notice of PHREEQC
# (https://water.usgs.gov/water-resources/software/PHREEQC/Phreeqc_UserRightsNotice.txt)
# allows the files to be used, copied and distributed, the notice kept and the
# authors and the USGS acknowledged. Each one's own source is cited as its
# header states it.
const _USGS_TAG = "v3.7.3"
const _USGS_PAGE = "https://github.com/phreeqc-dev/phreeqc3/tree/$(_USGS_TAG)/database " *
    "(PHREEQC: https://www.usgs.gov/software/phreeqc-version-3)"
const _USGS_LICENSE = "the U.S. Geological Survey User Rights Notice of PHREEQC: use, copy, " *
    "modification and distribution allowed, the notice kept and the authors and the USGS acknowledged"
const _PHREEQC_CITATION = "Parkhurst & Appelo (2013), Description of input and examples for " *
    "PHREEQC version 3, U.S. Geological Survey Techniques and Methods 6-A43, doi:10.3133/tm6A43"

_usgs_urls(file) = [
    "https://raw.githubusercontent.com/phreeqc-dev/phreeqc3/$(_USGS_TAG)/database/$file",
    "https://cdn.jsdelivr.net/gh/phreeqc-dev/phreeqc3@$(_USGS_TAG)/database/$file",
]

function _usgs(file, title, sha, source)
    return ThirdPartyDatabase(
        file, title, _usgs_urls(file), sha, _USGS_PAGE, _USGS_LICENSE,
        "$(_PHREEQC_CITATION); $source",
    )
end

"""
    THIRD_PARTY_DATABASES

The databases ChemistryLab obtains from their publishers, by file name. See
[`database_info`](@ref) for where each one currently resolves.
"""
const THIRD_PARTY_DATABASES = Dict(
    db.name => db for db in (
            _thermohub(
                "cemdata18-thermofun.json", "Cemdata18, ThermoFun format",
                "00cb56943cf22fc09c32e9c09e070bc1755ef68ac521dfb05a64782c92a5b756",
                "Lothenbach et al. (2019), Cem. Concr. Res. 115, 472-506, doi:10.1016/j.cemconres.2018.04.018",
            ),
            _thermohub(
                "psinagra-12-07-thermofun.json", "PSI/Nagra 12/07, ThermoFun format",
                "c20ad35614a89787b33ae87c3df70a3528d1b0bb0235e1c3c3ceb747558b41a4",
                "Thoenen, Hummel, Berner & Curti (2014), The PSI/Nagra Chemical Thermodynamic Database 12/07, PSI Bericht Nr. 14-04",
            ),
            _thermohub(
                "aq17-thermofun.json", "aq17, ThermoFun format",
                "8638d24a9bf5903b9cf261202502ecc4e4dce05272fdada440053ca508edf8f9",
                "Miron, Wagner, Kulik & Lothenbach (2017), Am. J. Sci. 317, 755-806, doi:10.2475/07.2017.01",
            ),
            _thermohub(
                "slop98-inorganic-thermofun.json", "SUPCRT slop98, inorganic species, ThermoFun format",
                "0d646eb8205a30a6565d0c933dac8ac7b3a9982f724a4eaf46114d066ecdfeb2",
                "the slop98 data file of SUPCRT92, Johnson, Oelkers & Helgeson (1992), Comput. Geosci. 18, 899-947, doi:10.1016/0098-3004(92)90029-Q",
            ),
            _thermohub(
                "slop98-organic-thermofun.json", "SUPCRT slop98, organic species, ThermoFun format",
                "d5ed815a162e20f46a45b35abb95faaf76f19d2fdb8382ba95198c618958ef7a",
                "the slop98 data file of SUPCRT92, Johnson, Oelkers & Helgeson (1992), Comput. Geosci. 18, 899-947, doi:10.1016/0098-3004(92)90029-Q",
            ),
            _usgs(
                "phreeqc.dat", "phreeqc.dat, the default database of PHREEQC",
                "3e819f36a78a134b9e53557e9fd9e640d3616b84282cf7095c60337c49c6357b",
                "its pressure and temperature data: Appelo, Parkhurst & Post (2014), Geochim. " *
                "Cosmochim. Acta 125, 49-67, doi:10.1016/j.gca.2013.10.003",
            ),
            _usgs(
                "llnl.dat", "llnl.dat, PHREEQC format of the LLNL thermo.com.V8.R6.230 data",
                "24d9266ff5c02aab84b2ba295567c8a5328d26f0a00fc69a43b6bb33b3efade2",
                "its data: thermo.com.V8.R6.230, prepared by J. Johnson at Lawrence Livermore " *
                "National Laboratory, converted to PHREEQC format by G. Anderson with D. Parkhurst",
            ),
            _usgs(
                "minteq.v4.dat", "minteq.v4.dat, PHREEQC format of MINTEQA2 version 4.02",
                "1d4dd3f14932ccc4236f862f56689be6b1d1bbaa14067b23277bfcfa7c9bea3b",
                "its data: MINTEQA2 version 4.02 of the U.S. Environmental Protection Agency, " *
                "translated as the PHREEQC release notes state",
            ),
            _usgs(
                "wateq4f.dat", "wateq4f.dat, PHREEQC format of the WATEQ4F database",
                "b12c4a9818c946a882c675c458499d76a71fa0ec0becc8e5b9e175391ffaeb39",
                "its data: Ball & Nordstrom (1991), User's manual for WATEQ4F, U.S. Geological " *
                "Survey Open-File Report 91-183, doi:10.3133/ofr91183",
            ),
            _usgs(
                "sit.dat", "sit.dat, PHREEQC format of ThermoChimie 9b0 with SIT",
                "427d6114ed3f3135054683882319a0852b593d2685f0ccc80855f10ae1c4b840",
                "its data: ThermoChimie, Giffaut et al. (2014), Appl. Geochem. 49, 225-236, " *
                "doi:10.1016/j.apgeochem.2014.05.007",
            ),
            _usgs(
                "pitzer.dat", "pitzer.dat, the Pitzer database of PHREEQC",
                "3895bc5caf3f843abbb6deade72c9515c2be1c2efa3dd57052d0bb10cba68996",
                "its pressure and temperature data: Appelo, Parkhurst & Post (2014), Geochim. " *
                "Cosmochim. Acta 125, 49-67, doi:10.1016/j.gca.2013.10.003",
            ),
            ThirdPartyDatabase(
                "CEMDATA18-31-03-2022-phaseVol.dat", "Cemdata18 in PHREEQC format, from Empa",
                String[],
                "bfbbfd399f5b5a025498493a18eb9b72e60c09a71e497bc51fa412de23f6279e",
                "https://www.empa.ch/web/s308/thermodynamic-data",
                "the terms of Empa's download page",
                "Lothenbach et al. (2019), Cem. Concr. Res. 115, 472-506, doi:10.1016/j.cemconres.2018.04.018",
            ),
        )
)

const DATABASE_DIR_VARIABLE = "CHEMISTRYLAB_DATABASE_DIR"

# Tests point the cache elsewhere; `nothing` means the Scratch space.
const _CACHE_OVERRIDE = Ref{Union{Nothing, String}}(nothing)

"""
    database_cache() -> String

The directory where ChemistryLab keeps the databases it has downloaded or built:
a Scratch space of the Julia depot, created on first use. Deleting it is safe;
the files are obtained again when next needed.
"""
database_cache() = something(_CACHE_OVERRIDE[], get_scratch!(@__MODULE__, "databases"))

# Obtaining a database is serialized. Independent calculations run on threads
# (the documentation computes six coupled trajectories at once), and each asks
# for its database: the memo and the announcement set below are a `Dict` and a
# `Set`, which concurrent insertions corrupt, and two threads building the same
# derived database would write the same file. The lock is reentrant because a
# derived database obtains its base through the same path.
const _DATABASE_LOCK = ReentrantLock()

# SHA-256 of a file, memoized on (path, size, mtime): the databases are read
# hundreds of times in a session and hashed once.
const _DIGESTS = Dict{Tuple{String, Int, Float64}, String}()
function _digest(path)
    st = stat(path)
    return lock(_DATABASE_LOCK) do
        get!(_DIGESTS, (abspath(path), Int(st.size), st.mtime)) do
            bytes2hex(open(sha256, path))
        end
    end
end

# Names whose foreign version has been announced once already this session.
const _ANNOUNCED = Set{String}()

# A temporary file beside `target`, so that moving it into place is atomic, and
# with a name of its own, so that two processes sharing a depot never write the
# same one.
_partial_path(target) = tempname(dirname(target); cleanup = false) * ".part"

function _local_candidates(dir, name)
    isdir(dir) || return String[]
    found = String[]
    direct = joinpath(dir, name)
    isfile(direct) && push!(found, direct)
    for entry in sort!(readdir(dir; join = true))
        isdir(entry) || continue
        candidate = joinpath(entry, name)
        isfile(candidate) && push!(found, candidate)
    end
    return found
end

# A user's own copy: the validated version if one of the candidates is it,
# otherwise the first candidate, announced once as a different version.
function _from_local_dir(db::ThirdPartyDatabase, dir)
    (dir === nothing || isempty(dir)) && return nothing
    candidates = _local_candidates(dir, db.name)
    isempty(candidates) && return nothing
    for c in candidates
        _digest(c) == db.sha256 && return c
    end
    chosen = first(candidates)
    if !(chosen in _ANNOUNCED)
        push!(_ANNOUNCED, chosen)
        @warn "`$(db.name)` in $(DATABASE_DIR_VARIABLE) is not the version ChemistryLab " *
            "was validated with; it is used, and results may differ from the documented ones." path =
            chosen found = _digest(chosen) expected = db.sha256
    end
    return chosen
end

# The cache holds either the validated version, or a version the user installed
# with `force = true`, whose digest is recorded beside it.
_accepted_path(cached) = cached * ".accepted-sha256"

function _cached(db::ThirdPartyDatabase, cache)
    cached = joinpath(cache, db.name)
    isfile(cached) || return nothing
    digest = _digest(cached)
    digest == db.sha256 && return cached
    accepted = _accepted_path(cached)
    isfile(accepted) && strip(read(accepted, String)) == digest && return cached
    # Damaged, or left by an interrupted run: remove it so the next step
    # obtains the file again.
    rm(cached; force = true)
    return nothing
end

function _download!(db::ThirdPartyDatabase, cache)
    mkpath(cache)
    target = joinpath(cache, db.name)
    attempts = String[]
    for url in db.urls
        tmp = _partial_path(target)
        try
            Downloads.download(url, tmp; timeout = 120)
            digest = bytes2hex(open(sha256, tmp))
            if digest == db.sha256
                mv(tmp, target; force = true)
                @info "ChemistryLab downloaded `$(db.name)` ($(db.title)) and verified " *
                    "its checksum." source = url license = db.license cite = db.citation stored =
                    target
                return target
            end
            push!(attempts, "$url: the content differs from the validated version (SHA-256 $digest)")
        catch err
            push!(attempts, "$url: $(sprint(showerror, err))")
        finally
            rm(tmp; force = true)
        end
    end
    throw(DatabaseUnavailable(db, attempts, cache))
end

"""
    DatabaseUnavailable <: Exception

Thrown when a database cannot be obtained. The message names the file, says
where it was looked for and what was tried, and lists the ways to provide it.
"""
struct DatabaseUnavailable <: Exception
    db::ThirdPartyDatabase
    attempts::Vector{String}
    cache::String
end

function Base.showerror(io::IO, e::DatabaseUnavailable)
    db = e.db
    dir = get(ENV, DATABASE_DIR_VARIABLE, "")
    print(io, "DatabaseUnavailable: ChemistryLab could not obtain `$(db.name)` ($(db.title)).\n")
    print(io, "  It is not in $(DATABASE_DIR_VARIABLE) ")
    print(io, isempty(dir) ? "(not set)" : "($dir)")
    print(io, ", nor in the cache ($(e.cache)).\n")
    if isempty(db.urls)
        print(io, "  ChemistryLab does not download it: its publisher distributes it through a web\n")
        print(io, "  page, under its own terms, so it has to be downloaded once by hand. To provide it, either:\n")
        print(io, "    - open $(db.page) in a web browser, download `$(db.name)`,\n")
        print(io, "      then call `install_database(\"<path of the downloaded file>\")`;\n")
    else
        print(io, "  Downloading it failed:\n")
        for a in e.attempts
            print(io, "    - ", a, "\n")
        end
        print(io, "  To provide it, either:\n")
        print(io, "    - restore network access and run the call again, or `fetch_databases()`;\n")
        print(io, "    - download it yourself from one of the addresses above, or from\n")
        print(io, "      $(db.page), then call `install_database(\"<path>\")`;\n")
    end
    print(io, "    - or set ENV[\"$(DATABASE_DIR_VARIABLE)\"] to a directory that holds it.\n")
    return print(io, "  Publisher's license: $(db.license). Cite: $(db.citation).")
end

"""
    database_path(name; download = true) -> String

Path to the database file `name` — a database obtained from its publisher (see
[`THIRD_PARTY_DATABASES`](@ref)) or one ChemistryLab builds from such a database
(see [`DERIVED_DATABASES`](@ref)) — obtaining or building it when needed.
[`datapath`](@ref) calls this for every database name, so code written
`datapath("cemdata18-thermofun.json")` needs no change.

With `download = false` nothing is downloaded, and a file absent from the local
directory and the cache throws [`DatabaseUnavailable`](@ref).
"""
function database_path(name::AbstractString; download::Bool = true)
    haskey(DERIVED_DATABASES, name) || haskey(THIRD_PARTY_DATABASES, name) || throw(
        ArgumentError(
            "`$name` is not a database ChemistryLab obtains. Known databases: " *
                join(sort!(vcat(collect(keys(THIRD_PARTY_DATABASES)), collect(keys(DERIVED_DATABASES)))), ", "),
        ),
    )
    haskey(DERIVED_DATABASES, name) && return _derived_path(DERIVED_DATABASES[name]; download)
    return _obtain(THIRD_PARTY_DATABASES[name], get(ENV, DATABASE_DIR_VARIABLE, nothing), database_cache(); download)
end

function _obtain(db::ThirdPartyDatabase, dir, cache; download::Bool = true)
    return lock(_DATABASE_LOCK) do
        local_copy = _from_local_dir(db, dir)
        local_copy === nothing || return local_copy
        cached = _cached(db, cache)
        cached === nothing || return cached
        (download && !isempty(db.urls)) || throw(DatabaseUnavailable(db, String[], cache))
        _download!(db, cache)
    end
end

is_database_name(name::AbstractString) =
    haskey(THIRD_PARTY_DATABASES, name) || haskey(DERIVED_DATABASES, name)

"""
    install_database(path; name = basename(path), force = false) -> String

Put a copy of a database the user downloaded into ChemistryLab's cache, and
return its path there. The file is identified by `name` and its SHA-256 must be
that of the version ChemistryLab was validated with; a different file is
refused unless `force = true`, in which case it is installed and a warning says
that results may differ from the documented ones.

This is how a database that ChemistryLab does not download is provided (one
whose publisher distributes it through a web page, under its own terms), and how
a machine without network access is prepared.
"""
function install_database(path::AbstractString; name::AbstractString = basename(path), force::Bool = false)
    db = get(THIRD_PARTY_DATABASES, name, nothing)
    db === nothing && throw(
        ArgumentError(
            "`$name` is not a database ChemistryLab obtains; pass `name = \"<database name>\"` " *
                "if the file was saved under another name. Known databases: " *
                join(sort!(collect(keys(THIRD_PARTY_DATABASES))), ", "),
        ),
    )
    isfile(path) || throw(ArgumentError("no such file: $path"))
    digest = bytes2hex(open(sha256, path))
    if digest != db.sha256 && !force
        throw(
            ArgumentError(
                "`$path` is not the version of `$name` ChemistryLab was validated with " *
                    "(SHA-256 $digest, expected $(db.sha256)). If it is a release you mean to " *
                    "use, call `install_database(path; force = true)`; results may then differ " *
                    "from the documented ones.",
            ),
        )
    end
    return lock(_DATABASE_LOCK) do
        cache = database_cache()
        mkpath(cache)
        target = joinpath(cache, name)
        tmp = _partial_path(target)
        try
            cp(path, tmp; force = true)
            mv(tmp, target; force = true)
        finally
            rm(tmp; force = true)
        end
        accepted = _accepted_path(target)
        if digest == db.sha256
            rm(accepted; force = true)
        else
            write(accepted, digest)
            @warn "`$name` installed in a version ChemistryLab was not validated with; results may " *
                "differ from the documented ones." found = digest expected = db.sha256
        end
        target
    end
end

"""
    fetch_databases(; derived = true) -> Vector{NamedTuple}

Obtain in advance every database that can be downloaded, and build the derived
ones, so that later work needs no network access: before a documentation build,
on a machine that will go offline, or in continuous integration. A database that
has to be downloaded by hand is reported as missing, not as an error.

Returns one `(; name, path, status)` per database, `status` being `:ready` or
`:missing`.
"""
function fetch_databases(; derived::Bool = true)
    out = NamedTuple{(:name, :path, :status), Tuple{String, Union{Nothing, String}, Symbol}}[]
    for name in sort!(collect(keys(THIRD_PARTY_DATABASES)))
        db = THIRD_PARTY_DATABASES[name]
        path = try
            _obtain(db, get(ENV, DATABASE_DIR_VARIABLE, nothing), database_cache(); download = !isempty(db.urls))
        catch err
            err isa DatabaseUnavailable && isempty(db.urls) || rethrow()
            nothing
        end
        push!(out, (; name, path, status = path === nothing ? :missing : :ready))
    end
    if derived
        for name in sort!(collect(keys(DERIVED_DATABASES)))
            push!(out, (; name, path = database_path(name), status = :ready))
        end
    end
    return out
end

"""
    database_info([io]) -> Vector{NamedTuple}

Where each database ChemistryLab knows currently resolves, without downloading
or building anything: `:local` (in `CHEMISTRYLAB_DATABASE_DIR`), `:cached`,
`:downloadable`, `:manual` (to be downloaded by hand) or, for a derived
database, `:buildable`. Printed as a table, each database with its source and
the license its publisher states, and returned as one
`(; name, status, path, source, license)` per database.
"""
function database_info(io::IO = stdout)
    dir = get(ENV, DATABASE_DIR_VARIABLE, nothing)
    cache = database_cache()
    rows = NamedTuple{(:name, :status, :path, :source, :license), Tuple{String, Symbol, Union{Nothing, String}, String, String}}[]
    for name in sort!(collect(keys(THIRD_PARTY_DATABASES)))
        db = THIRD_PARTY_DATABASES[name]
        local_copy, cached = lock(_DATABASE_LOCK) do
            local_copy = try
                _from_local_dir(db, dir)
            catch
                nothing
            end
            local_copy, (local_copy === nothing ? _cached(db, cache) : nothing)
        end
        status, path = if local_copy !== nothing
            (:local, local_copy)
        elseif cached !== nothing
            (:cached, cached)
        else
            (isempty(db.urls) ? :manual : :downloadable, nothing)
        end
        push!(rows, (; name, status, path, source = db.page, license = db.license))
    end
    for name in sort!(collect(keys(DERIVED_DATABASES)))
        d = DERIVED_DATABASES[name]
        path = joinpath(cache, "derived", name)
        status = isfile(path) ? :cached : :buildable
        push!(rows, (; name, status, path = isfile(path) ? path : nothing, source = "built from $(d.base) and $(d.description)", license = "that of $(d.base)"))
    end
    for r in rows
        println(io, rpad(r.name, 36), rpad(string(r.status), 14), r.source)
        println(io, " "^50, "license: ", r.license)
    end
    return rows
end
