# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

"""
    datapath(parts::AbstractString...) -> String

Absolute path to a data file of ChemistryLab. Called without argument, returns
the package's `data/` directory.

A thermodynamic database is named by its file name alone,
`datapath("cemdata18-thermofun.json")`, and resolved by
[`database_path`](@ref): from `ENV["CHEMISTRYLAB_DATABASE_DIR"]` if it is set
and holds the file, else from the package's cache, else downloaded from its
publisher on first use, or built from it for a derived database. Any other name
is a file under `data/`: the solid-solution and gel models, the published values
of `data/literature/`, the experimental datasets.

This is the recommended way to name a database, because it depends neither on
the working directory nor on where the file happens to be stored: a script
written this way runs identically from the package root, from an editor whose
REPL started elsewhere, inside a documentation build and on a machine that has
never used ChemistryLab before.

# Examples

```julia
substances = build_species(datapath("cemdata18-thermofun.json"))
ss_phases  = build_solid_solutions(datapath("solid_solutions.toml"), dict)
datapath("experimental", "README.md")     # subdirectories of data/ work too
```

See also [`database_info`](@ref), [`fetch_databases`](@ref),
[`read_thermofun_database`](@ref), [`build_species`](@ref).
"""
function datapath(parts::AbstractString...)
    length(parts) == 1 && is_database_name(only(parts)) && return database_path(only(parts))
    root = pkgdir(@__MODULE__)
    root === nothing && error(
        "cannot locate the ChemistryLab package directory, so its " *
            "`data/` files cannot be resolved. Pass an explicit path instead.",
    )
    return joinpath(root, "data", parts...)
end

"""
    resolve_data_path(path::AbstractString) -> String

Resolve `path` to an existing file, trying in order:

 1. `path` as given — relative to the working directory, or absolute;
 2. `datapath(path)` — the same relative path under `data/`;
 3. `datapath(basename(path))` — a data file or a database of that name,
    whatever directory prefix the caller wrote (a database is obtained as
    [`datapath`](@ref) describes);
 4. `joinpath(pkgdir(ChemistryLab), path)` — relative to the package root.

The working directory comes first, so a call that already resolves keeps
resolving to exactly the same file: the fallbacks can only turn a failure into a
success, never change an existing answer. And steps 2–3 can only succeed for a
name that *is* one of the package's data files or databases, so a mistyped path
to a file of one's own still fails loudly instead of silently loading something
else.

Throws `ArgumentError` listing the data files and databases when nothing
matches, and [`DatabaseUnavailable`](@ref) when the name is a database that
cannot be obtained.
"""
function resolve_data_path(path::AbstractString)
    isfile(path) && return String(path)

    for candidate in (datapath(path), datapath(basename(path)))
        isfile(candidate) && return candidate
    end

    root = pkgdir(@__MODULE__)
    if root !== nothing
        candidate = joinpath(root, path)
        isfile(candidate) && return candidate
    end

    files = try
        sort!(readdir(datapath()))
    catch
        String[]
    end
    databases = sort!(vcat(collect(keys(THIRD_PARTY_DATABASES)), collect(keys(DERIVED_DATABASES))))
    throw(
        ArgumentError(
            "no such data file: \"$path\". It was looked for relative to the " *
                "working directory ($(pwd())), then among ChemistryLab's data " *
                "files ($(datapath())) and databases. Data files: " *
                "$(join(files, ", ")). Databases: $(join(databases, ", ")). " *
                "Use `datapath(\"<name>\")` to name either from anywhere.",
        ),
    )
end

"""
    display_data_path(path::AbstractString) -> String

Short, machine-independent label for `path`, meant for the banner a reader sees
rather than for opening a file.

A database is shown by its file name (`cemdata18-thermofun.json`), wherever the
cache or the user's directory put it; another file inside the package is shown
relative to the package root (`data/solid_solutions.toml`); anything else is
shown unchanged. This
matters because these banners are captured verbatim into the documentation:
printing the absolute path would bake the build machine's directories
(`/home/runner/work/...`) into every page that loads a database.

Whether a file is inside the package is decided by a prefix test and not by
asking `relpath` and reading `..` off its answer. `relpath` compares two paths
component by component, and on Windows two paths on **different drives** share no
component at all: it returns a chain of `..` and the target, which starts with
`..` only by accident, or a path that starts with the drive letter and does not.
A GitHub Windows runner checks out on one drive and puts `tempdir()` on another,
so a file plainly outside the package could come back rewritten. The prefix test
has no such case.
"""
function display_data_path(path::AbstractString)
    # A database obtained from its publisher lives in the cache or in the user's
    # own directory, both machine-specific: it is shown by its name.
    name = basename(path)
    is_database_name(name) && !startswith(abspath(path), abspath(pwd())) && return name
    root = pkgdir(@__MODULE__)
    if root !== nothing
        absolute = abspath(path)
        sep = Base.Filesystem.path_separator
        prefix = endswith(root, sep) ? root : root * sep
        # The separator is part of the prefix on purpose: without it, a sibling
        # directory whose name merely begins with the package's would match.
        startswith(absolute, prefix) && return relpath(absolute, root)
    end
    return String(path)
end
