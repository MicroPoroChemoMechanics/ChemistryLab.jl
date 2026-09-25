# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

using JSON

# ── Aqueous species and real solids come from a shipped database ─────────────
#
# Not from numbers typed into a test file. A standard Gibbs energy recalled by
# hand drifts from the database it came from, and it drifts silently: the water
# written into four of these files read −237181 J/mol against slop98's
# −237183.0, and portlandite read −897010 against CEMDATA18's −897013.0. Neither
# is large enough to fail a test and both are wrong.
#
# It also removes a whole class of question from a reviewer's plate. A number in
# a test file has to be checked against something; a species pulled from the
# database the package ships is already the thing it would be checked against.
#
# What is still constructed by hand, and legitimately so, is the **surface
# species**: no database carries `XsOH2+`, and its standard energy is built to
# *be* the log K under test. That is the oracle, derived rather than looked up,
# and every file says so where it does it.
#
# `test/activities.jl` established this pattern; these files now follow it.

if !@isdefined(_REFERENCE_DATABASES)
    const _REFERENCE_DATABASES = Dict{Symbol, Dict{String, Any}}()
end

"""
    reference_species(symbol; db = :slop98) -> AbstractSpecies

A species from one of the databases the package ships, by its symbol.

`db` is `:slop98` (SUPCRT92/slop98 inorganic — the aqueous species) or
`:cemdata18` (the cement phases). Each database is read once per session.
"""
function reference_species(sym::AbstractString; db::Symbol = :slop98)
    file = db === :slop98 ? "slop98-inorganic-thermofun.json" :
        db === :cemdata18 ? "cemdata18-thermofun.json" :
        error("unknown reference database $db; use :slop98 or :cemdata18")
    cache = get!(_REFERENCE_DATABASES, db) do
        Dict{String, Any}(
            symbol(s) => s for s in build_species(datapath(file); verbose = false)
        )
    end
    haskey(cache, sym) || error("$sym is not in $file")
    return cache[sym]
end

"""
    reference_species(symbols; db = :slop98) -> Vector

Several at once, in the order given.
"""
reference_species(syms; db::Symbol = :slop98) =
    [reference_species(s; db) for s in syms]

"""
    moles_of_water(mass = 1.0u"kg") -> Real

How many moles of `H2O@` weigh `mass`, from the solvent species' own molar mass.

`55.5` is what this number is usually written as, and writing it is a small
mistake with a measurable consequence: `55.5 mol` of water weighs 0.99983 kg,
not 1 kg, so every molality computed against it is off by 1.7 × 10⁻⁴ relative.
Matching Reaktoro's water basis exactly is what took one of these comparisons
from 5.5e-5 to 4.1e-9.
"""
moles_of_water(mass = 1.0u"kg") =
    ustrip(us"mol", mass / reference_species("H2O@")[:M])

# ── Oracle fixtures ──────────────────────────────────────────────────────────

"""
    reference_oracle(name) -> NamedTuple

A cross-code fixture from `test/reference/<name>.json`, as a nested
`NamedTuple` so it reads like the Julia literal it replaces: `f.points[1].pH`.

These files are **generated**, by the scripts beside them, and were previously
pasted into the test files as Julia constants. They are data, and they now live
where data lives. Three things that cost time are gone with the paste:

  - a transcription step, which is the error class these generators exist to
    remove in the first place;
  - a float formatter inside each generator, because Python writes `6.0e-09`
    and `5e-06` where Runic — which gates CI — wants `6.0e-9` and `5.0e-6`, so
    a data file had to be emitted in the *formatter's* dialect;
  - a provenance header written as comments, unreadable by anything. Versions,
    the database and its md5 are now fields, and a test can assert on them.

Regenerate with the script named in the fixture's own `generator` field.
"""
function reference_oracle(name::AbstractString)
    path = joinpath(@__DIR__, "reference", name * ".json")
    isfile(path) || error("no oracle fixture at $path; regenerate it with its script")
    return _as_namedtuple(JSON.parsefile(path; dicttype = Dict{String, Any}))
end

"""
    gems_bdot() -> Float64

The B-dot a GEM-Selektor run of a CEMDATA18 Portland cement implies, identified
from the activity coefficients it printed for |z| = 1 and 2 at its ionic
strength (fixture `gems_cemdata18_portland`), the limiting-law slope taken from
the same two: about 0.0976. `test/aqueous_properties.jl` checks that it
predicts the coefficient of |z| = 3.
"""
function gems_bdot()
    g = reference_oracle("gems_cemdata18_portland")
    D = (log10(g.gamma.z1) - log10(g.gamma.z2)) / 3
    return (log10(g.gamma.z1) + D) / g.ionic_strength_mol_per_kg
end

_as_namedtuple(x::AbstractDict) =
    (; (Symbol(k) => _as_namedtuple(v) for (k, v) in x)...)
_as_namedtuple(x::AbstractVector) = [_as_namedtuple(v) for v in x]
_as_namedtuple(x) = x
