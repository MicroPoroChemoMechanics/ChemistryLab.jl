# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

using TOML

"""
    build_pitzer_parameters(path; format = :auto) -> PitzerParameters

Read a Pitzer interaction-parameter set from a TOML file (`format = :toml`) or
from the `PITZER` block of a PHREEQC-format database (`format = :phreeqc`);
`:auto` takes the TOML reader for a `.toml` file and the PHREEQC one otherwise.

Nothing about this is automatic: the file has to be named, which is the whole
point. `data/pitzer-reardon1990.toml` ships with the package and is reached as

```julia
p = build_pitzer_parameters(datapath("pitzer-reardon1990.toml"))
```

# File format

```toml
[meta]
name = "..."
temperature_K = 298.15
speciation = "dissociated"        # informational; read by the caller, not here

[[binary]]                        # one per cation-anion pair
cation = "Na+"
anion  = "Cl-"
beta0  = 0.0765
beta1  = 0.2664
beta2  = 0.0
Cphi   = 0.0013
origin = "fitted"                 # or "estimated:HSO4-"

[[theta]]                         # like-charge pair, either order
i = "Na+"
j = "K+"
value = -0.012

[[psi]]                           # triplet, any order
i = "Na+"
j = "Cl-"
k = "SO4-2"
value = 0.0014

[[lambda]]                        # neutral with an ion
neutral = "H2CO3@"
ion = "Na+"
value = 0.1
origin = "fitted"
```

Temperature terms (see [`PitzerParameters`](@ref)) are given as
`temperature = { beta0 = [A₁, A₂, A₃, A₄, A₅], Cphi = […] }` in a `[[binary]]`
entry and as `temperature = [A₁, …]` in a `[[theta]]`, `[[psi]]` or `[[lambda]]`
one, fewer than five meaning the rest are zero.

`beta1`, `beta2` and `Cphi` default to zero within a `[[binary]]` entry;
`beta0` does not, because an entry without it describes nothing. The shape
constants `alpha1`, `alpha1_22`, `alpha2` and `b` may be given under `[meta]`
and otherwise take the conventional values documented in
[`PitzerParameters`](@ref).

# What it refuses

A `[[binary]]` entry without `cation`, `anion` or `beta0`, and a duplicate pair
— silently keeping the last of two conflicting entries would be worse than
stopping. Species names are **not** checked against any database here: whether
the set covers the system is decided when the model is attached to it, which is
where the error can name the pairs that are missing.

# The PHREEQC format

```
PITZER
-B0
    Na+    Cl-    0.0765   A₁  A₂  A₃  A₄  A₅
-C0
    ...
-THETA
    Ca+2   Na+    0.07
-LAMDA
    CO2    Na+    0.1
-PSI
    Ca+2   Na+    Cl-     -0.007
```

`-B0`, `-B1`, `-B2` and `-C0` (`Cφ`) take a cation and an anion in either
order, `-THETA` two ions of like sign, `-LAMDA` a neutral species and an ion,
`-PSI` a triplet, each followed by `A₀` and up to five temperature coefficients
[ParkhurstAppelo2013](@cite). Species are named as PHREEQC writes them, with `@`
appended to a neutral species, which is how this package names them. The
identifiers this model has no counterpart for (`-ZETA`, `-ETA`, `-MU`, `-ALPHAS`,
`-APHI`, which the water model provides) are skipped, and a warning names them;
`-use_etheta` is the model's own keyword, `etheta`. **No database is stored in
this package**: this reads the file the caller has (`datapath("pitzer.dat")`
obtains PHREEQC's own), and every `origin` records its name and a truncated
SHA-256 of its contents.

See also: [`PitzerParameters`](@ref), [`PitzerActivityModel`](@ref),
[`pitzer_origin`](@ref).
"""
function build_pitzer_parameters(path::AbstractString; format::Symbol = :auto)
    format in (:auto, :toml, :phreeqc) || throw(
        ArgumentError("format must be :auto, :toml or :phreeqc; got :$format"),
    )
    fmt = format === :auto ? (lowercase(splitext(path)[2]) == ".toml" ? :toml : :phreeqc) : format
    return fmt === :toml ? _pitzer_from_toml(path) : _pitzer_from_phreeqc(path)
end

function _pitzer_from_toml(toml_file::AbstractString)
    path = resolve_data_path(toml_file)
    raw = TOML.parsefile(path)
    meta = get(raw, "meta", Dict{String, Any}())

    beta0 = Dict{Tuple{String, String}, Float64}()
    beta1 = Dict{Tuple{String, String}, Float64}()
    beta2 = Dict{Tuple{String, String}, Float64}()
    Cphi = Dict{Tuple{String, String}, Float64}()
    origin = Dict{Tuple{String, String}, String}()
    temperature = Dict{Symbol, Dict{Any, Vector{Float64}}}()
    tterm!(kind, key, v) = (get!(temperature, kind, Dict{Any, Vector{Float64}}())[key] = Float64.(v))

    for (i, e) in enumerate(get(raw, "binary", Any[]))
        for key in ("cation", "anion", "beta0")
            haskey(e, key) || error(
                "build_pitzer_parameters: [[binary]] entry $i in $path has no \"$key\"."
            )
        end
        k = (String(e["cation"]), String(e["anion"]))
        haskey(beta0, k) && error(
            "build_pitzer_parameters: $path gives the pair $(k[1])/$(k[2]) twice."
        )
        beta0[k] = float(e["beta0"])
        beta1[k] = float(get(e, "beta1", 0.0))
        beta2[k] = float(get(e, "beta2", 0.0))
        Cphi[k] = float(get(e, "Cphi", 0.0))
        origin[k] = String(get(e, "origin", "unrecorded"))
        for (name, v) in get(e, "temperature", Dict{String, Any}())
            tterm!(Symbol(name), k, v)
        end
    end

    theta = Dict{Tuple{String, String}, Float64}()
    for (i, e) in enumerate(get(raw, "theta", Any[]))
        haskey(e, "i") && haskey(e, "j") && haskey(e, "value") || error(
            "build_pitzer_parameters: [[theta]] entry $i in $path needs i, j and value."
        )
        theta[(String(e["i"]), String(e["j"]))] = float(e["value"])
        haskey(e, "temperature") && tterm!(:theta, (String(e["i"]), String(e["j"])), e["temperature"])
    end

    psi = Dict{Tuple{String, String, String}, Float64}()
    for (i, e) in enumerate(get(raw, "psi", Any[]))
        all(haskey(e, k) for k in ("i", "j", "k", "value")) || error(
            "build_pitzer_parameters: [[psi]] entry $i in $path needs i, j, k and value."
        )
        psi[(String(e["i"]), String(e["j"]), String(e["k"]))] = float(e["value"])
        haskey(e, "temperature") && tterm!(:psi, (String(e["i"]), String(e["j"]), String(e["k"])), e["temperature"])
    end

    lambda = Dict{Tuple{String, String}, Float64}()
    for (i, e) in enumerate(get(raw, "lambda", Any[]))
        all(haskey(e, k) for k in ("neutral", "ion", "value")) || error(
            "build_pitzer_parameters: [[lambda]] entry $i in $path needs " *
                "neutral, ion and value."
        )
        k = (String(e["neutral"]), String(e["ion"]))
        lambda[k] = float(e["value"])
        haskey(e, "origin") && (origin[k] = String(e["origin"]))
        haskey(e, "temperature") && tterm!(:lambda, k, e["temperature"])
    end

    return PitzerParameters(;
        beta0 = beta0, beta1 = beta1, beta2 = beta2, Cphi = Cphi,
        theta = theta, psi = psi, lambda = lambda,
        alpha1 = float(get(meta, "alpha1", 2.0)),
        alpha1_22 = float(get(meta, "alpha1_22", 1.4)),
        alpha2 = float(get(meta, "alpha2", 12.0)),
        b = float(get(meta, "b", 1.2)),
        origin = origin,
        temperature = temperature,
    )
end

# The options of a PITZER block, by their spellings, and those that set no
# coefficient this model has.
const _PITZER_OPTIONS = Dict(
    "-b0" => :beta0, "-b1" => :beta1, "-b2" => :beta2, "-c0" => :Cphi,
    "-theta" => :theta, "-lamda" => :lambda, "-lambda" => :lambda, "-psi" => :psi,
)
const _PITZER_SWITCHES = ("-use_etheta", "-macinnes", "-redox")

# The PITZER blocks of a PHREEQC database. The names are read by the rule of
# every reader (`_phreeqc_name`, then `_solute_symbol` for a neutral species,
# which gets `@`); a cation–anion pair is keyed cation first, a neutral species
# with an ion neutral first.
function _pitzer_from_phreeqc(path::AbstractString)
    isfile(path) || throw(ArgumentError("no such database: $path"))
    blocks = [b for b in phreeqc_blocks(path) if b.keyword == "PITZER"]
    isempty(blocks) && throw(ArgumentError("no PITZER block in $path."))
    src = _source_tag(path)

    tables = Dict(k => Dict{Any, Float64}() for k in (:beta0, :beta1, :beta2, :Cphi, :theta, :psi, :lambda))
    temperature = Dict{Symbol, Dict{Any, Vector{Float64}}}()
    skipped = Set{String}()
    pairs_seen = Set{Tuple{String, String}}()
    for b in blocks
        current = nothing
        for (_, line) in b.lines
            if startswith(line, '-')
                id = lowercase(first(split(line)))
                current = get(_PITZER_OPTIONS, id, nothing)
                current === nothing && !(id in _PITZER_SWITCHES) && push!(skipped, id)
                continue
            end
            current === nothing && continue
            fields = split(line)
            nsp = current === :psi ? 3 : 2
            length(fields) > nsp || throw(ArgumentError("$path: a $current line has no coefficient: \"$line\""))
            names = _phreeqc_name.(fields[1:nsp])
            coeffs = parse.(Float64, fields[(nsp + 1):end])
            length(coeffs) <= 6 || throw(ArgumentError("$path: more than six coefficients on \"$line\""))
            key = if current in _PITZER_PAIRS
                _cation_anion(names, current, path, line)
            elseif current === :lambda
                _neutral_first(names, path, line)
            else
                Tuple(_solute_symbol.(names))
            end
            tables[current][key] = coeffs[1]
            current in _PITZER_PAIRS && push!(pairs_seen, key)
            length(coeffs) > 1 && (get!(temperature, current, Dict{Any, Vector{Float64}}())[key] = coeffs[2:end])
        end
    end
    isempty(skipped) || @warn "build_pitzer_parameters: $(basename(path)) carries $(join(sort!(collect(skipped)), ", ")), which this model has no counterpart for; they are skipped." maxlog = 1

    # A pair described by β¹, β² or Cφ without β⁰ has β⁰ = 0, as PHREEQC takes it.
    for k in pairs_seen
        haskey(tables[:beta0], k) || (tables[:beta0][k] = 0.0)
    end
    P2, P3 = Tuple{String, String}, Tuple{String, String, String}
    conv2(d) = Dict{P2, Float64}(P2(k) => v for (k, v) in d)
    return PitzerParameters(;
        beta0 = conv2(tables[:beta0]), beta1 = conv2(tables[:beta1]), beta2 = conv2(tables[:beta2]),
        Cphi = conv2(tables[:Cphi]), theta = conv2(tables[:theta]),
        psi = Dict{P3, Float64}(P3(k) => v for (k, v) in tables[:psi]),
        lambda = conv2(tables[:lambda]),
        origin = Dict{P2, String}(k => src for k in keys(tables[:beta0])),
        temperature = temperature,
    )
end

const _PITZER_PAIRS = (:beta0, :beta1, :beta2, :Cphi)

# A cation–anion pair of a PITZER block, keyed cation first.
function _cation_anion(names, table, path, line)
    z = _name_charge.(names)
    z[1] * z[2] < 0 || throw(ArgumentError("$path: $table needs a cation and an anion: \"$line\""))
    return z[1] > 0 ? (names[1], names[2]) : (names[2], names[1])
end

# A neutral species with an ion (or another neutral species), keyed neutral first.
function _neutral_first(names, path, line)
    z = _name_charge.(names)
    count(iszero, z) >= 1 || throw(ArgumentError("$path: a lambda line needs a neutral species: \"$line\""))
    a, b = iszero(z[1]) ? (names[1], names[2]) : (names[2], names[1])
    return (_solute_symbol(a), _solute_symbol(b))
end
