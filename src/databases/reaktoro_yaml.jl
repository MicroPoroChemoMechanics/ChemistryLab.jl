# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

# ── Reading a database of Reaktoro ───────────────────────────────────────────
#
# A database of Reaktoro is a YAML file whose `Species` map gives, for each
# species, its elements, charge and aggregate state, and its standard
# thermodynamic model with the parameters of that model, in SI units. Those
# models of formation properties that the package computes are read: `HKF` (the
# aqueous species), `MaierKelley` and `HollandPowell` (a heat capacity
# a + bT + c/T² + d/√T, the latter with the modified Tait equation of state), and
# the solvent, water. The others are refused by name.
#
# Only the subset of YAML these files use is read: maps, by indentation, of
# scalars. A line of any other form is an error naming it.

# The maps of a YAML file, by indentation; a scalar is kept as a string.
function _yaml_maps(path)
    root = OrderedDict{String, Any}()
    stack = Tuple{Int, OrderedDict{String, Any}}[(-1, root)]
    for (n, raw) in enumerate(eachline(path))
        line = rstrip(first(split(raw, " #"; limit = 2)))
        s = strip(line)
        (isempty(s) || startswith(s, "#") || s == "---") && continue
        indent = length(line) - length(lstrip(line))
        m = match(r"^([^:]+?):\s*(.*)$", s)
        m === nothing && throw(ArgumentError("$(basename(path)):$n: `$s` is not a `key: value` line; only maps of scalars are read"))
        while first(last(stack)) >= indent
            pop!(stack)
        end
        key, value = String(strip(m.captures[1])), String(strip(m.captures[2]))
        parent = last(last(stack))
        if isempty(value)
            child = OrderedDict{String, Any}()
            parent[key] = child
            push!(stack, (indent, child))
        else
            # A list, in flow (`[a, b]`) or block (`- a`) style; a negative
            # number is a scalar.
            (startswith(value, "[") || value == "-" || startswith(value, "- ")) && throw(
                ArgumentError("$(basename(path)):$n: `$s` holds a list; only maps of scalars are read")
            )
            parent[key] = strip(value, ['"', '\''])
        end
    end
    return root
end

# `3:Fe 2:Al` as a composition.
function _reaktoro_elements(s::AbstractString)
    out = Dict{Symbol, Float64}()
    for w in split(s)
        m = match(r"^([\d.eE+-]+):(\S+)$", w)
        m === nothing && throw(ArgumentError("`$w` is not a `coefficient:element` pair"))
        out[Symbol(m.captures[2])] = parse(Float64, m.captures[1])
    end
    return out
end

const _REAKTORO_STATES = Dict(
    "Aqueous" => AS_AQUEOUS, "Solid" => AS_CRYSTAL, "Gas" => AS_GAS, "Liquid" => AS_LIQUID,
)

# The ChemistryLab symbol of a species of Reaktoro: a neutral solute gets `@` in
# place of `(aq)`, water is the solvent.
function _reaktoro_symbol(name, state, charge)
    state == AS_AQUEOUS || return name
    base = replace(name, r"\(aq\)$" => "")
    base == "H2O" && return "H2O@"
    return iszero(charge) ? base * "@" : base
end

"""
    read_reaktoro_database(path) -> (df_elements, df_substances, df_reactions)

Read a database of Reaktoro, a YAML file whose `Species` map gives each species
with its standard thermodynamic model, into the three tables
[`read_thermofun_database`](@ref) returns for a ThermoFun file. `df_substances`
has one row per species: `symbol` (ChemistryLab's), `name`, `formula`,
`aggregate_state`, `class`, `charge`, `atoms`, `standard_model` and its
`parameters` (in SI units, as the file gives them), and `gauge`; the file's
reactions are not read, and `df_reactions` is empty.

[`build_species`](@ref) builds from `df_substances` the models the package
computes: `HKF` (the revised Helgeson-Kirkham-Flowers equations of the aqueous
species, as `solute_hkf88_reaktoro` computes them), `MaierKelley` (a heat
capacity a + bT + c/T² and a constant volume), `HollandPowell` (a heat capacity
a + bT + c/T² + d/√T and the volume of the modified Tait equation of state of
[HollandPowell2011](@cite), integrated from the reference pressure; a constant
volume without `kappa0`), and water (`WaterHKF`, `WaterHGK`,
`WaterWagnerPruss`), the solvent of the package's own equation of state.
Reaktoro integrates the volume of `HollandPowell` from zero pressure, as the
article writes it: its energies exceed these by the molar volume times one bar,
a few J/mol. A species of `HollandPowell` with a Landau transition (`Tcr`,
`Smax`, `Vmax`), for which the article states no equations of its own, is used
at its reference state only. Any other model is left out, with a warning naming
it.

No database of Reaktoro is distributed with the package or downloaded by it: the
data of its files come from several sources, under their own terms, and a copy
is installed by hand ([`install_database`](@ref)).
"""
function read_reaktoro_database(path::AbstractString)
    file = basename(path)
    root = _yaml_maps(path)
    haskey(root, "Species") || throw(ArgumentError("$file: no `Species` map; not a database of Reaktoro"))
    rows = NamedTuple[]
    notes = String[]
    for (name, entry) in root["Species"]
        entry isa AbstractDict || continue
        model = get(entry, "StandardThermoModel", nothing)
        if !(model isa AbstractDict) || length(model) != 1
            push!(notes, "$file: $name has no single standard model; not read")
            continue
        end
        kind, p = only(model)
        state = get(_REAKTORO_STATES, get(entry, "AggregateState", ""), AS_UNDEF)
        charge = parse(Float64, get(entry, "Charge", "0"))
        sym = _reaktoro_symbol(name, state, charge)
        cls = state == AS_AQUEOUS ? (sym == "H2O@" ? SC_AQSOLVENT : SC_AQSOLUTE) : (state == AS_GAS ? SC_GASFLUID : SC_COMPONENT)
        parameters = Dict{String, Float64}()
        if p isa AbstractDict
            for (k, v) in p
                x = tryparse(Float64, v)
                x === nothing ? push!(notes, "$file: parameter `$k` of $name is not a number; not read") : (parameters[k] = x)
            end
        end
        push!(
            rows, (;
                symbol = sym, name = String(name), formula = String(get(entry, "Formula", name)),
                aggregate_state = state, class = cls, charge, atoms = _reaktoro_elements(get(entry, "Elements", "")),
                standard_model = String(kind), parameters, gauge = "formation from the elements",
            ),
        )
    end
    df_substances = DataFrame(Tables.dictrowtable(rows))
    for (key, value) in ("format" => "reaktoro", "path" => String(path), "source" => file, "notes" => notes)
        metadata!(df_substances, key, value; style = :note)
    end
    df_elements = DataFrame(symbol = sort!(unique(string(e) for r in rows for e in keys(r.atoms))))
    return df_elements, df_substances, DataFrame()
end

const _REAKTORO_MODELS = ("HKF", "MaierKelley", "HollandPowell", "WaterHKF", "WaterHGK", "WaterWagnerPruss")

# One species of a table of a database of Reaktoro.
function _reaktoro_species(row)
    atoms, charge, p, kind = row.atoms, row.charge, row.parameters, row.standard_model
    whole = all(isinteger, values(atoms))
    s = Species(
        whole ? Dict{Symbol, Int}(e => Int(n) for (e, n) in atoms) : atoms, isinteger(charge) ? Int(charge) : charge;
        name = row.name, symbol = row.symbol, aggregate_state = row.aggregate_state, class = row.class,
    )
    num(k) = get(p, k, 0.0)
    s.Tref = T_STANDARD_Q
    s.Pref = P_STANDARD_Q
    values0 = [
        :Cp⁰ => missing,
        :ΔₐH⁰ => num("Hf") * u"J/mol", :S⁰ => num("Sr") * u"J/(mol*K)", :ΔₐG⁰ => num("Gf") * u"J/mol",
        :V⁰ => haskey(p, "Vr") ? num("Vr") * u"m^3/mol" : missing,
    ]
    if kind == "HKF"
        s[:thermo_method] = "solute_hkf88_reaktoro"
        hkf = [Symbol(k) => num(k) for k in ("a1", "a2", "a3", "a4", "c1", "c2", "wref")]
        s[:thermo_params] = [hkf; :z => charge; values0; :T => s.Tref; :P => s.Pref]
    elseif kind in ("MaierKelley", "HollandPowell")
        s[:thermo_method] = "cp_ft_equation"
        cp = [
            :a₀ => num("a") * u"J/(mol*K)", :a₁ => num("b") * u"J/(mol*K^2)", :a₂ => num("c") * u"J*K/mol",
            :a₃ => num("d") * u"J/(mol*K^(1//2))",
        ]
        s[:thermo_params] = [cp; values0; :T => s.Tref; :P => s.Pref]
        tait = kind == "HollandPowell" && haskey(p, "kappa0")
        haskey(p, "Vr") && !tait && (s[:V_method] = "mv_constant")
        if kind == "HollandPowell" && any(k -> haskey(p, k), ("Tcr", "Smax", "Vmax"))
            s[:refused_method] = "HollandPowell with a Landau transition"
            _reference_state_only!(s)
        elseif tait
            missed = filter(k -> !haskey(p, k), ["Vr", "Sr", "alpha0", "kappa0p", "numatoms"])
            isempty(missed) || throw(ArgumentError("$(row.name): `HollandPowell` with `kappa0` needs $(join(missed, ", "))"))
            _holland_powell11!(s, p)
        end
    else
        row.symbol == "H2O@" || throw(ArgumentError("$(row.name): the model `$kind` is that of water"))
        _reaktoro_water!(s, p)
    end
    s[:gauge] = row.gauge
    return s
end

# The species of a table of a database of Reaktoro, those named (by the file's
# name or ChemistryLab's symbol) or all; a species of a model the package does
# not compute is left out with a warning.
function _build_reaktoro_species(df::AbstractDataFrame, names = nothing)
    out = Species[]
    refused = Tuple{String, String}[]
    for r in eachrow(df)
        names === nothing || r.name in names || r.symbol in names || continue
        if !(r.standard_model in _REAKTORO_MODELS)
            push!(refused, (r.name, r.standard_model))
            continue
        end
        push!(out, _reaktoro_species(r))
    end
    isempty(refused) || @warn "build_species: $(length(refused)) species left out, their standard models not computed by ChemistryLab: " *
        join(("$s (`$m`)" for (s, m) in refused), ", ")
    return out
end

# Water: the package's equation of state, referred to the triple point as
# SUPCRT92 refers it. The values the file gives there must be the same.
function _reaktoro_water!(s, p)
    tp = _WATER_TRIPLE_POINT
    for (key, ref) in (("Gtr", tp.G), ("Htr", tp.H), ("Str", tp.S), ("Ttr", tp.T))
        haskey(p, key) || continue
        v = p[key]
        abs(v - ref) <= 1.0e-6 * abs(ref) || throw(
            ArgumentError("the water of this database is referred to $key = $v, where the package's equation of state takes $ref")
        )
    end
    M = ustrip(us"kg/mol", s[:M])
    s[:V_method] = "water_eos"
    s[:thermo_params] = [
        :Cp⁰ => missing, :ΔₐH⁰ => missing, :S⁰ => missing, :ΔₐG⁰ => missing,
        :V⁰ => M / _hgk_density(T_STANDARD, P_STANDARD) * u"m^3/mol", :T => s.Tref, :P => s.Pref,
    ]
    return s
end

# ── Holland and Powell (2011): the modified Tait equation of state ───────────
#
# A mineral of the `HollandPowell` model with a bulk modulus κ₀ (`kappa0`), its
# first and second pressure derivatives κ₀′ and κ₀″ (`kappa0p`, `kappa0pp`), a
# thermal expansion α₀ (`alpha0`) and the number n of atoms of its formula
# (`numatoms`) has the volume of the modified Tait equation with a thermal
# pressure [HollandPowell2011; Eqs. 3, 11–13](@cite),
#
#     V/V₀ = 1 − a (1 − (1 + b (P − P_th))^(−c)),
#     a = (1 + κ₀′)/(1 + κ₀′ + κ₀κ₀″),  b = κ₀′/κ₀ − κ₀″/(1 + κ₀′),
#     c = (1 + κ₀′ + κ₀κ₀″)/(κ₀′² + κ₀′ − κ₀κ₀″),
#     P_th = α₀ κ₀ (θ/ξ₀) [1/(e^(θ/T) − 1) − 1/(e^(θ/T₀) − 1)],
#     ξ₀ = u₀² e^(u₀)/(e^(u₀) − 1)²,  u₀ = θ/T₀,
#
# θ the Einstein temperature of the mineral, 10636/(S°/n + 6.44) K with its
# entropy S° in J/(mol K) [HollandPowell2011; p. 346](@cite), both constants read
# from `data/literature/HollandPowell2011.json`. Without κ₀″ the article's own
# value −κ₀′/κ₀ is taken. The Gibbs energy gains ∫ V dP from the reference
# pressure. The article's Eq. 13 integrates from zero pressure, as Reaktoro
# does: at the reference state Reaktoro's energies exceed the file's by the
# volume times one bar, 3.7 J/mol for calcite, where these keep the file's.
const _HP11 = let q(name) = ustrip(literature_value("HollandPowell2011", name))
    (; θ = q("einstein_temperature_coefficient"), s = q("einstein_temperature_entropy_offset"))
end

# The parameters of the equation from those of the file.
function _tait_parameters(p, Tr)
    κ0, κ1 = p["kappa0"], p["kappa0p"]
    κ2 = get(p, "kappa0pp", -κ1 / κ0)
    θ = _HP11.θ / (p["Sr"] / p["numatoms"] + _HP11.s)
    u0 = θ / Tr
    return (;
        V0 = p["Vr"], a = (1 + κ1) / (1 + κ1 + κ0 * κ2), b = κ1 / κ0 - κ2 / (1 + κ1),
        c = (1 + κ1 + κ0 * κ2) / (κ1^2 + κ1 - κ0 * κ2), ακ = p["alpha0"] * κ0, θ, u0,
        ξ0 = u0^2 * exp(u0) / expm1(u0)^2,
    )
end

# ∫₀^P V dP of the modified Tait equation at T [HollandPowell2011; Eq. 13].
function _tait_integral(q, T, P)
    Pth = q.ακ * q.θ / q.ξ0 * (1 / expm1(q.θ / T) - 1 / expm1(q.u0))
    return q.V0 * ((1 - q.a) * P + q.a * ((1 - q.b * Pth)^(1 - q.c) - (1 + q.b * (P - Pth))^(1 - q.c)) / (q.b * (q.c - 1)))
end

function _holland_powell11!(s, p)
    Tr, Pr = ustrip(us"K", s.Tref), ustrip(us"Pa", s.Pref)
    q = _tait_parameters(p, Tr)
    return _add_excess!(s, (T, P) -> _tait_integral(q, T, P) - _tait_integral(q, T, Pr))
end
