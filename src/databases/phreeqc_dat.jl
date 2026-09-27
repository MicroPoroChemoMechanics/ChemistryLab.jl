# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

using JSON
using TOML

"""
    parse_reaction_stoich_cemdata(reaction_line::AbstractString) -> (Vector, String, String)

Parse a reaction line from CEMDATA format and extract stoichiometric information.

# Arguments

  - `reaction_line`: reaction string in CEMDATA format, optionally with a comment after '#'.

# Returns

  - `reactants`: vector of dictionaries with "symbol" and "coefficient" keys.
  - `modified_equation`: equation string with added "@" markers for aqueous species.
  - `comment`: extracted comment string (empty if none).

The function automatically adds "@" suffixes to aqueous species without explicit charges
(except for the first reactant).
"""
function parse_reaction_stoich_cemdata(reaction_line::AbstractString)
    # Split the equation and the comment
    equation_parts = split(reaction_line, '#')
    equation = strip(equation_parts[1])
    comment = length(equation_parts) > 1 ? strip(equation_parts[2]) : ""

    parts = split(equation, "=")
    if length(parts) != 2
        @warn "Equation format is incorrect: $equation"
        return [], equation, comment
    end

    modified_equation_parts = String[]
    reactants = []

    for (side_index, side) in enumerate(parts)
        tokens = split(side)
        modified_tokens = []
        for (token_index, token) in enumerate(tokens)
            if token == "+"
                push!(modified_tokens, token)
                continue
            end
            m = match(r"^([+-]?\d*\.?\d+)?(.+?)([\+\-]\d*\.?\d*)?$", token)
            if m !== nothing
                coeff_str = m.captures[1]
                base_symbol = m.captures[2]
                charge = m.captures[3]
                if isnothing(charge)
                    charge = ""
                end
                sp = base_symbol * charge
                coeff = if isnothing(coeff_str) || isempty(coeff_str)
                    1.0
                else
                    parse(Float64, coeff_str)
                end

                # Add "@" to aqueous species without charge
                if isempty(charge) && !(side_index == 1 && token_index == 1)
                    sp *= "@"
                    token = coeff_str === nothing ? sp : coeff_str * sp
                end

                push!(modified_tokens, token)

                # Add to reactants/products
                push!(
                    reactants,
                    Dict("symbol" => sp, "coefficient" => side_index == 1 ? -coeff : coeff),
                )
            else
                push!(modified_tokens, token)
            end
        end
        push!(modified_equation_parts, join(modified_tokens, " "))
    end
    modified_equation = join(modified_equation_parts, " = ")

    return reactants, modified_equation, comment
end

"""
    parse_float_array(line::AbstractString) -> Vector{Float64}

Parse a line containing space-separated floats, skipping the first token and any comments.

# Arguments

  - `line`: input string with format "keyword value1 value2 ...".

# Returns

  - Vector of successfully parsed Float64 values.

# Examples

```julia
julia> parse_float_array("-analytical_expression 1.5 2.3 4.7")
3-element Vector{Float64}:
 1.5
 2.3
 4.7

julia> parse_float_array("-log_K 5.2 # comment")
1-element Vector{Float64}:
 5.2
```
"""
function parse_float_array(line)
    parts = split(line)
    if length(parts) < 2
        return Float64[]
    end
    float_parts = Float64[]
    for part in parts[2:end]
        if !startswith(part, "#")
            try
                push!(float_parts, parse(Float64, part))
            catch e
                # @warn "Could not parse '$part' as Float64, skipping."
            end
        end
    end
    return float_parts
end

# The option of a PHREEQC data line, by the rules of the PHREEQC manual
# (Parkhurst & Appelo 2013, "Description of Data Input"): case does not matter,
# and each option has documented spellings, with or without the leading dash
# and abbreviated. Only the spellings the manual lists are accepted, so that an
# abbreviation shared by two options is never guessed.
const _PHREEQC_OPTIONS = Dict(
    :log_k => ("log_k", "logk", "-log_k", "-logk", "-l"),
    :analytic => ("-analytic", "analytic", "analytical_expression", "-analytical_expression", "a_e", "ae", "-a_e", "-ae", "-a"),
    :gamma => ("-gamma", "-g"),
    :vm => ("-vm",),
)
function _phreeqc_option(token::AbstractString)
    t = lowercase(token)
    for (opt, spellings) in _PHREEQC_OPTIONS
        t in spellings && return opt
    end
    return nothing
end

# The numbers following an option, up to a comment or a unit.
function _phreeqc_numbers(parts)
    out = Float64[]
    for part in parts
        startswith(part, "#") && break
        v = tryparse(Float64, part)
        v === nothing && break
        push!(out, v)
    end
    return out
end

"""
    parse_phases(dat_content::AbstractString) -> Dict{String,Any}

Extract phase information from PHREEQC .dat file content.

# Arguments

  - `dat_content`: full text content of a PHREEQC .dat file.

# Returns

  - Dictionary mapping phase names to their properties: `equation`, `reactants`,
    `logKr`, `analytical_expression` and `molar_volume`.

Parses the PHASES section. The options are recognized as PHREEQC recognizes
them, whatever their case and in every spelling the manual documents (`log_k`,
`-log_K`, `-l`; `-analytic`, `-analytical_expression`, `-a_e`; `-Vm`). `-Vm` is
the **molar volume of the phase**, in cm³/mol, and is stored as such, never as a
volume of reaction.
"""
function parse_phases(dat_content)
    phases = Dict{String, Any}()
    in_phases = false
    current_phase = nothing

    for line in eachline(IOBuffer(dat_content))
        line = strip(line)
        if startswith(line, "PHASES")
            in_phases = true
            continue
        elseif in_phases && !isempty(line) && !startswith(line, "#")
            parts = split(line)
            opt = _phreeqc_option(first(parts))
            if opt === nothing && !occursin("=", line)
                # A new phase: its name, alone on its line.
                current_phase = Dict{String, Any}("symbol" => first(parts))
                phases[first(parts)] = current_phase
            elseif opt === nothing && current_phase !== nothing
                reactants, equation, comment = parse_reaction_stoich_cemdata(line)
                current_phase["equation"] = equation
                current_phase["reactants"] = reactants
                isempty(comment) || (current_phase["comment"] = comment)
            elseif current_phase === nothing
                continue
            elseif opt === :log_k
                v = _phreeqc_numbers(parts[2:end])
                if isempty(v)
                    @warn "Could not parse log_K value for phase $(current_phase["symbol"]), skipping."
                else
                    current_phase["logKr"] = Dict("values" => [first(v)], "errors" => [2])
                end
            elseif opt === :analytic
                analytical_expression = _phreeqc_numbers(parts[2:end])
                # coef of log10 in .dat becomes a coef of log in .json
                if length(analytical_expression) > 3
                    analytical_expression[4] /= log(10)
                end
                current_phase["analytical_expression"] = analytical_expression
            elseif opt === :vm
                v = _phreeqc_numbers(parts[2:end])
                isempty(v) || (current_phase["molar_volume"] = first(v))
            end
        end
    end

    return phases
end

"""
    phreeqc_gamma_parameters(path) -> Dict{String, Tuple{Float64, Float64}}

The WATEQ activity-coefficient parameters `(å, b)` that the `-gamma` option of a
PHREEQC database gives its aqueous species, by species symbol, for
[`TruesdellJonesActivityModel`](@ref).

Every `SOLUTION_SPECIES` block is read, the master species as well as the
others. The species a reaction defines is the first one to the right of its
equal sign, as the PHREEQC manual requires, and its name becomes a ChemistryLab
symbol by the rule the readers of this package use: a name ending in a charge is
kept, any other gets `@` appended (`CO2` is `CO2@`), and `H2O` is the solvent.

A second `-gamma` for a species replaces the first, as PHREEQC reads the options
of a species in order. A `-gamma` line whose numbers do not parse is reported
with its line number and skipped.
"""
function phreeqc_gamma_parameters(path::AbstractString)
    params = Dict{String, Tuple{Float64, Float64}}()
    where_defined = Dict{String, Int}()
    in_species = false
    current = nothing
    for (lineno, raw) in enumerate(eachline(resolve_data_path(path)))
        line = strip(first(split(raw, '#')))
        isempty(line) && continue
        if occursin(r"^[A-Z][A-Z_0-9]*$", line)
            in_species = line == "SOLUTION_SPECIES"
            current = nothing
            continue
        end
        in_species || continue
        parts = split(line)
        opt = _phreeqc_option(first(parts))
        if opt === nothing && occursin("=", line)
            rhs = strip(split(line, "="; limit = 2)[2])
            first_product = first(split(rhs, r"\s+\+\s+"))
            name = replace(strip(first_product), r"^[0-9.]+" => "")
            current = _phreeqc_symbol(name)
        elseif opt === :gamma && current !== nothing
            v = _phreeqc_numbers(parts[2:end])
            if length(v) < 2
                @warn "phreeqc_gamma_parameters: line $lineno of $(basename(path)) has a -gamma that does not parse; skipped" line
                continue
            end
            # PHREEQC reads the options of a species in order, so a second
            # `-gamma` replaces the first; `phreeqc.dat` itself refines Na+ that
            # way ("halite solubility").
            haskey(params, current) &&
                @debug "phreeqc_gamma_parameters: line $lineno replaces the -gamma of $current given at line $(where_defined[current])"
            params[current] = (v[1], v[2])
            where_defined[current] = lineno
        end
    end
    return params
end

function _phreeqc_symbol(name::AbstractString)
    name == "e-" && return "Zz"
    name == "H2O" && return "H2O@"
    return occursin(r"[+-][0-9]*$", name) ? String(name) : String(name) * "@"
end

"""
    extract_primary_species(file_path::AbstractString) -> DataFrame

Extract primary aqueous species from a PHREEQC database file.

# Arguments

  - `file_path`: path to PHREEQC .dat file.

# Returns

  - DataFrame with columns: species, symbol, formula, aggregate_state, atoms, charge, gamma.

Parses the SOLUTION_SPECIES section to extract master species and their properties.
The "Zz" charge placeholder is handled specially. Gamma coefficients for activity
models are extracted from "-gamma" lines.
"""
function extract_primary_species(file_path)
    lines = readlines(resolve_data_path(file_path))

    start_idx = 0
    end_idx = 0
    in_primary_section = false

    for (i, line) in enumerate(lines)
        stripped = strip(line)

        if startswith(stripped, "SOLUTION_SPECIES")
            start_idx = i + 1
            in_primary_section = true

            while start_idx <= length(lines)
                next_line = strip(lines[start_idx])
                if startswith(next_line, "# PMATCH MASTER SPECIES") ||
                        occursin("=", next_line)
                    break
                end
                start_idx += 1
            end

            if startswith(strip(lines[start_idx]), "# PMATCH MASTER SPECIES")
                start_idx += 1
            end
        end

        if in_primary_section && startswith(stripped, "# PMATCH SECONDARY MASTER SPECIES")
            end_idx = i - 1
            break
        end
    end

    if in_primary_section && end_idx == 0
        end_idx = length(lines)
    end

    species_data = []

    for i in start_idx:end_idx
        line = strip(lines[i])

        if occursin("=", line)
            parts = split(line, "=")
            current_species = strip(parts[1])

            if current_species == "e-"
                symbol = "Zz"
            elseif occursin(r"[\+\-]\d*$", current_species)
                symbol = current_species
            else
                symbol = current_species * "@"
            end

            push!(species_data, (species = current_species, symbol = symbol, gamma = Float64[]))
        end

        if startswith(line, "-gamma") && !isempty(species_data)
            parts = split(line)
            gamma_values = Float64[]
            for val in parts[2:end]
                num = tryparse(Float64, val)
                if num !== nothing
                    push!(gamma_values, num)
                end
            end

            if !isempty(gamma_values)
                last_entry = species_data[end]
                species_data[end] = (
                    species = last_entry.species, symbol = last_entry.symbol, gamma = gamma_values,
                )
            end
        end
    end

    df = DataFrame(species_data)
    df.symbol = String.(df.symbol)
    df.formula .= df.symbol
    df.aggregate_state .= "AS_AQUEOUS"
    df.atoms .= parse_formula.(df.symbol)
    df.charge .= extract_charge.(df.symbol)
    df[df.symbol .== "Zz", :species] .= "Zz"
    df[df.symbol .== "Zz", :formula] .= "Zz"
    df[df.symbol .== "Zz", :charge] .= 1
    return df[sortperm(df.symbol .== "Zz"), :]
end
