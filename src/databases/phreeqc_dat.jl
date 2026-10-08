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

function _phreeqc_symbol(name::AbstractString)
    name == "e-" && return "Zz"
    name == "H2O" && return "H2O@"
    return occursin(r"[+-][0-9]*$", name) ? String(name) : String(name) * "@"
end
