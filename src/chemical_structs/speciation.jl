# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

using OrderedCollections

"""
    same_components(::AbstractVector{<:AbstractSpecies}) -> Function

Return the function to extract components from species vectors.

Returns `atoms_charge` for Species vectors, `oxides_charge` for CemSpecies vectors.
"""
same_components(::AbstractVector{<:AbstractSpecies}) = atoms_charge
same_components(::AbstractVector{<:CemSpecies}) = oxides_charge

"""
    item_order(::AbstractVector{<:AbstractSpecies}) -> Vector{Symbol}

Return the ordering vector for components.

Returns `ATOMIC_ORDER` for Species vectors, `OXIDE_ORDER` for CemSpecies vectors.
"""
item_order(::AbstractVector{<:AbstractSpecies}) = ATOMIC_ORDER
item_order(::AbstractVector{<:CemSpecies}) = OXIDE_ORDER

"""
    union_atoms(atom_dicts::AbstractVector{<:AbstractDict}, order_vec=ATOMIC_ORDER) -> Vector{Symbol}

Compute the union of all keys from dictionaries, sorted by a given order.

# Arguments

  - `atom_dicts`: vector of dictionaries (e.g., atomic compositions).
  - `order_vec`: ordering vector for sorting keys (default: ATOMIC_ORDER).

# Returns

  - Sorted vector of unique symbols appearing in any dictionary.

# Examples

```jldoctest
julia> d1 = OrderedDict(:H => 2, :O => 1);

julia> d2 = OrderedDict(:C => 1, :O => 2);

julia> union_atoms([d1, d2], ATOMIC_ORDER)
3-element Vector{Symbol}:
 :C
 :H
 :O
```
"""
function union_atoms(atom_vecs::AbstractVector{<:AbstractVector}, order_vec = ATOMIC_ORDER)
    function sortfunc(k)
        idx = findfirst(==(k), order_vec)
        return isnothing(idx) ? max(1, length(order_vec) - 1) : idx
    end
    return sort!(collect(union(atom_vecs...)); by = sortfunc)
end

union_atoms(atom_dicts::AbstractVector{<:AbstractDict}, order_vec = ATOMIC_ORDER) =
    union_atoms(collect.(keys.(atom_dicts)), order_vec)

union_atoms(species_list::AbstractVector{<:AbstractSpecies}, order_vec = ATOMIC_ORDER) =
    union_atoms(collect.(keys.(atoms.(species_list))), order_vec)

function idx_speciation(
        species_list, atoms_list::AbstractVector{Symbol};
        # EVERY state, written as `instances` rather than listed. The list used
        # to name all four members there were, so its intent was "do not filter
        # on state at all" -- and the day a fifth was added it silently became a
        # filter, excluding it. `instances` cannot go stale that way.
        aggregate_state = collect(instances(AggregateState)),
        # The classes are NOT `instances`, and the difference is deliberate:
        # `SC_SSENDMEMBER` is absent because a solid-solution end member enters a
        # system through its phase and not as a species of its own, and
        # `SC_SURFCOMPLEX` is absent for exactly the same reason -- a surface
        # species enters through its `SiteFamily`. Five of seven, on purpose.
        class = [SC_AQSOLUTE, SC_AQSOLVENT, SC_COMPONENT, SC_GASFLUID, SC_UNDEF],
        exclude_species = [],
        include_species = [],
    )
    mask1 = [all(k .∈ Ref(atoms_list)) for k in keys.(atoms.(species_list))]
    mask2 = ChemistryLab.aggregate_state.(species_list) .∈ Ref(aggregate_state)
    mask3 = ChemistryLab.class.(species_list) .∈ Ref(class)
    mask4 = (species_list .∉ Ref(exclude_species)) .&&
        (symbol.(species_list) .∉ Ref(exclude_species))
    mask5 = (species_list .∈ Ref(include_species)) .||
        (symbol.(species_list) .∈ Ref(include_species))
    return (mask1 .&& mask2 .&& mask3 .&& mask4) .|| mask5
end

"""
    speciation(species_list, atoms_list::AbstractVector{Symbol}; kwargs...) -> Vector
    speciation(species_list, seeds::AbstractVector{<:AbstractSpecies}; keep_seeds = true, kwargs...) -> Vector
    speciation(species_list, seed_symbols::AbstractVector{<:AbstractString}; keep_seeds = true, kwargs...) -> Vector

Filter a list of species to those whose atomic composition is a subset of the
given atoms, or of the atoms found in the seed species.

# Arguments

  - `species_list`: full list of species to filter.
  - `atoms_list`: vector of atom symbols that define the chemical space.
  - `seeds`: seed species whose union of atoms defines the chemical space.
  - `seed_symbols`: symbols of seed species, looked up in `species_list`; a
    symbol found there under no species raises an `ArgumentError`.

# Keyword arguments

  - `aggregate_state`: filter by aggregate state (default: all states).
  - `class`: filter by species class (default: `SC_AQSOLUTE`, `SC_AQSOLVENT`,
    `SC_COMPONENT`, `SC_GASFLUID` and `SC_UNDEF`). Solid-solution end-members
    and surface complexes are left out, because they enter a system through
    their phase or their site family, not as species of their own.
  - `exclude_species`: species (or symbols) to always exclude.
  - `include_species`: species (or symbols) to always include regardless of
    composition.
  - `keep_seeds`: whether the seeds are always included, whatever the filters
    (default), or go through `aggregate_state`, `class` and `exclude_species`
    like every other species. With the default, a solid seed is in the result
    even under `aggregate_state = [AS_AQUEOUS]`, which is what building a
    system from its phases needs.

# Returns

  - A new `Vector` holding the matching species, in the order of `species_list`.
    It is not a view: it can be filtered or extended without touching the list.

Two results join without duplicates with Base `union`, since species compare
by formula, aggregate state, class and symbol: `union(speciation(subs,
["C3S"]), speciation(subs, ["Portlandite"]))`.

# Examples

```julia
julia> all_species = [Species("H2O"), Species("OH-"), Species("H+"), Species("Ca+2"), Species("CaCO3")];

julia> speciation(all_species, [:H, :O])
3-element Vector{Species{Int64}}:
 H₂O {H₂O} [H₂O ◆ H2O]
 OH⁻ {OH⁻} [OH⁻ ◆ OH-]
 H⁺ {H⁺} [H⁺ ◆ H+]
```
"""
function speciation(
        species_list, atoms_list::AbstractVector{Symbol}; kwargs...
    )
    return species_list[idx_speciation(species_list, atoms_list; kwargs...)]
end

function speciation(
        species_list, seeds::AbstractVector{<:AbstractSpecies};
        include_species = [], keep_seeds::Bool = true, kwargs...
    )
    return speciation(
        species_list, union_atoms(seeds);
        include_species = keep_seeds ? vcat(include_species, seeds) : include_species, kwargs...
    )
end

function speciation(
        species_list, seed_symbols::AbstractVector{<:AbstractString}; kwargs...
    )
    symbols = symbol.(species_list)
    absent = setdiff(seed_symbols, symbols)
    isempty(absent) || throw(
        ArgumentError(
            "speciation: no species $(join(absent, ", ")) in the list. A seed that is " *
                "not there would leave its elements out of the chemical space without " *
                "a word; check the spelling or the database.",
        )
    )
    return speciation(species_list, species_list[symbols .∈ Ref(seed_symbols)]; kwargs...)
end
