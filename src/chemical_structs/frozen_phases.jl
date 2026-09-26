# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

# ── A solid solution, read and then frozen ───────────────────────────────────
#
# Some models need a solid whose composition an equilibrium has decided, but
# which must not take part in the next equilibrium. The surface model of the
# C-S-H is the case in point. Its sites bind calcium and alkalis, and the CSHQ
# solid solution holds the same calcium and alkalis in its end members: carrying
# both in one system counts them twice. The published surface models avoid this
# with a C-S-H of fixed composition that only carries the sites. What those
# models leave open is the composition and the amount of that C-S-H, which a
# first equilibrium with CSHQ can supply.
#
# So two stages. The first equilibrium decides the solid solution. The second
# takes it out of the system: its elements are set aside, except the alkalis,
# which the surface model describes on its own and which go back to the
# solution. Nothing else moves, and the conservation of every element across the
# two stages is exact.

"""
    solid_solution_totals(state::ChemicalState, name::AbstractString) -> NamedTuple

What the solid solution `name` of `state` holds:

  - `members`: the amount of each end member, mol, by symbol;
  - `amount`: their sum, mol;
  - `elements`: the moles of each element over the end members;
  - `mass`: the mass of the phase, kg, from the molar masses of the end members.

Ratios such as Ca/Si are quotients of `elements`. Each end member counts in its
own formula unit, so `amount` depends on how the database writes them, and
`elements` does not.

The element type of the result is that of the state, so a state carrying dual
numbers gives derivatives of these totals.
"""
function solid_solution_totals(state::ChemicalState, name::AbstractString)
    cs = state.system
    group = cs.ss_groups[_solid_solution_index(cs, name)]
    n = [ustrip(us"mol", state.n[i]) for i in group]
    sp = cs.species[group]
    members = OrderedDict(symbol(s) => x for (s, x) in zip(sp, n))
    els = unique(e for s in sp for e in keys(atoms(s)))
    elements = OrderedDict(e => sum(get(atoms(s), e, 0) * x for (s, x) in zip(sp, n)) for e in els)
    mass = sum(_molar_mass_si.(sp) .* n)
    return (; members, amount = sum(n), elements, mass)
end

function _solid_solution_index(cs::ChemicalSystem, name::AbstractString)
    names = cs.solid_solutions === nothing ? String[] : [ChemistryLab.name(ss) for ss in cs.solid_solutions]
    k = findfirst(==(name), names)
    k === nothing && throw(
        ArgumentError(
            "no solid solution \"$name\" in this system" *
                (isempty(names) ? ", which declares none." : "; it declares $(join(names, ", ")).")
        )
    )
    return k
end

"""
    freeze_solid_solution(state, name, target::ChemicalSystem;
                          release = (:Na, :K), buffer = nothing) -> NamedTuple

The first state of a second stage in which the solid solution `name` of `state`
no longer reacts: its elements are set aside, and `target` sees the rest.

Every species of `target` found in the system of `state`, under the same symbol
and with the same composition, starts at the amount `state` holds. The end
members of `name` are frozen, and `target` must not contain them. The elements
listed in `release` leave the frozen solid and return to the solution as their
monovalent cation, balanced by as much hydroxide: an end member written as
`((NaOH)2.5SiO2H2O)0.2` gives back 0.5 NaOH, and its silica and water stay
frozen. `release = ()` freezes everything.

Refused, with the names, rather than computed wrong:

  - a species that holds matter in `state` but is absent from `target` (its
    elements would vanish between the stages), or present under the same symbol
    with another composition;
  - an end member of `name` among the species of `target` (it could form again
    from the frozen elements, which are no longer in the budget);
  - a released element without its cation, or a solution without hydroxide, in
    `target`;
  - `buffer`, when given, absent from `state`. A frozen composition is a sound
    approximation while a coexisting phase fixes the activities of the elements
    it shares with the solid solution. For a C-S-H that phase is portlandite;
    without it the composition of the gel moves with the solution, and freezing
    it is not a model of anything. The check is at the first stage; whether the
    buffer survives the second is a result of that stage and is checked on it.

Returns `(state, frozen)`. `state` is the `ChemicalState` of `target`, at the
temperature and pressure of the first. `frozen` holds the `elements` and the
`mass` (kg) set aside, the `members` they came from, and the moles of each
`released` element. With `A₁` and `A₂` the element matrices of the two systems,
`A₂ n₂ + frozen = A₁ n₁` holds element by element; what is added to `state`
afterwards, a salt for instance, is the second stage's business.
"""
function freeze_solid_solution(
        state::ChemicalState, name::AbstractString, target::ChemicalSystem;
        release = (:Na, :K), buffer::Union{Nothing, AbstractString} = nothing,
    )
    cs = state.system
    totals = solid_solution_totals(state, name)
    members = collect(keys(totals.members))
    source = Dict(symbol(s) => i for (i, s) in enumerate(cs.species))
    into = Dict(symbol(s) => j for (j, s) in enumerate(target.species))

    if buffer !== nothing
        ib = get(source, buffer, nothing)
        ib === nothing && throw(ArgumentError("the buffer \"$buffer\" is not a species of the first system."))
        _certainly_zero(state.n[ib]) && throw(
            ArgumentError(
                "\"$buffer\" is absent from the first stage, so nothing fixes the composition " *
                    "of \"$name\" once the solution changes, and freezing it is not a model. " *
                    "Describe the uptake by the solid solution itself instead (for chloride in " *
                    "C-S-H, the CSHQ_Cl phase of the database cemdata18-chloride.json)."
            )
        )
    end

    back = filter(m -> haskey(into, m), members)
    isempty(back) || throw(
        ArgumentError(
            "the second system contains $(join(back, ", ")), end members of the frozen " *
                "\"$name\"; they could form again from elements no longer in the budget."
        )
    )
    lost = String[]
    for (i, s) in enumerate(cs.species)
        symbol(s) in members && continue
        j = get(into, symbol(s), nothing)
        if j === nothing
            _certainly_zero(state.n[i]) || push!(lost, symbol(s))
        elseif atoms_charge(target.species[j]) != atoms_charge(s)
            throw(
                ArgumentError(
                    "\"$(symbol(s))\" has composition $(phreeqc(formula(s))) in the first system " *
                        "and $(phreeqc(formula(target.species[j]))) in the second; a symbol is a " *
                        "label, and the amount cannot follow it."
                )
            )
        end
    end
    isempty(lost) || throw(
        ArgumentError(
            "$(join(lost, ", ")) hold matter in the first stage but are not species of the " *
                "second, so their elements would vanish between the stages. Add them to it."
        )
    )

    # The released elements: those the phase carries, decided on its formulas
    # and not on the amounts, so that the result has one structure whatever
    # the numbers.
    freed = [e for e in release if haskey(totals.elements, e)]
    cation = Dict(e => _ion_index(target, Dict(e => 1, :Zz => 1), "$(e)⁺") for e in freed)
    hydroxide = isempty(freed) ? 0 : _ion_index(target, Dict(:O => 1, :H => 1, :Zz => -1), "OH⁻")
    released = OrderedDict(e => totals.elements[e] for e in freed)

    n1 = [ustrip(us"mol", x) for x in state.n]
    n2 = zeros(eltype(n1), length(target.species))
    for (j, s) in enumerate(target.species)
        i = get(source, symbol(s), nothing)
        i === nothing || (n2[j] = n1[i])
    end
    for e in freed
        n2[cation[e]] += released[e]
        n2[hydroxide] += released[e]
    end

    elements = copy(totals.elements)
    mass = totals.mass
    for e in freed
        elements[e] -= released[e]
        elements[:O] -= released[e]
        elements[:H] -= released[e]
        mass -= released[e] * (_molar_mass_si(target.species[cation[e]]) + _molar_mass_si(target.species[hydroxide]))
    end
    st = ChemicalState(target, n2 .* u"mol"; T = temperature(state), P = pressure(state))
    return (; state = st, frozen = (; elements, mass, members = totals.members, released))
end

# The index of the species of `cs` with this exact composition and charge.
function _ion_index(cs::ChemicalSystem, composition, label)
    js = findall(s -> Dict(atoms_charge(s)) == composition, cs.species)
    isempty(js) && throw(
        ArgumentError("the second system has no $label to receive what the frozen phase releases.")
    )
    return only(js)
end

# An amount known to be zero. The solver returns an absent phase as an exact
# zero; a symbolic amount is never known to be zero.
_certainly_zero(x) = (z = iszero(ustrip(x)); z isa Bool && z)
