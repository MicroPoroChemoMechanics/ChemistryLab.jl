# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

using LinearAlgebra
using OrderedCollections

"""
    struct ChemicalSystem{T<:AbstractSpecies, R<:AbstractReaction, C, S, SS} <: AbstractVector{T}

An immutable, fully typed collection of chemical species and reactions
with derived index structures and stoichiometric matrices.

Immutability guarantees that all derived fields (`dict_species`, `dict_reactions`,
index vectors, `CSM`, `SM`) remain consistent with `species` and `reactions`
throughout the lifetime of the object. To modify the system, use `merge` to
construct a new `ChemicalSystem`.

# Fields

  - `species`: ordered list of all species.
  - `dict_species`: fast O(1) lookup by species symbol.
  - `idx_aqueous`, `idx_crystal`, `idx_gas`: indices by aggregate state.
  - `idx_solutes`, `idx_solvent`, `idx_components`, `idx_gasfluid`: indices by class.
  - `reactions`: ordered list of all reactions.
  - `dict_reactions`: fast O(1) lookup by reaction symbol.
  - `CSM`: canonical stoichiometric matrix.
  - `SM`: stoichiometric matrix with respect to primaries.
  - `solid_solutions`: `Nothing` when no solid solutions are present, or a concrete
    `Vector{<:AbstractSolidSolutionPhase}` describing each solid-solution phase and
    its end-members. Populated via the `solid_solutions` keyword constructor.
  - `ss_groups`: for each solid solution, the indices of its end-members in `species`.
  - `idx_ssendmembers`: union of all end-member indices (flattened `ss_groups`).
  - `idx_surface`: indices of species in `AS_SURFACE`, i.e. bound to a site.
  - `site_families`: `Nothing` when no surface is declared, or a concrete
    `Vector{<:SiteFamily}`. Populated through the `site_families` keyword.
  - `site_groups`: for each family, the indices of its members in `species`, the
    **free site first** — the order the site mixing and the solver's reference
    member both rely on.
  - `idx_kinetic`: indices of kinetic species (empty when none declared).
"""
struct ChemicalSystem{T <: AbstractSpecies, R <: AbstractReaction, C, S, SS, SF} <:
    AbstractVector{T}
    species::Vector{T}
    dict_species::Dict{String, T}               # fast O(1) lookup by symbol

    # Indices by aggregate_state
    idx_aqueous::Vector{Int}
    idx_crystal::Vector{Int}
    idx_gas::Vector{Int}
    idx_surface::Vector{Int}

    # Indices by class
    idx_solutes::Vector{Int}
    idx_solvent::Vector{Int}
    idx_components::Vector{Int}
    idx_gasfluid::Vector{Int}

    reactions::Vector{R}
    dict_reactions::Dict{String, R}             # fast O(1) lookup by reaction symbol

    CSM::C                                      # canonical stoichiometric matrix — typed for performance
    SM::S                                       # stoichiometric matrix w.r.t. primaries — typed for performance

    # Solid solutions — SS = Nothing (no SS) or Vector{<:AbstractSolidSolutionPhase}
    solid_solutions::SS
    ss_groups::Vector{Vector{Int}}              # per-SS end-member indices
    idx_ssendmembers::Vector{Int}               # all end-member indices (flattened)

    # Surface site families — SF = Nothing, or Vector{<:SiteFamily}
    site_families::SF
    site_groups::Vector{Vector{Int}}            # per-family member indices, free site first

    idx_kinetic::Vector{Int}                    # kinetic species indices (empty if none)
end

# ── Constructors ──────────────────────────────────────────────────────────────

# An empty collection means "no solid solutions", and has to be accepted as
# such. Julia gives an empty *filtered* comprehension the element type `Any`,
# so `[x for x in build_solid_solutions(...) if x.name == "CSHQ"]` that selects
# nothing is a `Vector{Any}`; annotating the keyword as
# `AbstractVector{<:AbstractSolidSolutionPhase}` turned that into a `TypeError`
# on the keyword rather than a system without solid solutions, which hides the
# real problem — that the phase was skipped — behind a type error.
function _normalize_solid_solutions(ss)
    ss === nothing && return nothing
    isempty(ss) && return nothing
    all(x -> x isa AbstractSolidSolutionPhase, ss) || throw(
        ArgumentError(
            "solid_solutions must hold `AbstractSolidSolutionPhase` values; got " *
                "element types $(unique(typeof.(ss))).",
        )
    )
    return collect(AbstractSolidSolutionPhase, ss)
end


"""
    _resolve_site_families(site_families, species, idx_surface) -> (families, groups)

Resolve declared [`SiteFamily`](@ref) objects against the species list.

Returns `(nothing, Vector{Int}[])` when none is declared, so a system without a
surface is byte-identical to what it was before surfaces existed.

Everything refused here is refused rather than discovered later, because the
alternative is a wrong number rather than an error:

  - a member that is not in the species list — the family would share a budget
    with a species the system cannot see;
  - **an `AS_SURFACE` species belonging to no family** — it would carry a site
    pseudo-element into the conservation matrix, and so a row, with nothing
    mixing on it. This is the surface counterpart of the check that already
    refuses two solid solutions sharing a composition;
  - two families sharing one pseudo-element — one symbol, one budget, one
    family, or the site balance silently merges them;
  - a member whose **identity** does not survive the lookup: the system's
    species of that name must have the family's formula, must carry the
    family's site symbol, and must be `AS_SURFACE`;
  - a species claimed by **two** families.

# Why the last two are not redundant with the third

They were once believed to be. The reasoning was that two families can only
share a member if they share a site symbol, that a species carries exactly one,
and that a shared symbol is refused above — so no input could reach the case.

The reasoning skips a step. Members are matched to the system by `symbol`, and a
symbol is a label: the species the lookup lands on need not be the species the
family validated. Two families with **different** site symbols and the same
member labels therefore resolved to the same indices, one family's conservation
row disappeared, and nothing said so. The same gap let a family hold `AS_SURFACE`
copies while the system kept the caller's unqualified originals, leaving
`idx_surface` empty with the site mixing still running.
"""
function _resolve_site_families(site_families, species, idx_surface)
    if site_families === nothing || isempty(site_families)
        isempty(idx_surface) || throw(
            ArgumentError(
                "species $(join([symbol(species[i]) for i in idx_surface], ", ")) are " *
                    "in AS_SURFACE but no site family is declared. A surface species " *
                    "carries a site pseudo-element, which becomes a conservation row; " *
                    "without a family nothing mixes on it. Pass `site_families = [...]`.",
            )
        )
        return nothing, Vector{Int}[]
    end

    all(f -> f isa SiteFamily, site_families) || throw(
        ArgumentError(
            "site_families must hold `SiteFamily` values; got element types " *
                "$(unique(typeof.(site_families))).",
        )
    )

    by_symbol = Dict(symbol(sp) => i for (i, sp) in enumerate(species))
    groups = Vector{Int}[]
    seen_sites = Dict{Symbol, String}()
    seen_members = Dict{Int, String}()

    for f in site_families
        if haskey(seen_sites, f.site)
            throw(
                ArgumentError(
                    "SiteFamily \"$(f.name)\" and \"$(seen_sites[f.site])\" both use " *
                        ":$(f.site). One pseudo-element is one site budget, so two " *
                        "families sharing it would share a conservation row.",
                )
            )
        end
        seen_sites[f.site] = f.name

        group = Int[]
        for sp in site_members(f)
            i = get(by_symbol, symbol(sp), nothing)
            i === nothing && throw(
                ArgumentError(
                    "SiteFamily \"$(f.name)\": member \"$(symbol(sp))\" is not in the " *
                        "species list. Add it to the species vector first.",
                )
            )

            # The label found a species. Whether it found the RIGHT one is a
            # separate question, and it was not being asked.
            #
            # A family validates its own members and qualifies copies of them as
            # `AS_SURFACE`; the system keeps whatever the caller passed, and the
            # two are matched by `symbol` alone. So a lookup can land on a
            # species that merely shares a name — a different formula, a
            # different site symbol, or the caller's unqualified original while
            # the family holds a qualified copy. Each produces a system that is
            # plausible and wrong: in the last case `idx_surface` comes back
            # empty, so the site mixing and the phase accounting disagree about
            # which species are on a surface, silently.
            got = species[i]
            formula(got) == formula(sp) || throw(
                ArgumentError(
                    "SiteFamily \"$(f.name)\": member \"$(symbol(sp))\" has formula " *
                        "$(phreeqc(formula(sp))), but the species of that name in this " *
                        "system has formula $(phreeqc(formula(got))). A symbol is a " *
                        "label; two different species must not share one.",
                )
            )
            get(atoms(got), f.site, 0) > 0 || throw(
                ArgumentError(
                    "SiteFamily \"$(f.name)\" claims \"$(symbol(sp))\", but that species " *
                        "carries no :$(f.site). A family owns the species carrying ITS " *
                        "site symbol; matching by name alone would give this family a " *
                        "member that belongs to another one, and leave its own " *
                        "conservation row out of the matrix.",
                )
            )
            aggregate_state(got) == AS_SURFACE || throw(
                ArgumentError(
                    "SiteFamily \"$(f.name)\": member \"$(symbol(sp))\" is " *
                        "$(aggregate_state(got)) in the species list. A family qualifies " *
                        "its own copies as AS_SURFACE, but the system keeps the species " *
                        "you passed, and `idx_surface` is built from those — so the two " *
                        "would disagree about what is on a surface. Declare it " *
                        "`aggregate_state = AS_SURFACE, class = SC_SURFCOMPLEX`.",
                )
            )

            owner = get(seen_members, i, nothing)
            owner === nothing || throw(
                ArgumentError(
                    "species \"$(symbol(got))\" is claimed by both SiteFamily " *
                        "\"$owner\" and \"$(f.name)\". One species belongs to one " *
                        "family: shared ownership means one of the two site budgets " *
                        "has no row of its own.",
                )
            )
            seen_members[i] = f.name
            push!(group, i)
        end
        push!(groups, group)
    end

    orphans = setdiff(idx_surface, keys(seen_members))
    isempty(orphans) || throw(
        ArgumentError(
            "species $(join([symbol(species[i]) for i in orphans], ", ")) are in " *
                "AS_SURFACE but belong to no declared site family.",
        )
    )

    return collect(SiteFamily, site_families), groups
end

"""
    _declared(ss) -> String

The name of the declaration a solid-solution phase is an instance of.

Falls back to the phase's own name for any `AbstractSolidSolutionPhase` that does
not carry the field, so a user-defined phase type keeps working.
"""
_declared(ss) = hasproperty(ss, :declared) ? getproperty(ss, :declared) : name(ss)

"""
    _instances(ss) -> Int

How many coexisting compositions a declaration asks for; 1 for a phase type that
does not carry the field.
"""
_instances(ss) = hasproperty(ss, :instances) ? Int(getproperty(ss, :instances)) : 1

"""
    _expand_instances(species, solid_solutions) -> (species, solid_solutions)

Give every declaration asking for `instances > 1` the extra species it needs.

A `SolidSolutionPhase` with `instances = k` becomes `k` phases: the declaration
itself, then `k-1` copies named `"\$name#2"`, `"\$name#3"`, … whose end-members
are copies of the originals under `"\$symbol#2"`, `"\$symbol#3"`, … Each copy
shares the whole property dictionary of the species it was made from, so the two
are one substance under two labels and carry byte-identical thermodynamic data.

Why the copies are needed at all: the composition vector has one entry per
species, so a species belongs to exactly one phase. Two coexisting compositions
of one substance — which is what the Gibbs minimum is inside a spinodal — cannot
be written down without the substance appearing twice.

Why it is done by duplicating the species rather than by letting two phases share
them: it keeps `ss_groups` **disjoint**, so the activity assembly, the mole
fraction fill and the optimality certificate need no change at all. Sharing an
index would make the last group written win.

The duplicated column leaves the row rank of the conservation matrix unchanged —
it is a copy of a column already there — so conservation is untouched, and the
copies start at zero amount, so a budget computed as `A n` is unchanged too.

Returns the inputs unchanged, and without copying, when no declaration asks for
more than one instance. That is every system the shipped data describes.
"""
function _expand_instances(species::AbstractVector, solid_solutions)
    solid_solutions === nothing && return species, solid_solutions
    phases = collect(solid_solutions)
    any(ss -> _instances(ss) > 1, phases) || return species, solid_solutions

    extra_species = eltype(species)[]
    expanded = AbstractSolidSolutionPhase[]
    known = Set(symbol(s) for s in species)
    # The copy is made from the species AS IT SITS IN THE LIST, not from the
    # phase's qualified end-member. The two differ in `class` -- an end-member is
    # requalified to `SC_SSENDMEMBER` when the phase is built, while the entry in
    # the species list keeps whatever the caller passed -- and copying the wrong
    # one would put the two instances in different class partitions of the same
    # system. `ChemicalSystem` then requalifies both, identically, downstream.
    in_list = Dict(symbol(s) => s for s in species)

    for ss in phases
        push!(expanded, ss)
        k = _instances(ss)
        k > 1 || continue
        for i in 2:k
            copies = map(end_members(ss)) do em
                src = get(in_list, symbol(em), nothing)
                src === nothing && error(
                    "SolidSolutionPhase \"$(name(ss))\": end-member " *
                        "\"$(symbol(em))\" is not in the species list, so its " *
                        "instance $i cannot be built. Add it to the species vector " *
                        "first.",
                )
                sym = "$(symbol(em))#$i"
                sym in known && error(
                    "SolidSolutionPhase \"$(name(ss))\" needs the symbol " *
                        "\"$sym\" for instance $i of its end-member " *
                        "\"$(symbol(em))\", and a species of that name is already " *
                        "in the system. Rename it, or declare fewer instances.",
                )
                push!(known, sym)
                cp = with_symbol(src, sym)
                push!(extra_species, cp)
                cp
            end
            # The convexity test has already run on the declaration, and it ran
            # the other way round: it REQUIRED a spinodal. Running it again here
            # would refuse the very phase this branch exists to build, so the
            # instance is constructed with the check off and `instances = 1` --
            # it is one composition; the declaration is what holds several.
            push!(
                expanded,
                SolidSolutionPhase(
                    "$(name(ss))#$i", copies; model = model(ss),
                    check_convexity = false, declared = _declared(ss),
                ),
            )
        end
    end
    return vcat(collect(species), extra_species), expanded
end

"""
    with_instances(cs, name => k, ...; T = 298.15) -> ChemicalSystem

`cs` with the solid solution `name` declared with `k` instances, every other
species and phase unchanged.

The system is rebuilt from its declarations: the copies of the end-members that
earlier instances needed are dropped, the phase is declared again with `k`
instances, and [`ChemicalSystem`](@ref) builds the copies that number asks for,
under the symbols `"\$symbol#2"`, … The primaries, and therefore the order of the
components of a budget, are those of `cs`, so a budget written for `cs` is a
budget for the result. [`with_instances(state, cs)`](@ref with_instances) maps a
state across.

`T` is the temperature the phase is checked at: more than one instance is
refused where its mixing energy is convex, as at declaration. A system with
kinetic species is refused, since a kinetic step works on a system it was given.
"""
function with_instances(cs::ChemicalSystem, changes::Pair{<:AbstractString, <:Integer}...; T::Real = T_STANDARD)
    isempty(cs.idx_kinetic) || throw(
        ArgumentError(
            "with_instances: this system has kinetic species, and a kinetic step " *
                "works on the system it was given; declare the instances before building it."
        ),
    )
    phases = cs.solid_solutions
    phases === nothing && throw(ArgumentError("with_instances: this system declares no solid solution."))
    declarations = [ss for ss in phases if name(ss) == _declared(ss)]
    wanted = Dict(String(first(c)) => Int(last(c)) for c in changes)
    for n in keys(wanted)
        any(ss -> name(ss) == n, declarations) || throw(
            ArgumentError(
                "with_instances: no solid solution named \"$n\"; this system declares " *
                    join(("\"" * name(ss) * "\"" for ss in declarations), ", ") * "."
            ),
        )
    end
    copies = Set(symbol(em) for ss in phases if name(ss) != _declared(ss) for em in end_members(ss))
    species = [sp for sp in cs.species if !(symbol(sp) in copies)]
    redeclared = map(declarations) do ss
        k = get(wanted, name(ss), nothing)
        k === nothing && return ss
        return SolidSolutionPhase(
            name(ss), end_members(ss); model = model(ss), instances = k, T = T,
            acknowledge_degenerate = true,
        )
    end
    return ChemicalSystem(species, cs.SM.primaries; solid_solutions = redeclared, site_families = cs.site_families)
end

"""
    _refuse_overlapping_solid_solutions(solid_solutions)

Refuse a system in which two declared solid solutions describe the same
substance, naming the pair.

The case this exists for is the calcium silicate hydrate. CEMDATA18
[Lothenbach2019](@cite) carries three descriptions of it — `CSHQ`, `CNASH_ss`
and the `ECSH` family — and they are three *models of one gel*, not three
phases. Declaring two of them counts the same hydrate twice: the calcium, the
silicon and the alkalis all enter the element balance once and come out
distributed over two phases that are supposed to be alternatives.

The overlap is exact rather than approximate, which is what makes it detectable
here: `KSiOH` (an end-member of `CSHQ`), `ECSH1-KSH` and `ECSH2-KSH` all carry
the formula `((KOH)2.5SiO2H2O)0.2`. Two end-members of two different declared
phases with the same composition are therefore the signature, and the check is
composition-based rather than name-based so that it does not depend on the
database's naming.

Two end-members of the SAME phase may of course share nothing — that is a
mixture — and a pure phase repeating a mixing phase's composition is a separate
question the rank test upstream already refuses.

**Composition cannot see every pair, so a second test reads the database's own
models.** `CSHQ` and `CNASH_ss` share no composition — no end-member of one is a
substance of the other — and until 0.25.1 the pair passed, the gel counted
twice without a word. `data/gel_models.toml` lists the end-member symbols of each
model of one gel, and two declared phases whose end-members belong to two
different models of the same gel are refused, naming both models. Matching is by
symbol, so it does not depend on the name a phase is declared under.

**Instances of one declaration are exempt**, and that exemption is the whole
reason `SolidSolutionPhase` carries a `declared` field. A miscibility gap is
represented by the same binary present twice, on purpose, as two coexisting
compositions — which is composition overlap of exactly the kind this function
refuses. Telling the deliberate case from the mistake cannot be done by
composition, since they look identical; it is done by provenance, and two phases
share a `declared` name only when `ChemicalSystem` itself made the second from
the first.
"""
function _refuse_overlapping_solid_solutions(solid_solutions)
    phases = collect(solid_solutions)
    length(phases) < 2 && return nothing
    for i in eachindex(phases), j in (i + 1):lastindex(phases)
        _declared(phases[i]) == _declared(phases[j]) && continue
        for a in end_members(phases[i]), b in end_members(phases[j])
            atoms(a) == atoms(b) || continue
            error(
                "solid solutions \"$(name(phases[i]))\" and " *
                    "\"$(name(phases[j]))\" both contain the composition " *
                    "$(unicode(a)) — as \"$(symbol(a))\" and \"$(symbol(b))\". Two " *
                    "declared phases sharing a composition describe the same " *
                    "substance twice, so its elements would be distributed over " *
                    "both. CEMDATA18's `CSHQ`, `CNASH_ss` and `ECSH` families are " *
                    "three models of one C-S-H gel: declare exactly one of them.",
            )
        end
    end
    _refuse_two_gel_models(phases)
    return nothing
end

# Below this Gibbs-energy difference, per formula unit of the smaller
# end-member, two proportional end-members are one substance: 0.1 RT, about
# 250 J/mol at 25 °C. The case it was set on differs by 5 J/mol; two datasets'
# values of one mineral that differ by tens of kJ/mol are not flagged, since one
# of them then always wins.
const _ONE_SUBSTANCE_RT = 0.1

"""
    _composition_ratio(a, b) -> Union{Float64, Nothing}

The factor `k` with `atoms(b) = k · atoms(a)`, or `nothing` when the two
compositions are not proportional.
"""
function _composition_ratio(a, b)
    A, B = atoms(a), atoms(b)
    keys(A) == keys(B) || return nothing
    isempty(A) && return nothing
    k = nothing
    for (el, na) in A
        r = Float64(B[el]) / Float64(na)
        k === nothing && (k = r)
        # CEMDATA18 writes rounded formulas (Al0.6666667 for 2/3), so "proportional"
        # has to allow for the rounding of the last printed digit.
        isapprox(r, k; rtol = 1.0e-6) || return nothing
    end
    return k
end

# The standard Gibbs energy of formation at 25 °C and 1 bar, in J/mol, or
# `nothing` for a species that does not carry one.
function _g298(s)
    return try
        ustrip(us"J/mol", s[:ΔₐG⁰](T = T_STANDARD_Q, P = P_STANDARD_Q; unit = true))
    catch
        nothing
    end
end

"""
    _warn_one_substance_two_phases(solid_solutions)

Warn when two declared solid solutions, one of them non-ideal, hold one
substance under two normalizations: an end-member of one whose composition is
`k` times that of an end-member of the other, with Gibbs energies that agree to
the same factor within `0.1 RT` per formula unit of the smaller, at 25 °C.

The case it was written for is CEMDATA18's aluminate sulfate. `ettringite03_ss`,
the SO4 end-member of the SO4/CO3 AFt binary, is ettringite divided by three,
5 J/mol apart, and GEM-Selektor's CEMDATA18 list declares the binary beside the
`ettringite` solid solution. With ideal mixing the two phases share the
substance and the certificate holds. With the published Redlich–Kister model
[RedlichKister1948](@cite) on the binary, measured on four cement pastes, the
certified search stopped short of the solution: moving the sulfate from one
phase to the other changes the Gibbs energy by next to nothing, a flat direction
the Newton cannot settle. Declaring the substance once certified all four.

A warning, not a refusal: with ideal mixing the double declaration is harmless,
and a page may keep it on purpose rather than decide in advance which
normalization the answer uses. Two end-members of exactly one composition are
refused earlier, by `_refuse_overlapping_solid_solutions`, and instances of one
declaration are exempt, as they are there.
"""
function _warn_one_substance_two_phases(solid_solutions)
    phases = collect(solid_solutions)
    length(phases) < 2 && return nothing
    RT = R_GAS * T_STANDARD
    for i in eachindex(phases), j in (i + 1):lastindex(phases)
        P, Q = phases[i], phases[j]
        _declared(P) == _declared(Q) && continue
        (model(P) isa IdealSolidSolutionModel && model(Q) isa IdealSolidSolutionModel) &&
            continue
        for a in end_members(P), b in end_members(Q)
            k = _composition_ratio(a, b)
            (k === nothing || isapprox(k, 1.0; rtol = 1.0e-6)) && continue
            ga, gb = _g298(a), _g298(b)
            (ga === nothing || gb === nothing) && continue
            gap = abs(gb - k * ga) / max(k, 1.0)
            gap < _ONE_SUBSTANCE_RT * RT || continue
            nonideal = model(P) isa IdealSolidSolutionModel ? name(Q) : name(P)
            @warn "solid solutions \"$(name(P))\" and \"$(name(Q))\" hold one " *
                "substance twice: \"$(symbol(b))\" is \"$(symbol(a))\" times " *
                "$(round(k; sigdigits = 4)), and their Gibbs energies agree to that " *
                "factor within $(round(gap; sigdigits = 2)) J/mol per formula unit. " *
                "With \"$(nonideal)\" non-ideal, moving the substance from one phase " *
                "to the other changes the Gibbs energy by next to nothing, and the " *
                "certified search can stop short of the solution. Declare it once: " *
                "keep the phase whose mixing the problem needs, and drop the other " *
                "description of the substance."
        end
    end
    return nothing
end

"""
    _warn_pure_phase_beside_nonideal(species, solid_solutions, idx_members)

Warn when a pure solid declared in the system repeats, up to a factor, an
end-member of a declared non-ideal solid solution, with Gibbs energies that agree
to that factor within `0.1 RT` per formula unit at 25 °C: `ettringite` declared
pure beside `AFt_SO4_CO3`, whose SO4 end-member `ettringite03_ss` is ettringite
divided by three. It is the case of [`_warn_one_substance_two_phases`](@ref
ChemistryLab._warn_one_substance_two_phases) with one of the two phases pure,
and the same flat direction: moving the substance between the pure phase and the
mixing phase changes the Gibbs energy by next to nothing.

`idx_members` are the indices, in `species`, of the end-members of the declared
solid solutions; every other crystalline species is a pure phase.
"""
function _warn_pure_phase_beside_nonideal(species, solid_solutions, idx_members)
    members = Set(idx_members)
    pure = [s for (i, s) in enumerate(species) if aggregate_state(s) == AS_CRYSTAL && !(i in members)]
    isempty(pure) && return nothing
    RT = R_GAS * T_STANDARD
    for P in solid_solutions
        model(P) isa IdealSolidSolutionModel && continue
        for a in end_members(P), b in pure
            symbol(a) == symbol(b) && continue
            k = _composition_ratio(a, b)
            k === nothing && continue
            ga, gb = _g298(a), _g298(b)
            (ga === nothing || gb === nothing) && continue
            gap = abs(gb - k * ga) / max(k, 1.0)
            gap < _ONE_SUBSTANCE_RT * RT || continue
            @warn "the pure phase \"$(symbol(b))\" and the solid solution \"$(name(P))\" " *
                "hold one substance twice: \"$(symbol(b))\" is \"$(symbol(a))\" times " *
                "$(round(k; sigdigits = 4)), and their Gibbs energies agree to that factor " *
                "within $(round(gap; sigdigits = 2)) J/mol per formula unit. With " *
                "\"$(name(P))\" non-ideal, moving the substance from one phase to the " *
                "other changes the Gibbs energy by next to nothing, and the certified " *
                "search can stop short of the solution. Declare it once: keep the phase " *
                "whose mixing the problem needs, and drop the other description."
        end
    end
    return nothing
end

const _GEL_MODELS_LOCK = ReentrantLock()
const _GEL_MODELS = Ref{Union{Nothing, Dict{String, NamedTuple{(:gel, :model), Tuple{String, String}}}}}(nothing)

"""
    _gel_models() -> Dict{String, NamedTuple{(:gel, :model)}}

The gel model each listed end-member symbol belongs to, read once from
`data/gel_models.toml`.
"""
function _gel_models()
    return lock(_GEL_MODELS_LOCK) do
        cached = _GEL_MODELS[]
        cached === nothing || return cached
        path = joinpath(pkgdir(@__MODULE__), "data", "gel_models.toml")
        include_dependency(path)
        table = Dict{String, NamedTuple{(:gel, :model), Tuple{String, String}}}()
        for entry in TOML.parsefile(path)["gel_model"], em in entry["end_members"]
            table[em] = (gel = entry["gel"], model = entry["model"])
        end
        _GEL_MODELS[] = table
        return table
    end
end

"""
    _refuse_two_gel_models(phases)

Refuse two declared solid solutions whose end-members belong to two different
models of one gel in `data/gel_models.toml` — `CSHQ` with `CNASH_ss`, which share
no composition and so pass the composition test of
`_refuse_overlapping_solid_solutions`. Instances of one declaration are exempt,
as they are there.
"""
function _refuse_two_gel_models(phases)
    models = _gel_models()
    # `#2`, `#3`: the symbols `ChemicalSystem` gives the copies of an instance.
    base(sym) = replace(String(sym), r"#\d+$" => "")
    tags = map(phases) do p
        unique(models[base(symbol(em))] for em in end_members(p) if haskey(models, base(symbol(em))))
    end
    for i in eachindex(phases), j in (i + 1):lastindex(phases)
        _declared(phases[i]) == _declared(phases[j]) && continue
        for a in tags[i], b in tags[j]
            a.gel == b.gel && a.model != b.model || continue
            error(
                "solid solutions \"$(name(phases[i]))\" and \"$(name(phases[j]))\" " *
                    "are two models of one $(a.gel) gel, `$(a.model)` and " *
                    "`$(b.model)` (data/gel_models.toml). Declared together they " *
                    "count the same hydrate twice: declare exactly one of them, and " *
                    "to compare the two, build one system with each.",
            )
        end
    end
    return nothing
end


"""
    _refuse_sites_on_mixing_hosts(site_families, solid_solutions)

Refuse a site family whose host is an end member of a declared solid solution
that carries an element the family binds, naming the family, the phase and the
elements.

The case is the C-S-H again. The published surface models of the C-S-H bind
calcium and alkalis on its silanol sites, and CSHQ holds its calcium and its
alkalis in its end members. Hosting the sites on CSHQ would let one mole of
calcium be counted twice, in the solid and on its surface; in the paste of
[Guo2018](@cite) the surface complexes hold about a fifth of the calcium of the
C-S-H. The elements a family binds are those of its complexes that its free
site does not carry, less hydrogen and oxygen, which a protonation or a
hydroxylation exchanges with the water. A family that binds nothing the phase
carries is accepted.

The way out is two stages: [`freeze_solid_solution`](@ref) sets the solid
solution aside after a first equilibrium, and the sites go on the frozen solid.

A family whose support names no host is held to the same rule unless the support
is declared `external = true`: its sites could otherwise sit on the solid
solution unnoticed, and the declaration is what says they belong to a solid
outside the system.
"""
function _refuse_sites_on_mixing_hosts(site_families, solid_solutions)
    (site_families === nothing || solid_solutions === nothing) && return nothing
    for f in site_families
        host = f.support.host
        f.support.external && continue
        free = keys(atoms(f.free_site))
        bound = setdiff(Set(e for c in f.complexes for e in keys(atoms(c))), free, (:H, :O))
        if host === nothing
            for ss in solid_solutions
                shared = sort!(collect(intersect(bound, Set(e for em in end_members(ss) for e in keys(atoms(em))))))
                isempty(shared) && continue
                throw(
                    ArgumentError(
                        "SiteFamily \"$(f.name)\" names no host and binds $(join(shared, ", ")), " *
                            "which the solid solution \"$(name(ss))\" also holds. If its sites " *
                            "are on that phase, those elements are counted twice, in the solid " *
                            "and on its surface. If they belong to a solid outside this system, " *
                            "such as a gel frozen by `freeze_solid_solution`, declare the support " *
                            "with `external = true`."
                    )
                )
            end
            continue
        end
        for ss in solid_solutions
            any(em -> symbol(em) == host, end_members(ss)) || continue
            shared = sort!(collect(intersect(bound, Set(e for em in end_members(ss) for e in keys(atoms(em))))))
            isempty(shared) && continue
            throw(
                ArgumentError(
                    "SiteFamily \"$(f.name)\" is hosted by \"$host\", an end member of the solid " *
                        "solution \"$(name(ss))\", and binds $(join(shared, ", ")), which that phase " *
                        "also holds: those elements would be counted twice, in the solid and on its " *
                        "surface. Freeze the solid solution after a first equilibrium " *
                        "(`freeze_solid_solution`) and put the sites on the frozen solid."
                )
            )
        end
    end
    return nothing
end


# The standard energies of a database are counted from a zero it chooses: the
# elements in their reference states for a database of formation properties
# (ThermoFun), its master species for a database of reactions (PHREEQC), each
# master at zero at every temperature. Within one database the choice does not
# change the equilibrium; beside the species of another it is meaningless, a
# master species of PHREEQC being at zero where its formation energy from the
# elements is hundreds of kilojoules. A species read from a database carries its
# zero under `:gauge`, and a system holding two is refused. A species built by
# hand carries none and is not checked.
function _refuse_mixed_gauges(species)
    seen = Dict{String, String}()
    for s in species
        haskey(s, :gauge) || continue
        g = s[:gauge]
        haskey(seen, g) || (seen[g] = symbol(s))
    end
    length(seen) > 1 || return nothing
    throw(
        ArgumentError(
            "the species come from databases whose energies are counted from different zeros: " *
                join(("$(v) (and others) from $(k)" for (k, v) in seen), "; ") * ". " *
                "Equilibria are meaningful only among the species of one of them; build the system " *
                "from one database, or give the others energies counted from the same zero.",
        )
    )
end

"""
    ChemicalSystem(species, primaries=species; kinetic_species, solid_solutions) -> ChemicalSystem

Construct a fully typed `ChemicalSystem` from a vector of species,
an optional vector of primary species, optional kinetic species with rates,
and optional solid-solution phases.

All derived fields are computed once at construction time and remain
consistent for the lifetime of the object.

# Arguments

  - `species`: vector of `AbstractSpecies`.
  - `primaries`: subset used as independent components (default: all species).
  - `kinetic_species`: `nothing` (default) or a dictionary / vector of pairs mapping
    each kinetic species (by name `String` or `Species` object) to its rate function.
    Rate functions must be callable as `(T, P, t, n, lna, n_initial) → Real [mol/s]`
    (see [`KineticFunc`](@ref)). The rate is given per mole of kinetic species
    (stoichiometric coefficient = 1); the constructor corrects by `1/|νₖ|` automatically.
    When provided, the nullspace N of the stoichiometric matrix is diagonalized so that
    each kinetic species appears in exactly one reaction. Those reactions are stored in
    the `reactions` field with their rate attached via `rxn[:rate]`.
  - `solid_solutions`: vector of [`SolidSolutionPhase`](@ref) (default: `nothing`).
    When provided, end-members must already appear in `species` (matched by symbol) and
    must carry `aggregate_state = AS_CRYSTAL` and `class = SC_SSENDMEMBER`.

# Examples
```jldoctest
julia> sp = [
           Species("H2O"; aggregate_state=AS_AQUEOUS, class=SC_AQSOLVENT),
           Species("Na+"; aggregate_state=AS_AQUEOUS, class=SC_AQSOLUTE),
       ];

julia> cs = ChemicalSystem(sp);

julia> length(cs)
2

julia> cs["H2O"] == sp[1]
true
```

```jldoctest
julia> em1 = Species("AFm1"; aggregate_state=AS_CRYSTAL, class=SC_SSENDMEMBER);

julia> em2 = Species("AFm2"; aggregate_state=AS_CRYSTAL, class=SC_SSENDMEMBER);

julia> ss = SolidSolutionPhase("AFm", [em1, em2]);

julia> cs = ChemicalSystem([em1, em2]; solid_solutions=[ss]);

julia> cs.ss_groups
1-element Vector{Vector{Int64}}:
 [1, 2]

julia> cs.idx_ssendmembers
2-element Vector{Int64}:
 1
 2
```
"""
function ChemicalSystem(
        species::AbstractVector{T},
        primaries::AbstractVector{<:AbstractSpecies} = species;
        kinetic_species = nothing,
        solid_solutions::Union{Nothing, AbstractVector} = nothing,
        site_families::Union{Nothing, AbstractVector} = nothing,
    ) where {T <: AbstractSpecies}
    solid_solutions = _normalize_solid_solutions(solid_solutions)
    species, solid_solutions = _expand_instances(species, solid_solutions)
    _refuse_mixed_gauges(species)
    idx(f) = findall(f, species)
    # Extract kinetic species keys for StoichMatrix construction
    kin_keys = if isnothing(kinetic_species)
        nothing
    elseif kinetic_species isa AbstractDict
        collect(keys(kinetic_species))
    else
        # Vector of Pairs
        [first(p) for p in kinetic_species]
    end
    CSM = CanonicalStoichMatrix(species)
    SM = StoichMatrix(species, primaries; kinetic_species = kin_keys)

    # Build kinetic reactions from diagonalized nullspace
    idx_kinetic = isnothing(kin_keys) ? Int[] : _resolve_kinetic_indices(kin_keys, SM.species)
    kin_reactions = if isempty(idx_kinetic)
        Reaction[]
    else
        all_rxns = reactions(SM)
        kin_pairs = if kinetic_species isa AbstractDict
            collect(pairs(kinetic_species))
        else
            kinetic_species
        end
        rxn_list = Reaction[]
        for (name_or_sp, rate_fn) in kin_pairs
            sp_idx = _resolve_kinetic_indices([name_or_sp], SM.species)[1]
            # Find the unique reaction where this kinetic species has a non-zero coefficient
            matching = filter(all_rxns) do r
                any(!iszero(ν) && s == SM.species[sp_idx] for (s, ν) in r)
            end
            isempty(matching) && throw(
                ArgumentError(
                    "No reaction found for kinetic species \"$(symbol(SM.species[sp_idx]))\"."
                )
            )
            rxn = first(matching)
            # Get the stoichiometric coefficient and correct the rate
            νk = sum(ν for (s, ν) in rxn if s == SM.species[sp_idx]; init = 0)
            abs_νk = abs(νk)
            corrected_rate = isone(abs_νk) ? rate_fn :
                (T, P, t, n, lna, n0) -> rate_fn(T, P, t, n, lna, n0) / abs_νk
            rxn[:rate] = corrected_rate
            push!(rxn_list, rxn)
        end
        rxn_list
    end
    R = isempty(kin_reactions) ? AbstractReaction : eltype(kin_reactions)

    idx_surface = idx(s -> aggregate_state(s) == AS_SURFACE)
    sf, site_groups = _resolve_site_families(site_families, species, idx_surface)

    if isnothing(solid_solutions)
        return ChemicalSystem{T, R, typeof(CSM), typeof(SM), Nothing, typeof(sf)}(
            collect(T, species),
            Dict{String, T}(symbol(s) => s for s in species),
            idx(s -> aggregate_state(s) == AS_AQUEOUS),
            idx(s -> aggregate_state(s) == AS_CRYSTAL),
            idx(s -> aggregate_state(s) == AS_GAS),
            idx_surface,
            idx(s -> class(s) == SC_AQSOLUTE),
            idx(s -> class(s) == SC_AQSOLVENT),
            idx(s -> class(s) == SC_COMPONENT),
            idx(s -> class(s) == SC_GASFLUID),
            collect(R, kin_reactions),
            Dict{String, R}(symbol(r) => r for r in kin_reactions),
            CSM,
            SM,
            nothing,
            Vector{Int}[],
            Int[],
            sf,
            site_groups,
            idx_kinetic,
        )
    else
        ss_groups = map(solid_solutions) do ss
            map(end_members(ss)) do em
                idx_em = findfirst(s -> symbol(s) == symbol(em), species)
                idx_em === nothing &&
                    error(
                    "SolidSolutionPhase \"$(name(ss))\": end-member \"$(symbol(em))\" " *
                        "not found in the species list. Add it to the species vector first.",
                )
                idx_em
            end
        end
        idx_ssendmembers = isempty(ss_groups) ? Int[] : vcat(ss_groups...)
        _refuse_overlapping_solid_solutions(solid_solutions)
        _warn_one_substance_two_phases(solid_solutions)
        _warn_pure_phase_beside_nonideal(species, solid_solutions, idx_ssendmembers)
        _refuse_sites_on_mixing_hosts(sf, solid_solutions)
        ss = collect(solid_solutions)

        return ChemicalSystem{T, R, typeof(CSM), typeof(SM), typeof(ss), typeof(sf)}(
            collect(T, species),
            Dict{String, T}(symbol(s) => s for s in species),
            idx(s -> aggregate_state(s) == AS_AQUEOUS),
            idx(s -> aggregate_state(s) == AS_CRYSTAL),
            idx(s -> aggregate_state(s) == AS_GAS),
            idx_surface,
            idx(s -> class(s) == SC_AQSOLUTE),
            idx(s -> class(s) == SC_AQSOLVENT),
            idx(s -> class(s) == SC_COMPONENT),
            idx(s -> class(s) == SC_GASFLUID),
            collect(R, kin_reactions),
            Dict{String, R}(symbol(r) => r for r in kin_reactions),
            CSM,
            SM,
            ss,
            ss_groups,
            idx_ssendmembers,
            sf,
            site_groups,
            idx_kinetic,
        )
    end
end

"""
    ChemicalSystem(species, primaries::AbstractVector{<:AbstractString}; kinetic_species, solid_solutions) -> ChemicalSystem

Convenience constructor that resolves primary species from their symbol strings.

The components must span the species: every species is written as a combination
of them, and that combination **is** the conservation law the equilibrium enforces
for it. A list that cannot express a species is refused by name — see
[`StoichMatrix`](@ref).

# Examples
```jldoctest
julia> sp = [
           Species("H2O";  aggregate_state=AS_AQUEOUS),
           Species("NaCl"; aggregate_state=AS_CRYSTAL),
       ];

julia> cs = ChemicalSystem(sp, ["H2O", "NaCl"]);

julia> symbol.(cs.SM.primaries)
2-element Vector{String}:
 "H2O"
 "NaCl"
```

Water alone would not do. Sodium and chlorine would have nowhere to be conserved,
and `ChemicalSystem(sp, ["H2O"])` raises an `ArgumentError` naming `NaCl` rather
than projecting it onto `H2O` — a projection that would balance arithmetically
and let the solver make salt out of water.
"""
function ChemicalSystem(
        species::AbstractVector{T},
        primaries::AbstractVector{<:AbstractString};
        kinetic_species = nothing,
        solid_solutions::Union{Nothing, AbstractVector} = nothing,
        site_families = nothing,
    ) where {T <: AbstractSpecies}
    # Resolve string symbols to species objects, preserving order
    primaries_species = species[symbol.(species) .∈ Ref(primaries)]
    return ChemicalSystem(
        species,
        primaries_species;
        kinetic_species = kinetic_species,
        solid_solutions = solid_solutions,
        site_families = site_families,
    )
end

# ── Solid solution accessor ───────────────────────────────────────────────────

"""
    solid_solutions(cs::ChemicalSystem) -> Nothing | Vector{<:AbstractSolidSolutionPhase}

Return the registered solid-solution phases, or `nothing` if none were declared.

# Examples
```jldoctest
julia> em1 = Species("Em1"; aggregate_state=AS_CRYSTAL, class=SC_SSENDMEMBER);

julia> em2 = Species("Em2"; aggregate_state=AS_CRYSTAL, class=SC_SSENDMEMBER);

julia> cs = ChemicalSystem(
           [em1, em2];
           solid_solutions=[SolidSolutionPhase("SS", [em1, em2])],
       );

julia> solid_solutions(cs) isa Vector
true

julia> length(solid_solutions(cs))
1
```
"""
solid_solutions(cs::ChemicalSystem) = cs.solid_solutions

"""
    site_families(cs::ChemicalSystem) -> Union{Nothing, Vector{<:SiteFamily}}

The surface site families declared on `cs`, or `nothing` when it has no surface.

Mirrors [`solid_solutions`](@ref), and for the same reason: a family is a named
group of species with its own mixing, and the system has to carry the
declaration because the species alone do not say which budget they share.
"""
site_families(cs::ChemicalSystem) = cs.site_families

"""
    kinetic_species(cs::ChemicalSystem) -> SubArray

Return a view of the kinetic species declared at construction time.
Empty when no kinetic species were declared.
"""
kinetic_species(cs::ChemicalSystem) = @view cs.species[cs.idx_kinetic]

# ── AbstractVector interface ──────────────────────────────────────────────────

"""
    Base.size(cs::ChemicalSystem) -> Tuple

Return the size of the underlying species vector.

# Examples
```jldoctest
julia> cs = ChemicalSystem([Species("H2O"; aggregate_state=AS_AQUEOUS)]);

julia> size(cs)
(1,)
```
"""
Base.size(cs::ChemicalSystem) = size(cs.species)

"""
    Base.getindex(cs::ChemicalSystem, i::Int) -> AbstractSpecies

Return the species at position `i`.

# Examples
```jldoctest
julia> cs = ChemicalSystem([Species("H2O"; aggregate_state=AS_AQUEOUS)]);

julia> cs[1] == Species("H2O"; aggregate_state=AS_AQUEOUS)
true
```
"""
Base.getindex(cs::ChemicalSystem, i::Int) = cs.species[i]

"""
    Base.getindex(cs::ChemicalSystem, i::AbstractString) -> AbstractSpecies

Return the species whose symbol matches `i`. Runs in O(1) via `dict_species`.

# Examples
```jldoctest
julia> cs = ChemicalSystem([Species("H2O"; aggregate_state=AS_AQUEOUS)]);

julia> cs["H2O"] == Species("H2O"; aggregate_state=AS_AQUEOUS)
true
```
"""
Base.getindex(cs::ChemicalSystem, i::AbstractString) = cs.dict_species[i]

# ── Reaction accessor ─────────────────────────────────────────────────────────

"""
    get_reaction(cs::ChemicalSystem, sym::AbstractString) -> AbstractReaction

Return the reaction identified by symbol `sym`. Runs in O(1) via `dict_reactions`.

# Examples
```julia
cs = ChemicalSystem(
    [Species("H2O"; aggregate_state=AS_AQUEOUS)];
);
get_reaction(cs, "some_rxn")  # returns the Reaction with that symbol
```
"""
get_reaction(cs::ChemicalSystem, sym::AbstractString) = cs.dict_reactions[sym]

# ── Merge ─────────────────────────────────────────────────────────────────────

"""
    Base.merge(cs1::ChemicalSystem, cs2::ChemicalSystem) -> ChemicalSystem

Construct a new `ChemicalSystem` from the union of two systems.

Species are unioned by symbol — duplicates from `cs2` are discarded — and the
solid solutions and site families by name. `CSM`, `SM` and the reactions are
built from scratch from the full species list. Primaries are taken as the union
of both systems' primaries, filtered to those actually present in the merged
species list. A solid solution declared with several instances is declared
again, and its copies are built anew. Kinetic species are not carried over: a
kinetic system is to be declared as such.

In case of a conflict of symbol or name, `cs1` takes priority over `cs2`.
The return type is inferred from the merged collections and may differ from
`typeof(cs1)` or `typeof(cs2)` if they contain different concrete types.

# Examples
```jldoctest
julia> cs1 = ChemicalSystem(
           [Species("H2O"; aggregate_state=AS_AQUEOUS, class=SC_AQSOLVENT),
            Species("H+";  aggregate_state=AS_AQUEOUS, class=SC_AQSOLUTE)],
       );

julia> cs2 = ChemicalSystem(
           [Species("OH-"; aggregate_state=AS_AQUEOUS, class=SC_AQSOLUTE),
            Species("H2O"; aggregate_state=AS_AQUEOUS, class=SC_AQSOLVENT)],
       );

julia> cs = merge(cs1, cs2);

julia> length(cs)
3
```
"""
function Base.merge(cs1::ChemicalSystem, cs2::ChemicalSystem)
    # The phases as declared, and the copies of end-members that the instances
    # of a phase declared several times added to the species: the declarations
    # are carried over and the copies built again, as `with_instances` does.
    declarations(cs) = cs.solid_solutions === nothing ? AbstractSolidSolutionPhase[] :
        [ss for ss in cs.solid_solutions if name(ss) == _declared(ss)]
    copies(cs) = cs.solid_solutions === nothing ? Set{String}() :
        Set(symbol(em) for ss in cs.solid_solutions if name(ss) != _declared(ss) for em in end_members(ss))
    families(cs) = cs.site_families === nothing ? SiteFamily[] : collect(SiteFamily, cs.site_families)

    # Species: cs1 first, then those of cs2 it does not hold (cs1 wins on conflict).
    existing_symbols = Set(symbol.(cs1.species))
    extra_species = filter(s -> symbol(s) ∉ existing_symbols, cs2.species)
    dropped = union(copies(cs1), copies(cs2))
    all_species = filter(s -> symbol(s) ∉ dropped, vcat(cs1.species, extra_species))

    # Union of primaries: cs1 first, then new ones from cs2, among the species kept.
    existing_primary_symbols = Set(symbol.(cs1.SM.primaries))
    extra_primaries = filter(p -> symbol(p) ∉ existing_primary_symbols, cs2.SM.primaries)
    all_species_symbols = Set(symbol.(all_species))
    all_primaries = filter(p -> symbol(p) ∈ all_species_symbols, vcat(cs1.SM.primaries, extra_primaries))

    # Solid solutions and site families by name, cs1 winning on conflict.
    ss1, sf1 = declarations(cs1), families(cs1)
    names1, fnames1 = Set(name.(ss1)), Set(f.name for f in sf1)
    ss = vcat(ss1, filter(p -> name(p) ∉ names1, declarations(cs2)))
    sf = vcat(sf1, filter(f -> f.name ∉ fnames1, families(cs2)))
    return ChemicalSystem(
        all_species, all_primaries;
        solid_solutions = isempty(ss) ? nothing : ss, site_families = isempty(sf) ? nothing : sf,
    )
end

"""
    Base.merge(css::ChemicalSystem...) -> ChemicalSystem

Construct a new `ChemicalSystem` from the union of an arbitrary number of systems,
processed left-to-right. Earlier systems take priority over later ones
in case of symbol conflicts.

# Examples
```jldoctest
julia> cs1 = ChemicalSystem([Species("H2O"; aggregate_state=AS_AQUEOUS, class=SC_AQSOLVENT)]);

julia> cs2 = ChemicalSystem([Species("H+"; aggregate_state=AS_AQUEOUS, class=SC_AQSOLUTE)]);

julia> cs3 = ChemicalSystem([Species("OH-"; aggregate_state=AS_AQUEOUS, class=SC_AQSOLUTE)]);

julia> cs = merge(cs1, cs2, cs3);

julia> length(cs)
3
```
"""
Base.merge(css::ChemicalSystem...) = reduce(merge, css)

# ── Views by aggregate state ──────────────────────────────────────────────────

"""
    aqueous(cs::ChemicalSystem) -> SubArray

Return a view of all aqueous species.

# Examples
```jldoctest
julia> cs = ChemicalSystem([
           Species("H2O";  aggregate_state=AS_AQUEOUS),
           Species("NaCl"; aggregate_state=AS_CRYSTAL),
       ]);

julia> length(aqueous(cs))
1

julia> aggregate_state(aqueous(cs)[1]) == AS_AQUEOUS
true
```
"""
aqueous(cs::ChemicalSystem) = @view cs.species[cs.idx_aqueous]

"""
    crystal(cs::ChemicalSystem) -> SubArray

Return a view of all crystalline species.

# Examples
```jldoctest
julia> cs = ChemicalSystem([
           Species("H2O";  aggregate_state=AS_AQUEOUS),
           Species("NaCl"; aggregate_state=AS_CRYSTAL),
       ]);

julia> aggregate_state(crystal(cs)[1]) == AS_CRYSTAL
true
```
"""
crystal(cs::ChemicalSystem) = @view cs.species[cs.idx_crystal]

"""
    surface(cs::ChemicalSystem) -> SubArray

A view of the species bound to a surface site, i.e. those in `AS_SURFACE`.

They are deliberately **not** in `crystal(cs)`: `idx_crystal` is read by the
solver and by the start repair as "a pure mineral phase", which a site occupancy
is not. They are nonetheless counted in the *solid* compartment of a
[`ChemicalState`](@ref), because that is where their matter is.
"""
surface(cs::ChemicalSystem) = @view cs.species[cs.idx_surface]

"""
    gas(cs::ChemicalSystem) -> SubArray

Return a view of all gas-phase species.

# Examples
```jldoctest
julia> cs = ChemicalSystem([Species("CO2"; aggregate_state=AS_GAS)]);

julia> aggregate_state(gas(cs)[1]) == AS_GAS
true
```
"""
gas(cs::ChemicalSystem) = @view cs.species[cs.idx_gas]

# ── Views by class ────────────────────────────────────────────────────────────

"""
    solutes(cs::ChemicalSystem) -> SubArray

Return a view of all aqueous solute species.

# Examples
```jldoctest
julia> cs = ChemicalSystem([Species("Na+"; aggregate_state=AS_AQUEOUS, class=SC_AQSOLUTE)]);

julia> class(solutes(cs)[1]) == SC_AQSOLUTE
true
```
"""
solutes(cs::ChemicalSystem) = @view cs.species[cs.idx_solutes]

"""
    solvent(cs::ChemicalSystem) -> AbstractSpecies

Return the unique solvent species directly (not a view),
since a chemical system contains at most one solvent.

# Examples
```jldoctest
julia> cs = ChemicalSystem([Species("H2O"; aggregate_state=AS_AQUEOUS, class=SC_AQSOLVENT)]);

julia> class(solvent(cs)) == SC_AQSOLVENT
true
```
"""
solvent(cs::ChemicalSystem) = cs.species[cs.idx_solvent][1]  # unique element, return directly

"""
    components(cs::ChemicalSystem) -> SubArray

Return a view of all component species.

# Examples
```jldoctest
julia> cs = ChemicalSystem([Species("SiO2"; aggregate_state=AS_CRYSTAL, class=SC_COMPONENT)]);

julia> class(components(cs)[1]) == SC_COMPONENT
true
```
"""
components(cs::ChemicalSystem) = @view cs.species[cs.idx_components]

"""
    gasfluid(cs::ChemicalSystem) -> SubArray

Return a view of all gas/fluid species.

# Examples
```jldoctest
julia> cs = ChemicalSystem([Species("CO2"; aggregate_state=AS_GAS, class=SC_GASFLUID)]);

julia> class(gasfluid(cs)[1]) == SC_GASFLUID
true
```
"""
gasfluid(cs::ChemicalSystem) = @view cs.species[cs.idx_gasfluid]
