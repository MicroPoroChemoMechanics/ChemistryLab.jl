# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

# ── The surface potential as an unknown of the solve ─────────────────────────
#
# A [`DiffuseLayer`](@ref) written as an activity coefficient is exact, and the
# solver cannot always reach it. The inner loop of the dual Newton recovers a
# mixing phase's composition from its own stationarity with `lnγ` read at the
# previous iterate — a fixed point, which contracts only while the activity
# moves less than the composition does. `electrostatic_stiffness` is that
# factor, and above roughly five it does not.
#
# Carrying `Ψ` as an unknown removes the fixed point instead of taming it. The
# activity model is then TOLD its potential, so during the inner solve `lnγ` is
# a constant and the loop converges in one pass, exactly as it does for ideal
# mixing; the coupling moves into the OUTER Newton, which has a Jacobian for it.
#
# The machinery is the one `_constraint_blocks` already provides — `nq` extra
# unknowns, `hq` for an activity model that depends on them, `cq` for the
# equations that close them — and the adiabatic constraint is the template:
# there the unknown is a temperature the activity model must be rebuilt at,
# here it is a potential it must be evaluated at.

"""
    needs_potential_unknown(model) -> Bool

Whether a site mixing model's electrostatic term must be carried as an unknown
of the solve rather than evaluated from the composition.

`true` for [`DiffuseLayer`](@ref) and `false` for everything else, including
[`ConstantCapacitance`](@ref) — whose potential is linear in the composition,
so the fixed point contracts and eliminating it costs nothing.
"""
needs_potential_unknown(::AbstractSiteMixingModel) = false
needs_potential_unknown(::DiffuseLayer) = true
needs_potential_unknown(m::ConstantCapacitance) = needs_potential_unknown(m.base)

"""
    _potential_families(cs) -> Vector{Int}

Which of `cs`'s site families need a potential unknown, by index.

`cs.site_families` is `nothing` — not an empty vector — for a system that
declares no surface, which is most of them, so the guard is the first line
rather than an afterthought.
"""
function _potential_families(cs::ChemicalSystem)
    fams = cs.site_families
    fams === nothing && return Int[]
    return [k for (k, f) in enumerate(fams) if needs_potential_unknown(f.model)]
end

"""
    _potential_supports(cs) -> Vector{Vector{Int}}

The distinct supports of `cs` that need a potential unknown, each as the list of
family indices on it.

One unknown per **surface**, not per family: ferrihydrite's strong and weak
sites are two families on one oxide and share one `Ψ`. See
[`support_group`](@ref).

Refuses a support whose families disagree — one asking for a potential and
another not, or two declaring different areas for the same surface. Both are
incoherent rather than merely unusual, and both would otherwise produce a
number.
"""
function _potential_supports(cs::ChemicalSystem)
    groups = support_group(cs)
    isempty(groups) && return Vector{Vector{Int}}()
    out = Vector{Vector{Int}}()
    seen = Set{Int}()
    for grp in groups
        first(grp) in seen && continue
        push!(seen, first(grp))
        wants = [needs_potential_unknown(cs.site_families[g].model) for g in grp]
        any(wants) || continue
        name = cs.site_families[first(grp)].support.name
        all(wants) || throw(
            ArgumentError(
                "the families on support \"$name\" disagree about the surface " *
                    "potential: $(count(wants)) of $(length(grp)) declare a model " *
                    "that needs one. A surface carries one potential, so either " *
                    "all of its families describe it or none does."
            ),
        )
        areas = unique(cs.site_families[g].model.area for g in grp)
        listed = join(areas, ", ")
        length(areas) == 1 || throw(
            ArgumentError(
                "the families on support \"$name\" declare different areas for it " *
                    "($listed m²). One surface has one area, and the charge " *
                    "density it carries depends on it."
            ),
        )
        push!(out, grp)
    end
    return out
end

"""
    _surface_potential_blocks(des, state, p, n0) -> NamedTuple or nothing

The parameter block the **system** contributes, as opposed to the one its
constraint does: one unknown `ψ̃ = FΨ/RT` per surface that needs it, with the
Gouy-Chapman closure as its equation.

`nothing` when no surface needs one, which is every system without a diffuse
layer — so nothing pays for the possibility.

# The closure, and why its residual is already dimensionless

```math
c_k(n, \\tilde\\psi) \\;=\\; \\tilde\\psi_k
  \\;-\\; 2\\operatorname{asinh}\\!\\left(\\frac{\\sigma_k(n)}{\\kappa\\sqrt{I(n)}}\\right)
```

`σ_k` is summed over every member of every family on the k-th surface. `ψ̃` is
in units of `RT` by construction, which is what the stationarity rows are in, so
no scaling is needed — unlike the enthalpy residual of an adiabatic solve, which
is `10⁵ J` and has to be divided by `RT` before it can sit in the same Newton
system.

`q0 = 0` starts the solve on an uncharged surface, which is the composition a
solve without electrostatics would return and therefore the natural cold start.
"""
function _surface_potential_blocks(des, state, p, n0)
    cs = state.system
    supports = _potential_supports(cs)
    isempty(supports) && return nothing

    models = [cs.site_families[first(grp)].model for grp in supports]
    members = [
        reduce(vcat, (cs.site_groups[g] for g in grp); init = Int[]) for grp in supports
    ]
    charges = [
        Float64[charge(sp) for g in grp for sp in site_members(cs.site_families[g])]
            for grp in supports
    ]
    solvent = isempty(cs.idx_solvent) ? 0 : only(cs.idx_solvent)
    ions = [i for i in cs.idx_solutes if !iszero(charge(cs.species[i]))]
    ion_z = Float64[charge(cs.species[i]) for i in ions]
    M_w = iszero(solvent) ? 1.0 : ustrip(us"kg/mol", cs.species[solvent][:M])
    nq = length(supports)

    # The activity model, evaluated at the potential the solve currently holds.
    hq = (x, q, params) ->
    des.lna(x, merge(params, (ψ_site = _scatter(cs, supports, q),)))

    cq = function (x, q, params)
        T = hasproperty(params, :T) ? params.T : 298.15
        I = _aqueous_ionic_strength(x, ions, ion_z, solvent, M_w)
        return [
            q[k] - diffuse_layer_potential(
                models[k], charges[k], [x[i] for i in members[k]], I, T,
            ) for k in 1:nq
        ]
    end

    return (;
        nq = nq, gq = nothing, hq = hq, cq = cq,
        Aq = zeros(Float64, size(des.A, 1), nq),
        q0 = zeros(Float64, nq), qscale = ones(Float64, nq),
        apply = (T, P, q) -> (T, P),
        supports = supports,
    )
end

"""
    _scatter(cs, supports, q) -> Vector

Each surface's potential placed at every family standing on it, `nothing`
elsewhere, so the activity kernel can index by family without knowing which
families share a surface.
"""
function _scatter(cs::ChemicalSystem, supports::Vector{Vector{Int}}, q)
    out = Vector{Any}(nothing, length(cs.site_families))
    @inbounds for (k, grp) in enumerate(supports), f in grp
        out[f] = q[k]
    end
    return out
end

"""
    _compose_blocks(a, b) -> NamedTuple

One parameter block from two: the constraint's unknowns first, the system's
after. Each side keeps its own indices into `q`, so neither has to know the
other exists.

# The one combination that is refused

Both sides may declare `hq`, the callback that rebuilds the activity model at
the current unknowns — an adiabatic solve does, because the Debye-Hückel
coefficients are functions of temperature, and a diffuse layer does, because it
is told its potential. Composing two of them would mean threading one's
parameter merge through the other's, and the shipped `hq` closures are opaque
to that.

Rather than compose them wrongly, this refuses and says so. A diffuse layer
under an adiabatic constraint is not a case anyone has yet; when someone does,
the fix is to have each block contribute the `params` it merges instead of a
whole `hq`, and this error message is where they will start.
"""
function _compose_blocks(a, b)
    b === nothing && return a
    na = a.nq
    if a.hq !== nothing && b.hq !== nothing
        throw(
            ArgumentError(
                "this constraint makes the activity model depend on an unknown of " *
                    "its own, and so does the surface potential of a `DiffuseLayer`. " *
                    "Composing the two is not implemented: each would have to merge " *
                    "its parameters through the other's. Use a `ConstantCapacitance`, " *
                    "or drop the constraint, or open an issue with the case."
            )
        )
    end
    return (;
        nq = na + b.nq,
        gq = a.gq === nothing ? nothing : (q, params) -> a.gq(q[1:na], params),
        hq = a.hq !== nothing ? (x, q, params) -> a.hq(x, q[1:na], params) :
            b.hq === nothing ? nothing :
            (x, q, params) -> b.hq(x, q[(na + 1):end], params),
        cq = (x, q, params) -> vcat(
            a.cq === nothing ? Float64[] : a.cq(x, q[1:na], params),
            b.cq(x, q[(na + 1):end], params),
        ),
        Aq = hcat(a.Aq, b.Aq),
        q0 = vcat(a.q0, b.q0),
        qscale = vcat(a.qscale, b.qscale),
        apply = (T, P, q) -> a.apply(T, P, q[1:na]),
    )
end
