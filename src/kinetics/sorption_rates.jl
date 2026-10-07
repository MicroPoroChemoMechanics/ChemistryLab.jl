# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

# ── sorption_rate: a slow surface ────────────────────────────────────────────

"""
    sorption_rate(k, cs, rxn) -> KineticFunc

The rate of a reaction that takes one surface site from one state to another,
written as an elementary step whose reverse is fixed by the equilibrium
constant of the reaction:

```math
r = k\\, n_s \\prod_{i} a_i^{\\nu_i}\\,(1 - \\Omega),
\\qquad
\\ln\\Omega = \\sum_j \\nu_j \\left(\\ln a_j + \\frac{\\Delta_a G^\\circ_j}{RT}\\right)
```

in mol/s. `n_s` is the amount of the site species among the reactants, the
product runs over the other reactants with their coefficients `ν_i`, and `Ω` is
the saturation ratio of the reaction as written, over every participant,
products counted positive. `k` is in s⁻¹ (the activities are dimensionless): a
number, or an [`AbstractFunc`](@ref) of the temperature such as
[`arrhenius_rate_constant`](@ref).

# Why this form

With site fractions for the activities of the two site states, `x_s = n_s/N`
and `x_c = n_c/N` over the family's `N` sites, the second term is
`(k/K) n_c Π a_p^{ν_p}`, the products other than the site: the law is the
difference of a forward and a backward mass-action rate whose ratio is the
equilibrium constant `K` of the database. Three properties follow, whatever `k`:

  - the rate vanishes where the equilibrium solver puts the surface, `Ω = 1`,
    because both read the same activities and the same standard energies;
  - it has the sign of the affinity `𝒜 = −RT ln Ω`, so `𝒜 r ≥ 0`;
  - two reactions closing a cycle on one family have rate constants whose
    ratios multiply to the product of their `K`, which holds by construction.

A non-ideal site family (an exchange convention, a surface potential) puts its
activity coefficients in `Ω`, hence in the backward term only. Thermodynamics
fixes the ratio of the two terms, not how a non-ideality is shared between
them; this is the choice made here.

# Requirements, each checked

  - Exactly one reactant and one product occupy a surface site, each with a
    coefficient of one, and both belong to the same site family.
  - The site reactant is the species the reaction controls, which
    [`KineticReaction`](@ref) takes as the first reactant when the reaction
    has no crystal: write it first.
  - Every participant is a species of `cs` and carries a standard Gibbs energy.

Used in a [`KineticsProblem`](@ref), the reaction makes the whole site family
kinetic, the free site included; each of its states then needs a reaction that
forms it.

# Example

```julia
# ≡XOH + Ca+2 → ≡XOHCa+2, the free site first
rxn = Reaction(OrderedDict(free => 1, ca => 1), OrderedDict(bound => 1);
               symbol = "Ca on X", equal_sign = '→')
kr = KineticReaction(cs, rxn, sorption_rate(1.0e-3, cs, rxn))
```

See also [`transition_state`](@ref) for minerals.
"""
function sorption_rate(k, cs::ChemicalSystem, rxn::AbstractReaction)
    k isa Union{Real, AbstractFunc} || throw(
        ArgumentError("sorption_rate: `k` must be a number or an AbstractFunc of the temperature.")
    )
    k isa Real && k < 0 && throw(ArgumentError("sorption_rate: `k` must be non-negative, got $k."))

    surface = Set(cs.idx_surface)
    index_of(sp) = begin
        i = findfirst(s -> s == sp, cs.species)
        i === nothing && throw(
            ArgumentError("sorption_rate: \"$(symbol(sp))\" is not a species of the system.")
        )
        haskey(cs.species[i], :ΔₐG⁰) || throw(
            ArgumentError(
                "sorption_rate: \"$(symbol(sp))\" carries no standard Gibbs energy, " *
                    "which the saturation ratio of the reaction needs."
            )
        )
        i
    end
    reactants = [(index_of(sp), Float64(ν)) for (sp, ν) in rxn.reactants]
    products = [(index_of(sp), Float64(ν)) for (sp, ν) in rxn.products]

    site_r = [(i, ν) for (i, ν) in reactants if i in surface]
    site_p = [(i, ν) for (i, ν) in products if i in surface]
    (length(site_r) == 1 && length(site_p) == 1) || throw(
        ArgumentError(
            "sorption_rate: the reaction must take one site state to another, " *
                "one surface species on each side; it has $(length(site_r)) among its " *
                "reactants and $(length(site_p)) among its products."
        )
    )
    (i_s, ν_s), (i_c, ν_c) = only(site_r), only(site_p)
    (ν_s == 1 && ν_c == 1) || throw(
        ArgumentError(
            "sorption_rate: the two site states must have a coefficient of one " *
                "(got $ν_s and $ν_c): one site moves from one state to the other."
        )
    )
    fam = findfirst(grp -> i_s in grp, cs.site_groups)
    (fam !== nothing && i_c in cs.site_groups[fam]) || throw(
        ArgumentError(
            "sorption_rate: \"$(symbol(cs.species[i_s]))\" and " *
                "\"$(symbol(cs.species[i_c]))\" must be members of one site family."
        )
    )
    _find_mineral_idx(cs, rxn) == i_s || throw(
        ArgumentError(
            "sorption_rate: the reaction would be controlled by " *
                "\"$(symbol(cs.species[_find_mineral_idx(cs, rxn)]))\", not by the site " *
                "\"$(symbol(cs.species[i_s]))\". Write the site species first among the reactants."
        )
    )

    key(i) = _rate_lookup_key(cs, cs.species[i])
    site_key = key(i_s)
    forward = [(key(i), ν) for (i, ν) in reactants if i != i_s]
    saturation = vcat(
        [(key(i), -ν, cs.species[i][:ΔₐG⁰]) for (i, ν) in reactants],
        [(key(i), ν, cs.species[i][:ΔₐG⁰]) for (i, ν) in products],
    )

    f = (T, P, _t, n, lna, _n0) -> begin
        n_site = n[site_key]
        ln_fwd = log(max(n_site, oneunit(n_site) * 1.0e-300)) +
            sum((ν * lna[s] for (s, ν) in forward); init = zero(T))
        ln_Ω = sum(
            ν * (lna[s] + g(; T = T, P = P, unit = false) / (R_GAS * T))
                for (s, ν, g) in saturation
        )
        kv = k isa Real ? k : k(; T = T)
        # k e^{ln_fwd} (1 − Ω), with 1 − Ω as −expm1(ln Ω): exact near equilibrium.
        return -kv * exp(ln_fwd) * expm1(ln_Ω)
    end

    refs = (T = 298.15u"K", P = 1.0e5u"Pa")
    return KineticFunc(f, refs, u"mol/s")
end
