# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

# ── A site family from the reactions that define it ──────────────────────────
#
# A published surface model is a list of reactions and their constants, written
# the way PHREEQC writes them: `Hfo_wOH + Zn+2 = Hfo_wOZn+ + H+`, log K = −1.99.
# A `SiteFamily` wants species with standard energies. The bridge between the
# two is one line of arithmetic per reaction, and every page and test that built
# a surface by hand carried its own copy of it, with the one trap it has: the
# complex must carry the energies of the aqueous species its reaction consumes,
# or it certifies with none of it formed.

"""
    site_family(name, reactions, aqueous; master, site, capacity, support,
                model = IdealSiteMixing(), T = 298.15u"K", P = 1.0e5u"Pa")
        -> SiteFamily

A site family built from its surface reactions, written as PHREEQC writes them,
and their constants.

`reactions` holds `equation => log_K` pairs, or the [`SorptionReaction`](@ref)s
[`read_sorption_model`](@ref) returns. Each reaction forms **one** surface
complex from **the** free site and aqueous species, with a coefficient of one on
both. `master` is the prefix naming the surface species in the equations
(`"Hfo_w"`), and it is replaced by `site`, one of the [`SITE_SYMBOLS`](@ref)
(`"Xw"`), so that `Hfo_wOZn+` becomes the species `XwOZn+`.

`aqueous` supplies the aqueous species the reactions name, as a vector of species
or a dictionary by symbol. `H2O` in an equation is the solvent, looked up as
`H2O@` when the collection has no `H2O`.

The standard energies are built so that every reaction has its constant, with
the free site at zero and the aqueous species at their own `ΔₐG⁰`:

```math
\\Delta_a G^\\circ_{\\text{complex}}(T, P) = -RT \\ln 10 \\, \\log K
  + \\sum_{\\text{reactants}} \\nu_i \\Delta_a G^\\circ_i(T, P)
  - \\sum_{\\text{other products}} \\nu_i \\Delta_a G^\\circ_i(T, P) .
```

The free site at zero is a gauge while the site budget is fixed, and it is the
required value under `SITES_FOLLOW_HOST`, where the free sites are counted as part
of the host; [`host_coupling_bias`](@ref) explains why.

A published surface constant carries no reaction enthalpy, and `log K` is then
held at every temperature, as PHREEQC holds a constant given without one: the
energy of each complex follows the aqueous species of its reaction (ASSUMED:
`ΔᵣH = 0`). The keywords `T` and `P` are accepted for compatibility and do not
change the answer.

# Examples

```julia
dat = read_sorption_model(datapath("phreeqc.dat"))
weak = [r for r in dat.surfaces["Hfo_w"].reactions if !haskey(r.stoichiometry, "Hfo_sOH")]
family = site_family("Hfo_w", weak, aqueous; master = "Hfo_w", site = "Xw",
                     capacity = TotalSiteAmount(2.0e-4u"mol"), support)
```

See also: [`SiteFamily`](@ref), [`read_sorption_model`](@ref).
"""
function site_family(
        name::AbstractString, reactions, aqueous;
        master::AbstractString, site::AbstractString,
        capacity::AbstractSiteCapacity, support::SurfaceSupport,
        model::AbstractSiteMixingModel = IdealSiteMixing(),
        T = T_STANDARD_Q, P = P_STANDARD_Q,
    )
    is_site_symbol(Symbol(site)) || throw(
        ArgumentError("site must be one of SITE_SYMBOLS, such as \"Xw\"; got \"$site\".")
    )
    lookup = aqueous isa AbstractDict ? aqueous : Dict(symbol(s) => s for s in aqueous)
    Tq = T isa Real ? T * u"K" : T
    Pq = P isa Real ? P * u"Pa" : P
    surface(sym) = startswith(sym, master)
    rename(sym) = site * sym[(lastindex(master) + 1):end]

    free = nothing
    complexes = AbstractSpecies[]
    for r in reactions
        equation, logK = _equation_and_logk(r)
        stoich = _parse_sorption_stoichiometry(equation)
        surf = [(sp, ν) for (sp, ν) in stoich if surface(sp)]
        reactant = [sp for (sp, ν) in surf if ν < 0]
        product = [sp for (sp, ν) in surf if ν > 0]
        (
            length(reactant) == 1 && length(product) == 1 &&
                stoich[only(reactant)] == -1 && stoich[only(product)] == 1
        ) || throw(
            ArgumentError(
                "site_family \"$name\": \"$equation\" must turn one free site " *
                    "into one complex, each with a coefficient of one; species " *
                    "starting with \"$master\" are read as surface species."
            ),
        )
        free === nothing && (free = only(reactant))
        only(reactant) == free || throw(
            ArgumentError(
                "site_family \"$name\": \"$equation\" starts from " *
                    "\"$(only(reactant))\" where the others start from \"$free\"; " *
                    "one family has one free site."
            ),
        )
        # −RT ln10 log K minus the energies of the other participants at the
        # same (T, P): reactants have ν < 0, so they add.
        terms = Tuple{Float64, Any}[
            (Float64(ν), _reaction_species(lookup, sp, name)[:ΔₐG⁰]) for (sp, ν) in stoich if !surface(sp)
        ]
        lnK = log(10) * logK
        G = (T, P) -> -R_GAS * T * lnK -
            sum((ν * g(; T = T, P = P, unit = false) for (ν, g) in terms); init = zero(T * lnK))
        c = Species(rename(only(product)); aggregate_state = AS_SURFACE, class = SC_SURFCOMPLEX)
        c[:ΔₐG⁰] = NumericFunc(G, (:T, :P), (T = Tq, P = Pq), u"J/mol")
        push!(complexes, c)
    end
    free === nothing && throw(ArgumentError("site_family \"$name\": no reaction given."))
    f = Species(rename(free); aggregate_state = AS_SURFACE, class = SC_SURFCOMPLEX)
    f[:ΔₐG⁰] = SymbolicFunc(0.0u"J/mol")
    return SiteFamily(name, f, complexes; capacity, support, model)
end

# A constant in the number type it is given: a log K being fitted carries its
# derivative into the energy of its complex.
_equation_and_logk(r::SorptionReaction) = (r.equation, float(value(r.log_K)))
_equation_and_logk(r::Pair) = (String(first(r)), float(last(r)))

function _reaction_species(lookup, sym, name)
    haskey(lookup, sym) && return lookup[sym]
    sym == "H2O" && haskey(lookup, "H2O@") && return lookup["H2O@"]
    throw(
        ArgumentError(
            "site_family \"$name\": no aqueous species \"$sym\" among the ones given."
        )
    )
end
