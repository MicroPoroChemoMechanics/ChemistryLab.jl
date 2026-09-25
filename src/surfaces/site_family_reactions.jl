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
\\Delta_a G^\\circ_{\\text{complex}} = -RT \\ln 10 \\, \\log K
  + \\sum_{\\text{reactants}} \\nu_i \\Delta_a G^\\circ_i
  - \\sum_{\\text{other products}} \\nu_i \\Delta_a G^\\circ_i .
```

The free site at zero is a gauge while the site budget is fixed, and only then:
under `SITES_FOLLOW_HOST` the reference of the free site reaches the host's
solubility, and [`host_coupling_bias`](@ref) says what it must be.

The energies are constants, evaluated at `T` and `P`: a published log K carries
no temperature dependence, and neither does the family built from it.

# Examples

```julia
dat = read_sorption_model(joinpath(pkgdir(ChemistryLab), "test", "reference", "phreeqc.dat"))
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
        T = 298.15u"K", P = 1.0e5u"Pa",
    )
    is_site_symbol(Symbol(site)) || throw(
        ArgumentError("site must be one of SITE_SYMBOLS, such as \"Xw\"; got \"$site\".")
    )
    lookup = aqueous isa AbstractDict ? aqueous : Dict(symbol(s) => s for s in aqueous)
    Tq = T isa Real ? T * u"K" : T
    Pq = P isa Real ? P * u"Pa" : P
    RT = R_GAS * ustrip(us"K", Tq)
    energy(s) = ustrip(us"J/mol", s[:ΔₐG⁰](T = Tq, P = Pq; unit = true))
    aqueous_energy(sym) = energy(_reaction_species(lookup, sym, name))
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
        G = -RT * log(10) * logK
        for (sp, ν) in stoich
            surface(sp) && continue
            G -= ν * aqueous_energy(sp)      # reactants have ν < 0, so they add
        end
        c = Species(rename(only(product)); aggregate_state = AS_SURFACE, class = SC_SURFCOMPLEX)
        c[:ΔₐG⁰] = SymbolicFunc(G * u"J/mol")
        push!(complexes, c)
    end
    free === nothing && throw(ArgumentError("site_family \"$name\": no reaction given."))
    f = Species(rename(free); aggregate_state = AS_SURFACE, class = SC_SURFCOMPLEX)
    f[:ΔₐG⁰] = SymbolicFunc(0.0u"J/mol")
    return SiteFamily(name, f, complexes; capacity, support, model)
end

_equation_and_logk(r::SorptionReaction) = (r.equation, Float64(value(r.log_K)))
_equation_and_logk(r::Pair) = (String(first(r)), Float64(last(r)))

function _reaction_species(lookup, sym, name)
    haskey(lookup, sym) && return lookup[sym]
    sym == "H2O" && haskey(lookup, "H2O@") && return lookup["H2O@"]
    throw(
        ArgumentError(
            "site_family \"$name\": no aqueous species \"$sym\" among the ones given."
        )
    )
end
