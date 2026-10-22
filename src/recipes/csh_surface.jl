# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2026 Jean-François Barthélémy (Cerema, UMR MCD)

# The surface of the C-S-H of a paste, after Guo et al. (2018): the route of
# docs/src/examples/chloride_binding_blended.md (route A) applied to the
# certified equilibrium of a recipe. The gel is frozen at its composition, the
# sites of the published model are carried on it, and the paste is solved again
# with them; a Donnan layer may be added on the diffuse layer.

# The Debye length [m] of a solution of ionic strength `I` [mol/kg] at `T` [K],
# from the permittivity of water and the constants of the library.
function _debye_length(I, T)
    ε = water_relative_permittivity(T, P_STANDARD) * VACUUM_PERMITTIVITY
    return sqrt(ε * R_GAS * T / (2 * FARADAY^2 * I * 1000))
end

# The aqueous species a surface reaction written as PHREEQC writes it reads,
# every term of the equation that is not a surface species.
_reaction_solutes(eq::AbstractString, master::AbstractString) =
    [String(t) for side in split(eq, "=") for t in split(side, " + ") if !startswith(strip(t), master)] .|> strip

# The thickness [m] of the Donnan layer `with_csh_surface` is asked for: a layer,
# a length, a number of meters, or `:debye` for one Debye length `debye`.
_donnan_thickness(d::DonnanLayer, debye) = d.thickness
_donnan_thickness(d::Real, debye) = float(d)
_donnan_thickness(d::UnionAbstractQuantity, debye) = ustrip(us"m", d)
function _donnan_thickness(d::Symbol, debye)
    d === :debye || throw(ArgumentError("with_csh_surface: `donnan = :$d`; give a thickness, a `DonnanLayer` or `:debye`."))
    return debye
end

"""
    with_csh_surface(rs::RecipeState; donnan = nothing, water = :taken)
        -> NamedTuple

The paste of `rs` with the surface of its C-S-H, as the model of
[Guo2018](@citet) writes it: its site density and its specific surface, read
from its data file and referred to the silicon of the gel (the molar mass per
silicon of their C-S-H), its surface complexes with the ions the system holds
(those whose ion it lacks are left out), and a diffuse layer
([`DiffuseLayer`](@ref)). The route is that of
[Chloride binding in CEM III/A and CEM III/B](@ref sec-example-chloride-blended), route A:

 1. the solid solution `CSHQ` of `rs` is frozen at its composition
    ([`freeze_solid_solution`](@ref)), its alkalis returned to the solution,
    portlandite fixing the activities it shares with the gel;
 2. the sites are put on a support outside the system, the frozen gel;
 3. the paste is solved again with them ([`equilibrate_certified`](@ref)).

`donnan`, a [`DonnanLayer`](@ref) or a thickness (m, or `:debye` for one Debye
length of the solution), adds the layer of water on the diffuse layer
([`equilibrate_donnan`](@ref), its water `water = :taken` from the pore solution
of a closed paste, or `:added`). It is refused when the pore water is a film
thinner than the layer or than the Debye length, the case of a hardened paste
at the area of [Guo2018](@citet): a Donnan layer assumes a free solution beyond
it, and in such a film the diffuse layers of the walls overlap (a layer made
thin enough to fit has to hold the charge of the surface in so little water
that it has no fixed point). It is meant for a suspension.

Two numbers belong to the model rather than to the paste, and are worth reading
before the rest: the published model counts calcium on the surface on top of a
gel whose measured Ca/Si includes it, so that portlandite gives some up; and
the area makes the pore water a film of a nanometer or less on the gel, which
`film` reports, against the Debye length.

Refused, with the reason: a paste without `CSHQ` in its system, without C-S-H,
or without portlandite.

Returns `(; state, certificate, family, area, sites, film, debye, potential,
complexes, bound, frozen, donnan)`: the state of the second stage and its
certificate, the site family, the area (m²) and the sites (mol), the thickness
of the pore water on the area and the Debye length (m), the potential of the
surface (V), the moles of each surface species, the moles of each element held
by the complexes, what the freezing set aside, and, with `donnan`, the result
of [`equilibrate_donnan`](@ref) and the contents of the layer
([`diffuse_layer_contents`](@ref)).
"""
function with_csh_surface(rs::RecipeState; donnan = nothing, water::Symbol = :taken)
    eq1 = rs.state
    cs1 = eq1.system
    model = rs.model
    sss = cs1.solid_solutions === nothing ? [] : collect(cs1.solid_solutions)
    k = findfirst(ss -> ss.name == "CSHQ", sss)
    k === nothing && throw(
        ArgumentError(
            "with_csh_surface: the system has no solid solution CSHQ, the C-S-H the surface model is " *
                "carried on; choose a phase list with it."
        )
    )
    gel = solid_solution_totals(eq1, "CSHQ")
    Si = get(gel.elements, :Si, 0.0)
    # below a nanomole of silicon the gel is a trace the solver floors, not a gel
    Si > 1.0e-9 || throw(ArgumentError("with_csh_surface: the paste holds no C-S-H at this point."))

    # The sites and the area of Guo et al., referred to the silicon of the gel.
    guo(q) = literature_value("Guo2018", q)
    t = literature_table("Guo2018", "csh_formula")
    si = t.coefficient[findfirst(==("SiO2"), t.oxide)]
    M_per_si = sum(c / si * ustrip(us"g/mol", Species(ox)[:M]) for (ox, c) in zip(t.oxide, t.coefficient))
    n_sites = ustrip(us"mol/g", guo("csh_site_density")) * M_per_si * Si
    area = ustrip(us"m^2/g", guo("csh_specific_surface_area")) * M_per_si * Si

    # The complexes whose ions the system holds.
    members = Set(symbol(m) for m in sss[k].end_members)
    sp2 = [s for s in cs1.species if !(symbol(s) in members)]
    aqueous = [s for s in sp2 if aggregate_state(s) == AS_AQUEOUS]
    have = Set(symbol.(aqueous))
    rows = literature_table("Guo2018", "surface_reactions_phreeqc")
    reactions = [r => lk for (r, lk) in zip(rows.reaction, rows.log_K) if all(in(have), _reaction_solutes(r, "Csh_w"))]
    family = site_family(
        "Csh_w", reactions, aqueous;
        master = "Csh_w", site = "Xw", capacity = TotalSiteAmount(n_sites * u"mol"),
        support = SurfaceSupport("C-S-H", FixedSurfaceArea(area); external = true), model = DiffuseLayer(; area),
    )
    sites = vcat([family.free_site], family.complexes)
    others = [ss for (i, ss) in enumerate(sss) if i != k]
    cs2 = ChemicalSystem(
        vcat(sp2, sites), vcat(collect(cs1.SM.primaries), [family.free_site]);
        solid_solutions = isempty(others) ? nothing : others, site_families = [family],
    )
    first_stage = freeze_solid_solution(eq1, "CSHQ", cs2; buffer = "Portlandite")
    n = ustrip.(us"mol", first_stage.state.n)
    n[findfirst(s -> symbol(s) == symbol(family.free_site), cs2.species)] = n_sites
    start = ChemicalState(cs2, n .* u"mol"; T = temperature(eq1), P = pressure(eq1))
    eq2, cert = equilibrate_certified(start; model)

    # What the surface holds, and what it does to the pore water.
    nn = ustrip.(us"mol", eq2.n)
    idx = Dict(symbol(s) => i for (i, s) in enumerate(cs2.species))
    complexes = [symbol(s) => nn[idx[symbol(s)]] for s in sites]
    # the elements the complexes take from the solution: not the site's own,
    # nor oxygen and hydrogen, and not a trace (an element the paste lacks)
    own = Set(keys(atoms(family.free_site)))
    bound = Dict{Symbol, Float64}()
    for s in family.complexes, (el, c) in atoms(s)
        (el in own || el in (:O, :H)) && continue
        bound[el] = get(bound, el, 0.0) + c * nn[idx[symbol(s)]]
    end
    filter!(((el, x),) -> x > 1.0e-12 * n_sites, bound)
    T = ustrip(us"K", temperature(eq2))
    I = ionic_strength(eq2)
    w = nn[idx["H2O@"]] * ustrip(us"kg/mol", cs2.species[idx["H2O@"]][:M])
    film = w / 1000 / area
    debye = _debye_length(I, T)
    z = Float64[charge(s) for s in sites]
    psi = diffuse_layer_potential(family.model, z, [nn[idx[symbol(s)]] for s in sites], I, T) * R_GAS * T / FARADAY

    layer = nothing
    if donnan !== nothing
        thickness = _donnan_thickness(donnan, debye)
        film > max(thickness, debye) || throw(
            ArgumentError(
                "with_csh_surface: the pore water is a film of $(round(1.0e9film; sigdigits = 2)) nm on the " *
                    "$(round(area; sigdigits = 3)) m² of the C-S-H, against a Donnan layer of " *
                    "$(round(1.0e9thickness; sigdigits = 2)) nm and a Debye length of $(round(1.0e9debye; sigdigits = 2)) nm: " *
                    "a Donnan layer assumes a free solution beyond it, which the film does not hold."
            )
        )
        d = DonnanLayer(; thickness)
        r = equilibrate_donnan(eq2, d; model, water)
        layer = (; result = r, contents = diffuse_layer_contents(r.state, d))
    end
    return (;
        state = eq2, certificate = cert, family, area, sites = n_sites, film, debye, potential = psi,
        complexes, bound, frozen = first_stage.frozen, donnan = layer,
    )
end
