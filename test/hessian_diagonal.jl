# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

using ChemistryLab
using DynamicQuantities
using ForwardDiff
using LinearAlgebra
using Random
using Test

include("reference_species.jl")

# ── The diagonal of the Gibbs Hessian, against the Jacobian it replaces ───────
#
# `_lna_diagonal!` gives `∂ln aᵢ/∂nᵢ` without forming the Jacobian of the
# activities; the reference is that Jacobian, by forward mode, of which only the
# diagonal is read. One system holds every kind of term the diagonal is built
# from: ions and a neutral species, a gas phase, a solid solution with an excess
# term, and four site families, one per electrostatic model and one without.
# The coefficients of the Truesdell–Jones and LLNL models below are arbitrary:
# what is tested is an identity, which holds whatever they are.

const RT25_HD = ChemistryLab.R_GAS * 298.15
_g0_hd(v) = SymbolicFunc(v * u"J/mol")

function _hd_surface(name, g; planes = nothing)
    s = Species(name; aggregate_state = AS_SURFACE, class = SC_SURFCOMPLEX)
    s[:ΔₐG⁰] = _g0_hd(g)
    return planes === nothing ? s : with_plane_charges(s, planes...)
end

# A family `X` of three members, free, protonated and deprotonated, on a support
# of its own.
function _hd_family(X, model; planes = false)
    pl(c) = planes ? (c, 0.0, 0.0) : nothing
    free = _hd_surface("$(X)OH", 0.0; planes = pl(0.0))
    prot = _hd_surface("$(X)OH2+", -RT25_HD * log(10.0^7.3); planes = pl(1.0))
    depr = _hd_surface("$(X)O-", -RT25_HD * log(10.0^-8.9); planes = pl(-1.0))
    support = SurfaceSupport("s$(X)", nothing, FixedSurfaceArea(60.0))
    return SiteFamily(X, free, [prot, depr]; capacity = TotalSiteAmount(2.0e-3), support, model)
end

function _hd_system(; gas_eos::Bool)
    cem = Dict(symbol(s) => s for s in build_species(datapath("cemdata18-thermofun.json"); verbose = false))
    aqueous = [cem[s] for s in ("H2O@", "H+", "OH-", "Na+", "Cl-", "Ca+2", "CO2@", "HCO3-", "CO3-2", "Ca(CO3)@")]
    gases = gas_eos ? [peng_robinson(cem["CO2"]), peng_robinson(cem["N2"])] : [cem["CO2"], cem["N2"]]
    afm = SolidSolutionPhase("AFm", [cem["monosulphate12"], cem["C4AClH10"]]; model = RedlichKisterModel(a0 = 3000.0, a1 = 400.0))
    others = [cem["SO4-2"], cem["AlO2-"]]
    families = [
        _hd_family("Xa", IdealSiteMixing()),
        _hd_family("Xb", ConstantCapacitance(; C = 1.2, area = 60.0)),
        _hd_family("Xc", DiffuseLayer(; area = 60.0)),
        _hd_family("Xd", ChargePlanes(; area = 60.0, C1 = 1.1); planes = true),
    ]
    members = reduce(vcat, [vcat([f.free_site], f.complexes) for f in families])
    species = vcat(aqueous, others, gases, end_members(afm), members)
    primaries = vcat(
        ["H2O@", "H+", "Na+", "Cl-", "Ca+2", "CO2@", "SO4-2", "AlO2-", "N2"],
        [symbol(f.free_site) for f in families],
    )
    return ChemicalSystem(species, primaries; solid_solutions = [afm], site_families = families)
end

# The state: every amount drawn around a pore solution, two of them at zero so
# that the activity floor takes over and their rows go flat.
function _hd_amounts(cs, rng)
    at(s) = findfirst(sp -> symbol(sp) == s, cs.species)
    n = 1.0e-3 .* exp.(2 .* randn(rng, length(cs.species)))
    n[at("H2O@")] = 55.5
    n[at("Na+")] = 0.3
    n[at("Cl-")] = 0.25
    n[at("Ca(CO3)@")] = 0.0
    n[at("XbO-")] = 0.0
    return n
end

function _hd_params(cs, n; T = 298.15)
    st = ChemicalState(cs, n .* u"mol"; T = 298.15u"K")
    p = ChemistryLab._build_params(st)
    return merge(p, (T = T,))
end

_jacobian_diagonal(lna, n, p) = diag(ForwardDiff.jacobian(x -> lna(x, p), n))

@testsection "the diagonal of the Gibbs Hessian, without the Jacobian" begin
    rng = Xoshiro(20261103)
    models = (
        HKFActivityModel(),
        DaviesActivityModel(),
        TruesdellJonesActivityModel(; parameters = Dict("Na+" => (4.0, 0.075), "Cl-" => (3.5, 0.015), "Ca+2" => (5.0, 0.165))),
        LLNLActivityModel(;
            temperatures = [0.0, 25.0, 60.0], A = [0.49, 0.51, 0.55], B = [0.325, 0.329, 0.335],
            Bdot = [0.039, 0.041, 0.044], co2 = [0.1, 0.01, 1.0e-4, -2.0e-3, 1.0e-6],
            sizes = Dict("Na+" => 4.0, "Cl-" => 3.0, "Ca+2" => 6.0), co2_species = ["CO2@"],
        ),
    )

    @testset "equal to the diagonal of the Jacobian, $(nameof(typeof(m))), $(eos ? "real" : "ideal") gas" for m in models, eos in (false, true)
        cs = _hd_system(; gas_eos = eos)
        lna = activity_model(cs, m)
        @test ChemistryLab._has_lna_diagonal(lna)
        for _ in 1:5
            n = _hd_amounts(cs, rng)
            p = _hd_params(cs, n)
            ref = _jacobian_diagonal(lna, n, p)
            d = ChemistryLab._lna_diagonal!(similar(n), lna, n, p)
            @test d ≈ ref rtol = 1.0e-12
            # the rows the floor flattens, exactly
            @test findall(iszero, d) == findall(iszero, ref)
        end
    end

    @testset "derivatives through the parameters, nested under their own tag" begin
        # The temperature carried as a dual number: the diagonal and its
        # derivative with respect to `T`, against those of the Jacobian.
        cs = _hd_system(; gas_eos = true)
        lna = activity_model(cs, HKFActivityModel())
        n = _hd_amounts(cs, rng)
        mine(T) = ChemistryLab._lna_diagonal!(zeros(typeof(T), length(n)), lna, n, _hd_params(cs, n; T))
        theirs(T) = _jacobian_diagonal(lna, n, _hd_params(cs, n; T))
        @test ForwardDiff.derivative(mine, 310.0) ≈ ForwardDiff.derivative(theirs, 310.0) rtol = 1.0e-10
    end

    @testset "the hook the interior point reads" begin
        cs = _hd_system(; gas_eos = false)
        μ = build_potentials(cs, HKFActivityModel())
        n = _hd_amounts(cs, rng)
        p = _hd_params(cs, n)
        hdiag = ChemistryLab._hessian_diagonal_hook(μ, p)
        hf = hdiag(similar(n), n)
        # the diagonal of ∂μ/∂n, floored at zero as the interior point floors it
        @test hf ≈ max.(diag(ForwardDiff.jacobian(x -> μ(x, p), n)), 0.0) rtol = 1.0e-12
        @test all(>=(0), hf)
        # and none for an activity closure it has no method for, nor for
        # anything but potentials
        @test ChemistryLab._hessian_diagonal_hook(ChemistryLab._Potentials((n, p) -> n), p) === nothing
        @test ChemistryLab._hessian_diagonal_hook((n, p) -> n, p) === nothing
    end

    @testset "the Hessian the active set reads" begin
        cs = _hd_system(; gas_eos = false)
        μ = build_potentials(cs, HKFActivityModel())
        n = _hd_amounts(cs, rng)
        p = _hd_params(cs, n)
        ns = length(n)
        H = ChemistryLab._hessian_hook(μ, p)(zeros(ns, ns), n)
        J = ForwardDiff.jacobian(x -> μ(x, p), n)
        at(s) = findfirst(sp -> symbol(sp) == s, cs.species)
        # the diagonal, exact
        @test diag(H) ≈ diag(J) rtol = 1.0e-12
        # off it, the coupling of an ideal mixture: exact within the ideal gas,
        # whose members' log activities are their log mole fractions
        gas = cs.idx_gas
        @test H[gas, gas] ≈ J[gas, gas] rtol = 1.0e-12
        @test H[at("CO2"), at("N2")] ≈ -1 / sum(n[gas])
        # between a solute and the solvent, and nowhere else across phases
        w = at("H2O@")
        @test H[at("Na+"), w] == H[w, at("Na+")] == -1 / n[w]
        @test iszero(H[at("Na+"), at("Cl-")])
        @test iszero(H[at("Na+"), gas[1]])
        # an amount on the floor moves nothing: its column is zero, as the
        # Jacobian's is
        for s in ("Ca(CO3)@", "XbO-")
            @test iszero(H[:, at(s)]) && iszero(J[:, at(s)])
        end
        # with the temperature a dual number, under its own tag
        mine(T) = diag(ChemistryLab._hessian_hook(μ, _hd_params(cs, n; T))(zeros(typeof(T), ns, ns), n))
        theirs(T) = diag(ForwardDiff.jacobian(x -> μ(x, _hd_params(cs, n; T)), n))
        @test ForwardDiff.derivative(mine, 310.0) ≈ ForwardDiff.derivative(theirs, 310.0) rtol = 1.0e-10
        # none without a diagonal to build it on
        @test ChemistryLab._hessian_hook(ChemistryLab._Potentials((n, p) -> n), p) === nothing
        @test ChemistryLab._hessian_hook((n, p) -> n, p) === nothing
    end
end
