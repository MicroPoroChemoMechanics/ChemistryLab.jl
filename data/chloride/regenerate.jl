# Build `data/cemdata18-chloride.json`: CEMDATA18 plus a chloride end member of
# CSHQ, fitted on the sorption tests of Hirao et al. (2005).
#
#   julia --project=docs data/chloride/regenerate.jl
#
# The file this writes is VENDORED, and this script exists so that it can be
# reproduced and audited rather than trusted. CEMDATA18 is copied through
# verbatim; one substance, `CSHQ-Cl`, is appended. Its Gibbs energy is the one
# number fitted here (data/chloride/member.jl writes the end member and models
# the tests). Two candidates are fitted, and the one the file ships is the one
# whose binding grows with the Ca/Si of the gel, as measured; the other is
# recorded with the reason it was rejected. README.md in this directory records
# the choice of the data and the sets that were examined and rejected.

import ChemistryLab
using ChemistryLab, DynamicQuantities, OptimaSolver, Printf

const JSON = ChemistryLab.JSON
const HERE = @__DIR__
const SRC = joinpath(pkgdir(ChemistryLab), "data", "cemdata18-thermofun.json")
const OUT = joinpath(pkgdir(ChemistryLab), "data", "cemdata18-chloride.json")

include(joinpath(HERE, "member.jl"))
using .ChlorideMember: AfterNaSiOH, CalciumChloride, SHIPPED, T_REF, calcium_trend, calibration_points,
    composition, forward, member_entry, member_formula, member_species, member_symbol, standard_properties

# The activity model of the CEMDATA18 pages: the Debye-Hückel limiting law with a
# B-dot term, the B-dot identified from the activity coefficients a GEM-Selektor
# run printed for |z| = 1 and 2 (test/reference/gems_cemdata18_portland.json).
function gems_model()
    g = JSON.parsefile(joinpath(pkgdir(ChemistryLab), "test", "reference", "gems_cemdata18_portland.json"))
    lg1, lg2 = log10(g["gamma"]["z1"]), log10(g["gamma"]["z2"])
    return HKFActivityModel(å = 0.0, Ḃ = (lg1 + (lg1 - lg2) / 3) / g["ionic_strength_mol_per_kg"], Kₙ = 0.0)
end

# A golden-section search, which is all a single bracketed parameter needs.
function golden(f, a, b; tol = 1.0)
    φ = (sqrt(5) - 1) / 2
    c, d = b - φ * (b - a), a + φ * (b - a)
    fc, fd = f(c), f(d)
    while b - a > tol
        if fc < fd
            b, d, fd = d, c, fc
            c = b - φ * (b - a)
            fc = f(c)
        else
            a, c, fc = c, d, fd
            d = a + φ * (b - a)
            fd = f(d)
        end
    end
    return (a + b) / 2
end

# The brackets, from a scan of each candidate: below them the gel takes all of
# the chloride, above them none.
bracket(::AfterNaSiOH) = (-20.0e3, 0.0)
bracket(::CalciumChloride) = (-15.0e3, 5.0e3)

function fit(c, subs, db, model, water_fraction)
    f = forward(c, subs, db; model, water_fraction)
    y = calibration_points().bound
    ssr(δ) = sum(abs2, f(δ).bound .- y)
    δ = golden(ssr, bracket(c)...; tol = 1.0)
    r = f(δ)
    return (; δ, f, r, rmse = sqrt(sum(abs2, r.bound .- y) / length(y)))
end

function main()
    db = JSON.parsefile(SRC)
    subs = build_species(SRC; verbose = false)
    byname = Dict(symbol(s) => s for s in subs)
    model = gems_model()
    hg(q) = ustrip(literature_value("HongGlasser1999", q))
    pts = calibration_points()
    w = hg("water_content_high_CaSi")

    # Both candidates, on the same data, and the Ca/Si trend each then gives.
    fits = Dict()
    for c in (AfterNaSiOH(), CalciumChloride())
        ft = fit(c, subs, db, model, w)
        trend = calcium_trend(subs, member_species(member_entry(c, db, byname, ft.δ)); model)
        fits[c] = (; ft..., trend)
        @printf(
            "%-16s δ = %8.0f J/mol  rmse %.4f mmol/g  held per Si at Ca/Si %.2f: %.4f, at %.2f: %.4f\n",
            member_symbol(c), ft.δ, ft.rmse, trend.ratio_low, trend.low, trend.ratio_high, trend.high
        )
    end
    rejected = fits[AfterNaSiOH()]
    rejected.trend.high < rejected.trend.low ||
        error("the candidate written after NaSiOH no longer inverts the Ca/Si trend; revisit the choice.")
    shipped = fits[SHIPPED]
    shipped.trend.high > shipped.trend.low ||
        error("the shipped candidate does not bind more at the higher Ca/Si; nothing written.")

    # The water of the dried gel is not reported. Hong & Glasser adopt Taylor's
    # 18 % for gels of Ca/Si 1.5 and 1.8 and measured 14 % below; the fit is made
    # at the first and repeated at the second, and the difference is part of the
    # uncertainty.
    alt = fit(SHIPPED, subs, db, model, hg("water_content_low_CaSi"))
    δ = shipped.δ
    residuals = shipped.r.bound .- pts.bound

    # The uncertainty, through the constant K = exp(-δ/RT) so that the relative
    # step of `identifiability` is defined: σ(δ) = RT σ(ln K).
    RT = R_GAS * T_REF
    id = identifiability(
        K -> shipped.f(-RT * log(only(K))).bound, [exp(-δ / RT)];
        observed = pts.bound, relstep = 0.01, names = ["K_formation"],
    )
    σ_fit = RT * only(id.stderr)
    σ_water = abs(alt.δ - δ)
    σ = sqrt(σ_fit^2 + σ_water^2)
    p = standard_properties(SHIPPED, byname, δ)

    @printf("shipped: δ = %.0f J/mol (water 18 %%), %.0f J/mol (water 14 %%)\n", δ, alt.δ)
    @printf("σ(δ): %.0f J/mol from the residuals, %.0f from the water, %.0f combined\n", σ_fit, σ_water, σ)
    @printf("gel Ca/Si at portlandite saturation: %.4f\n", shipped.r.r)
    println("  c (mol/L)   bound, measured   model   held by the end member (mmol/g)")
    for k in eachindex(pts.chloride)
        @printf("  %8.4f   %8.4f   %8.4f   %8.4f\n", pts.chloride[k], pts.bound[k], shipped.r.bound[k], shipped.r.held[k])
    end

    sym = member_symbol(SHIPPED)
    existing = Set(String(s["symbol"]) for s in db["substances"])
    sym in existing && error("$sym would overwrite a CEMDATA18 substance; refusing.")
    reaction = join(["$(ν) $(s)" for (s, ν) in composition(SHIPPED)], " + ")
    provenance = Dict(
        "kind" => "fitted",
        "source" => Dict("key" => "Hirao2005", "doi" => literature("Hirao2005").source["doi"]),
        "data" => "Fig. 5, bound chloride per gram of hydrated alite, the points up to 1 mol/L (data/literature/Hirao2005.json)",
        "parameter" => "delta_J_per_mol: the Gibbs energy of forming $sym from $reaction at 293.15 K and 1 bar",
        "delta_J_per_mol" => δ,
        "delta_uncertainty_J_per_mol" => σ,
        "delta_uncertainty_from_residuals" => σ_fit,
        "delta_uncertainty_from_gel_water" => σ_water,
        "gel_water_fraction" => w,
        "chloride_mol_per_L" => pts.chloride,
        "bound_measured_mmol_per_g" => pts.bound,
        "bound_model_mmol_per_g" => shipped.r.bound,
        "residuals_mmol_per_g" => residuals,
        "rmse_mmol_per_g" => shipped.rmse,
        "calcium_trend_held_per_Si" => [shipped.trend.low, shipped.trend.high],
        "calcium_trend_Ca_Si" => [shipped.trend.ratio_low, shipped.trend.ratio_high],
        "estimated" => "S and V: those of the ions it is made of (no entropy or volume of dissolution); heat capacity zero; H from G and S",
        "rejected_candidate" => Dict(
            "formula" => member_formula(AfterNaSiOH()),
            "delta_J_per_mol" => rejected.δ,
            "rmse_mmol_per_g" => rejected.rmse,
            "calcium_trend_held_per_Si" => [rejected.trend.low, rejected.trend.high],
            "calcium_trend_Ca_Si" => [rejected.trend.ratio_low, rejected.trend.ratio_high],
            "reason" => "binds less chloride at the higher Ca/Si, against the measured trend; the silica it carries is favored by the silica activity of a low-Ca/Si gel",
        ),
        "activity_model" => "Debye-Hueckel limiting law with B-dot, a = 0, B-dot identified from test/reference/gems_cemdata18_portland.json",
        "generator" => "data/chloride/regenerate.jl",
    )
    push!(db["substances"], member_entry(SHIPPED, db, byname, δ; provenance))
    db["thermodataset"] = "cemdata18-chloride"
    db["chloride_extension"] = Dict(
        "base" => "cemdata18 (Lothenbach et al. 2019, doi:10.1016/j.cemconres.2018.04.018)",
        "added_substances" => [sym],
        "solid_solution" => "CSHQ_Cl in data/solid_solutions.toml",
        "note" => "CEMDATA18 entries are copied unchanged; nothing is overwritten. The added end member is fitted, not measured: see its chloride_provenance.",
    )
    tmp = OUT * ".tmp"
    open(tmp, "w") do io
        JSON.print(io, db, 2)
    end
    mv(tmp, OUT; force = true)
    println("wrote $OUT ($(length(db["substances"])) substances)")
    return nothing
end

abspath(PROGRAM_FILE) == (@__FILE__) && main()
