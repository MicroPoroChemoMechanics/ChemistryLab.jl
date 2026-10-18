# The generated script: it parses, its compute part is sections 1 to 5, and the
# values the model holds are the ones the code passes, read back from the
# parsed expressions rather than from the text the generator wrote.

using Dates

m0 = CLS.default_model()
g = CLS.generate(m0; version = v"0.1.0", date = Date(2026, 10, 18))

"Every call to `f` in an expression tree."
function calls(ex, f::Symbol, found = Expr[])
    ex isa Expr || return found
    ex.head === :call && ex.args[1] === f && push!(found, ex)
    foreach(a -> calls(a, f, found), ex.args)
    return found
end
"The value of keyword `name` in a call, written after `;` or not."
function kwarg(call::Expr, name::Symbol)
    params = [a for a in call.args if a isa Expr && a.head === :parameters]
    kws = vcat((p.args for p in params)..., call.args)
    found = [a for a in kws if a isa Expr && a.head === :kw && a.args[1] === name]
    isempty(found) && return any(a -> a === name, kws) ? name : nothing   # `; w_b` stands for `w_b = w_b`
    return only(found).args[2]
end
"Every assignment to `name`."
assignments(ex, name, found = Any[]) = (ex isa Expr && (ex.head === :(=) && ex.args[1] === name && push!(found, ex.args[2]); foreach(a -> assignments(a, name, found), ex.args)); found)
parses(src) = !any(a -> a isa Expr && a.head in (:error, :incomplete), Meta.parseall(src).args)
compute(m) = Meta.parseall(CLS.generate(m).compute)
with(f) = (m = CLS.model_from_dict(CLS.to_dict(m0)); f(m); m)

@testset "the default script" begin
    @test parses(g.source)
    @test startswith(g.source, g.compute)
    @test occursin("## 5. Calculation", g.compute) && !occursin("## 6.", g.compute)
    @test all(occursin("## $k. ", g.source) for k in 1:7)
    @test occursin("Written by ChemistryLab Studio 0.1.0 on 2026-10-18.", g.source)
    @test occursin("julia --project=@chemistrylab-studio cem_i_52_5_n_paste.jl", g.source)
    c = Meta.parseall(g.compute)
    @test only(assignments(c, :w_b)) == m0.mix.w_b
    @test only(calls(c, :Recipe)).args[end] == :(cem_i_52_5_n => 1.0)
    @test kwarg(only(calls(c, :Recipe)), :w_b) === :w_b
    extents = only(calls(c, :with_extents))
    @test Dict(e.args[2] => e.args[3] for e in extents.args[3].args[2:end]) ==
        Dict(k => v.value for (k, v) in m0.materials[1].extents)
    @test only(calls(c, :phase_list_system)).args[2] == m0.system.phase_list
    @test only(calls(c, :cemdata18_activity_model)).args[2] == QuoteNode(:KOH)
    d = only(calls(c, :diagnose_system))
    @test kwarg(d, :t) === nothing                                   # constant laws read no age
    @test kwarg(d, :reactants) == :(phase_list($(m0.system.phase_list)).reactants)
    @test kwarg(only(calls(c, :equilibrate_certified)), :t) === nothing
    @test occursin("masses = phase_masses(state)", g.source) && occursin("bar(", g.source)
end

@testset "materials" begin
    logistic = CLS.LogisticExtentChoice(0.8, 3.0, 1.5, 0.0)
    m = with(
        m -> (
            m.materials = CLS.MaterialChoice[
                CLS.TemplateMaterial(m0.materials[1].template, 0.7, ["C3S" => logistic, "C3A" => CLS.ParrottKillohExtentChoice(400.0, true)], nothing),
                CLS.TemplateMaterial("S1 slag (Durdzinski 2017)", 0.3, Pair{String, CLS.ExtentChoice}[], CLS.ConstantExtentChoice(0.4)),
            ]
        )
    )
    c = compute(m)
    @test only(calls(c, :Recipe)).args[(end - 1):end] == [:(cem_i_52_5_n => 0.7), :(s1_slag => 0.3)]
    @test kwarg(only(calls(c, :LogisticExtent)), :half_time) == 3.0
    pk = only(calls(c, :ParrottKillohExtent))
    @test pk.args[end] == "C3A" && kwarg(pk, :blaine) == 400.0 && kwarg(pk, :α_max) == :(powers_alpha_max(w_b))
    slag = only(e for e in calls(c, :with_extents) if e.args[3] == :s1_slag)
    @test kwarg(slag, :material_extent) == 0.4
    # an age is read by the solve and the diagnosis once a law depends on it
    @test kwarg(only(calls(c, :equilibrate_certified)), :t) == 28.0
    @test kwarg(only(calls(c, :diagnose_system)), :t) == 28.0
    # no Powers cap, a table, a capped law
    m = with(
        m -> (
            m.materials[1] = CLS.TemplateMaterial(
                m0.materials[1].template, 1.0, [
                    "C3A" => CLS.ParrottKillohExtentChoice(385.0, false), "C4AF" => CLS.TabulatedExtentChoice([1.0, 7.0], [0.2, 0.5]),
                    "C2S" => CLS.CappedExtentChoice(logistic, 0.6),
                ], nothing
            )
        )
    )
    c = compute(m)
    @test kwarg(only(calls(c, :ParrottKillohExtent)), :α_max) === nothing
    @test only(calls(c, :TabulatedExtent)).args[2:3] == [:([1.0, 7.0]), :([0.2, 0.5])]
    @test only(calls(c, :CappedExtent)).args[end] == 0.6
    # an oxide analysis by Bogue's calculation, and as oxides with a density
    rows = [CLS.AmountRow("CaO", 64.0, "%"), CLS.AmountRow("SiO2", 21.0, "%")]
    c = compute(with(m -> (m.materials[1] = CLS.OxideMaterial("my cement", "cement", 1.0, rows, "bogue", "Anh", 0.0, Pair{String, CLS.ExtentChoice}[], nothing))))
    @test only(calls(c, :bogue)).args[(end - 1):end] == [:my_cement_oxides, :cs] && kwarg(only(calls(c, :bogue)), :sulfate) == "Anh"
    @test only(assignments(c, :my_cement_oxides)) == :(Dict("CaO" => 0.64, "SiO2" => 0.21))
    c = compute(with(m -> (m.materials[1] = CLS.OxideMaterial("slag", "scm", 1.0, rows, "oxides", "Gp", 2.9, Pair{String, CLS.ExtentChoice}[], nothing))))
    om = only(calls(c, :oxide_material))
    @test om.args[(end - 1):end] == ["slag", :slag_oxides] && kwarg(om, :density) == 2.9 && kwarg(om, :kind) == QuoteNode(:scm)
    c = compute(with(m -> (m.materials[1] = CLS.OxideMaterial("slag", "scm", 1.0, rows, "oxides", "Gp", 0.0, Pair{String, CLS.ExtentChoice}[], nothing))))
    @test kwarg(only(calls(c, :oxide_material)), :density) === nothing
    # a mineral composition
    c = compute(with(m -> (m.materials[1] = CLS.MineralMaterial("clinker", "cement", 1.0, [CLS.AmountRow("C3S", 65.0, "%")], Pair{String, CLS.ExtentChoice}[], nothing))))
    mc = only(calls(c, :MineralConstituent))
    @test mc.args[end] == :(cs["C3S"]) && kwarg(mc, :mass_fraction) == 0.65
end

@testset "database, system, activity, mix, constraint" begin
    c = compute(with(m -> (m.database = CLS.ImportedDatabase("/data/phreeqc.dat", "phreeqc"))))
    @test only(calls(c, :import_database)).args[end] == "/data/phreeqc.dat"
    @test kwarg(only(calls(c, :import_database)), :format) == QuoteNode(:phreeqc)
    sys = CLS.PhaseListSystem(m0.system.phase_list, ["C4AClH10"], ["straetlingite"], "CNASH_ss", ["CH4@"])
    p = only(calls(compute(with(m -> (m.system = sys))), :phase_list_system))
    @test kwarg(p, :add) == :(["C4AClH10"]) && kwarg(p, :remove) == :(["straetlingite"])
    @test kwarg(p, :replace) == :(Dict("CSHQ" => "CNASH_ss")) && kwarg(p, :exclude_aqueous) == :(["CH4@"])
    c = compute(with(m -> (m.system = CLS.CustomSystem(["C3S", "Portlandite"], ["CSHQ"], ["H2@"]))))
    @test only(calls(c, :ChemicalSystem)).args[(end - 1):end] == [:species, :CEMDATA_PRIMARIES]
    @test kwarg(only(calls(c, :diagnose_system)), :reactants) === nothing      # no list, no reactants of a list
    @test only(calls(compute(with(m -> (m.activity = CLS.Cemdata18Activity("NaOH", false)))), :cemdata18_activity_model)).args[end] == QuoteNode(:NaOH)
    @test kwarg(only(calls(compute(with(m -> (m.activity = CLS.DaviesActivity(0.25)))), :DaviesActivityModel)), :b) == 0.25
    @test kwarg(only(calls(compute(with(m -> (m.activity = CLS.PitzerActivity("zero")))), :PitzerActivityModel)), :missing_pairs) == QuoteNode(:zero)
    r = only(calls(compute(with(m -> (m.mix = CLS.Mix(0.4, 35.0, 250.0, 2.0)))), :Recipe))
    # a unit macro carries a line number: compare as text
    @test string(kwarg(r, :P)) == string(:(2.0u"bar")) && string(kwarg(r, :binder_mass)) == string(:(250.0u"g"))
    c = compute(with(m -> (m.constraint = CLS.FixedpHChoice(12.0))))
    @test kwarg(only(calls(c, :equilibrate_certified)), :constraint) == :(FixedpH(12.0))
    c = compute(with(m -> (m.constraint = CLS.SaturatedCuringChoice())))
    @test only(assignments(c, :fresh)) == :(budget(recipe, cs).state)
    @test kwarg(only(calls(c, :equilibrate_certified)), :constraint) == :(SaturatedCuring(; reference = fresh))
    c = compute(with(m -> (m.constraint = CLS.FixedActivityChoice("CO2@", 1.0e-5))))
    @test kwarg(only(calls(c, :equilibrate_certified)), :constraint) == :(FixedActivity("CO2@", 1.0e-5))
end

@testset "calculations" begin
    logistic = ["C3S" => CLS.LogisticExtentChoice(0.8, 3.0, 1.5, 0.0)]
    timed(m) = (m.materials[1] = CLS.TemplateMaterial(m0.materials[1].template, 1.0, logistic, nothing))
    c = compute(with(m -> (timed(m); m.calculation = CLS.HydrationPath([1.0, 7.0]))))
    @test only(assignments(c, :ages)) == :([1.0, 7.0])
    @test only(calls(c, :hydrate)).args[(end - 2):end] == [:recipe, :cs, :ages]
    @test kwarg(only(calls(c, :diagnose_system)), :t) == 1.0          # the first age
    c = compute(with(m -> (timed(m); m.calculation = CLS.HydrationPath(Float64[]))))
    @test kwarg(only(calls(c, :diagnose_system)), :t) === nothing     # no age yet
    c = compute(with(m -> (timed(m); m.calculation = CLS.HydrationPath([1.0, 7.0]); m.constraint = CLS.FixedpHChoice(12.0))))
    @test kwarg(only(calls(c, :hydrate)), :constraint) == :(FixedpH(12.0))
    c = compute(with(m -> (m.calculation = CLS.BlendPath("S1 slag (Durdzinski 2017)", 0.4, [0.0, 0.3], 90.0))))
    @test only(assignments(c, :fractions)) == :([0.0, 0.3])
    @test kwarg(only(calls(c, :blend)), :t) === nothing                # constant laws: no age
    c = compute(with(m -> (timed(m); m.calculation = CLS.BlendPath("S1 slag (Durdzinski 2017)", 0.4, [0.0, 0.3], 90.0))))
    @test kwarg(only(calls(c, :blend)), :t) == 90.0
    c = compute(with(m -> (m.calculation = CLS.CarbonationPath([0.0, 0.5], 28.0))))
    @test only(calls(c, :carbonate)).args[end] == :([0.0, 0.5])
    c = compute(with(m -> (m.calculation = CLS.SaltPath("NaCl", [0.0, 0.1], 28.0); m.constraint = CLS.FixedpHChoice(12.0))))
    @test only(calls(c, :add_salt)).args[end - 1] == "NaCl"
    @test kwarg(only(calls(c, :equilibrate_certified)), :constraint) == :(FixedpH(12.0))   # on the paste the process starts from
    @test only(calls(compute(with(m -> (m.calculation = CLS.TitrationPath("HCl", [0.1], 28.0)))), :titrate)).args[end - 1] == "HCl"
    l = only(calls(compute(with(m -> (m.calculation = CLS.LeachingPath(5, 1000.0, 28.0)))), :leach))
    @test l.args[end] == 5 && string(kwarg(l, :renewal)) == string(:(1000.0u"g"))
    c = compute(with(m -> (m.calculation = CLS.AttackPath([CLS.AmountRow("Na2SO4", 50.0, "mmol/L"), CLS.AmountRow("MgSO4", 0.01, "mol/L")], 250.0, 4, 28.0))))
    @test only(assignments(c, :V)) == 250.0
    @test kwarg(only(calls(c, :leach)), :solution) === :solution
    @test 0.05 in [x.args[2] for x in calls(c, :.*) if x.args[2] isa Real]
    @test 0.01 in [x.args[2] for x in calls(c, :.*) if x.args[2] isa Real]
    c = compute(with(m -> (m.calculation = CLS.TemperatureSweep([5.0, 40.0], 28.0))))
    @test only(calls(c, :ProcessResult)).args[2:end] == [QuoteNode(:temperature), :temperatures, :states]
    c = compute(with(m -> (m.calculation = CLS.WaterSweep([0.4, 0.6], 28.0))))
    @test only(calls(c, :ProcessResult)).args[2] == QuoteNode(:w_b)
    # a path is tabulated and drawn as a path
    s = CLS.generate(with(m -> (m.calculation = CLS.CarbonationPath([0.0, 0.5], 28.0)))).source
    @test occursin("process_table(path)", s) && occursin("areaplot(", s)
end

@testset "kinetics" begin
    law(kind; kw...) = CLS.choice_from_dict(CLS.RateChoice, Dict{String, Any}("kind" => kind, (string(k) => v for (k, v) in kw)...))
    run(rates; kw...) = CLS.KineticsCalc(get(kw, :horizon, 28.0), get(kw, :ages, [1.0, 28.0]), rates, get(kw, :cal, CLS.IsothermalChoice()), 1.0e-6)
    pk = law("parrott_killoh_avrami")
    m = with(m -> (m.calculation = run(["C3S" => pk, "Gp" => law("first_order"; k = 0.5)])))
    @test isempty(CLS.problems(m))
    c = compute(m)
    @test !isempty(calls(c, :KineticsProblem)) && only(calls(c, :integrate)).args[2] === :problem
    p = only(calls(c, :parrott_killoh_avrami))
    @test p.args[(end - 1):end] == [:PK84_PARAMS_C3S, "C3S"]
    @test kwarg(p, :blaine) == 385.0 && kwarg(p, :w_c) === :w_c && kwarg(p, :α_max) == :(powers_alpha_max(w_c))
    @test length(assignments(c, :w_c)) == 1
    @test only(assignments(c, :horizon)) == 28.0 && only(assignments(c, :ages)).args[end] == :([1.0, 28.0])
    # the first-order law, in mol/s from a constant per day
    f = only(calls(c, :KineticFunc))
    @test f.args[2].args[2] isa Expr && occursin("0.5 / 86400", string(f.args[2]))
    @test isempty(calls(c, :SemiAdiabaticCalorimeter)) && kwarg(only(calls(c, :KineticsProblem)), :calorimeter) === nothing
    # no clinker law that reads the water, no `w_c`
    c = compute(with(m -> (m.calculation = run(["C3S" => CLS.PKAvramiRate(0.0, false, false)]))))
    p = only(calls(c, :parrott_killoh_avrami))
    @test isempty(assignments(c, :w_c)) && kwarg(p, :blaine) === nothing && kwarg(p, :w_c) === nothing
    # the laws of Waller, the user's own, a semi-adiabatic cell
    rates = ["C3S" => law("waller"), "C2S" => law("waller"; preset = "silica_fume", tau = 40.0, blaine = 2000.0), "Gp" => law("custom_rate")]
    s = CLS.generate(with(m -> (m.calculation = run(rates; cal = CLS.SemiAdiabaticChoice(900.0, 0.3, 20.0))))).source
    @test parses(s)
    c = Meta.parseall(s)
    w = calls(c, :waller)
    @test w[1].args[(end - 1):end] == [:WALLER_PARAMS_FLY_ASH, "C3S"]
    @test string(w[2].args[end - 1]) == string(:(merge(WALLER_PARAMS_SILICA_FUME, (τ = 40.0u"d",)))) && kwarg(w[2], :blaine) == 2000.0
    @test only(assignments(c, :key)) == "Gp"
    cal = only(calls(c, :SemiAdiabaticCalorimeter))
    @test string(kwarg(cal, :Cp)) == string(:(900.0u"J/K")) && kwarg(cal, :T0) === :T
    @test kwarg(only(calls(c, :KineticsProblem)), :calorimeter) === :calorimeter
    @test occursin("degrees_of_hydration(solution, problem; times = ages .* 86400)", s)
    @test occursin("_degrees.png", s)
    # what keeps a kinetic run from running
    codes(m) = [p.code for p in CLS.problems(m)]
    @test codes(with(m -> (m.calculation = run(Pair{String, CLS.RateChoice}[])))) == ["no_rate"]
    @test codes(with(m -> (m.calculation = run(["C3S" => pk]; ages = [1.0, 90.0])))) == ["ages_beyond_horizon"]
    @test codes(with(m -> (m.calculation = run(["C3S" => pk]; ages = Float64[])))) == ["empty_path"]
    @test codes(with(m -> (m.calculation = run(["Gp" => pk])))) == ["pk_rate_phase"]
    @test codes(with(m -> (m.calculation = run(["Gp" => CLS.CustomRate("n[key] *")])))) == ["rate_syntax"]
    @test codes(with(m -> (m.calculation = run(["Gp" => CLS.CustomRate("  ")])))) == ["rate_syntax"]
    @test codes(with(m -> (m.calculation = run(["C3S" => pk]); m.constraint = CLS.FixedpHChoice(12.0)))) == ["kinetics_constraint"]
    @test codes(with(m -> (m.calculation = run(["C3S" => CLS.FirstOrderRate(-1.0)])))) == ["out_of_range"]
    # the laws of the degrees of reaction are read at the start of the run
    logistic = ["C2S" => CLS.LogisticExtentChoice(0.8, 3.0, 1.5, 0.0)]
    c = compute(with(m -> (m.materials[1] = CLS.TemplateMaterial(m0.materials[1].template, 1.0, logistic, nothing); m.calculation = run(["C3S" => pk]))))
    @test kwarg(only(calls(c, :diagnose_system)), :t) == 0.0
end

@testset "names and numbers" begin
    @test CLS.material_name("S1 slag (Durdzinski 2017)", String[]) == "s1_slag"
    @test CLS.material_name("PC (Durdzinski 2017), Bogue", ["pc"]) == "pc_2"
    @test CLS.material_name("PC (x)", ["pc", "pc_2"]) == "pc_3"
    @test CLS.material_name("42 (x)", String[]) == "material_42"
    @test CLS.material_name("Recipe", String[]) == "recipe_material"
    @test CLS.material_name("scm", String[]) == "scm_material"
    @test CLS.material_name("(x)", String[]) == "material_"
    @test CLS.julia_number(1) == "1.0"
    @test CLS.julia_number(0.45) == "0.45"
    @test CLS.script_stem("CEM I 52.5 N, at 28 d!") == "cem_i_52_5_n_at_28_d"
    @test CLS.script_stem("…") == "chemistrylab_study"
    @test textwidth(CLS.heading(3, "Chemical system")) == CLS.RULE_WIDTH
    @test CLS.fractions_code([CLS.AmountRow("CaO", 50, "%")]) == "Dict(\"CaO\" => 0.5)"
end
