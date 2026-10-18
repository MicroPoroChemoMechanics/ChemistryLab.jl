# The diagnostics of a system and a recipe: each finding raised by a system or a
# recipe built to show it, and absent where it does not apply.

using ChemistryLab, DynamicQuantities, Test

@testset "Diagnostics of a system" begin
    substances = build_species(datapath("cemdata18-thermofun.json"); verbose = false)
    db = Dict(symbol(s) => s for s in substances)
    list = "Portland paste (Lothenbach and Winnefeld 2006)"
    cs = phase_list_system(list, substances)
    codes(ds) = [d.code for d in ds]
    only_code(ds, c) = only(d for d in ds if d.code === c)

    # the published list alone, without a recipe or a temperature: nothing to say
    @test isempty(diagnose_system(cs))

    # a cement whose constituents the list holds, at 20 °C: the clinker phases are
    # read below the 25 °C their data start at, as reactants that dissolve; the
    # other reactants of the list are such only when named
    cem = Recipe(material_template("CEM I 52.5 N (Gruyaert 2010), Bogue", cs) => 1.0; w_b = 0.5)
    ds = diagnose_system(cs; recipe = cem, reactants = phase_list(list).reactants)
    @test codes(ds) == [:reactant_temperature_range]
    d = only(ds)
    @test d.level === :info && "C3S" in d.subjects && isempty(d.fix)
    alone = only_code(diagnose_system(cs; recipe = cem), :temperature_range)
    @test !isempty(alone.subjects) && all(in(phase_list(list).reactants), alone.subjects)

    # minor oxides whose elements no species carries are set aside
    pc = Recipe(material_template("PC (Durdzinski 2017), Bogue", cs) => 1.0; w_b = 0.5)
    d = only_code(diagnose_system(cs; recipe = pc), :not_in_system)
    @test d.level === :warning && "minor oxides: TiO2" in d.subjects

    # at 80 °C, hydrates are read beyond their data; the reactants are not reported then
    hot = Recipe(material_template("CEM I 52.5 N (Gruyaert 2010), Bogue", cs) => 1.0; w_b = 0.5, T = T_ZERO_CELSIUS_Q + 80u"K")
    d = only_code(diagnose_system(cs; recipe = hot), :temperature_range)
    @test d.level === :warning && "ettringite" in d.subjects && "C4AH13" in d.subjects
    @test :reactant_temperature_range ∉ codes(diagnose_system(cs; recipe = hot))
    @test codes(diagnose_system(cs; T = 80u"K" + T_ZERO_CELSIUS_Q)) == [:temperature_range]

    # the slow Al-Si hydrogarnet: a warning at ambient temperature, information above
    with_garnet = phase_list_system(list, substances; add = ["C3AS0.84H4.32"])
    d = only_code(diagnose_system(with_garnet), :slow_hydrogarnet)
    @test d.level === :warning && d.subjects == ["C3AS0.84H4.32"]
    @test only_code(diagnose_system(with_garnet; T = 333.15), :slow_hydrogarnet).level === :info

    # ettringite twice, as itself and as a third of itself
    twice = phase_list_system(list, substances; add = ["ettringite03_ss"])
    @test only_code(diagnose_system(twice), :declared_twice).subjects == ["ettringite", "ettringite03_ss"]

    # a crystal without a molar volume or an enthalpy of formation
    bare = Species("CaO"; aggregate_state = AS_CRYSTAL, symbol = "bare_lime")
    small = ChemicalSystem(vcat(speciation(substances, ["Portlandite"]; aggregate_state = [AS_AQUEOUS]), [bare]), CEMDATA_PRIMARIES)
    ds = diagnose_system(small)
    @test only_code(ds, :missing_volume).subjects == ["bare_lime"]
    @test only_code(ds, :missing_enthalpy).subjects == ["bare_lime"]

    # SIT: the pairs of ions without a coefficient, written as pairs
    sit = SITActivityModel(; parameters = build_sit_parameters(datapath("sit.dat")))
    d = only_code(diagnose_system(cs; model = sit), :missing_sit_pairs)
    @test d.level === :warning && all(occursin("–", p) for p in d.subjects)
    @test length(d.subjects) == length(missing_epsilon_pairs(cs, sit))
    @test isempty(diagnose_system(cs; model = cemdata18_activity_model(:KOH)))

    # a system beyond the size given
    @test only_code(diagnose_system(cs; large = 10), :large_system).level === :info
    @test isempty(diagnose_system(cs; large = length(cs.species)))

    # a diagnostic prints as a sentence
    @test startswith(sprint(show, d), "warning (missing_sit_pairs): ")
    @test sprint(show, SystemDiagnostic(:info, :x, "said.", "", String[])) == "info (x): said."
end
