# The extensions whose dependencies the main suite does not load, each tested in
# this environment of their own. SymbolicNumericIntegration brings some fifty
# packages, ModelingToolkitBase among them, that the suite has no use for; the
# coverage job of the CI runs this file after the suite.
#     julia --project=test/extensions test/extensions/runtests.jl
using ChemistryLab, DynamicQuantities, SymbolicNumericIntegration, Test
import SymPyPythonCall

@testset "SymbolicNumericIntegrationExt: a model from its heat capacity" begin
    @test Base.get_extension(ChemistryLab, :SymbolicNumericIntegrationExt) !== nothing

    # H = ∫Cp dT, S = ∫Cp/T dT and G = −∫S dT, integrated symbolically: their
    # differences between two temperatures are those of the closed forms.
    name = :test_linear_model
    factories = add_thermo_model(name, :(c₀ + c₁ * T))
    try
        @test haskey(ChemistryLab.THERMO_MODELS, name) && haskey(ChemistryLab.THERMO_FACTORIES, name)
        @test issetequal(keys(factories), (:Cp, :H, :S, :G))
        c₀, c₁, T₁, T₂ = 2.0, 0.5, 300.0, 400.0
        at(k, T) = factories[k](; c₀, c₁)(; T)
        @test at(:Cp, T₁) ≈ c₀ + c₁ * T₁
        @test at(:H, T₂) - at(:H, T₁) ≈ c₀ * (T₂ - T₁) + c₁ * (T₂^2 - T₁^2) / 2
        @test at(:S, T₂) - at(:S, T₁) ≈ c₀ * log(T₂ / T₁) + c₁ * (T₂ - T₁)
        F(T) = c₀ * (T * log(T) - T) + c₁ * T^2 / 2      # ∫S dT
        @test at(:G, T₂) - at(:G, T₁) ≈ -(F(T₂) - F(T₁))
    finally
        delete!(ChemistryLab.THERMO_MODELS, name)
        delete!(ChemistryLab.THERMO_FACTORIES, name)
    end

    # With units: those given, and the functions' own completed as J/mol/K for
    # Cp and S and J/mol for H and G.
    name = :test_linear_model_units
    factories = add_thermo_model(name, :(a + b * T), [:T => u"K", :a => u"J/mol/K", :b => u"J/(mol*K^2)"])
    try
        @test ChemistryLab.get_unit(factories[:Cp]) == u"J/mol/K"
        @test ChemistryLab.get_unit(factories[:S]) == u"J/mol/K"
        @test ChemistryLab.get_unit(factories[:H]) == u"J/mol"
        @test ChemistryLab.get_unit(factories[:G]) == u"J/mol"
        @test ChemistryLab.get_unit(factories[:H], :b) == u"J/(mol*K^2)"
    finally
        delete!(ChemistryLab.THERMO_MODELS, name)
        delete!(ChemistryLab.THERMO_FACTORIES, name)
    end
end

@testset "Reactions with SymPy coefficients" begin
    # Stoichiometric coefficients may be SymPy symbols. Their matrix rank came
    # back as a `Sym`, and `Int(x)` on it fails on SymPyPythonCall, so building
    # such a reaction threw. The reference is the combustion balanced by hand:
    # CₙH₂ₙ₊₂ + (3n+1)/2 O₂ = (n+1) H₂O + n CO₂.
    n = SymPyPythonCall.symbols("n", real = true)
    a = SymPyPythonCall.symbols("a", real = true)
    @test ChemistryLab._float_rank([a 1; 1 a]) isa Int
    alkane = Species(:C => n, :H => 2n + 2)
    O₂, H₂O, CO₂ = Species.(split("O₂ H₂O CO₂"))
    r = Reaction([alkane, O₂], [H₂O, CO₂])
    @test SymPyPythonCall.simplify(r[O₂] + (3n + 1) / 2) == 0
    @test SymPyPythonCall.simplify(r[H₂O] - (n + 1)) == 0
    @test SymPyPythonCall.simplify(r[CO₂] - n) == 0
end
