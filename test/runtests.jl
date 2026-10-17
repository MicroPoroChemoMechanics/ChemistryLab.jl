using ChemistryLab
using Aqua
using DynamicQuantities
using ForwardDiff
using OptimaSolver
using OrdinaryDiffEq
using OrderedCollections
using Test
using TimerOutputs

# The suite runs in parts, so that CI can run them as parallel jobs.
# `Pkg.test(; test_args = ["equilibrium"])` runs the files of one part, several
# parts may be named, and with no argument every file runs, as it always has.
# A unit is a test file, or the name of a section that includes none (Aqua).
# Every unit belongs to exactly one part, which the last testset checks.
include("test_parts.jl")

function rununit(unit::AbstractString)
    isempty(ARGS) && return true
    for part in ARGS
        haskey(TEST_PARTS, part) ||
            error("unknown test part \"$(part)\"; the parts are $(join(sort!(collect(keys(TEST_PARTS))), ", "))")
    end
    return any(part -> unit in TEST_PARTS[part], ARGS)
end

# A file of a section: run if its part is selected, and timed on its own, so
# that the timer printed at the end gives the cost of every file, which is what
# the parts are balanced on.
function testfile(file::AbstractString)
    rununit(file) || return nothing
    t0 = time()
    @timeit file include(file)
    println("⏱ ", lpad(round(time() - t0; digits = 1), 7), " s  ", file)
    flush(stdout)
    return nothing
end

# In the sections of this file, `include("x.jl")` becomes `testfile("x.jl")`,
# and a section that includes no file is one unit, named after the section. The
# test files use `@testsection` too: theirs are left as they were, run whenever
# their file runs — rewritten or gated, they would be skipped silently.
_route(ex) = ex
function _route(ex::Expr)
    if ex.head === :call && ex.args[1] === :include && length(ex.args) == 2 && ex.args[2] isa String
        return :(testfile($(ex.args[2])))
    end
    return Expr(ex.head, map(_route, ex.args)...)
end
_includes(ex) = false
_includes(ex::Expr) = (ex.head === :call && ex.args[1] === :testfile) || any(_includes, ex.args)

macro testsection(str, block)
    here = __source__.file === nothing ? "" : basename(String(__source__.file))
    toplevel = here == "runtests.jl"
    body = toplevel ? _route(block) : block
    gated = toplevel && !_includes(body) ? :(rununit($str) && $body) : body
    return quote
        @timeit "$($(esc(str)))" begin
            @testset "$($(esc(str)))" begin
                $(esc(gated))
            end
        end
    end
end

reset_timer!()

"""
    aqua_persistent_tasks(m::Module) -> Bool

`Aqua.test_persistent_tasks` in a child process with default bounds checking.

The check precompiles a synthetic package that depends on `m` and waits for a
sentinel written from that package's module body. Julia writes no
precompilation cache when `--check-bounds` is forced, so the body never runs,
the sentinel never appears, and the check fails for a reason that says nothing
about `m`. `Pkg.test()` forces the flag by default.

Running it in a child with `--check-bounds=auto` keeps both halves: the test
suite proper still runs with bounds checking on, and the check still runs.
"""
function aqua_persistent_tasks(m::Module)
    code = """
    using Aqua, $(nameof(m)), Test
    @testset "persistent_tasks" begin
        Aqua.test_persistent_tasks($(nameof(m)))
    end
    """
    cmd = `$(first(Base.julia_cmd())) --check-bounds=auto --startup-file=no --project=$(Base.active_project()) -e $code`
    return success(run(ignorestatus(cmd)))
end

@testsection "Aqua" begin
    Aqua.test_all(
        ChemistryLab;
        # No `piracies` exemption any more. `thermo_factories.jl` used to extend
        # some thirty `Base` math functions to `DynamicQuantities.Quantity`,
        # declared here rather than hidden; those methods are gone, and this
        # check is what keeps them from coming back.
        # Run apart, see `aqua_persistent_tasks` above.
        persistent_tasks = false,
    )
    @test aqua_persistent_tasks(ChemistryLab)
end

@testsection "Construction tests" begin

    include("species.jl")
    include("formulas.jl")
    include("stoich_matrices.jl")
    include("databases.jl")
    include("databases_remote.jl")
    include("thermofun_methods.jl")
    include("literature.jl")
    include("zeolites.jl")
    include("asr_extension.jl")
    include("cshq_chloride.jl")
    include("parsing_utils.jl")
    include("reactions.jl")
    include("surfaces.jl")
    include("size_classes.jl")

end

@testsection "Speciation tests" begin
    include("speciation.jl")
end

@testsection "System and State tests" begin
    include("chemical_systems.jl")
    include("oxide_budget.jl")
    include("glass_thermochemistry.jl")
    include("chemical_states.jl")
    include("site_families.jl")
    include("test_volume_fractions.jl")
end

@testsection "Thermodynamics tests" begin
    include("thermodynamics.jl")
    include("hkf.jl")
    include("cemdata18_reference.jl")
    include("duan2016_reference.jl")
    include("lothenbach2010_reference.jl")
    include("lothenbach2008_reference.jl")
end

@testsection "Equilibrium tests" begin
    include("activities.jl")
    include("aqueous_properties.jl")
    include("redox.jl")
    include("solid_solutions.jl")
    include("cashplus.jl")
    include("surface_complexation.jl")
    include("ion_exchange.jl")
    include("diffuse_layer.jl")
    include("surface_temperature.jl")
    include("charge_planes.jl")
    include("csh_surface.jl")
    include("sit.jl")
    include("phreeqc_databases.jl")
    include("reaction_formats.jl")
    include("reaktoro_yaml.jl")
    include("provenance.jl")
    include("identifiability.jl")
    include("ignition_loss.jl")
    include("thermogram.jl")
    include("sorption_models.jl")
    include("claysor.jl")
    include("equilibrium_reference.jl")
    include("published_values.jl")
    include("atkins1992_reference.jl")
    include("limestone_blending_reference.jl")
    include("chloride_binding_reference.jl")
    include("chloride_blended_reference.jl")
    include("hong_glasser1999_reference.jl")
    include("test_dual_solver.jl")
    include("certified_equilibrium.jl")
    include("aqueous_inversion.jl")
    include("trace_sensitivity.jl")
    include("equilibrium_constraints.jl")
    include("pressure.jl")
    include("water_eos_reference.jl")
    include("validation_high_temperature.jl")
    include("leaching.jl")
    include("sulfate_attack.jl")
    include("seawater.jl")
    include("seawater_flushing.jl")
    include("delayed_ettringite.jl")
    include("hemicarbonate.jl")
    include("chloride_afm.jl")
    include("real_gases.jl")
    include("capillary.jl")
    include("pitzer.jl")
    include("recipes.jl")
    include("validation_lw2006.jl")
    include("validation_deweerdt2011.jl")
    include("validation_shi2016.jl")
    include("validation_lothenbach2008.jl")
    include("validation_lothenbach2008_temperature.jl")
    include("validation_lothenbach2014.jl")
    include("validation_snellings2022_heat.jl")
    include("validation_chloride_temperature.jl")
    include("validation_logk_temperature.jl")
    include("validation_lhopital2016_alkali.jl")
    include("validation_lhopital2016_aluminum.jl")
    include("validation_wiebe_gaddy1940.jl")
    include("validation_scholer2017.jl")
    include("validation_deschner2012.jl")
    include("validation_deschner2013.jl")
end

@testsection "Utils tests" begin
    include("utils.jl")
    include("scripts.jl")
    include("docs_nomenclature.jl")
end

@testsection "Kinetics tests" begin
    include("kinetics/test_rate_models.jl")
    include("kinetics/test_calorimetry.jl")
    include("kinetics/test_postprocessing.jl")
    include("kinetics/test_feasibility.jl")
    include("kinetics/test_calibration.jl")
    include("kinetics/test_implicit_step.jl")
    include("kinetics/test_pore_humidity.jl")
    include("kinetics/test_surface_coupling.jl")
    include("kinetics/test_sorption_kinetics.jl")
    include("kinetics/test_rate_gates.jl")
    include("kinetics/test_seed_populations.jl")
    include("kinetics/test_covering_layers.jl")
    include("kinetics/test_parameter_ad.jl")
    include("kinetics/test_rhs_speciation.jl")
    include("kinetics/test_compound_energy_kinetics.jl")
    include("kinetics/test_blended_kinetics.jl")
    include("kinetics/test_quaternary_kinetics.jl")
    include("kinetics/test_ternary_kinetics.jl")
    include("kinetics/test_slag_temperature.jl")
    include("kinetics/test_slag_temperature_pastes.jl")
    include("kinetics/test_palandri_kharaka.jl")
    include("kinetics/test_glass_dissolution.jl")
    include("coupling_reference.jl")
end

# Last: loading Ipopt registers one more starting point for the certified search,
# which no test above may depend on.
@testsection "Ipopt extension" begin
    include("ipopt_extension.jl")
end

# Every unit of this file belongs to exactly one part, and every unit a part
# names exists: a file added to a section and to no part would run in no CI job.
# Read from the source, so that every part checks it, not only a full run.
@testset "Test parts" begin
    units = String[]
    for sec in eachmatch(r"^@testsection \"([^\"]+)\" begin(.*?)^end"ms, read(@__FILE__, String))
        files = [m.captures[1] for m in eachmatch(r"^\s*include\(\"([^\"]+)\"\)"m, sec.captures[2])]
        append!(units, isempty(files) ? [sec.captures[1]] : files)
    end
    listed = reduce(vcat, collect(values(TEST_PARTS)))
    @test allunique(listed)
    @test sort(listed) == sort(units)
end

print_timer()
println()
