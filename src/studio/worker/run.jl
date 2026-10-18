# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2026 Jean-François Barthélémy (Cerema, UMR MCD)

# Running a piece of a script: in a fresh module, so that one run leaves nothing
# behind for the next, with what it prints and what it warns kept for the page.

"""
    run_in_module(source, name) -> NamedTuple

Run `source` in a fresh module, its output and its log messages written to one
text: `(; ok, error, output, elapsed, mod)`. A failure of the code is reported,
with its line; an interrupt (the page's Stop) propagates, so that the caller
knows the run was stopped rather than failed.
"""
function run_in_module(source::AbstractString, name::Symbol)
    mod = Module(name)
    failure = Ref("")
    log = tempname()
    t0 = time()
    open(log, "w") do io
        redirect_stdio(; stdout = io, stderr = io) do
            with_logger(ConsoleLogger(io)) do
                try
                    Base.include_string(mod, source, "studio_script.jl")
                catch err
                    # `include_string` wraps the interrupt of a Stop in a LoadError
                    (err isa LoadError ? err.error : err) isa InterruptException && rethrow()
                    # the script's failure is the answer the page shows, with its line
                    failure[] = sprint(showerror, err)
                end
            end
        end
    end
    elapsed = time() - t0
    output = read(log, String)
    rm(log; force = true)
    return (; ok = isempty(failure[]), error = failure[], output, elapsed, mod)
end

# A global of the module, read in the world the script ran in.
global_of(mod::Module, name::Symbol) =
    Base.invokelatest(isdefined, mod, name) ? Base.invokelatest(getfield, mod, name) : nothing

"""
    run_compute(source) -> Dict{String, Any}

Run sections 1 to 5 of a studio script, and return whether they succeeded, the
error, what they printed, the time they took and the results read from
`study`.
"""
function run_compute(source::AbstractString)
    r = run_in_module(source, :StudioRun)
    study = r.ok ? global_of(r.mod, :study) : nothing
    results = study === nothing ? nothing : Base.invokelatest(extract, study)
    return Dict{String, Any}(
        "ok" => r.ok, "error" => r.error, "output" => r.output, "elapsed" => r.elapsed, "results" => results,
    )
end

"""
    run_diagnostics(source) -> Dict{String, Any}

Run the setup of a studio script (database, system, materials, model) ending
with `diagnostics = diagnose_system(…)`, and return the findings, with the
warnings ChemistryLab printed on the way.
"""
function run_diagnostics(source::AbstractString)
    r = run_in_module(source, :StudioDiagnostics)
    found = r.ok ? global_of(r.mod, :diagnostics) : nothing
    return Dict{String, Any}(
        "ok" => r.ok, "error" => r.error, "output" => r.output,
        "diagnostics" => found === nothing ? Any[] : [
                Dict{String, Any}(
                    "level" => String(d.level), "code" => String(d.code), "message" => d.message,
                    "fix" => d.fix, "subjects" => d.subjects,
                ) for d in found
            ],
    )
end
