# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2026 Jean-François Barthélémy (Cerema, UMR MCD)

# Running the compute sections of a script: in a fresh module, so that one run
# leaves nothing behind for the next, with what it prints kept for the page.

"""
    run_compute(source) -> Dict{String, Any}

Run `source` (sections 1 to 5 of a studio script) in a fresh module, and return
whether it succeeded, its error, what it printed, the time it took and the
results read from its `study`. An interrupt (the page's Stop) propagates, so
that the caller knows the run was stopped rather than failed.
"""
function run_compute(source::AbstractString)
    m = Module(:StudioRun)
    failure = Ref("")
    log = tempname()
    t0 = time()
    open(log, "w") do io
        redirect_stdio(; stdout = io, stderr = io) do
            try
                Base.include_string(m, source, "studio_script.jl")
            catch err
                # `include_string` wraps the interrupt of a Stop in a LoadError
                (err isa LoadError ? err.error : err) isa InterruptException && rethrow()
                # the script's failure is the answer the page shows, with its line
                failure[] = sprint(showerror, err)
            end
        end
    end
    elapsed = time() - t0
    output = read(log, String)
    rm(log; force = true)
    ok = isempty(failure[])
    # the script's globals belong to a later world than this function's
    defined = ok && Base.invokelatest(isdefined, m, :study)
    results = defined ? Base.invokelatest(extract, Base.invokelatest(getfield, m, :study)) : nothing
    return Dict{String, Any}(
        "ok" => ok, "error" => failure[], "output" => output, "elapsed" => elapsed, "results" => results,
    )
end
