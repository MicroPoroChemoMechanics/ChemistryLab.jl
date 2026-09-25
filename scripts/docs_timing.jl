# =============================================================================
#  docs_timing.jl — where the documentation build spends its time
#
#  The executed blocks of the documentation pages are run here the way
#  Documenter runs them: one module per block name (a fresh one for each unnamed
#  block), pages in the order of `docs/pages.jl`, and every page in one process,
#  so that compilation is paid once as it is in the real build. The time of each
#  block is appended to a CSV as soon as its page finishes, so that a run which
#  is interrupted keeps what it had measured.
#
#    julia --project=docs scripts/docs_timing.jl                  # every page
#    julia --project=docs scripts/docs_timing.jl cem4 theory/     # pages whose
#                                                                 # path matches
#
#  The CSV (`page,block,name,seconds,error`) goes to `docs_timing.csv` in the
#  current directory, or to the path in `CHEMLAB_TIMING_CSV`, and a summary
#  sorted by cost is printed at the end.
#
#  What is not reproduced: rendering, the VitePress build and the doctests.
#  What is measured is the executed blocks, where nearly all of the build goes.
#  It is one heavy process, like the build itself, and is not to be run beside
#  another.
# =============================================================================

using Printf

# The packages `docs/make.jl` loads, in the same order, BEFORE any page runs.
# Which extensions are active, and in which order they registered, decides the
# route a certified solve takes: measured, a map that let the pages load Ipopt
# themselves, late, put a different back end first and timed a different build.
using ChemistryLab
using Optimization, OptimizationIpopt
using OptimaSolver
using OrdinaryDiffEq

ENV["GKSwstype"] = "100"          # headless GR, as in docs/make.jl
using Plots
Plots.gr()
Plots.default(; fontfamily = "sans-serif")   # the font docs/make.jl pins

const ROOT = dirname(@__DIR__)
const SRC = joinpath(ROOT, "docs", "src")

include(joinpath(ROOT, "docs", "pages.jl"))

leaves(node::AbstractString) = [node]
leaves(node::Pair) = leaves(node.second)
leaves(node::AbstractVector) = reduce(vcat, leaves.(node); init = String[])

"""
    executed_blocks(path) -> Vector{Tuple{Union{String,Nothing},String}}

The `@example` and `@setup` blocks of a page, in order, with their names
(`nothing` for an unnamed block).
"""
function executed_blocks(path)
    text = read(path, String)
    return [
        (
            m.captures[2] === nothing || isempty(strip(m.captures[2])) ? nothing :
                String(strip(m.captures[2])), String(m.captures[3]),
        )
            for m in eachmatch(r"^```@(example|setup)(?:[ \t]+([^\n]*))?\n(.*?)^```"ms, text)
    ]
end

# A module as Documenter builds one for a named block: a bare `Module` has no
# `include` or `eval` of its own, and a page that includes a script needs both.
function sandbox(name)
    mod = Module(name)
    Core.eval(mod, :(include(x) = Base.include($mod, abspath(x))))
    Core.eval(mod, :(eval(x) = Core.eval($mod, x)))
    return mod
end

function time_page(rel, io)
    blocks = executed_blocks(joinpath(SRC, rel))
    isempty(blocks) && return 0.0
    modules = Dict{String, Module}()
    total = 0.0
    dir = mktempdir()
    for (k, (name, code)) in enumerate(blocks)
        mod = name === nothing ? sandbox(gensym("ex")) :
            get!(() -> sandbox(Symbol("ex_", name)), modules, name)
        err = ""
        t = @elapsed try
            cd(dir) do
                redirect_stdout(devnull) do
                    Base.CoreLogging.with_logger(Base.CoreLogging.NullLogger()) do
                        Base.include_string(mod, code, "$rel:block$k")
                    end
                end
            end
        catch e
            err = replace(first(split(sprint(showerror, e), '\n')), "," => ";")
        end
        total += t
        @printf(io, "%s,%d,%s,%.3f,%s\n", rel, k, something(name, ""), t, err)
    end
    flush(io)
    rm(dir; recursive = true, force = true)
    return total
end

function main(patterns = ARGS)
    selected = [p for p in leaves(pages) if isempty(patterns) || any(occursin(q, p) for q in patterns)]
    csv = get(ENV, "CHEMLAB_TIMING_CSV", "docs_timing.csv")
    totals = Pair{String, Float64}[]
    open(csv, "w") do io
        println(io, "page,block,name,seconds,error")
        for rel in selected
            t = time_page(rel, io)
            t > 0 && push!(totals, rel => t)
            t > 0 && (@printf("%9.1f s  %s\n", t, rel); flush(stdout))
        end
    end
    println("\nsorted by cost:")
    for (rel, t) in sort(totals; by = last, rev = true)
        @printf("%9.1f s  %s\n", t, rel)
    end
    @printf("%9.1f s  TOTAL over %d pages\n", sum(last, totals; init = 0.0), length(totals))
    println("per-block detail in ", csv)
    return
end

abspath(PROGRAM_FILE) == (@__FILE__) && main()
