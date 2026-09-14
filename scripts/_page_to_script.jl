# =============================================================================
#  _page_to_script.jl — turn a documentation page into a runnable script
#
#  The binder-family pages carry their calculation in a dozen small blocks
#  separated by prose, each showing an intermediate result. That is good for a
#  reader and useless for someone who wants to RUN the calculation, so the blocks
#  are extracted here into one file per page.
#
#  Every block of those pages shares a single name, which is what makes this
#  sound: Documenter gives one module per block name, so the blocks of a page
#  already run as one program in one scope, and concatenating them changes
#  nothing. `_page_blocks` refuses a page where that is not true rather than
#  emitting a script that would run differently from the page.
#
#  The page stays the source. The scripts are generated, and `test/scripts.jl`
#  fails if any of them has drifted from its page — so they cannot rot. That test
#  compares PARSED EXPRESSIONS, not bytes, so a generated file can be formatted
#  like the rest of the repository without the comparison becoming a formatting
#  check.
#
#  Regenerate, then format — in that order, because regenerating overwrites the
#  formatting:
#    julia --project=. scripts/_page_to_script.jl
#    julia --project=@runic -e 'using Runic; Runic.main(["--inplace", "scripts"])'
#
#  Runic is pinned to 1.10 by .github/workflows/Format.yml; use the same version
#  or the two will disagree about what "formatted" means.
# =============================================================================

const PAGES = [
    "examples/cem1_from_clinker.md",
    "examples/cem2_blended.md",
    "examples/cem3_slag.md",
    "examples/cem4_pozzolanic.md",
    "examples/cem5_composite.md",
    "examples/cement_wc_ratio.md",
    "examples/miscibility_gap.md",
]

"""
    _page_blocks(path) -> (name, Vector{String})

The executed blocks of a documentation page, in order, with the name they share.

Throws when the page mixes block names: those run in different modules and
concatenating them would produce a program the page never ran.
"""
function _page_blocks(path::AbstractString)
    text = read(path, String)
    names, codes = String[], String[]
    for m in eachmatch(r"^```@(example|setup)(?:[ \t]+([^\n]*))?\n(.*?)^```"ms, text)
        push!(names, m.captures[2] === nothing ? "" : strip(m.captures[2]))
        push!(codes, m.captures[3])
    end
    uniq = unique(names)
    length(uniq) == 1 && !isempty(only(uniq)) || error(
        "$(basename(path)) mixes block names $(uniq): its blocks run in " *
            "different modules and cannot be concatenated into one script."
    )
    return only(uniq), codes
end

"""
    _script_text(page, relpath) -> String

The generated script for one page, header included.
"""
function _script_text(path::AbstractString, rel::AbstractString)
    _, codes = _page_blocks(path)
    io = IOBuffer()
    println(io, "# " * "="^77)
    println(io, "#  GENERATED FILE — do not edit.")
    println(io, "#")
    println(io, "#  Extracted from docs/src/$rel, whose executed blocks it reproduces in")
    println(io, "#  order. Edit the page, then regenerate:")
    println(io, "#")
    println(io, "#      julia --project=. scripts/_page_to_script.jl")
    println(io, "#")
    println(io, "#  `test/scripts.jl` fails if this file and its page disagree.")
    println(io, "#")
    println(io, "#  Usage:")
    println(io, "#      julia --project=docs scripts/$(replace(basename(rel), ".md" => ".jl"))")
    println(io, "#")
    println(io, "#  No `Pkg.activate` here, deliberately: the active project is global process")
    println(io, "#  state and this file is meant to be `include`d as well as run.")
    println(io, "#")
    println(io, "#  Formatted with Runic like the rest of the repository. `test/scripts.jl`")
    println(io, "#  compares SYNTAX TREES rather than text, so formatting is free to differ")
    println(io, "#  from the page while a real divergence still fails the suite.")
    println(io, "# " * "="^77)
    for code in codes
        println(io)
        print(io, code)
    end
    return String(take!(io))
end

function main(root = dirname(@__DIR__))
    for rel in PAGES
        path = joinpath(root, "docs", "src", rel)
        out = joinpath(root, "scripts", replace(basename(rel), ".md" => ".jl"))
        write(out, _script_text(path, rel))
        println("  wrote ", relpath(out, root))
    end
    return
end

# `@__FILE__` needs the parentheses: bare, the parser reads `&& main()` as
# more macro arguments.
abspath(PROGRAM_FILE) == (@__FILE__) && main()
