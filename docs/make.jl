using ChemistryLab
using Optimization, OptimizationIpopt  # load extension OptimizationIpoptExt
using OptimaSolver                     # load extension OptimaSolverExt
using OrdinaryDiffEq                  # load extension KineticsOrdinaryDiffEqExt
using Documenter
using Logging
using DocumenterCitations
# VitePress renders the site from the Markdown that Documenter emits, and
# typesets every formula — chemical equations included — at build time into
# static SVG, so no MathJax bundle is ever fetched by the reader's browser.
using DocumenterVitepress
using PrettyTables

include("pages.jl")

# ── Guard: the precomputed trajectories must not be stale ────────────────────
#
# The heavy coupled runs are computed once by `scripts/precompute_docs.jl` and
# read back by the pages. That is a cache, and a cache that can go stale without
# saying so is not a cache but a false claim: the documentation would keep
# showing the trajectory of a solver that no longer exists.
#
# So each file records the commit it was produced at AND the version of the
# optimizer that produced it, and this compares both against what this build
# resolves. It refuses the build rather than warning, because a warning in a
# three-hundred-line build log is not read.
let
    dir = joinpath(@__DIR__, "src", "assets", "precomputed")
    files = isdir(dir) ? filter(endswith(".csv"), readdir(dir)) : String[]
    if !isempty(files)
        # What the stored results depend on: the solver, the kinetics, and the
        # two scripts that drive them. Documentation and tests are excluded --
        # they cannot change a trajectory.
        watched = ["src", "scripts/ionic_hydration.jl", "scripts/hydration_calibration.jl",
            "scripts/precompute_docs.jl"]
        code_commit = try
            readchomp(`git -C $(dirname(@__DIR__)) log -1 --format=%H -- $watched`)
        catch
            ""
        end
        # The solver is source too, and it is NOT in this repository, so a commit
        # cannot speak for it. 0.5.3 added a stability test to the certificate,
        # which changes how many instants a replay reports as proved; a stored
        # result produced under an earlier one is stale in exactly the same sense.
        solver_version = try
            string(pkgversion(OptimaSolver))
        catch
            ""
        end
        stale = String[]
        for f in files
            stored = ""
            stored_solver = ""
            for line in eachline(joinpath(dir, f))
                startswith(line, "#") || break
                m = match(r"^#\s*commit:\s*(\S+)", line)
                m === nothing || (stored = m.captures[1])
                m = match(r"^#\s*OptimaSolver version:\s*(\S+)", line)
                m === nothing || (stored_solver = m.captures[1])
            end
            if !isempty(stored) && !isempty(code_commit) &&
                    !startswith(code_commit, stored)
                push!(stale, "$f (produced at $stored)")
            elseif !isempty(solver_version) && stored_solver != solver_version
                push!(
                    stale,
                    "$f (produced with OptimaSolver " *
                        (isempty(stored_solver) ? "unrecorded" : stored_solver) *
                        ", this build resolves $solver_version)",
                )
            end
        end
        isempty(stale) || error(
            "precomputed results are older than the code that produces them.\n" *
                "  last commit touching the package or the driving scripts: " *
                "$(code_commit[1:min(end, 8)])\n" *
                "  OptimaSolver resolved by this build: $solver_version\n" *
                join("  stale: " .* stale, "\n") *
                "\n\nRegenerate them with\n" *
                "    julia --project=docs scripts/precompute_docs.jl\n" *
                "or, if the change cannot have moved a trajectory, re-run it anyway: " *
                "a stored result that no longer matches its source is worse than " *
                "a slow build.",
        )
    end
end

# ── Guard: no page may set the plot font ─────────────────────────────────────
#
# Documenter runs every `@example` block in ONE process, so `default(fontfamily
# = ...)` on any page silently applies to every page built after it. When the
# chosen font lacks a glyph a label needs — Computer Modern has no Unicode
# sub/superscripts, and the labels here use them freely — GR prints `GKS: glyph
# missing from current font` for each one and hunts for a fallback. That is
# instant locally, where fontconfig is warm, but not in a CI container, and it
# buries the build log besides.
#
# The font is therefore chosen once, above, as one that carries those glyphs,
# and this refuses any page that tries to change it. Two seconds of scanning
# instead of a buried log and a slow fallback.
let
    offenders = String[]
    for (root, _, files) in walkdir(joinpath(@__DIR__, "src")),
            f in filter(endswith(".md"), files)

        path = joinpath(root, f)
        for (i, line) in enumerate(eachline(path))
            startswith(strip(line), "#") && continue        # a comment, not a call
            occursin(r"\bfontfamily\s*=", line) || continue
            push!(offenders, "  $(relpath(path, @__DIR__)):$i  $(strip(line))")
        end
    end
    isempty(offenders) || error(
        "a page sets `fontfamily`. Documenter shares one process, so this " *
        "applies to every page built after it; a font missing a label's glyph " *
        "then sends GR into a fallback that is slow in CI. Set the font only " *
        "in docs/make.jl; a page may still set `framestyle`, `grid` and the " *
        "like.\n" * join(offenders, "\n"),
    )
end

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"); style = :authoryear)

DocMeta.setdocmeta!(
    ChemistryLab,
    :DocTestSetup,
    :(using ChemistryLab, DynamicQuantities, OrderedCollections, Symbolics);
    recursive=true,
)

ENV["FORCE_COLOR"] = "true"
ENV["COLUMNS"] = "200"
ENV["LINES"] = "100"
ENV["GKSwstype"] = "100"   # headless GR backend — prevents Plots from hanging in doc builds

# Composite panels crop their outer labels unless the margins are explicit — a
# big enough canvas is not enough on its own. Setting the defaults once here
# covers every `@example` block, which run inside this process.
using Plots
Plots.gr()
Plots.default(;
    left_margin = 6Plots.mm,
    bottom_margin = 6Plots.mm,
    right_margin = 4Plots.mm,
    top_margin = 3Plots.mm,
    # The font is set HERE and nowhere else, and that is not a style preference.
    #
    # Documenter runs every `@example` block in one process, so a page calling
    # `default(fontfamily = ...)` changes the font for every page built after it.
    # `examples/cem1_solid_solutions.md` did exactly that with "Computer Modern",
    # which has no Unicode sub/superscripts.
    #
    # Whenever a label needs a glyph the font lacks, GR prints `GKS: glyph
    # missing from current font` and falls back. Locally that fallback is
    # instant, because fontconfig's cache is warm; in a CI container it is not,
    # and a documentation build spent its whole two-hour budget emitting those
    # lines and was canceled by the timeout. So this is not cosmetic.
    #
    # Two defenses, because either alone is fragile: the font is pinned here to
    # one that carries the glyphs, and no plot label in `docs/src` uses Unicode
    # sub/superscripts at all. A page may set `framestyle`, `grid` and the like;
    # it must not set the font.
    fontfamily = "sans-serif",
)

# ── Stopgap: CitationSiteNode in the Markdown writer ─────────────────────────
#
# GUARDED, and the guard is what makes one file serve both versions.
#
# `CitationSiteNode` exists only in DocumenterCitations 1.5, while
# DocumenterVitepress declares a weak dependency `DocumenterCitations = "1 - 1.4"`
# — an upstream cap, so the docs environment resolves to 1.4 and the name is not
# there. Referring to it unconditionally is an `UndefVarError` raised before the
# first page is built, which is exactly how this was found. On 1.4 the citations
# need no help at all: the node is what 1.5 introduced.
#
# When DocumenterVitepress raises its cap, widening the bound in
# `docs/Project.toml` is the only edit needed — this method starts applying again
# on its own.
if isdefined(DocumenterCitations, :CitationSiteNode)
# DocumenterCitations 1.5 wraps every expanded citation in a `CitationSiteNode`,
# whose only purpose is to give the citation an HTML anchor so the bibliography
# can link back to it. Its own docstring calls it "transparent in any output
# format other than HTML", and both the LaTeX writer and MDFlatten implement it
# as "render my children".
#
# DocumenterVitepress 0.3.5 ships a DocumenterCitations extension, but it covers
# only `BibliographyNode`. With no method for `CitationSiteNode`, the writer
# falls through to its generic branch, which prints `Markdown.plain(element)` —
# so every one of this manual's citations came out as the literal text
# `DocumenterCitations.CitationSiteNode("kachanov1992-cite-1")`.
#
# The same one-line treatment as the other non-HTML writers. Remove this once
# DocumenterVitepress covers the node upstream.
    function DocumenterVitepress.render(
            io::IO,
            mime::MIME"text/plain",
            node::Documenter.MarkdownAST.Node,
            ::DocumenterCitations.CitationSiteNode,
            page,
            doc;
            kwargs...,
        )
        return DocumenterVitepress.render(
            io, mime, node, node.children, page, doc; kwargs...,
        )
    end
end

# ── Stopgap: heading anchors that contain LaTeX ──────────────────────────────
# DocumenterVitepress builds each heading as `## <text> {#<slug>}`, where the
# slug is Documenter's anchor label passed through its own
# `sanitized_anchor_label` — whose comment says "vitepress doesn't like special
# markdown characters in the id slug", but which only strips `[ ] ( ) *`.
#
# A heading such as `## The COD tensor ``\boldsymbol{B}``` yields the slug
# `The-COD-tensor-\boldsymbol{B}`. VitePress's `{#...}` parser rejects the
# backslash and the braces, so it treats the whole suffix as *text*: the heading
# renders as "The COD tensor {#The-COD-tensor-\boldsymbol{B}}", the formula is
# dropped, and the same garbage lands in the "On this page" outline. Twenty-three
# headings across seven pages were affected.
#
# Stripping those characters from the slug is safe here: nothing links to those
# anchors (checked against every `](…#…)` destination in the generated Markdown),
# and this narrows to headings only, leaving docstring anchors — which legitimately
# carry braces, are emitted as raw `<a id=…>`, and *are* linked to — untouched.
#
# Remove once `sanitized_anchor_label` covers these characters upstream.
function DocumenterVitepress.render(
        io::IO,
        mime::MIME"text/plain",
        node::Documenter.MarkdownAST.Node,
        header::Documenter.AnchoredHeader,
        page,
        doc;
        kwargs...,
    )
    anchor = header.anchor
    label = DocumenterVitepress.sanitized_anchor_label(anchor)
    id = replace(replace(label, r"[\\{}]" => ""), " " => "-")
    heading = first(node.children)
    println(io)
    print(io, "#"^(heading.element.level), " ")
    heading_iob = IOBuffer()
    DocumenterVitepress.render(heading_iob, mime, node, heading.children, page, doc; kwargs...)
    print(io, rstrip(String(take!(heading_iob))))
    print(io, " {#$(id)}")
    if haskey(kwargs, :inventory)
        item = DocumenterVitepress.InventoryItem(
            name = anchor.id,
            domain = "std",
            role = "label",
            dispname = DocumenterVitepress._get_inventory_dispname(
                anchor.id, Documenter.MDFlatten.mdflatten(anchor.node)
            ),
            priority = -1,
            uri = DocumenterVitepress._get_inventory_uri(doc, page, id),
        )
        push!(kwargs[:inventory], item)
    end
    println(io)
    return nothing
end

# ── Stopgap: ordered lists start at 2, and swallow their first item ──────────
# DocumenterVitepress numbers ordered-list items with `bullet(i) = "$(i+1). "`,
# but `enumerate` is already 1-based: every ordered list in the manual came out
# numbered from 2. It also emits no blank line before the list.
#
# Together those two produce the damage seen on the References page. A list whose
# first marker is `2.` cannot interrupt a paragraph — CommonMark allows that only
# for a list starting at `1.` — so, with no blank line to separate them, the first
# entry was absorbed into the preceding prose as plain text and the list began at
# `3.`. Eighteen ordered lists across the manual were affected; not one of them
# started at 1.
#
# Remove once the numbering is fixed upstream.
function DocumenterVitepress.render(
        io::IO,
        mime::MIME"text/plain",
        node::Documenter.MarkdownAST.Node,
        list::Documenter.MarkdownAST.List,
        page,
        doc;
        kwargs...,
    )
    bullet(i) = list.type === :ordered ? "$(i). " : "- "
    println(io)
    iob = IOBuffer()
    for (i, item) in enumerate(node.children)
        DocumenterVitepress.render(
            iob, mime, item, item.children, page, doc; prenewline = false, kwargs...
        )
        eachline = split(String(take!(iob)), '\n')
        # Continuation lines must line up with the text, i.e. under the marker's
        # full width. Upstream hard-codes two spaces, which fits `- ` but not
        # `1. `: a display equation inside an ordered item fell out of the list,
        # splitting it in two and restarting the numbering.
        pad = " "^length(bullet(i))
        eachline[2:end] .= pad .* eachline[2:end]
        final_string = join(eachline, '\n')
        endswith(final_string, '\n') || (final_string *= "\n")
        print(io, bullet(i))
        print(io, final_string)
    end
    return nothing
end

# ── Per-block timing, so a slow build says what is slow ──────────────────────
#
# Under `JULIA_DEBUG=Documenter` Documenter announces each block it is about to
# evaluate, but without timing — so a three-hour build names three hundred
# blocks and does not say which one spent the three hours. Diagnosing that from
# outside is guesswork, and guesswork on this has already cost several rounds.
#
# This logger passes every message through untouched and, each time a new block
# starts, reports on stderr how long the PREVIOUS one took. Blocks under the
# threshold stay silent, so the log gains a line only where there is something
# to see. `println` rather than `@info`, deliberately: emitting a log record
# from inside a log handler re-enters the handler.
const SLOW_BLOCK_SECONDS = 5.0

struct BlockTimer{L <: AbstractLogger} <: AbstractLogger
    inner::L
    t0::Base.RefValue{Float64}
    label::Base.RefValue{String}
    total::Base.RefValue{Float64}
end

Logging.min_enabled_level(l::BlockTimer) = Logging.min_enabled_level(l.inner)
Logging.shouldlog(l::BlockTimer, args...) = Logging.shouldlog(l.inner, args...)
Logging.catch_exceptions(l::BlockTimer) = Logging.catch_exceptions(l.inner)

function Logging.handle_message(
        l::BlockTimer, level, message, _module, group, id, file, line; kwargs...,
    )
    msg = string(message)
    if occursin("Evaluating ", msg) && occursin("block:", msg)
        now = time()
        dt = now - l.t0[]
        l.total[] += dt
        if dt >= SLOW_BLOCK_SECONDS
            println(
                stderr,
                "⏱  previous block took ", round(dt; digits = 1), " s",
                "  (running total ", round(l.total[] / 60; digits = 1), " min)",
                "  — ", l.label[],
            )
            flush(stderr)
        end
        l.t0[] = now
        body = replace(msg, r"^.*?block:\s*"s => "")
        l.label[] = first(split(strip(body), '\n'))
    end
    return Logging.handle_message(
        l.inner, level, message, _module, group, id, file, line; kwargs...,
    )
end

Logging.with_logger(
    BlockTimer(Logging.current_logger(), Ref(time()), Ref("start"), Ref(0.0)),
) do

makedocs(;
    # `clean = false` lets pages deleted from the source survive in `build/`
    # and go on being deployed. Nothing writes there before `makedocs`.
    modules=[ChemistryLab],
    remotes=nothing,
    authors="Jean-François Barthélémy and Anthony Soive",
    sitename="ChemistryLab.jl",
    # The favicon and the logo are picked up automatically from `docs/src/assets`,
    # and the sidebar is derived from `pages`, so neither needs declaring here.
    # mhchem now loads in `docs/src/.vitepress/mathjax-plugin.ts` instead of in a
    # `mathengine`, because the typesetting happens at build time.
    format=DocumenterVitepress.MarkdownVitepress(;
        repo = "https://github.com/MicroPoroChemoMechanics/ChemistryLab.jl",
        devbranch = "main",
        devurl = "dev",
        deploy_url = "https://MicroPoroChemoMechanics.github.io/ChemistryLab.jl",
        description = "Aqueous and cement chemistry in Julia: speciation, equilibrium and kinetics",
    ),
    pages=pages,
    plugins=[bib],
    warnonly=[:docs_block],
    draft=false,
)

end  # Logging.with_logger

# DocumenterVitepress writes a real directory per version rather than the
# symlinks Documenter used, so it needs its own `deploydocs`.
DocumenterVitepress.deploydocs(;
    repo         = "github.com/MicroPoroChemoMechanics/ChemistryLab.jl.git",
    target       = joinpath(@__DIR__, "build"),
    branch       = "gh-pages",
    devbranch    = "main",
    push_preview = false,
)
