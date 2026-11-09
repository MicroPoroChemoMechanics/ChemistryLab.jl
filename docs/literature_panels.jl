# The published data a page reads, shown where it reads them.
#
# Every value a page takes from `data/literature/` goes through
# `ChemistryLab.literature(key)`, which keeps what it read in a cache. Each
# executed block of the documentation starts with that cache empty, so that
# what is in it once the block has run is exactly what the block read; under the
# block, each source read for the first time on the page is then shown, folded:
# a card for the source (citation, DOI, where in it the values are read, how
# they were transcribed and checked), its quantities and its tables, with their
# units, the kind of each value and, for a table read off a figure, how it was
# digitized. No page has to ask for it, so no page can forget it.
#
# The hook is a Documenter expander placed just before the ones that run
# `@eval`, `@example`, `@repl` and `@setup` blocks (Documenter's `Selectors` extension
# point, the one DocumenterCitations uses for its bibliography blocks): it
# empties the cache, hands the block to Documenter's own expander, and inserts
# the panels after it. Included by `make.jl`, after `source_links.jl`, whose
# `source_url` it reuses.

using JSON
using Documenter
import Documenter: Selectors
const _Expanders = Documenter.Expanders
const _MarkdownAST = Documenter.MarkdownAST

# ── HTML pieces ──────────────────────────────────────────────────────────────

"The HTML text of `x`, with the five characters that HTML reads escaped."
_esc(x) = replace(string(x), '&' => "&amp;", '<' => "&lt;", '>' => "&gt;", '"' => "&quot;", '\'' => "&#39;")
_esc(::Nothing) = ""

"A table cell: a missing value is a dash, a number is set in figures, anything else is text."
_cell(::Nothing) = """<td class="na">—</td>"""
_cell(x::Real) = """<td class="num">$(_esc(x))</td>"""
_cell(x) = """<td>$(_esc(x))</td>"""

"An HTML table of `header` (already rendered) over `rows` (vectors of rendered cells)."
function _html_table(header, rows; class = "")
    head = join(("<th>$h</th>" for h in header))
    body = join(("<tr>" * join(r) * "</tr>" for r in rows), "\n")
    return """<div class="data-table-wrap"><table class="data-table $class"><thead><tr>$head</tr></thead><tbody>
    $body
    </tbody></table></div>"""
end

"A unit as it reads after a number or in a column head: a pure number has none."
_unit(u) = u === nothing || u == "1" ? "" : """ <span class="unit">$(_esc(u))</span>"""
_head(c, u) = "$(_esc(c))" * (u === nothing || u == "1" ? "" : """ <span class="unit">[$(_esc(u))]</span>""")

# Every link of a panel opens a new tab, as VitePress does for the external links
# of a page. It must: the sites of the MPCM packages share one origin, and
# VitePress's router takes a link of the same origin without a `target` for a
# page of the current site, which it cannot find (a 404 until the page is
# reloaded).
"A link to the DOI `d`, or nothing when there is none."
_doi_link(::Nothing) = ""
_doi_link(d) = """ <a target="_blank" rel="noreferrer" href="https://doi.org/$(_esc(d))">doi:$(_esc(d))</a>"""

const _KIND_WORDS = Dict(
    "measured" => "measured", "fitted" => "fitted", "estimated" => "estimated",
    "placeholder" => "placeholder", "unstated" => "provenance unstated",
)
# Most values are published as they are; a badge marks the others only.
_kind(k) = k === nothing || k == "published" ? "" :
    """<span class="badge kind-$(_esc(k))">$(_esc(get(_KIND_WORDS, k, k)))</span>"""

"An uncertainty after a value: `± u` for a number, the words of a relative or unit one."
_uncertainty(::Nothing, unit) = ""
_uncertainty(u::Real, unit) = """ <span class="unc">± $(_esc(u))</span>"""
function _uncertainty(u::AbstractDict, unit)
    haskey(u, "relative") && return """ <span class="unc">± $(_esc(100 * u["relative"])) %</span>"""
    haskey(u, "value") && return """ <span class="unc">± $(_esc(u["value"]))$(_unit(get(u, "unit", unit)))</span>"""
    return ""
end
_uncertainty(u, unit) = """ <span class="unc">± $(_esc(u))</span>"""

# ── Sparklines ───────────────────────────────────────────────────────────────

const _SPARK_MAX_POINTS = 300
const _SPARK_MAX_SERIES = 4

"""
    _sparkline(rows) -> String

An inline SVG of a numerical table: its first column as abscissa, up to four of
the next numerical columns as curves, each scaled to the box. An abscissa spread
over decades (a time, a concentration) is read on a logarithmic scale. Empty
when the table has no numerical abscissa or no numerical column after it.
"""
function _sparkline(rows::AbstractVector; w = 160, h = 40, pad = 3)
    ## The cells are JSON values: `Real` here means an Int or a Float64 of the file.
    length(rows) < 2 && return ""
    all(r -> !isempty(r) && r[1] isa Real, rows) || return ""
    ncol = minimum(length, rows)
    ys = [j for j in 2:ncol if all(r -> r[j] isa Real || r[j] === nothing, rows) && count(r -> r[j] isa Real, rows) ≥ 2]
    isempty(ys) && return ""
    ys = ys[1:min(end, _SPARK_MAX_SERIES)]
    step = max(1, cld(length(rows), _SPARK_MAX_POINTS))
    sample = rows[1:step:end]
    x = [float(r[1]) for r in sample]
    xpos = filter(>(0), x)
    if all(≥(0), x) && !isempty(xpos) && maximum(xpos) > 100 * minimum(xpos)
        x₀ = minimum(xpos)
        x = log1p.(x ./ x₀)
    end
    x0, x1 = extrema(x)
    sx = x1 > x0 ? (w - 2pad) / (x1 - x0) : 0.0
    lines = map(enumerate(ys)) do (k, j)
        pts = [(xi, float(r[j])) for (xi, r) in zip(x, sample) if r[j] isa Real]
        y0, y1 = extrema(last.(pts))
        sy = y1 > y0 ? (h - 2pad) / (y1 - y0) : 0.0
        path = join(("$(round(pad + sx * (a - x0); digits = 1)),$(round(h - pad - sy * (b - y0); digits = 1))" for (a, b) in pts), " ")
        """<polyline class="s$k" points="$path"/>"""
    end
    return """<svg class="sparkline" viewBox="0 0 $w $h" width="$w" height="$h" aria-hidden="true">$(join(lines))</svg>"""
end

# ── One source ───────────────────────────────────────────────────────────────

const _MAX_ROWS_SHOWN = 200

"The error a reading off a figure adds, in words."
_reading(::Nothing) = "reading error not quantified"
_reading(u) = "reading error" * _uncertainty(u, nothing)

"The digitization of a table read off a figure, as a short note."
_digitization(::Nothing) = ""
function _digitization(d)
    reading = _reading(get(d, "uncertainty", nothing))
    return """<p class="digitization"><span class="badge digit">read off a figure</span>
    $(_esc(get(d, "method", ""))), $(_esc(get(d, "tool", ""))); $reading.
    $(haskey(d, "note") ? "<span class=\"note\">$(_esc(d["note"]))</span>" : "")</p>"""
end

"A table of a source, folded, with its units, kind, location and sparkline."
function _table(name, t)
    cols, units, rows = t["columns"], get(t, "units", fill(nothing, length(t["columns"]))), t["rows"]
    shown = rows[1:min(end, _MAX_ROWS_SHOWN)]
    more = length(rows) > length(shown) ?
        """<p class="more">The first $(length(shown)) of $(length(rows)) rows; the file holds them all.</p>""" : ""
    tbl = _html_table([_head(c, u) for (c, u) in zip(cols, units)], [[_cell(x) for x in r] for r in shown]; class = "series")
    return """<details class="data-series"><summary><code>$(_esc(name))</code>: $(length(rows)) rows $(_kind(get(t, "kind", nothing))) $(_sparkline(rows))</summary>
    <p class="where">$(_esc(get(t, "location", "")))</p>
    $(_digitization(get(t, "digitization", nothing)))
    $tbl
    $more
    </details>"""
end

"""
    literature_panel(key) -> String

The file `data/literature/<key>.json` as a folded HTML panel: a card for its
source, the table of its quantities, its tables folded one by one, its notes,
and a link to the file at the commit the site is built from.
"""
function literature_panel(key::AbstractString)
    rel = "data/literature/$key.json"
    file = joinpath(REPO_ROOT, rel)
    isfile(file) || throw(ArgumentError("literature_panel: no file $file"))
    d = JSON.parsefile(file)
    get(d, "schema", nothing) == "chemistrylab-literature/1" ||
        throw(ArgumentError("literature_panel: $rel is not of schema chemistrylab-literature/1"))
    src, tr = d["source"], d["transcription"]
    checked = get(tr, "checked_against_source", false) ?
        """<span class="badge ok">checked against the source</span>""" :
        """<span class="badge warn">not checked against the source</span>"""
    card = """<div class="data-source">
    <p class="citation">$(_esc(src["citation"]))$(_doi_link(get(src, "doi", nothing)))</p>
    <p><span class="label">Read in</span> $(_esc(get(src, "location", ""))) $checked</p>
    <p><span class="label">Transcribed from</span> $(_esc(get(tr, "from", "")))</p>
    $(get(tr, "note", nothing) === nothing ? "" : "<p class=\"note\">$(_esc(tr["note"]))</p>")
    </div>"""
    qs = get(d, "quantities", Dict())
    quantities = isempty(qs) ? "" : "<h4>Quantities</h4>\n" * _html_table(
            ["name", "value", "description", "where"],
            [
                [
                    "<td><code>$(_esc(k))</code></td>",
                    """<td class="num">$(_esc(q["value"]))$(_unit(get(q, "unit", nothing)))$(_uncertainty(get(q, "uncertainty", nothing), get(q, "unit", nothing)))</td>""",
                    "<td>$(_esc(get(q, "description", ""))) $(_kind(get(q, "kind", nothing)))</td>",
                    """<td class="where">$(_esc(get(q, "location", "")))</td>""",
                ] for (k, q) in qs
            ]
        )
    ts = get(d, "tables", Dict())
    tables = isempty(ts) ? "" : "<h4>Tables</h4>\n" * join((_table(k, t) for (k, t) in ts), "\n")
    ns = get(d, "notes", [])
    notes = isempty(ns) ? "" : "<h4>Notes</h4>\n<ul>" * join(("<li>$(_esc(n))</li>" for n in ns)) * "</ul>"
    link = """<p class="data-link">The file <a target="_blank" rel="noreferrer" href="$(source_url(rel))"><code>$rel</code></a>
    on GitHub, at the commit this page was built from.</p>"""
    count = join(
        filter(
            !isempty, [
                isempty(qs) ? "" : "$(length(qs)) quantities",
                isempty(ts) ? "" : "$(length(ts)) tables",
            ]
        ), ", "
    )
    return """<details class="data-file"><summary><span class="data-icon" aria-hidden="true"></span>Data read above: <code>$(_esc(key))</code>
    <span class="data-cite">$(_esc(_short_citation(src["citation"])))</span>
    <span class="data-count">$count</span></summary><div class="data-file-body">
    $card
    $quantities
    $tables
    $notes
    $link
    </div></details>"""
end

"""
Authors and year of a citation, for the bar of a folded panel: what precedes the
title, the first author followed by *et al.* beyond three.
"""
function _short_citation(c::AbstractString)
    m = match(r"^(.*?)\s*(\(\d{4}[a-z]?\))", c)
    m === nothing && return ""
    authors = split(m.captures[1], r"(?<=\.),\s+(?:& )?")
    names = length(authors) > 3 ? first(authors) * " et al." : m.captures[1]
    return "$names $(m.captures[2])"
end

# ── The hook into Documenter ─────────────────────────────────────────────────

"The sources already shown on each page, so that each is shown once, where first read."
const _SHOWN = Dict{String, Set{String}}()

"""
    run_and_show_literature(expander, node, page, doc)

Run `node` with Documenter's `expander`, the literature cache emptied before,
then insert after it a panel for each source it read that the page has not yet
shown.
"""
function run_and_show_literature(expander, node, page, doc)
    lock(() -> empty!(ChemistryLab._LITERATURE_CACHE), ChemistryLab._LITERATURE_LOCK)
    Selectors.runner(expander, node, page, doc)
    read = lock(() -> sort!(collect(keys(ChemistryLab._LITERATURE_CACHE))), ChemistryLab._LITERATURE_LOCK)
    shown = get!(Set{String}, _SHOWN, page.source)
    last = node
    for key in read
        key in shown && continue
        push!(shown, key)
        panel = _MarkdownAST.Node(Documenter.RawNode(:html, literature_panel(key)))
        _MarkdownAST.insert_after!(last, panel)
        last = panel
    end
    return nothing
end

# One expander per kind of executed block, each ordered just before
# Documenter's own, whose matcher it reuses. They are nested pipelines, so that a
# block inside an admonition or a list item is handled too.
for (name, expander, order) in (
        (:LiteratureEvalBlocks, _Expanders.EvalBlocks, 4.9),
        (:LiteratureExampleBlocks, _Expanders.ExampleBlocks, 7.9),
        (:LiteratureREPLBlocks, _Expanders.REPLBlocks, 8.9),
        (:LiteratureSetupBlocks, _Expanders.SetupBlocks, 9.9),
    )
    @eval begin
        abstract type $name <: _Expanders.NestedExpanderPipeline end
        Selectors.order(::Type{$name}) = $order
        Selectors.matcher(::Type{$name}, node, page, doc) = Selectors.matcher($expander, node, page, doc)
        Selectors.runner(::Type{$name}, node, page, doc) = run_and_show_literature($expander, node, page, doc)
    end
end
