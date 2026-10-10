# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

# The plot defaults of the documentation, set once before any page runs, by
# `docs/make.jl` and by `scripts/docs_timing.jl`. Documenter runs every
# `@example` block in one process: a default a page sets applies to every page
# built after it, so a figure would look different in the full build and in a
# shard that does not build that page first. Set here, every page is drawn the
# same way, whatever the order. A page may still set `framestyle`, `grid` and
# the like for its own figures; it must not set the font.
#
# Composite panels crop their outer labels unless the margins are explicit — a
# big enough canvas is not enough on its own.
Plots.default(;
    left_margin = 6Plots.mm,
    bottom_margin = 6Plots.mm,
    right_margin = 4Plots.mm,
    top_margin = 3Plots.mm,
    framestyle = :box,
    grid = false,
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
