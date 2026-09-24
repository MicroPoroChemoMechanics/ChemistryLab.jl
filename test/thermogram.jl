# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

using ChemistryLab
using ChemistryLab: value
using DynamicQuantities
using LinearAlgebra
using Test

include("reference_species.jl")

@testsection "a thermogram, and the windows it takes to have one" begin

    c18 = Dict(
        symbol(s) => s for s in
            build_species(datapath("cemdata18-thermofun.json"); verbose = false)
    )
    aq = Dict(
        symbol(s) => s for s in
            build_species(datapath("slop98-inorganic-thermofun.json"); verbose = false)
    )
    solids = [c18[k] for k in ("Portlandite", "Gp", "Cal")]
    aqueous = [aq[k] for k in ("H2O@", "H+", "Ca+2", "SO4-2", "CO3-2")]
    cs = ChemicalSystem(
        vcat(aqueous, solids),
        [aq["H2O@"], aq["H+"], aq["Ca+2"], aq["SO4-2"], aq["CO3-2"]],
    )
    idx = Dict(symbol(sp) => k for (k, sp) in enumerate(cs.species))
    n = Any[fill(1.0e-12u"mol", length(cs.species))...]
    n[idx["H2O@"]] = moles_of_water() * u"mol"
    n[idx["Portlandite"]] = 1.0u"mol"
    n[idx["Gp"]] = 2.0u"mol"
    n[idx["Cal"]] = 3.0u"mol"
    st = ChemicalState(cs, n)

    # DECLARED PLACEHOLDERS, and they print as such. These are not literature
    # values and nothing here pretends they are: they exist to exercise the
    # forward path and to be the starting point an identification moves away
    # from. The real ones come from a publication or from a measured curve.
    placeholder(phase, T, w; releases = :water) = DecompositionWindow(
        phase, T, w;
        releases, kind = PROV_PLACEHOLDER, source = "illustrative, not measured",
    )
    windows = [
        placeholder("Gp", 400.0, 15.0),            # gypsum, low
        placeholder("Portlandite", 720.0, 12.0),   # portlandite, its own step
        placeholder("Cal", 950.0, 20.0; releases = :carbon_dioxide),
    ]
    grid = range(300.0, 1200.0; length = 451)

    @testset "the logistic is a release, with the properties of one" begin
        w = windows[2]
        @test released_fraction(w, value(w.midpoint)) ≈ 0.5
        @test released_fraction(w, value(w.midpoint) - 10 * value(w.width)) < 1.0e-4
        @test released_fraction(w, value(w.midpoint) + 10 * value(w.width)) > 1 - 1.0e-4
        # And the rate integrates to one, which is what makes it a fraction —
        # a peak that did not would put mass into or out of nothing.
        dT = step(grid)
        @test sum(released_rate(w, t) for t in grid) * dT ≈ 1.0 atol = 1.0e-4
        # Symmetric about the midpoint, and peaking there.
        @test released_rate(w, value(w.midpoint)) ≈ 1 / (4 * value(w.width))
    end

    @testset "it integrates to what ignition_loss says" begin
        tg = thermogram(st, windows; temperatures = grid)
        loss_total = ustrip(us"kg", ignition_loss(st).total)
        @test tg.loss[end] ≈ loss_total rtol = 1.0e-4

        # A LOGISTIC HAS INFINITE TAILS, so the curve does not start at exactly
        # zero: at 300 K the gypsum window is already 1.3 ‰ through. That is the
        # model, not an error, and asserting the identity is a stronger test
        # than a tolerance would be — it says the leak IS the tail.
        water = Dict(p.first => ustrip(us"kg", p.second) for p in bound_water_per_phase(st))
        co2 = Dict(
            k => ustrip(us"kg", v) for (k, v) in ChemistryLab._co2_per_phase(st)
        )
        tail = sum(
            (w.releases === :water ? water : co2)[w.phase] *
                released_fraction(w, first(grid)) for w in windows
        )
        @test tg.loss[1] ≈ tail rtol = 1.0e-10
        @test tail / loss_total < 1.0e-3            # and it is small, measured

        # The mass on the pan starts at the SOLID mass less whatever has already
        # gone, and ends at the residue.
        m0 = ustrip(us"kg", mass(st).solid)
        @test tg.mass[1] ≈ m0 - tg.loss[1] rtol = 1.0e-12
        @test tg.mass[end] ≈ m0 - loss_total rtol = 1.0e-4
        # DTG integrates to the same total, which is the consistency between the
        # curve and its derivative and is not automatic once three peaks are
        # summed.
        @test sum(tg.dtg) * step(grid) ≈ loss_total rtol = 1.0e-3
        # Three peaks, and they are where they were put.
        peaks = [grid[argmax([abs(released_rate(w, t)) for t in grid])] for w in windows]
        @test peaks ≈ [value(w.midpoint) for w in windows] atol = 2 * step(grid)
    end

    @testset "a phase with no window is named, not silently dropped" begin
        short = windows[1:2]                       # calcite left out
        missing_ones = phases_without_windows(st, short)
        @test ("Cal" => :carbon_dioxide) in missing_ones
        @test !(("Gp" => :water) in missing_ones)

        # COVERAGE IS PER PHASE **AND** PRODUCT. A carbonated hydrate carries
        # hydrogen and carbon, releases both, and at different temperatures —
        # so a phase covered for its water and not for its carbon dioxide is
        # half covered, and counting per phase alone would call it done.
        half = [
            placeholder("Gp", 400.0, 15.0),
            placeholder("Portlandite", 720.0, 12.0),
            placeholder("Cal", 500.0, 10.0),      # water window on a carbonate
        ]
        @test ("Cal" => :carbon_dioxide) in phases_without_windows(st, half)
        @test "Cal" in windows_without_phases(st, half)   # calcite releases no water
        # And the curve is short by exactly the carbon dioxide, which is the
        # silent failure this guards against.
        tg = thermogram(st, short; temperatures = grid)
        co2 = ustrip(us"kg", ignition_loss(st).carbon_dioxide)
        @test ustrip(us"kg", ignition_loss(st).total) - tg.loss[end] ≈ co2 rtol = 1.0e-3
        @test isempty(phases_without_windows(st, windows))

        # AND THE MIRROR, which is the one that catches a typo: a window on a
        # phase the state has nothing to release from contributes nothing and
        # raises nothing, so the only symptom would be a missing peak.
        typo = [placeholder("Portlandit", 720.0, 12.0)]      # one letter short
        @test windows_without_phases(st, typo) == ["Portlandit"]
        @test thermogram(st, typo; temperatures = grid).loss[end] == 0
        @test isempty(windows_without_phases(st, windows))
        # A window that names the right phase but the wrong product is caught
        # too: portlandite releases no carbon dioxide.
        wrong = [placeholder("Portlandite", 720.0, 12.0; releases = :carbon_dioxide)]
        @test windows_without_phases(st, wrong) == ["Portlandite"]
    end

    @testset "a phase that goes in stages" begin
        # Gypsum really does lose its two waters in two steps,
        # CaSO4·2H2O -> CaSO4·0.5H2O -> CaSO4, so a window per stage with the
        # fractions splitting the release is the ordinary case, not an exotic
        # one.
        staged = [
            DecompositionWindow(
                "Gp", 380.0, 8.0; fraction = 0.75,
                kind = PROV_PLACEHOLDER, source = "illustrative"
            ),
            DecompositionWindow(
                "Gp", 430.0, 8.0; fraction = 0.25,
                kind = PROV_PLACEHOLDER, source = "illustrative"
            ),
            windows[2], windows[3],
        ]
        tg = thermogram(st, staged; temperatures = grid)
        @test tg.loss[end] ≈ ustrip(us"kg", ignition_loss(st).total) rtol = 1.0e-4
        @test occursin("75.0 %", sprint(show, staged[1]))

        # AND TWO WINDOWS THAT EACH CLAIM ALL OF IT ARE REFUSED. Without the
        # check they would release the phase's mass twice, and the only symptom
        # would be a curve integrating to more than `ignition_loss` — a silent
        # doubling rather than an error.
        doubled = [windows[1], windows[1], windows[2], windows[3]]
        err = try
            thermogram(st, doubled; temperatures = grid)
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("summing to 2.0", sprint(showerror, err))
        @test occursin("twice", sprint(showerror, err))

        # A fraction outside (0, 1] is refused at construction.
        @test_throws ArgumentError DecompositionWindow("Gp", 400.0, 15.0; fraction = 0.0)
        @test_throws ArgumentError DecompositionWindow("Gp", 400.0, 15.0; fraction = 1.5)

        # And a round trip through the parameters keeps the split.
        θs, _ = window_parameters(staged)
        @test [w.fraction for w in with_window_parameters(staged, θs)] ==
            [w.fraction for w in staged]
    end

    @testset "what a window claims about itself" begin
        @test provenance(windows[1].midpoint) === PROV_PLACEHOLDER
        @test !is_evidence(windows[1].midpoint)
        @test occursin("placeholder", sprint(show, windows[1]))
        r = provenance_report([w.midpoint for w in windows])
        @test r.weakest === PROV_PLACEHOLDER
        @test !r.all_evidence

        # A bare number claims nothing rather than claiming to be measured.
        bare = DecompositionWindow("Gp", 400.0, 15.0)
        @test provenance(bare.midpoint) === PROV_UNSTATED

        @test_throws ArgumentError DecompositionWindow("Gp", 400.0, 0.0)
        @test_throws ArgumentError DecompositionWindow("Gp", 400.0, -5.0)
        @test_throws ArgumentError DecompositionWindow("Gp", 400.0, 15.0; releases = :argon)
    end

    @testset "parameters pack and unpack" begin
        θ, names = window_parameters(windows)
        @test length(θ) == 6
        @test names[1] == "T½(Gp)" && names[2] == "w(Gp)"
        @test θ[1] ≈ 400.0 && θ[2] ≈ 15.0
        back = with_window_parameters(windows, θ)
        @test [value(w.midpoint) for w in back] ≈ [value(w.midpoint) for w in windows]
        @test [w.releases for w in back] == [w.releases for w in windows]
        # Out of an optimizer, so fitted — which is what the round trip means.
        @test provenance(back[1].midpoint) === PROV_FITTED
        @test_throws ArgumentError with_window_parameters(windows, θ[1:3])
    end

    @testset "the windows are identified from the curve" begin
        # THE POINT OF THE WHOLE CHAIN. The windows are not consequences of the
        # formulas, so either a publication supplies them or a measurement does.
        # This is the second route, on a curve generated from known windows so
        # the answer is known: a real one would be a measured thermogram.
        truth = windows
        target = thermogram(st, truth; temperatures = grid).dtg

        forward(θ) = thermogram(
            st, with_window_parameters(truth, θ; kind = PROV_UNSTATED);
            temperatures = grid,
        ).dtg

        # Started well away from the answer — 40 K and 40 % off on every window.
        θ0 = [
            v * f for (v, f) in zip(
                    window_parameters(truth)[1],
                    [1.1, 1.4, 0.945, 0.6, 1.042, 1.4],
                )
        ]
        @test maximum(abs.(θ0 .- window_parameters(truth)[1])) > 30

        # Gauss-Newton in log space, on the same sensitivity matrix
        # `identifiability` uses — which is the point: the machinery that says
        # whether a parameter is determined is the machinery that determines it.
        θ = copy(θ0)
        for _ in 1:40
            J = log_sensitivity(forward, θ; relstep = 1.0e-3)
            r = forward(θ) .- target
            δ = -(J \ r)
            θ = θ .* exp.(clamp.(δ, -0.3, 0.3))
        end

        truth_θ = window_parameters(truth)[1]
        @test maximum(abs.(θ .- truth_θ) ./ truth_θ) < 1.0e-3

        # And what the curve actually determined, rather than what came out.
        id = identifiability(
            forward, θ; observed = target, names = window_parameters(truth)[2],
        )
        @test id.rank == 6                       # three separated peaks: all six
        @test id.rmse < 1.0e-9 * maximum(abs, target)
        out = as_traced(id, θ; source = "a synthetic thermogram")
        @test all(t -> provenance(t) === PROV_FITTED, out)
        @test all(!is_evidence, out)             # fitted, never evidence
        @info "TGA windows identified" worst_rel = maximum(abs.(θ .- truth_θ) ./ truth_θ) rank = id.rank
    end

    @testset "overlapping peaks are where identifiability earns its keep" begin
        # Two phases releasing in the same window is the real difficulty of a
        # thermogram — C-S-H, AFt and AFm all go below 200 °C — and a fit that
        # reported four numbers there would be reporting two.
        overlapped = [
            placeholder("Gp", 400.0, 15.0),
            placeholder("Portlandite", 405.0, 15.0),
        ]
        forward(θ) = thermogram(
            st, with_window_parameters(overlapped, θ; kind = PROV_UNSTATED);
            temperatures = grid,
        ).dtg
        θ, names = window_parameters(overlapped)
        id = identifiability(forward, θ; names = names)
        # Two of four, and not pinned to exactly two on purpose: the two
        # qualifying ratios are 10.1 and 11.2, close enough that a different
        # machine could swap which is the larger and answer one instead. What is
        # structural is that a single peak's worth of parameters is visible.
        @test id.rank <= 2                       # not four independent numbers
        @test id.condition > 20
        @info "overlapping windows" rank = id.rank condition = id.condition
    end
end
