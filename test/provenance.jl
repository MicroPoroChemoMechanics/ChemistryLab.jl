# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

using ChemistryLab
using ChemistryLab: value, source
using Test

@testsection "a number that says where it came from" begin

    @testset "the ordering is the claim" begin
        # Ordered from the weakest claim to the strongest, and `weakest` is what
        # a derived quantity can honestly say about itself.
        @test PROV_UNSTATED < PROV_PLACEHOLDER < PROV_ESTIMATED <
            PROV_FITTED < PROV_PUBLISHED < PROV_MEASURED

        # An UNSTATED value is the weakest of all, deliberately: a number that
        # forgot to say where it came from must never strengthen a result.
        @test weakest(Traced(1.0, PROV_MEASURED), Traced(2.0)) === PROV_UNSTATED
        @test weakest(
            Traced(1.0, PROV_MEASURED, "this work"),
            Traced(2.0, PROV_PUBLISHED, "a paper"),
        ) === PROV_PUBLISHED
        @test weakest() === PROV_UNSTATED
        @test weakest(Traced(1.0, PROV_FITTED), Traced(2.0, PROV_PLACEHOLDER)) ===
            PROV_PLACEHOLDER
    end

    @testset "what counts as evidence, and what does not" begin
        @test is_evidence(Traced(1.0, PROV_MEASURED, "this work"))
        @test is_evidence(Traced(1.0, PROV_PUBLISHED, "a paper"))
        # A FITTED value is not evidence here. It may be excellent; whether it
        # is depends on the data, the model and the parameter's identifiability,
        # and this predicate cannot make that judgement.
        @test !is_evidence(Traced(1.0, PROV_FITTED, "a calibration"))
        @test !is_evidence(Traced(1.0, PROV_ESTIMATED, "an analog"))
        @test !is_evidence(Traced(1.0, PROV_PLACEHOLDER, "pending a measurement"))
        @test !is_evidence(Traced(1.0))
        # A bare number claims nothing rather than claiming to be measured.
        @test !is_evidence(1.0)
        @test provenance(1.0) === PROV_UNSTATED
        @test source(1.0) == ""
    end

    @testset "it is not a Real, and that is the design" begin
        # Flowing silently into arithmetic is exactly how provenance gets lost,
        # so it does not: unwrapping is an act a reader can see.
        t = Traced(2.0, PROV_MEASURED, "this work")
        @test !(t isa Real)
        @test_throws MethodError t + 1.0
        @test value(t) + 1.0 == 3.0
        @test value(3.5) == 3.5            # a plain number unwraps to itself
    end

    @testset "it prints what it is" begin
        @test occursin("measured", sprint(show, Traced(1.0, PROV_MEASURED)))
        @test occursin("this work", sprint(show, Traced(1.0, PROV_MEASURED, "this work")))
        @test occursin("placeholder", sprint(show, Traced(1.0, PROV_PLACEHOLDER)))
        # with no source, no stray separator
        @test !occursin(": ]", sprint(show, Traced(1.0, PROV_MEASURED)))
    end

    @testset "a report beside a result" begin
        # A table where nine coefficients are published and one is a placeholder
        # is a different object from one where all ten are published, and the
        # difference shows in none of the numbers.
        nine = [Traced(float(i), PROV_PUBLISHED, "a paper") for i in 1:9]
        r = provenance_report(vcat(nine, [Traced(0.0, PROV_PLACEHOLDER, "pending")]))
        @test r.total == 10
        @test r.counts[PROV_PUBLISHED] == 9
        @test r.counts[PROV_PLACEHOLDER] == 1
        @test r.weakest === PROV_PLACEHOLDER
        @test !r.all_evidence

        clean = provenance_report(nine)
        @test clean.all_evidence
        @test clean.weakest === PROV_PUBLISHED

        empty_r = provenance_report(Traced{Float64}[])
        @test empty_r.total == 0
        @test empty_r.weakest === PROV_UNSTATED
        @test !empty_r.all_evidence       # nothing is not everything
    end

    @testset "it carries any value type" begin
        # Including a Dual, which is what identifying a parameter needs, and a
        # dimensioned quantity, which is what the rest of this package speaks.
        @test value(Traced(1 // 2, PROV_MEASURED)) == 1 // 2
        @test provenance(Traced([1.0, 2.0], PROV_FITTED, "a fit")) === PROV_FITTED
    end
end
