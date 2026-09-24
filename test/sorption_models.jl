# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

using ChemistryLab
using ChemistryLab: value, source
using Test

# A PHREEQC sorption model written here on purpose. The point is to test the
# READER, and a published compilation is not this package's to carry — ClaySor
# 2023 is CC-BY-4.0 and freely available from its own Zenodo deposit, which is
# where a user gets it. The shape below is ClaySor's, down to the `error:` and
# `ref:` tags, because those are what a reader has to keep.
const _TINY_SORPTION = """
# A sorption model, in the shape a published one has.
# Consistent with some aqueous database, which is the point of saying so.

EXCHANGE_MASTER_SPECIES

Mntx\tMntx-\t # IonExchangeSite_mont

SURFACE_MASTER_SPECIES

Mnt_s\tMnt_sOH\t # EdgeSite_S_mont
Mnt_w\tMnt_wOH\t # EdgeSite_W2_mont

EXCHANGE_SPECIES

Mntx- = Mntx-\t# ionExchangeSite
\t-log_K\t0.0

Mntx- + Na+ = MntxNa
\t-log_k = 0.0

Ca+2 + 2 MntxNa = Mntx2Ca + 2 Na+ # Ca exchange error: 0.27 ref: Curti:2023:rep:
\t-log_K\t0.396

SURFACE_SPECIES

Mnt_sOH = Mnt_sOH\t# EdgeSite_S_mont
\t-log_K\t0.0

H+ + Mnt_sOH = Mnt_sOH2+ # Protonated strong site ref: Bradbury_Baeyens:1997:pap:
\t-log_K\t4.5

Mnt_sOH = Mnt_sO- + H+ # Deprotonated strong site ref: Bradbury_Baeyens:1997:pap:
\t-log_K\t-7.9

Ni+2 + Mnt_wOH = Mnt_wONi+ + H+ # Ni on the weak site error: 0.11 ref: Somebody:2020:pap:
\t-log_K\t-2.5

SOLUTION_SPECIES

H2O = H2O
\tlog_k 0.0
"""

_write_tiny() = let p = joinpath(mktempdir(), "tiny-sorption.dat")
    write(p, _TINY_SORPTION)
    p
end

@testsection "reading a published sorption model" begin

    @testset "what it finds" begin
        m = read_sorption_model(_write_tiny())
        @test sort(collect(keys(m.surfaces))) == ["Mnt_s", "Mnt_w"]
        @test collect(keys(m.exchangers)) == ["Mntx"]
        @test m.surfaces["Mnt_s"].reference == "Mnt_sOH"
        @test m.exchangers["Mntx"].reference == "Mntx-"
        @test occursin("EdgeSite_S_mont", m.surfaces["Mnt_s"].comment)
        # The block ends at the next keyword: nothing from SOLUTION_SPECIES.
        @test !any(
            r -> occursin("H2O = H2O", r.equation), m.surfaces["Mnt_s"].reactions,
        )
        @test occursin("sha256", m.source)
        # The header is kept because that is where a model says which aqueous
        # database its constants belong to.
        @test occursin("Consistent with some aqueous database", m.header)
        @test occursin("SorptionModel", sprint(show, m))
    end

    @testset "a reaction keeps its constant, its source and its uncertainty" begin
        m = read_sorption_model(_write_tiny())
        prot = only(
            r for r in m.surfaces["Mnt_s"].reactions
                if r.stoichiometry == Dict("H+" => -1, "Mnt_sOH" => -1, "Mnt_sOH2+" => 1)
        )
        @test value(prot.log_K) ≈ 4.5
        @test provenance(prot.log_K) === PROV_PUBLISHED
        @test source(prot.log_K) == "Bradbury_Baeyens:1997:pap:"
        # This entry states no error, and `nothing` is not zero: a constant that
        # says nothing about how well it is known is not a constant known
        # exactly.
        @test uncertainty(prot.log_K) === nothing

        ni = only(r for r in m.surfaces["Mnt_w"].reactions if haskey(r.stoichiometry, "Ni+2"))
        @test uncertainty(ni.log_K) ≈ 0.11
        @test source(ni.log_K) == "Somebody:2020:pap:"
        @test occursin("Ni on the weak site", ni.comment)
        # the tags are taken OUT of the plain comment rather than left in it
        @test !occursin("error:", ni.comment)
        @test !occursin("ref:", ni.comment)
    end

    @testset "a charge is not a separator" begin
        # `Ca+2 + 2 MntxNa = Mntx2Ca + 2 Na+` split on every `+` yields "Ca",
        # "2" and "2 MntxNa" — three terms, none of them the calcium ion. The
        # separator is a plus with whitespace on both sides.
        m = read_sorption_model(_write_tiny())
        ca = only(
            r for r in m.exchangers["Mntx"].reactions if haskey(r.stoichiometry, "Ca+2")
        )
        @test ca.stoichiometry ==
            Dict("Ca+2" => -1, "MntxNa" => -2, "Mntx2Ca" => 1, "Na+" => 2)
        @test value(ca.log_K) ≈ 0.396
        @test uncertainty(ca.log_K) ≈ 0.27
        @test source(ca.log_K) == "Curti:2023:rep:"
        # And `-log_k = 0.0`, the other spelling the same file uses.
        @test any(
            r -> haskey(r.stoichiometry, "MntxNa") && !haskey(r.stoichiometry, "Ca+2") &&
                value(r.log_K) == 0.0,
            m.exchangers["Mntx"].reactions,
        )
    end

    @testset "what the compilation is, as a whole" begin
        m = read_sorption_model(_write_tiny())
        r = provenance_report(log_constants(m))
        @test r.total == length(log_constants(m))
        @test r.weakest === PROV_PUBLISHED
        @test r.all_evidence
        # The number that matters when deciding what a model supports: how many
        # of its constants say nothing about how well they are known.
        @test r.without_uncertainty > 0
        @test r.without_uncertainty < r.total
    end

    @testset "taking a subset, by naming it" begin
        m = read_sorption_model(_write_tiny())
        # A compilation covering thirty elements is not something to import
        # whole; saying which part was taken is part of saying what was
        # reproduced.
        @test length(reactions_involving(m, "Ni+2")) == 1
        @test length(reactions_involving(m, "Ca+2")) == 1
        @test isempty(reactions_involving(m, "Am+3"))
        @test length(reactions_involving(m, "H+")) >= 2
    end

    @testset "refusals" begin
        @test_throws ArgumentError read_sorption_model("/nonexistent/path.dat")
        # A database with no sorption blocks is refused by name rather than
        # returning an empty model that would silently sorb nothing.
        plain = joinpath(mktempdir(), "plain.dat")
        write(plain, "SOLUTION_SPECIES\n    H2O = H2O\n    log_k 0.0\n")
        @test_throws ArgumentError read_sorption_model(plain)
    end
end
