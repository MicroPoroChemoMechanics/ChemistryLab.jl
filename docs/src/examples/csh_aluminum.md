# [Aluminum in the C-S-H: two ways that were tried, and why neither is shipped](@id ex-csh-aluminum)

!!! info "Before this page"
    [Solid solutions in a calculation](@ref sec-tutorial-solid-solutions), for
    how a gel is declared and what its end-members can hold, and [Aluminum
    uptake by C-S-H](@ref sec-validation-aluminum-uptake), which computes the
    syntheses used below with the two gels of Cemdata18.

!!! warning "What this page fits, and what it ships"
    The two extensions below are ChemistryLab's own, each with one energy fitted
    on published measurements. Neither is part of Cemdata18, of a published model
    or of the package's databases: they live in `scripts/csh_aluminum.jl`, and
    this page records what they give, so that the question of a C-S-H with
    aluminum at a high Ca/Si is documented rather than answered by a model that
    has not earned it.

The calcium silicate hydrate of a hydrated cement holds aluminum. In the pastes
of [DeWeerdt2011](@citet), measured by electron microprobe, the gel of a
Portland cement has a Ca/Si of 1.8 ± 0.1 and an Al/Si of 0.06 ± 0.01; with fly
ash it moves to about 1.4 and 0.13 after 140 days. A model of the gel for such
binders has to hold both the calcium of a gel beside portlandite and some
aluminum, and none of the gels ChemistryLab ships does. CSHQ [Kulik2011](@cite)
and the CASH+ family [Kulik2022, Miron2022a](@cite) hold no aluminum; CNASH_ss
[Myers2014](@cite) holds it, and its Ca/Si stops at 1.5, that of its
calcium-richest end-member, its gel reaching 1.16 in the pastes of De Weerdt et
al. where 1.8 was measured ([the validation page](@ref sec-validation-blended-gels)).

## 1. What the literature offers

The authors of CASH+ describe how aluminum would enter their model, as an
aluminate on the bridging site of the silicate chain, and cite the extension as
in preparation [Kulik2022; ref. 25](@cite). [Yan2022](@citet) compare their
measurements with a CASH+ model that holds aluminum, cited in turn as in
preparation. A search of Crossref on 2026-10-06 finds no publication of it, and
its parameters are not available, so that ChemistryLab can neither use it nor
recompute the curves of Yan et al.

The measurements are another matter. [LHopital2015](@citet) and
[LHopital2016a](@citet) synthesized C-A-S-H at 20 °C at a Ca/Si from 0.6 to 1.6
and an Al/Si up to 0.33, without alkali and, for the first, at a Ca/Si of 1.0 in
0.5 M KOH; [LHopital2016b](@citet) at an Al/Si of 0.05 in KOH and NaOH solutions
up to 0.5 M; [Yan2022](@citet) at a Ca/Si of 1.0 in NaOH and KOH up to 1 M.
Each reports the gel by mass balance, its Al/Si included, and the solution.
Two other sources were examined and not used: [Haas2015](@citet) give a model of
their own, by surface reactions, on C-A-S-H made by adding a C-S-H to a calcium
aluminate solution for a week, the kind of short experiment that L'Hôpital et
al. found to put more aluminum in the gel than equilibrium does; and
[Roosz2018](@citet) measured two C-A-S-H by calorimetry, whose aluminum, an
Al/Si of 0.04 and 0.016, is too little to fix the energy of an aluminum
end-member.

```@example cah
using ChemistryLab, DynamicQuantities, Logging, OptimaSolver, Printf, Plots
default(framestyle = :box, grid = false)
include(joinpath(pkgdir(ChemistryLab), "scripts", "csh_aluminum.jl"))
nothing # hide
```

## 2. What the syntheses say before any model

Two readings of the measurements decide how a model can be fitted on them.

The first concerns the other solids. Beyond an Al/Si of about 0.05 the
syntheses hold strätlingite and katoite beside the gel, and the Al/Si of the
gel stays near 0.1, which reads like the limit a model should reproduce. The
saturation indices L'Hôpital et al. compute from their own solutions say
otherwise: those solutions are undersaturated with respect to the strätlingite
they hold, by 1.2 to 3.4 log units, and to the katoite, by 6 to 11 (their
Appendix C, recorded in `data/literature/LHopital2015.json`). These phases are
not at equilibrium with the solution, and the Al/Si the gel keeps beside them is
not an equilibrium datum. What is one is the pair the gel and its solution
form: the Al/Si of the gel, measured by mass balance, against the aluminum
dissolved beside it. That pair, the uptake isotherm, is what the models below
are fitted and tested on, with the aluminum-bearing pure phases left out of the
calculation.

The second concerns the shape of that isotherm. Without alkali, at a Ca/Si of
1.0 and below an Al/Si of 0.05, where the gel holds all the aluminum:

```@example cah
g = literature_table("LHopital2015", "gel_composition")
p = literature_table("LHopital2015", "pore_solution")
pts = [(g.Al_Si[i], ustrip(us"mmol/L", p.concentration[j])) for i in eachindex(g.Al_Si)
       for j in eachindex(p.element)
       if !ismissing(g.Al_Si[i]) && g.Al_Si_target[i] <= 0.05 && p.Al_Si_target[j] == g.Al_Si_target[i] &&
          p.time_d[j] == g.time_d[i] && p.element[j] == "Al" && p.qualifier[j] == "measured"]
x, y = log10.(first.(pts)), log10.(last.(pts))
slope = sum((x .- sum(x) / length(x)) .* (y .- sum(y) / length(y))) / sum(abs2, x .- sum(x) / length(x))
@printf("%d samples, slope of log(Al dissolved) against log(gel Al/Si): %.2f\n", length(pts), slope)
fig = scatter(first.(pts), last.(pts); xscale = :log10, yscale = :log10, color = :steelblue, label = "measured",
              xlabel = "Al/Si of the gel", ylabel = "dissolved aluminum (mmol/L)",
              title = "C-A-S-H at Ca/Si 1.0, no alkali, 20 °C", legend = :topleft, size = (680, 420),
              left_margin = 6Plots.mm, bottom_margin = 6Plots.mm)
xs = [0.01, 0.06]
plot!(fig, xs, 0.4 .* xs; color = :gray, linestyle = :dash, label = "proportional")
savefig(fig, "cah-isotherm.svg"); nothing # hide
```

![The dissolved aluminum against the Al/Si of the gel, measured, on logarithmic axes, with a line of slope one.](cah-isotherm.svg)

The dissolved aluminum is close to proportional to the aluminum of the gel, the
slope on logarithmic axes being 1.3. In an ideal mixture an end-member carrying
``\nu`` atoms of aluminum per formula unit gives a dissolved aluminum
proportional to the power ``1/\nu`` of its mole fraction, a slope of ``1/\nu``:
the isotherm asks for about one aluminum per formula unit. 5CA, the aluminum
end-member of CNASH_ss, carries 0.25 as Cemdata18 writes it, and taken as an
ideal member it would give a slope of 4.

## 3. Two ways, one energy each

**Way A** extends CASH+, whose gel is a mixture on two sites: the bridging
tetrahedron of the silicate chain (a silicate, a vacancy or a calcium) and the
interlayer cation (a proton, a calcium, a sodium or a potassium), every
combination a compound [Kulik2022](@cite). It adds the aluminate the authors
describe for the bridging site, AlO(OH)₂⁻, which carries the same charge as the
silicate it replaces; their Table 1 lists Al(OH)₄⁻ instead, which differs by a
molecule of water and changes the water of the gel, not its solubility in
dilute solutions. Each new compound is the silicate compound of the same
interlayer occupant with its bridging silicate exchanged:

```math
\mathrm{TSvh} + \mathrm{Al(OH)_3} \;\longrightarrow\; \mathrm{TAvh} + \mathrm{SiO_2} + \mathrm{H_2O},
```

gibbsite and amorphous silica of Cemdata18 on either side, so that TAvh is
Ca₂Si₂AlO₁₁H₇, TACh Ca₃Si₂AlO₁₃H₉, and so on with sodium and potassium. Its
Gibbs energy is that of the reaction's left-hand side minus the right-hand side's,
plus one energy ``\delta``, the same for every interlayer occupant; its entropy,
heat capacity and volume follow the reaction unchanged, and the aluminate is
given no site interaction of its own. One aluminum per compound: the isotherm
of section 2 is what the formalism gives.

**Way B** keeps CSHQ, the ideal mixture Cemdata18 ships, and adds one
end-member: four formula units of 5CA, (CaO)₅(SiO₂)₄(Al₂O₃)₀.₅(H₂O)₆.₅, which
carries one aluminum, its entropy, heat capacity and volume four times those of
5CA, its Gibbs energy four times that of 5CA plus ``\delta``.

Each way has one energy to fit, ``\delta``, and it is fitted on the syntheses of
L'Hôpital et al. (2016a) without alkali, at a Ca/Si from 0.6 to 1.6: the
dissolved aluminum the model gives for the gel at its measured Al/Si, against
the dissolved aluminum measured, in decimal logarithm. A synthesis whose
aluminum was below the detection limit enters only by how far the model exceeds
that limit.

```@example cah
fit = csh_al_samples("fit")
δs = Dict(:A => [8.0, 9.0, 9.8, 11.0, 12.0], :B => [-10.0, -5.0, 0.0, 5.0, 10.0])   # kJ/mol
scan = Dict(way => [csh_al_misfit(csh_al_isotherm(way, 1000δ, fit)) for δ in δs[way]] for way in (:A, :B))
println("way   δ (kJ/mol)   rms log   bias    above a detection limit")
for way in (:A, :B), (δ, m) in zip(δs[way], scan[way])
    @printf("%-5s %8.1f   %8.2f   %+6.2f   %6.2f\n", way, δ, m.rms, m.bias, m.excess)
end
best = csh_al_isotherm(:A, 9.8e3, fit)
println("\nway A at 9.8 kJ/mol, by the Ca/Si of the synthesis: bias (log), or the excess over a detection limit")
for ca in sort(unique(x.Ca_Si for x in best))
    m = csh_al_misfit([x for x in best if x.Ca_Si == ca])
    @printf("Ca/Si %.1f: %s\n", ca, m.n > 0 ? @sprintf("bias %+.2f on %d", m.bias, m.n) : @sprintf("%.2f above the limit", m.excess))
end
```

```@example cah
fig = plot(δs[:A], [m.rms for m in scan[:A]]; marker = :circle, color = :steelblue, linewidth = 2,
           label = "way A, CASH+ with the aluminate", xlabel = "exchange energy δ (kJ/mol)",
           ylabel = "rms of log(computed / measured Al)", title = "fit on L'Hôpital et al. 2016a, no alkali",
           legend = :top, size = (680, 420), left_margin = 6Plots.mm, bottom_margin = 6Plots.mm)
plot!(fig, δs[:B], [m.rms for m in scan[:B]]; marker = :diamond, color = :darkorange, linewidth = 2,
      label = "way B, CSHQ with 4 × 5CA")
savefig(fig, "cah-fit.svg"); nothing # hide
```

![The root mean square of the logarithmic misfit of the dissolved aluminum against the exchange energy, for the two ways.](cah-fit.svg)

Way A is fitted at ``\delta = 9.8`` kJ/mol, where the fifteen measured
aluminum concentrations are reproduced within 0.28 in decimal logarithm, a
factor of 1.9, without a trend in the Ca/Si from 1.0 to 1.6. Way B does less
than half as well at its best, ``\delta = 0``. Both fail on the syntheses of
Ca/Si 0.6 and 0.8, where the measured aluminum is below
the detection limit of 0.0037 mmol/L: way A puts about 200 times more in
solution there, way B more still, its calcium-rich member having no place in a
gel beside amorphous silica. A second energy for way A, one for the aluminate
balanced by a proton and one for that balanced by calcium, was tried to reach
those syntheses as well; the values it needs make the mixture of the two sites
non-convex, so that the gel would unmix, and most equilibria no longer certify.
It was not pursued.

## 4. On data the fit has not seen

The energies are now held, and the isotherm computed for syntheses in alkali
hydroxide solutions, from three series none of the fits used:

```@example cah
δbest = Dict(:A => 9.8, :B => 0.0)
sets = ["lh16b" => "L'Hôpital et al. 2016b", "lh15" => "L'Hôpital et al. 2015", "yan" => "Yan et al. 2022"]
# The equilibria that do not certify warn; they are counted in the table and
# left out of the misfit.
diagnostics = IOBuffer() # hide
val = with_logger(ConsoleLogger(diagnostics)) do # hide
val = Dict((way, set) => csh_al_isotherm(way, 1000δbest[way], csh_al_samples(set)) for way in (:A, :B), (set, _) in sets)
end # hide
println("series                     way   n    rms log   bias   certified")
for (set, label) in sets, way in (:A, :B)
    r = val[(way, set)]
    m = csh_al_misfit(r)
    @printf("%-26s %-4s %3d   %7.2f  %+6.2f   %d/%d\n", label, way, m.n, m.rms, m.bias, count(x -> x.certified, r), length(r))
end
```

```@example cah
println("L'Hôpital et al. 2016b, way A, by the Ca/Si of the synthesis")
println("Ca/Si   n   bias (log)")
lh16b = val[(:A, "lh16b")]
for ca in sort(unique(x.Ca_Si for x in lh16b))
    m = csh_al_misfit([x for x in lh16b if x.Ca_Si == ca])
    @printf("%4.1f  %3d   %+6.2f\n", ca, m.n, m.bias)
end
```

```@example cah
fig = plot([1.0e-3, 10.0], [1.0e-3, 10.0]; xscale = :log10, yscale = :log10, color = :gray, linestyle = :dash,
           label = "", xlabel = "measured aluminum (mmol/L)", ylabel = "computed aluminum (mmol/L)",
           title = "way A, δ = 9.8 kJ/mol, alkaline syntheses", legend = :topleft, size = (680, 480),
           left_margin = 6Plots.mm, bottom_margin = 6Plots.mm)
for ((set, label), (color, shape)) in zip(sets, [(:steelblue, :circle), (:darkorange, :diamond), (:seagreen, :utriangle)])
    shown = [x for x in val[(:A, set)] if x.certified && x.qualifier == "measured"]
    scatter!(fig, [x.Al for x in shown], [x.model for x in shown]; color, markershape = shape, label)
end
savefig(fig, "cah-validation.svg"); nothing # hide
```

![The aluminum computed by way A against the aluminum measured in the alkaline syntheses of three series, with the line of equality.](cah-validation.svg)

Way A carries over to the alkaline syntheses at a Ca/Si of 1.0: the series of
Yan et al., in NaOH and KOH up to 1 M, within about 0.4 in decimal logarithm,
and those of L'Hôpital et al. (2016b) at that Ca/Si with a bias of −0.07. Away from it,
it does not: at a higher Ca/Si it puts more aluminum in solution than was
measured, increasingly so up to an order of magnitude at 1.6, and at a lower one
the excess of the alkali-free syntheses returns. Way B is further off on the two
larger series and about as close on the four syntheses of L'Hôpital et al.
(2015). The
interlayer occupants sodium and potassium were given the same ``\delta`` as the
proton and the calcium, an assumption these series test rather than fit. The
equilibria that do not certify, eight of 87 in the first series under way A, are
counted in the table and left out of the misfit.

## 5. In a cement paste

The test that matters for a cement is a paste, where the gel shares the
aluminum with the AFm, AFt and hydrogarnet phases. The pastes of De Weerdt et
al. at 140 days, computed as on [the validation page](@ref ex-validation-blended)
with the gel of way A in place of CSHQ, against the gel they measured:

```@example cah
include(joinpath(pkgdir(ChemistryLab), "scripts", "de_weerdt_2011.jl"))
cashplus = build_species(datapath("cemdata18-cashplus.json"); verbose = false)
base = phase_list_system(DW11_PHASES, cashplus; replace = Dict("CSHQ" => "CASH+NK"), exclude_aqueous = DW11_CASHPLUS_EXCLUDED)
gelA = csh_al_gel_A(9.8e3, ["v", "C", "N", "K"])
cs_paste = ChemicalSystem(
    vcat(base.species, [m for m in end_members(gelA) if !(symbol(m) in symbol.(base.species))]), CEMDATA_PRIMARIES;
    solid_solutions = vcat([p for p in base.solid_solutions if ChemistryLab.name(p) != "CASH+NK"], [gelA]),
)
paste_model = cemdata18_activity_model(:KOH)
println("paste       gel Ca/Si   gel Al/Si   measured Ca/Si, Al/Si   certified")
measured = Dict("OPC" => "1.8, 0.06", "OPC-L" => "1.8, 0.06", "OPC-FA" => "1.4, 0.13", "OPC-FA-L" => "1.4, 0.13")
for mix in DW11_MIXES
    rs = hydrate(dw11_recipe(mix), cs_paste, [140.0]; model = paste_model)
    e = solid_solution_totals(rs[1].state, "C-S-H").elements
    @printf("%-10s %9.2f %11.3f   %20s     %s\n", mix, e[:Ca] / e[:Si], e[:Al] / e[:Si], measured[mix], rs[1].certificate.optimal)
end
```

The gel of way A reaches a Ca/Si of about 1.5, where the measured one is 1.8
without fly ash and 1.4 with it; that is the reach of CASH+ itself, not of the
aluminate. Its aluminum is right in one paste of the four, the CEM II/B-V
(OPC-FA), twice the measured one in the CEM I, and nearly nothing in the two
pastes with limestone, where the gel was measured to hold 0.06 and 0.13 and the
calculation puts the aluminum in the other aluminate hydrates. The energy fitted
on syntheses does not carry over to the pastes, where the aluminum of the gel is
set by its competition with those hydrates.

## 6. What is proposed, and why

Neither way is shipped. Way A is the better one by every test: it reproduces
the uptake isotherm of the syntheses at a Ca/Si from 1.0 to 1.6 with one
energy, and carries over to alkaline syntheses at a Ca/Si of 1.0. It fails
below a Ca/Si of 1.0 and in alkaline solutions above it, and in cement pastes
it gives the gel a wrong amount of aluminum in three of four. A fitted
end-member in a database looks exactly like a published one, and this one would
be used precisely where it fails: in pastes, at a high Ca/Si, in alkaline pore
solutions.

What a calculation can do meanwhile is what [the aluminum-uptake page](@ref
sec-validation-aluminum-uptake) concludes: CNASH_ss where the aluminum held by
the gel matters, at a Ca/Si that cannot exceed about 1.2; CSHQ or CASH+NK where
the calcium of a Portland cement's gel matters, the aluminum then going to the
hydrates. What would settle the question is the extension of CASH+ its authors
have announced, fitted together with the energies of the aluminate hydrates,
and measurements of the gel's isotherm at a high Ca/Si in alkaline solutions,
which none of the series above reaches beyond 1.6. `scripts/csh_aluminum.jl`
reproduces everything on this page, and the transcriptions it reads,
`data/literature/LHopital2015.json` and `data/literature/Yan2022.json`, record
how each value was read and checked.

## Where to go next

[Aluminum uptake by C-S-H](@ref sec-validation-aluminum-uptake) and [Alkali
uptake by C-A-S-H](@ref sec-validation-alkali-uptake) compute the syntheses with
the gels the package ships, and [The CASH+ model](@ref ex-cashplus-csh) the gel
way A extends.
