# [Validation against a silica fume shotcrete paste (Lothenbach et al. 2014)](@id sec-validation-silica-fume)

!!! info "Before this page"
    [Validation against a measured paste](@ref ex-validation), for the
    Portland paste whose phase list this page reuses, and [the CASH+
    model](@ref ex-cashplus-csh), for the C-S-H model it tests.

[Lothenbach2014](@citet) followed a low-pH shotcrete cement, ESDRED, made for
contact with the clay barrier of a repository: 60 % CEM I 42.5 N and 40 % silica
fume, with 4.8 g of an aluminum sulfate set accelerator and 1.2 g of a
polycarboxylate superplasticizer per 100 g of binder, at a water/binder ratio
of 0.5, sealed at 20 °C for 3.5 years. They analyzed the pore solution from one
hour on and followed the clinker and the silica fume by ²⁹Si NMR.
[Miron2022b](@citet) computed this paste to test the alkali end-members of the
CASH+NK model of the C-S-H on a binder whose gel ends at a low Ca/Si. This page
computes it at the measured degrees of reaction, with the C-S-H as CSHQ and as
CASH+NK, and compares the pore solutions.

```@example esdred
using ChemistryLab, DynamicQuantities, Logging, OptimaSolver, Printf, Plots
default(framestyle = :box, grid = false)
include(joinpath(pkgdir(ChemistryLab), "scripts", "esdred_2014.jl"))
nothing # hide
```

## 1. The binder and how far it has reacted

The silicates of the clinker and the silica fume are followed by NMR: the
clinker's share of the silicon, and the silica fume left unreacted, which falls
from 40 to 10 g per 100 g of binder in 3.5 years
([Lothenbach2014; Table 2](@cite)). They react here at those degrees; the
aluminate and the ferrite, which the NMR does not see, by the law of
[ParrottKilloh1984](@citet) with the constants of this cement.

```@example esdred
nmr = es14_nmr_extents()
println("  days   alite and belite   silica fume")
for d in (1, 7, 28, 360, 1310)
    @printf("%6d   %16.2f %13.2f\n", d, extent(nmr.silicates, d), extent(nmr.silica_fume, d))
end
```

The binder is entered as [Lothenbach2014](@citet) describe it: the normative
phases of the CEM I (their Table 1, which sums to 100.31 g and is scaled to
100), the silica fume by its oxides, and the accelerator as an addition, its
aluminum, sulfate and alkalis dissolved from the start. [Miron2022b](@citet)
describe ESDRED as 40 % CEM I and 60 % silica fume; the paper they cite, and its
Table 2, give 60 and 40, which is used here. The assumptions the script makes
are listed at its head, `scripts/esdred_2014.jl`.

## 2. The formate of the accelerator

The accelerator brings organic carbon, which [Lothenbach2014](@citet) identify
as formate: 202 mM in the mixing water at the start, of which they estimate from
the charge balance of each pore solution how much is still dissolved, 140 mM
after one day and about 80 mM after a year (footnote a of their Table 3).
Formate is an anion, and at these concentrations it is the main one. The
measured solution after 28 days balances its charge only with it:

```@example esdred
d = 28.0
cations = es14_measured(d, "K") + es14_measured(d, "Na") + 2es14_measured(d, "Ca")
anions = es14_measured(d, "OH-") + 2es14_measured(d, "S") + es14_measured(d, "formate")
@printf("at %d days: cations %.0f meq/L; anions %.0f meq/L, of which formate %.0f\n",
        d, cations, anions, es14_measured(d, "formate"))
```

Cemdata18 [Lothenbach2019](@cite) has no formate, and [Miron2022b](@citet) left
it out, naming it as a possible cause of the gaps they found. The formate of the
SUPCRT [Johnson1992](@cite) organic database can be added. It then has to stay
formate, as it does in the paste: the reduced species of sulfur and iron of the
database are left out, so that nothing in the system can oxidize it. But the
model has no solid that takes formate up, where [Lothenbach2014](@citet)
estimate that the solids take a third to three-fifths of it. The two
calculations below are therefore two bounds: without the formate, and with all
of it in solution.

```@example esdred
model = cemdata18_activity_model(:KOH)
days = es14_days()
diagnostics = IOBuffer() # hide
runs = with_logger(ConsoleLogger(diagnostics)) do # hide
runs = Dict((gel, f) => hydrate(es14_recipe(; formate = f), es14_system(gel; formate = f), days; model)
            for gel in (:CSHQ, :CASHNK), f in (false, true))
end # hide
for gel in (:CSHQ, :CASHNK), f in (false, true)
    @printf("%-7s %-15s certified at %d of %d ages\n", gel, f ? "with formate" : "without formate",
            count(rs -> rs.certificate.optimal, runs[(gel, f)].states), length(days))
end
```

## 3. The pore solutions

```@example esdred
gel_name(gel) = gel === :CSHQ ? "CSHQ" : "CASH+NK"
println("  days  element  measured   CSHQ: no formate  formate   CASH+NK: no formate  formate")
for (k, d) in enumerate(days), e in ("K", "Na", "Ca", "pH")
    m = es14_measured(d, e)
    v = [e == "pH" ? pH(runs[(g, f)].states[k].state, model) : pore_solution_mmol(runs[(g, f)].states[k].state)[e]
         for g in (:CSHQ, :CASHNK) for f in (false, true)]
    @printf("%6.0f  %-7s %9.3g %15.3g %9.3g %18.3g %9.3g\n", d, e, m, v...)
end
```

```@example esdred
panels = map(("K", "Ca", "pH")) do e
    # The legend sits above the curves of calcium, its axis raised to make room.
    p = plot(; xscale = :log10, title = e == "pH" ? "pH" : "$e (mmol/L or mmol/kg)", xlabel = "time (days)",
             legend = e == "Ca" ? :top : false, ylims = e == "Ca" ? (0, 240) : :auto)
    for (gel, color) in ((:CSHQ, :darkorange), (:CASHNK, :steelblue)), (f, style) in ((false, :solid), (true, :dash))
        y = [e == "pH" ? pH(rs.state, model) : pore_solution_mmol(rs.state)[e] for rs in runs[(gel, f)].states]
        plot!(p, days, y; color, linestyle = style, linewidth = 2,
              label = "$(gel_name(gel)), " * (f ? "all formate" : "no formate"))
    end
    scatter!(p, days, [es14_measured(d, e) for d in days]; color = :black, label = "measured")
    p
end
fig = plot(panels...; layout = (1, 3), size = (1000, 380), left_margin = 6Plots.mm, bottom_margin = 7Plots.mm)
savefig(fig, "esdred-pore-solution.svg"); nothing # hide
```

![The potassium, the calcium and the pH of the pore solution against time, measured, and computed with CSHQ and CASH+NK, without the formate and with all of it in solution.](esdred-pore-solution.svg)

## 4. What the comparison says

Without the formate, CASH+NK follows the fall of the alkalis as the silica fume
reacts and decalcifies the gel, from two weeks on: the potassium comes down to
8 mmol/kg after 3.5 years, where 11 mmol/L were measured, while CSHQ keeps 96.
That is the result [Miron2022b](@citet) report for their fine-tuned model. In
the first days CASH+NK takes too much of it, 43 to 59 mmol/kg where 139 to 236
were measured, and the sodium stays below the measurement throughout with both
gels. The calcium is reproduced at 28 and 56 days by CASH+NK, 13 mmol/kg against
13 and 20 mmol/L, and by neither later, where the measured solution holds 27 to
29 and the computed ones 0.6 to 5. Between two weeks and two months the pH falls
below the measured one, by half a unit to a unit with CSHQ and by one to two with
CASH+NK; later CASH+NK comes back within a unit of it and CSHQ goes above it.

With all the formate in solution the calcium and the alkalis rise past the
measured ones, the formate having to be balanced by cations, and the pH falls
further. After the first weeks the measured calcium lies between the two
bounds, as it should if the solids take up part of the formate. No model of the
package, nor of Cemdata18, describes that uptake by the hydrates, and a
calculation of this paste is only as good as what replaces it. What this page
settles is narrower: the late alkalis of a silica fume binder are followed by
CASH+NK and not by CSHQ, and its calcium and pH are not followed by either
without an account of the formate of its accelerator.

## Where to go next

[The CASH+ model](@ref ex-cashplus-csh) computes the Portland cement paste
[Miron2022b](@citet) also used, and [the CEM IV page](@ref ex-cem4-pozzolanic) a
pozzolanic binder with CSHQ and CNASH_ss.
