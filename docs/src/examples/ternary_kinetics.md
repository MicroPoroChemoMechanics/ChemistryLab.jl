# [CEM II/B-V and CEM II/B-M (V-LL), with their CEM I, integrated in time](@id ex-ternary-kinetics)

!!! info "Before this page"
    [Validation against measured blended pastes](@ref ex-validation-blended),
    the same four pastes at the extents their authors measured, and
    [CEM I 52.5 R with slag, fly ash and limestone, integrated in time](@ref ex-quaternary-kinetics),
    the same path on another blend.

[DeWeerdt2011](@citet) blended a clinker interground with gypsum with 5 %
limestone powder, 35 % siliceous fly ash, or 30 % fly ash and 5 % limestone, and
followed the pastes at w/b = 0.5 and 20 °C for six months by X-ray diffraction:
the clinker phases left, and the portlandite and the ettringite formed (their
Table 7). [The validation page](@ref ex-validation-blended) computes these
pastes at the degrees of reaction they measured. Here the clinker reacts under
a rate law, and only the fly ash follows the measurement: what the pastes form
over six months is then a prediction of the law and of the thermodynamics
together, which nothing below was fitted to.

## 1. The laws

The materials are the templates of the validation page: the clinker by its
Rietveld phases and its minor oxides, the fly ash by its crystals, inert, and
its glass, the limestone by its calcite. Each clinker phase dissolves under the
law of [ParrottKilloh1984](@citet) with the parameters
[Lavergne2018](@citet) tabulate, at the Blaine fineness of the cement (Table 1)
and at the water/clinker ratio of the paste. The glass of the fly ash reacts as
the authors measured it by image analysis: the fit printed on their Fig. 7,
in percent of the fly ash at `t` days,

```@example ternary
using ChemistryLab, DynamicQuantities, OptimaSolver, Printf
using Logging # hide
include(joinpath(pkgdir(ChemistryLab), "scripts", "deweerdt2011_kinetics.jl"))

a, b, c = dw11_value("fly_ash_fit_a"), dw11_value("fly_ash_fit_b"), dw11_value("fly_ash_fit_c")
@printf("y = %g + %g ln(t + %g)\n", a, b, c)
```

written on the degree of reaction `α` of the glass rather than on the time, as
the Waller law is: with `s` the glass's share of the fly ash, `y = 100 s α`,
`t + c = exp((y − a)/b)`, and the rate is `dα/dt = b exp(−(y − a)/b) / (100 s)`
per day (`dw11k_fly_ash_law`). The limestone, the gypsum and everything else are
at equilibrium.

ASSUMED: the minor oxides of the clinker (free lime, the alkalis, magnesia and
the sulfate of the clinker itself) enter the equilibrium at the mixing. The
validation page releases them with the clinker as a whole.

## 2. The four pastes

```@example ternary
setup = dw11k_setup()
day = 86400.0
diagnostics = IOBuffer() # hide
# `Any`: the type of a run is long enough that a dictionary specialized on it
# takes minutes to compile.
runs = with_logger(ConsoleLogger(diagnostics)) do # hide
runs = Dict{String, Any}(mix => dw11k_run(setup, mix) for mix in DW11_MIXES)
end # hide
occursin("re-speciation failed", String(take!(diagnostics))) && error("a re-speciation failed") # hide
for mix in DW11_MIXES
    @printf("%-9s %s, %3d steps\n", mix, runs[mix].sol.retcode, length(runs[mix].sol.t))
end
```

Neither law reads the equilibrium partition, so the trajectory does not depend
on it: every number below is computed on the certified replay of each run
([`speciated_states`](@ref)), in percent of the solids of the paste, every solid
with its water and the unreacted part of the binder included.

## 3. The clinker phases

```@example ternary
days = [1, 7, 28, 90, 180]
states = with_logger(ConsoleLogger(diagnostics)) do # hide
states = Dict{String, Any}(mix => dw11k_replay(runs[mix], days) for mix in DW11_MIXES)
end # hide
contents = Dict(mix => dw11k_contents(runs[mix], states[mix]) for mix in DW11_MIXES)
occursin("could not be certified", String(take!(diagnostics))) && error("an instant was not certified") # hide
cell(mix, k, d, q) = (m = dw11_phase_content(mix, d, q); m === nothing ?
    @sprintf("%5.1f/  – ", contents[mix][k][q]) : @sprintf("%5.1f/%4.1f", contents[mix][k][q], m))
for mix in ("OPC", "OPC-FA")
    println(mix, ", wt.% computed/measured")
    println("  days     C3S        C2S        C3A        C4AF     OPC reacted")
    for (k, d) in enumerate(days)
        println(lpad(d, 6), "  ", join((cell(mix, k, d, q) for q in ("C3S", "C2S", "C3A", "C4AF", "opc_reacted")), "  "))
    end
end
```

The last column is the table's own measure of the clinker reacted: one minus the
sum of the four phases over that sum at the mixing, the same for the model and
for the measurement.

## 4. Portlandite and ettringite

```@example ternary
for mix in DW11_MIXES
    println(mix, ", wt.% computed/measured")
    println("  days   portlandite   ettringite")
    for (k, d) in enumerate(days)
        println(lpad(d, 6), "     ", join((cell(mix, k, d, q) for q in ("portlandite", "ettringite")), "    "))
    end
end
```

Where the aluminum and the sulfate are at six months, in mmol per 100 g of
binder:

```@example ternary
phases = ("ettringite", "monosulphate12", "C4AH13", "monocarbonate", "hemicarbonate",
          "straetlingite", "C3AH6", "hydrotalcite")
println(rpad("", 16), join((lpad(mix, 10) for mix in DW11_MIXES)))
for name in phases
    println(rpad(name, 16), join((@sprintf("%10.2f", 1000 * dw11k_amount(states[mix][end], name)) for mix in DW11_MIXES)))
end
```

## 5. What the comparison says

**The clinker as a whole reacts at about the measured rate, its phases do not.**
The table's measure, the clinker reacted, is reproduced within five points from
seven days to six months in the four pastes, and seven to ten points low at one
day. It hides two errors in opposite directions. The law dissolves the alite too
slowly at first, 25 wt.% of the solids of the plain cement at one day where 15
are left, and then lets it stall, at 6.7 wt.% at six months where the
measurement finds 1.4; and it dissolves the belite from the first day, 8.9
against 17.1 wt.% at seven days, where the pastes leave it barely touched for a
month. The parameters are not this clinker's: fitting them on the plain paste
and testing them on the three others is the next step.

**Without fly ash the portlandite is the measured one.** In the plain and in the
limestone cement it is within two points of the measurement at every age, the
measurement being good to ±1 wt.%. The ettringite of the limestone cement is
within three points, and the limestone holds it as the authors found, as
monocarbonate takes the aluminum; in the plain cement it falls below the
measurement after a month, the minimization turning part of it into
monosulfate, which the authors also found, but more of it.

**With fly ash the portlandite and, without limestone, the ettringite are
lost.** The glass, at the degree of reaction the authors measured, consumes the
portlandite at equilibrium: at six months a fifth of the measured amount is
left in the paste without limestone, under half in the paste with it. This is
the finding of [the validation page](@ref ex-validation-blended) at measured
extents, which the kinetics carries in time: the C-S-H of this phase list,
CSHQ, takes no aluminum and stays at the calcium-to-silicon ratio of a gel in
equilibrium with portlandite, where the authors measured a gel poorer in calcium
holding aluminum. Without limestone, the aluminum of the glass then goes to the
AFm phases and to hydrogarnet, and the ettringite is gone from 28 days, where
the paste keeps 7 wt.%; with limestone it goes to monocarbonate, and the
ettringite stays. The two other gels the package ships do not fix it: at
measured extents `CASH+NK` behaves as `CSHQ`, and `CNASH_ss` takes the aluminum
at a calcium-to-silicon ratio a Portland cement's gel does not have
([the validation page](@ref sec-validation-blended-gels)).

## Where to go next

[Validation against measured blended pastes](@ref ex-validation-blended)
separates the thermodynamics from the kinetics on the same pastes, and
[Kinetics under partial equilibrium](@ref sec-theory-pe-kinetics) writes the
formulation these runs integrate.
