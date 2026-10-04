# [CEM I 52.5 N and slag pastes, integrated in time](@id ex-blended-slag-kinetics)

!!! info "The same pastes"
    [CEM I 52.5 N and slag in an isothermal calorimeter, read off the states](@ref sec-example-isothermal)
    and [the bound water of the same pastes](@ref sec-example-tga) read them at
    the degrees of hydration the image analysis measured;
    [the clinker kinetics page](@ref sec-clinker-kinetics) introduces the
    Parrott–Killoh law.

The two pages before this one take the degrees of hydration of the pastes of
[Gruyaert2010](@citet) as measured and compute the states they lead to. Here the
degrees are not given: the pastes are integrated from the mixing, the four
clinker phases dissolving under the Parrott–Killoh law, the slag glass under the
Waller law, and the equilibrium solved along the way. The comparison with the
article is then a test of the rate laws, and only one number is fitted to it.

The pastes, the system and the activity model are those of the two pages before,
built by `scripts/gruyaert2010.jl`; `scripts/gruyaert2010_kinetics.jl` adds the
slag glass as a species ([`glass_species`](@ref)), the recipe of the unreacted
mix, the rate laws and the run.

## 1. The slag's characteristic time, from one measurement

The Waller law gives a degree of reaction that follows a sigmoid in the logarithm
of time,

```math
\alpha(t) = \frac{1}{1 + (\tau/t)^{n}} ,
```

at the reference temperature and fineness of the law, which are 20 °C and the
slag's Blaine fineness here. No parameter set for a slag is shipped: the one the
package attributed to [Waller1999](@citet) is not in the thesis. The exponent
`n` is taken from the fly-ash set of [Lavergne2018](@citet), an assumption, and
`τ` is solved from one measurement, the 72 % of slag the image analysis found at
28 months in the paste with 50 % slag.

```@example g10k
using ChemistryLab, DynamicQuantities, OptimaSolver, Printf
using Logging # hide
include(joinpath(pkgdir(ChemistryLab), "scripts", "gruyaert2010_kinetics.jl"))

G10(table, column; where...) = only(getproperty(literature_table("Gruyaert2010", table; where...), column))
αs_measured(days, sb) = G10("hydration_degree_slag", :alpha_slag; age_days = days * u"d", slag_to_binder = sb) / 100
τ = gruyaert_slag_time()
α_waller(t) = 1 / (1 + (ustrip(us"d", τ) / t)^WALLER_PARAMS_FLY_ASH.n)
@printf("τ = %.0f days\n", ustrip(us"d", τ))
for (days, sb) in ((852, 0.5), (2, 0.85), (852, 0.85))
    @printf("slag %.2f, %3d days: law %4.2f, measured %4.2f\n", sb, days, α_waller(days), αs_measured(days, sb))
end
```

The first line is the calibration and reproduces its point. The second is a check
it was not set on, the paste with 85 % slag at 2 days, and the law meets it. The
third it misses by a factor of nearly two: in that paste the cement is a sixth of
the binder, and the portlandite that activates a slag runs out, which a law
written in the slag's own degree of reaction cannot see.

## 2. The cement alone

The plain paste, integrated over the 1018 days of the last bound-water
measurement. The rate law of Parrott and Killoh carries its own water/cement
factor ([`pk_wc_factor`](@ref)), which slows a phase once its degree passes
`1.333 w/c`, here 0.67.

```@example g10k
cs0 = gruyaert_system()
diagnostics = IOBuffer() # hide
plain = with_logger(ConsoleLogger(diagnostics)) do # hide
plain = gruyaert_run(cs0; slag = 0.0)
end # hide
occursin("re-speciation failed", String(take!(diagnostics))) && error("a re-speciation failed") # hide
clinker = ["C3S", "C2S", "C3A", "C4AF"]
day = 86400.0
αc_measured(days, sb) = G10("hydration_degree", :alpha_cement; age_days = days * u"d", slag_to_binder = sb) / 100
@printf("%d accepted steps, retcode %s\n", length(plain.sol.t), plain.sol.retcode)
for days in (2, 852)
    @printf("cement at %3d days: computed %4.2f, image analysis %4.2f\n", days,
            gruyaert_degree(plain, clinker, days * day), αc_measured(days, 0.0))
end
```

Neither law reads the equilibrium partition, so the trajectory does not depend
on it: the run solves it at each accepted step only to report it, with the
interior point. Every number below is computed on the certified replay of the
run ([`speciated_states`](@ref)).

The heat is read off the certified replay of the run ([`heat_release`](@ref)),
from the unreacted mix, and compared with the reaction degree of Table 6 times
the total heat of Table 2.

```@example g10k
Q_total = ustrip(us"J/g", G10("isothermal_total_heat", :heat_infinity; slag_to_binder = 0.0))
r(days) = G10("reaction_degree", :r; age_days = days * u"d", slag_to_binder = 0.0) / 100
t, Q, _ = with_logger(ConsoleLogger(diagnostics)) do # hide
t, Q, _ = heat_release(plain.sol, plain.kp; times = [2.0, 7.0] .* day, reference = plain.kp.initial_state)
end # hide
occursin("could not be certified", String(take!(diagnostics))) && error("an instant was not certified") # hide
for (days, q) in zip((2, 7), Q)
    @printf("heat at %d days: computed %5.1f J/g, measured %5.1f J/g\n", days, q / 100, r(days) * Q_total)
end
```

Nothing in this was fitted. The heat at 2 days is the measured one to 0.1 %, the
degree of hydration 13 % below the image analysis; at 7 days the heat is 14 %
below. At 28 months the cement is at 85 % against 74 %: the water/cement factor
brings the law down from the 94 % it reaches without it, and the article's 74 %
is the ultimate degree it computes for this w/c.

## 3. The blends

With slag the cement has more water to hydrate in, `w/c = w/b / (1 − s)`, 1.0 for
50 % slag and 3.3 for 85 %, and the water/cement factor no longer acts.

```@example g10k
slag_sp = gruyaert_slag_species()
cs1 = gruyaert_kinetic_system(slag_sp)
blends = with_logger(ConsoleLogger(diagnostics)) do # hide
# `Any`: the type of a run is long enough that a dictionary specialized on it
# takes minutes to compile.
blends = Dict{Float64, Any}(sb => gruyaert_run(cs1; slag = sb, τ_slag = τ, slag_species = slag_sp) for sb in (0.5, 0.85))
end # hide
occursin("re-speciation failed", String(take!(diagnostics))) && error("a re-speciation failed") # hide
# Not every age was measured for every paste: a dash where Table 5 has nothing.
measured(f, days, sb) = try
    @sprintf("%4.2f", f(days, sb))
catch
    "  — "
end
println("slag  days   cement: computed  measured   slag: computed  measured")
for sb in (0.5, 0.85), days in (2, 852)
    @printf("%4.2f  %4d           %4.2f      %s           %4.2f      %s\n", sb, days,
            gruyaert_degree(blends[sb], clinker, days * day), measured(αc_measured, days, sb),
            gruyaert_degree(blends[sb], ["BFS"], days * day), measured(αs_measured, days, sb))
end
```

The cement of the blends reaches 94 % at 28 months, the image analysis 94 % and
91 %: more than in the plain paste, as measured, and for the reason the factor
encodes. Two numbers are missed, and they say what the laws lack. At 2 days the
cement of the paste with 85 % slag is at 29 %, the law at 47 %: the
Parrott–Killoh constants know nothing of the slag around the grains. And the slag
of that paste stops at 39 %, where the law takes it to 72 %, for the reason of
Section 1.

## 4. Bound water over 1018 days

The bound water of the replayed states against the thermogravimetry of Fig. 7,
batch b, whose materials are those of the calorimetry. The calculation counts
every hydrogen of the solids and the thermobalance what leaves above 105 °C once
the sample has been dried, so the computed water is the larger, as on
[the page of the same pastes at measured degrees](@ref sec-example-tga).

```@example g10k
wb_measured(sb) = (tb = literature_table("Gruyaert2010", "bound_water"; batch = "b", slag_to_binder = sb);
                   (ustrip.(us"d", tb.age_days), ustrip.(tb.bound_water)))
g(q) = ustrip(uconvert(us"g", q))
println("slag   days  computed  measured  ratio")
for (sb, run) in ((0.0, plain), (0.5, blends[0.5]), (0.85, blends[0.85]))
    ages, wb = wb_measured(sb)
    keep = [i for i in eachindex(ages) if ages[i] >= 1 && ages[i] <= 1019]
    states = with_logger(ConsoleLogger(diagnostics)) do # hide
    states = speciated_states(run.sol, run.kp; times = ages[keep] .* day)
    end # hide
    occursin("could not be certified", String(take!(diagnostics))) && error("an instant was not certified") # hide
    for (k, st) in zip(keep, states)
        c = g(bound_water(st))
        @printf("%4.2f  %5.0f   %6.2f    %6.2f   %4.2f\n", sb, ages[k], c, wb[k], c / wb[k])
    end
end
```

For the plain paste and the paste with 50 % slag the ratio stays between 1.1 and
1.8 over the 1018 days, the gap of definition the page at measured degrees finds
at its four points. The paste with 85 % slag is different: the ratio is 1.03 to
1.13 up to a month, then grows to 2.6, the measured water holding near 10 g
while the computed one keeps rising with a slag the law goes on reacting. It is
the miss of Section 1, read on a second observable.

## What these runs do not do

- **The heat of a blend.** A glass has no enthalpy of formation in any database,
  so the heat of the blends is not computed; the calorimetry page derives the
  enthalpy the glass would need, and fed back here it would reproduce the
  measured heat by construction.
- **Temperature.** The slag's rate carries the apparent activation energy the
  article fits for it at the paste's cement-to-binder ratio (its Eq. 3), and the
  clinker the energies of [Lavergne2018](@citet); at 20 °C, the reference
  temperature of both laws, neither acts.
- **The activation of the slag by portlandite**, and the early slowdown of the
  cement in a paste that is mostly slag: no published law of either is wired, and
  Sections 1 and 3 show where their absence is measured.
