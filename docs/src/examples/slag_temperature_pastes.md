# [CEM I 52.5 R with slag and limestone, integrated in time at 5, 20 and 40 °C](@id ex-slag-temperature-pastes)

!!! info "Before this page"
    [CEM I 52.5 R with slag and limestone at 5, 20 and 40 °C](@ref ex-slag-temperature),
    where the law of the slag used here is fitted on its degree of reaction.

The pastes of [Snellings2022](@citet), a CEM I 52.5 R with 40 % slag and 10 %
limestone, were followed at 5, 20 and 40 °C for six months: the degrees of
reaction of their clinker and slag, which [the previous page](@ref ex-slag-temperature)
compares with the rate laws, and what the hydration formed, the bound water and
the portlandite by thermogravimetry (their Fig. 8), the portlandite, the
ettringite, the hydrotalcite and the carboaluminates by X-ray diffraction
(their Fig. 10). Here the pastes at w/b 0.5 are integrated from the mixing at
each temperature, the clinker and the slag under their laws and everything else
at equilibrium, and what they form is compared with those measurements, which
nothing was fitted on.

## 1. The pastes and their laws

The three materials are built from [Snellings2022; Table 1](@cite) (the
templates `PC (Snellings 2022)`, `slag (Snellings 2022)` and
`limestone (Snellings 2022)` of `data/recipe_templates.toml`): the cement by its
Rietveld phases, the slag by its calcite, its quartz and its glass, found by
difference from its analysis, the limestone by its calcite, quartz and dolomite.
The quartz is inert. The four clinker phases dissolve under the law of
[ParrottKilloh1984](@citet) with the published constants and activation
energies, the glass of the slag under the Waller law [Waller1999](@cite) with
the constants fitted on the previous page for w/b 0.5, at the temperature of the
paste. Everything else is at equilibrium from the mixing
(`scripts/snellings2022_pastes.jl`).

ASSUMED: the minor oxides of the cement, those of its analysis that its
Rietveld phases do not hold, magnesia among them, enter the equilibrium at the
mixing; the slag's sulfur, which its analysis does not report, is left out.

```@example slag-pastes
using ChemistryLab, DynamicQuantities, OptimaSolver, Printf
using Logging # hide
include(joinpath(pkgdir(ChemistryLab), "scripts", "snellings2022_pastes.jl"))

setup = sn22p_setup()
θ = sn22_slag_fit(; per_wb = (:α_max,), τ₀ = 10.0).θ[0.5]
@printf("slag: τ = %.2f d at 20 °C, n = %.3f, Ea = %.1f kJ/mol, ceiling %.3f\n", θ.τ, θ.n, θ.Ea, θ.α_max)
```

## 2. The three runs

```@example slag-pastes
temperatures = (5.0, 20.0, 40.0)
days = [1, 2, 7, 28, 90, 180]
diagnostics = IOBuffer() # hide
runs, obs = with_logger(ConsoleLogger(diagnostics)) do # hide
runs = Dict{Float64, Any}(T => sn22p_run(setup, T, θ) for T in temperatures)
obs = Dict{Float64, Any}(T => sn22p_observables(runs[T], days) for T in temperatures)
runs, obs # hide
end # hide
messages = String(take!(diagnostics)) # hide
occursin("re-speciation failed", messages) && error("a re-speciation failed") # hide
occursin("could not be certified", messages) && error("an instant was not certified") # hide
for T in temperatures
    @printf("%2.0f °C: %s, %d steps\n", T, runs[T].sol.retcode, length(runs[T].sol.t))
end
```

Neither law reads the equilibrium partition, so the trajectories do not depend
on it; every number below is computed on the certified replay of each run
([`speciated_states`](@ref)), in g per 100 g of binder, the thermogravimetry
over the mass at 550 °C as [Snellings2022](@citet) report it.

## 3. Bound water and portlandite

```@example slag-pastes
row(key, name, technique, T) = join((@sprintf("%5.1f/%5.1f", getproperty(o, key), m) for (o, m) in zip(obs[T], sn22p_measured(name, technique, T))), " ")
println("computed/measured        ", join((lpad("$d d", 11) for d in days)))
for (key, name, technique) in ((:bound_water, "bound water", "TGA"), (:portlandite_tga, "portlandite", "TGA"), (:portlandite, "portlandite", "XRD"))
    for T in temperatures
        @printf("%-12s %-4s %2.0f °C  %s\n", name, technique, T, row(key, name, technique, T))
    end
end
```

**The bound water and the portlandite follow the measurement.** The bound water
is within two points of it from 28 days on at the three temperatures, and at
seven days at 5 and 20 °C. Before, it is high, at 5 °C by four and a half points
at one day, where the calculation holds all its ettringite from the first day
and the paste half of it (Section 4). The portlandite is within 1.4 points of
the thermogravimetry and 1.8 of the diffraction at every age and temperature.
Its trend with the temperature is not the measured one: the thermogravimetry
finds less portlandite at six months the hotter the paste, 6.9, 6.1 and 5.4 %,
where the calculation finds a little more, 6.4, 6.5 and 6.8 %.

## 4. Ettringite, hydrotalcite and the carboaluminates

```@example slag-pastes
for (key, name) in ((:ettringite, "ettringite"), (:hydrotalcite, "hydrotalcite"), (:hemicarbonate, "hemicarbonate"), (:monocarbonate, "monocarbonate"))
    for T in temperatures
        @printf("%-13s %2.0f °C  %s\n", name, T, row(key, name, "XRD", T))
    end
end
for T in temperatures
    measured = sn22p_measured("hemicarbonate", "XRD", T) .+ sn22p_measured("monocarbonate", "XRD", T)
    computed = [o.hemicarbonate + o.monocarbonate for o in obs[T]]
    @printf("both          %2.0f °C  %s\n", T, join((@sprintf("%5.1f/%5.1f", c, m) for (c, m) in zip(computed, measured)), " "))
end
```

The sulfate of the paste bounds its ettringite. All of it held in ettringite
makes, per 100 g of binder:

```@example slag-pastes
st = runs[20.0].kp.initial_state
n_S = sum(ustrip(us"mol", st.n[i]) * get(atoms(s), :S, 0) for (i, s) in enumerate(st.system.species))
ett = st.system.dict_species["ettringite"]
@printf("%.1f g of ettringite at most\n", n_S / get(atoms(ett), :S, 0) * ustrip(us"g/mol", ett[:M]))
```

**The ettringite is that bound, at every temperature.** From the first day all
the sulfate of the cement is in ettringite, and it stays there: below the
temperature where monosulfate takes over (see
[the hydrates of two cements from 0 to 60 °C](@ref ex-hydrates-temperature)),
the calcite turns the aluminum the
sulfate leaves into carboaluminates. The diffraction finds as much at 40 °C,
but more at 20 °C, 12 %, and more still at 5 °C, 14 %: more than the sulfate
of the cement can make. Either the slag gives sulfur, which its analysis does
not report and this calculation leaves out, or the diffraction counts more
ettringite than there is; the measurement alone does not tell which.

**The carboaluminate is monocarbonate here, hemicarbonate in the paste.** With
calcite in excess, monocarbonate is the stable carboaluminate of Cemdata18
[Lothenbach2019](@cite), and the calculation forms it alone; the diffraction
finds mostly hemicarbonate, which forms first and turns into monocarbonate
slowly, as [Snellings2022](@citet) observe. From 28 days on the calculation
forms more carboaluminate than the two measured together, the last three lines
of the table: at six months half as much again at 5 °C and three times as much
at 40 °C, where the measured ettringite holds aluminum the calculation puts in
the AFm.

**The hydrotalcite is there from the first day.** The magnesia of the cement's
minor oxides and of the dolomite, at equilibrium from the mixing, forms it at
once, where the diffraction finds none before seven days at 20 °C and before
28 days at 5 °C; [Snellings2022](@citet) tie its formation to the reaction of
the slag. At three and six months the calculation is within 2.3 points of the
measurement.

## 5. What the comparison says

The rate laws of the previous page carry the temperature into the bound water
and the portlandite of the pastes within two points, from four weeks on; the
hydrates that hold the aluminum and the sulfate do not follow as well. The
ettringite of the calculation is fixed by the sulfate of the cement where the
measurement finds up to half as much again, the carboaluminate is the stable
monocarbonate where the paste holds hemicarbonate, and the hydrotalcite forms
from magnesia the calculation releases at the mixing. Each points to an input
that [Snellings2022](@citet) do not give, the sulfur of the slag and where the
magnesia of the cement sits, or to a hydrate the paste has not yet reached.

## Where to go next

[CEM I 52.5 N HTS and CEM II/A-L 42.5 R from 0 to 60 °C](@ref ex-hydrates-temperature)
follows the same hydrates at equilibrium up to the temperature where
monosulfate replaces ettringite and monocarbonate.
