# [Calculation of thermodynamic properties of calcite dissolution](@id ex-from-scratch)

!!! info "Before this page"
    [Getting started](@ref sec-quickstart), whose calcite example this page
    rebuilds without a database.

This example is equivalent to quickstart except that no database is needed.

Let us recall the example which consisted of studying the dissolution of calcite in water.

$\ce{CaCO3 <=> Ca^2+ + CO3^2-}$


It is possible to calculate the thermodynamic properties of the reaction, in particular the solubility constant of the reaction ($\ln K$) which is related to the Gibbs free energy of the reaction ($\Delta_r G^°$). This solubility constant is a function of temperature and the calculation is performed at a reference temperature of 298 K and at a pressure of 1 Atm using the following equation:

$\Delta_r G^° = - RT\ln K$

where $\Delta_r G^°$ is deduced from the Gibbs energies of formation ($\Delta_f {G_i}^°$) of the other chemical species involved in the reaction:

$\Delta_r G^° = \sum_i \nu_i \Delta_f {G_i}^°$

### First step: species declaration

The first step, therefore, is to construct each of the species present in the reaction. This can be done with the [`Species`](@ref) function. Its simplest use is as follows:

```@example example1
using ChemistryLab

calcite = Species("CaCO3", aggregate_state=AS_CRYSTAL, class=SC_COMPONENT)
Ca²⁺ = Species("Ca+2", aggregate_state=AS_AQUEOUS, class=SC_COMPONENT)
CO₃²⁻ = Species("CO3-2", aggregate_state=AS_AQUEOUS, class=SC_COMPONENT)
CO₃²⁻
```

The created object contains a certain amount of information whose properties can be entered *a posteriori* (or during the construction of the [`Species`](@ref)).

!!! note "Calculation of Molar Mass"
    It can be noted that during the construction of the species, a calculation of the molar mass is systematically performed.

### Second step: calculation of the thermodynamic properties of each species

For each species, it is possible to assign thermodynamic properties, such as the Gibbs energy of formation or the heat capacity. These data come from a database, here [Thermoddem](https://thermoddem.brgm.fr) [Blanc2012](@cite), whose values for the three species at 25 °C are transcribed in `data/literature/Blanc2012.json`:

```@example example1
using DynamicQuantities, Printf

thermoddem(sp) = literature_row("Blanc2012", "individual_properties", sp)
@printf("%-8s %12s %12s %12s %12s %12s\n", "", "ΔfG° kJ/mol", "ΔfH° kJ/mol", "S° J/K/mol", "Cp° J/K/mol", "V° cm³/mol")
for sp in ("Calcite", "Ca+2", "CO3-2")
    r = thermoddem(sp)
    @printf("%-8s %12.3f %12.3f %12.2f %12.2f %12.2f\n", sp, ustrip(us"kJ/mol", r.dfG), ustrip(us"kJ/mol", r.dfH),
            ustrip(us"J/(mol*K)", r.S), ustrip(us"J/(mol*K)", r.Cp), ustrip(us"cm^3/mol", r.V))
end
```

#### Thermodynamic properties of formation

The first step involves associating the values of the thermodynamic properties of formation for each species. The values are read from `data/literature/Blanc2012.json` rather than typed, so that the page cannot drift from its source. Typing them by hand, as `83.47u"J/K/mol"` for the heat capacity, gives the same entry. For calcite:

```@example example1
cal = thermoddem("Calcite")
th_prop_0_calcite = Dict(:Cp⁰ => cal.Cp, :ΔₐH⁰ => cal.dfH, :S⁰ => cal.S, :ΔₐG⁰ => cal.dfG, :V⁰ => cal.V)
```

#### Heat capacity, enthalpy and free energy as a function of temperature

The second step is to describe the evolution of heat capacity as a function of temperature for each species. For calcite, the database expresses the heat capacity as a function of temperature, $C_p = a + bT + cT^{-2}$, whose coefficients are transcribed in the same file. A reference temperature can then be defined in order to construct thermodynamic functions, such as heat capacity, entropy, enthalpy and free enthalpy. For calcite, this can be done as follows:

```@example example1
mk = literature_row("Blanc2012", "maier_kelley", "Calcite")
params_Cp_calcite = Dict(:a₀ => mk.a, :a₁ => mk.b, :a₂ => mk.c)
T_ref = Dict(:T => 298.15u"K")
params_calcite = merge(th_prop_0_calcite, params_Cp_calcite, T_ref)
dtf_calcite = build_thermo_functions(:cp_ft_equation, params_calcite)
```

!!! tip "Expression of Cp as a function of temperature"
    In the Thermoddem database, the expression for heat capacity as a function of temperature is written as:

    $C_p(T) = a + b * T + c * T^{-2}$

    However, the expression given in Thermoddem for this species is a simplification of the following more complete function:
    
    $a_0 + a_1 * T + a_2 * T^{-2} + a_3 * T^{-0.5} + a_4 * T^2 + a_5 * T^3 + a_6 * T^4 + a_7 * T^{-3} + a_8 * T^{-1} + a_9 * T^{0.5} + a_{10} * log(T)$

    This expression is implemented in ChemistryLab and can be called by passing `:cp_ft_equation` as an argument.

Calling the `thermo_functions` (here `build_thermo_functions`) allows the calculation of the expressions for the different thermodynamic properties as a function of temperature, according to the following expressions:

$\Delta_a {H^°}_T = \int_{T_{ref}}^T C_p(\tau) d\tau + \Delta_f {H^°}$

${S^°}_T = \int_{T_{ref}}^T \frac{C_p(\tau)}{\tau} d\tau + {S^°}$

$\Delta_a {G^°}_T = \int_{T_{ref}}^T C_p(\tau) d\tau - T * \int_{T_{ref}}^T \frac{C_p(\tau)}{\tau} d\tau - (T - T_{ref}){S^°}_{T_{ref}} + \Delta_f G^°$

where $\Delta_a {H^°}_T$ and $\Delta_a {G^°}_T$ are the apparent enthalpy and free energy (Gibbs) at T.

The expressions for the thermodynamic properties of calcite can be added to the species `calcite` as follows. The symbolic expression is computed for each property, and the last line shows the one of the Gibbs energy of calcite:

```@example example1
calcite.Cp⁰ = dtf_calcite[:Cp⁰]
calcite.ΔₐH⁰ = dtf_calcite[:ΔₐH⁰]
calcite.S⁰ = dtf_calcite[:S⁰]
calcite.ΔₐG⁰ = dtf_calcite[:ΔₐG⁰]
```

```@example example1
using Plots

p1 = plot(xlabel="Temperature [°C]", ylabel="ΔₐG⁰ [J.mol⁻¹]", title="Gibbs energy of calcite \nas a function of temperature")
plot!(p1, θ -> calcite.ΔₐG⁰(T = 273.15+θ), 0:0.1:100, label="ΔₐG⁰ of calcite")
```


Similarly, we can provide the properties of the species $\ce{Ca^2+}$ and $\ce{CO3^2-}$, printed above. Unlike that of calcite, the heat capacity of an aqueous ion follows the Helgeson-Kirkham-Flowers equations of state (HKF [Helgeson1981](@cite)). ChemistryLab has them (`:solute_hkf88_reaktoro`), but they take the HKF coefficients of each ion, which are not among the values transcribed here; we therefore hold the heat capacity at its value at 25 °C, -26.38 and -276.88 J mol⁻¹ K⁻¹ respectively.


```@example example1
ca = thermoddem("Ca+2")
th_prop_0_Ca²⁺ = Dict(:Cp⁰ => ca.Cp, :ΔₐH⁰ => ca.dfH, :S⁰ => ca.S, :ΔₐG⁰ => ca.dfG, :V⁰ => ca.V)
params_Cp_Ca²⁺ = Dict(:a₀ => ca.Cp)   # the heat capacity at 25 °C, held constant
params_Ca²⁺ = merge(th_prop_0_Ca²⁺, params_Cp_Ca²⁺, T_ref)
dtf_Ca²⁺ = build_thermo_functions(:cp_ft_equation, params_Ca²⁺)
Ca²⁺.ΔₐG⁰ = dtf_Ca²⁺[:ΔₐG⁰]

co3 = thermoddem("CO3-2")
th_prop_0_CO₃²⁻ = Dict(:Cp⁰ => co3.Cp, :ΔₐH⁰ => co3.dfH, :S⁰ => co3.S, :ΔₐG⁰ => co3.dfG, :V⁰ => co3.V)
params_Cp_CO₃²⁻ = Dict(:a₀ => co3.Cp)
params_CO₃²⁻ = merge(th_prop_0_CO₃²⁻, params_Cp_CO₃²⁻, T_ref)
dtf_CO₃²⁻ = build_thermo_functions(:cp_ft_equation, params_CO₃²⁻)
CO₃²⁻.ΔₐG⁰ = dtf_CO₃²⁻[:ΔₐG⁰]
```

!!! warning "A constant heat capacity for the ions"
    The expressions for enthalpies and free energies remain temperature-dependent thanks to the integration performed on Cp. Within the temperature and pressure ranges typically considered in ChemistryLab (0–100 °C, 1 atm), assuming a constant Cp has little impact on the solubility product.

### Third step: writing the reaction

Now we can write the dissolution/precipitation reaction of calcite.

```@example example1
r = Reaction([calcite, Ca²⁺, CO₃²⁻]; equal_sign='↔')
```

During the construction of this reaction, the thermodynamic properties of the reaction are calculated.

$RT \; ln(K) = - \Delta_r G^° = - \sum_i \nu_i  \Delta_f {G^°}_i$

```@example example1
using Plots

p1 = plot(xlabel="Temperature [°C]", ylabel="pKs", title="Solubility product (pKs) of calcite \nas a function of temperature")
plot!(p1, θ -> r.ΔᵣG⁰(T = 273.15+θ) / (R_GAS * (273.15+θ)) / log(10), 0:0.1:100, label="pKs")
```

---

## Notes and next steps

This example demonstrates the **manual** workflow: create species, attach thermodynamic data from an external source, build a reaction and evaluate its temperature-dependent properties.

In practice, loading species from a database (see [Importing thermodynamic databases](@ref sec-importing-databases)) is faster and less error-prone:

```julia
all_species = build_species(datapath("cemdata18-thermofun.json"))
species = speciation(all_species, split("Cal H2O@");
              aggregate_state=[AS_AQUEOUS], exclude_species=split("H2@ O2@ CH4@"))
```

## Where to go next

The equilibrium state of the same system, with its species amounts, its pH and
the saturation index of calcite, follows from a `ChemicalSystem` built on these
species and passed to `equilibrate`; the tutorial
[Chemical Equilibrium](@ref sec-equilibrium) takes that workflow from end to end.