# =============================================================================
#  snellings2022_kinetics.jl — a slag-limestone cement at 5, 20 and 40 °C
#
#  The ternary cement of Snellings et al. (2022), 50 % of a CEM I 52.5 R, 40 %
#  blast-furnace slag and 10 % limestone, at w/b 0.4, 0.5 and 0.6, cured at 5,
#  20 and 40 °C: the degrees of reaction of the clinker and of the slag the
#  authors measured by X-ray diffraction (their Fig. 6), against the laws of the
#  package at the temperature of each paste. The clinker reacts under the law of
#  Parrott and Killoh with its published constants and activation energies,
#  nothing fitted. No constant of the Waller law is published for a slag: its
#  time, its exponent and its activation energy are fitted here on the slag of
#  each w/b, and the activation energy compared with the one the authors fit
#  with another law. The page docs/src/examples/slag_temperature.md runs it.
# =============================================================================

using ChemistryLab, DynamicQuantities
using ForwardDiff: ForwardDiff
using LinearAlgebra: Diagonal, diag

const SN22 = "Snellings2022"
sn22_value(name) = ustrip(literature_value(SN22, name))

# The four clinker phases of Table 1 under the names of the rate laws, the two
# polymorphs of C3A together.
const SN22_CLINKER = ("C3S" => ("Alite",), "C2S" => ("Belite",), "C3A" => ("C3A ortho", "C3A cubic"), "C4AF" => ("C4AF",))
const SN22_PK = Dict(
    "C3S" => PK84_PARAMS_C3S, "C2S" => PK84_PARAMS_C2S,
    "C3A" => PK84_PARAMS_C3A, "C4AF" => PK84_PARAMS_C4AF,
)

# 0 °C is 273.15 K, by the definition of the Celsius scale.
_sn22_kelvin(T_C) = (T_C + 273.15) * u"K"

"""
    sn22_degrees(constituent) -> NamedTuple

The degrees of reaction of Fig. 6, `constituent` "clinker" or "slag": the
columns `temperature_C`, `w_b`, `age` (days) and `degree_percent`, as plain
numbers, one row per marker.
"""
function sn22_degrees(constituent)
    t = literature_table(SN22, "degree_of_reaction"; constituent)
    return (
        temperature_C = Float64.(ustrip.(t.temperature_C)), w_b = Float64.(ustrip.(t.w_b)),
        age = Float64.(ustrip.(us"d", t.age)), degree_percent = Float64.(ustrip.(t.degree_percent)),
    )
end

"""
    sn22_clinker_masses() -> Dict

The four clinker phases of the Portland cement, percent of its mass (Table 1,
XRD-Rietveld), the two polymorphs of C3A added.
"""
function sn22_clinker_masses()
    t = literature_table(SN22, "phase_composition"; material = "PC")
    pct(name) = ustrip(only(t.percent[t.phase .== name]))
    return Dict(ph => sum(pct, names) for (ph, names) in SN22_CLINKER)
end

"""
    sn22_clinker_degree(T_C, w_b, days) -> Vector

The degree of reaction of the clinker, percent, at `days`, at `T_C` °C and
water/binder `w_b`: each phase under the law of Parrott and Killoh with the
constants and activation energies of Lavergne et al. (2018), at the Blaine
fineness of the cement and its water/cement ratio, and the four weighted by
their masses, as the authors weight them (one minus the sum of the phases left
over their sum at the mixing).
"""
function sn22_clinker_degree(T_C, w_b, days)
    m = sn22_clinker_masses()
    blaine = literature_value(SN22, "blaine_pc")
    w_c = w_b / (sn22_value("pc_percent") / 100)
    T = _sn22_kelvin(T_C)
    e = Dict(ph => ParrottKillohExtent(ph; T, blaine, w_c) for ph in keys(m))
    total = sum(values(m))
    return [100 * sum(m[ph] * extent(e[ph], d) for ph in keys(m)) / total for d in days]
end

"""
    sn22_waller_degree(θ, T_C, days) -> Vector

The degree of reaction of the slag, percent, under the Waller law at the
constant temperature `T_C` °C: `θ = (τ, n, Ea, α_max)`, `τ` in days at the
reference temperature of the law, `Ea` in kJ/mol, `α_max` the ceiling of
[`waller`](@ref). At a constant temperature that law integrates to
`α = α_max/(1 + (α_max τ/(A t))ⁿ)`, `A` the Arrhenius factor from the
reference temperature: the rate is `A` times the one at the reference, so the
time is scaled by `A` (the test checks the identity on the law itself).
"""
function sn22_waller_degree(θ, T_C, days)
    τ, n, Ea, α_max = θ
    T_ref = ustrip(us"K", literature_value("Lavergne2018", "T_ref"))
    T = ustrip(us"K", _sn22_kelvin(T_C))
    A = exp(-Ea * 1000 / R_GAS * (1 / T - 1 / T_ref))
    return [100 * α_max / (1 + (α_max * τ / (A * t))^n) for t in days]
end

# The constants the fit starts from: the exponent and the activation energy of
# the fly-ash set of Lavergne et al. (2018), the only ones the law has, and no
# ceiling.
_sn22_waller_start(τ₀) = (
    τ = float(τ₀), n = float(WALLER_PARAMS_FLY_ASH.n),
    Ea = ustrip(us"kJ/mol", WALLER_PARAMS_FLY_ASH.Ea), α_max = 1.0,
)

"""
    sn22_slag_fit(w_b = (0.4, 0.5, 0.6); shared = (:τ, :n, :Ea), per_wb = (), τ₀ = 30.0)
        -> NamedTuple

The Waller law of the slag fitted on its degrees of reaction at the w/b of
`w_b` (Fig. 6b: three temperatures, six ages each), the constants of `shared`
common to them, those of `per_wb` one per w/b, the others those the fit starts
from (the fly-ash set, no ceiling): Levenberg–Marquardt on their logarithms
with the exact Jacobian (`ForwardDiff`), from `τ₀` days.

Returns `θ` (a dictionary from each w/b to its four constants at the fit), the
names and values of the free constants (`names`, `values`), the
root-mean-square misfit in points of degree (`rms`), the model and the
measurement at each marker with its temperature, w/b and age, and the
[`identifiability`](@ref) of the free constants at the fit for a measurement
good to the 10 points the authors give.
"""
function sn22_slag_fit(w_b = (0.4, 0.5, 0.6); shared = (:τ, :n, :Ea), per_wb = (), τ₀ = 30.0)
    w_bs = Tuple(w_b)
    d = sn22_degrees("slag")
    rows = findall(in(w_bs), d.w_b)
    measured = d.degree_percent[rows]
    start = _sn22_waller_start(τ₀)
    names = vcat([string(k) for k in shared], [string(k, " (w/b ", x, ")") for x in w_bs for k in per_wb])
    function constants(θ, x)
        q = merge(start, NamedTuple{shared}(Tuple(θ[1:length(shared)])))
        off = length(shared) + (findfirst(==(x), w_bs) - 1) * length(per_wb)
        q = merge(q, NamedTuple{per_wb}(Tuple(θ[(off + 1):(off + length(per_wb))])))
        return (q.τ, q.n, q.Ea, q.α_max)
    end
    model(θ) = [only(sn22_waller_degree(constants(θ, d.w_b[i]), d.temperature_C[i], (d.age[i],))) for i in rows]
    residual(z) = model(exp.(z)) .- measured
    z = log.(vcat([getproperty(start, k) for k in shared], [getproperty(start, k) for _ in w_bs for k in per_wb]))
    r = residual(z)
    f = sum(abs2, r)
    λ = 1.0e-2
    for _ in 1:300
        J = ForwardDiff.jacobian(residual, z)
        A, gr = J' * J, J' * r
        accepted = false
        for _ in 1:25
            zn = z .- (A + λ * Diagonal(diag(A) .+ 1.0e-12)) \ gr
            rn = residual(zn)
            fn = sum(abs2, rn)
            if fn < f
                z, r, f, accepted = zn, rn, fn, true
                λ = max(λ / 3, 1.0e-9)
                break
            end
            λ *= 4
        end
        accepted || break
    end
    θ = exp.(z)
    id = identifiability(model, θ; observed = measured, noise = sn22_value("slag_degree_precision"), names)
    return (;
        θ = Dict(x => NamedTuple{(:τ, :n, :Ea, :α_max)}(constants(θ, x)) for x in w_bs),
        names, values = θ, rms = sqrt(f / length(rows)), model = model(θ), measured,
        temperature_C = d.temperature_C[rows], w_b = d.w_b[rows], age = d.age[rows], identifiability = id,
    )
end

"""
    sn22_slag_law(θ; name = "slag") -> KineticFunc

The Waller law of the package ([`waller`](@ref)) with the four constants of
`θ`, one entry of the `θ` of [`sn22_slag_fit`](@ref) (`τ` in days at the
reference temperature, `Ea` in kJ/mol): the law whose isothermal integral
[`sn22_waller_degree`](@ref) writes in closed form.
"""
sn22_slag_law(θ; name = "slag") = waller(
    merge(WALLER_PARAMS_FLY_ASH, (τ = θ.τ * u"d", n = θ.n, Ea = θ.Ea * u"kJ/mol")), name; α_max = θ.α_max,
)

"""
    sn22_durdzinski(θ) -> NamedTuple

The law of `θ` (one entry of the `θ` of [`sn22_slag_fit`](@ref)) against the
two slags of the RILEM round robin (Durdziński et al. 2017, Table 5): the
degrees of reaction measured by image analysis of electron micrographs
(laboratories B and E) on the sealed pastes, 40 % slag at w/b 0.4, the w/b of
the `θ` to give. ASSUMED: the pastes cured at 20 °C, a temperature the round
robin does not state.

Returns the columns `material`, `lab`, `age` (days), `measured` and `model`
(percent).
"""
function sn22_durdzinski(θ)
    t = literature_table("Durdzinski2017", "degree_of_reaction"; technique = "SEM-IA", curing = "sealed")
    keep = findall(m -> m in ("S1", "S2"), t.material)
    age = Float64.(ustrip.(t.age_days[keep]))
    T_C = ustrip(us"K", literature_value("Lavergne2018", "T_ref")) - 273.15
    model = sn22_waller_degree((θ.τ, θ.n, θ.Ea, θ.α_max), T_C, age)
    return (
        material = t.material[keep], lab = t.lab[keep], age,
        measured = Float64.(ustrip.(t.degree_percent[keep])), model,
    )
end
