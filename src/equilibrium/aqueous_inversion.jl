# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

# ── the aqueous solutes from their potentials, through the ionic strength ─────
#
# The dual Newton of OptimaSolver recovers the amounts of a phase from the
# potentials its members have to meet, `hᵢ(x) = cᵢ`. Its default does so by
# sweeps, `wᵢ ← wᵢ + cᵢ − hᵢ`, which assumes `∂hᵢ/∂wᵢ = 1` and nothing else. The
# Debye–Hückel family couples the solutes through the ionic strength: `∂h/∂w`
# has a part of rank one, and where it is strong the sweeps cycle instead of
# converging. Measured on the coupled hydration of a CEM I over three hours, 94 %
# of 790 000 inversions ended unconverged; on the cement pastes of a thesis, 84 %
# where the model has a solution, and every one where the limiting law has none.
#
# For a model whose coefficients depend on the composition through the ionic
# strength alone, the solutes are explicit at a given `I`,
#
#     ln mᵢ(I) = cᵢ − ln γᵢ(I),
#
# and `I = ½ Σ zᵢ² mᵢ(I)` is one equation in one unknown. That is how PHREEQC
# treats the ionic strength, as an unknown of its own [ParkhurstAppelo2013](@cite).
# It is solved here exactly, on `s = ln I`, or found to have no root, which is the
# runaway of a model past its range said in so many words. Measured on the same
# inversions, the root reproduces `hᵢ = cᵢ` to 1e-12 where the sweeps did not
# converge, and the cases the sweeps cycled on without a solution have no root.

"""
    _aqueous_form(model, cs, members) -> Union{Nothing, NamedTuple}

What the inversion of the aqueous phase needs from an activity model: the charge
and the `log₁₀ γ(I)` of each member, the molar mass of the solvent and the
Debye–Hückel parameters at given conditions. `nothing` for a model whose
coefficients depend on more than the ionic strength (Pitzer, SIT), which the
sweeps of OptimaSolver then handle.
"""
_aqueous_form(model, cs, members) = nothing

function _aqueous_form(::DiluteSolutionModel, cs::ChemicalSystem, members)
    M_w = ustrip(us"kg/mol", cs.species[only(cs.idx_solvent)][:M])
    return (; kind = :dilute, ln_c_solvent = log(1.0 / M_w))
end

function _ionic_form(cs, members, log10γ, AB)
    jw = only(cs.idx_solvent)
    M_w = ustrip(us"kg/mol", cs.species[jw][:M])
    z = [Int(charge(cs.species[i])) for i in members]
    return (; kind = :ionic, z, M_w, log10γ, AB)
end

function _aqueous_form(model::HKFActivityModel, cs::ChemicalSystem, members)
    å = _promoted([iszero(charge(cs.species[i])) ? 0.0 : _hkf_lookup_å(cs.species[i], model) for i in members])
    Kn = _promoted([iszero(charge(cs.species[i])) ? _setschenow(cs.species[i], model) : 0.0 for i in members])
    log10γ = (t, z, I, sqrtI, A, B) -> iszero(z) ? _log10γ_neutral(model, I, Kn[t]) :
        _log10γ_ion(model, z, å[t], I, sqrtI, A, B)
    AB = p -> (model.temperature_dependent && hasproperty(p, :T) && hasproperty(p, :P)) ?
        (ab = hkf_debye_huckel_params(p.T, p.P); (ab.A, ab.B)) : (model.A, model.B)
    return _ionic_form(cs, members, log10γ, AB)
end

function _aqueous_form(model::DaviesActivityModel, cs::ChemicalSystem, members)
    log10γ = (t, z, I, sqrtI, A, B) -> iszero(z) ? _log10γ_neutral(model, I) :
        _log10γ_ion(model, z, 0.0, I, sqrtI, A, 0.0)
    AB = p -> (model.temperature_dependent && hasproperty(p, :T) && hasproperty(p, :P)) ?
        (hkf_debye_huckel_params(p.T, p.P).A, 0.0) : (model.A, 0.0)
    return _ionic_form(cs, members, log10γ, AB)
end

function _aqueous_form(model::TruesdellJonesActivityModel, cs::ChemicalSystem, members)
    par = [get(model.parameters, symbol(cs.species[i]), nothing) for i in members]
    log10γ = (t, z, I, sqrtI, A, B) -> _truesdell_jones(par[t], z, I, sqrtI, A, B)
    AB = p -> (model.temperature_dependent && hasproperty(p, :T) && hasproperty(p, :P)) ?
        (ab = hkf_debye_huckel_params(p.T, p.P); (ab.A, ab.B)) : (model.A, model.B)
    return _ionic_form(cs, members, log10γ, AB)
end

"""
    _aqueous_inverter(des, pq = nothing) -> Union{Nothing, Function}

The `invert` of OptimaSolver's `SolutionPhase` for the aqueous phase of `des`:
`(c, ref, w, q, params) -> log-amounts`, or `nothing` from the call when no ionic
strength solves the system. `pq(q, params)` gives the parameters the activity
model sees when a constraint makes them unknowns (the temperature of an
adiabatic solve); without it they are `params`. For SIT and Pitzer, whose
coefficients depend on more than the ionic strength, it is Newton's method on the
solutes (`_newton_inverter`). Returns `nothing` itself for a model neither covers.
"""
function _aqueous_inverter(des::DualEquilibriumSolver, pq = nothing)
    form = _aqueous_form(des.model, des.system, des.idx_aq)
    form === nothing && return _newton_inverter(des, pq)
    jref = des.j_solvent
    # Built inside a scope that solves on values (`_STRIP_TAGS`), it returns
    # values, as the activity closure it stands for does (`_scoped_lna`).
    tags = _STRIP_TAGS[]
    return function (c, ref, w, q, params)
        p = pq === nothing ? params : pq(q, params)
        out = _invert_aqueous(form, c, ref, w, p, jref)
        return (out === nothing || isempty(tags)) ? out : _strip_tags(out, tags)
    end
end

function _invert_aqueous(form, c, ref, w, p, jref)
    if form.kind === :dilute
        lnw = log(ref)
        return [j == jref ? w[j] : c[j] + lnw - form.ln_c_solvent for j in eachindex(c)]
    end
    z = form.z
    A, B = form.AB(p)
    ϵ = p.ϵ
    ln10 = log(10.0)
    lndenom = log(ref * form.M_w)
    ions = [t for t in eachindex(z) if t != jref && z[t] != 0]
    lnγ(t, I) = ln10 * form.log10γ(t, z[t], I, sqrt(I + ϵ), A, B)
    # `ln(½ zₜ²) + cₜ`, the part of `ln(½ zₜ² mₜ)` that does not depend on `I`.
    a0 = [log(abs(z[t])^2 / 2) + c[t] for t in ions]
    # `ln I` of the composition the potentials give at an ionic strength `exp(s)`,
    # minus `s`: its root is the self-consistent ionic strength. One pass, a
    # running log-sum-exp, and nothing allocated: it is evaluated some twenty
    # times per inversion, an inversion per trial step of the outer Newton.
    F(s) = begin
        I = exp(s)
        M = a0[1] - lnγ(ions[1], I)
        acc = one(M)
        @inbounds for k in 2:length(ions)
            x = a0[k] - lnγ(ions[k], I)
            if x > M
                acc = acc * exp(M - x) + 1
                M = x
            else
                acc += exp(x - M)
            end
        end
        M + log(acc) - s
    end
    I = if isempty(ions)
        zero(eltype(c))
    else
        # Values for the search, and the derivative in `s` taken on `F` itself
        # before its values are read: stripping first would strip the derivative.
        Fv = s -> _plain(F(s))
        dFv = s -> _plain(ForwardDiff.derivative(F, s))
        # Below the dilute limit, where the solutes are those of zero ionic
        # strength: the Debye–Hückel coefficients only raise them from there, so
        # the root lies above.
        s_start = max(log(1.0e-30), Fv(log(1.0e-30)) + log(1.0e-30) - 1.0)
        sv = _ionic_strength_root(Fv, dFv, s_start)
        sv === nothing && return nothing
        # Found on the values. Whatever carries dual numbers (the potentials, the
        # temperature through A and B, a parameter of the model) shows in the
        # type of `F`; two Newton steps in that type then give the root its
        # derivatives, every order of a nested differentiation, by the
        # implicit-function theorem.
        f1 = F(sv)
        s = sv + zero(f1)
        if !(f1 isa AbstractFloat)
            for _ in 1:2
                s -= F(s) / ForwardDiff.derivative(F, s)
            end
        end
        exp(s)
    end
    return [j == jref ? w[j] : c[j] - lnγ(j, I) + lndenom for j in eachindex(c)]
end

"""
    _ionic_strength_root(F, dF, s_start; smax = log(1e4), hmax = 4) -> Union{Nothing, Float64}

The first root of `F(s) = ln S(eˢ) − s` above `s_start`, where `F` crosses from
positive to negative: the self-consistent ionic strength `eˢ` of a solution whose
solutes follow from it, on the branch connected to the dilute limit. `nothing`
when that branch holds none.

The first one, and not the one nearest the composition the solve holds. The
limiting law has up to three: past its range it can dip below zero over a
narrow interval only, two roots close together, and the coefficient `Ḃ I`
brings a third back far above, at thousands of mol/kg for a CEM I pore
solution with sodium chloride. The physical answer is the dilute branch; taking
the root nearest a start made the answer depend on the path, and a search
stepping up from a start sitting on the root itself, with `F` at `+2e-16`,
stepped over a dip `0.15` wide in `ln I`. The steps are Newton's from the left,
which a convex `F` keeps short of its first root, and a dip a step spans shows
as `F′` changing sign there and is found whatever its width.

A dip that stays above zero ends the dilute branch without a root, and so does
`smax`, 1e4 mol/kg, a molality no aqueous model is stated for: the potentials
then hold no solution on it, which is the answer. Searching on past such a dip
reached the third root, and from the potentials of a cold start the solve then
began at a composition three orders of magnitude off, which it never left: the
start of a chloride-loaded paste failed where recovering the solutes one by one
from that start had certified it at once. `s_start` lies below the root
(`_invert_aqueous`); a start already past it, with potentials that hold more
than `smax` at zero ionic strength, is searched downwards, under the same
ceiling.
"""
function _ionic_strength_root(F, dF, s_start; smax = log(1.0e4), hmax = 4.0)
    a = s_start
    Fa = F(a)
    isfinite(Fa) || return nothing
    if !(Fa > 0)                                     # already past it: below
        r = _bracketed_root(F, dF, a - 60.0, a)
        return r <= smax ? r : nothing
    end
    da = dF(a)
    while a < smax
        # Newton's step towards the root while `F` falls, which from the left of
        # a convex `F` never passes it; at most `hmax`, and that much where `F`
        # does not fall.
        step = da < 0 ? min(-Fa / da, hmax) : hmax
        b = min(a + step, smax)
        b - a <= 4 * eps(max(1.0, abs(a))) && return b
        Fb = F(b)
        isfinite(Fb) || return nothing
        Fb <= 0 && return _bracketed_root(F, dF, a, b)
        db = dF(b)
        if da < 0 && db > 0
            # A minimum inside: below zero, it holds the first root; above, the
            # dilute branch ends there. Located by the secant on `F′`, kept inside
            # its bracket (Illinois): a few evaluations, where under the limiting
            # law past its range nine searches in ten end on such a dip.
            m = _bracketed_zero(dF, a, da, b, db)
            return F(m) <= 0 ? _bracketed_root(F, dF, a, m) : nothing
        end
        a, Fa, da = b, Fb, db
    end
    return nothing
end

# A zero of `f` in `[a, b]`, `fa < 0 < fb`, by the Illinois variant of regula
# falsi: superlinear, and never outside the bracket.
function _bracketed_zero(f, a, fa, b, fb)
    side = 0
    m = (a * fb - b * fa) / (fb - fa)
    fm = f(m)
    k = 1
    while !(k == 60 || fm == 0 || b - a <= 1.0e-10 * max(1.0, abs(a)))
        if fm < 0
            a, fa = m, fm
            side == -1 && (fb /= 2)
            side = -1
        else
            b, fb = m, fm
            side == 1 && (fa /= 2)
            side = 1
        end
        m = (a * fb - b * fa) / (fb - fa)
        fm = f(m)
        k += 1
    end
    return m
end

# The root of `F` in `[a, b]`, `F(a) > 0 ≥ F(b)`, by Newton's method kept inside
# the bracket, `dF` its derivative.
function _bracketed_root(F, dF, a, b)
    lo, hi = a, b
    F(lo) > 0 || return lo
    s = (lo + hi) / 2
    sn, k, done = s, 0, false
    while !done
        k += 1
        fs = F(s)
        fs > 0 ? (lo = s) : (hi = s)
        sn = s - fs / dF(s)
        (isfinite(sn) && lo < sn < hi) || (sn = (lo + hi) / 2)
        done = k == 200 || abs(sn - s) <= 4 * eps(max(1.0, abs(s))) ||
            hi - lo <= 4 * eps(max(1.0, abs(lo)))
        s = sn
    end
    return sn
end

# ── models of more than the ionic strength: Newton on the solutes ─────────────
#
# SIT adds to the Debye–Hückel term of an ion a sum over its counter-ions,
# `Σₖ ε(i,k) mₖ`, and Pitzer adds pair and triplet terms: their coefficients
# depend on the composition through more than `I`, and no single equation gives
# the solutes back. Their potentials are met by Newton's method on the
# log-amounts, `hᵢ(w) = cᵢ`, with the model's own Jacobian, exact by forward mode:
# the identity, plus the rank-one Debye–Hückel term, plus the interactions, where
# OptimaSolver's sweeps assume the identity alone. The start is the better of the
# composition the solve holds and the one the model's Debye–Hückel part gives
# through the ionic strength (`_invert_aqueous`). The globalization is that of
# OptimaSolver's Newton on a mixing phase: a step of at most 30 in any
# log-amount, halved until the squared residual falls.
#
# A solute whose amount lies below the activity floor has an activity that no
# longer moves with it, so no Newton step can place it. It is held out of the
# iteration and placed afterwards from its activity coefficient at the
# composition found, which its own vanishing amount does not change.

"""
    _newton_predictor_form(model, cs, members) -> Union{Nothing, NamedTuple}

The Debye–Hückel part of a model the Newton inversion covers, as an
`_aqueous_form`, which gives that inversion its start; `nothing` for a model it
does not cover.
"""
_newton_predictor_form(model, cs, members) = nothing

function _newton_predictor_form(model::SITActivityModel, cs::ChemicalSystem, members)
    b = model.b
    log10γ = (t, z, I, sqrtI, A, B) -> iszero(z) ? zero(sqrtI) : -z^2 * A * sqrtI / (1 + b * sqrtI)
    AB = p -> (model.temperature_dependent && hasproperty(p, :T) && hasproperty(p, :P)) ?
        (hkf_debye_huckel_params(p.T, p.P).A, 0.0) : (model.A, 0.0)
    return _ionic_form(cs, members, log10γ, AB)
end

function _newton_predictor_form(model::PitzerActivityModel, cs::ChemicalSystem, members)
    bp = model.parameters.b
    # `ln γᵢ = zᵢ² f^γ`, with `f^γ` on the osmotic basis `A_φ = A ln 10 / 3`.
    fγ(A, sqrtI) = -(A / 3) * (sqrtI / (1 + bp * sqrtI) + 2 / bp * log1p(bp * sqrtI))
    log10γ = (t, z, I, sqrtI, A, B) -> z^2 * fγ(A, sqrtI)
    AB = p -> (model.temperature_dependent && hasproperty(p, :T) && hasproperty(p, :P)) ?
        (hkf_debye_huckel_params(p.T, p.P).A, 0.0) : (_DH_A_25C, 0.0)
    return _ionic_form(cs, members, log10γ, AB)
end

"""
    _newton_inverter(des, pq = nothing) -> Union{Nothing, Function}

The `invert` of the aqueous phase for SIT and Pitzer: the solutes recovered from
their potentials by Newton's method on their log-amounts, with the model's exact
Jacobian. `nothing` from the call when the iteration finds no composition, which
OptimaSolver reads as it reads the inversion through the ionic strength: an
iterate is then swept, and a trial of its line search passed over. Returns
`nothing` itself for a model it does not cover.
"""
function _newton_inverter(des::DualEquilibriumSolver, pq = nothing)
    pred = _newton_predictor_form(des.model, des.system, des.idx_aq)
    pred === nothing && return nothing
    aq, jref, lna = des.idx_aq, des.j_solvent, des.lna
    ns = length(des.system.species)
    tags = _STRIP_TAGS[]
    return function (c, ref, w, q, params)
        p = pq === nothing ? params : pq(q, params)
        out = _invert_aqueous_newton(lna, pred, ns, aq, jref, c, ref, w, p)
        return (out === nothing || isempty(tags)) ? out : _strip_tags(out, tags)
    end
end

# The log-amount below which a solute is placed rather than solved for: a
# hundred thousand times the activity floor, in moles.
const _NEWTON_TRACE_LOG = log(1.0e-25)

function _invert_aqueous_newton(lna, pred, ns, aq, jref, c, ref, w, p; maxit::Int = 50, tol = nothing)
    cv = Float64[_plain(x) for x in c]
    wv = Float64[clamp(_plain(x), -700.0, 20.0) for x in w]
    live = [t for t in eachindex(cv) if t != jref && isfinite(cv[t])]
    lndenom = log(_plain(ref) * pred.M_w)

    # The amounts of the phase at log-amounts `z` of the members `idx` (positions
    # in `live`), every other live member at `held`: the species outside the
    # phase do not enter its activities, and a dead member is held at zero.
    function amounts(z, idx, held)
        T = promote_type(eltype(z), eltype(held), typeof(ref))
        x = ones(T, ns)
        for (t, i) in enumerate(aq)
            x[i] = t == jref ? ref : zero(T)
        end
        for (k, t) in enumerate(live)
            x[aq[t]] = exp(held[k])
        end
        for (k, j) in enumerate(idx)
            x[aq[live[j]]] = exp(z[k])
        end
        return x
    end
    function residual(z, idx, held)
        h = lna(amounts(z, idx, held), p)
        return [h[aq[live[j]]] - c[live[j]] for j in idx]
    end

    # Two starts: the composition the solve holds, and the Debye–Hückel part's.
    z = wv[live]
    all_idx = collect(eachindex(live))
    worst(zz) = maximum(abs, Float64[_plain(r) for r in residual(zz[all_idx], all_idx, zz)]; init = 0.0)
    wp = _invert_aqueous(pred, cv, _plain(ref), wv, p, jref)
    if wp !== nothing
        zp = Float64[clamp(_plain(wp[t]), -700.0, 20.0) for t in live]
        worst(zp) < worst(z) && (z = zp)
    end

    # At the rounding of the potentials, which are tens to hundreds.
    tol = something(tol, 1.0e-12 * max(1.0, maximum(abs, view(cv, live); init = 0.0)))
    for _ in 1:3
        sig = [k for k in eachindex(live) if z[k] > _NEWTON_TRACE_LOG]
        zs = z[sig]
        F(zz) = Float64[_plain(r) for r in residual(zz, sig, z)]
        J(zz) = _plain.(ForwardDiff.jacobian(zz -> residual(zz, sig, z), zz))
        Fz = F(zs)
        φ = sum(abs2, Fz)
        converged = isempty(sig) || maximum(abs, Fz) <= tol
        it = 0
        while !converged && it < maxit
            it += 1
            δ = qr(J(zs), ColumnNorm()) \ (-Fz)
            all(isfinite, δ) || return nothing
            α = min(1.0, 30.0 / max(maximum(abs, δ), eps()))
            accepted = false
            for _ in 1:30
                zt = clamp.(zs .+ α .* δ, -700.0, 20.0)
                Ft = F(zt)
                φt = sum(abs2, Ft)
                if isfinite(φt) && φt <= (1 - 1.0e-4 * α) * φ
                    zs, Fz, φ, accepted = zt, Ft, φt, true
                    break
                end
                α /= 2
            end
            converged = maximum(abs, Fz) <= tol
            if !accepted
                # A residual at the rounding of the potentials stops the
                # halving before it stops the iteration.
                converged = maximum(abs, Fz) <= 100 * tol
                break
            end
        end
        converged || return nothing
        z[sig] = zs
        # The traces, from their activity coefficient at that composition: the
        # model's log-activity at a test amount above the floor, minus the log
        # of that molality.
        traces = [k for k in eachindex(live) if !(k in sig)]
        isempty(traces) && break
        probe = copy(z)
        for k in traces
            probe[k] = log(1.0e-20)
        end
        h = lna(amounts(Float64[], Int[], probe), p)
        moved = false
        for k in traces
            lnγ = _plain(h[aq[live[k]]]) - (log(1.0e-20) - lndenom)
            z[k] = clamp(cv[live[k]] - lnγ + lndenom, -700.0, 20.0)
            z[k] > _NEWTON_TRACE_LOG && (moved = true)
        end
        moved || break
    end

    # Found on the values. Whatever carries dual numbers (the potentials, the
    # temperature, a parameter of the model) shows in the type of the residual;
    # a Newton step in that type, with the Jacobian on the values, gives the
    # composition its first derivatives, by the implicit-function theorem. That
    # is what OptimaSolver differentiates through an inversion (the Jacobian of
    # its outer residual where an inversion has not converged); a nested
    # derivative of an equilibrium is lifted by its tangent at the answer, which
    # does not go through the inversion. A Jacobian taken in the dual type, which
    # a nested lift would need, nests its tag inside the caller's and ForwardDiff
    # cannot order the two.
    sig = [k for k in eachindex(live) if z[k] > _NEWTON_TRACE_LOG]
    R0 = residual(z[sig], sig, z)
    T = promote_type(eltype(R0), eltype(w), typeof(ref))
    zd = Vector{T}(z)
    if !(T <: AbstractFloat) && !isempty(sig)
        Jv = _plain.(ForwardDiff.jacobian(zz -> residual(zz, sig, z), z[sig]))
        zs = zd[sig]
        for _ in 1:2
            zs = zs - Jv \ residual(zs, sig, zd)
        end
        zd[sig] = zs
        traces = [k for k in eachindex(live) if !(k in sig)]
        if !isempty(traces)
            probe = copy(zd)
            for k in traces
                probe[k] = log(1.0e-20)
            end
            h = lna(amounts(T[], Int[], probe), p)
            for k in traces
                lnγ = h[aq[live[k]]] - (log(1.0e-20) - log(ref * pred.M_w))
                zd[k] = c[live[k]] - lnγ + log(ref * pred.M_w)
            end
        end
    end
    out = Vector{T}(undef, length(c))
    for t in eachindex(c)
        out[t] = w[t]
    end
    for (k, t) in enumerate(live)
        out[t] = zd[k]
    end
    return out
end
