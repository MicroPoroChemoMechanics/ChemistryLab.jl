# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

# ── Physical constants, in one place and taken from the library ──────────────
#
# `DynamicQuantities.Constants` carries the CODATA values WITH THEIR DIMENSIONS,
# and that is where they come from here. Writing `8.31446261815324` into a
# formula is wrong in three ways at once: it is a magic number a reader has to
# recognize, it silently asserts a unit system, and the next occurrence of it
# drifts from the first. There were nine of them in this package.
#
# Both forms are provided because both are needed and the difference matters.
# The DIMENSIONAL constant is the one to use when a quantity carries units; the
# stripped one is for the inner loops, where an activity model or a rate law
# works on plain `Float64` (and on `ForwardDiff.Dual`) and a `Quantity` would
# cost an allocation per evaluation. Deriving the second from the first, rather
# than typing it again, is what keeps them from drifting apart.

using DynamicQuantities
using DynamicQuantities: Constants as _DQConstants

"""
    R_GAS_Q

The molar gas constant with its dimensions, `8.314462618… J/(mol·K)`, taken from
`DynamicQuantities.Constants.R` (CODATA).

Use this wherever the surrounding quantities carry units. For the plain number
used inside an activity model or a rate law, see [`R_GAS`](@ref).
"""
const R_GAS_Q = _DQConstants.R

"""
    R_GAS

The molar gas constant as a plain `Float64` in `J/(mol·K)`, **derived from**
[`R_GAS_Q`](@ref) rather than written again.

This is the form the inner loops want: an activity model is evaluated once per
species per Newton iteration, on `Float64` and on `ForwardDiff.Dual`, and a
dimensional `Quantity` there costs an allocation each time for a factor whose
units are known at the call site anyway.
"""
const R_GAS = ustrip(us"J/(mol*K)", R_GAS_Q)

"""
    FARADAY_Q

The Faraday constant with its dimensions, `96485.332… C/mol`, from
`DynamicQuantities.Constants.F` (CODATA, exact by the 2019 SI redefinition).
"""
const FARADAY_Q = _DQConstants.F

"""
    FARADAY

The Faraday constant as a plain `Float64` in `C/mol`, derived from
[`FARADAY_Q`](@ref). Used by the Nernst relation between `pe` and `Eh`.
"""
const FARADAY = ustrip(us"C/mol", FARADAY_Q)

"""
    RT_over_F(T)

`RT/F` in volts at temperature `T` in kelvin — the natural scale of an
electrochemical potential, and the factor in Nernst's relation. At 25 °C,
`RT ln(10)/F` is the familiar `0.05916 V` per unit of `pe`.
"""
RT_over_F(T::Real) = R_GAS * T / FARADAY
