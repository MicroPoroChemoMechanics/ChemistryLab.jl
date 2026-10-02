# =============================================================================
#  validation_common.jl — what the validation pastes share
#
#  Included by the scripts of the published pastes (lothenbach_winnefeld_2006.jl,
#  de_weerdt_2011.jl), so that a pore solution and a budget are read the same way
#  on every one of them.
# =============================================================================

using ChemistryLab
using DynamicQuantities
using OrderedCollections

# Included by every script that needs it, so that each runs on its own; defined
# once per session, since `Pkg.test` warns about every method redefined.
if !@isdefined(pore_solution_mmol)
    """
        pore_solution_mmol(state) -> Dict

    The dissolved elements of an equilibrium, and the hydroxide, in mmol per kg of
    water, the elements summed over every aqueous species that holds them.
    """
    function pore_solution_mmol(state)
        cs = state.system
        n = ustrip.(us"mol", state.n)
        iw = only(cs.idx_solvent)
        kg = n[iw] * ustrip(us"kg/mol", cs.species[iw][:M])
        out = Dict{String, Float64}()
        for i in cs.idx_solutes, (el, k) in atoms(cs.species[i])
            (el === :H || el === :O) && continue
            out[String(el)] = get(out, String(el), 0.0) + 1000 * k * n[i] / kg
        end
        out["OH-"] = 1000 * n[findfirst(s -> symbol(s) == "OH-", cs.species)] / kg
        return out
    end

    """
        budget_elements(cs, b) -> OrderedDict

    A budget of `cs`, written in its primaries, as amounts of the elements (and of
    the charge, `Zz`), in mol: what a second code is given to replay it.
    """
    function budget_elements(cs, b)
        content(p, el) = el === :Zz ? Float64(charge(p)) : Float64(get(atoms(p), el, 0))
        return OrderedDict(
            String(el) => sum(content(p, el) * x for (p, x) in zip(cs.SM.primaries, b)) for el in cs.CSM.primaries
        )
    end
end
