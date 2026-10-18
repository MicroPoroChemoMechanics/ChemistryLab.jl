# [The library of formulations](@id studio-library)

The *Start* tab of the studio offers formulations published in papers, each
loaded with one click: every tab is then filled, and the calculation can be run
at once, changed, and saved under a name of your own. Nothing is downloaded:
the formulations ship with the studio, and the compositions they use are those
ChemistryLab holds, with their sources, in its material templates.

Only formulations published in a paper are offered. Where an entry takes a
proportion from its paper (a replacement, a water-to-binder ratio), the tests
of the studio check it against the data ChemistryLab holds for that paper. The
degrees of reaction are illustrative unless the entry says otherwise: set them
for the age you are interested in, or let rate laws decide them.

| Formulation | What it computes | Sources |
| --- | --- | --- |
| CEM I 52.5 N paste | one equilibrium at 28 days, w/b 0.5 | [Gruyaert2010](@cite), [LothenbachWinnefeld2006](@cite) |
| Portland cement paste, RILEM round robin | one equilibrium; its volumes are unknown, three minor oxides having no sourced density | [Durdzinski2017](@cite), [LothenbachWinnefeld2006](@cite) |
| Portland cement with 40 % slag, RILEM round robin | one equilibrium of the paste PC-S1, w/b 0.4 | [Durdzinski2017](@cite), [LothenbachWinnefeld2006](@cite) |
| Portland cement with 30 % siliceous fly ash, RILEM round robin | one equilibrium of the paste PC-SFA, w/b 0.4 | [Durdzinski2017](@cite), [LothenbachWinnefeld2006](@cite) |
| Slag cements: the replacement of a CEM I 52.5 N by slag | a path along the slag-to-binder ratios of the paper, 0 to 0.85 | [Gruyaert2010](@cite), [LothenbachWinnefeld2006](@cite) |
| Hydration of a CEM I 52.5 N in time, by rate laws | a kinetic run over 28 days, the clinker phases by the law of [ParrottKilloh1984](@citet) with the parameters of [Lavergne2018](@citet) | [Gruyaert2010](@cite), [Lavergne2018](@cite), [LothenbachWinnefeld2006](@cite) |
| Carbonation of a CEM I 52.5 N paste | a path of CO₂ added, up to 1 mol per 100 g of binder | [Gruyaert2010](@cite), [LothenbachWinnefeld2006](@cite) |
| Leaching of a CEM I 52.5 N paste by pure water | twelve renewals of the pore solution by 5 kg of water | [Gruyaert2010](@cite), [LothenbachWinnefeld2006](@cite) |

The phases of every entry are those [LothenbachWinnefeld2006](@citet) list for
a Portland paste, in the thermodynamic database Cemdata18
[Lothenbach2019](@cite).
