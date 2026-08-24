# Control battery results, 2026-08-03

Run with `python run_controls.py`. Receptor `prepared/receptor.pdbqt`, box
30 x 30 x 30 Å at (21.52, −7.70, 23.55), exhaustiveness 32, seed 42,
num_modes 9, energy_range 3. Records in `validation-run/controls/`.

## The numbers

| Control | Score | Position RMSD | Shape RMSD | Verdict |
|---|---|---|---|---|
| 1. Crystal glucose | −5.905 | **0.518 Å** | 0.128 Å | PASS |
| 2. Re-embedded glucose | −5.775 | **5.358 Å** | 1.265 Å | FAIL |
| 3. Isomaltose, the gate | −8.084 | **2.895 Å** | 2.150 Å | **FAIL** |

Position RMSD is symmetry-corrected with **no superposition** (RDKit
`CalcRMS`): did the pose land in the right place. Shape RMSD superimposes
first (`GetBestRMS`): did it find the right shape, wherever it ended up. Only
the first is a redocking control. An earlier version of `run_controls.py` used
`GetBestRMS` for the verdict, which would have passed a pose 5 Å away with
perfect internal geometry. Corrected before any verdict was reported.

**Cross-validation of the two implementations.** On control 1, Talanai's
standard-library lower bound and RDKit's symmetry-corrected value agree exactly
at 0.518 Å, and both reproduce the published `validation.json` figure of
0.519 Å. Three independent routes to the same number.

## Every pose is in the right pocket

Minimum heavy-atom distance from each rank-1 pose to each catalytic residue:

| Pose | Asp215 | Glu277 | Asp352 |
|---|---|---|---|
| 1. Crystal glucose | 2.77 Å | 2.87 Å | 2.75 Å |
| 2. Re-embedded glucose | 3.23 Å | 2.75 Å | 3.68 Å |
| 3. Isomaltose | 2.83 Å | 2.69 Å | 2.83 Å |

All three are hydrogen-bonded to the entire catalytic triad. **The box is
right, the site is right, and nothing is docking into a decoy pocket.**

## What actually failed

Not the pocket. The **orientation within** it.

Control 2 is the finding. The same molecule, glucose, into the same receptor,
same box, same settings, same seed:

- fed the **crystal coordinates**, it recovers the crystal pose at 0.518 Å
- built through **the study's own ligand pipeline** from SMILES via ETKDGv3 +
  MMFF94 + Meeko, the exact route every screened compound took, it lands
  **5.358 Å away**, still contacting all three catalytic residues, scoring
  **0.13 kcal/mol** differently

A 5.4 Å RMSD on a 12 heavy-atom sugar that still hydrogen-bonds the whole
triad means the ring is placed in a different orientation in the same pocket,
not in a different site.

So: **the scoring function cannot distinguish the crystallographic glucose
orientation from an alternative 5.4 Å away.** They differ by less than the
engine's own run-to-run noise. The published 0.519 Å validation lands on the
correct one because the correct one was supplied as the starting conformation,
and Vina perturbs from its input rather than sampling wholly independently of
it.

Control 3 tells the same story one step harder. Isomaltose, cross-docked from
3AXH, is 2.895 Å from its crystal pose while hydrogen-bonded to all three
catalytic residues at 2.69 to 2.83 Å. Right pocket, right contacts, placement
outside the 2.0 Å bar.

## What this does and does not undermine

**Does not undermine:** that the compounds occupy the catalytic site. Every
control lands there and engages the triad. Binding-affinity rankings are
driven by pocket occupancy and burial, and nothing here says the site is wrong.

**Does undermine:** any claim about *which specific residues a given compound
contacts and how*. The PLIP interaction data describes orientations that this
protocol demonstrably cannot pin down to better than several angstrom. Two
orientations differing by 5.4 Å will produce different hydrogen-bond
inventories, and the scoring function ranks them 0.13 kcal/mol apart.

**Also undermines:** the published redocking validation as evidence about the
screening pipeline. It validates docking a crystal-coordinate file. Every
screened compound came from the SMILES pipeline instead, and that route fails
to recover a known pose for the simplest ligand available.

## Status: the screen did not run

The gate failed, so per section 4 of the protocol nothing downstream was
started. The threshold was not relaxed to get through it; a gate that moves
when it fails is not a gate.

## Open options, none of them taken unilaterally

1. **Report the protocol as unvalidated at 2.0 Å for pose claims**, keep the
   affinity work, and drop or heavily caveat the interaction analysis.
2. **Argue a stated allowance in the methods** for cross-docking a
   disaccharide from a 1.80 Å structure into a 1.60 Å one. Defensible in
   writing, but it must be argued out loud, not applied silently.
3. **Test whether the orientation degeneracy is a sampling problem** by
   re-running control 2 at much higher exhaustiveness and across seeds. If a
   harder search recovers the crystal orientation, the fix is search effort
   rather than a limit of the scoring function. This is the cheapest
   informative next experiment and it is a few minutes per run.
4. **Restrict pose claims to the ensemble level**: report which residues are
   contacted across the top N poses rather than by the single rank-1 pose.

Option 3 should be run before any of the others is chosen, because it decides
whether this is a fixable sampling limit or a scoring-function limit.
