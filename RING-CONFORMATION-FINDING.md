# The ligand preparation freezes a single ring conformation

**2026-08-03.** Root cause of the control-2 failure, traced from a failed
redock to a single line of code.

## The causal chain, each step verified

1. `_prep_assets_local.py:41` calls `AllChem.EmbedMolecule` and generates
   **one** conformer, seeded, from ETKDGv3.
2. Line 43 runs `MMFFOptimizeMolecule` with 400 iterations. That is a **local**
   minimisation and cannot cross a ring-flip barrier, which for a pyranose is
   on the order of 10 kcal/mol.
3. Meeko writes ring bonds as rigid. This is correct; PDBQT has no other option.
4. **AutoDock Vina cannot change ring geometry at all.**

So whichever ring conformation ETKDG happened to produce on that single attempt
is frozen for the entire study.

## Demonstrated for glucose

| Molecule | Ring torsions | Signs | Mean abs | Verdict |
|---|---|---|---|---|
| Crystal `glc_ligand.pdbqt` | 53.9, −62.9, 64.8, −59.2, 50.7, −46.9 | `+-+-+-` | 56.4° | textbook chair |
| Pipeline `glc_reembedded.pdbqt` | 29.1, 25.6, −57.2, 25.7, 33.2, −62.5 | `++-++-` | 38.9° | twist-boat |
| Docked pose, exhaustiveness 512 | 29.1, 25.6, −57.2, 25.6, 33.2, −62.5 | `++-++-` | 38.9° | **unchanged** |

The docked pose carries the input torsions to one decimal place. Vina did not
alter the ring, because it cannot. The crystal pose was therefore unreachable,
which is exactly what the escalation showed: rank-1 RMSD pinned at 5.35 to
5.36 Å across exhaustiveness 128 and 512 and three seeds, scores varying by
0.002 kcal/mol. A converged search on an unreachable target.

## How far it reaches

Audit of every non-aromatic six-ring in the prepared ligands
(`audit_ring_conformations.py`). No docking, pure geometry on the files that
were actually docked.

**Genuine defects.** These are saturated pyranose or cyclohexane rings with no
sp2 carbon, so they have no excuse not to be chairs:

| Compound | Rings flagged | Note |
|---|---|---|
| **Isovitexin** | 1 of 1 | sugar ring, signs `++-++-`, 38.4°, the same pattern as the broken glucose |
| **Acarbose (reference)** | 2 of 3 | at least one is a sugar ring. **Every "beats acarbose" claim is measured against this molecule** |
| Spinosin | 1 of 2 | sugar ring. This is the file *I* prepared on 2026-08-03, same pipeline, same defect |
| Isomaltose | 1 of 2 | sugar ring. Also mine, and it is the intended gate ligand |

**Clean:** Rutin (2 of 2 chairs), Betulinic acid (4 of 4), Vitexin (1 of 1).
**Not applicable:** Quercetin, Kaempferol, Luteolin have no saturated six-ring;
aromatic rings are planar by construction and cannot suffer this.

**Probable false positives, stated as such.** The heuristic flags any ring that
is not an alternating-sign chair near 55°, and a ring containing an sp2 carbon
*cannot* be a chair:

- **Oleanolic acid**, 1 of 5 flagged at 30.0° with alternating signs. Oleanane
  triterpenes carry a Δ12-13 double bond, so one half-chair is expected and
  correct. Almost certainly fine.
- **Ursolic acid**, 5 of 5 flagged, 34 to 43°, with non-alternating sign
  patterns. The Δ12 double bond explains one half-chair, not five distorted
  rings. This looks like a genuine defect, and the contrast with oleanolic
  acid, a near-identical scaffold that came out with four clean chairs, is
  what a single-conformer embedding lottery would look like. **Needs a proper
  check before being called either way.**

## What this does and does not mean

**Does not mean** the compounds do not occupy the catalytic site. Every control
pose hydrogen-bonds the full triad. Nor does it invalidate the whole affinity
ranking on its own.

**Does mean** that for the affected ligands the docked score was computed on a
strained, non-native conformer, and that pose-level claims about them describe
geometry the search was never able to explore.

**The sharpest consequence:** acarbose. It is the yardstick for every
comparison in the study, it is the most conformationally complex ligand in the
set, and two of its three saturated rings are non-chair.

## The fix, one line

```python
# current, _prep_assets_local.py:41
AllChem.EmbedMolecule(mol, params)

# what it should be
cids = AllChem.EmbedMultipleConfs(mol, numConfs=50, params=params)
AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=400)
# then keep the lowest-energy conformer
```

Generating one conformer and minimising locally is a lottery on ring pucker.
Generating many and keeping the lowest-energy one is the standard practice this
pipeline was missing.

## Status

Nothing has been re-prepared or re-docked. This is a section-3 item under the
delegation charter: it changes what the study can conclude, so the decision is
the author's.

Ordered by value:

1. **Verify the acarbose rings properly**, distinguishing its cyclohexene
   (legitimately a half-chair) from its pyranose rings (which must be chairs).
   Cheap, and it decides how serious this is.
2. **Settle ursolic acid** the same way.
3. **Fix the pipeline** to embed many conformers and keep the lowest.
4. **Re-prepare and re-dock** whatever turns out to be genuinely affected, and
   report the deltas beside the published values.
