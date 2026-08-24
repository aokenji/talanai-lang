# Validation redock of the 10 screened compounds

**Revision 2, 2026-08-03.** Revision 1 (2026-07-31) was written before the
adversarial expert review and **must not be executed**: its sole go/no-go gate
was a control that cannot fail. See section 4 for what replaced it.

This document is written to be executed by a session with none of the
originating conversation in front of it.

Companion documents: `REVIEW-FINDINGS.md` (all 17 surviving review findings),
`BUILD-SPEC.md` section 12 (the finding and its decomposition).

---

## 0. Read this first: what "no going back" actually means

The redock is **not destructive**. It writes new files into
`D:\THESIS_VSC\talanai-lang\validation-run\` and nothing else.

- `compounds.js`, `protocol.js` and everything under `D:\BALAKATDBV2` are
  **read only**. Nothing here writes there.
- The defended thesis numbers stay where they are, including the Zenodo deposit
  `10.5281/zenodo.20384660`, which is immutable regardless.
- The TalanaiDock assets are read only. No receptor or ligand is regenerated in
  place.

What is irreversible is **knowing**. Once these numbers exist, continuing to
present the unvalidated ranking without mentioning them stops being an
oversight and becomes a choice.

The honest reading is that this strengthens the work. A study that validates
its own protocol, finds a problem and reports it is worth more than one that
never looked.

---

## 1. Why this run is happening

**The screening receptor carried no hydrogen-bond donors or acceptors at all.**

The screen passed the raw `receptor_clean.pdb` directly to Vina. Vina reads
AutoDock atom types from columns 77-78, which in a PDBQT hold `OA`, `NA`, `SA`
and `HD`, and in a PDB hold plain element symbols. Plain `O` and `N` neither
donate nor accept. Vina accepts a `.pdb` receptor silently, with no warning.

Measured with `--score_only`, glucose at its crystal coordinates:

| Receptor | Ligand-receptor energy |
|---|---|
| `prepared/receptor.pdbqt` (Meeko) | **−5.732** kcal/mol |
| `receptor_clean.pdb` (the screen) | **−0.067** kcal/mol |
| both, under `--weight_hydrogen 0` | −0.067 and −0.067 |

The third row is the proof: switch off the hydrogen-bond term and the two
receptors become identical, so the raw receptor's H-bond term is exactly zero.
For this ligand that term is 98.8 percent of the interaction energy. The screen
ranked candidates on shape and hydrophobic burial alone.

Type census: `prepared/receptor.pdbqt` has 923 OA, 14 NA, 13 SA, 1061 HD.
`receptor_clean.pdb` has 923 O, 804 N, 13 S and zero hydrogens.

**The mechanism is atom typing, not charges.** Vina has no electrostatic term,
so Gasteiger charges are irrelevant to it. Missing polar hydrogens account for
roughly half the effect; the rest is `O` to `OA` and `N` to `NA` acceptor
typing, which adding hydrogens alone would not fix. Any wording that blames
charges is wrong and should be corrected wherever it appears.

The redock failure is the symptom. The `--score_only` decomposition is the
diagnosis, it takes one minute, and it belongs in the September deck.

### The 2x2 that isolated it

| Case | Receptor | Box | Exh | Score | RMSD |
|---|---|---|---|---|---|
| B calibration | Meeko | 18 Å | 16 | −5.885 | 0.513 |
| D search only | Meeko | 30 Å | 8 | −5.905 | 0.518 |
| C prep only | raw | 18 Å | 16 | −4.208 | 5.592 |
| A screening | raw | 30 Å | 8 | −4.024 | 6.164 |

Case D shows the 30 Å box at exhaustiveness 8 recovers the crystal pose
perfectly well. The box was never the problem.

Case B agrees with the published `validation.json` (0.519 Å, −5.897) to 0.006 Å
and 0.012 kcal/mol. Note this is **agreement between two runs**, not a
calibration of the RMSD implementation, and `validate_redock.py` used
`--energy_range 4` against this document's 3.

### What the consequence is, stated carefully

Re-docking the original five on the prepared receptor at exhaustiveness 8:

| Compound | Thesis (raw) | Prepared | Δ |
|---|---|---|---|
| Rutin | −8.857 | −10.620 | −1.763 |
| Betulinic acid | −8.290 | −9.644 | −1.354 |
| Quercetin | −7.503 | −8.889 | −1.386 |
| Kaempferol | −7.479 | −8.240 | −0.761 |
| Oleanolic acid | −6.922 | −7.792 | −0.870 |
| Acarbose | −6.660 | −8.576 | −1.916 |

Three seeds on the compounds whose order changed: acarbose mean −8.730 (spread
0.362), kaempferol −8.342 (0.265), oleanolic acid −7.809 (0.032).

**The defensible claim is "3 or 4 of 5, with kaempferol indistinguishable from
acarbose at this sampling".** Not "3 of 5". Revision 1 said "acarbose's worst
seed beats kaempferol's best"; that is not derivable from a mean and a range,
and the nine per-seed values were never written to disk. With n=3 against n=3 a
permutation test tops out at p=0.05, so complete separation is the ceiling of
this design, not a significant result. Oleanolic acid's demotion (0.921 gap
against a 0.032 spread) survives any reading. Kaempferol's does not. Quercetin
at −8.889 against acarbose's −8.730 is a 0.16 gap from a single seed and is
also unsettled.

**Do not write "the validated receptor".** The redock validates pose recovery
for one small sugar. No enrichment test, no decoys and no known actives were
run, so neither receptor is validated for *ranking*. Write "the receptor that
passes the redocking control".

### The Δ column is not a single-variable measurement

`reference/reference.json` records acarbose at **−5.932** on the raw receptor
under the exact thesis configuration, against a published −6.660. That is 0.728
kcal/mol, thirty-five times the largest deviation among the five compounds
(0.021), and acarbose is the one compound absent from the reproducibility
table. The thesis's own acarbose input file does not exist on disk; the
available `acarbose.pdbqt` was regenerated from SMILES. So the −1.916 Δ changed
the receptor **and** the ligand file. Like for like it is −2.644.

---

## 2. Toolchain, verified 2026-07-31

| Component | Version | Location |
|---|---|---|
| AutoDock Vina | **1.2.7** | `...\docking_assets\vina.exe` |
| Python (system) | 3.12.10 | `python` on PATH |
| Python (TalanaiDock bundled) | 3.12.10 | `...\dist\TalanaiDock\python\python.exe` |
| RDKit | 2026.03.3 | both interpreters |
| Meeko | 0.7.1 | both interpreters |
| PLIP | present | bundled interpreter only |
| Talanai | 0.1.0 | `D:\THESIS_VSC\talanai-lang` |

Asset root:
`D:\THESIS_VSC\BALAKATDBV2\talanaidock\dist\TalanaiDock\app\docking_assets`

Talanai's validator stays standard library only and `tal check` must keep
working with nothing installed. RDKit and Meeko are used here as **upstream
preparation and cross-check** tools, run separately, exactly as the original
pipeline used them.

### Input checksums (SHA-256, first 16 hex)

Verify before running. A mismatch means the inputs changed and the comparison
to the thesis numbers is no longer clean.

| File | SHA-256 (16) | Bytes |
|---|---|---|
| `vina.exe` | `E0C4B2715E0C1A74` | 1233920 |
| `docking_data\receptor_clean.pdb` | `BE48D1A17CBF7790` | 396634 |
| `docking_data\prepared\receptor.pdbqt` | `87068D49698B6519` | 477576 |
| `docking_data\prepared\receptor_H.pdb` | `AD102AF2ED78C79A` | 761643 |
| `docking_data\grid_center.json` | `F6DA6416C147F14E` | 149 |
| `validation\glc_ligand.pdbqt` | `22C7165D9DDB1D6A` | 1788 |
| `validation\glc_3a4a_601.pdb` | `AD7D342694349CED` | 989 |
| `validation\glc_crystal_ref.sdf` | `233FF304162EF521` | 2123 |
| `reference\acarbose.pdbqt` | `0DC1261C7265A359` | 6267 |

Ligands: quercetin `EEB05B8DD0063F6B`, rutin `2353B3EEF3C6D689`, kaempferol
`4B1C5E91F596F06A`, betulinic `D8453D141A887A09`, oleanolic `9701B4C62764D5CD`,
luteolin `A0B055B5C6AE79F2`, vitexin `EAD50138728D9E77`, isovitexin
`72AB44640EFD7F2E`, ursolic `40076B5892E5BC19`.

---

## 3. The compound set

`compounds.js` as of 2026-07-31 01:20: **10 screened, 32 candidates, 1
reference.** Torsion counts must be added to every results table, because
sampling adequacy scales with them and the review found this is where the
comparison breaks.

| # | Compound | Thesis (raw) | Formula | Heavy | Prepared file |
|---|---|---|---|---|---|
| 1 | Rutin | −8.857 | C27H30O16 | 43 | `rutin.pdbqt` |
| 2 | Betulinic acid | −8.290 | C30H48O3 | 33 | `betulinic.pdbqt` |
| 3 | Ursolic acid | −8.105 | C30H48O3 | 33 | `ursolic.pdbqt` |
| 4 | Isovitexin | −8.024 | C21H20O10 | 31 | `isovitexin.pdbqt` |
| 5 | Spinosin | −7.859 | C28H32O15 | 43 | **see note** |
| 6 | Luteolin | −7.733 | C15H10O6 | 21 | `luteolin.pdbqt` |
| 7 | Quercetin | −7.503 | C15H10O7 | 22 | `quercetin.pdbqt` |
| 8 | Kaempferol | −7.479 | C15H10O6 | 21 | `kaempferol.pdbqt` |
| 9 | Vitexin | −7.469 | C21H20O10 | 31 | `vitexin.pdbqt` |
| 10 | Oleanolic acid | −6.922 | C30H48O3 | 33 | `oleanolic.pdbqt` |
| ref | Acarbose | −6.660 | C25H43NO18 | 44 | `acarbose.pdbqt`, **not the thesis file** |

**Spinosin note, and it is more serious than revision 1 assumed.** Investigated
2026-08-03. `compounds.js` records `pubchemCID: 155692` and −7.859, but **no
spinosin structure file exists anywhere**: not in `D:\BALAKATDBV2`, not in the
frozen TalanaiDock checkout (commit `053de17`, 2026-06-18), and not in the
`ziziphus-docking` Hugging Face Space, checked both locally and against the
live Space file tree via its API.

`interactions.js:11-13` states plainly that "Ursolic Acid and Spinosin come
from the TalanaiDock Vina + PLIP pipeline (July 2026)" and that "the two runs
report different distance conventions". Every local snapshot predates July, so
the run that produced −7.859 used a newer app or Space version than any copy on
this machine, or was an ephemeral bring-your-own run whose input was never
persisted. **The origin is unverified and the input file is unrecoverable.**

Two consequences:

1. By this project's own rule R105, −7.859 is not currently citable. It has no
   findable run behind it.
2. **Ursolic acid and Spinosin may not be comparable to the other eight.** They
   went through a different pipeline run, and the code itself says the two runs
   use different distance conventions. That is an R601-class question about the
   published set, independent of the receptor problem, and it should be settled
   before September.

A freshly prepared file now exists at `validation-inputs\spinosin.pdbqt`
(SHA-256 `af4133c660d57773...`, 43 heavy atoms, C28H32O15 confirmed, TORSDOF
16, SMILES from PubChem CID 155692, built through the exact
`_prep_assets_local.py` route: ETKDGv3 seed `0xF00D`, MMFF94 maxIters 400,
Meeko 0.7.1). The methods must state it was prepared on 2026-08-03 and not by
the original pipeline.

**Acarbose note.** `reference/acarbose.pdbqt` was regenerated from SMILES. The
file the thesis docked does not exist. State this in the methods.

---

## 4. The control battery, which replaces revision 1's single gate

Revision 1 gated everything on redocking `glc_ligand.pdbqt`. The review
established that this control **cannot fail**: of its six declared torsions,
five rotate only a hydroxyl hydrogen, so heavy-atom pose recovery is a
rigid-body problem over 12 atoms; the file carries crystal coordinates verbatim
and so never exercises the ligand pipeline; and the box is centred on its own
centroid. A gate that cannot fail is not a gate.

Three controls now run, in this order. **Control 3 is the go/no-go.**

### Control 1: crystal glucose, retained for continuity

Redock `glc_ligand.pdbqt` into the prepared receptor, screening box,
exhaustiveness as chosen in section 6, seed 42. Expect roughly 0.52 Å. Its only
purpose is continuity with `validation.json`. **It does not gate anything.**

### Control 2: re-embedded glucose, tests the ligand pipeline

Same, but with glucose built from SMILES through ETKDGv3 + MMFF94 + Meeko, the
route every screened compound went through, with no crystal coordinates used.
File: `validation-inputs\glc_reembedded.pdbqt` (SHA-256 `d28c32a3e6f59ffc...`,
12 heavy atoms, TORSDOF 6, SMILES verified against PubChem CID 79025 by both
canonical SMILES and InChI).

**Be precise about what this does and does not test.** Verified 2026-08-03: the
re-embedded file has *identical* topology to the crystal-derived one, same
rotatable bonds, same TORSDOF 6, because Meeko assigns rotatable bonds from
connectivity rather than coordinates. And Vina randomises ligand position,
orientation and torsions at the start of each Monte Carlo run, so the input
coordinates matter little to the search itself.

So Control 2 is **not** a harder search. What it tests is whether the ligand
pipeline produces a chemically correct molecule: right protonation state, right
AutoDock atom types, right rotatable-bond assignment, no distorted ring
geometry. Given that mis-typing the *receptor* is the entire finding of this
study, checking that the same toolchain types a *ligand* correctly is worth one
run. Do not describe it as testing search difficulty.

The torsional-difficulty gap is Control 3's job, not this one's.

### Control 3: isomaltose from 3AXH, cross-docked. The real gate.

**Corrected 2026-08-03. Revision 2 originally specified maltose from PDB 3AJ7.
That was wrong and could not have run.** 3AJ7 contains no sugar at all, only
one calcium ion and 608 waters. The error was mine: I took "the 3A4A companion
structure with maltose" from a review suggestion and never checked the ligand
records.

All five PDB entries of this enzyme (UniProt P53051) were then checked:

| PDB | Res. | Ligand | Intact disaccharide? |
|---|---|---|---|
| 3AJ7 | 1.30 Å | none | no |
| 3A4A | 1.60 Å | GLC 601, 12 atoms | no, one ring only |
| 3A47 | 1.59 Å | none | no |
| 3AXI | 1.40 Å | GLC 601, 12 atoms | no, **despite its title saying "in complex with maltose"** |
| 3AXH | 1.80 Å | GLC-GLC, 23 atoms, O6-C1 linked | **yes, isomaltose** |

**No structure of this enzyme contains intact maltose.** Yamamoto 2010's own
abstract explains why: in the maltose-soaked crystal only the nonreducing-end
glucose gave density, with incomplete density at the reducing end. That
structure is 3A4A itself, which is why 3A4A's ligand is a lone glucose.

Note 3AXI as a trap for anyone working from titles rather than coordinates: it
is *titled* as a maltose complex and contains only glucose.

**The substitute is scientifically better than the original plan.** 3A4A is an
*isomaltase*, an oligo-1,6-glucosidase. Isomaltose is its cognate α-1,6
substrate; maltose is α-1,4 and is the poorer substrate. Validating on the
enzyme's actual substrate is more defensible than validating on maltose, so
the correction improves the control rather than merely rescuing it.

**The control:** isomaltose, prepared from SMILES through the standard ligand
route, cross-docked from 3AXH into 3A4A. 23 heavy atoms against glucose's 12,
and the flexibility gap this control exists to close is real rather than
nominal. **Require rank-1 RMSD under 2.0 Å against the 3AXH isomaltose pose
superposed into 3A4A's frame. If it fails, stop and report. Do not dock the
compounds.**

**Superposition is mandatory and is not optional bookkeeping.** 3AJ7 and 3A4A
happen to share a crystallographic frame (backbone RMSD 0.110 Å over 586 CA),
but **3AXH does not**: its isomaltose sits 61 Å from the 3A4A box centre
unaligned. After Kabsch superposition onto 3A4A it lands **2.676 Å from the
box centre**, comfortably inside the 30 Å box and 6.5 to 8.8 Å from the
catalytic triad. So the existing box is correct for this control, but only
against a superposed reference. Comparing against raw 3AXH coordinates would
produce a 61 Å RMSD and a false failure.

Caveats to state in the methods: 3AXH is 1.80 Å against 3A4A's 1.60 Å, and it
comes from a later paper (Yamamoto et al. 2011) than the structure being
docked into. Cross-docking between the two is the point, since it tests induced
fit, but the resolution difference belongs in the record.

A `maltose.pdbqt` was also prepared (23 heavy atoms, TORSDOF 12 against
glucose's 6) before the structure problem was found. It is retained in
`validation-inputs\` but **cannot serve as a validation control**, because no
crystal reference exists to measure an RMSD against. Do not promote it to a
gate.

### RMSD, reported honestly for every control

Report both:
1. Talanai's `tal control` value, which is a **lower bound** when atom
   correspondence cannot be established, labelled as such.
2. The RDKit symmetry-corrected value, mirroring `validate_redock.py`, against
   `glc_crystal_ref.sdf` or the equivalent reference.

A lower bound above the threshold is conclusive for failure. A lower bound
below the threshold is **not** a pass. Only the symmetry-corrected number can
confirm a pass.

---

## 5. The run

### 5.1 Parameters

Fixed, and not negotiable if the comparison is to mean anything:

- Receptor: `prepared\receptor.pdbqt`
- Box centre `21.52, -7.70, 23.55`, size `30 x 30 x 30 Å`, unchanged from the
  thesis and shown adequate by case D
- Engine Vina 1.2.7, `num_modes 9`, `energy_range 3`, `cpu 4`
- Seeds **42, 43, 44**, reported as mean **and** every individual value
- Ligand preparation unchanged, and the existing prepared files reused
  byte-for-byte so that only the receptor differs

### 5.2 Acarbose gets converged separately

Acarbose has 26 active torsions against 3 for oleanolic acid, and it anchors
every claim in the study while being its least converged number. Its three-seed
spread was 0.362.

Escalate acarbose alone: exhaustiveness 32, then 64, then 128, three seeds
each, until its spread approaches the rigid ligands' 0.03. Report the
convergence curve. Note also that Vina divides the intermolecular term by
(1 + 0.05846 x N_rot), which is 2.52 for acarbose against 1.29 for kaempferol,
so under-sampling and the torsion penalty both push toward "everything beats
acarbose".

### 5.3 Also re-dock acarbose on the raw receptor

Seeds 42, 43, 44, thesis configuration, using `reference/acarbose.pdbqt`.
Publish the mean beside both −6.660 and −5.932. Until that exists, the Δ column
must be dropped or labelled "receptor plus unexplained".

### 5.4 Keep the poses

Revision 1 discarded them, which made "3 of 5 beat acarbose" a score-only claim
with no evidence any pose touched the active site. For every rank-1 pose,
record and publish the minimum heavy-atom distance to Asp215, Glu277 and
Asp352. A compound that scores well without contacting the catalytic residues
does not support an active-site inhibition claim.

### 5.5 Outputs

Everything under `D:\THESIS_VSC\talanai-lang\validation-run\`:

- `run-record.json` per ligand per seed: full command, input checksums, score,
  pose file path, Talanai version, RMSD method
- `results.csv`: thesis value, each seed value, mean, spread (defined as range,
  stated explicitly), Δ, heavy atoms, **torsion count**, ligand efficiency,
  minimum distance to each catalytic residue, and whether it beats acarbose
- `methods.txt`, generated from the parameters that actually ran
- `validated.tal`, which must pass `tal check`
- every pose as `.pdbqt`

### 5.6 Not part of this run

PLIP re-analysis of contacts. The current interaction data derives from poses
on the raw receptor and inherits the same question, but re-running it is a
larger change and a separate decision.

---

## 6. Decisions, settled 2026-07-31

**A. Exhaustiveness 32** for the main run. R401's rule of thumb for a 27,000 Å³
box, and the seed spreads at exhaustiveness 8 (0.27 to 0.36 on flexible
ligands, 0.03 on the rigid triterpene) are direct evidence that 8 under-samples.
Roughly 3.2 hours for 11 ligands across 3 seeds on this i3, plus the acarbose
escalation in 5.2.

Report the already-measured exhaustiveness 8 numbers alongside, since 32
changes the receptor **and** the sampling at once and the two effects should
stay separable.

**B. Prepare Spinosin** through the original route so all 10 are covered, with
the methods stating it was prepared later and not by the original pipeline.

---

## 7. Safety rules

1. **Never write to `D:\BALAKATDBV2`.** Read only, always.
2. **Never edit the thesis numbers.** New numbers live beside them, labelled.
3. **Verify the section 2 checksums first.** If one differs, stop and say so.
4. **Report failures as findings.** A control that does not pass is a result.
5. **No single-seed results, and publish every per-seed value**, not just the
   mean and a spread. Revision 1 broke this rule in the sentence the whole run
   exists to support.
6. **State the RMSD method every time**, and whether it is exact or a bound.
7. **Write a run record for every run, at the time.** Rule R105 says a number
   without a findable run is not citable, and the project has already violated
   its own rule once.

---

## 8. Resuming in a fresh session

Say: **"Run the validation redock from REDOCK-PROTOCOL.md."**

Order: verify checksums, run controls 1 and 2 for information, run control 3 as
the gate, stop if it fails, then dock 11 ligands across 3 seeds, escalate
acarbose separately, keep every pose, and write the outputs in 5.5.

Report the table with thesis values beside the new values, both rankings, the
torsion counts and the catalytic-residue distances. Do not touch
`D:\BALAKATDBV2`.
