# Validation-input provenance

Generated 2026-08-03. This directory holds ligand files prepared for a docking
validation run against the reproducible pipeline documented in
`D:\THESIS_VSC\BALAKATDBV2\talanaidock\dist\TalanaiDock\app\_prep_assets_local.py`
(read-only source; not modified).

Tool versions (verified on both interpreters: system Python and the bundled
`D:\THESIS_VSC\BALAKATDBV2\talanaidock\dist\TalanaiDock\python\python.exe`):

- Python 3.12.10
- RDKit 2026.03.3
- Meeko 0.7.1

Preparation route (copied verbatim from `_prep_assets_local.py`'s `prepare_one()`,
same random seed `0xF00D`):

```python
mol = Chem.MolFromSmiles(smiles)
n_heavy = mol.GetNumAtoms()              # heavy-atom count, pre-AddHs
mol = Chem.AddHs(mol)
params = AllChem.ETKDGv3()
params.randomSeed = 0xF00D
params.maxIterations = 2000
AllChem.EmbedMolecule(mol, params)
AllChem.MMFFOptimizeMolecule(mol, maxIters=400)   # MMFF94, bounded
setups = MoleculePreparation().prepare(mol)
pdbqt, ok, err = PDBQTWriterLegacy.write_string(setups[0])
```

Run with:

```
python prep_ligand.py
```

(script kept at
`C:\Users\himaru\AppData\Local\Temp\claude\D--THESIS-VSC-BALAKATDBV2\dc2ea37e-f6ff-49f2-93e2-eedaf9df7839\scratchpad\prep_ligand.py`
during this session; it embeds the same `prepare_one()` body unchanged, targeting
the two SMILES below, and writes the two files in this directory).

---

## spinosin.pdbqt

- **SMILES (isomeric, PubChem CID 155692):**
  `COC1=C(C(=C2C(=C1)OC(=CC2=O)C3=CC=C(C=C3)O)O)[C@H]4[C@@H]([C@H]([C@@H]([C@H](O4)CO)O)O)O[C@H]5[C@@H]([C@H]([C@@H]([C@H](O5)CO)O)O)O`
- **Source:** `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/155692/property/IsomericSMILES/TXT`
  (fetched twice, both requests returned the same string). CID 155692 is the
  `pubchemCID` recorded for Spinosin at `D:\BALAKATDBV2\src\data\compounds.js:153`.
- **Molecular formula (RDKit-computed from the embedded 3D structure):** C28H32O15,
  matching the formula in compounds.js.
- **Heavy atom count:** 43 (RDKit `mol.GetNumAtoms()` pre-AddHs; also matches the
  `PUBCHEM_HEAVY_ATOM_COUNT` field of the existing `D:\BALAKATDBV2\public\mol\cid_155692.sdf`).
- **Embedding:** ETKDGv3 succeeded on the first attempt (seed `0xF00D`); MMFF94
  optimization converged (`MMFFOptimizeMolecule` returned 0, i.e. no strain
  warning). Geometry not flagged as strained.
- **Active torsions (TORSDOF in the output PDBQT):** 16.
- **SHA-256:** `af4133c660d57773eec0caac682b68567e28009231112438d722908579fa3a6d`
- **Date generated:** 2026-08-03.

### Investigation notes on Spinosin's docking origin

- `compounds.js` records `pubchemCID: 155692`, `bindingAffinity: -7.859`, status
  `'screened'`.
- No file named anything like `spinosin.*` (.pdbqt/.sdf/.mol2/.smi/.pdb) exists
  anywhere under `D:\BALAKATDBV2`, `D:\THESIS_VSC\BALAKATDBV2` (including the
  read-only `talanaidock\dist` tree), or `D:\THESIS_VSC\BALAKATDBV2\ziziphus-docking`
  (the Hugging Face Space source checkout). The live HF Space file tree
  (`huggingface.co/api/spaces/aokenji/ziziphus-docking/tree/main/docking_data/ligands`,
  fetched during this session) also lists only the same 15 ligands already on
  disk; no `spinosin.pdbqt` there either.
- `docking_assets/docking_data/compounds.tsv` in the read-only TalanaiDock dist
  app lists exactly 15 compounds (matching the 15 `.pdbqt` files already in
  `ligands/`) and Spinosin is not one of them. `SOURCES.md` in the same folder
  documents the same 15-compound set. Ursolic Acid **is** in this set
  (`ursolic.pdbqt` exists), so it was prepared through this local pipeline;
  Spinosin was not.
- `D:\BALAKATDBV2\src\data\interactions.js:11-13` carries a code comment: "Ursolic
  Acid and Spinosin come from the TalanaiDock Vina + PLIP pipeline (July 2026)
  ... Contact COUNTS are comparable across batches; bond STRENGTH labels are
  not, because the two runs report different distance conventions." This is a
  direct, specific claim in the codebase (not my inference) that Spinosin was
  docked through TalanaiDock's Vina+PLIP pipeline, separately and later than the
  original 15-compound batch.
- However, the local `talanaidock\dist` app is a shallow git checkout frozen at
  commit `053de17` (2026-06-18, confirmed via `git log -1 --format=%cI`), and
  the local `ziziphus-docking` HF Space checkout's ligand files are dated
  2026-06-12 to 06-14, both **before** the "July 2026" date in the interactions.js
  comment. Neither snapshot, nor the live HF Space, contains a Spinosin ligand
  file or compounds.tsv entry.
- **Verdict:** the docking almost certainly ran through the TalanaiDock Vina+PLIP
  pipeline per the explicit code comment, using a newer version of TalanaiDock
  or the HF Space than any copy present on this machine (or an ephemeral
  BYO-ligand run whose intermediate `.pdbqt` was never persisted back to any
  repo I could read). I could not locate the actual ligand file or job artifact
  used for that run on this machine or on the live HF Space, so the exact SMILES
  and prep parameters used for the original -7.859 kcal/mol result remain
  **unverified** (the codebase's own claim is the best evidence available; I am
  not able to independently confirm it from a surviving file). The
  `spinosin.pdbqt` produced here is a fresh, from-scratch preparation for the
  validation run, not a recovery of whatever file produced the original score.

---

## glc_reembedded.pdbqt

- **SMILES (alpha-D-glucopyranose):** `OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O`
  (task-supplied SMILES).
- **Verification against PubChem CID 79025:** PubChem's IsomericSMILES for CID
  79025 is `C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O)O)O)O)O`, written in a
  different atom order. RDKit canonicalization of both strings produces the
  identical canonical SMILES (`OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O`) and
  identical InChI
  (`InChI=1S/C6H12O6/c7-1-2-3(8)4(9)5(10)6(11)12-2/h2-11H,1H2/t2-,3-,4+,5-,6+/m1/s1`).
  **The task-given SMILES is confirmed correct for alpha-D-glucopyranose /
  CID 79025.** No disagreement.
- **No crystal coordinates used:** built purely from the SMILES string through
  ETKDGv3 embedding + MMFF94 optimization, independent of the 3A4A `GLC 601`
  crystallographic coordinates.
- **Heavy atom count:** 12 (6 C + 6 O), as expected for C6H12O6.
- **Embedding:** ETKDGv3 succeeded (seed `0xF00D`); MMFF94 converged (return 0).
- **Active torsions (TORSDOF):** 6.
- **SHA-256:** `d28c32a3e6f59ffcd33d42b0cf5256378e719f90a40c5b20ebb67bf87ad3968d`
- **Date generated:** 2026-08-03.

### Torsion comparison against the crystal-derived `glc_ligand.pdbqt`

Read from
`D:\THESIS_VSC\BALAKATDBV2\talanaidock\dist\TalanaiDock\app\docking_assets\docking_data\validation\glc_ligand.pdbqt`
(read-only, not modified).

- **TORSDOF in the crystal file: 6.** Same numeric torsion count as the
  re-embedded file, 6.
- **The review's claim, checked branch by branch:** in PDBQT/AutoDock BRANCH
  semantics, a torsion's rotation axis is the bond connecting the parent atom to
  the branch's first (root) atom; that root atom sits on the axis and does not
  move, only atoms further downstream in the branch sweep through space. Walking
  the six `BRANCH` records in `glc_ligand.pdbqt`:
  - `BRANCH 1 7` (ring C1 to hydroxyl O): only the downstream H (atom 8) moves.
  - `BRANCH 2 9` (ring C2 to hydroxyl O): only H (atom 10) moves.
  - `BRANCH 3 11` (ring C3 to hydroxyl O): only H (atom 12) moves.
  - `BRANCH 4 16` (ring C4 to hydroxyl O): only H (atom 17) moves.
  - `BRANCH 5 13` (ring C5 to exocyclic C6): the C6 carbon (atom 13) is on-axis
    and does not move, but the nested branch hanging off it (O6 = atom 14, plus
    its H) is off-axis and does move. **This is the one torsion that displaces
    a heavy atom (O6).**
  - Nested `BRANCH 13 14` (C6 to O6): only the terminal H (atom 15) moves; O6
    itself is on-axis for this inner rotation.
  - That is 6 BRANCH records total, 5 of which (C1-OH, C2-OH, C3-OH, C4-OH, and
    the inner C6-O6-H) rotate only a hydroxyl hydrogen with the oxygen fixed on
    the axis, and exactly 1 (the outer C5-C6 bond) moves a heavy atom (O6).
- **Verdict: the review's claim is TRUE**, verified directly from the file's
  BRANCH/TORSDOF records. TORSDOF is 6, and only the C5-C6 (O6) torsion
  repositions a heavy atom; the other five only reorient a terminal hydroxyl
  hydrogen.
- **The freshly re-embedded file has the identical topology** (same rotatable
  bonds: 4 ring-hydroxyl C-OH bonds, the exocyclic C5-C6 bond, and the C6-O6
  bond), because it is the same molecule and RDKit/Meeko assign rotatable bonds
  from connectivity, not coordinates. So the torsion *count* and *which bonds
  are rotatable* are unchanged by re-embedding. What changes is the *starting
  geometry*: the crystal file starts already in the bound pose (Vina's search
  has nothing left to discover, including for the one heavy-atom-moving
  torsion), while `glc_reembedded.pdbqt` starts from a generic gas-phase
  RDKit/MMFF conformation and must actually be translated, rotated, and
  torsionally adjusted by Vina to reach any docked pose. That is the intended
  fix for the "control exercises no search" concern; it is not visible in the
  TORSDOF number itself.

---

## Control 3 prerequisite: maltose cross-docking, investigated 2026-08-03

Task: prepare the flexibility-matched positive control for Control 3 of
`REDOCK-PROTOCOL.md` section 4, using PDB 3AJ7 as the source of the maltose
ligand. Verified before proceeding, as instructed, rather than assumed.

### Verification finding: PDB 3AJ7 does not contain maltose

- Fetched from `https://files.rcsb.org/download/3AJ7.pdb`, saved as
  `3AJ7.pdb` (SHA-256 `c8c70ea215cb84cedfedc1825e26935dba875316c4c04529a09acc84ee1a1647`,
  478629 bytes).
- Header: `HEADER HYDROLASE 26-MAY-10 3AJ7`; `TITLE CRYSTAL STRUCTURE OF
  ISOMALTASE FROM SACCHAROMYCES CEREVISIAE`; resolution **1.30 A** (`REMARK 2`).
  Citation (`JRNL`): K. Yamamoto, H. Miyake, M. Kusunoki, S. Osaki, "Crystal
  structures of isomaltase from Saccharomyces cerevisiae and in complex with
  its competitive inhibitor maltose," FEBS J. 277, 4205 to 4214 (2010), DOI
  `10.1111/j.1742-4658.2010.07810.x`, PMID 20812985. This matches the citation
  given in the task.
- Direct census of every `HET`/`HETNAM`/`FORMUL` record and every unique
  `HETATM` residue name in the file: exactly two heteroatom species, `CA`
  (1 calcium ion) and `HOH` (608 waters). No sugar of any kind, no maltose, no
  glucose, nothing else.
- Per the hard rule in the task ("if 3AJ7 does not contain maltose, STOP and
  report what it does contain; do not substitute a different structure or
  ligand"), the extraction pipeline stops here. **No `maltose_3aj7_ref.pdb` was
  created**, because there is nothing in 3AJ7 to extract. Producing that file
  anyway would mean inventing coordinates, which "never fabricate a value"
  forbids.

### What the cited paper's own structures actually contain

The FEBS J 2010 paper's title promises two structures: an apo form and "in
complex with its competitive inhibitor maltose." Its own abstract (fetched
from `https://pubmed.ncbi.nlm.nih.gov/20812985/`) gives the two resolutions as
1.30 A and 1.60 A, and states: "An electron density corresponding to a
nonreducing end glucose residue was observed in the active site of isomaltase
in complex with maltose; however, only incomplete density was observed for the
reducing end." That is the paper itself describing partial hydrolysis or
disorder of the soaked maltose during crystallization.

Matching resolutions directly against the deposited entries (both fetched and
inspected here, not assumed from the abstract):

- **3AJ7**: 1.30 A, apo. Matches the paper's apo branch.
- **3A4A**: 1.60 A, contains one `HET GLC 601` (alpha-D-glucopyranose, 12 heavy
  atoms, all-atom occupancy 1.00, no alternate conformation). Matches the
  paper's "maltose complex" branch and the "nonreducing-end glucose only"
  description: one sugar ring was ordered enough to deposit; the other half of
  the maltose was never modeled at all (not present even at partial
  occupancy).

So **the entry that is actually the paper's maltose-soaked structure is 3A4A,
not 3AJ7**, and even 3A4A does not contain intact maltose: it contains the same
12-heavy-atom monosaccharide already used as this project's Control 1/2
glucose ligand (`GLC 601`, the source of `glc_ligand.pdbqt`). The task's
framing of 3AJ7 as "the 3A4A companion structure... isomaltase in complex with
maltose" has the two PDB IDs' roles reversed, and even the corrected pairing
(3A4A as the maltose-soaked entry) does not yield an intact maltose ligand.

### Broader search: no isomaltase structure in the PDB contains intact maltose

Queried RCSB for every entry of this exact protein (UniProt `P53051`,
oligo-1,6-glucosidase / isomaltase, S. cerevisiae): exactly 5 entries exist,
`3A47`, `3A4A`, `3AJ7`, `3AXH`, `3AXI`. Each was fetched and its
`HET`/`HETNAM`/`FORMUL`/`LINK` records inspected directly:

| PDB | Resolution | Paper | Ligand(s) besides Ca2+/water | Disaccharide intact? |
|---|---|---|---|---|
| 3AJ7 | 1.30 A | Yamamoto 2010 FEBS J | none | no ligand at all |
| 3A4A | 1.60 A | Yamamoto 2010 FEBS J | `GLC 601` (12 heavy atoms) | no, single ring only (this is the paper's "maltose complex," but only the nonreducing-end glucose was ordered) |
| 3A47 | 1.59 A | unpublished ("TO BE PUBLISHED") | none | no ligand at all |
| 3AXI | 1.40 A | Yamamoto 2011 J. Biosci. Bioeng. | `GLC 601` (12 heavy atoms) | no, single ring only, despite the entry title "in complex with maltose" |
| 3AXH | 1.80 A | Yamamoto 2011 J. Biosci. Bioeng. | `GLC B 1` + `GLC B 2` (23 heavy atoms total), linked `O6(B1)-C1(B2)` at 1.74 A via a `LINK` record | **yes, but this is isomaltose (alpha-1,6 linkage), not maltose (alpha-1,4 linkage)** |

**Conclusion: no deposited structure of this enzyme contains an intact,
fully-ordered maltose ligand anywhere in the PDB.** The only entry with a
genuine, fully-resolved disaccharide bound is 3AXH, and it is isomaltose, a
constitutional isomer of maltose with a different glycosidic linkage
regiochemistry (alpha-1,6 vs alpha-1,4), from a different, later paper
(Yamamoto, Miyake, Kusunoki, Osaki, "Steric hindrance by 2 amino acid residues
determines the substrate specificity of isomaltase from Saccharomyces
cerevisiae," J. Biosci. Bioeng. 112, 545 to 550 (2011), PMID 21925939, DOI
`10.1016/j.jbiosc.2011.08.016`), not the FEBS J 2010 paper the task cites for
3AJ7.

Per the hard rule, 3AXH's isomaltose was **not** used to build a stand-in
`maltose_3aj7_ref.pdb` and is **not** offered here as Control 3. It is reported
only as an existing, verified fact about what the PDB contains for this
enzyme, for whoever revises the protocol.

Raw files saved, all fetched directly from `files.rcsb.org`, unmodified:

| File | Source URL | SHA-256 | Bytes |
|---|---|---|---|
| `3AJ7.pdb` | `https://files.rcsb.org/download/3AJ7.pdb` | `c8c70ea215cb84cedfedc1825e26935dba875316c4c04529a09acc84ee1a1647` | 478629 |
| `3A4A.pdb` | `https://files.rcsb.org/download/3A4A.pdb` | `4a7c75be6d33eee3360b1ba455555cd21c312117c6e11ef8eb43c7544f6f2eca` | 477738 |
| `3AXI.pdb` | `https://files.rcsb.org/download/3AXI.pdb` | `38870a291a449748b5fc6ebaada772f2daee46b64839d364d8a4bbbd97f46b64` | 480330 |
| `3AXH.pdb` | `https://files.rcsb.org/download/3AXH.pdb` | `c6c568299e4f80fd25ce2227e173a09bf40cd8371775b16d143e453557ed9e6f` | 466398 |
| `3A47.pdb` | `https://files.rcsb.org/download/3A47.pdb` | `d1bbcbdf26a2e3362ec24bacdd74c191a5b72ff0aa20b6701bb98f88279dd61e` | 468828 |

`3A4A.pdb` here is a fresh RCSB fetch used only to inspect its ligand content
for this investigation; it was not used to modify anything in the read-only
`docking_assets` tree under `D:\BALAKATDBV2` or the TalanaiDock dist app, and
nothing in that tree was written to.

### `maltose.pdbqt`: prepared, but explicitly unvalidated against any crystal reference

Built from SMILES through the identical route in `_prep_assets_local.py`'s
`prepare_one()`, same random seed `0xF00D`:

```python
mol = Chem.MolFromSmiles(smiles)
n_heavy = mol.GetNumAtoms()              # heavy-atom count, pre-AddHs
mol = Chem.AddHs(mol)
params = AllChem.ETKDGv3()
params.randomSeed = 0xF00D
params.maxIterations = 2000
AllChem.EmbedMolecule(mol, params)               # returned 0: success
AllChem.MMFFOptimizeMolecule(mol, maxIters=400)  # returned 0: converged
setups = MoleculePreparation().prepare(mol)
pdbqt, ok, err = PDBQTWriterLegacy.write_string(setups[0])
```

- **SMILES (isomeric, PubChem CID 6255, maltose):**
  `C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@@H]2[C@H](O[C@H]([C@@H]([C@H]2O)O)O)CO)O)O)O)O`
- **Source:** `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/6255/property/IsomericSMILES/TXT`
  (fetched twice, both requests returned the identical string).
- **Cross-check:** PubChem `MolecularFormula` for CID 6255 = `C12H22O11`,
  matching maltose; PubChem `InChI`
  (`InChI=1S/C12H22O11/c13-1-3-5(15)6(16)9(19)12(22-3)23-10-4(2-14)21-11(20)8(18)7(10)17/h3-20H,1-2H2/t3-,4-,5-,6+,7-,8-,9-,10-,11-,12-/m1/s1`)
  and `IUPACName` are both consistent with a glucopyranosyl-glucose
  disaccharide.
- **Heavy atom count:** 23 (RDKit `mol.GetNumAtoms()` pre-AddHs), matching the
  expected count for C12H22O11 (12 C + 11 O = 23).
- **Embedding:** ETKDGv3 succeeded on the first attempt (seed `0xF00D`); MMFF94
  optimization converged (`MMFFOptimizeMolecule` returned 0).
- **Active torsions (TORSDOF):** 12.
- **SHA-256:** `ee7f5148a73dd17ae8e64676da6ccd75d1637397702b0aefe51f9c79d5ef83d2`
- **Bytes:** 3234
- **Date generated:** 2026-08-03.

**Flexibility comparison against glucose:** TORSDOF is 12 for maltose versus 6
for both `glc_ligand.pdbqt` (crystal-derived) and `glc_reembedded.pdbqt`
(documented above), so maltose carries exactly double the rotatable-bond count
of glucose in this pipeline's accounting, and its two rings are additionally
joined by a rotatable glycosidic bond that has no counterpart in a
monosaccharide. That is the flexibility gap Control 3 is meant to exercise.

**Why this file alone does not satisfy Control 3.** Control 3 as written
requires cross-docking this ligand into 3A4A and comparing the rank-1 pose
against maltose's crystallographic position in 3AJ7, gated at RMSD under
2.0 A. That crystallographic ground truth does not exist: not in 3AJ7 (no
ligand at all), and not anywhere else for this enzyme (the only genuine
disaccharide-bound structure, 3AXH, is isomaltose, a different molecule).
`maltose.pdbqt` can still be docked into the 3A4A box and scored, but there is
no pose to validate it against, so no RMSD gate can be computed and Control
3's go/no-go verdict cannot be produced from it alone. It is provided as ready
material in case the protocol is revised (for example, to accept isomaltose
from 3AXH as the flexibility-matched disaccharide control instead of maltose,
or to source a maltose-bound structure from a homologous GH13 enzyme), not as
a working Control 3 input.

### Cross-docking geometry check: performed on what is actually available

Since 3AJ7 contains no maltose to place in the 3A4A box, the check specified in
the task cannot be run for maltose itself. Two adjacent, genuinely informative
checks were run instead, using only `numpy` and standard-library PDB parsing
(no other installs):

**1. Are 3AJ7 and 3A4A, the two real structures from the cited paper, already
in the same crystallographic frame?** Both are isomorphous crystals, space
group C2 (`C 1 2 1`), with very similar unit cells (3AJ7: 95.516, 115.538,
61.756 A, beta 91.05 deg; 3A4A: 95.428, 115.404, 61.606 A, beta 91.19 deg).
Direct, unaligned comparison of the three catalytic-residue CA positions gives:

| Residue | 3AJ7 to 3A4A distance, no transformation applied |
|---|---|
| Asp215 | 0.441 A |
| Glu277 | 0.430 A |
| Asp352 | 0.551 A |

A full Kabsch superposition on all 586 shared backbone CA atoms gives RMSD
**0.1103 A**. Verdict: yes, 3AJ7 and 3A4A are already effectively in the same
frame; no meaningful superposition is needed between these two specifically.

**2. Bonus/reference only, not a maltose control: where would 3AXH's
isomaltose fall if cross-docked into 3A4A's box?** 3AXH does **not** share
3A4A's raw frame: its catalytic-residue CA coordinates are approximately
mirrored (for example, Asp215 CA sits at `[-19.31, -2.393, -27.159]` in 3AXH
versus `[19.425, -3.039, 27.246]` in 3A4A), consistent with a different
origin/axis choice within the same space group rather than a different fold.
A full Kabsch superposition on 586 shared CA atoms gives backbone RMSD
**0.1202 A** (an excellent fit once properly aligned). Under that
superposition, the isomaltose centroid (23 atoms, chain B, residues 1 and 2)
maps to `[19.841, -8.264, 21.544]` in 3A4A's frame. The existing 3A4A box is
centered at `[21.52, -7.70, 23.55]` with 30 A sides (15 A half-width per
axis):

- Offset of the isomaltose centroid from the box center:
  `[-1.679, -0.564, -2.006]`, magnitude **2.676 A**, well inside the box on
  all three axes.
- Distances from this centroid to the catalytic triad in 3A4A: Asp215
  7.744 A, Glu277 8.821 A, Asp352 6.454 A, consistent with active-site
  binding.
- Without superposition (raw, unaligned coordinates), the same centroid sits
  **61.065 A** from the box center: nowhere near it.

This shows that if the protocol were revised to substitute isomaltose for
maltose, no box change would be needed, but the coordinate-frame mismatch
between 3AXH and 3A4A is real, and a superposition step would be mandatory
first; skipping it would silently place the ligand about 60 A from where it
belongs.

### Recommendation for the protocol authors, not acted on here

Control 3 as specified cannot run: its named source structure (3AJ7) has no
ligand, and no isomaltase structure in the PDB has intact maltose. Options to
consider, none chosen here:

1. Substitute 3AXH's isomaltose as the flexibility-matched disaccharide
   control, accepting that it tests a different linkage (alpha-1,6) than the
   compounds actually screened.
2. Find a maltose-bound structure of a homologous GH13 alpha-glucosidase from
   a different organism and use that as the crystallographic reference,
   accepting cross-species RMSD as a weaker comparison.
3. Drop the crystallographic RMSD gate for maltose specifically and substitute
   a different pass criterion (for example, active-site contact or pose
   plausibility only, no RMSD threshold).
4. Treat the review's original concern as adequately addressed by Control 2
   alone and drop Control 3.

This is a finding to route back through the same review process that produced
revision 2 of `REDOCK-PROTOCOL.md`, not a decision made in this session.

---

## Files in this directory

| File | SHA-256 | Heavy atoms | TORSDOF | Notes |
|---|---|---|---|---|
| `spinosin.pdbqt` | `af4133c660d57773eec0caac682b68567e28009231112438d722908579fa3a6d` | 43 | 16 | |
| `glc_reembedded.pdbqt` | `d28c32a3e6f59ffcd33d42b0cf5256378e719f90a40c5b20ebb67bf87ad3968d` | 12 | 6 | |
| `maltose.pdbqt` | `ee7f5148a73dd17ae8e64676da6ccd75d1637397702b0aefe51f9c79d5ef83d2` | 23 | 12 | prepared 2026-08-03; NOT validated against any crystal reference, no maltose-bound structure of this enzyme exists (see above) |
| `3AJ7.pdb` | `c8c70ea215cb84cedfedc1825e26935dba875316c4c04529a09acc84ee1a1647` | n/a (raw structure) | n/a | raw RCSB fetch; apo, contains no maltose (the verification finding) |
| `3A4A.pdb` | `4a7c75be6d33eee3360b1ba455555cd21c312117c6e11ef8eb43c7544f6f2eca` | n/a (raw structure) | n/a | raw RCSB fetch; this is the paper's actual "maltose complex" entry, but contains only `GLC 601` (12 atoms); fetched fresh for inspection only, the read-only receptor asset tree was not touched |
| `3AXI.pdb` | `38870a291a449748b5fc6ebaada772f2daee46b64839d364d8a4bbbd97f46b64` | n/a (raw structure) | n/a | raw RCSB fetch; titled "in complex with maltose" but contains only `GLC 601` (12 atoms), from a different (2011) paper |
| `3AXH.pdb` | `c6c568299e4f80fd25ce2227e173a09bf40cd8371775b16d143e453557ed9e6f` | n/a (raw structure) | n/a | raw RCSB fetch; contains isomaltose (`GLC`-`GLC`, 23 atoms, alpha-1,6 linked), the only intact disaccharide bound to this enzyme anywhere in the PDB; this is NOT maltose |
| `3A47.pdb` | `d1bbcbdf26a2e3362ec24bacdd74c191a5b72ff0aa20b6701bb98f88279dd61e` | n/a (raw structure) | n/a | raw RCSB fetch; apo, unpublished structure |
| `isomaltose_3axh_raw.pdb` | `b488526f558267ed7a9fff1feb39ec0da207a09b593571b6252e2a3ada335398` | 23 | n/a | raw extraction, 3AXH native coordinates, see below |
| `isomaltose_3axh_ref.pdb` | `bcf82d0be8d599b5bf184aa28545c4f34fb5f6f1ca72a4ce52da1eee3e4a3464` | 23 | n/a | Kabsch-superposed into 3A4A frame, the Control 3 RMSD reference target, see below |
| `isomaltose.pdbqt` | `b64d8856a2ea4d49c46453f3d1e0807baf9b76878ef8a51d8e93d7696266699c` | 23 | 12 | prepared 2026-08-03 from SMILES (PubChem CID 439193), see below |

`maltose_3aj7_ref.pdb` was **not created**: 3AJ7 has no maltose ligand to
extract, and no substitute ligand or structure was used in its place, per the
task's hard rule. **Control 3 has since been amended by the protocol owner to
use isomaltose cross-docked from 3AXH instead of maltose from 3AJ7; the
isomaltose files above are that amended control's inputs, documented in full
below.**

---

## Control 3, amended 2026-08-03: isomaltose from 3AXH, cross-docked into 3A4A

`REDOCK-PROTOCOL.md` section 4 was amended by the protocol owner after the
maltose investigation above. Control 3 now uses isomaltose from PDB 3AXH
(already downloaded and inspected during the maltose search) as the
flexibility-matched disaccharide, cross-docked into 3A4A. This section
documents the three files that implement it.

### `isomaltose_3axh_raw.pdb`: extraction and independent verification

Extracted directly from the already-downloaded `3AXH.pdb`
(`https://files.rcsb.org/download/3AXH.pdb`, SHA-256
`c6c568299e4f80fd25ce2227e173a09bf40cd8371775b16d143e453557ed9e6f`), every
`HETATM` record with residue name `GLC` and chain `B`, coordinates
unmodified.

- **Chain and residue numbers:** chain `B`, residues `1` and `2` (two `GLC`
  groups, no shared/ambiguous numbering with the protein chain, which is
  chain `A`).
- **Heavy atom count:** 23 total, 12 on residue 1 (`C1`-`C6`, `O1`-`O6`) and
  11 on residue 2 (`C1`-`C6`, `O2`-`O6`, no `O1`). Matches the expected count
  for C12H22O11.
- **Occupancy and alternate conformations:** every one of the 23 atoms has
  occupancy `1.00` and a blank alternate-location indicator (column 17 in the
  PDB record). No partial occupancy, no altloc `A`/`B` splitting, on either
  residue.
- **The alpha-1,6 linkage, confirmed from coordinates, not from naming or the
  deposited `LINK` record:** the distance between `O6` of `GLC B 1` and `C1`
  of `GLC B 2`, computed directly from their `x, y, z` columns, is **1.744 A**.
  This is a bonded distance (the file's own `LINK` record separately and
  independently states 1.74 A, so the two agree, but the 1.744 A value here
  was computed from the coordinates themselves, not copied from that record).
  Residue 1 retains a free `O1` (an unlinked anomeric hydroxyl, the reducing
  end of the disaccharide), while residue 2 has no `O1` at all (its `C1`, the
  anomeric carbon, is instead the atom bonded through the glycosidic oxygen to
  residue 1's `O6`). A glucose whose C1 bonds through an oxygen to another
  glucose's O6 is, by definition, an alpha-1,6 linkage: this is isomaltose,
  verified from the atom connectivity itself rather than assumed from the
  `GLC`/`GLC` residue naming.
- **SHA-256:** `b488526f558267ed7a9fff1feb39ec0da207a09b593571b6252e2a3ada335398`
- **Bytes:** 2448
- **Date generated:** 2026-08-03.

### `isomaltose_3axh_ref.pdb`: superposed into 3A4A's frame, the RMSD reference

Every RMSD Control 3 reports will be measured against this file, so it carries
the full superposition provenance as `REMARK` lines in the file itself, not
only here.

**Superposition.** A Kabsch alignment of the protein backbone CA atoms was
computed between 3AXH and 3A4A (chain A, matching residue numbers, no
alternate-location atoms other than blank/`A`): **586 shared CA atoms**,
backbone RMSD after superposition **0.1202 A**. (Note: this file uses the
3AXH-to-3A4A transform, 0.1202 A. Do not confuse it with the separate
3AJ7-to-3A4A backbone RMSD of 0.1103 A reported earlier in this document,
which answers a different question, whether 3AJ7 and 3A4A share a frame, and
is not used anywhere in building this reference file.)

The resulting rotation and translation were applied to all 23 atoms of
`isomaltose_3axh_raw.pdb`. This was re-run independently for this deliverable
(not copied from the earlier investigation's numbers) and reproduced the same
586-atom count and the same 0.1202 A backbone RMSD, confirming the transform
is deterministic and stable.

- **Ligand centroid in 3A4A's frame:** `[19.841, -8.264, 21.544]`.
- **Offset from the 3A4A box center** (`21.52, -7.70, 23.55`, 30 A cube, 15 A
  half-width per axis), independently re-verified: `[-1.679, -0.564, -2.006]`,
  magnitude **2.676 A**. This reproduces the 2.676 A figure from the earlier
  investigation exactly. Inside the box on all three axes.
- **Per-residue minimum heavy-atom distance** (every isomaltose heavy atom
  against every heavy atom of the named 3A4A residue, chain A; this replaces
  the earlier centroid-to-CA range of 6.5 to 8.8 A with the actual closest
  contact per residue):

  | Catalytic residue | Minimum heavy-atom distance | Closest atom pair |
  |---|---|---|
  | Asp215 | **2.724 A** | receptor `OD2` to isomaltose `GLC B 2` atom `O6` |
  | Glu277 | **3.045 A** | receptor `OE2` to isomaltose `GLC B 2` atom `O2` |
  | Asp352 | **2.379 A** | receptor `OD2` to isomaltose `GLC B 2` atom `O2` |

  All three minimum distances are in a plausible hydrogen-bonding range
  (2.4 to 3.0 A), consistent with the isomaltose sitting in genuine contact
  with the catalytic triad in 3A4A's frame, not merely somewhere inside the
  box.
- **REMARK lines in the file itself** record: the number of backbone CA atoms
  used (586), the backbone RMSD (0.1202 A), the ligand centroid in 3A4A
  coordinates, the box center and offset, and the three per-residue minimum
  distances above, so the file is self-describing.
- **SHA-256:** `bcf82d0be8d599b5bf184aa28545c4f34fb5f6f1ca72a4ce52da1eee3e4a3464`
- **Bytes:** 2885
- **Date generated:** 2026-08-03.
- This file is a reference/target for RMSD scoring only. It is not meant to be
  docked itself.

### `isomaltose.pdbqt`: built from SMILES, the actual Control 3 ligand input

Built through the identical `_prep_assets_local.py` route used for every
other ligand in this directory, same seed `0xF00D`:

```python
mol = Chem.MolFromSmiles(smiles)
n_heavy = mol.GetNumAtoms()
mol = Chem.AddHs(mol)
params = AllChem.ETKDGv3()
params.randomSeed = 0xF00D
params.maxIterations = 2000
AllChem.EmbedMolecule(mol, params)               # returned 0: success
AllChem.MMFFOptimizeMolecule(mol, maxIters=400)  # returned 0: converged
setups = MoleculePreparation().prepare(mol)
pdbqt, ok, err = PDBQTWriterLegacy.write_string(setups[0])
```

- **SMILES (isomeric, PubChem CID 439193, isomaltose), confirmed:**
  `C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)OC[C@@H]2[C@H]([C@@H]([C@H](C(O2)O)O)O)O)O)O)O)O`
- **Source:** `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/439193/property/IsomericSMILES/TXT`
  (fetched twice, identical both times).
- **CID 439193 confirmed as isomaltose:** PubChem synonyms for this CID
  include `ISOMALTOSE`, CAS `499-40-1`, and `6-O-alpha-D-glucopyranosyl-D-glucopyranose`
  (the "6-O" naming is the alpha-1,6 linkage). `MolecularFormula` = `C12H22O11`.
  `InChI` = `InChI=1S/C12H22O11/c13-1-3-5(14)8(17)10(19)12(23-3)21-2-4-6(15)7(16)9(18)11(20)22-4/h3-20H,1-2H2/t3-,4-,5-,6-,7+,8+,9-,10-,11?,12+/m1/s1`.
- **Confirmed NOT maltose:** for contrast, PubChem CID 6255 (maltose)
  synonyms include `4-O-alpha-D-glucopyranosyl-` and
  `alpha-D-glucopyranosyl-(1->4)-beta-D-glucopyranose` (the "4-O" / "1->4"
  naming is the alpha-1,4 linkage). RDKit was used to parse both SMILES
  independently and compare: same molecular formula (`C12H22O11`) but
  **different canonical SMILES and different InChI** (the connectivity layer
  differs at the glycosidic bridge: isomaltose's InChI routes the bridge
  through an exocyclic CH2 carbon, matching a 6-position/primary-carbon
  linkage, while maltose's routes it directly to a ring carbon, matching a
  4-position/ring linkage). The two molecules are confirmed structurally
  distinct, and the isomaltose SMILES used here reproduces PubChem's own
  InChI for CID 439193 exactly when passed through RDKit.
- **Heavy atom count:** 23 (RDKit `mol.GetNumAtoms()` pre-AddHs), matching
  C12H22O11 (12 C + 11 O).
- **Ring check:** RDKit's ring perception on the parsed (pre-embedding)
  molecule finds exactly two 6-membered rings (two pyranose rings), as
  expected for a glucose-glucose disaccharide.
- **Embedding:** ETKDGv3 succeeded on the first attempt (seed `0xF00D`); MMFF94
  optimization converged (`MMFFOptimizeMolecule` returned 0).
- **Active torsions (TORSDOF):** 12.
- **SHA-256:** `b64d8856a2ea4d49c46453f3d1e0807baf9b76878ef8a51d8e93d7696266699c`
- **Bytes:** 3231
- **Date generated:** 2026-08-03.

**TORSDOF comparison across the sugar series:**

| Ligand | Heavy atoms | TORSDOF |
|---|---|---|
| glucose (`glc_ligand.pdbqt`, crystal, and `glc_reembedded.pdbqt`, re-embedded) | 12 | 6 |
| maltose (`maltose.pdbqt`, alpha-1,4, not validated, see above) | 23 | 12 |
| **isomaltose (`isomaltose.pdbqt`, alpha-1,6, the actual Control 3 ligand)** | **23** | **12** |

Isomaltose and maltose are constitutional isomers (same formula, same heavy
atom count, one extra glycosidic linkage each versus a monosaccharide) and, in
this pipeline's accounting, carry the identical TORSDOF of 12 versus glucose's
6, exactly double. The flexibility-matching argument made for maltose in the
original task applies identically to isomaltose: this is not a coincidence of
which disaccharide was chosen, it is a property of being a disaccharide built
from two flexible pyranose rings joined by a rotatable glycosidic bond.

### Pre-docking baseline RMSD, using `talanai.control.rmsd`

Computed with `D:\THESIS_VSC\talanai-lang` added to `sys.path` and the
project's own comparator, not a bespoke calculation:

```python
sys.path.insert(0, r"D:\THESIS_VSC\talanai-lang")
from talanai import control
pose = control.read_heavy_atoms("isomaltose.pdbqt")
reference = control.read_heavy_atoms("isomaltose_3axh_ref.pdb")
value, method, matched, exact = control.rmsd(pose, reference)
```

- **Atoms:** 23 in each file, all 23 matched.
- **Method:** `LOWER BOUND, nearest same-element atom, correspondence
  unresolved (atom names differ between pose and reference)`. This is
  `talanai.control`'s case 3, triggered because the Meeko-prepared PDBQT names
  every atom generically (`C`, `O`) while the crystallographic reference keeps
  the PDB chemical-component names (`C1`, `O6`, and so on); the two files also
  do not list atoms in the same order, so neither the by-name nor the
  by-file-order fast paths apply.
- **RMSD (lower bound): 26.406 A.**
- **Exact:** `False` (this is a lower bound on the true correspondence-based
  RMSD, not an exact value; `talanai.control`'s own documentation is explicit
  that a lower bound is only conclusive when it exceeds a threshold, which is
  not the use made of it here).

This is the intended result, not an error: `isomaltose.pdbqt` is a generic
MMFF94 conformer with no relationship to the receptor's coordinate frame
(RDKit embeds a molecule around its own centroid, unrelated to any binding
site), while `isomaltose_3axh_ref.pdb` sits at the real, superposed bound-pose
coordinates in 3A4A's frame, roughly 20 A from the origin. The two are
expected to be far apart before any docking search runs. This number is the
**pre-docking baseline**: once Control 3 is actually run (Vina cross-docking
`isomaltose.pdbqt` into 3A4A), the rank-1 pose's RMSD against this same
`isomaltose_3axh_ref.pdb` should drop sharply, and the size of that drop is the
evidence that the search, not the starting geometry, is what produces the
docked pose. A rank-1 RMSD that is not dramatically smaller than 26.406 A
would itself be a finding worth reporting.

### Why isomaltose and not maltose

*(Written for the thesis methods section; every figure below was independently
verified in this session, not taken from a secondary source.)*

The original protocol specified maltose from PDB 3AJ7 as the flexibility-
matched positive control, cross-docked into 3A4A, on the premise that 3AJ7 is
"the 3A4A companion structure... isomaltase in complex with maltose" from
Yamamoto, Miyake, Kusunoki and Osaki, *FEBS Journal* 277, 4205 to 4214 (2010),
DOI `10.1111/j.1742-4658.2010.07810.x`, PMID 20812985.

Direct inspection of 3AJ7 (fetched from RCSB, `HET`/`HETATM` records checked
by hand) showed it contains **no ligand at all** besides one calcium ion and
608 waters. This alone falsified the premise and, per this project's rule
against substituting a different structure or ligand without saying so, the
extraction pipeline stopped rather than silently continuing.

A five-entry survey followed. Every PDB structure of this exact protein
(UniProt `P53051`, oligo-1,6-glucosidase / isomaltase, EC 3.2.1.10,
*Saccharomyces cerevisiae*) was retrieved and its ligand content checked
directly:

| PDB | Resolution | Paper | Ligand(s) besides Ca2+/water | Disaccharide intact? |
|---|---|---|---|---|
| 3AJ7 | 1.30 A | Yamamoto 2010, *FEBS J.* 277:4205 | none | no ligand at all |
| 3A4A | 1.60 A | Yamamoto 2010, *FEBS J.* 277:4205 | `GLC 601` (12 heavy atoms) | no, single ring only |
| 3A47 | 1.59 A | unpublished | none | no ligand at all |
| 3AXI | 1.40 A | Yamamoto 2011, *J. Biosci. Bioeng.* 112:545 | `GLC 601` (12 heavy atoms) | no, single ring only |
| 3AXH | 1.80 A | Yamamoto 2011, *J. Biosci. Bioeng.* 112:545 | `GLC B 1` + `GLC B 2` (23 heavy atoms, alpha-1,6 linked) | **yes: isomaltose** |

**The Yamamoto 2010 density explanation.** The FEBS J 2010 paper's own
abstract states that a maltose-soaked crystal showed "an electron density
corresponding to a nonreducing end glucose residue... in the active site...
however, only incomplete density was observed for the reducing end," which
matches 3A4A exactly (1.60 A, a single ordered `GLC 601`, the second glucose
of the soaked maltose never modeled at all, not even at partial occupancy).
So the paper's own maltose-soaked structure never actually contained intact,
fully-ordered maltose; only a single hydrolysis or disorder product, glucose,
survived to be deposited. No structure of this enzyme, across either paper,
contains intact maltose.

**The 3AXI title trap.** 3AXI's own RCSB title reads "Crystal structure of
isomaltase in complex with maltose," which reads as exactly the structure the
original protocol wanted. Direct inspection shows it, too, contains only a
single `GLC 601` (12 heavy atoms), not an intact disaccharide: the title
describes the soaking experiment, not what the crystal preserved. Taking PDB
entry titles at face value here would have produced a control built on a
monosaccharide while believing it to be a flexible disaccharide, silently
reintroducing the exact "rigid 12-atom sugar" problem Control 3 exists to
avoid. Only 3AXH, titled "in complex with isomaltose," turned out to contain a
genuine, fully-ordered, alpha-1,6-linked 23-atom disaccharide, independently
confirmed from its coordinates (the 1.744 A `O6`-`C1` bond and the missing
`O1` on the linked residue, both above), not from its title or residue naming.

**The cognate-substrate argument.** Isomaltase is formally oligo-1,6-
glucosidase (EC 3.2.1.10): its designated activity is hydrolyzing alpha-1,6
glycosidic bonds. Isomaltose (alpha-D-glucopyranosyl-(1->6)-D-glucose) is
exactly that bond, making it the enzyme's cognate substrate. Maltose
(alpha-D-glucopyranosyl-(1->4)-beta-D-glucopyranose, PubChem CID 6255) carries
the alpha-1,4 bond instead, and the FEBS J 2010 paper itself describes maltose
as isomaltase's "competitive inhibitor," not its substrate, consistent with
alpha-1,4 being the non-preferred linkage for this enzyme. Validating the
docking protocol's pose recovery on isomaltose therefore tests the protocol
against a ligand the receptor is naturally built to bind and process, which is
a better-justified positive control than validating it against a linkage the
enzyme merely tolerates as an inhibitor, independent of the fact that
isomaltose is also the only one of the two disaccharides actually available,
intact, anywhere in the PDB for this protein.

Citations:

- K. Yamamoto, H. Miyake, M. Kusunoki, S. Osaki. "Crystal structures of
  isomaltase from Saccharomyces cerevisiae and in complex with its
  competitive inhibitor maltose." *FEBS Journal* 277(19), 4205 to 4214 (2010).
  DOI `10.1111/j.1742-4658.2010.07810.x`. PMID 20812985. PDB entries 3AJ7,
  3A4A (and 3A47, unpublished, same UniProt entry).
- K. Yamamoto, H. Miyake, M. Kusunoki, S. Osaki. "Steric hindrance by 2 amino
  acid residues determines the substrate specificity of isomaltase from
  Saccharomyces cerevisiae." *Journal of Bioscience and Bioengineering*
  112(6), 545 to 550 (2011). DOI `10.1016/j.jbiosc.2011.08.016`.
  PMID 21925939. PDB entries 3AXH, 3AXI.
