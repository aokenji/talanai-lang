# Run records index: the receptor-preparation finding

This file exists because rule R105 ("a live endpoint reading is not a citable
artefact; someone else has to be able to find the run this number came from")
is violated by the project's own central finding otherwise. The receptor-
preparation conclusion in BUILD-SPEC.md section 12 and REDOCK-PROTOCOL.md
section 1 rests on four redocking cases and a six-ligand re-dock. Only one of
those five runs has anything on disk: this directory's
`redock_screening_seed42.json`, for case "A screening", and even that record's
`rank1_rmsd_A` value was computed by a method later withdrawn (see that file's
`superseded_by` field, added 2026-08-03; the corrected value for that run is
6.164 A, matching the index below).

**None of the numbers in this file has an independent run record beyond that
one exception.** They are transcribed here, accurately and without alteration,
from BUILD-SPEC.md section 12 and the identical table in REDOCK-PROTOCOL.md
section 1, which is as close to primary as this project currently has for
cases B, C and D and for the six-ligand re-dock. Every one of them must be
regenerated with `tal control` (which now writes `rmsd_method` and
`talanai_version` into every record it produces, see `control.py`) before it
is cited in anything that leaves this repository.

---

## The four redocking cases

Fixed across all four: box centre (21.52, -7.70, 23.55); ligand GLC 601 from
`validation\glc_ligand.pdbqt`; reference `validation\glc_3a4a_601.pdb`; Vina
1.2.7; num_modes 9; energy_range 3; cpu 4; seed 42. Only the receptor
preparation, box size and exhaustiveness vary between cases, which is the
whole point of running all four (BUILD-SPEC.md section 12; REDOCK-PROTOCOL.md
section 1).

| Case | Receptor file | Box size | Exhaustiveness | Seed | Score (kcal/mol) | RMSD (A) | RMSD method |
|---|---|---|---|---|---|---|---|
| B calibration | `prepared\receptor.pdbqt` (Meeko) | 18 A | 16 | 42 | -5.885 | 0.513 | LOWER BOUND, nearest same-element atom, correspondence unresolved |
| D search only | `prepared\receptor.pdbqt` (Meeko) | 30 A | 8 | 42 | -5.905 | 0.518 | LOWER BOUND, nearest same-element atom, correspondence unresolved |
| C prep only | `receptor_clean.pdb` (raw) | 18 A | 16 | 42 | -4.208 | 5.592 | LOWER BOUND, nearest same-element atom, correspondence unresolved |
| A screening | `receptor_clean.pdb` (raw) | 30 A | 8 | 42 | -4.024 | 6.164 | LOWER BOUND, nearest same-element atom, correspondence unresolved |

Source: BUILD-SPEC.md section 12, item 1 ("R106, blocking, and now
measured"), cross-checked against the identical table in REDOCK-PROTOCOL.md
section 1. The RMSD method for all four is the lower bound described in both
documents: "the docking input names every atom C or O while the crystal
reference uses C1..C6 / O1..O6", which is `control.py`'s `rmsd()` case 3, the
only one of its three correspondence strategies that applies to this ligand
pairing. A lower bound above the 2.0 A threshold is conclusive for failure; a
lower bound below it is not by itself a pass (`control.py` docstring, and
`rmsd()`'s own comments).

**Run record status: NOT CAPTURED AT THE TIME.** Only case A ("screening") has
a JSON artefact (`redock_screening_seed42.json`, this directory), and that
artefact's `rank1_rmsd_A` field was produced by an atom-pairing method that no
longer exists in `control.py` (see its `superseded_by` field). Cases B, C and
D have no run record of any kind on disk: no command line, no pose file, no
checksum, nothing to audit against. Regenerate all four with `tal control`
before citing any of these four numbers in a publication.

---

## The six-ligand re-dock ("Consequence, exploratory")

Re-docking the thesis screening set on the receptor that passes its own
redocking control (the Meeko-prepared receptor of cases B and D), same box
(21.52, -7.70, 23.55; 30 x 30 x 30 A), same exhaustiveness (8), same seed (42):

| Compound | Thesis (raw) kcal/mol | Prepared kcal/mol | Delta |
|---|---|---|---|
| Rutin | -8.857 | -10.620 | -1.763 |
| Betulinic acid | -8.290 | -9.644 | -1.354 |
| Quercetin | -7.503 | -8.889 | -1.386 |
| Kaempferol | -7.479 | -8.240 | -0.761 |
| Oleanolic acid | -6.922 | -7.792 | -0.870 |
| Acarbose | -6.660 | -8.576 | -1.916 |

Source: BUILD-SPEC.md section 12, item 2 ("Consequence, exploratory"); the
identical table appears in REDOCK-PROTOCOL.md section 1.

**Run record status: NOT CAPTURED AT THE TIME.** No `validation-run\`
directory exists anywhere under this repository, so none of these six scores
has a command line, pose file or checksum behind it. Both source documents
already flag this table as provisional rather than final: "Single seed,
exhaustiveness 8, so this is a signal to investigate, not a number to
publish" (BUILD-SPEC.md section 12, item 2). Do not cite any of these six
numbers in a publication until the three-seed, exhaustiveness-32 re-dock
described in REDOCK-PROTOCOL.md is actually run and each ligand's own run
record is written.

---

## What this file is, and is not

This is an audit index, not a run record. It exists so that a reader who
follows R105 to its logical conclusion lands on an accurate accounting of
what was and was not captured, rather than on silence or on a single stale
artefact standing in for five runs. It does not substitute for regenerating
the missing JSON artefacts, and every number above should be treated as
"reported in the build documents, not independently verifiable from a stored
run" until that regeneration happens.

Compiled 2026-08-03 from BUILD-SPEC.md section 12 and REDOCK-PROTOCOL.md
section 1 only. No value in this file was computed or re-derived; all are
transcribed and cross-checked for consistency between those two documents.
