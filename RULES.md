# The rules

This is the product. Everything else in Talanai exists so this table can run.

Three severities:

- **REFUSE** the experiment does not run. Not a preference: the result would
  not be interpretable.
- **WARN** it runs, and you should know this before you present it.
- **RECORD** neither wrong nor optional. A fact that must appear in the run
  record and the generated methods paragraph.

A check that cannot be performed reports RECORD with the word **UNVERIFIED**.
It never reports a pass it did not earn.

---

## R0xx  the file itself

| ID | Severity | Rule | Reason |
|---|---|---|---|
| R000 | REFUSE | The file could not be read | A setting outside any block belongs to nothing |
| R001 | REFUSE | A required block is missing | Without `study`, `receptor`, `site`, `control`, `ligands` and `dock` the experiment is sketched, not described |
| R002 | WARN | Unknown block | A misspelt block would otherwise be silently ignored, which is worse than an error |
| R003 | WARN | Unknown key | Same reasoning one level down. A typo must never read as a default |
| R004 | REFUSE | A single-value setting is given twice | Two exhaustiveness values means nobody knows which one ran |

## R1xx  the control, which unlocks everything else

| ID | Severity | Rule | Reason |
|---|---|---|---|
| R101 | REFUSE | The redocking control has no recorded result | A protocol that cannot reproduce a known answer has no business predicting an unknown one. This is the positive-control lane on a blot |
| R102 | REFUSE | The control failed its threshold | A failed control is information. The pocket, the preparation or the search is wrong |
| R103 | WARN | The threshold is looser than 2.0 Å | 2.0 Å is the accepted bar. A looser one has to be argued for in text |
| R104 | WARN | The control ran under a different box or a higher exhaustiveness than the screen | Focused-box redocking is standard and asks an easier question than the screen does. Legitimate, but a reader will otherwise assume the RMSD validates the screening search |
| R105 | REFUSE | The control result has no source | A live endpoint reading is not a citable artefact. Someone else has to be able to find the run |
| R106 | REFUSE | The control validated a different receptor preparation | **The one that bites.** A control validates the protocol it ran under. Change the receptor preparation and the thing that passed at 2.0 Å is not the thing that produced the affinities |

## R2xx  words that mean two things

| ID | Severity | Rule | Reason |
|---|---|---|---|
| R201 | REFUSE | Bare `rmsd` | This project uses RMSD for two unrelated quantities: redocking accuracy, and per-compound pose-cluster spread (1.61 to 2.00 Å). Quoting one where the other was meant is a scientific error, not a typo. Write `redock_rmsd` or `pose_cluster_rmsd` |

## R3xx  geometry

| ID | Severity | Rule | Reason |
|---|---|---|---|
| R301 | REFUSE | The box is not fully specified | A box is three coordinates and three lengths |
| R302 | WARN | The box declares no residues to enclose | An unanchored box is the commonest way to produce confident nonsense: the search runs perfectly, in the wrong pocket |
| R303 | REFUSE / RECORD | The box must enclose the named catalytic residues | Verified geometrically when the structure is on disk, UNVERIFIED when it is not |
| R304 | WARN | A named residue is not in the structure | Numbering or chain mismatch |
| R305 | WARN | The box is large enough to be blind docking, undeclared | Blind docking is a legitimate method with a different interpretation and a far larger search requirement. State it |
| R307 | REFUSE / RECORD | A ligand's saturated six-rings must be chairs | Vina holds ring bonds rigid, so whatever ring conformation the preparation produced is what docks, scores and ends up in the results. A boat or twist-boat is strained and usually non-physiological, and docking cannot repair it. Rings containing an sp2 atom are exempt, since a double bond forces a half chair. **2026-08-04: this study's own preparation embedded one conformer and minimised locally, which cannot cross a ring-flip barrier. Three ligands were affected including the reference compound, and correcting the reference alone moved another compound from "loses to acarbose" to "tied with it"** |
| R306 | WARN / RECORD | The rank-1 pose must actually contact the named catalytic residues, not merely share a box with them | R303 checks the box; this checks the pose. UNVERIFIED when no pose file is on disk to check, WARN when the nearest named residue is past the 4.0 Å contact cutoff. A pose that never reaches the catalytic residues does not support an active-site inhibition claim, whatever the box around it encloses |

## R4xx  the search

| ID | Severity | Rule | Reason |
|---|---|---|---|
| R401 | WARN | Search budget is low for the box volume | Vina's default exhaustiveness of 8 was chosen for small boxes. Effort has to scale with the volume searched. Rule of thumb, not a published standard |
| R402 | RECORD | `modes` or `energy_range` not recorded | The pose ensemble cannot be described without them |

## R5xx  the search is stochastic

| ID | Severity | Rule | Reason |
|---|---|---|---|
| R501 | WARN | Single seed | Vina searches; where it starts affects where it lands. Repeating one seed reproduces the answer exactly and says nothing about convergence. A fixed seed is a scale that reads the same number every time and can still be off |

## R6xx  preparation is part of the measurement

| ID | Severity | Rule | Reason |
|---|---|---|---|
| R601 | REFUSE | Ligands and reference were prepared differently | June 2026, this project: TalanaiDock reproduced Quercetin at −7.525 only once Meeko preparation was removed. A Vina score belongs to a pair of prepared files, not to a molecule, so two preparations are two assays |
| R602 | REFUSE | Receptor preparation not declared | Waters, hetatms, hydrogens and charges all move the score. It is the most commonly omitted item in published docking methods |
| R603 | WARN | Ligand protonation state not declared | The structure drawn in a paper is not what exists at pH 7.4 |
| R604 | REFUSE / WARN / PASS | A prepared file does not match its recorded checksum; or none is recorded | September 2026, this project: re-preparing quercetin under the recorded recipe (“rdkit ETKDGv3 + MMFF94 with ring-aware conformer selection”) scored −8.380 against a published −8.818, while the published FILE scored −8.877 under the same receptor, box, seed and exhaustiveness. The engine reproduces; the conformer does not. Rutin, whose saturated rings let ring-awareness constrain the choice, reproduced to 0.013. A recipe alone is a description; a recipe plus a digest is a protocol |

## R7xx  reading the results

| ID | Severity | Rule | Reason |
|---|---|---|---|
| R701 | WARN | Ranking on raw affinity only | Vina scores rise with heavy-atom count, so ranking on the raw number partly ranks by molecular weight |
| R702 | RECORD | The strongest binder is also the heaviest compound | The size confound is live and a reviewer will say so. Recording it is the answer prepared in advance |

## R8xx  what the study may claim

| ID | Severity | Rule | Reason |
|---|---|---|---|
| R801 | REFUSE | Claim scope is undeclared or species-level | These compounds are documented across the genus, in *Z. jujuba* and *Z. mauritiana*. They were not isolated from *Z. talanai* material |
| R803 | WARN | The claim reads as demonstrated, not predicted | Docking predicts binding. "Inhibits", "proves", "outperforms" turn a prediction into a finding with no new evidence |
| R804 | REFUSE | A `ziziphus/` library declares a species scope | Every experiment importing it would inherit a species-level claim from genus-level evidence, quietly and at scale. The boundary is in the namespace name, and this rule keeps it there |

## R9xx  the receptor's provenance

| ID | Severity | Rule | Reason |
|---|---|---|---|
| R901 | WARN | Surrogate receptor undisclosed | 3A4A is *S. cerevisiae* isomaltase standing in for the human intestinal enzyme. Standard and defensible; burying it is what causes trouble |
| R902 | WARN | Resolution missing, or worse than 2.5 Å | Resolution bounds what the contact analysis can claim about side-chain placement |

## RAxx  replication and reporting

| ID | Severity | Rule | Reason |
|---|---|---|---|
| RA01 | RECORD | Replication kind must be stated | Same-seed replication proves the pipeline is faithful. Different-seed replication tests whether the search converged. Both are worth having and neither substitutes for the other |
| RA02 | WARN | Replication deviates beyond run-to-run reproducibility | This checks run-to-run numerical reproducibility, a project convention (0.1 kcal/mol), not a published Vina scoring-precision figure; Vina's own accuracy against experiment is reported at 2 to 3 kcal/mol. Agreement inside the convention is the point. Disagreement outside it means something changed |
