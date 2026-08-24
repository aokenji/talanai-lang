# Talanai: build specification

**Status:** v0.1.0 built and green, 2026-07-31.
**Target:** innovation convention, 3 to 5 September 2026. Feature freeze 29 August.
**Thesis defence:** around December 2026. No collision.

> This document was reconstructed on 2026-07-31 after the original, written in
> a read-only side session, never reached disk. Where the original resurfaces
> and disagrees, the original wins on rule IDs, vocabulary and command names,
> since those were designed first.

---

## 1. What this is

Talanai is a domain-specific language for describing molecular docking
experiments, plus a validator that refuses the ones that could not mean what
they claim to mean.

A `.tal` file is one experiment written down completely: the receptor and how
it was prepared, the pocket, the positive control, the ligands, the search
settings, and the results.

**The validator is the product.** Docking is the easy part. Anyone can dock.
Almost nothing in this field stops you from docking badly.

- Command: `tal`
- Files: `.tal`
- Shared library namespace: `ziziphus/`, genus level only, never species level
- Home: `D:\THESIS_VSC\talanai-lang`, touching nothing in `D:\BALAKATDBV2`

## 2. Why, and the claim to make

A docking study today usually cannot be repeated from its own methods section.
The box coordinates, the seed, the exhaustiveness, the ligand protonation and
the receptor preparation are routinely absent.

**The claim is reproducibility, not accessibility.**

This distinction decides whether the pitch survives a judge who knows the
field. "Docking made easy for non-programmers" has incumbents who will be
named at you within thirty seconds: PyRx, CB-Dock2, DockingPie, SwissDock,
Galaxy, KNIME. A portable, checkable, self-refusing experiment record has
none.

The sentence to say:

> A docking study today cannot be repeated from its own methods section. We
> built a format where the experiment is a single file another lab can
> execute, and which refuses to run when the protocol has not passed its own
> positive control.

Accessibility is the second sentence, never the first.

## 3. The design rule

**Hide plumbing. Never hide science.**

| Hidden, because nobody's result depends on it | Visible and short, because it *is* the result |
|---|---|
| PDB to PDBQT conversion, atom typing | Which pocket, and how the box was placed |
| Invoking obabel, Meeko, RDKit, Vina | Receptor preparation: waters, hetatms, H, charges |
| Threads, temp files, retries, caching | Ligand protonation state at pH |
| Parsing output, writing CSV and figures | Exhaustiveness, seeds, replicate count |
| Install locations, file paths | Scoring function, and what a ranking ranked on |
| | Whether the protocol passed a redocking control |

A single `dock` keyword that silently picks preparation, box and protonation
would be worse than the status quo: it would hand unexamined choices to people
who by construction cannot examine them, and print them to three decimals.

## 4. Non-goals

Deliberately not built, so that nobody has to maintain them:

- A compiler, lexer, parser generator, or AST
- A package manager, a formatter, a language server
- Turing completeness, user-defined functions, control flow beyond iteration
- Any network path. The runner uses **local Vina**, never the Hugging Face Space
- A GUI. TalanaiHub and TalanaiDock already are the GUI
- Automatic fixes. The validator explains; the scientist decides

## 5. The `.tal` grammar

Three syntax rules, and no others.

```
# a hash starts a comment, to end of line

receptor 3A4A                 a line at column 0 opens a block, with an argument
  resolution   1.6 A          an indented line is a key and its value
  compound     Rutin          a key may repeat; repeats keep their order
  prepare_note SOURCES.md says this is the file the
               thesis actually ran                    a line indented DEEPER
                                                      than its key continues it
```

No punctuation, no nesting, no quoting rules, no declared types. Units are
written as a person writes them (`1.6 A`, `-8.857`, `pH 7.4`) and read with
the number.

An unknown key is a **warning**, never a silent ignore. A misspelt setting
that quietly reverted to a default would be the worst possible failure mode,
so the reader refuses to be helpful about it.

### The eleven words

| Word | Holds |
|---|---|
| `study` | what this is, who ran it, and what it claims |
| `receptor` | the protein, and exactly how it was prepared |
| `site` | where in the protein you are looking |
| `control` | the redocking control that unlocks docking |
| `ligands` | what is being screened, and how it was prepared |
| `reference` | the comparator drug |
| `dock` | the search settings |
| `rank` | how results are ordered |
| `report` | what to emit |
| `results` | what was obtained |
| `replication` | a re-run, and what kind it was |

A twelfth word, `library`, appears only in `ziziphus/` files. Those are a
different kind of file, judged against different rules, because they are not
experiments.

Full key list: `FORMAT.md`.

## 6. The rules

Three severities:

- **REFUSE** the experiment does not run. Not a preference: the result would
  not be interpretable.
- **WARN** it runs, and you should know this before presenting it.
- **RECORD** neither wrong nor optional. A fact that must appear in the run
  record and the generated methods paragraph.

A check that cannot be performed reports RECORD with the word **UNVERIFIED**.
It never reports a pass it did not earn.

| ID | Sev | Rule | Reason, in biology |
|---|---|---|---|
| R000 | REFUSE | File unreadable | A setting outside any block belongs to nothing |
| R001 | REFUSE | Required block missing | Without study, receptor, site, control, ligands, dock the experiment is sketched, not described |
| R002 | WARN | Unknown block | A misspelt block would otherwise be silently ignored |
| R003 | WARN | Unknown key | A typo must never read as a default |
| R004 | REFUSE | Single-value setting given twice | Two exhaustiveness values means nobody knows which ran |
| R101 | REFUSE | Control has no recorded result | A protocol that cannot reproduce a known answer has no business predicting an unknown one. This is the positive-control lane on a blot |
| R102 | REFUSE | Control failed its threshold | A failed control is information. The pocket, the prep or the search is wrong |
| R103 | WARN | Threshold looser than 2.0 Å | 2.0 Å is the accepted bar; looser has to be argued in text |
| R104 | WARN | Control used a different box or higher exhaustiveness than the screen | Focused-box redocking asks an easier question than the screen. Legitimate, but a reader will assume the RMSD validates the screening search |
| R105 | REFUSE | Control result has no source | A live endpoint reading is not a citable artefact |
| R106 | REFUSE | Control validated a different receptor preparation | **The one that bites.** Change the preparation and the thing that passed at 2.0 Å is not the thing that produced the affinities |
| R201 | REFUSE | Bare `rmsd` | This project uses RMSD for redocking accuracy *and* for pose-cluster spread (1.61 to 2.00 Å). Quoting one for the other is a scientific error, not a typo |
| R301 | REFUSE | Box not fully specified | A box is three coordinates and three lengths |
| R302 | WARN | Box declares no residues to enclose | An unanchored box is the commonest way to produce confident nonsense: a perfect search in the wrong pocket |
| R303 | REFUSE / RECORD | Box must enclose the named catalytic residues | Verified against the PDB when on disk, UNVERIFIED when not |
| R304 | WARN | Named residue absent from the structure | Numbering or chain mismatch |
| R305 | WARN | Box large enough to be blind docking, undeclared | Blind docking is a different method with a different interpretation and a far larger search requirement |
| R306 | WARN / RECORD | Rank-1 pose does not contact the named catalytic residues | R303 checks the box; this checks the pose. UNVERIFIED with no pose file on disk, WARN past the 4.0 Å contact cutoff. A pose that never reaches the site does not support an active-site claim |
| R401 | WARN | Search budget low for box volume | Vina's default exhaustiveness of 8 was chosen for small boxes. Rule of thumb, not a standard |
| R402 | RECORD | `modes` or `energy_range` unrecorded | The pose ensemble cannot be described without them |
| R501 | WARN | Single seed | Vina searches; where it starts affects where it lands. One seed repeated is reproducibility, not convergence. A fixed seed is a scale that reads the same every time and can still be off |
| R601 | REFUSE | Ligands and reference prepared differently | June 2026, this project: Quercetin returned −7.525 only once Meeko prep was removed. A Vina score belongs to a pair of prepared files, so two preparations are two assays |
| R602 | REFUSE | Receptor preparation undeclared | The most commonly omitted item in published docking methods |
| R603 | WARN | Ligand protonation undeclared | The structure drawn in a paper is not what exists at pH 7.4 |
| R701 | WARN | Ranking on raw affinity only | Vina scores rise with heavy-atom count, so this partly ranks by molecular weight |
| R702 | RECORD | Strongest binder is also the heaviest compound | The size confound is live. Recording it is the answer prepared in advance |
| R801 | REFUSE | Claim scope undeclared or species-level | These compounds are documented across the genus, in *Z. jujuba* and *Z. mauritiana*, not isolated from *Z. talanai* |
| R803 | WARN | Claim reads as demonstrated, not predicted | "Inhibits", "proves", "outperforms" turn a prediction into a finding with no new evidence |
| R804 | REFUSE | A `ziziphus/` library declares a species scope | Every experiment importing it would inherit a species-level claim, quietly and at scale |
| R901 | WARN | Surrogate receptor undisclosed | 3A4A is *S. cerevisiae* isomaltase standing in for the human enzyme. Standard; burying it is what causes trouble |
| R902 | WARN | Resolution missing or worse than 2.5 Å | It bounds what contact analysis may claim about side chains |
| RA01 | RECORD | Replication kind must be stated | Same-seed proves the pipeline is faithful. Different-seed tests convergence. Neither substitutes for the other |
| RA02 | WARN | Replication deviates beyond run-to-run reproducibility | A project convention (0.1 kcal/mol) for pipeline reproducibility, not a published Vina scoring-precision figure; Vina's own accuracy against experiment is reported at 2 to 3 kcal/mol. Agreement inside the convention is the point; outside it, something changed |

Authoritative version with full prose: `RULES.md`.

## 7. Architecture

```
tal.cmd, tal.py        launcher; TALANAI_PYTHON selects the interpreter
talanai/
  parse.py             the .tal reader. Not a parser generator
  model.py             the eleven words, and typed access
  rules.py             THE PRODUCT
  chem.py              formulae, heavy atoms, ligand efficiency
  pdb.py               enough PDB reading to check the box encloses the site
  engine.py            offline local Vina, content-addressed cache
  report.py            console output, CSV, generated methods paragraph
  cli.py               the five commands
examples/              the defended thesis screen as a .tal
ziziphus/              shared, genus-level compound sets
tests/                 unit tests, and acceptance against compounds.js
```

### The five commands

| Command | Does | Needs Vina |
|---|---|---|
| `tal check` | Validate. Returns instantly | no |
| `tal explain` | Validate, plus both rankings and the generated methods text | no |
| `tal control` | Run the redocking control that unlocks docking | yes |
| `tal run` | Dock, only if the file validates. `--dry-run` prints the plan | yes |
| `tal report` | Write `results.csv`, `methods.txt`, run record | no |

Exit codes: `0` accepted, `1` refused, `2` usage error.

### The cache key

Content-addressed on everything that can change the number: receptor file
identity and preparation, box, ligand file identity and preparation, engine
version, exhaustiveness, modes, energy range, seed. Change any of them and the
old result is not reused, because it is not the same measurement. This is what
makes an interrupted screen resumable rather than restartable.

## 8. Constraints

1. **Python standard library only.** No pip, no bindings. Vina is invoked as
   the executable it already is.
2. **Must run on the interpreter bundled with TalanaiDock.** `TALANAI_PYTHON`
   selects it; the launcher falls back to `python`.
3. **Fully offline.** Local Vina at
   `...\TalanaiDock\app\docking_assets\vina.exe`. No Hugging Face Space, no
   network path anywhere in the runner.
4. **Light enough for an i3 with 8 GB.** One ligand at a time by default.
5. **`tal check` needs no docking engine and must return instantly.**
   Validation has to be usable on a laptop with nothing installed, which is
   most of the people this is for.
6. **Read-only toward `D:\BALAKATDBV2`.** Talanai reads `compounds.js`,
   `protocol.js` and the TalanaiDock assets. It writes nothing there.

## 9. Acceptance criteria

**A. The thesis numbers, read not typed.** `tests/acceptance.py` reads
affinities and formulae out of `D:\BALAKATDBV2\src\data\compounds.js` at
runtime and requires the example file to agree to the digit. The numbers are
never written into the test, so it fails the moment the two records drift.

| Compound | kcal/mol | Heavy atoms | Per atom |
|---|---|---|---|
| Rutin | −8.857 | 43 | 0.206 |
| Betulinic acid | −8.290 | 33 | 0.251 |
| Quercetin | −7.503 | 22 | 0.341 |
| Kaempferol | −7.479 | 21 | 0.356 |
| Oleanolic acid | −6.922 | 33 | 0.210 |
| Acarbose (reference) | −6.660 | 44 | 0.151 |

**B. The ranking inverts.** Rutin wins on raw score and finishes last of the
five per heavy atom; Kaempferol leads. Everything still beats acarbose on both
measures, so the conclusion holds. The test asserts the two rankings disagree,
because that disagreement is the reason the fourth column exists.

**C. Rules are individually tested.** Each rule has a minimal `.tal` that
triggers it, so a failure names the rule rather than "the example changed".

**D. `tal check` returns in well under a second** with no engine installed.

## 10. The stage demo, three minutes

Offline, cached, with a recorded parachute. Nothing docks live.

| Time | Beat |
|---|---|
| 0:00 | A real published methods paragraph on screen. Ask the room to rerun it. Four numbers are missing: box, seed, exhaustiveness, preparation |
| 0:30 | The `.tal` file. Read it aloud. It is a methods section that happens to execute |
| 1:15 | `tal check` on a good file: ACCEPTED |
| 1:35 | Delete the control line. `tal check` again: **REFUSED**, docking locked. "You cannot screen with an unvalidated protocol" |
| 2:00 | `tal explain`: the ranking inverts under ligand efficiency. The tool argues with its own author |
| 2:40 | `tal report`: the methods paragraph, generated from the parameters that ran. No transcription step, so no drift |
| 3:00 | Stop |

**Parachute:** a screen recording of the full sequence, on the laptop, playable
without a terminal. Convention wifi does not exist as far as this demo is
concerned. Cached results ship in the repository so nothing has to compute.

## 11. Schedule

| Week | Work |
|---|---|
| Aug 3 | Close R106: redock GLC 601 on the **screening** receptor and record it. Wire `tal control` to run it |
| Aug 10 | Three-seed run on the five compounds, to answer convergence. `tal init`. A second worked example on a different target, to prove nothing is hardcoded to 3A4A |
| Aug 17 | The opening-slide paper. Slides and poster. First full rehearsal |
| Aug 24 | Rehearse to time. Record the parachute video |
| **Aug 29** | **Feature freeze.** Bug fixes only after this |
| Sept 3 to 5 | Convention |

**Ken's own task, not delegable:** find a published docking paper whose methods
you genuinely cannot rerun. That paper is the opening slide, and no tool can
produce it for you.

## 12. Open scientific items

1. **R106, blocking, and now measured.** `tal control` was run on 2026-07-31
   and the picture is sharper than a disclosure issue.

   The recorded control (`docking_data/validation/validation.json`, 2026-06-16,
   RMSD 0.519 Å) used the **Meeko-prepared** receptor, an 18 Å box and
   exhaustiveness 16. The screen used the **raw cleaned** receptor, a 30 Å box
   and exhaustiveness 8. Those three differences were decomposed:

   | Case | Receptor | Box | Exh | Score | RMSD |
   |---|---|---|---|---|---|
   | B calibration | Meeko | 18 Å | 16 | −5.885 | **0.513** |
   | D search only | Meeko | 30 Å | 8 | −5.905 | **0.518** |
   | C prep only | raw | 18 Å | 16 | −4.208 | **5.592** |
   | A screening | raw | 30 Å | 8 | −4.024 | **6.164** |

   Case B reproduces the published 0.519 Å and −5.897, which calibrates the
   RMSD implementation. Case D shows the 30 Å screening box and exhaustiveness
   8 recover the crystal pose perfectly well. **The box was never the problem.
   The receptor preparation is.** The raw path fails its own positive control
   by a wide margin at every search setting tried.

   RMSD here is a lower bound, not symmetry-corrected, because the docking
   input names every atom `C` or `O` while the crystal reference uses
   `C1..C6` / `O1..O6`. A lower bound above the threshold is conclusive for
   failure, so the 5.6 to 6.2 Å results stand.

2. **Consequence, exploratory.** Re-docking the thesis set on the receptor
   that passes its control, same box, same exhaustiveness, same seed:

   | Compound | Thesis (raw) | Prepared | Δ |
   |---|---|---|---|
   | Rutin | −8.857 | −10.620 | −1.763 |
   | Betulinic acid | −8.290 | −9.644 | −1.354 |
   | Quercetin | −7.503 | −8.889 | −1.386 |
   | Kaempferol | −7.479 | −8.240 | −0.761 |
   | Oleanolic acid | −6.922 | −7.792 | −0.870 |
   | **Acarbose** | −6.660 | **−8.576** | **−1.916** |

   Everything binds more strongly on the prepared receptor, but acarbose gains
   most and overtakes kaempferol and oleanolic acid. On the validated receptor
   the result is **3 of 5 beat acarbose**, not 5 of 5. Rutin remains the
   strongest binder by a clear margin in both.

   Single seed, exhaustiveness 8, so this is a signal to investigate, not a
   number to publish. A three-seed check on the affected compounds is the
   deciding test.
2. **R401, R501.** Exhaustiveness 8 on a 30 Å cube, single seed. The 2026-07-31
   replication was same-seed, demonstrating pipeline fidelity rather than
   search convergence.
3. **R702.** The size confound is live and should be answered in the thesis
   before it is asked at a panel.
