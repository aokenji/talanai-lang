# Talanai

A language for describing molecular docking experiments, and a validator that
refuses the ones that could not mean what they claim to mean.

A `.tal` file is one experiment written down completely: the receptor and how
it was prepared, the pocket, the positive control, the ligands, the search
settings, and the results. The validator is the product. Docking is the easy
part.

**v0.1.0**, 2026-07-31. Python standard library only, offline, no installs.

## Why

A docking study today usually cannot be repeated from its own methods section.
The box coordinates, the seed, the exhaustiveness, the ligand protonation and
the receptor preparation are routinely missing. A `.tal` file is one page that
another lab can execute, and that will not run when the protocol has not passed
its own positive control.

The claim is reproducibility, not accessibility. Easier docking already exists:
PyRx, CB-Dock2, DockingPie, SwissDock, Galaxy. A portable, checkable,
self-refusing experiment record does not.

## Try it

```
tal check   examples/alpha-glucosidase.tal
tal explain examples/alpha-glucosidase.tal
tal run     examples/alpha-glucosidase.tal --dry-run
tal report  examples/alpha-glucosidase.tal

python tests/test_rules.py
python tests/acceptance.py
```

On the thesis protocol as it stood in May 2026, `check` **refused**, on R106.
After the three corrections of August 2026 it **accepts**, with one warning and
two recorded items. Both outcomes are the tool working: the refusal is what it
is for, and the acceptance is what closing the finding looks like. The findings
and their closures are listed below rather than quietly dropped.

## The five commands

| Command | Does | Needs Vina |
|---|---|---|
| `tal check` | Validate. Returns instantly | no |
| `tal explain` | Validate, plus both rankings and the generated methods text | no |
| `tal control` | Run the redocking control that unlocks docking | yes |
| `tal run` | Dock, only if the file validates. `--dry-run` prints the plan | yes |
| `tal report` | Write `results.csv`, `methods.txt` and the run record | no |

`check` never touches the docking engine. Validation has to stay usable on a
laptop with nothing installed, which is most of the people this is for.

**Timing, stated accurately.** `check` returns in about 0.3 s on its own. When
the experiment points at prepared ligand files and RDKit is installed, R307
reads their ring geometry and the command takes about 3 s, almost all of it
RDKit's import. Without RDKit, R307 reports UNVERIFIED and the command stays
sub-second. It never requires Vina.

## What it refuses

Thirty rules across three severities, each with the real incident behind it.
The full table is in [RULES.md](RULES.md); the format is in
[FORMAT.md](FORMAT.md). The ones that matter most:

- **R101 / R105 / R106** the control. Docking stays locked until a redocking
  control is recorded, has a citable source, and was run on the same receptor
  preparation as the screen.
- **R601** anything compared must share a preparation. A Vina score belongs to
  a pair of prepared files, not to a molecule.
- **R201** `rmsd` alone is banned, because this project uses the word for two
  different quantities.
- **R801 / R804** a claim may not exceed its evidence, and the `ziziphus/`
  namespace is genus-level by construction.
- **R303** the box must actually contain the catalytic residues, verified
  against the PDB when it is on disk and reported UNVERIFIED when it is not.
- **R306** the box is not the pose. R303 only checks that the search box
  encloses the catalytic residues; R306 checks that the rank-1 pose itself
  comes within 4.0 Å of one, WARN when it does not and UNVERIFIED when there
  is no pose file to check. A screen can pass R303 while every reported pose
  sits in a non-catalytic pocket, which is the exact failure this rule set
  otherwise only names in prose.

## Layout

```
tal.cmd, tal.py        launcher
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
tests/                 unit tests, and the acceptance test against compounds.js
```

## Acceptance test

`tests/acceptance.py` reads the affinities out of
`D:\BALAKATDBV2\src\data\compounds.js` at runtime and checks the example file
agrees to the digit, formulae included. Nothing is typed into the test, so it
breaks the moment the two records drift apart. It also asserts that the
ranking inverts under ligand efficiency, because that inversion is the reason
the fourth column exists.

Talanai only ever reads from `D:\BALAKATDBV2`. It writes nothing there.

## What v0 does not do

Deliberately absent, so nobody has to maintain them: a compiler, a lexer, a
parser generator, an AST, a package manager, a language server. The runner
executes local Vina offline and never contacts the Hugging Face Space.

`tal control` does not yet run the control itself; it locates Vina and tells
you what to record. That is the next thing to build.

## What the validator found, and what closed it

Found by running the validator against the real files, not by reading them.
Every blocking finding below is now closed. They are kept on the page because a
tool that only ever shows its clean state is not evidence of anything.

1. **R106, was blocking, closed 2026-08-04.** The redocking control in
   `docking_data/validation/validation.json` had been run on the
   **Meeko-prepared** receptor, in an 18 Å box, at exhaustiveness 16. The screen
   ran on the **raw cleaned** receptor (`receptor_clean.pdb`, which `SOURCES.md`
   calls the file "used as the thesis ran it"), in a 30 Å box, at exhaustiveness
   8. The control therefore validated a protocol that was not the screening
   protocol, in exactly the dimension this project had already learned changes
   the number. Closed by moving the screen onto the Meeko-prepared receptor,
   one of the three variables in the 2026-08-04 correction.
2. **R401, R501, closed 2026-08-04.** Exhaustiveness 8 on a 30 Å cube, single
   seed, with a same-seed replication that demonstrated pipeline fidelity rather
   than search convergence. Closed by re-running at exhaustiveness 32 across
   seeds 42/43/44 and carrying the per-compound spread as data.
3. **R702, standing, and not a defect.** Ligand efficiency inverts the ranking.
   Rutin wins on raw score and finishes last of the five per heavy atom;
   Kaempferol leads. On raw score three of the five beat acarbose and two do
   not; per heavy atom all five do. The panel will ask.

Still open on the current file, and correctly so:

- **WARN R104.** The control ran in a tighter box than the screen. Focused-box
  redocking is standard and asks an easier question, so it is flagged rather
  than refused: no reader should assume the RMSD validates the screening search.
- **RECORD R306.** Whether the pose contacts the catalytic site is UNVERIFIED.
  The redocking gate is self-docking biased, and an independently built
  conformer of the same molecule misses by 5.9 Å. This protocol ranks
  compounds; it does not place them, and no pose-level claim rests on it.
- **RECORD RA01.** The replication used different seeds, which is the only way
  replicating a deterministic pipeline measures anything at all.

## Reconcile

`BUILD-SPEC.md` from the read-only side session is still not on disk. Where it
differs from what is here, it wins: its rule IDs, grammar, vocabulary and
command names were designed first. The naming, layout and rule numbering in
this repository were chosen independently and should be treated as provisional
until the two are merged.
