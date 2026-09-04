# The `.tal` format

A `.tal` file is one docking experiment, written down completely. It is meant
to be read by a biologist and executed by a computer, in that order of
priority.

## Three syntax rules, and that is all

```
# a hash starts a comment, to end of line

receptor 3A4A                 a line at column 0 opens a block, with an argument
  resolution   1.6 A          an indented line is a key and its value
  compound     Rutin          a key may repeat; repeats keep their order
  prepare_note SOURCES.md says this is the file the
               thesis actually ran                    a line indented DEEPER
                                                      than its key continues it
```

There is no punctuation, no nesting, no quoting rules, no types to remember.
Units are written as you would write them (`1.6 A`, `-8.857`, `pH 7.4`) and
read with the number.

## The eleven words

An experiment is built from eleven block words and nothing else.

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

Library files in `ziziphus/` use one further word, `library`, and are a
different kind of file: a shared, genus-level compound set. They are checked
against different rules, because they are not experiments.

## The smallest file that runs

```
study "screen"
  claim_scope genus Ziziphus

receptor 3A4A
  prepare raw cleaned, waters removed

site
  center 21.52 -7.70 23.55
  size 30 30 30
  must_enclose Asp215 Glu277 Asp352

control
  redock GLC 601
  measure redock_rmsd
  require under 2.0 A
  result 0.519 A
  source validation.json
  prepare raw cleaned, waters removed

ligands
  prepare rdkit + meeko
  protonate pH 7.4
  compound Quercetin C15H10O7

dock
  engine vina 1.2.7
  exhaustiveness 32
  seeds 42 43 44
```

Note what is compulsory and why. `claim_scope`, because a screen has to say
what it is a claim about. `prepare` on both the receptor and the ligands,
because preparation is part of the measurement. `source` on the control,
because a number with no artefact behind it is a claim rather than a result.
And `prepare` on the control, because a control validates the protocol it was
actually run under.

## Keys by block

**study** author, date, doi, claim_scope, claim, note
**receptor** name, organism, surrogate_for, resolution, file, prepare, prepare_note, note
**site** center, size, anchored_on, anchor_note, must_enclose, blind, note
**control** redock, measure, require, result, source, prepare, box_size, exhaustiveness, receptor_file, result_note, note
**ligands** source, dir, prepare, protonate, compound, note
**reference** formula, prepare, protonate, role, note
**dock** engine, exhaustiveness, modes, energy_range, seeds, cpu, vina_path, note
**rank** by
**results** affinity, pose_cluster_rmsd, pose_file
**replication** date, pipeline, kind, rerun
**library** scope, source, compound, reported_in, note

An unknown key is a warning, never a silent ignore. A misspelt setting that
quietly reverted to a default would be the worst possible failure mode, so the
reader refuses to be helpful about it.

## Compound lines

```
compound  Rutin  C27H30O16  rutin.pdbqt
```

Name, then optional formula, then optional prepared structure file. The
formula is what makes the ligand-efficiency column possible; without it that
column reads `-` rather than guessing.

## Checksums

`checksum` records what a preparation recipe actually produced. One line per
file, repeatable, allowed in `receptor`, `ligands`, `reference` and `control`:

    checksum      quercetin.pdbqt 3f1a9c...

The filename is matched on its basename, so a file named in one block and
checksummed in another is still covered, and the same file does not need
pinning twice.

Generate them with `tal checksum <file>`, which digests every prepared file the
experiment names and prints the lines to paste back.

R604 asks for these because a recipe does not determine its own output. Under
the recipe this study records, a re-prepared quercetin docked half a kilocalorie
weaker than the published one, on the same receptor, box, seed and
exhaustiveness. The engine reproduces perfectly; the conformer does not. Where
ring-aware selection applies it constrains the choice and the number comes back
(rutin, to 0.013). Where it does not, the fallback is lowest gas-phase energy,
which has nothing to do with docking.

If a checksum is recorded and the file on disk disagrees, R604 REFUSES. Do not
update the checksum to match whatever is there now: that is the check deleting
itself.

## Pose file

```
results
  pose_file  redock_screening_seed42.pdbqt
```

An optional key on `results`, pointing at the rank-1 docked pose (PDB or
PDBQT), resolved the same way `receptor.file` is: relative to the `.tal`
file itself. It exists for R306: R303 already verifies that the search BOX
encloses the residues named in `site.must_enclose`, but nothing checked that
the resulting POSE actually came near them until this key gave R306
something to read. Without it, R306 reports UNVERIFIED rather than a pass it
did not earn, the same convention R303 uses when no receptor file is on
disk.
