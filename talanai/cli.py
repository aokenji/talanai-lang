"""
The six commands.

    tal check    <file>            validate. No engine needed. Instant.
    tal explain  <file>            validate, plus rankings and the methods text
    tal checksum <file>            digest the prepared files, for R604
    tal control  <file>            run the redocking control that unlocks docking
    tal run      <file>            dock, only if the file validates
    tal report   <file>            write csv, methods and the run record

Exit codes: 0 accepted, 1 refused, 2 usage error. `check` touching no engine
is a design requirement, not an accident: validation has to stay usable on a
laptop with nothing installed.
"""

from __future__ import annotations

import os
import sys

from . import chem, engine, model, parse, pdb, report, rules
from . import control as control_runner

USAGE = __doc__


def load(path):
    if not os.path.isfile(path):
        print("  no such file: %s" % path)
        return None
    document = parse.parse_file(path, known_keys=model.KNOWN_KEYS)
    return model.Experiment(document)


def _validate(experiment, show_passes=True):
    findings = rules.run_all(experiment)
    report.print_findings(findings, experiment.path, experiment.title,
                          show_passes=show_passes)
    return findings


# --------------------------------------------------------------------------
def command_check(path, args):
    experiment = load(path)
    if experiment is None:
        return 2
    findings = _validate(experiment, show_passes="--quiet" not in args)
    return report.summarise(findings, rules.docking_is_locked(findings),
                            experiment.kind)


def command_explain(path, args):
    experiment = load(path)
    if experiment is None:
        return 2
    findings = _validate(experiment)
    report.print_ranking(experiment)
    print("  " + "-" * report.WIDTH)
    print("  Methods, generated from the parameters above")
    print("")
    for line in report.methods_paragraph(experiment).splitlines():
        print("    %s" % line)
    print("")
    return report.summarise(findings, rules.docking_is_locked(findings),
                            experiment.kind)


def _resolve(base, hint):
    if not hint:
        return None
    return hint if os.path.isabs(hint) else os.path.join(base, hint)


def command_control(path, args):
    """
    Run the redocking control. This is the command that unlocks `tal run`.

    It deliberately does not accept a number typed by hand. The point of a
    control is that the file records a run which happened, with an artefact
    behind it that someone else can find.

    By default it runs the control under the SCREENING protocol: the same
    receptor preparation, box and exhaustiveness the experiment uses. That is
    what makes the result mean something about this screen rather than about a
    neighbouring one. Pass --as-recorded to reproduce the control block's own
    settings instead, which is useful for comparing the two.
    """
    experiment = load(path)
    if experiment is None:
        return 2
    control = experiment.control
    if control is None:
        print("  no control block in %s" % path)
        return 1

    base = os.path.dirname(os.path.abspath(path))
    # Default mode runs against the SCREENING receptor, which is the whole
    # point: the control has to validate the protocol that produced the
    # affinities. control.receptor_file records what the recorded control used
    # and is only picked up under --as-recorded.
    receptor = _resolve(base, experiment.receptor.one("file")
                        if experiment.receptor else None)
    ligand = _resolve(base, control.one("ligand_file"))
    reference = _resolve(base, control.one("reference_file")) or ligand

    as_recorded = "--as-recorded" in args
    center, size = experiment.box()
    exhaustiveness = experiment.exhaustiveness() or 8
    if as_recorded:
        control_size = chem.numbers(control.one("box_size", ""))
        if len(control_size) == 3:
            size = control_size
        control_exhaustiveness = chem.as_float(control.one("exhaustiveness"))
        if control_exhaustiveness:
            exhaustiveness = control_exhaustiveness
        if control.one("receptor_file"):
            receptor = _resolve(base, control.one("receptor_file"))

    seed = chem.as_float(control.one("seed"))
    if seed is None:
        seeds = experiment.seeds()
        seed = seeds[0] if seeds else 42

    missing = [label for label, value in (
        ("control.ligand_file (the co-crystallised ligand, prepared)", ligand),
        ("a reference structure (control.reference_file)", reference),
        ("a receptor (receptor.file or control.receptor_file)", receptor))
        if not value or not os.path.isfile(value)]
    if missing or not (center and size):
        print("")
        print("  Cannot run the control. Missing:")
        for item in missing:
            print("    %s" % item)
        if not (center and size):
            print("    a fully specified box")
        print("")
        return 1

    try:
        vina = engine.find_vina(experiment)
    except engine.EngineMissing as error:
        print("")
        print("  %s" % error)
        print("")
        return 1

    outdir = os.path.join(report.default_output_dir(experiment), "control")
    mode = "as recorded" if as_recorded else "under the SCREENING protocol"
    print("")
    print("  Redocking control, %s" % mode)
    print("    vina           %s" % vina)
    print("    receptor       %s" % os.path.basename(receptor))
    print("    ligand         %s" % os.path.basename(ligand))
    print("    reference      %s" % os.path.basename(reference))
    print("    box            %.0f x %.0f x %.0f A at (%.2f, %.2f, %.2f)"
          % (size[0], size[1], size[2], center[0], center[1], center[2]))
    print("    exhaustiveness %d, seed %d" % (exhaustiveness, seed))
    print("")
    print("  running ...")

    record = control_runner.run(
        vina, receptor, ligand, reference, center, size, exhaustiveness,
        int(seed), outdir,
        modes=int(chem.as_float(experiment.dock.one("modes")) or 9),
        energy_range=int(chem.as_float(experiment.dock.one("energy_range")) or 3),
        label="screening" if not as_recorded else "as-recorded")

    name = "redock_%s_seed%d.json" % (
        "screening" if not as_recorded else "as-recorded", int(seed))
    artefact = control_runner.store(record, outdir, name)

    if record.get("error"):
        print("")
        print("  FAILED: %s" % record["error"].strip().splitlines()[-1][:200])
        print("  record: %s" % artefact)
        return 1

    value = record.get("rank1_rmsd_A")
    exact = record.get("rmsd_is_exact")
    limit = experiment.control_limit() or model.ACCEPTED_REDOCK_LIMIT
    print("")
    print("    rank 1 score   %.3f kcal/mol" % record["rank1_score_kcal_mol"])
    if value is None:
        print("    rank 1 RMSD    not computable: %s" % record["rmsd_method"])
    else:
        print("    rank 1 RMSD    %.3f A%s" % (value, "" if exact else "  (lower bound)"))
        print("    method         %s" % record["rmsd_method"])
        print("    best pose      %.3f A" % record["best_pose_rmsd_A"])
    print("    took           %.0f s" % record["seconds"])
    print("")
    print("  record written to %s" % artefact)
    print("")

    if value is None:
        print("  NO VERDICT. The RMSD could not be computed, so this run")
        print("  neither passes nor fails. Fix the atom correspondence first.")
        print("")
        return 1

    if not exact:
        # A lower bound is conclusive in one direction only.
        if value >= limit:
            print("  DEFINITIVELY DID NOT PASS under %.1f A." % limit)
            print("  This is a lower bound, so no atom assignment could do")
            print("  better. The pose is not in the crystal position.")
            print("")
            return 1
        print("  INCONCLUSIVE. The lower bound is under %.1f A, but a lower" % limit)
        print("  bound below a threshold is not a pass, only an absence of")
        print("  proof of failure. A symmetry-corrected RMSD is needed to")
        print("  settle it, and this tool cannot compute one without RDKit.")
        print("")
        return 1

    if value < limit:
        print("  PASSED under %.1f A. Record it in %s:" % (limit, os.path.basename(path)))
        print("")
        print("      control")
        print("        result   %.3f A" % value)
        print("        source   %s" % os.path.relpath(artefact, base))
        print("        prepare  %s" % (experiment.receptor.one("prepare") or ""))
        print("")
        return 0

    print("  DID NOT PASS under %.1f A." % limit)
    print("  That is a result, not a failure of the tool. Under this protocol")
    print("  the search does not recover the crystal pose.")
    print("")
    return 1


def command_run(path, args):
    experiment = load(path)
    if experiment is None:
        return 2
    findings = rules.run_all(experiment)
    locked = rules.docking_is_locked(findings)
    if locked:
        report.print_findings([f for f in findings if f.level == rules.REFUSE],
                              experiment.path, experiment.title)
        print("  " + "-" * report.WIDTH)
        print("  REFUSED. Docking is locked until these are resolved.")
        print("")
        return 1

    jobs = engine.plan(experiment)
    dry = "--dry-run" in args
    print("")
    print("  %d job(s): %d compound(s) x %d seed(s)"
          % (len(jobs), len(experiment.compounds()), len(experiment.seeds()) or 1))

    missing = engine.missing_inputs(jobs)
    if missing:
        print("")
        print("  Missing input files:")
        for item in missing:
            print("    %s" % item)
        print("")
        print("  Declare 'file' on the receptor block and 'dir' on the ligands")
        print("  block, pointing at prepared structures.")
        if not dry:
            return 1

    try:
        vina = engine.find_vina(experiment)
    except engine.EngineMissing as error:
        print("")
        print("  %s" % error)
        if not dry:
            return 1
        vina = "vina"

    workdir = report.default_output_dir(experiment)
    if dry:
        print("")
        print("  Dry run. Commands that would execute:")
        print("")
        for job in jobs[:6]:
            print("    %s --config %s --out %s"
                  % (os.path.basename(vina),
                     os.path.join(workdir, "%s_seed%s_%s.conf"
                                  % (job.ligand, job.seed, job.fingerprint())),
                     os.path.join(workdir, "%s_seed%s_out.pdbqt"
                                  % (job.ligand, job.seed))))
        if len(jobs) > 6:
            print("    ... and %d more" % (len(jobs) - 6))
        print("")
        print("  Cache keys are content-addressed, so a re-run after an")
        print("  interruption resumes rather than restarting.")
        print("")
        return 0

    done, reused = 0, 0
    for job in jobs:
        record = engine.cached(experiment, job)
        if record:
            reused += 1
        else:
            record = engine.execute(experiment, job, vina, workdir)
            engine.store(experiment, job, record)
        affinity = record.get("affinity")
        print("    %-24s %s" % (job.label(),
                                "%.3f kcal/mol" % affinity if affinity
                                else "FAILED"))
        done += 1
    print("")
    print("  %d run, %d reused from cache" % (done - reused, reused))
    print("")
    return 0


def command_report(path, args):
    experiment = load(path)
    if experiment is None:
        return 2
    findings = rules.run_all(experiment)
    outdir = report.default_output_dir(experiment)
    os.makedirs(outdir, exist_ok=True)
    written = []
    if experiment.affinities():
        written.append(report.write_csv(experiment, os.path.join(outdir, "results.csv")))
    written.append(report.write_methods(experiment, os.path.join(outdir, "methods.txt")))
    print("")
    for item in written:
        print("  wrote %s" % item)
    if rules.docking_is_locked(findings):
        print("")
        print("  Note: this experiment is REFUSED. The methods text says so.")
    print("")
    return 0


def command_checksum(path, args):
    """
    Print a checksum line for every prepared file this experiment names, ready
    to paste back into the .tal.

    R604 asks for these because a preparation recipe does not determine its own
    output. Re-preparing quercetin from its SMILES under the recorded recipe
    produced a conformer worth half a kilocalorie less than the published one,
    on the same receptor, box and seed. A digest is what turns "prepared this
    way" into "prepared into exactly this file".

        tal checksum examples/alpha-glucosidase.tal
    """
    experiment = load(path)
    if experiment is None:
        return 2

    named = rules.structure_files(experiment)
    if not named:
        print("")
        print("  %s names no prepared structure files." % os.path.basename(path))
        print("")
        return 0

    base = os.path.dirname(os.path.abspath(path))
    print("")
    print("  CHECKSUMS for %s" % os.path.basename(path))
    print("  " + "=" * 74)

    found, missing = [], []
    for what, raw in named:
        resolved = pdb.resolve(raw, base)
        if resolved:
            found.append((what, raw, rules._digest(resolved)))
        else:
            missing.append((what, raw))

    if missing:
        print("")
        print("  NOT ON THIS MACHINE, so no digest can be taken:")
        for what, raw in missing:
            print("    %-22s %s" % (what, raw))
        print("")
        print("  A checksum can only be recorded where the file is. Run this")
        print("  on the machine that holds the prepared files.")

    if not found:
        print("")
        return 1

    by_block = {}
    for what, raw, digest in found:
        block = ("receptor" if what == "receptor"
                 else "control" if what.startswith("control ")
                 else "ligands")
        by_block.setdefault(block, []).append((os.path.basename(raw), digest))

    print("")
    print("  Paste these into the matching blocks:")
    for block in ("receptor", "ligands", "control"):
        if block not in by_block:
            continue
        print("")
        print("  %s" % block)
        for name, digest in by_block[block]:
            print("    checksum      %s %s" % (name, digest))

    print("")
    print("  %d file(s) digested with %s.%s"
          % (len(found), rules.CHECKSUM_ALGORITHM,
             " %d could not be reached." % len(missing) if missing else ""))
    print("")
    return 0


COMMANDS = {
    "check": command_check,
    "checksum": command_checksum,
    "explain": command_explain,
    "control": command_control,
    "run": command_run,
    "report": command_report,
}


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] in ("-h", "--help", "help"):
        print(USAGE)
        return 2
    command = argv[0]
    if command not in COMMANDS:
        print("  unknown command: %s" % command)
        print(USAGE)
        return 2
    if len(argv) < 2:
        print("  %s needs a .tal file" % command)
        return 2
    return COMMANDS[command](argv[1], argv[2:])
