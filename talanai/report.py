"""
Output. Console findings, the ranking table, CSV, and the methods paragraph.

The methods paragraph matters more than it looks. It is generated from the
same file the run used, so the sentence in the paper and the parameters in the
experiment cannot drift apart. That is the difference between a methods
section and a description of one.

ASCII only, deliberately: this has to render in the Windows console on the
machine that runs the docking.
"""

from __future__ import annotations

import csv
import os

from . import chem, rules

WIDTH = 74


def _wrap(text, indent):
    """Wrap to WIDTH without importing textwrap's whole machinery."""
    out, line = [], ""
    for word in str(text).split():
        if len(line) + len(word) + 1 > WIDTH - len(indent):
            out.append(line)
            line = word
        else:
            line = (line + " " + word).strip()
    if line:
        out.append(line)
    return out


def print_findings(findings, path, title, show_passes=True):
    print("")
    print("  talanai  %s" % path)
    print("  %s" % title)
    print("  " + "-" * WIDTH)
    ordered = sorted(findings, key=lambda f: (rules.SEVERITY_ORDER[f.level], f.code))
    for finding in ordered:
        if finding.level == rules.PASS and not show_passes:
            continue
        location = " (line %d)" % finding.line if finding.line else ""
        print("  %-7s %s  %s%s" % (finding.level, finding.code,
                                   finding.title, location))
        for raw in str(finding.detail).split("\n"):
            if not raw.strip():
                continue
            for line in _wrap(raw, " " * 12):
                print("            %s" % line)
        if finding.fix:
            wrapped = _wrap(finding.fix, " " * 15)
            print("            -> %s" % wrapped[0])
            for line in wrapped[1:]:
                print("               %s" % line)
        print("")


def print_ranking(experiment):
    affinities = experiment.affinities()
    if not affinities:
        return
    compounds = {c.name: c for c in experiment.compounds()}
    rows = []
    for name, score in affinities.items():
        compound = compounds.get(name)
        atoms = compound.heavy_atoms if compound else 0
        efficiency = chem.ligand_efficiency(score, compound.formula) if compound and compound.formula else None
        rows.append((name, score, atoms, efficiency))

    print("  " + "-" * WIDTH)
    print("  Ranking, two ways")
    print("")
    print("    %-18s %10s %7s %10s" % ("compound", "kcal/mol", "heavy", "per atom"))
    for name, score, atoms, efficiency in sorted(rows, key=lambda r: r[1]):
        print("    %-18s %10.3f %7s %10s"
              % (name.replace("_", " "), score, atoms or "-",
                 "%.3f" % efficiency if efficiency else "-"))
    print("")

    scored = [r for r in rows if r[3]]
    if len(scored) >= 2:
        best_raw = min(scored, key=lambda r: r[1])[0]
        best_efficiency = max(scored, key=lambda r: r[3])[0]
        if best_raw != best_efficiency:
            print("    Strongest raw score:     %s" % best_raw.replace("_", " "))
            print("    Most efficient per atom: %s" % best_efficiency.replace("_", " "))
            print("    The ranking changes with the measure. Report both.")
            print("")


def summarise(findings, locked, kind="experiment"):
    refusals = [f for f in findings if f.level == rules.REFUSE]
    warnings = [f for f in findings if f.level == rules.WARN]
    records = [f for f in findings if f.level == rules.RECORD]
    library = kind == "library"
    print("  " + "-" * WIDTH)
    if locked:
        print("  REFUSED. %d blocking, %d warning, %d recorded."
              % (len(refusals), len(warnings), len(records)))
        print("  This library cannot be imported." if library else
              "  Docking is locked. Nothing would run from this file.")
    else:
        print("  ACCEPTED. %d warning, %d recorded." % (len(warnings), len(records)))
        print("  This library is importable." if library else
              "  This experiment is runnable.")
    print("")
    return 1 if locked else 0


# --------------------------------------------------------------------------
# Artefacts
# --------------------------------------------------------------------------
def write_csv(experiment, path):
    affinities = experiment.affinities()
    compounds = {c.name: c for c in experiment.compounds()}
    reruns = experiment.reruns()
    rows = []
    for name, score in sorted(affinities.items(), key=lambda kv: kv[1]):
        compound = compounds.get(name)
        formula = compound.formula if compound else ""
        atoms = compound.heavy_atoms if compound else ""
        efficiency = (chem.ligand_efficiency(score, formula) if formula else None)
        rows.append({
            "compound": name.replace("_", " "),
            "formula": formula or "",
            "heavy_atoms": atoms,
            "affinity_kcal_per_mol": "%.3f" % score,
            "ligand_efficiency": "%.3f" % efficiency if efficiency else "",
            "rerun_kcal_per_mol": "%.3f" % reruns[name] if name in reruns else "",
        })
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    return path


def methods_paragraph(experiment):
    """
    The methods section, written from the parameters that actually ran.

    Only states what the file records. Where the file records nothing, the
    sentence says so rather than filling in a plausible default.
    """
    exp = experiment
    receptor = exp.receptor
    center, size = exp.box()
    sentences = []

    organism = receptor.one("organism") if receptor else None
    surrogate = receptor.one("surrogate_for") if receptor else None
    resolution = receptor.one("resolution") if receptor else None
    target = "%s (PDB %s%s)" % (
        receptor.one("name") or "the receptor",
        receptor.arg,
        ", %s" % resolution if resolution else "")
    if organism:
        target += ", %s" % organism
        if surrogate:
            target += " used as a structural surrogate for %s" % surrogate

    sentences.append(
        "Molecular docking was performed against %s." % target)

    if receptor and receptor.one("prepare"):
        sentences.append(
            "The receptor was prepared %s." % receptor.one("prepare"))
    if exp.ligands and exp.ligands.one("prepare"):
        protonation = exp.ligands.one("protonate")
        sentences.append(
            "Ligands were prepared with %s%s."
            % (exp.ligands.one("prepare"),
               " at %s" % protonation if protonation else ""))

    if center and size:
        anchor = exp.site.one("anchored_on") if exp.site else None
        sentences.append(
            "The search box measured %.0f x %.0f x %.0f A centred at "
            "(%.2f, %.2f, %.2f)%s."
            % (size[0], size[1], size[2], center[0], center[1], center[2],
               ", on %s" % anchor if anchor else ""))
    if exp.site and exp.site.one("must_enclose"):
        sentences.append(
            "The box encloses %s." % exp.site.one("must_enclose").replace(" ", ", "))

    if exp.dock:
        seeds = exp.seeds()
        detail = []
        if exp.dock.one("engine"):
            detail.append("AutoDock %s" % exp.dock.one("engine").replace("vina", "Vina"))
        if exp.exhaustiveness():
            detail.append("exhaustiveness %d" % exp.exhaustiveness())
        if exp.dock.one("modes"):
            detail.append("%s modes" % exp.dock.one("modes"))
        if exp.dock.one("energy_range"):
            detail.append("energy range %s kcal/mol" % exp.dock.one("energy_range"))
        if seeds:
            detail.append("seed%s %s" % ("s" if len(seeds) > 1 else "",
                                         ", ".join(str(s) for s in seeds)))
        sentences.append("Docking used %s." % "; ".join(detail))
        if len(seeds) == 1:
            sentences.append(
                "A single seed was used, so reported values are reproducible "
                "rather than replicate means.")

    got, limit = exp.control_result(), exp.control_limit()
    if got is not None and limit is not None:
        ligand = exp.control.one("redock") or "the co-crystallised ligand"
        control_prep = (exp.control.one("prepare") or "").strip()
        screen_prep = (receptor.one("prepare") or "").strip() if receptor else ""
        differs = (control_prep and screen_prep
                   and control_prep.lower() != screen_prep.lower())
        control_box = exp.control.one("box_size")
        control_exhaustiveness = exp.control.one("exhaustiveness")

        if differs:
            # The paragraph must not assert what the validator refuses. R106.
            qualifier = "on a receptor prepared %s" % control_prep
            if control_box:
                qualifier += ", in a %s A box" % control_box.split()[0]
            if control_exhaustiveness:
                qualifier += " at exhaustiveness %s" % control_exhaustiveness
            sentences.append(
                "Redocking of %s %s gave redock RMSD %.2f A, within the "
                "accepted %.1f A criterion." % (ligand, qualifier, got, limit))
            sentences.append(
                "This control was run under a different receptor preparation "
                "from the screen, which used %s, and it therefore validates "
                "pose recovery rather than the screening protocol itself."
                % screen_prep)
        else:
            sentences.append(
                "The protocol was validated by redocking %s, giving redock "
                "RMSD %.2f A against the accepted %.1f A criterion."
                % (ligand, got, limit))
    else:
        sentences.append(
            "The redocking control has not yet been recorded, so these values "
            "are not presented as validated.")

    scope = exp.study.one("claim_scope") if exp.study else None
    if scope:
        sentences.append(
            "Compounds were selected at the level of %s; they are candidate "
            "compounds documented across the genus rather than isolated from "
            "a single species, and binding is predicted in silico rather than "
            "demonstrated." % scope)

    if exp.replication and exp.reruns():
        original, rerun = exp.affinities(), exp.reruns()
        shared = [n for n in rerun if n in original]
        if shared:
            worst = max(abs(rerun[n] - original[n]) for n in shared)
            kind = (exp.replication.one("kind") or "").replace("_", "-")
            sentences.append(
                "A %s replication of %d compounds through %s agreed to within "
                "%.3f kcal/mol."
                % (kind or "repeat", len(shared),
                   exp.replication.one("pipeline") or "the same pipeline", worst))

    paragraph = " ".join(sentences)
    return "\n".join(_wrap(paragraph, ""))


def write_methods(experiment, path):
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(methods_paragraph(experiment))
        handle.write("\n")
    return path


def default_output_dir(experiment):
    base = os.path.dirname(os.path.abspath(experiment.path))
    stem = os.path.splitext(os.path.basename(experiment.path))[0]
    return os.path.join(base, stem + ".out")
