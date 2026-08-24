#!/usr/bin/env python3
"""
The control battery from REDOCK-PROTOCOL.md section 4.

Runs three controls, in order, and decides whether the validation screen is
allowed to proceed. Control 3 is the gate; controls 1 and 2 are information.

    python run_controls.py

WHY RDKIT APPEARS HERE AND NOT IN talanai/
    Talanai's own rmsd() is standard-library only and returns a LOWER BOUND
    when atom correspondence cannot be established, which is the case for
    every ligand here: the docking inputs name atoms C and O while the crystal
    references use C1, O6 and so on.

    A lower bound above the threshold is conclusive for FAILURE. A lower bound
    below the threshold is NOT a pass, only an absence of proof of failure. A
    gate therefore cannot be decided on it. RDKit's GetBestRMS does the
    symmetry-corrected graph matching that settles it, mirroring TalanaiDock's
    own validate_redock.py. Both numbers are recorded for every control.
"""

import json
import os
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from talanai import control as tal_control   # noqa: E402

ASSETS = (r"D:\THESIS_VSC\BALAKATDBV2\talanaidock\dist\TalanaiDock\app"
          r"\docking_assets")
DATA = os.path.join(ASSETS, "docking_data")
VINA = os.path.join(ASSETS, "vina.exe")
RECEPTOR = os.path.join(DATA, "prepared", "receptor.pdbqt")

HERE = os.path.dirname(os.path.abspath(__file__))
INPUTS = os.path.join(HERE, "validation-inputs")
OUT = os.path.join(HERE, "validation-run", "controls")

CENTER = ["21.52", "-7.7", "23.55"]
BOX = ["30", "30", "30"]
EXHAUSTIVENESS = "32"
SEED = "42"
THRESHOLD = 2.0

CONTROLS = [
    {
        "id": "1_crystal_glucose",
        "gate": False,
        "why": "continuity with validation.json. Cannot fail; five of its six "
               "torsions rotate only a hydroxyl hydrogen.",
        "ligand": os.path.join(DATA, "validation", "glc_ligand.pdbqt"),
        "reference": os.path.join(DATA, "validation", "glc_3a4a_601.pdb"),
    },
    {
        "id": "2_reembedded_glucose",
        "gate": False,
        "why": "tests that the ligand pipeline produces a chemically correct "
               "molecule. Not a harder search.",
        "ligand": os.path.join(INPUTS, "glc_reembedded.pdbqt"),
        "reference": os.path.join(DATA, "validation", "glc_3a4a_601.pdb"),
    },
    {
        "id": "3_isomaltose_crossdock",
        "gate": True,
        "why": "THE GATE. Cognate alpha-1,6 substrate of this isomaltase, 23 "
               "heavy atoms, cross-docked from 3AXH so it also tests induced "
               "fit. Reference is superposed into 3A4A's frame.",
        "ligand": os.path.join(INPUTS, "isomaltose.pdbqt"),
        "reference": os.path.join(INPUTS, "isomaltose_3axh_ref.pdb"),
    },
]


def symmetry_corrected_rmsd(pose_pdbqt, reference_pdb):
    """
    RDKit symmetry-corrected heavy-atom RMSD, no superposition.

    Returns (value, note). Returns (None, reason) rather than a number when it
    cannot be computed, because a gate must never be decided on a guess.
    """
    try:
        from rdkit import Chem, RDLogger
        from rdkit.Chem import AllChem
        RDLogger.DisableLog("rdApp.*")
    except ImportError as exc:
        return None, "rdkit unavailable: %s" % exc

    try:
        from meeko import PDBQTMolecule, RDKitMolCreate
        pose = RDKitMolCreate.from_pdbqt_mol(
            PDBQTMolecule.from_file(pose_pdbqt, skip_typing=True))[0]
        pose = Chem.RemoveHs(pose)
    except Exception as exc:
        return None, "could not read pose: %s: %s" % (type(exc).__name__, exc)

    ref = Chem.MolFromPDBFile(reference_pdb, removeHs=True, sanitize=False)
    if ref is None:
        return None, "could not read reference %s" % os.path.basename(reference_pdb)

    if pose.GetNumAtoms() != ref.GetNumAtoms():
        return None, "atom counts differ (%d pose, %d reference)" % (
            pose.GetNumAtoms(), ref.GetNumAtoms())

    # TWO DIFFERENT QUESTIONS, AND ONLY ONE OF THEM IS THE CONTROL.
    #
    #   CalcRMS      symmetry-corrected, NO superposition. "Did the pose land
    #                in the right PLACE?" This is the redocking control.
    #   GetBestRMS   symmetry-corrected, superimposes first. "Did the pose find
    #                the right SHAPE, wherever it ended up?" Diagnostic only.
    #
    # An earlier version of this file used GetBestRMS for the verdict. That was
    # wrong: a pose 5 A away with perfect internal geometry scores near zero
    # under GetBestRMS and would have been passed as a valid redock.
    from rdkit.Chem import rdMolAlign
    try:
        placed = rdMolAlign.CalcRMS(pose, ref, maxMatches=1000000)
    except Exception as exc:
        return None, "CalcRMS failed: %s: %s" % (type(exc).__name__, exc)
    try:
        shape = AllChem.GetBestRMS(pose, ref, maxMatches=1000000)
    except Exception:
        shape = None
    return placed, ("rdkit CalcRMS, symmetry-corrected, NO superposition"
                    + ("; shape-only GetBestRMS %.3f A" % shape
                       if shape is not None else ""))


def dock(spec):
    os.makedirs(OUT, exist_ok=True)
    pose_out = os.path.join(OUT, spec["id"] + "_pose.pdbqt")
    command = [VINA, "--receptor", RECEPTOR, "--ligand", spec["ligand"],
               "--center_x", CENTER[0], "--center_y", CENTER[1],
               "--center_z", CENTER[2],
               "--size_x", BOX[0], "--size_y", BOX[1], "--size_z", BOX[2],
               "--exhaustiveness", EXHAUSTIVENESS, "--num_modes", "9",
               "--energy_range", "3", "--seed", SEED, "--cpu", "4",
               "--out", pose_out]
    started = time.time()
    done = subprocess.run(command, capture_output=True, text=True, timeout=7200)
    elapsed = time.time() - started

    record = {
        "control": spec["id"],
        "is_gate": spec["gate"],
        "purpose": spec["why"],
        "receptor": RECEPTOR,
        "ligand": spec["ligand"],
        "reference": spec["reference"],
        "box_center": CENTER, "box_size": BOX,
        "exhaustiveness": int(EXHAUSTIVENESS), "seed": int(SEED),
        "command": " ".join(command),
        "returncode": done.returncode,
        "seconds": round(elapsed, 1),
        "pose_file": pose_out if os.path.isfile(pose_out) else None,
        "talanai_version": "0.1.0",
    }
    if done.returncode != 0:
        record["error"] = (done.stderr or done.stdout)[-1500:]
        return record

    record["rank1_score_kcal_mol"] = tal_control.parse_score(done.stdout)

    crystal = tal_control.read_heavy_atoms(spec["reference"])
    pose = tal_control.read_heavy_atoms(pose_out, model=1)
    bound, method, matched, exact = tal_control.rmsd(pose, crystal)
    record["talanai_rmsd_A"] = round(bound, 3) if bound is not None else None
    record["talanai_rmsd_method"] = method
    record["talanai_rmsd_is_exact"] = exact

    value, note = symmetry_corrected_rmsd(pose_out, spec["reference"])
    record["rdkit_rmsd_A"] = round(value, 3) if value is not None else None
    record["rdkit_rmsd_note"] = note
    return record


def verdict(record):
    """Only an exact or symmetry-corrected number may confirm a pass."""
    bound = record.get("talanai_rmsd_A")
    exact = record.get("talanai_rmsd_is_exact")
    rdkit = record.get("rdkit_rmsd_A")

    if bound is not None and bound >= THRESHOLD:
        return "FAIL", "lower bound %.3f A exceeds %.1f A, conclusive" % (bound, THRESHOLD)
    if rdkit is not None:
        return ("PASS" if rdkit < THRESHOLD else "FAIL",
                "symmetry-corrected %.3f A against %.1f A" % (rdkit, THRESHOLD))
    if exact and bound is not None:
        return ("PASS" if bound < THRESHOLD else "FAIL",
                "exact %.3f A against %.1f A" % (bound, THRESHOLD))
    return "INCONCLUSIVE", ("only a lower bound (%s) is available and it is under "
                           "the threshold, which is not a pass" % bound)


def main():
    os.makedirs(OUT, exist_ok=True)
    print("")
    print("  Control battery, REDOCK-PROTOCOL.md section 4")
    print("  receptor %s" % os.path.basename(RECEPTOR))
    print("  box %s A at (%s), exhaustiveness %s, seed %s"
          % (BOX[0], ", ".join(CENTER), EXHAUSTIVENESS, SEED))
    print("  " + "-" * 68)

    results = []
    for spec in CONTROLS:
        if not os.path.isfile(spec["ligand"]):
            print("  %-26s MISSING LIGAND %s" % (spec["id"], spec["ligand"]))
            continue
        print("  %-26s running ..." % spec["id"], flush=True)
        record = dock(spec)
        state, why = verdict(record)
        record["verdict"] = state
        record["verdict_reason"] = why
        results.append(record)

        with open(os.path.join(OUT, spec["id"] + ".json"), "w",
                  encoding="utf-8") as handle:
            json.dump(record, handle, indent=2, sort_keys=True)

        if record.get("error"):
            print("      FAILED TO RUN: %s" % record["error"].splitlines()[-1][:120])
            continue
        print("      score        %.3f kcal/mol" % record["rank1_score_kcal_mol"])
        print("      talanai RMSD %s A%s"
              % (record["talanai_rmsd_A"], "" if record["talanai_rmsd_is_exact"]
                 else "  (lower bound)"))
        print("      rdkit RMSD   %s A" % (record["rdkit_rmsd_A"]
                                           if record["rdkit_rmsd_A"] is not None
                                           else "not computed: " + record["rdkit_rmsd_note"][:60]))
        print("      %-6s %s%s" % (state, why, "   <-- GATE" if spec["gate"] else ""))
        print("      %.0f s" % record["seconds"])
        print("")

    print("  " + "-" * 68)
    gate = [r for r in results if r["is_gate"]]
    if not gate:
        print("  NO GATE RAN. The screen must not proceed.")
        return 1
    state = gate[0]["verdict"]
    if state == "PASS":
        print("  GATE PASSED. The screen may proceed.")
        return 0
    print("  GATE %s. The screen must not proceed." % state)
    print("  %s" % gate[0]["verdict_reason"])
    return 1


if __name__ == "__main__":
    sys.exit(main())
