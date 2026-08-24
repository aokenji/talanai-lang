#!/usr/bin/env python3
"""
Re-run the control gate with a ring-aware isomaltose.

WHY
    Control 3, the gate, failed at 2.895 A against a 2.0 A bar. But the
    isomaltose used was built with the single-conformer method, and
    audit_rings_rigorous.py flagged one of its two saturated rings as a boat.

    Vina cannot fix a ring. So the gate may have failed for exactly the reason
    acarbose scored badly, and it has to be retested with a correctly prepared
    ligand before the protocol is judged.

    This does not re-open the self-docking-bias finding, which stands on the
    glucose experiments. It only asks whether the isomaltose gate is passable.
"""

import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from talanai import control as tal_control      # noqa: E402
import prep_ringaware as ringaware              # noqa: E402

ASSETS = (r"D:\THESIS_VSC\BALAKATDBV2\talanaidock\dist\TalanaiDock\app"
          r"\docking_assets")
DATA = os.path.join(ASSETS, "docking_data")
VINA = os.path.join(ASSETS, "vina.exe")
RECEPTOR = os.path.join(DATA, "prepared", "receptor.pdbqt")

INPUTS = os.path.join(HERE, "validation-inputs")
REFERENCE = os.path.join(INPUTS, "isomaltose_3axh_ref.pdb")
OUT = os.path.join(HERE, "validation-run", "gate-retest")

CENTER = ["21.52", "-7.7", "23.55"]
BOX = ["30", "30", "30"]
THRESHOLD = 2.0

# isomaltose, PubChem CID 439193, alpha-1,6 linked
ISOMALTOSE = ("OC[C@H]1O[C@H](OC[C@H]2O[C@H](O)[C@H](O)[C@@H](O)[C@@H]2O)"
              "[C@H](O)[C@@H](O)[C@@H]1O")


def rmsd_no_superposition(pose_file, reference):
    from rdkit import Chem, RDLogger
    from rdkit.Chem import rdMolAlign
    from meeko import PDBQTMolecule, RDKitMolCreate
    RDLogger.DisableLog("rdApp.*")
    ref = Chem.MolFromPDBFile(reference, removeHs=True, sanitize=False)
    mol = Chem.RemoveHs(RDKitMolCreate.from_pdbqt_mol(
        PDBQTMolecule.from_file(pose_file, skip_typing=True))[0])
    if ref is None or mol.GetNumAtoms() != ref.GetNumAtoms():
        return None, "atom counts differ or reference unreadable"
    return [round(rdMolAlign.CalcRMS(mol, ref, prbId=c, refId=0,
                                     maxMatches=1000000), 3)
            for c in range(mol.GetNumConformers())], "rdkit CalcRMS, no superposition"


def main():
    os.makedirs(OUT, exist_ok=True)
    print("")
    print("  Gate retest: isomaltose, ring-aware preparation")
    print("  Previous gate result: 2.895 A, FAIL, with a single-conformer ligand")
    print("  whose audit showed 1 of 2 saturated rings in the boat family.")
    print("  " + "-" * 68)

    result = ringaware.prepare_ring_aware(ISOMALTOSE, 300)
    ligand = os.path.join(INPUTS, "isomaltose_ringaware.pdbqt")
    with open(ligand, "w", encoding="utf-8") as handle:
        handle.write(result["pdbqt"])

    thetas = result["thetas"]
    ok_rings = ringaware.all_chairs(thetas) if thetas else False
    print("  preparation  %d of %d conformers all-chair, picked #%d"
          % (result["n_qualifying"], result["n_conformers"], result["chosen"] + 1))
    print("  ring theta   %s  ->  %s"
          % (", ".join("%.1f" % t for t in thetas),
             "ALL CHAIRS" if ok_rings else "STILL DEFECTIVE"))
    print("  energy       %.2f (global best %.2f)"
          % (result["energy"], result["best_energy"]))
    print("")

    if not ok_rings:
        print("  Rings still wrong. Not docking. Fix preparation first.")
        return 1

    pose = os.path.join(OUT, "isomaltose_ringaware_e32_s42.pdbqt")
    command = [VINA, "--receptor", RECEPTOR, "--ligand", ligand,
               "--center_x", CENTER[0], "--center_y", CENTER[1],
               "--center_z", CENTER[2],
               "--size_x", BOX[0], "--size_y", BOX[1], "--size_z", BOX[2],
               "--exhaustiveness", "32", "--num_modes", "9",
               "--energy_range", "3", "--seed", "42", "--cpu", "4",
               "--out", pose]
    print("  docking ...", flush=True)
    started = time.time()
    done = subprocess.run(command, capture_output=True, text=True, timeout=21600)
    elapsed = time.time() - started

    if done.returncode != 0 or not os.path.isfile(pose):
        print("  FAILED: %s" % (done.stderr or done.stdout)[-200:])
        return 1

    score = tal_control.parse_score(done.stdout)
    values, method = rmsd_no_superposition(pose, REFERENCE)
    record = {
        "control": "3_isomaltose_crossdock_ringaware", "is_gate": True,
        "ligand": ligand, "reference": REFERENCE, "receptor": RECEPTOR,
        "preparation": "ring-aware: lowest energy among all-chair conformers",
        "ring_theta": thetas, "exhaustiveness": 32, "seed": 42,
        "rank1_score_kcal_mol": score, "per_pose_rmsd_A": values,
        "rmsd_method": method, "threshold_A": THRESHOLD,
        "previous_gate_rmsd_A": 2.895,
        "seconds": round(elapsed, 1), "command": " ".join(command),
        "talanai_version": "0.1.0",
    }
    if values:
        record["rank1_rmsd_A"] = values[0]
        record["best_of_ensemble_rmsd_A"] = min(values)
    with open(os.path.join(OUT, "gate_retest.json"), "w", encoding="utf-8") as h:
        json.dump(record, h, indent=2, sort_keys=True)

    print("")
    print("  score        %.3f kcal/mol" % score)
    if not values:
        print("  RMSD         not computable: %s" % method)
        return 1
    print("  rank1 RMSD   %.3f A   (was 2.895 A)" % values[0])
    print("  best of 9    %.3f A" % min(values))
    print("  took         %.0f s" % elapsed)
    print("")
    print("  " + "-" * 68)
    if values[0] < THRESHOLD:
        print("  GATE PASSES at %.3f A. The earlier failure was the ligand's"
              % values[0])
        print("  ring geometry, not the protocol. The screen may proceed.")
        return 0
    print("  GATE STILL FAILS at %.3f A against %.1f A." % (values[0], THRESHOLD))
    print("  Ring geometry was not the cause. The protocol itself does not")
    print("  recover a cross-docked cognate substrate.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
