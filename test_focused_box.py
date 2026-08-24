#!/usr/bin/env python3
"""
Is the 30 A box the reason independently built ligands never recover?

WHAT IS ESTABLISHED
    Only crystal-derived coordinates recover their pose in this protocol:

      glucose, crystal-derived            0.518 A   PASS
      glucose, SMILES, twist-boat         5.358 A   fail
      glucose, SMILES, correct chair      5.919 A   fail
      isomaltose, SMILES, 1 boat ring     2.895 A   fail
      isomaltose, SMILES, all chairs      9.773 A   fail

    Ring geometry is not the cause. Search budget is not the cause; glucose
    was identical at exhaustiveness 32, 128 and 512.

THE UNTESTED VARIABLE
    Every one of those independent-ligand runs used the 30 A SCREENING box.
    The only configuration that ever passed, the published validation, used an
    18 A FOCUSED box, and it also used crystal coordinates. So box size and
    ligand origin have never been separated.

    A 30 A cube on this enzyme is close to blind docking. If a focused box
    recovers the pose from an independently built ligand, the protocol is fine
    and the box was simply too large for pose work. If it still fails, the
    protocol genuinely cannot recover a pose it was not handed.

    This is the last cheap experiment that distinguishes those.
"""

import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from talanai import control as tal_control   # noqa: E402

ASSETS = (r"D:\THESIS_VSC\BALAKATDBV2\talanaidock\dist\TalanaiDock\app"
          r"\docking_assets")
DATA = os.path.join(ASSETS, "docking_data")
VINA = os.path.join(ASSETS, "vina.exe")
RECEPTOR = os.path.join(DATA, "prepared", "receptor.pdbqt")
INPUTS = os.path.join(HERE, "validation-inputs")
OUT = os.path.join(HERE, "validation-run", "focused-box")

CENTER = ["21.52", "-7.7", "23.55"]
THRESHOLD = 2.0

CASES = [
    ("glucose crystal-derived", os.path.join(DATA, "validation", "glc_ligand.pdbqt"),
     os.path.join(DATA, "validation", "glc_3a4a_601.pdb")),
    ("glucose SMILES, chair", os.path.join(INPUTS, "glc_multiconf.pdbqt"),
     os.path.join(DATA, "validation", "glc_3a4a_601.pdb")),
    ("isomaltose SMILES, chairs", os.path.join(INPUTS, "isomaltose_ringaware.pdbqt"),
     os.path.join(INPUTS, "isomaltose_3axh_ref.pdb")),
]

BOXES = [18, 24, 30]


def rmsd_values(pose_file, reference):
    from rdkit import Chem, RDLogger
    from rdkit.Chem import rdMolAlign
    from meeko import PDBQTMolecule, RDKitMolCreate
    RDLogger.DisableLog("rdApp.*")
    ref = Chem.MolFromPDBFile(reference, removeHs=True, sanitize=False)
    mol = Chem.RemoveHs(RDKitMolCreate.from_pdbqt_mol(
        PDBQTMolecule.from_file(pose_file, skip_typing=True))[0])
    if ref is None or mol.GetNumAtoms() != ref.GetNumAtoms():
        return None
    return [round(rdMolAlign.CalcRMS(mol, ref, prbId=c, refId=0,
                                     maxMatches=1000000), 3)
            for c in range(mol.GetNumConformers())]


def main():
    os.makedirs(OUT, exist_ok=True)
    print("")
    print("  Does a focused box recover an independently built ligand?")
    print("  Receptor %s, exhaustiveness 32, seed 42."
          % os.path.basename(RECEPTOR))
    print("  " + "-" * 68)
    print("  %-26s %5s %9s %9s %9s"
          % ("ligand", "box", "score", "rank1", "best/9"))

    for label, ligand, reference in CASES:
        if not os.path.isfile(ligand):
            print("  %-26s MISSING" % label)
            continue
        for box in BOXES:
            tag = "%s_box%d" % (os.path.splitext(os.path.basename(ligand))[0], box)
            pose = os.path.join(OUT, tag + ".pdbqt")
            command = [VINA, "--receptor", RECEPTOR, "--ligand", ligand,
                       "--center_x", CENTER[0], "--center_y", CENTER[1],
                       "--center_z", CENTER[2],
                       "--size_x", str(box), "--size_y", str(box),
                       "--size_z", str(box),
                       "--exhaustiveness", "32", "--num_modes", "9",
                       "--energy_range", "3", "--seed", "42", "--cpu", "4",
                       "--out", pose]
            started = time.time()
            done = subprocess.run(command, capture_output=True, text=True,
                                  timeout=21600)
            if done.returncode != 0 or not os.path.isfile(pose):
                print("  %-26s %5d   FAILED" % (label, box), flush=True)
                continue
            score = tal_control.parse_score(done.stdout)
            values = rmsd_values(pose, reference)
            record = {"ligand": ligand, "reference": reference, "box": box,
                      "exhaustiveness": 32, "seed": 42,
                      "rank1_score_kcal_mol": score, "per_pose_rmsd_A": values,
                      "rmsd_method": "rdkit CalcRMS, no superposition",
                      "seconds": round(time.time() - started, 1),
                      "command": " ".join(command), "talanai_version": "0.1.0"}
            with open(os.path.join(OUT, tag + ".json"), "w",
                      encoding="utf-8") as handle:
                json.dump(record, handle, indent=2, sort_keys=True)
            if not values:
                print("  %-26s %5d %9.3f    rmsd n/a" % (label, box, score),
                      flush=True)
                continue
            flag = "  PASS" if values[0] < THRESHOLD else ""
            print("  %-26s %5d %9.3f %9.3f %9.3f%s"
                  % (label, box, score, values[0], min(values), flag), flush=True)
        print("")

    print("  " + "-" * 68)
    print("  If the SMILES-built ligands pass at 18 A but not at 30 A, the box")
    print("  was too large for pose work and the protocol is sound. If they")
    print("  fail at every box size while the crystal-derived one passes at")
    print("  every box size, the protocol cannot recover a pose it was not")
    print("  handed, and pose-level claims cannot be supported at all.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
