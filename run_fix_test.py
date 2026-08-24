#!/usr/bin/env python3
"""
Does properly preparing the ligand fix the redocking control?

Apples to apples against control 2: same receptor, same box, same
exhaustiveness 32, same seeds. The ONLY thing that changed is that the ligand
was built by embedding 50 conformers and keeping the lowest-energy one,
instead of embedding a single conformer and minimising it locally.

  original prep, exh 32, seed 42  ->  5.358 A   (twist-boat ring, FAIL)

If the fixed preparation lands under 2.0 A, then the diagnosis is proven and
the repair is a one-line change to the preparation pipeline.
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
REFERENCE = os.path.join(DATA, "validation", "glc_3a4a_601.pdb")

HERE = os.path.dirname(os.path.abspath(__file__))
LIGAND = os.path.join(HERE, "validation-inputs", "glc_multiconf.pdbqt")
OUT = os.path.join(HERE, "validation-run", "fixtest")

CENTER = ["21.52", "-7.7", "23.55"]
BOX = ["30", "30", "30"]
THRESHOLD = 2.0
# The escalation was previously run on the twist-boat file only. The correct
# chair has to be tested at the same effort before search can be ruled in or
# out: the earlier exhaustiveness 32 run tested the ligand, not the budget.
PLAN = [(32, 42), (32, 43), (32, 44), (128, 42), (128, 43), (512, 42)]


def per_pose_rmsd(pose_file, reference):
    from rdkit import Chem, RDLogger
    from rdkit.Chem import rdMolAlign
    from meeko import PDBQTMolecule, RDKitMolCreate
    RDLogger.DisableLog("rdApp.*")
    ref = Chem.MolFromPDBFile(reference, removeHs=True, sanitize=False)
    mol = Chem.RemoveHs(RDKitMolCreate.from_pdbqt_mol(
        PDBQTMolecule.from_file(pose_file, skip_typing=True))[0])
    return [round(rdMolAlign.CalcRMS(mol, ref, prbId=c, refId=0,
                                     maxMatches=1000000), 3)
            for c in range(mol.GetNumConformers())]


def main():
    os.makedirs(OUT, exist_ok=True)
    if not os.path.isfile(LIGAND):
        print("  missing %s. Run prep_multiconf.py first." % LIGAND)
        return 2

    print("")
    print("  Correct 4C1 chair ligand, escalating search effort")
    print("  receptor %s, box %s A" % (os.path.basename(RECEPTOR), BOX[0]))
    print("  twist-boat prep gave 5.358 A at exh 32 and 5.36 A at exh 512")
    print("  crystal-coordinate input gives 0.518 A at every effort")
    print("  " + "-" * 66)
    print("  %5s %5s %9s %9s %9s %5s"
          % ("exh", "seed", "score", "rank1", "best/9", "at#"))

    rows = []
    for EXHAUSTIVENESS, seed in PLAN:
        pose_out = os.path.join(OUT, "glc_multiconf_e%d_s%d.pdbqt"
                                % (EXHAUSTIVENESS, seed))
        command = [VINA, "--receptor", RECEPTOR, "--ligand", LIGAND,
                   "--center_x", CENTER[0], "--center_y", CENTER[1],
                   "--center_z", CENTER[2],
                   "--size_x", BOX[0], "--size_y", BOX[1], "--size_z", BOX[2],
                   "--exhaustiveness", str(EXHAUSTIVENESS), "--num_modes", "9",
                   "--energy_range", "3", "--seed", str(seed), "--cpu", "4",
                   "--out", pose_out]
        started = time.time()
        done = subprocess.run(command, capture_output=True, text=True, timeout=7200)
        elapsed = time.time() - started
        if done.returncode != 0 or not os.path.isfile(pose_out):
            print("  %5d   FAILED %s" % (seed, (done.stderr or "")[-60:]))
            continue

        score = tal_control.parse_score(done.stdout)
        values = per_pose_rmsd(pose_out, REFERENCE)
        best = min(values)
        record = {
            "ligand": LIGAND, "preparation": "EmbedMultipleConfs 50, lowest energy",
            "receptor": RECEPTOR, "box_center": CENTER, "box_size": BOX,
            "exhaustiveness": EXHAUSTIVENESS, "seed": seed,
            "rank1_score_kcal_mol": score, "per_pose_rmsd_A": values,
            "rank1_rmsd_A": values[0], "best_of_ensemble_rmsd_A": best,
            "best_pose_index": values.index(best) + 1,
            "rmsd_method": "rdkit CalcRMS, symmetry-corrected, no superposition",
            "threshold_A": THRESHOLD, "seconds": round(elapsed, 1),
            "command": " ".join(command), "talanai_version": "0.1.0",
            "comparison": "original single-conformer prep gave 5.358 A rank-1",
        }
        with open(os.path.join(OUT, "glc_multiconf_e%d_s%d.json"
                               % (EXHAUSTIVENESS, seed)), "w",
                  encoding="utf-8") as handle:
            json.dump(record, handle, indent=2, sort_keys=True)
        rows.append(record)
        print("  %5d %5d %9.3f %9.3f %9.3f %5d"
              % (EXHAUSTIVENESS, seed, score, values[0], best,
                 record["best_pose_index"]), flush=True)

    print("  " + "-" * 66)
    if not rows:
        print("  No runs completed.")
        return 1
    passed = [r for r in rows if r["rank1_rmsd_A"] < THRESHOLD]
    print("")
    if len(passed) == len(rows):
        print("  FIXED. Rank-1 recovers the crystal pose in all %d seeds,"
              % len(rows))
        print("  against 5.358 A for the single-conformer preparation.")
        print("  The failure was ligand preparation, not the receptor, not the")
        print("  box, and not the search budget.")
    elif passed:
        print("  PARTIALLY FIXED: %d of %d seeds recover the pose."
              % (len(passed), len(rows)))
    else:
        print("  NOT FIXED by preparation alone. Rank-1 still misses at %.1f A."
              % THRESHOLD)
        print("  The ring conformation was necessary but not sufficient.")
    print("")
    return 0


if __name__ == "__main__":
    sys.exit(main())
