#!/usr/bin/env python3
"""
Re-prepare and re-dock the three compounds with defective rings.

TARGETS, from audit_rings_rigorous.py
    Ursolic acid   3 of 5 saturated rings in the boat family, published -8.105
    Isovitexin     its single sugar ring, published -8.024
    Acarbose       one sugar ring, the reference for every comparison

    Vina holds rings rigid, so none of these could have relaxed during docking.

TWO ARMS, because they answer different questions

    ARM A  raw receptor, exhaustiveness 8, seed 42. Exactly the thesis
           configuration, run with the OLD and the NEW preparation. Everything
           is held constant except ligand preparation, so the difference is
           attributable to preparation alone and is directly comparable with
           the published numbers.

    ARM B  Meeko-prepared receptor, exhaustiveness 32, three seeds, new
           preparation only. Not comparable with the published numbers, since
           the receptor and the search both differ, but it is the configuration
           a corrected study would actually use.

Nothing here writes to D:\\BALAKATDBV2. Published values are never overwritten.
"""

import json
import os
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from talanai import control as tal_control          # noqa: E402
from prep_multiconf import prepare, ring_signature  # noqa: E402

ASSETS = (r"D:\THESIS_VSC\BALAKATDBV2\talanaidock\dist\TalanaiDock\app"
          r"\docking_assets")
DATA = os.path.join(ASSETS, "docking_data")
VINA = os.path.join(ASSETS, "vina.exe")
RAW = os.path.join(DATA, "receptor_clean.pdb")
MEEKO = os.path.join(DATA, "prepared", "receptor.pdbqt")

HERE = os.path.dirname(os.path.abspath(__file__))
INPUTS = os.path.join(HERE, "validation-inputs", "reprepared")
OUT = os.path.join(HERE, "validation-run", "reprep")

CENTER = ["21.52", "-7.7", "23.55"]
BOX = ["30", "30", "30"]

TARGETS = [
    {"name": "Ursolic acid", "stem": "ursolic", "published": -8.105,
     "old": os.path.join(DATA, "ligands", "ursolic.pdbqt"),
     "smiles": "CC1CCC2(C(=O)O)CCC3(C)C(=CCC4C5(C)CCC(O)C(C)(C)C5CCC43C)C2C1C"},
    {"name": "Isovitexin", "stem": "isovitexin", "published": -8.024,
     "old": os.path.join(DATA, "ligands", "isovitexin.pdbqt"),
     "smiles": "O=c1cc(-c2ccc(O)cc2)oc2cc(O)c(C3OC(CO)C(O)C(O)C3O)c(O)c12"},
    {"name": "Acarbose (ref)", "stem": "acarbose", "published": -6.660,
     "old": os.path.join(DATA, "reference", "acarbose.pdbqt"),
     "smiles": ("C[C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@@H]2[C@H](O[C@@H]"
                "([C@@H]([C@H]2O)O)O[C@H]([C@@H](CO)O)[C@@H]([C@H](C=O)O)O)CO)"
                "O)O)N[C@H]3C=C([C@H]([C@@H]([C@H]3O)O)O)CO")},
]


def dock(receptor, ligand, exhaustiveness, seed, tag):
    os.makedirs(OUT, exist_ok=True)
    pose = os.path.join(OUT, "%s.pdbqt" % tag)
    command = [VINA, "--receptor", receptor, "--ligand", ligand,
               "--center_x", CENTER[0], "--center_y", CENTER[1],
               "--center_z", CENTER[2],
               "--size_x", BOX[0], "--size_y", BOX[1], "--size_z", BOX[2],
               "--exhaustiveness", str(exhaustiveness), "--num_modes", "9",
               "--energy_range", "3", "--seed", str(seed), "--cpu", "4",
               "--out", pose]
    started = time.time()
    done = subprocess.run(command, capture_output=True, text=True, timeout=21600)
    record = {"tag": tag, "receptor": receptor, "ligand": ligand,
              "exhaustiveness": exhaustiveness, "seed": seed,
              "seconds": round(time.time() - started, 1),
              "returncode": done.returncode, "command": " ".join(command),
              "pose_file": pose if os.path.isfile(pose) else None,
              "talanai_version": "0.1.0"}
    if done.returncode != 0:
        record["error"] = (done.stderr or done.stdout)[-800:]
        return None, record
    score = tal_control.parse_score(done.stdout)
    record["rank1_score_kcal_mol"] = score
    with open(os.path.join(OUT, tag + ".json"), "w", encoding="utf-8") as handle:
        json.dump(record, handle, indent=2, sort_keys=True)
    return score, record


def main():
    os.makedirs(INPUTS, exist_ok=True)
    os.makedirs(OUT, exist_ok=True)

    print("")
    print("  STEP 1  re-prepare with a 50-conformer search")
    print("  " + "-" * 72)
    for target in TARGETS:
        pdbqt, mol, chosen, energies = prepare(target["smiles"], 50)
        path = os.path.join(INPUTS, target["stem"] + ".pdbqt")
        with open(path, "w", encoding="utf-8") as handle:
            handle.write(pdbqt)
        target["new"] = path
        rings = ring_signature(mol)
        boats = sum(1 for _s, mean_abs in rings if mean_abs < 45.0)
        print("  %-16s conformer %2d of %d, energy %8.2f (worst %8.2f)"
              % (target["name"], chosen + 1, len(energies),
                 min(energies), max(energies)))
        print("  %-16s %d six-ring(s), %d below chair amplitude"
              % ("", len(rings), boats))
    print("")

    print("  STEP 2  ARM A: thesis configuration, raw receptor, exh 8, seed 42")
    print("  Only the ligand preparation differs. Comparable with published.")
    print("  " + "-" * 72)
    print("  %-16s %10s %10s %10s %10s"
          % ("compound", "published", "old prep", "new prep", "delta"))
    arm_a = []
    for target in TARGETS:
        old, _r = dock(RAW, target["old"], 8, 42, target["stem"] + "_raw_old")
        new, _r = dock(RAW, target["new"], 8, 42, target["stem"] + "_raw_new")
        if old is None or new is None:
            print("  %-16s FAILED" % target["name"])
            continue
        arm_a.append((target["name"], target["published"], old, new))
        print("  %-16s %10.3f %10.3f %10.3f %+10.3f"
              % (target["name"], target["published"], old, new, new - old),
              flush=True)
    print("")

    print("  STEP 3  ARM B: Meeko receptor, exh 32, seeds 42/43/44, new prep")
    print("  The configuration a corrected study would use.")
    print("  " + "-" * 72)
    print("  %-16s %9s %9s %9s %10s %8s"
          % ("compound", "seed 42", "seed 43", "seed 44", "mean", "spread"))
    for target in TARGETS:
        scores = []
        for seed in (42, 43, 44):
            score, _r = dock(MEEKO, target["new"], 32, seed,
                             "%s_meeko_new_s%d" % (target["stem"], seed))
            scores.append(score)
        clean = [s for s in scores if s is not None]
        if not clean:
            print("  %-16s all runs failed" % target["name"])
            continue
        mean = sum(clean) / len(clean)
        print("  %-16s %9s %9s %9s %10.3f %8.3f"
              % (target["name"],
                 *["%.3f" % s if s is not None else "-" for s in scores],
                 mean, max(clean) - min(clean)), flush=True)

    print("")
    print("  " + "-" * 72)
    print("  Arm A isolates ligand preparation. Arm B is not comparable with")
    print("  the published values: receptor and search both differ.")
    print("  No published number was modified.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
