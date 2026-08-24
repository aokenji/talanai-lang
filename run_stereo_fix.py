#!/usr/bin/env python3
"""
Stereochemistry correction: re-dock the three compounds whose 2026-08-04
corrected-screen SMILES had unspecified stereocenters.

WHAT WAS WRONG, found 2026-08-11 during defence-panel review
    Ursolic Acid, Vitexin and Isovitexin were docked with SMILES carrying
    ZERO of their stereocenters specified (10, 5 and 5 respectively; verified
    with RDKit AssignStereochemistry against the real PubChem structures).
    RDKit's conformer generator therefore picked essentially arbitrary
    configurations at every one of those centers. What was docked as "Ursolic
    Acid" etc. may not have been the natural-product stereoisomer at all.

    Rutin, Betulinic Acid, Spinosin, Oleanolic Acid and Acarbose were checked
    at the same time and are fully specified; they are unaffected and are not
    re-run here. Luteolin, Quercetin and Kaempferol have no stereocenters.

FIX
    SMILES re-pulled from PubChem by the CID already recorded in
    src/data/compounds.js (not re-guessed), verified fully stereo-specified
    with RDKit before use. Everything else held IDENTICAL to
    run_full_screen.py: same receptor, same box, same exhaustiveness 32,
    same three seeds 42/43/44, same ring-aware conformer selection.

    python run_stereo_fix.py
"""

import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from talanai import control as tal_control       # noqa: E402
import prep_ringaware as ringaware               # noqa: E402

ASSETS = (r"D:\THESIS_VSC\BALAKATDBV2\talanaidock\dist\TalanaiDock\app"
          r"\docking_assets")
DATA = os.path.join(ASSETS, "docking_data")
VINA = os.path.join(ASSETS, "vina.exe")
RECEPTOR = os.path.join(DATA, "prepared", "receptor.pdbqt")

INPUTS = os.path.join(HERE, "validation-inputs", "stereo-fix")
OUT = os.path.join(HERE, "validation-run", "stereo-fix")

CENTER = ["21.52", "-7.7", "23.55"]
BOX = ["30", "30", "30"]
EXHAUSTIVENESS = 32
SEEDS = (42, 43, 44)

ACARBOSE = -8.706  # the corrected-screen reference, unaffected by this fix

# name, stem, pubchem_cid (as used in compounds.js), formula, conformers,
# fully-stereo-specified SMILES (fetched from PubChem 2026-08-11, verified
# with RDKit AssignStereochemistry before being placed here)
COMPOUNDS = [
    ("Ursolic Acid", "ursolic_v2", 64945, "C30H48O3", 150,
     "C[C@@H]1CC[C@@]2(CC[C@@]3(C(=CC[C@H]4[C@]3(CC[C@@H]5[C@@]4(CC[C@@H]"
     "(C5(C)C)O)C)C)[C@@H]2[C@H]1C)C)C(=O)O"),
    ("Vitexin", "vitexin_v2", 5280441, "C21H20O10", 150,
     "C1=CC(=CC=C1C2=CC(=O)C3=C(O2)C(=C(C=C3O)O)[C@H]4[C@@H]([C@H]([C@@H]"
     "([C@H](O4)CO)O)O)O)O"),
    ("Isovitexin", "isovitexin_v2", 162350, "C21H20O10", 150,
     "C1=CC(=CC=C1C2=CC(=O)C3=C(O2)C=C(C(=C3O)[C@H]4[C@@H]([C@H]([C@@H]"
     "([C@H](O4)CO)O)O)O)O)O"),
]


def dock(ligand, seed, tag):
    pose = os.path.join(OUT, tag + ".pdbqt")
    command = [VINA, "--receptor", RECEPTOR, "--ligand", ligand,
               "--center_x", CENTER[0], "--center_y", CENTER[1],
               "--center_z", CENTER[2],
               "--size_x", BOX[0], "--size_y", BOX[1], "--size_z", BOX[2],
               "--exhaustiveness", str(EXHAUSTIVENESS), "--num_modes", "9",
               "--energy_range", "3", "--seed", str(seed), "--cpu", "4",
               "--out", pose]
    started = time.time()
    done = subprocess.run(command, capture_output=True, text=True, timeout=21600)
    record = {"tag": tag, "ligand": ligand, "receptor": RECEPTOR,
              "preparation": "ring-aware, fully stereo-specified SMILES (fixed 2026-08-11)",
              "box_center": CENTER, "box_size": BOX,
              "exhaustiveness": EXHAUSTIVENESS, "seed": seed,
              "seconds": round(time.time() - started, 1),
              "returncode": done.returncode, "command": " ".join(command),
              "talanai_version": "0.1.0",
              "not_validated": ("the redocking gate does not pass for this "
                                "configuration; ranking only, no pose claims"),
              "pose_file": pose if os.path.isfile(pose) else None}
    if done.returncode != 0:
        record["error"] = (done.stderr or done.stdout)[-800:]
    else:
        record["rank1_score_kcal_mol"] = tal_control.parse_score(done.stdout)
    with open(os.path.join(OUT, tag + ".json"), "w", encoding="utf-8") as handle:
        json.dump(record, handle, indent=2, sort_keys=True)
    return record.get("rank1_score_kcal_mol")


# Prior (broken) values, for the printed before/after comparison only.
PRIOR = {"Ursolic Acid": -8.717, "Vitexin": -9.490, "Isovitexin": -9.356}
PRIOR_SPREAD = {"Ursolic Acid": 0.035, "Vitexin": 0.043, "Isovitexin": 0.029}


def main():
    os.makedirs(INPUTS, exist_ok=True)
    os.makedirs(OUT, exist_ok=True)

    print("")
    print("  STEP 1  ring-aware prep, fully stereo-specified structures")
    print("  " + "-" * 74)
    ready = []
    for name, stem, cid, formula, confs, smiles in COMPOUNDS:
        result = ringaware.prepare_ring_aware(smiles, confs)
        thetas = result["thetas"]
        ok = ringaware.all_chairs(thetas) if thetas else True
        path = os.path.join(INPUTS, stem + ".pdbqt")
        with open(path, "w", encoding="utf-8") as handle:
            handle.write(result["pdbqt"])
        note = ("no saturated ring" if not thetas
                else ("all chairs" if ok else "*** STILL DEFECTIVE ***"))
        print("  %-16s %3d/%3d all-chair, %s" % (name, result["n_qualifying"],
                                                 result["n_conformers"], note))
        ready.append((name, stem, cid, formula, path))

    print("")
    print("  STEP 2  docking %d ligands x %d seeds at exhaustiveness %d"
          % (len(ready), len(SEEDS), EXHAUSTIVENESS))
    print("  " + "-" * 74)

    results = []
    for name, stem, cid, formula, path in ready:
        scores = []
        for seed in SEEDS:
            score = dock(path, seed, "%s_seed%d" % (stem, seed))
            if score is None:
                print("  %-16s seed %d FAILED" % (name, seed))
            else:
                scores.append(score)
                print("  %-16s seed %d  %8.3f" % (name, seed, score))
        if scores:
            mean = sum(scores) / len(scores)
            spread = max(scores) - min(scores)
            results.append({"name": name, "mean": round(mean, 3),
                            "spread": round(spread, 3)})

    print("")
    print("  RESULT  before / after, against Acarbose %.3f" % ACARBOSE)
    print("  " + "-" * 74)
    print("  %-16s %14s %14s   %s" % ("compound", "prior (broken)",
                                      "corrected", "verdict change"))
    for r in results:
        name = r["name"]
        half = r["spread"] / 2.0

        def verdict(mean, spread):
            h = spread / 2.0
            if mean + h < ACARBOSE - 0.005:
                return "beats"
            if mean - h > ACARBOSE + 0.005:
                return "loses to"
            return "indistinguishable from"

        prior_verdict = verdict(PRIOR[name], PRIOR_SPREAD[name])
        new_verdict = verdict(r["mean"], r["spread"])
        change = ("SAME (%s)" % new_verdict if prior_verdict == new_verdict
                  else "CHANGED: %s -> %s" % (prior_verdict, new_verdict))
        print("  %-16s %14.3f %14.3f   %s" % (name, PRIOR[name], r["mean"],
                                              change))

    summary = {"fix": "stereochemistry correction",
               "date": "2026-08-11",
               "protocol": {"receptor": RECEPTOR, "box_center": CENTER,
                            "box_size": BOX, "exhaustiveness": EXHAUSTIVENESS,
                            "seeds": list(SEEDS)},
               "reference_acarbose": ACARBOSE,
               "results": results, "prior": PRIOR}
    with open(os.path.join(OUT, "summary.json"), "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
    print("")
    print("  wrote %s" % os.path.join(OUT, "summary.json"))


if __name__ == "__main__":
    main()
