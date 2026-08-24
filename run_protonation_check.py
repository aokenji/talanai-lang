#!/usr/bin/env python3
"""
Ligand protonation-state check: re-dock the three carboxylic-acid triterpenes
(betulinic acid, oleanolic acid, ursolic acid) with their C-28 carboxylic
acid deprotonated (COO-) instead of neutral (COOH), and compare against the
live site's current values.

WHY THIS IS BEING CHECKED, flagged 2026-08-11 during defence-panel review
    Betulinic, oleanolic and ursolic acid each carry one carboxylic acid
    group. The pipeline (this file's SMILES source, run_full_screen.py /
    run_stereo_fix.py) has always written it neutral, "C(=O)O", with no
    pH-assignment step anywhere in prep_ringaware.py or the Meeko call.
    Standard AutoDock Vina practice calls for assigning the dominant
    protonation state at a target pH before docking, since Vina scores one
    static structure and cannot represent an ionization equilibrium.
    A generic carboxylic acid pKa (~4.5-5) means these compounds should be
    >99% deprotonated at both physiological pH (7.4) and the ~pH 6.8 most
    yeast alpha-glucosidase inhibition assays in this project's own
    literature review use.

WHAT THIS SCRIPT DOES
    Deprotonates the C-28 carboxylic acid ("C(=O)O" -> "C(=O)[O-]") in each
    of the three current, fully-stereo-specified SMILES already used for the
    live site (betulinic/oleanolic unchanged since the original screen;
    ursolic acid uses its 2026-08-11 stereochemistry-corrected SMILES).
    Everything else held IDENTICAL to run_full_screen.py / run_stereo_fix.py:
    same receptor, same box, same exhaustiveness 32, same three seeds
    42/43/44, same ring-aware conformer selection.

    This is a CHECK, not an applied fix. Nothing here touches compounds.js
    or the live site. Results are for review before any decision on whether
    to actually change the pipeline's default protonation state.

    python run_protonation_check.py
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

INPUTS = os.path.join(HERE, "validation-inputs", "protonation-check")
OUT = os.path.join(HERE, "validation-run", "protonation-check")

CENTER = ["21.52", "-7.7", "23.55"]
BOX = ["30", "30", "30"]
EXHAUSTIVENESS = 32
SEEDS = (42, 43, 44)

ACARBOSE = -8.706  # unaffected by this check, no carboxylic acid on acarbose

# name, stem, pubchem_cid, formula (anion), conformers,
# SMILES = the live site's current neutral SMILES with C(=O)O -> C(=O)[O-]
# at the C-28 carboxylic acid only (verified: appears exactly once per
# molecule, the only other oxygens present are hydroxyls, not carbonyls).
COMPOUNDS = [
    ("Betulinic Acid", "betulinic_anion", 64971, "C30H47O3-", 150,
     "C=C(C)[C@@H]1CC[C@]2(C(=O)[O-])CC[C@]3(C)[C@H](CC[C@@H]4[C@@]5(C)CC"
     "[C@H](O)C(C)(C)[C@@H]5CC[C@]43C)[C@@H]12"),
    ("Oleanolic Acid", "oleanolic_anion", 10494, "C30H47O3-", 150,
     "CC1(C)CC[C@]2(C(=O)[O-])CC[C@]3(C)C(=CC[C@@H]4[C@@]5(C)CC[C@H](O)"
     "C(C)(C)[C@@H]5CC[C@]43C)[C@@H]2C1"),
    ("Ursolic Acid", "ursolic_anion", 64945, "C30H47O3-", 150,
     "C[C@@H]1CC[C@@]2(CC[C@@]3(C(=CC[C@H]4[C@]3(CC[C@@H]5[C@@]4(CC[C@@H]"
     "(C5(C)C)O)C)C)[C@@H]2[C@H]1C)C)C(=O)[O-]"),
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
              "preparation": "ring-aware, C-28 carboxylate anion (deprotonated)",
              "box_center": CENTER, "box_size": BOX,
              "exhaustiveness": EXHAUSTIVENESS, "seed": seed,
              "seconds": round(time.time() - started, 1),
              "returncode": done.returncode, "command": " ".join(command),
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


# Current live (protonated COOH) values, for the printed before/after only.
PRIOR = {"Betulinic Acid": -9.720, "Oleanolic Acid": -7.919, "Ursolic Acid": -8.651}
PRIOR_SPREAD = {"Betulinic Acid": 0.026, "Oleanolic Acid": 0.022, "Ursolic Acid": 0.038}


def main():
    os.makedirs(INPUTS, exist_ok=True)
    os.makedirs(OUT, exist_ok=True)

    print("")
    print("  STEP 1  ring-aware prep, C-28 carboxylate anion")
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
    print("  RESULT  protonated (live) / deprotonated (anion), vs Acarbose %.3f"
          % ACARBOSE)
    print("  " + "-" * 74)
    print("  %-16s %16s %16s   %s" % ("compound", "protonated (live)",
                                      "deprotonated", "verdict change"))
    for r in results:
        name = r["name"]

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
        print("  %-16s %16.3f %16.3f   %s" % (name, PRIOR[name], r["mean"],
                                              change))

    summary = {"check": "ligand protonation state (C-28 carboxylate)",
               "date": "2026-08-13",
               "protocol": {"receptor": RECEPTOR, "box_center": CENTER,
                            "box_size": BOX, "exhaustiveness": EXHAUSTIVENESS,
                            "seeds": list(SEEDS)},
               "reference_acarbose": ACARBOSE,
               "results": results,
               "prior_protonated_live": PRIOR,
               "prior_protonated_live_spread": PRIOR_SPREAD,
               "note": "This is a CHECK. compounds.js and the live site were "
                       "not modified by this run."}
    with open(os.path.join(OUT, "summary.json"), "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
    print("")
    print("  wrote %s" % os.path.join(OUT, "summary.json"))


if __name__ == "__main__":
    main()
