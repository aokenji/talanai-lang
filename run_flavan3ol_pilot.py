#!/usr/bin/env python3
"""
Flavan-3-ol pilot: dock catechin and epicatechin, the only flavan-3-ols in the
candidate library, on the identical protocol used for every other compound.

WHY THIS IS BEING RUN, flagged 2026-08-13 during defence-panel review
    The only measured phytochemistry that actually exists for Ziziphus talanai
    (Reyes et al. 2018, class-level screening) reports glycosides, condensed
    tannins and leucoanthocyanins. Condensed tannins are proanthocyanidin
    oligomers of flavan-3-ols; leucoanthocyanidins are flavan-3,4-diols. None
    of the ten compounds in the validated screen is a flavan-3-ol, a
    leucoanthocyanidin, or a proanthocyanidin - the one compound class the
    plant's own evidence actually points to is the class that was never
    docked. Catechin and Epicatechin are already curated in the candidate
    library (ids 14, 15) with status 'candidate', bindingAffinity: null.

WHAT THIS SCRIPT DOES
    Docks catechin and epicatechin only. SMILES pulled fresh from PubChem by
    CID (9064, 72276, the CIDs already recorded in compounds.js), verified
    with RDKit before use: both parse to C15H14O6, both have exactly 2 of 2
    stereocenters assigned, and their InChIKeys are genuinely distinct
    (PFTAWBLQPZVEMU-DZGCQCFKSA-N vs PFTAWBLQPZVEMU-UKRRQHHQSA-N) - this is
    the exact pair an earlier audit found byte-identical SMILES for in the
    SEPARATE HF Space backend (compounds.tsv); confirmed here that the site's
    own candidate-library CIDs are correct and genuinely distinct.

    Everything else held IDENTICAL to every other corrected run this project
    has made: same receptor, same box, same exhaustiveness 32, same three
    seeds 42/43/44, same ring-aware conformer selection (both compounds have
    one saturated pyran ring, the C-ring, so ring-awareness applies).

    THIS IS A PILOT, NOT AN APPLIED RESULT. Nothing here touches compounds.js
    or the live site. A real, documented caveat applies to any strong result:
    Vina's score divides its conformation-dependent term by
    (1 + 0.05846 * n_rotatable_bonds), so molecules are not compared on a
    torsion-neutral scale, and ligand efficiency (|score| / heavy atoms) does
    NOT remove this confound - it inverts it, favoring small molecules
    instead of favoring large ones (see REVIEW-FINDINGS.md Finding 11).
    Catechin and epicatechin are small (2 rotatable bonds, 22 heavy atoms)
    relative to the existing screen (3-16 rotatable bonds), so a strong raw
    score here is expected partly BECAUSE of Vina's own normalization, not
    purely because of genuine binding strength. Report the raw score, the
    rotatable-bond count, and this caveat together - never the raw score
    alone as a new headline.

    python run_flavan3ol_pilot.py
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

INPUTS = os.path.join(HERE, "validation-inputs", "flavan3ol-pilot")
OUT = os.path.join(HERE, "validation-run", "flavan3ol-pilot")

CENTER = ["21.52", "-7.7", "23.55"]
BOX = ["30", "30", "30"]
EXHAUSTIVENESS = 32
SEEDS = (42, 43, 44)

ACARBOSE = -8.706

# name, stem, pubchem_cid, formula, conformers, rotatable_bonds (for the LE
# caveat), fully-stereo-specified SMILES (fetched from PubChem 2026-08-13,
# verified with RDKit: formula, 2/2 stereocenters, distinct InChIKeys)
COMPOUNDS = [
    ("Catechin", "catechin", 9064, "C15H14O6", 150, 2,
     "C1[C@@H]([C@H](OC2=CC(=CC(=C21)O)O)C3=CC(=C(C=C3)O)O)O"),
    ("Epicatechin", "epicatechin", 72276, "C15H14O6", 150, 2,
     "C1[C@H]([C@H](OC2=CC(=CC(=C21)O)O)C3=CC(=C(C=C3)O)O)O"),
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
              "preparation": "ring-aware, fully stereo-specified SMILES from PubChem CID",
              "box_center": CENTER, "box_size": BOX,
              "exhaustiveness": EXHAUSTIVENESS, "seed": seed,
              "seconds": round(time.time() - started, 1),
              "returncode": done.returncode, "command": " ".join(command),
              "not_validated": ("the redocking gate does not pass for this "
                                "configuration; ranking only, no pose claims. "
                                "This is a pilot screen, not an applied result."),
              "pose_file": pose if os.path.isfile(pose) else None}
    if done.returncode != 0:
        record["error"] = (done.stderr or done.stdout)[-800:]
    else:
        record["rank1_score_kcal_mol"] = tal_control.parse_score(done.stdout)
    with open(os.path.join(OUT, tag + ".json"), "w", encoding="utf-8") as handle:
        json.dump(record, handle, indent=2, sort_keys=True)
    return record.get("rank1_score_kcal_mol")


def main():
    os.makedirs(INPUTS, exist_ok=True)
    os.makedirs(OUT, exist_ok=True)

    print("")
    print("  STEP 1  ring-aware prep, flavan-3-ol pilot")
    print("  " + "-" * 74)
    ready = []
    for name, stem, cid, formula, confs, rotb, smiles in COMPOUNDS:
        result = ringaware.prepare_ring_aware(smiles, confs)
        thetas = result["thetas"]
        ok = ringaware.all_chairs(thetas) if thetas else True
        path = os.path.join(INPUTS, stem + ".pdbqt")
        with open(path, "w", encoding="utf-8") as handle:
            handle.write(result["pdbqt"])
        note = ("no saturated ring" if not thetas
                else ("all chairs" if ok else "*** STILL DEFECTIVE ***"))
        print("  %-14s %3d/%3d all-chair, %s" % (name, result["n_qualifying"],
                                                 result["n_conformers"], note))
        ready.append((name, stem, cid, formula, rotb, path))

    print("")
    print("  STEP 2  docking %d ligands x %d seeds at exhaustiveness %d"
          % (len(ready), len(SEEDS), EXHAUSTIVENESS))
    print("  " + "-" * 74)

    results = []
    for name, stem, cid, formula, rotb, path in ready:
        scores = []
        for seed in SEEDS:
            score = dock(path, seed, "%s_seed%d" % (stem, seed))
            if score is None:
                print("  %-14s seed %d FAILED" % (name, seed))
            else:
                scores.append(score)
                print("  %-14s seed %d  %8.3f" % (name, seed, score))
        if scores:
            mean = sum(scores) / len(scores)
            spread = max(scores) - min(scores)
            results.append({"name": name, "mean": round(mean, 3),
                            "spread": round(spread, 3), "rotatable_bonds": rotb})

    print("")
    print("  RESULT  vs Acarbose %.3f (rotatable-bond count shown - see LE caveat)"
          % ACARBOSE)
    print("  " + "-" * 74)
    print("  %-14s %10s %10s %6s   %s" % ("compound", "mean", "spread", "rotb", "verdict"))
    for r in results:
        h = r["spread"] / 2.0
        if r["mean"] + h < ACARBOSE - 0.005:
            verdict = "beats"
        elif r["mean"] - h > ACARBOSE + 0.005:
            verdict = "loses to"
        else:
            verdict = "indistinguishable from"
        print("  %-14s %10.3f %10.3f %6d   %s" % (r["name"], r["mean"], r["spread"],
                                                    r["rotatable_bonds"], verdict))

    summary = {"pilot": "flavan-3-ol screen (catechin, epicatechin)",
               "date": "2026-08-13",
               "status": "PILOT - not applied to compounds.js or the live site",
               "protocol": {"receptor": RECEPTOR, "box_center": CENTER,
                            "box_size": BOX, "exhaustiveness": EXHAUSTIVENESS,
                            "seeds": list(SEEDS)},
               "reference_acarbose": ACARBOSE,
               "results": results,
               "ligand_efficiency_caveat": (
                   "Vina's score divides its conformation-dependent term by "
                   "(1 + 0.05846 * rotatable_bonds). Catechin and epicatechin "
                   "have 2 rotatable bonds each, far fewer than most of the "
                   "existing 10-compound screen (3-16). A strong raw score here "
                   "is expected partly because of this normalization favoring "
                   "small, rigid molecules, not purely because of binding "
                   "strength. Ligand efficiency (|score|/heavy atoms) does not "
                   "remove this confound, it inverts it. Report rotatable-bond "
                   "count and this caveat alongside any score from this run.")}
    with open(os.path.join(OUT, "summary.json"), "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
    print("")
    print("  wrote %s" % os.path.join(OUT, "summary.json"))


if __name__ == "__main__":
    main()
