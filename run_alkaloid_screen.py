#!/usr/bin/env python3
"""
The cyclopeptide alkaloid arm: five genus-Ziziphus macrocyclic alkaloids.

WHY THIS EXISTS
    The screened set is seven flavonoids and three triterpenes. That is exactly
    the compound profile the network-pharmacology literature is criticised for:
    quercetin, kaempferol, luteolin and rutin are in almost every plant, so a
    reviewer reads the set as generic.

    Cyclopeptide alkaloids are different. They are a chemotaxonomic signature of
    Ziziphus and Rhamnaceae rather than a ubiquitous class, so the genus-level
    inclusion argument is materially STRONGER for them than for the flavonoids.
    They are also the class implicated in the only Ziziphus study that measured
    a cholinesterase endpoint and identified compounds (Foyet et al. 2019).

PROTOCOL
    Identical to run_full_screen.py in every respect: same receptor with Meeko
    hydrogen-bond atom typing, same 30 A box on the same centre, exhaustiveness
    32, seeds 42/43/44, ring-aware preparation. Nothing is changed, so these
    scores are directly comparable to the corrected screen.

WHAT THIS IS NOT
    A validated protocol. The redocking gate does not pass; see
    VALIDATION-FINDINGS.md section 4c. Ranking only. No pose or interaction
    claim may be drawn from this run.

NOTE ON MACROCYCLES
    These ligands are 13- to 14-membered macrocycles. ETKDGv3 handles them with
    macrocycle torsion preferences, but embedding is harder than for the
    flavonoids and conformer yield may be low. The saturated-six-ring chair test
    that R307 enforces mostly does not apply here, because the macrocycle is not
    a six-membered saturated ring. Where a ligand has no saturated six-ring the
    check correctly reports "no saturated ring" and passes.

    SMILES were fetched from PubChem by CID, not typed from memory.

    python run_alkaloid_screen.py
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

INPUTS = os.path.join(HERE, "validation-inputs", "alkaloids")
OUT = os.path.join(HERE, "validation-run", "alkaloids")

CENTER = ["21.52", "-7.7", "23.55"]
BOX = ["30", "30", "30"]
EXHAUSTIVENESS = 32
SEEDS = (42, 43, 44)

# name, stem, PubChem CID, formula, conformers, SMILES (from PubChem PUG REST)
COMPOUNDS = [
    ("Sanjoinine A", "sanjoinine_a", 11731186, "C31H42N4O4", 200,
     "CC(C)C[C@H]1C(=O)N/C=C\\C2=CC=C(C=C2)O[C@H]([C@@H](C(=O)N1)NC(=O)"
     "[C@H](CC3=CC=CC=C3)N(C)C)C(C)C"),
    ("Mauritine A", "mauritine_a", 11353668, "C32H41N5O5", 200,
     "C[C@@H](C(=O)N[C@@H](C(C)C)C(=O)N1CC[C@H]2[C@H]1C(=O)N[C@H](C(=O)N/C=C\\"
     "C3=CC=C(O2)C=C3)CC4=CC=CC=C4)N(C)C"),
    ("Frangulanine", "frangulanine", 163078670, "C28H44N4O4", 200,
     "CC[C@H](C)[C@@H](C(=O)N[C@H]1[C@H](OC2=CC=C(C=C2)C=CNC(=O)[C@@H](NC1=O)"
     "CC(C)C)C(C)C)N(C)C"),
    ("Zizyphine A", "zizyphine_a", 6324833, "C33H49N5O6", 250,
     "CC[C@H](C)[C@@H](C(=O)N1CC[C@H]2[C@H]1C(=O)N3CCC[C@H]3C(=O)N/C=C\\C4="
     "C(C=C(O2)C=C4)OC)NC(=O)[C@H]([C@@H](C)CC)N(C)C"),
    ("Amphibine D", "amphibine_d", 5318120, "C36H49N5O5", 250,
     "CCC(C)C1C(=O)N/C=C\\C2=CC=C(C=C2)OC3CCN(C3C(=O)N1)C(=O)C(C(C)CC)NC(=O)"
     "C(CC4=CC=CC=C4)N(C)C"),
]

# The corrected-screen reference, for context in the printed table.
ACARBOSE = -8.706


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
              "preparation": "ring-aware: lowest energy among all-chair conformers",
              "box_center": CENTER, "box_size": BOX,
              "exhaustiveness": EXHAUSTIVENESS, "seed": seed,
              "seconds": round(time.time() - started, 1),
              "returncode": done.returncode, "command": " ".join(command),
              "talanai_version": "0.1.0",
              "arm": "cyclopeptide alkaloids",
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


def main():
    os.makedirs(INPUTS, exist_ok=True)
    os.makedirs(OUT, exist_ok=True)

    print("")
    print("  STEP 1  ring-aware preparation of %d macrocyclic alkaloids"
          % len(COMPOUNDS))
    print("  " + "-" * 74)
    ready, skipped = [], []
    for name, stem, cid, formula, confs, smiles in COMPOUNDS:
        started = time.time()
        try:
            result = ringaware.prepare_ring_aware(smiles, confs)
        except Exception as exc:
            print("  %-16s PREPARATION FAILED: %s" % (name, str(exc)[:50]))
            skipped.append((name, "preparation failed: %s" % str(exc)[:80]))
            continue
        thetas = result["thetas"]
        ok = ringaware.all_chairs(thetas) if thetas else True
        path = os.path.join(INPUTS, stem + ".pdbqt")
        with open(path, "w", encoding="utf-8") as handle:
            handle.write(result["pdbqt"])
        note = ("no saturated six-ring" if not thetas
                else ("all chairs" if ok else "*** STILL DEFECTIVE ***"))
        print("  %-16s %3d/%3d qualifying, %-22s %5.1fs"
              % (name, result["n_qualifying"], result["n_conformers"], note,
                 time.time() - started))
        if ok:
            ready.append((name, stem, cid, formula, path))
        else:
            skipped.append((name, "rings still defective"))

    if skipped:
        print("")
        print("  SKIPPED:")
        for name, why in skipped:
            print("    %-16s %s" % (name, why))

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
            results.append({"name": name, "stem": stem, "pubchem_cid": cid,
                            "formula": formula, "seeds": list(SEEDS),
                            "scores": scores, "mean": round(mean, 3),
                            "spread": round(spread, 3)})

    results.sort(key=lambda r: r["mean"])

    print("")
    print("  RESULT  ranked, against the corrected-screen Acarbose of %.3f"
          % ACARBOSE)
    print("  " + "-" * 74)
    print("  %-16s %9s %8s   %s" % ("compound", "mean", "spread", "vs Acarbose"))
    for r in results:
        # Same separation rule as the site: ranges must not overlap.
        half = r["spread"] / 2.0
        if r["mean"] + half < ACARBOSE - 0.005:
            verdict = "beats"
        elif r["mean"] - half > ACARBOSE + 0.005:
            verdict = "loses to"
        else:
            verdict = "indistinguishable from"
        print("  %-16s %9.3f %8.3f   %s" % (r["name"], r["mean"], r["spread"],
                                            verdict))

    summary = {"arm": "cyclopeptide alkaloids",
               "protocol": {"receptor": RECEPTOR, "box_center": CENTER,
                            "box_size": BOX, "exhaustiveness": EXHAUSTIVENESS,
                            "seeds": list(SEEDS),
                            "preparation": "ring-aware"},
               "reference_acarbose": ACARBOSE,
               "not_validated": ("redocking gate does not pass; ranking only, "
                                 "no pose or interaction claims"),
               "results": results, "skipped": skipped}
    with open(os.path.join(OUT, "summary.json"), "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
    print("")
    print("  wrote %s" % os.path.join(OUT, "summary.json"))


if __name__ == "__main__":
    main()
