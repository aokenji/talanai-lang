#!/usr/bin/env python3
"""
Independent replication of the five example compounds plus the reference.

WHAT THIS IS
    A re-run of the CURRENT corrected protocol, from SMILES, with a fresh
    ligand preparation, to answer one question: does the published number come
    back when the pipeline is run again?

    Same receptor, same box, same exhaustiveness 32, same ligand preparation.
    The variable is the SEARCH: seeds 45/46/47 instead of the published
    42/43/44.

    This is deliberate. prep_ringaware embeds with a fixed seed (0xF00D) and
    Vina is deterministic given identical inputs and seed, so a same-seed re-run
    would reproduce the published file byte for byte and prove nothing. The
    language says so itself: RA01 records that same-seed replication "does not
    measure search convergence". Different seeds ask the question worth asking,
    which is whether the ranking survives starting the search somewhere else.

WHY IT EXISTS
    Until now there was no replication of the FINAL protocol. The published
    values are a composite: run_full_screen.py (2026-08-04) for compounds the
    later corrections did not touch, and run_protonation_check.py (2026-08-13)
    for the C-28 acids. Both are primary runs. Neither checks the other.

    examples/alpha-glucosidase.tal previously carried a replication block whose
    numbers replicated the PRE-correction screen. It was removed rather than
    updated, because inventing replacement numbers would be fabrication. This
    script produces the real ones.

WHAT IT IS NOT
    A validated protocol. The redocking gate does not pass: this configuration
    cannot recover a crystal pose from an independently built ligand. See
    VALIDATION-FINDINGS.md section 4c. These are RANKINGS. No pose-level or
    interaction claim may be drawn from them.

SMILES PROVENANCE, per compound, deliberately not uniform
    Rutin, Quercetin, Kaempferol, Acarbose
        neutral SMILES exactly as in run_full_screen.py. The 2026-08-11 and
        2026-08-13 corrections did not touch these.
    Betulinic Acid, Oleanolic Acid
        C-28 carboxylate ANION SMILES exactly as in run_protonation_check.py,
        which is the form the published values were measured on.

    Mixing those two sources is correct and is the whole subtlety here: the
    published dataset is itself a composite, so a replication has to reproduce
    the composite, not one of its halves.
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

INPUTS = os.path.join(HERE, "validation-inputs", "replication")
OUT = os.path.join(HERE, "validation-run", "replication")

CENTER = ["21.52", "-7.7", "23.55"]
BOX = ["30", "30", "30"]
EXHAUSTIVENESS = 32
SEEDS = (45, 46, 47)   # NEW seeds: this is a convergence test, not a replay

# name, stem, published (src/data/compounds.js), formula, conformers, SMILES
COMPOUNDS = [
    ("Rutin", "rutin", -10.317, "C27H30O16", 200,
     "C[C@@H]1O[C@@H](OC[C@H]2O[C@@H](Oc3c(-c4ccc(O)c(O)c4)oc4cc(O)cc(O)c4c3=O)"
     "[C@H](O)[C@@H](O)[C@@H]2O)[C@H](O)[C@H](O)[C@H]1O"),
    ("Betulinic Acid", "betulinic", -9.084, "C30H47O3-", 150,
     "C=C(C)[C@@H]1CC[C@]2(C(=O)[O-])CC[C@]3(C)[C@H](CC[C@@H]4[C@@]5(C)CC"
     "[C@H](O)C(C)(C)[C@@H]5CC[C@]43C)[C@@H]12"),
    ("Quercetin", "quercetin", -8.818, "C15H10O7", 100,
     "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12"),
    ("Kaempferol", "kaempferol", -8.510, "C15H10O6", 100,
     "O=c1c(O)c(-c2ccc(O)cc2)oc2cc(O)cc(O)c12"),
    ("Oleanolic Acid", "oleanolic", -7.593, "C30H47O3-", 150,
     "CC1(C)CC[C@]2(C(=O)[O-])CC[C@]3(C)C(=CC[C@@H]4[C@@]5(C)CC[C@H](O)"
     "C(C)(C)[C@@H]5CC[C@]43C)[C@@H]2C1"),
    ("Acarbose", "acarbose", -8.706, "C25H43NO18", 300,
     "C[C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@@H]2[C@H](O[C@@H]([C@@H]([C@H]2O)"
     "O)O[C@H]([C@@H](CO)O)[C@@H]([C@H](C=O)O)O)CO)O)O)N[C@H]3C=C([C@H]([C@@H]"
     "([C@H]3O)O)O)CO"),
]

# The project's own convention for run-to-run agreement. Not a published Vina
# precision figure, and deliberately not presented as one.
REPRODUCIBILITY_LIMIT = 0.10


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
              "preparation": "ring-aware, identical to published; only the search seed differs",
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


def main():
    os.makedirs(INPUTS, exist_ok=True)
    os.makedirs(OUT, exist_ok=True)

    print("")
    print("  Replication of the corrected protocol")
    print("  exhaustiveness %d, seeds %s (convergence test vs published 42/43/44)"
          % (EXHAUSTIVENESS, "/".join(str(s) for s in SEEDS)))
    print("  " + "-" * 68)

    ready = []
    for name, stem, published, formula, confs, smiles in COMPOUNDS:
        result = ringaware.prepare_ring_aware(smiles, confs)
        thetas = result.get("thetas") if isinstance(result, dict) else None
        chairs = ringaware.all_chairs(thetas) if thetas else True
        path = os.path.join(INPUTS, stem + ".pdbqt")
        block = result.get("pdbqt") if isinstance(result, dict) else result
        with open(path, "w", encoding="utf-8") as handle:
            handle.write(block)
        ready.append((name, stem, published, path))
        print("    prepared %-16s %s" % (name, "chairs ok" if chairs else "RING NOT CHAIR"))

    print("")
    rows = []
    for name, stem, published, path in ready:
        scores = []
        for seed in SEEDS:
            score = dock(path, seed, "%s_s%d" % (stem, seed))
            if score is None:
                print("    %-16s seed %d FAILED" % (name, seed))
                continue
            scores.append(score)
            print("    %-16s seed %d  %8.3f" % (name, seed, score))
        if not scores:
            continue
        mean = sum(scores) / len(scores)
        spread = max(scores) - min(scores)
        delta = mean - published
        rows.append({"name": name, "scores": scores, "mean": round(mean, 3),
                     "spread": round(spread, 3), "published": published,
                     "delta": round(delta, 3),
                     "within_reproducibility": abs(delta) <= REPRODUCIBILITY_LIMIT})
        print("    %-16s mean %8.3f  published %8.3f  delta %+.3f  %s"
              % (name, mean, published, delta,
                 "ok" if abs(delta) <= REPRODUCIBILITY_LIMIT else "OUTSIDE CONVENTION"))
        print("")

    summary = {
        "what": "independent replication of the corrected protocol",
        "date": time.strftime("%Y-%m-%d"),
        "pipeline": "run_replication.py, vina 1.2.7, 3A4A prepared receptor",
        "kind": "different_seed",
        "varied": ("search seeds only: 45/46/47 against the published 42/43/44. "
                   "Receptor, box, exhaustiveness and ligand preparation are "
                   "identical, so any deviation is search convergence, not drift."),
        "reproducibility_limit_kcal_mol": REPRODUCIBILITY_LIMIT,
        "reproducibility_limit_note": ("this project's convention for run-to-run "
                                       "agreement, not a published Vina precision figure"),
        "protocol": {"box_center": CENTER, "box_size": BOX,
                     "exhaustiveness": EXHAUSTIVENESS, "seeds": list(SEEDS),
                     "receptor": RECEPTOR},
        "not_validated": ("the redocking gate does not pass; these are rankings, "
                          "not placements, and no interaction claim follows"),
        "results": rows,
    }
    with open(os.path.join(OUT, "summary.json"), "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)

    worst = max(rows, key=lambda r: abs(r["delta"])) if rows else None
    print("  " + "-" * 68)
    if worst:
        print("  worst deviation: %s %+.3f kcal/mol" % (worst["name"], worst["delta"]))
        print("  all within %.2f convention: %s"
              % (REPRODUCIBILITY_LIMIT, all(r["within_reproducibility"] for r in rows)))
    print("  Written to %s" % OUT)
    print("")


if __name__ == "__main__":
    main()
