#!/usr/bin/env python3
"""
The full corrected screen: all ten screened compounds plus the reference.

WHAT THIS IS
    One internally consistent dataset in which
      - the receptor carries AutoDock hydrogen-bond atom types (Meeko), rather
        than the raw PDB whose H-bond term was measured at exactly zero
      - every ligand's saturated rings are chairs, verified by Cremer-Pople,
        selected ring-aware rather than by energy alone
      - sampling is exhaustiveness 32 across three seeds, which was shown to
        converge to spreads of 0.01 to 0.04 kcal/mol

WHAT THIS IS NOT
    A validated protocol. The redocking gate does not pass: this configuration
    cannot recover a crystal pose from an independently built ligand. See
    VALIDATION-FINDINGS.md section 4c. These numbers are a RANKING, and the
    study's own conclusion is that this protocol ranks but does not place.

    Pose-level and interaction claims are not supported by this run and must
    not be derived from it.

BOX
    Kept at 30 A, matching the published screen, so the ranking stays
    comparable. The 18 A focused box helps pose recovery but pose claims are
    out of scope here, and changing two things at once is the error this whole
    investigation exists to document.

    python run_full_screen.py
"""

import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from talanai import chem, control as tal_control   # noqa: E402
import prep_ringaware as ringaware                 # noqa: E402

ASSETS = (r"D:\THESIS_VSC\BALAKATDBV2\talanaidock\dist\TalanaiDock\app"
          r"\docking_assets")
DATA = os.path.join(ASSETS, "docking_data")
VINA = os.path.join(ASSETS, "vina.exe")
RECEPTOR = os.path.join(DATA, "prepared", "receptor.pdbqt")

INPUTS = os.path.join(HERE, "validation-inputs", "screen")
OUT = os.path.join(HERE, "validation-run", "full-screen")

CENTER = ["21.52", "-7.7", "23.55"]
BOX = ["30", "30", "30"]
EXHAUSTIVENESS = 32
SEEDS = (42, 43, 44)

# name, stem, published affinity, formula, SMILES, conformers to generate
COMPOUNDS = [
    ("Rutin", "rutin", -8.857, "C27H30O16", 200,
     "C[C@@H]1O[C@@H](OC[C@H]2O[C@@H](Oc3c(-c4ccc(O)c(O)c4)oc4cc(O)cc(O)c4c3=O)"
     "[C@H](O)[C@@H](O)[C@@H]2O)[C@H](O)[C@H](O)[C@H]1O"),
    ("Betulinic Acid", "betulinic", -8.290, "C30H48O3", 150,
     "C=C(C)[C@@H]1CC[C@]2(C(=O)O)CC[C@]3(C)[C@H](CC[C@@H]4[C@@]5(C)CC[C@H](O)"
     "C(C)(C)[C@@H]5CC[C@]43C)[C@@H]12"),
    ("Ursolic Acid", "ursolic", -8.105, "C30H48O3", 150,
     "CC1CCC2(C(=O)O)CCC3(C)C(=CCC4C5(C)CCC(O)C(C)(C)C5CCC43C)C2C1C"),
    ("Isovitexin", "isovitexin", -8.024, "C21H20O10", 150,
     "O=c1cc(-c2ccc(O)cc2)oc2cc(O)c(C3OC(CO)C(O)C(O)C3O)c(O)c12"),
    ("Spinosin", "spinosin", -7.859, "C28H32O15", 300,
     "COC1=C(C(=C2C(=C1)OC(=CC2=O)C3=CC=C(C=C3)O)O)[C@H]4[C@@H]([C@H]([C@@H]"
     "([C@H](O4)CO)O)O)O[C@H]5[C@@H]([C@H]([C@@H]([C@H](O5)CO)O)O)O"),
    ("Luteolin", "luteolin", -7.733, "C15H10O6", 100,
     "O=c1cc(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12"),
    ("Quercetin", "quercetin", -7.503, "C15H10O7", 100,
     "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12"),
    ("Kaempferol", "kaempferol", -7.479, "C15H10O6", 100,
     "O=c1c(O)c(-c2ccc(O)cc2)oc2cc(O)cc(O)c12"),
    ("Vitexin", "vitexin", -7.469, "C21H20O10", 150,
     "O=c1cc(-c2ccc(O)cc2)oc2c(C3OC(CO)C(O)C(O)C3O)c(O)cc(O)c12"),
    ("Oleanolic Acid", "oleanolic", -6.922, "C30H48O3", 150,
     "CC1(C)CC[C@]2(C(=O)O)CC[C@]3(C)C(=CC[C@@H]4[C@@]5(C)CC[C@H](O)C(C)(C)"
     "[C@@H]5CC[C@]43C)[C@@H]2C1"),
    ("Acarbose", "acarbose", -6.660, "C25H43NO18", 300,
     "C[C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@@H]2[C@H](O[C@@H]([C@@H]([C@H]2O)"
     "O)O[C@H]([C@@H](CO)O)[C@@H]([C@H](C=O)O)O)CO)O)O)N[C@H]3C=C([C@H]([C@@H]"
     "([C@H]3O)O)O)CO"),
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
              "preparation": "ring-aware: lowest energy among all-chair conformers",
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


def main():
    os.makedirs(INPUTS, exist_ok=True)
    os.makedirs(OUT, exist_ok=True)

    print("")
    print("  STEP 1  ring-aware preparation of all %d ligands" % len(COMPOUNDS))
    print("  " + "-" * 74)
    ready, skipped = [], []
    for name, stem, published, formula, confs, smiles in COMPOUNDS:
        try:
            result = ringaware.prepare_ring_aware(smiles, confs)
        except Exception as exc:
            print("  %-16s PREPARATION FAILED: %s" % (name, exc))
            skipped.append((name, "preparation failed"))
            continue
        thetas = result["thetas"]
        ok = ringaware.all_chairs(thetas) if thetas else True
        path = os.path.join(INPUTS, stem + ".pdbqt")
        with open(path, "w", encoding="utf-8") as handle:
            handle.write(result["pdbqt"])
        note = ("no saturated ring" if not thetas
                else ("all chairs" if ok else "*** STILL DEFECTIVE ***"))
        print("  %-16s %3d/%3d all-chair, %s" % (name, result["n_qualifying"],
                                                 result["n_conformers"], note))
        if ok:
            ready.append((name, stem, published, formula, path))
        else:
            skipped.append((name, "rings still defective"))

    if skipped:
        print("")
        print("  NOT DOCKED, because a bad ring cannot be fixed by docking:")
        for name, why in skipped:
            print("    %-16s %s" % (name, why))

    print("")
    print("  STEP 2  docking, Meeko receptor, %s A box, exh %d, seeds %s"
          % (BOX[0], EXHAUSTIVENESS, ", ".join(str(s) for s in SEEDS)))
    print("  " + "-" * 74)
    print("  %-16s %8s %8s %8s %9s %7s %9s"
          % ("compound", "seed42", "seed43", "seed44", "mean", "spread", "published"))

    rows = []
    for name, stem, published, formula, path in ready:
        scores = [dock(path, seed, "%s_s%d" % (stem, seed)) for seed in SEEDS]
        clean = [s for s in scores if s is not None]
        if not clean:
            print("  %-16s all runs failed" % name)
            continue
        mean = sum(clean) / len(clean)
        rows.append({"name": name, "published": published, "formula": formula,
                     "scores": scores, "mean": mean,
                     "lo": min(clean), "hi": max(clean),
                     "spread": max(clean) - min(clean),
                     "heavy": chem.heavy_atoms(formula)})
        print("  %-16s %8s %8s %8s %9.3f %7.3f %9.3f"
              % (name, *["%.3f" % s if s is not None else "-" for s in scores],
                 mean, max(clean) - min(clean), published), flush=True)

    if not rows:
        print("  nothing docked")
        return 1

    reference = next((r for r in rows if r["name"] == "Acarbose"), None)

    print("")
    print("  STEP 3  ranking, with size correction")
    print("  " + "-" * 74)
    print("  %-16s %9s %6s %9s %s"
          % ("compound", "mean", "heavy", "per atom", "vs acarbose"))
    for row in sorted(rows, key=lambda r: r["mean"]):
        efficiency = abs(row["mean"]) / row["heavy"] if row["heavy"] else 0
        verdict = ""
        if reference and row["name"] != "Acarbose":
            # Ranges must separate before a winner is named. Means alone are
            # not enough; this project has been caught on that twice.
            if row["hi"] < reference["lo"]:
                verdict = "beats"
            elif row["lo"] > reference["hi"]:
                verdict = "loses to"
            else:
                verdict = "INDISTINGUISHABLE"
        print("  %-16s %9.3f %6d %9.3f  %s"
              % (row["name"], row["mean"], row["heavy"], efficiency, verdict))

    summary = {"generated": "2026-08-04", "receptor": RECEPTOR,
               "box_center": CENTER, "box_size": BOX,
               "exhaustiveness": EXHAUSTIVENESS, "seeds": list(SEEDS),
               "preparation": "ring-aware, all saturated rings verified chairs",
               "validated": False,
               "validation_note": ("the redocking gate does not pass for this "
                                   "configuration. Ranking only. No pose or "
                                   "interaction claims may be derived from this "
                                   "run. See VALIDATION-FINDINGS.md section 4c."),
               "results": rows}
    with open(os.path.join(OUT, "summary.json"), "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)

    print("")
    print("  " + "-" * 74)
    print("  Written to %s" % OUT)
    print("  NOT a validated protocol. Ranking only. No pose claims.")
    print("  Published values are unchanged and were never written to.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
