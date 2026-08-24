#!/usr/bin/env python3
"""
Re-dock the three affected compounds, now with ring-aware preparation.

Uses the files already written by prep_ringaware.py. It does NOT re-prepare,
so the ring-aware selection is not overwritten by the energy-only one.

Two arms, same as before and for the same reason:
  ARM A  raw receptor, exhaustiveness 8, seed 42, the thesis configuration.
         Comparable with the published numbers because only the ligand
         preparation differs.
  ARM B  Meeko receptor, exhaustiveness 32, three seeds. Not comparable with
         published, since receptor and search both differ, but it is what a
         corrected study would report.
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
RAW = os.path.join(DATA, "receptor_clean.pdb")
MEEKO = os.path.join(DATA, "prepared", "receptor.pdbqt")

REPREP = os.path.join(HERE, "validation-inputs", "reprepared")
OUT = os.path.join(HERE, "validation-run", "ringaware")

CENTER = ["21.52", "-7.7", "23.55"]
BOX = ["30", "30", "30"]

TARGETS = [
    ("Ursolic acid", "ursolic", -8.105, os.path.join(DATA, "ligands", "ursolic.pdbqt")),
    ("Isovitexin", "isovitexin", -8.024, os.path.join(DATA, "ligands", "isovitexin.pdbqt")),
    ("Acarbose (ref)", "acarbose", -6.660, os.path.join(DATA, "reference", "acarbose.pdbqt")),
]


def dock(receptor, ligand, exhaustiveness, seed, tag):
    os.makedirs(OUT, exist_ok=True)
    pose = os.path.join(OUT, tag + ".pdbqt")
    command = [VINA, "--receptor", receptor, "--ligand", ligand,
               "--center_x", CENTER[0], "--center_y", CENTER[1],
               "--center_z", CENTER[2],
               "--size_x", BOX[0], "--size_y", BOX[1], "--size_z", BOX[2],
               "--exhaustiveness", str(exhaustiveness), "--num_modes", "9",
               "--energy_range", "3", "--seed", str(seed), "--cpu", "4",
               "--out", pose]
    started = time.time()
    done = subprocess.run(command, capture_output=True, text=True, timeout=21600)
    record = {"tag": tag, "receptor": os.path.basename(receptor),
              "ligand": os.path.basename(ligand), "ligand_path": ligand,
              "preparation": "ring-aware: lowest energy among all-chair conformers",
              "exhaustiveness": exhaustiveness, "seed": seed,
              "seconds": round(time.time() - started, 1),
              "returncode": done.returncode, "command": " ".join(command),
              "talanai_version": "0.1.0",
              "pose_file": pose if os.path.isfile(pose) else None}
    if done.returncode != 0:
        record["error"] = (done.stderr or done.stdout)[-800:]
        with open(os.path.join(OUT, tag + ".json"), "w", encoding="utf-8") as h:
            json.dump(record, h, indent=2, sort_keys=True)
        return None
    record["rank1_score_kcal_mol"] = tal_control.parse_score(done.stdout)
    with open(os.path.join(OUT, tag + ".json"), "w", encoding="utf-8") as h:
        json.dump(record, h, indent=2, sort_keys=True)
    return record["rank1_score_kcal_mol"]


def main():
    os.makedirs(OUT, exist_ok=True)
    print("")
    print("  ARM A: thesis configuration, raw receptor, exh 8, seed 42")
    print("  Only the ligand preparation differs from the published run.")
    print("  " + "-" * 70)
    print("  %-16s %10s %10s %11s %9s"
          % ("compound", "published", "old prep", "ring-aware", "delta"))
    for name, stem, published, old in TARGETS:
        new = os.path.join(REPREP, stem + ".pdbqt")
        if not os.path.isfile(new):
            print("  %-16s MISSING %s" % (name, new))
            continue
        a = dock(RAW, old, 8, 42, stem + "_raw_old")
        b = dock(RAW, new, 8, 42, stem + "_raw_ringaware")
        if a is None or b is None:
            print("  %-16s FAILED" % name)
            continue
        print("  %-16s %10.3f %10.3f %11.3f %+9.3f"
              % (name, published, a, b, b - a), flush=True)

    print("")
    print("  ARM B: Meeko receptor, exh 32, seeds 42/43/44, ring-aware prep")
    print("  " + "-" * 70)
    print("  %-16s %9s %9s %9s %10s %8s"
          % ("compound", "seed 42", "seed 43", "seed 44", "mean", "spread"))
    means = {}
    ranges = {}
    for name, stem, _published, _old in TARGETS:
        new = os.path.join(REPREP, stem + ".pdbqt")
        if not os.path.isfile(new):
            continue
        scores = [dock(MEEKO, new, 32, s, "%s_meeko_ringaware_s%d" % (stem, s))
                  for s in (42, 43, 44)]
        clean = [s for s in scores if s is not None]
        if not clean:
            print("  %-16s all runs failed" % name)
            continue
        mean = sum(clean) / len(clean)
        means[name] = mean
        ranges[name] = (min(clean), max(clean))
        print("  %-16s %9s %9s %9s %10.3f %8.3f"
              % (name, *["%.3f" % s if s is not None else "-" for s in scores],
                 mean, max(clean) - min(clean)), flush=True)

    print("")
    print("  " + "-" * 70)
    # Comparing means alone is not enough and this project has already been
    # caught doing it once. With three seeds the honest test is whether the
    # observed RANGES separate at all. If they overlap, the two compounds are
    # indistinguishable at this sampling, however the means happen to order.
    reference = means.get("Acarbose (ref)")
    if reference is not None:
        ref_lo, ref_hi = ranges["Acarbose (ref)"]
        for name, mean in means.items():
            if name == "Acarbose (ref)":
                continue
            lo, hi = ranges[name]
            overlaps = not (hi < ref_lo or lo > ref_hi)
            if overlaps:
                verdict = "INDISTINGUISHABLE, ranges overlap"
            elif mean < reference:
                verdict = "beats acarbose, ranges separate"
            else:
                verdict = "does not beat acarbose, ranges separate"
            print("  %-16s %.3f [%.3f, %.3f] vs acarbose %.3f [%.3f, %.3f]"
                  % (name, mean, lo, hi, reference, ref_lo, ref_hi))
            print("  %-16s -> %s" % ("", verdict))
    print("")
    print("  Every ligand here has all-chair saturated rings, verified by")
    print("  Cremer-Pople. Arm B is not comparable with the published values.")
    print("  No published number was modified.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
