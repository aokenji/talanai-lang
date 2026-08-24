#!/usr/bin/env python3
"""
Is the ranking a property of chemistry, or a property of AutoDock Vina?

THE QUESTION
    The corrected screen ranks eleven compounds using Vina's scoring function.
    If a different, independently parameterised function ranks them the same
    way, the ordering is more likely to reflect something real. If it does not,
    the ordering is a property of one scoring function and should be presented
    that way.

    This is the check a reviewer means by "is your ranking robust?", and it is
    cheap because the poses already exist.

TWO TESTS, WHICH ANSWER DIFFERENT THINGS
    RESCORE   score the SAME Vina poses with Vinardo. Isolates the scoring
              function: same geometry, different function. Seconds per ligand.
    REDOCK    dock again from scratch with Vinardo driving the search. Fully
              independent ranking, but slow, so seed 42 only.

    Vinardo (Quiroga and Villarreal, PLoS ONE 2016) is a reparameterised
    scoring function shipped with Vina 1.2.x, so this needs no extra software.

    python run_consensus.py            rescore only
    python run_consensus.py --redock   rescore, then re-dock with Vinardo
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
VINA = os.path.join(ASSETS, "vina.exe")
RECEPTOR = os.path.join(ASSETS, "docking_data", "prepared", "receptor.pdbqt")

SCREEN = os.path.join(HERE, "validation-run", "full-screen")
INPUTS = os.path.join(HERE, "validation-inputs", "screen")
OUT = os.path.join(HERE, "validation-run", "consensus")

CENTER = ["21.52", "-7.7", "23.55"]
BOX = ["30", "30", "30"]

STEMS = [
    ("Rutin", "rutin"), ("Betulinic Acid", "betulinic"),
    ("Ursolic Acid", "ursolic"), ("Isovitexin", "isovitexin"),
    ("Spinosin", "spinosin"), ("Luteolin", "luteolin"),
    ("Quercetin", "quercetin"), ("Kaempferol", "kaempferol"),
    ("Vitexin", "vitexin"), ("Oleanolic Acid", "oleanolic"),
    ("Acarbose", "acarbose"),
]


def spearman(a, b):
    """Rank correlation, standard library only. Ties get average ranks."""
    def ranks(values):
        order = sorted(range(len(values)), key=lambda i: values[i])
        out = [0.0] * len(values)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
                j += 1
            average = (i + j) / 2.0 + 1
            for k in range(i, j + 1):
                out[order[k]] = average
            i = j + 1
        return out

    ra, rb = ranks(a), ranks(b)
    n = len(a)
    mean_a, mean_b = sum(ra) / n, sum(rb) / n
    num = sum((ra[i] - mean_a) * (rb[i] - mean_b) for i in range(n))
    den_a = sum((ra[i] - mean_a) ** 2 for i in range(n)) ** 0.5
    den_b = sum((rb[i] - mean_b) ** 2 for i in range(n)) ** 0.5
    return num / (den_a * den_b) if den_a and den_b else 0.0


def extract_first_model(pose_file):
    """
    Write MODEL 1 of a Vina output file to its own single-model PDBQT.

    Vina's --score_only refuses a multi-MODEL ligand file, which is exactly
    what its own --out produces. The MODEL and ENDMDL tags are dropped and only
    the rank-1 block is kept, since that is the pose the screen reported.
    """
    kept, inside, seen = [], False, False
    with open(pose_file, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("MODEL"):
                if seen:
                    break
                inside, seen = True, True
                continue
            if line.startswith("ENDMDL"):
                if inside:
                    break
                continue
            if inside or not seen:
                kept.append(line)
    if not kept:
        return None
    single = os.path.join(OUT, "_rank1_" + os.path.basename(pose_file))
    with open(single, "w", encoding="utf-8") as handle:
        handle.writelines(kept)
    return single


def rescore(pose_file):
    """Vinardo score of an existing pose. Returns None if it cannot be read."""
    single = extract_first_model(pose_file)
    if single is None:
        return None, "could not extract a model from the pose file"
    pose_file = single
    command = [VINA, "--receptor", RECEPTOR, "--ligand", pose_file,
               "--scoring", "vinardo", "--score_only"]
    done = subprocess.run(command, capture_output=True, text=True, timeout=900)
    if done.returncode != 0:
        # Some builds want a box even for score_only.
        command = [VINA, "--receptor", RECEPTOR, "--ligand", pose_file,
                   "--scoring", "vinardo", "--score_only",
                   "--center_x", CENTER[0], "--center_y", CENTER[1],
                   "--center_z", CENTER[2],
                   "--size_x", BOX[0], "--size_y", BOX[1], "--size_z", BOX[2]]
        done = subprocess.run(command, capture_output=True, text=True, timeout=900)
        if done.returncode != 0:
            return None, (done.stderr or done.stdout)[-200:]
    for line in done.stdout.splitlines():
        low = line.lower()
        if "estimated free energy" in low or low.strip().startswith("affinity"):
            for token in line.replace(":", " ").split():
                try:
                    return float(token), None
                except ValueError:
                    continue
    return None, "could not parse a score"


def redock_vinardo(ligand, seed, tag):
    pose = os.path.join(OUT, tag + ".pdbqt")
    command = [VINA, "--receptor", RECEPTOR, "--ligand", ligand,
               "--scoring", "vinardo",
               "--center_x", CENTER[0], "--center_y", CENTER[1],
               "--center_z", CENTER[2],
               "--size_x", BOX[0], "--size_y", BOX[1], "--size_z", BOX[2],
               "--exhaustiveness", "32", "--num_modes", "9",
               "--energy_range", "3", "--seed", str(seed), "--cpu", "4",
               "--out", pose]
    started = time.time()
    done = subprocess.run(command, capture_output=True, text=True, timeout=21600)
    if done.returncode != 0:
        return None
    score = tal_control.parse_score(done.stdout)
    with open(os.path.join(OUT, tag + ".json"), "w", encoding="utf-8") as handle:
        json.dump({"tag": tag, "scoring": "vinardo", "ligand": ligand,
                   "seed": seed, "exhaustiveness": 32,
                   "rank1_score_kcal_mol": score,
                   "seconds": round(time.time() - started, 1),
                   "command": " ".join(command), "talanai_version": "0.1.0"},
                  handle, indent=2, sort_keys=True)
    return score


def main():
    os.makedirs(OUT, exist_ok=True)
    summary_path = os.path.join(SCREEN, "summary.json")
    if not os.path.isfile(summary_path):
        print("  run run_full_screen.py first")
        return 2
    with open(summary_path, encoding="utf-8") as handle:
        screen = json.load(handle)
    vina_means = {r["name"]: r["mean"] for r in screen["results"]}

    print("")
    print("  TEST 1  rescore the same poses with Vinardo")
    print("  Same geometry, different scoring function.")
    print("  " + "-" * 66)
    print("  %-16s %10s %10s" % ("compound", "vina", "vinardo"))

    rows = []
    for name, stem in STEMS:
        pose = os.path.join(SCREEN, "%s_s42.pdbqt" % stem)
        if not os.path.isfile(pose) or name not in vina_means:
            print("  %-16s missing pose or screen entry" % name)
            continue
        value, error = rescore(pose)
        if value is None:
            print("  %-16s rescore failed: %s" % (name, error))
            continue
        rows.append({"name": name, "vina": vina_means[name], "vinardo": value})
        print("  %-16s %10.3f %10.3f" % (name, vina_means[name], value), flush=True)

    if len(rows) < 3:
        print("  too few results to correlate")
        return 1

    rho = spearman([r["vina"] for r in rows], [r["vinardo"] for r in rows])
    print("")
    print("  Spearman rank correlation, Vina against Vinardo: %.3f" % rho)
    print("")
    print("  %-16s %6s %8s" % ("compound", "vina#", "vinardo#"))
    by_vina = sorted(rows, key=lambda r: r["vina"])
    by_vinardo = sorted(rows, key=lambda r: r["vinardo"])
    vina_rank = {r["name"]: i + 1 for i, r in enumerate(by_vina)}
    vinardo_rank = {r["name"]: i + 1 for i, r in enumerate(by_vinardo)}
    for r in by_vina:
        moved = vinardo_rank[r["name"]] - vina_rank[r["name"]]
        flag = "" if moved == 0 else ("  %+d" % moved)
        print("  %-16s %6d %8d%s"
              % (r["name"], vina_rank[r["name"]], vinardo_rank[r["name"]], flag))

    reference = next((r for r in rows if r["name"] == "Acarbose"), None)
    if reference:
        beats_vina = [r["name"] for r in rows
                      if r["name"] != "Acarbose" and r["vina"] < reference["vina"]]
        beats_vinardo = [r["name"] for r in rows
                         if r["name"] != "Acarbose" and r["vinardo"] < reference["vinardo"]]
        print("")
        print("  beats acarbose under Vina    (%d): %s"
              % (len(beats_vina), ", ".join(beats_vina)))
        print("  beats acarbose under Vinardo (%d): %s"
              % (len(beats_vinardo), ", ".join(beats_vinardo)))
        disagree = set(beats_vina) ^ set(beats_vinardo)
        if disagree:
            print("  THE TWO FUNCTIONS DISAGREE ABOUT: %s" % ", ".join(sorted(disagree)))
        else:
            print("  The two functions agree on which compounds beat acarbose.")

    with open(os.path.join(OUT, "rescore.json"), "w", encoding="utf-8") as handle:
        json.dump({"spearman_vina_vs_vinardo": rho, "rows": rows,
                   "note": ("Vinardo rescoring of Vina poses. Same geometry, "
                            "different function. Does not re-search."),
                   "generated": "2026-08-04"}, handle, indent=2, sort_keys=True)

    if "--redock" in sys.argv:
        print("")
        print("  TEST 2  re-dock with Vinardo driving the search, seed 42")
        print("  " + "-" * 66)
        print("  %-16s %10s %10s" % ("compound", "vina mean", "vinardo dock"))
        redocked = []
        for name, stem in STEMS:
            ligand = os.path.join(INPUTS, stem + ".pdbqt")
            if not os.path.isfile(ligand) or name not in vina_means:
                continue
            score = redock_vinardo(ligand, 42, stem + "_vinardo_s42")
            if score is None:
                print("  %-16s FAILED" % name)
                continue
            redocked.append({"name": name, "vina": vina_means[name],
                             "vinardo_dock": score})
            print("  %-16s %10.3f %10.3f" % (name, vina_means[name], score),
                  flush=True)
        if len(redocked) >= 3:
            rho2 = spearman([r["vina"] for r in redocked],
                            [r["vinardo_dock"] for r in redocked])
            print("")
            print("  Spearman, Vina mean against independent Vinardo dock: %.3f" % rho2)
            with open(os.path.join(OUT, "redock.json"), "w", encoding="utf-8") as h:
                json.dump({"spearman": rho2, "rows": redocked,
                           "generated": "2026-08-04"}, h, indent=2, sort_keys=True)

    print("")
    print("  " + "-" * 66)
    print("  Agreement between scoring functions is evidence the ordering is")
    print("  not an artefact of one function. It is NOT evidence the ordering")
    print("  is experimentally correct. Only measured data can say that.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
