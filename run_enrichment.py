#!/usr/bin/env python3
"""
Dock the actives and their property-matched decoys, for the enrichment panel.

    python make_decoys.py --library <a SMILES file>     (run this first)
    python run_enrichment.py
    python analyze_enrichment.py

RESUMABLE ON PURPOSE
    This is hours of work, and a laptop that sleeps mid-run must not lose it.
    Every ligand writes its own JSON record on completion and a rerun skips
    anything already recorded. Interrupt it freely.

WHY ONE SEED, WHEN EVERY OTHER RUN IN THIS PROJECT USES THREE
    Because the comparison has to be symmetric, and that matters more here
    than precision does.

    The published affinities are means of seeds 42/43/44. Docking 200 decoys
    on three seeds would be roughly 600 runs and several days. Docking decoys
    on ONE seed and comparing them against three-seed active means would be
    worse than slow, it would be WRONG: run_acarbose_seeds.py measured a
    search-failure rate of about 1 in 12 at exhaustiveness 32, and a miss
    makes a ligand look weaker than it is. Applied to decoys only, that noise
    pushes decoys down the ranking and inflates enrichment. The bias would run
    in exactly the direction that flatters the result.

    So both sides get one seed, seed 42, and the actives are re-docked here
    rather than read from compounds.js. Vina is deterministic on identical
    inputs, so the actives' single-seed numbers are reproducible, and
    analyze_enrichment.py additionally reports the three-seed means as a
    secondary check that the conclusion does not turn on the choice.

    Exhaustiveness stays at 32. If this is too slow, cut the decoy COUNT.
    Lowering exhaustiveness would benchmark a protocol the thesis did not use,
    which is worse than not benchmarking at all.

PREPARATION IS IDENTICAL ON BOTH SIDES
    Same ring-aware conformer selection, same conformer count, same meeko
    settings for actives and decoys. Vina holds rings rigid, so a decoy handed
    a twist-boat is a decoy that was set up to lose. Each record carries
    whether its rings came out as chairs, and the analysis reports the result
    with and without the ones that did not, so that bias stays visible instead
    of being argued about.

NOTHING HERE IS A POSE CLAIM
    The redocking gate fails for this protocol. Enrichment asks only whether
    the SCORE separates actives from lookalikes. Every record carries
    "validated": false, like every other run record in this project.
"""

import argparse
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

INPUTS = os.path.join(HERE, "validation-inputs", "enrichment")
OUT = os.path.join(HERE, "validation-run", "enrichment")
DECOYS = os.path.join(INPUTS, "decoys.json")

CENTER = ["21.52", "-7.7", "23.55"]
BOX = ["30", "30", "30"]
EXHAUSTIVENESS = 32
SEED = 42
CONFORMERS = 150


def slug(text):
    keep = [c.lower() if c.isalnum() else "_" for c in text]
    return "".join(keep).strip("_")


def prepare(smiles, stem, confs):
    """Ring-aware prep. Returns (pdbqt path, all-chair flag, note)."""
    path = os.path.join(INPUTS, "ligands", stem + ".pdbqt")
    meta = path + ".json"
    if os.path.isfile(path) and os.path.isfile(meta):
        with open(meta, encoding="utf-8") as handle:
            cached = json.load(handle)
        return path, cached["all_chairs"], cached["note"]

    os.makedirs(os.path.dirname(path), exist_ok=True)
    result = ringaware.prepare_ring_aware(smiles, confs)
    thetas = result["thetas"]
    chairs = ringaware.all_chairs(thetas) if thetas else True
    note = ("no saturated six-ring" if not thetas
            else ("all chairs" if chairs else "ring(s) not chairs"))
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(result["pdbqt"])
    with open(meta, "w", encoding="utf-8") as handle:
        json.dump({"smiles": smiles, "all_chairs": chairs, "note": note,
                   "thetas": thetas, "n_qualifying": result["n_qualifying"],
                   "n_conformers": result["n_conformers"]},
                  handle, indent=2, sort_keys=True)
    return path, chairs, note


def dock(ligand, tag, role, name, extra):
    record_path = os.path.join(OUT, tag + ".json")
    if os.path.isfile(record_path):
        with open(record_path, encoding="utf-8") as handle:
            return json.load(handle), True

    pose = os.path.join(OUT, "poses", tag + ".pdbqt")
    os.makedirs(os.path.dirname(pose), exist_ok=True)
    command = [VINA, "--receptor", RECEPTOR, "--ligand", ligand,
               "--center_x", CENTER[0], "--center_y", CENTER[1],
               "--center_z", CENTER[2],
               "--size_x", BOX[0], "--size_y", BOX[1], "--size_z", BOX[2],
               "--exhaustiveness", str(EXHAUSTIVENESS), "--num_modes", "9",
               "--energy_range", "3", "--seed", str(SEED), "--cpu", "4",
               "--out", pose]
    started = time.time()
    try:
        done = subprocess.run(command, capture_output=True, text=True,
                              timeout=21600)
    except subprocess.TimeoutExpired:
        record = {"tag": tag, "role": role, "name": name, "error": "timeout",
                  "seconds": round(time.time() - started, 1)}
        with open(record_path, "w", encoding="utf-8") as handle:
            json.dump(record, handle, indent=2, sort_keys=True)
        return record, False

    record = {
        "tag": tag, "role": role, "name": name, "ligand": ligand,
        "receptor": RECEPTOR, "box_center": CENTER, "box_size": BOX,
        "exhaustiveness": EXHAUSTIVENESS, "seed": SEED,
        "preparation": "ring-aware conformer selection, meeko",
        "conformers_embedded": CONFORMERS,
        "seconds": round(time.time() - started, 1),
        "returncode": done.returncode, "command": " ".join(command),
        "validated": False,
        "scope": ("Enrichment run. The redocking gate does not pass for this "
                  "protocol, so this is a claim about score separation only, "
                  "never about pose."),
        "pose_file": pose if os.path.isfile(pose) else None,
    }
    record.update(extra)
    if done.returncode != 0:
        record["error"] = (done.stderr or done.stdout)[-800:]
    else:
        record["rank1_score_kcal_mol"] = tal_control.parse_score(done.stdout)

    with open(record_path, "w", encoding="utf-8") as handle:
        json.dump(record, handle, indent=2, sort_keys=True)
    return record, False


def build_worklist():
    with open(DECOYS, encoding="utf-8") as handle:
        manifest = json.load(handle)

    work = []
    for active in manifest["actives"]:
        work.append({"role": "active", "name": active["name"],
                     "smiles": active["smiles"],
                     "tag": "active__" + slug(active["name"]),
                     "extra": {"pubchem_cid": active["cid"],
                               "n_yeast_ic50_sources":
                                   active["n_yeast_ic50_sources"]}})
    for active_name, decoys in sorted(manifest["decoys"].items()):
        for i, decoy in enumerate(decoys, 1):
            work.append({"role": "decoy", "name": "%s decoy %d" % (active_name, i),
                         "smiles": decoy["smiles"],
                         "tag": "decoy__%s__%03d" % (slug(active_name), i),
                         "extra": {"matched_to": active_name,
                                   "library_id": decoy["library_id"]}})
    return manifest, work


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--limit", type=int, default=0,
                        help="stop after this many NEW dockings (for a dry run)")
    parser.add_argument("--prep-only", action="store_true",
                        help="prepare ligands and stop, no docking")
    args = parser.parse_args()

    if not os.path.isfile(DECOYS):
        raise SystemExit("no %s - run make_decoys.py first" % DECOYS)
    os.makedirs(OUT, exist_ok=True)

    manifest, work = build_worklist()
    actives = [w for w in work if w["role"] == "active"]
    decoys = [w for w in work if w["role"] == "decoy"]

    print("")
    print("  ENRICHMENT RUN")
    print("  " + "=" * 74)
    print("  %d actives, %d decoys, %d dockings total"
          % (len(actives), len(decoys), len(work)))
    print("  exhaustiveness %d, seed %d, one seed on BOTH sides for symmetry"
          % (EXHAUSTIVENESS, SEED))
    already = sum(1 for w in work if os.path.isfile(os.path.join(OUT, w["tag"] + ".json")))
    print("  %d already recorded, %d to go" % (already, len(work) - already))
    print("")

    if args.prep_only:
        print("  PREPARING EVERY LIGAND  (--prep-only)")
        print("  " + "-" * 74)
        not_chairs = []
        for i, item in enumerate(work, 1):
            _path, chairs, _note = prepare(item["smiles"], item["tag"], CONFORMERS)
            if not chairs:
                not_chairs.append(item["tag"])
            if i % 25 == 0 or i == len(work):
                print("    prepared %d / %d" % (i, len(work)))
        print("")
        print("    %d of %d could not be given all-chair rings."
              % (len(not_chairs), len(work)))
        return

    # Preparation is LAZY, done immediately before each docking, so that
    # --limit 3 costs three preparations rather than a hundred and eighty.
    # prepare() caches to disk, so nothing is repeated across runs.
    print("  DOCKING  (ligands are prepared as they come up)")
    print("  " + "-" * 74)
    done_now, elapsed, not_chairs = 0, 0.0, []
    for i, item in enumerate(work, 1):
        record_path = os.path.join(OUT, item["tag"] + ".json")
        if os.path.isfile(record_path):
            continue

        path, chairs, note = prepare(item["smiles"], item["tag"], CONFORMERS)
        item["extra"]["all_chairs"] = chairs
        item["extra"]["ring_note"] = note
        if not chairs:
            not_chairs.append(item["tag"])

        record, _cached = dock(path, item["tag"], item["role"], item["name"],
                               item["extra"])
        done_now += 1
        elapsed += record.get("seconds", 0)
        score = record.get("rank1_score_kcal_mol")
        remaining = len(work) - i
        eta = (elapsed / done_now) * remaining / 3600
        print("    [%4d/%4d] %-34s %8s  %5.0fs   eta %4.1f h  %s"
              % (i, len(work), item["tag"][:34],
                 ("%.3f" % score) if score is not None else "FAILED",
                 record.get("seconds", 0), eta, note))
        if args.limit and done_now >= args.limit:
            print("")
            print("    --limit of %d reached. Rerun without it to continue;"
                  % args.limit)
            print("    everything already recorded is skipped.")
            break

    print("")
    if not_chairs:
        print("  %d ligand(s) this session could not be given all-chair rings."
              % len(not_chairs))
        print("  They are docked anyway. analyze_enrichment.py reports the")
        print("  result with and without them, because a ligand handed a")
        print("  strained ring was set up to lose.")
        print("")
    if done_now:
        print("  %d dockings this session, %.1f s each on average."
              % (done_now, elapsed / done_now))
        print("  Projected for the remaining %d: %.1f hours."
              % (len(work) - done_now,
                 (elapsed / done_now) * (len(work) - done_now) / 3600))
    print("  Records in %s" % OUT)
    print("  Next: python analyze_enrichment.py")
    print("")


if __name__ == "__main__":
    main()
