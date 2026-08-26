#!/usr/bin/env python3
"""
Acarbose seed sweep: how often does the search miss?

WHY
    The 2026-08-25 replication returned acarbose at -8.720, -8.708 and -8.449 on
    seeds 45/46/47, a spread of 0.271 against the 0.010 recorded on the
    published seeds 42/43/44. Reading that as "acarbose has high variance" is
    wrong, and the per-seed numbers say so:

        42 -8.700   43 -8.710   44 -8.708
        45 -8.720   46 -8.708   47 -8.449

    Five of six agree inside 0.020. One missed by 0.26. That is not scatter, it
    is an OCCASIONAL SEARCH FAILURE: exhaustiveness 32 usually finds the basin
    for this ligand and sometimes does not. The distinction matters, because a
    failure rate can be measured and quoted, while "high variance" cannot be
    acted on.

    Acarbose is the reference every "beats acarbose" verdict is measured
    against, and it is the largest, most flexible ligand in the set at 46 heavy
    atoms. If the search misses once in six, margins near the reference are less
    trustworthy than the published spread of 0.010 suggests.

WHAT THIS DOES
    Six more seeds, 48 to 53, identical protocol and preparation. Combined with
    the six already on disk that gives n = 12, enough to say something about how
    often the miss happens rather than that it happened once.

WHAT IT DOES NOT DO
    Change any published value. This is a measurement of search behaviour, not a
    re-screen. compounds.js is not touched.
"""

import json
import os
import statistics
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import run_replication as R                      # noqa: E402

NEW_SEEDS = (48, 49, 50, 51, 52, 53)
OUT = os.path.join(HERE, "validation-run", "acarbose-seeds")
EXISTING = [
    os.path.join(HERE, "validation-run", "full-screen"),
    os.path.join(HERE, "validation-run", "replication"),
]


def collect_existing():
    found = {}
    for folder in EXISTING:
        if not os.path.isdir(folder):
            continue
        for name in os.listdir(folder):
            if not (name.startswith("acarbose_s") and name.endswith(".json")):
                continue
            with open(os.path.join(folder, name), encoding="utf-8") as handle:
                record = json.load(handle)
            seed = record.get("seed")
            score = record.get("rank1_score_kcal_mol")
            if seed is not None and score is not None:
                found[int(seed)] = float(score)
    return found


def main():
    os.makedirs(OUT, exist_ok=True)
    R.OUT = OUT   # write new records here, leave the earlier runs untouched

    name, stem, published, formula, confs, smiles = [
        c for c in R.COMPOUNDS if c[0] == "Acarbose"][0]

    print("")
    print("  Acarbose seed sweep, exhaustiveness %d" % R.EXHAUSTIVENESS)
    print("  " + "-" * 68)

    prepared = R.ringaware.prepare_ring_aware(smiles, confs)
    path = os.path.join(R.INPUTS, "acarbose.pdbqt")
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(prepared["pdbqt"])
    print("    prepared, rule: %s" % prepared["rule"])
    print("")

    scores = collect_existing()
    print("    already on disk: %s" % ", ".join(
        "%d=%.3f" % (s, scores[s]) for s in sorted(scores)))
    print("")

    for seed in NEW_SEEDS:
        started = time.time()
        score = R.dock(path, seed, "acarbose_s%d" % seed)
        if score is None:
            print("    seed %d FAILED" % seed)
            continue
        scores[seed] = float(score)
        print("    seed %d  %8.3f   (%.0f s)" % (seed, score, time.time() - started))

    values = [scores[s] for s in sorted(scores)]
    best = min(values)                      # most negative = best score
    # A "miss" is a run that lands materially above the basin the others agree on.
    # 0.10 is this project's run-to-run convention, used here as the threshold for
    # "did not reach the same answer", not as a precision claim about Vina.
    misses = [s for s in sorted(scores) if scores[s] - best > R.REPRODUCIBILITY_LIMIT]
    basin = [scores[s] for s in sorted(scores) if s not in misses]

    summary = {
        "what": "acarbose seed sweep at exhaustiveness 32",
        "date": time.strftime("%Y-%m-%d"),
        "question": ("is the wide replication spread scatter, or an occasional "
                     "search failure?"),
        "n_seeds": len(values),
        "seeds": {str(s): scores[s] for s in sorted(scores)},
        "best_kcal_mol": round(best, 3),
        "full_spread_kcal_mol": round(max(values) - min(values), 3),
        "basin_n": len(basin),
        "basin_mean_kcal_mol": round(statistics.fmean(basin), 3),
        "basin_spread_kcal_mol": round(max(basin) - min(basin), 3),
        "miss_seeds": misses,
        "miss_rate": "%d/%d" % (len(misses), len(values)),
        "miss_threshold_kcal_mol": R.REPRODUCIBILITY_LIMIT,
        "miss_threshold_note": ("this project's run-to-run convention, used as a "
                                "did-not-reach-the-same-answer threshold, not a "
                                "published Vina precision figure"),
        "published_kcal_mol": published,
        "published_spread_kcal_mol": 0.010,
        "not_validated": ("the redocking gate does not pass; this measures search "
                          "behaviour, not pose accuracy"),
        "does_not_change": ("no published value is altered by this run; "
                            "src/data/compounds.js is untouched"),
    }
    with open(os.path.join(OUT, "summary.json"), "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)

    print("")
    print("  " + "-" * 68)
    print("  n = %d seeds" % len(values))
    print("  basin: %d seeds, mean %.3f, spread %.3f"
          % (len(basin), statistics.fmean(basin), max(basin) - min(basin)))
    print("  misses: %s  (rate %s, threshold %.2f kcal/mol)"
          % (misses or "none", summary["miss_rate"], R.REPRODUCIBILITY_LIMIT))
    print("  full spread across all seeds: %.3f" % (max(values) - min(values)))
    print("  Written to %s" % OUT)
    print("")


if __name__ == "__main__":
    main()
