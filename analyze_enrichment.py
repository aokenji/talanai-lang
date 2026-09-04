#!/usr/bin/env python3
"""
Enrichment metrics for the docking protocol: can it put known actives above
property-matched lookalikes?

    python analyze_enrichment.py
    python analyze_enrichment.py --selftest    (checks the maths, no data needed)

WHAT IS REPORTED, AND WHY THESE
    ROC AUC       probability a random active outranks a random decoy. 0.5 is
                  chance. Robust, but it weights the whole list equally, and
                  nobody screens the whole list.
    BEDROC a=20   the same question weighted towards the TOP of the list,
                  which is the only part anyone acts on. a=20 corresponds to
                  roughly the first 8 percent carrying most of the weight
                  (Truchon and Bayly 2007).
    EF at x%      how many times more actives appear in the top x% than chance
                  would give. The number a reader recognises.

    Reported together on purpose. A protocol can have a respectable AUC and
    still be useless for a shortlist, and BEDROC is what catches that.

WHY THERE IS NO EF AT 1 PERCENT HERE
    One percent of about 200 molecules is TWO molecules. An EF computed on two
    slots moves in enormous steps and carries no information at this scale.
    Quoting it would be numerology. EF is reported at 5 and 10 percent, and
    the absolute counts are printed beside every ratio so the reader can see
    exactly how few molecules the number rests on.

WHY EVERY NUMBER CARRIES A CONFIDENCE INTERVAL
    There are about ten actives. That is a small set, and every metric here is
    correspondingly noisy. The intervals come from bootstrapping the actives
    with replacement, 10000 times. If an interval includes the chance value,
    the honest reading is that this benchmark could not tell, NOT that the
    protocol failed and NOT that it passed.

TIES
    Vina scores can tie. Ties are given average ranks, which is the neutral
    treatment; breaking them by input order would let the order of the file
    decide the result.

THIS IS NOT A POSE CLAIM
    The redocking gate fails for this protocol. Enrichment measures score
    separation only. A good enrichment result would say the ranking is
    informative; it would say nothing whatever about where a ligand sits.
"""

import argparse
import glob
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.join(HERE, "validation-run", "enrichment")
INPUTS = os.path.join(HERE, "validation-inputs", "enrichment")
OUT = os.path.join(RUNS, "enrichment.json")

BEDROC_ALPHA = 20.0
FRACTIONS = (0.05, 0.10, 0.20)
BOOTSTRAP = 10000


# ── ranking ─────────────────────────────────────────────────────────────────

def average_ranks(scores):
    """
    1-indexed ranks, best (most negative) first, ties sharing the average rank.
    """
    order = sorted(range(len(scores)), key=lambda i: scores[i])
    ranks = [0.0] * len(scores)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and scores[order[j + 1]] == scores[order[i]]:
            j += 1
        shared = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[order[k]] = shared
        i = j + 1
    return ranks


# ── metrics ─────────────────────────────────────────────────────────────────

def roc_auc(active_ranks, n_total, n_active):
    """Mann-Whitney U over ranks; identical to the ROC area."""
    n_decoy = n_total - n_active
    if n_active == 0 or n_decoy == 0:
        return None
    rank_sum = sum(active_ranks)
    u = n_active * n_decoy + n_active * (n_active + 1) / 2.0 - rank_sum
    return u / (n_active * n_decoy)


def bedroc(active_ranks, n_total, n_active, alpha=BEDROC_ALPHA):
    """Truchon and Bayly (2007), equation 36."""
    if n_active == 0 or n_active == n_total:
        return None
    ra = n_active / float(n_total)
    total = sum(math.exp(-alpha * r / n_total) for r in active_ranks)
    rie_denominator = (1.0 - math.exp(-alpha)) / (n_total * (math.exp(alpha / n_total) - 1.0))
    rie = (total / n_active) / rie_denominator
    factor = (ra * math.sinh(alpha / 2.0)
              / (math.cosh(alpha / 2.0) - math.cosh(alpha / 2.0 - alpha * ra)))
    return rie * factor + 1.0 / (1.0 - math.exp(alpha * (1.0 - ra)))


def enrichment_factor(active_ranks, n_total, n_active, fraction):
    cutoff = max(1, int(round(fraction * n_total)))
    hits = sum(1 for r in active_ranks if r <= cutoff)
    expected = n_active * cutoff / float(n_total)
    return {
        "fraction": fraction,
        "top_n": cutoff,
        "actives_found": hits,
        "actives_expected_by_chance": round(expected, 2),
        "ef": round(hits / expected, 3) if expected > 0 else None,
    }


def percentile(values, p):
    if not values:
        return None
    ordered = sorted(values)
    k = (len(ordered) - 1) * p
    lo, hi = int(math.floor(k)), int(math.ceil(k))
    if lo == hi:
        return ordered[lo]
    return ordered[lo] + (ordered[hi] - ordered[lo]) * (k - lo)


def bootstrap_ci(scores, labels, metric, rounds=BOOTSTRAP):
    """
    Resample the ACTIVES with replacement, keeping the decoys fixed. The
    actives are the scarce, noisy side; the decoys are not what the interval
    is about.
    """
    active_idx = [i for i, a in enumerate(labels) if a]
    decoy_scores = [scores[i] for i, a in enumerate(labels) if not a]
    if not active_idx or not decoy_scores:
        return None, None
    n_active = len(active_idx)
    n_total = n_active + len(decoy_scores)

    # Deterministic LCG: a fixed sequence, so the interval is reproducible.
    state = 0xF00D
    draws = []
    for _ in range(rounds):
        sample = []
        for _ in range(n_active):
            state = (1103515245 * state + 12345) & 0x7FFFFFFF
            sample.append(scores[active_idx[state % n_active]])
        merged = sample + decoy_scores
        flags = [True] * n_active + [False] * len(decoy_scores)
        ranks = average_ranks(merged)
        ar = [ranks[i] for i, f in enumerate(flags) if f]
        value = metric(ar, n_total, n_active)
        if value is not None:
            draws.append(value)
    return percentile(draws, 0.025), percentile(draws, 0.975)


# ── data ────────────────────────────────────────────────────────────────────

def load_records():
    rows = []
    for path in sorted(glob.glob(os.path.join(RUNS, "*.json"))):
        if os.path.basename(path) == "enrichment.json":
            continue
        with open(path, encoding="utf-8") as handle:
            record = json.load(handle)
        score = record.get("rank1_score_kcal_mol")
        if score is None:
            continue
        rows.append({
            "tag": record.get("tag"),
            "role": record.get("role"),
            "name": record.get("name"),
            "score": float(score),
            "all_chairs": record.get("all_chairs", True),
            "matched_to": record.get("matched_to"),
        })
    return rows


def evaluate(rows, label):
    scores = [r["score"] for r in rows]
    labels = [r["role"] == "active" for r in rows]
    n_total, n_active = len(rows), sum(labels)
    if n_active == 0 or n_active == n_total:
        return None

    ranks = average_ranks(scores)
    active_ranks = [ranks[i] for i, a in enumerate(labels) if a]

    auc = roc_auc(active_ranks, n_total, n_active)
    bed = bedroc(active_ranks, n_total, n_active)
    auc_lo, auc_hi = bootstrap_ci(scores, labels, roc_auc)
    bed_lo, bed_hi = bootstrap_ci(scores, labels, bedroc)

    return {
        "set": label,
        "n_total": n_total,
        "n_active": n_active,
        "n_decoy": n_total - n_active,
        "decoys_per_active": round((n_total - n_active) / float(n_active), 1),
        "roc_auc": round(auc, 4),
        "roc_auc_ci95": [round(auc_lo, 4), round(auc_hi, 4)] if auc_lo is not None else None,
        "roc_auc_chance": 0.5,
        "bedroc_alpha20": round(bed, 4),
        "bedroc_ci95": [round(bed_lo, 4), round(bed_hi, 4)] if bed_lo is not None else None,
        "enrichment_factors": [
            enrichment_factor(active_ranks, n_total, n_active, f) for f in FRACTIONS],
        "best_ranked_active": min(active_ranks),
        "worst_ranked_active": max(active_ranks),
    }


def verdict(result):
    """Say what the interval supports, and nothing more."""
    if result is None:
        return "no result"
    ci = result.get("roc_auc_ci95")
    if not ci:
        return "no interval"
    if ci[0] > 0.5:
        return ("SEPARATES. The 95% interval for ROC AUC lies entirely above "
                "chance, so on this set the score does put actives above "
                "matched decoys.")
    if ci[1] < 0.5:
        return ("ANTI-CORRELATED. The interval lies entirely BELOW chance. "
                "The score is actively preferring the decoys, which is a "
                "result and must be reported as one.")
    return ("COULD NOT TELL. The interval spans chance. With this few actives "
            "the benchmark does not have the resolution to decide, and that "
            "is the finding. It is not a pass and it is not a failure.")


# ── self test ───────────────────────────────────────────────────────────────

def selftest():
    print("")
    print("  SELF TEST of the metric implementations")
    print("  " + "-" * 74)
    n_active, n_decoy = 10, 190
    n_total = n_active + n_decoy
    ok = True

    perfect = list(range(1, n_active + 1))
    worst = list(range(n_decoy + 1, n_total + 1))
    middle = [n_decoy / 2.0 + i for i in range(n_active)]

    for name, ranks, want_auc in (("perfect", perfect, 1.0),
                                  ("worst", worst, 0.0)):
        auc = roc_auc(ranks, n_total, n_active)
        bed = bedroc(ranks, n_total, n_active)
        good = abs(auc - want_auc) < 1e-9
        ok = ok and good
        print("    %-8s ROC AUC %.4f (want %.1f) %s   BEDROC %.4f"
              % (name, auc, want_auc, "ok" if good else "FAIL", bed))

    auc_mid = roc_auc(middle, n_total, n_active)
    near_half = abs(auc_mid - 0.5) < 0.05
    ok = ok and near_half
    print("    %-8s ROC AUC %.4f (want ~0.5) %s"
          % ("middle", auc_mid, "ok" if near_half else "FAIL"))

    bed_perfect = bedroc(perfect, n_total, n_active)
    bed_worst = bedroc(worst, n_total, n_active)
    ordered = bed_perfect > 0.9 and bed_worst < 0.1
    ok = ok and ordered
    print("    BEDROC bounds: perfect %.4f, worst %.4f %s"
          % (bed_perfect, bed_worst, "ok" if ordered else "FAIL"))

    ranks = average_ranks([-9.0, -8.0, -8.0, -7.0])
    tied = ranks == [1.0, 2.5, 2.5, 4.0]
    ok = ok and tied
    print("    average ranks on a tie: %s %s" % (ranks, "ok" if tied else "FAIL"))

    ef = enrichment_factor(perfect, n_total, n_active, 0.05)
    ef_ok = ef["actives_found"] == n_active and ef["top_n"] == 10
    ok = ok and ef_ok
    print("    EF 5%% on a perfect ranking: %s of %s in top %d, EF %.2f %s"
          % (ef["actives_found"], n_active, ef["top_n"], ef["ef"],
             "ok" if ef_ok else "FAIL"))

    print("")
    print("    %s" % ("ALL CHECKS PASSED" if ok else "*** SOME CHECKS FAILED ***"))
    print("")
    return 0 if ok else 1


# ── main ────────────────────────────────────────────────────────────────────

def show(result):
    if result is None:
        print("    not computable for this subset")
        return
    print("    %d molecules: %d actives, %d decoys (%.1f per active)"
          % (result["n_total"], result["n_active"], result["n_decoy"],
             result["decoys_per_active"]))
    ci = result["roc_auc_ci95"]
    print("    ROC AUC        %.4f   95%% CI [%.4f, %.4f]   chance 0.5"
          % (result["roc_auc"], ci[0], ci[1]))
    bci = result["bedroc_ci95"]
    print("    BEDROC a=20    %.4f   95%% CI [%.4f, %.4f]"
          % (result["bedroc_alpha20"], bci[0], bci[1]))
    for ef in result["enrichment_factors"]:
        print("    EF %3d%%        %.2f   (%d of %d actives in the top %d, "
              "chance would give %.2f)"
              % (ef["fraction"] * 100, ef["ef"], ef["actives_found"],
                 result["n_active"], ef["top_n"],
                 ef["actives_expected_by_chance"]))
    print("")
    print("    %s" % verdict(result))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selftest", action="store_true")
    args = parser.parse_args()

    if args.selftest:
        sys.exit(selftest())

    if not os.path.isdir(RUNS):
        raise SystemExit("no %s - run run_enrichment.py first" % RUNS)

    rows = load_records()
    if not rows:
        raise SystemExit("no completed run records in %s" % RUNS)

    n_active = sum(1 for r in rows if r["role"] == "active")
    print("")
    print("  ENRICHMENT RESULT")
    print("  " + "=" * 74)
    print("  %d scored molecules, %d actives" % (len(rows), n_active))

    print("")
    print("  ALL LIGANDS")
    print("  " + "-" * 74)
    overall = evaluate(rows, "all")
    show(overall)

    chairs = [r for r in rows if r["all_chairs"]]
    dropped = len(rows) - len(chairs)
    chair_only = None
    print("")
    print("  ALL-CHAIR LIGANDS ONLY  (%d dropped)" % dropped)
    print("  " + "-" * 74)
    if dropped == 0:
        print("    every ligand got all-chair rings; identical to the set above")
        chair_only = overall
    else:
        chair_only = evaluate(chairs, "all_chairs_only")
        show(chair_only)
        print("")
        print("    A ligand Vina was handed with a strained ring was set up to")
        print("    lose. If the two blocks disagree, the preparation is doing")
        print("    part of the discriminating and that has to be said out loud.")

    payload = {
        "generated_from": RUNS,
        "bedroc_alpha": BEDROC_ALPHA,
        "bootstrap_rounds": BOOTSTRAP,
        "ties": "average ranks",
        "results": [r for r in (overall, chair_only) if r],
        "scope": ("Score separation only. The redocking gate does not pass "
                  "for this protocol, so nothing here supports a pose or "
                  "interaction claim."),
        "small_sample_caveat": (
            "About ten actives. Every interval here is wide, and an interval "
            "spanning 0.5 means the benchmark could not tell, not that the "
            "protocol failed."),
        "no_ef1_percent": (
            "EF at 1 percent is omitted deliberately: one percent of this set "
            "is two molecules, and a ratio computed on two slots carries no "
            "information."),
    }
    with open(OUT, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)

    print("")
    print("  written to %s" % OUT)
    print("")


if __name__ == "__main__":
    main()
