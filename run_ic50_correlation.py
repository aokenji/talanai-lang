#!/usr/bin/env python3
"""
Does the docking ranking track measured potency?

THE TEST THAT MATTERS
    Redocking validates pose prediction. Consensus validates robustness to the
    scoring function. Neither says whether the ranking is RIGHT. The only thing
    that does is comparing predictions against measured inhibition.

    So: rank correlation between predicted binding score and published IC50
    against yeast alpha-glucosidase, the assay that matches PDB 3A4A.

    SIGN CONVENTION, which is easy to get backwards and was, on the first
    version of this script:

        a potent compound has a VERY NEGATIVE score and a LOW IC50
        a weak compound has a LESS NEGATIVE score and a HIGH IC50

    Both quantities are therefore low for potent compounds and high for weak
    ones, so a correct prediction gives a POSITIVE Spearman. Positive is good.
    Zero means the ranking carries no information about potency. Negative means
    the ranking is anti-predictive: it prefers the compounds that measure
    weaker.

THE HONESTY PROBLEM, AND HOW IT IS HANDLED
    Published IC50s for the same compound against the same nominal enzyme vary
    by 60 to 90 fold between labs. There is no single "true" value to correlate
    against. Picking one would be choosing the answer.

    So the correlation is computed three ways, using the LOWEST, the GEOMETRIC
    MEAN, and the HIGHEST reported value for each compound. If the three agree,
    the conclusion is robust to that choice. If they disagree, the experimental
    literature is too inconsistent to validate against, and that is itself the
    finding.

    python run_ic50_correlation.py
"""

import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
IC50 = os.path.join(HERE, "validation-run", "ic50", "ic50.json")
SCREEN = os.path.join(HERE, "validation-run", "full-screen", "summary.json")
CONSENSUS = os.path.join(HERE, "validation-run", "consensus", "redock.json")
OUT = os.path.join(HERE, "validation-run", "ic50")


def spearman(a, b):
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
    ma, mb = sum(ra) / n, sum(rb) / n
    num = sum((ra[i] - ma) * (rb[i] - mb) for i in range(n))
    da = sum((ra[i] - ma) ** 2 for i in range(n)) ** 0.5
    db = sum((rb[i] - mb) ** 2 for i in range(n)) ** 0.5
    return num / (da * db) if da and db else 0.0


def micromolar(entry):
    """IC50 in uM, or None. Only accepts values already in uM."""
    unit = str(entry.get("ic50_unit", "")).lower()
    value = entry.get("ic50_value")
    if value is None:
        return None
    if unit in ("um", "µm", "umol/l", "micromolar"):
        try:
            return float(value)
        except (TypeError, ValueError):
            return None
    return None


def main():
    for path in (IC50, SCREEN):
        if not os.path.isfile(path):
            print("  missing %s" % path)
            return 2

    with open(IC50, encoding="utf-8") as handle:
        lit = json.load(handle)
    with open(SCREEN, encoding="utf-8") as handle:
        screen = json.load(handle)
    vina = {r["name"]: r["mean"] for r in screen["results"]}

    vinardo = {}
    if os.path.isfile(CONSENSUS):
        with open(CONSENSUS, encoding="utf-8") as handle:
            vinardo = {r["name"]: r["vinardo_dock"] for r in json.load(handle)["rows"]}

    aliases = {"betulinic_acid": "Betulinic Acid", "ursolic_acid": "Ursolic Acid",
               "oleanolic_acid": "Oleanolic Acid"}

    rows = []
    for key, block in lit.get("compounds", {}).items():
        name = aliases.get(key.lower().replace(" ", "_"), key)
        if name not in vina:
            name = key.replace("_", " ").title()
        if name not in vina:
            continue
        values = [v for v in (micromolar(e) for e in block.get("yeast_values", []))
                  if v is not None and v > 0]
        if not values:
            continue
        geo = math.exp(sum(math.log(v) for v in values) / len(values))
        rows.append({"name": name, "n_sources": len(values),
                     "lo": min(values), "geo": geo, "hi": max(values),
                     "fold_spread": max(values) / min(values),
                     "vina": vina[name], "vinardo": vinardo.get(name)})

    if len(rows) < 4:
        print("  only %d compounds have usable data; too few to correlate"
              % len(rows))
        return 1

    rows.sort(key=lambda r: r["geo"])

    print("")
    print("  Predicted score against measured yeast alpha-glucosidase IC50")
    print("  Potent = very negative score AND low IC50, so both are low for a")
    print("  good compound: a correct prediction gives a POSITIVE Spearman.")
    print("  Zero means no information. Negative means anti-predictive.")
    print("  " + "-" * 76)
    print("  %-16s %8s %9s %9s %9s %8s %9s"
          % ("compound", "n", "IC50 lo", "IC50 geo", "IC50 hi", "fold", "vina"))
    for r in rows:
        print("  %-16s %8d %9.2f %9.2f %9.2f %7.0fx %9.3f"
              % (r["name"], r["n_sources"], r["lo"], r["geo"], r["hi"],
                 r["fold_spread"], r["vina"]))

    print("")
    print("  " + "-" * 76)
    print("  Spearman, score against log IC50, computed three ways")
    print("")
    results = {}
    for label, key in (("lowest reported", "lo"), ("geometric mean", "geo"),
                       ("highest reported", "hi")):
        ic50 = [math.log10(r[key]) for r in rows]
        rho_vina = spearman([r["vina"] for r in rows], ic50)
        results["vina_" + key] = rho_vina
        line = "  %-18s  Vina %+.3f" % (label, rho_vina)
        if all(r["vinardo"] is not None for r in rows):
            rho_vinardo = spearman([r["vinardo"] for r in rows], ic50)
            results["vinardo_" + key] = rho_vinardo
            line += "     Vinardo %+.3f" % rho_vinardo
        print(line)

    vina_values = [v for k, v in results.items() if k.startswith("vina_")]
    swing = max(vina_values) - min(vina_values)
    flips = min(vina_values) < 0 < max(vina_values)

    print("")
    print("  " + "-" * 76)
    print("  Correlation swings by %.3f depending only on which published"
          % swing)
    print("  value is chosen for each compound.")
    print("")
    for label, values in (("Vina", vina_values),
                          ("Vinardo", [v for k, v in results.items()
                                       if k.startswith("vinardo_")])):
        if not values:
            continue
        lo, hi = min(values), max(values)
        if lo < 0 < hi:
            verdict = "CHANGES SIGN with the literature value chosen"
        elif hi < -0.3:
            verdict = "ANTI-PREDICTIVE: it prefers the weaker compounds"
        elif abs(lo) < 0.25 and abs(hi) < 0.25:
            verdict = "NO PREDICTIVE VALUE: indistinguishable from chance"
        elif lo > 0.5:
            verdict = "predictive"
        else:
            verdict = "weakly predictive at best"
        print("  %-8s %+.3f to %+.3f   ->  %s" % (label, lo, hi, verdict))

    if flips:
        print("")
        print("  A correlation that changes sign purely by choosing a different")
        print("  but equally citable published value cannot validate anything.")

    with open(os.path.join(OUT, "correlation.json"), "w", encoding="utf-8") as h:
        json.dump({"generated": "2026-08-05", "rows": rows,
                   "spearman": results, "swing": swing, "changes_sign": flips,
                   "method": ("Spearman between predicted score and log10 IC50 "
                              "in uM, yeast alpha-glucosidase only, computed "
                              "with the lowest, geometric-mean and highest "
                              "published value per compound"),
                   "caveat": ("Published IC50s vary by up to 90-fold for the "
                              "same compound and enzyme. No single true value "
                              "exists to correlate against.")},
                  h, indent=2, sort_keys=True)
    print("")
    print("  written to %s" % os.path.join(OUT, "correlation.json"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
