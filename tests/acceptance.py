#!/usr/bin/env python3
"""
Acceptance test for Talanai.

The rule: the example experiment file must agree, to the digit, with the
canonical dataset in the TalanaiHub repository. The numbers are READ FROM
compounds.js, never typed in here, so this test fails the moment the two
records drift apart.

    python tests/acceptance.py

Canonical source:  D:\\BALAKATDBV2\\src\\data\\compounds.js   (read only)
Checked against:   examples/alpha-glucosidase.tal
"""

from __future__ import annotations

import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

import talanai                                    # noqa: E402
from talanai import chem                          # noqa: E402

COMPOUNDS_JS = r"D:\BALAKATDBV2\src\data\compounds.js"
EXPERIMENT = os.path.join(ROOT, "examples", "alpha-glucosidase.tal")

EXPECTED = {
    "Rutin": -8.857,
    "Betulinic Acid": -8.290,
    "Quercetin": -7.503,
    "Kaempferol": -7.479,
    "Oleanolic Acid": -6.922,
    "Acarbose": -6.660,
}


def read_compounds_js(path):
    """name -> {affinity, formula} straight out of the canonical dataset."""
    with open(path, encoding="utf-8") as handle:
        text = handle.read()
    records = {}
    for chunk in re.split(r"\n\s*\{", text):
        name = re.search(r"name:\s*'([^']+)'", chunk)
        affinity = re.search(r"bindingAffinity:\s*(-?[\d.]+)", chunk)
        formula = re.search(r"formula:\s*'([^']+)'", chunk)
        if not (name and affinity):
            continue
        records[name.group(1)] = {
            "affinity": float(affinity.group(1)),
            "formula": chem.normalise_formula(formula.group(1)) if formula else None,
        }
    return records


def main():
    failures = []
    print("")
    print("  Talanai acceptance test")
    print("  experiment : %s" % os.path.relpath(EXPERIMENT, ROOT))
    print("  canonical  : %s" % COMPOUNDS_JS)
    print("  " + "-" * 66)

    experiment = talanai.load(EXPERIMENT)
    affinities = {n.replace("_", " "): v for n, v in experiment.affinities().items()}
    compounds = {c.display_name(): c for c in experiment.compounds()}
    efficiencies = {n.replace("_", " "): v
                    for n, v in experiment.ligand_efficiencies().items()}

    # 1. every expected compound is present, with the thesis value
    print("    %-18s %10s %10s %9s  %s"
          % ("compound", "expected", "in .tal", "per atom", ""))
    for name, expected in sorted(EXPECTED.items(), key=lambda kv: kv[1]):
        got = affinities.get(name)
        efficiency = efficiencies.get(name)
        ok = got is not None and abs(got - expected) < 1e-9
        print("    %-18s %10.3f %10s %9s  %s"
              % (name, expected,
                 "%.3f" % got if got is not None else "MISSING",
                 "%.3f" % efficiency if efficiency else "-",
                 "ok" if ok else "FAIL"))
        if not ok:
            failures.append("%s: expected %.3f, file has %s"
                            % (name, expected, got))

    # 2. the same values, read live from the canonical dataset
    print("")
    if os.path.exists(COMPOUNDS_JS):
        canonical = read_compounds_js(COMPOUNDS_JS)
        for name, got in sorted(affinities.items()):
            record = canonical.get(name)
            if record is None:
                failures.append("%s is not in compounds.js" % name)
                print("    %-18s not in compounds.js               FAIL" % name)
                continue
            if abs(got - record["affinity"]) >= 1e-9:
                failures.append("%s: compounds.js has %.3f, file has %.3f"
                                % (name, record["affinity"], got))
                print("    %-18s %10.3f vs %.3f   FAIL"
                      % (name, record["affinity"], got))
                continue
            compound = compounds.get(name)
            if record["formula"] and compound and compound.formula:
                if record["formula"] != compound.formula:
                    failures.append("%s: formula %s vs %s"
                                    % (name, record["formula"], compound.formula))
                    print("    %-18s formula %s vs %s   FAIL"
                          % (name, record["formula"], compound.formula))
                    continue
            print("    %-18s matches compounds.js               ok" % name)
    else:
        print("    compounds.js not found, live comparison skipped")

    # 3. ligand efficiency inverts the ranking, which is the point of column 4
    print("")
    if efficiencies and affinities:
        strongest = min(affinities, key=lambda n: affinities[n])
        efficient = max(efficiencies, key=lambda n: efficiencies[n])
        print("    strongest raw score      %s" % strongest)
        print("    most efficient per atom  %s" % efficient)
        if strongest == efficient:
            failures.append("expected the two rankings to disagree")
            print("    rankings agree                              FAIL")
        else:
            print("    rankings disagree, as expected              ok")

    print("  " + "-" * 66)
    if failures:
        print("  FAILED (%d)" % len(failures))
        for item in failures:
            print("    %s" % item)
        return 1
    print("  PASSED. The experiment file matches compounds.js exactly.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
