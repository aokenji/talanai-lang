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

# The compounds this experiment file must cover. Their values are deliberately
# NOT listed here: they are read from compounds.js at run time, which is the
# whole point of this test. A literal table of affinities sitting beside the
# canonical dataset is a second copy of it, and the second copy is the one that
# goes stale.
#
# This file used to pin the May 2026 affinities (Rutin -8.857, Acarbose -6.660,
# and so on) and went on asserting them after three corrections had superseded
# them, while the docstring above claimed the numbers were never typed in here.
# The withdrawn values are preserved on purpose, in TalanaiHub's REPRODUCIBILITY
# record and in the Zenodo deposit. They do not belong in a live assertion.
REQUIRED = [
    "Rutin",
    "Betulinic Acid",
    "Quercetin",
    "Kaempferol",
    "Oleanolic Acid",
    "Acarbose",
]


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
    incomplete = False
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

    # 1. coverage: every required compound is present in the experiment file.
    #    Values are not checked here. Section 2 checks them, against the
    #    canonical dataset rather than against a copy of it.
    print("    %-18s %10s %9s  %s"
          % ("compound", "in .tal", "per atom", ""))
    for name in REQUIRED:
        got = affinities.get(name)
        efficiency = efficiencies.get(name)
        print("    %-18s %10s %9s  %s"
              % (name,
                 "%.3f" % got if got is not None else "MISSING",
                 "%.3f" % efficiency if efficiency else "-",
                 "ok" if got is not None else "FAIL"))
        if got is None:
            failures.append("%s is required but absent from the experiment file"
                            % name)

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
        incomplete = True
        print("    compounds.js is not reachable at that path.")
        print("    Coverage was checked. No affinity was verified.")

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
    if incomplete:
        # Not a pass. Nothing was compared against the canonical dataset, and a
        # test that verifies nothing must not print the word that means it did.
        print("  INCOMPLETE. Coverage passed, but the canonical dataset was out")
        print("  of reach, so no affinity was checked against it.")
        return 2
    print("  PASSED. The experiment file matches compounds.js exactly.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
