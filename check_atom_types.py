#!/usr/bin/env python3
"""
What does Vina actually see differently between the two glucose files?

RULED OUT SO FAR
    chemical identity   identical canonical SMILES
    ring conformation   both 4C1 chairs (Cremer-Pople theta 13.5 vs 2.1)
    search budget       identical result at exhaustiveness 32, 128 and 512

REMAINING
    Vina reads AutoDock atom types from columns 78-79 of the PDBQT and nothing
    else about chemistry. Its hydrogen-bond term depends entirely on those
    types: OA and NA accept, HD donates, C and A do not participate. If the two
    files carry different type counts, or a different torsion tree, then they
    are different molecules AS FAR AS VINA IS CONCERNED, regardless of being
    the same molecule to a chemist.

    That would explain a different global optimum with the same everything else.
"""

import os
import sys
from collections import Counter

DATA = (r"D:\THESIS_VSC\BALAKATDBV2\talanaidock\dist\TalanaiDock\app"
        r"\docking_assets\docking_data")
HERE = os.path.dirname(os.path.abspath(__file__))

FILES = [
    ("crystal input", os.path.join(DATA, "validation", "glc_ligand.pdbqt")),
    ("pipeline single conf", os.path.join(HERE, "validation-inputs", "glc_reembedded.pdbqt")),
    ("pipeline 50 confs", os.path.join(HERE, "validation-inputs", "glc_multiconf.pdbqt")),
]


def parse(path):
    types = Counter()
    branches = 0
    torsdof = None
    atoms = 0
    with open(path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(("ATOM", "HETATM")):
                atoms += 1
                types[line[77:79].strip()] += 1
            elif line.startswith("BRANCH"):
                branches += 1
            elif line.startswith("TORSDOF"):
                torsdof = line.split()[1]
    return types, branches, torsdof, atoms


def main():
    print("")
    print("  What AutoDock Vina sees in each glucose file")
    print("  Vina's hydrogen-bond term depends only on these atom types:")
    print("  OA and NA accept, HD donates, C and A do not participate.")
    print("  " + "-" * 68)

    parsed = []
    for label, path in FILES:
        if not os.path.isfile(path):
            print("  %-22s MISSING" % label)
            continue
        types, branches, torsdof, atoms = parse(path)
        parsed.append((label, types, branches, torsdof, atoms))
        summary = "  ".join("%s=%d" % (k, v) for k, v in sorted(types.items()))
        print("  %-22s atoms %2d  BRANCH %d  TORSDOF %s"
              % (label, atoms, branches, torsdof))
        print("  %-22s types  %s" % ("", summary))
        print("")

    print("  " + "-" * 68)
    if len(parsed) < 2:
        return 0

    ref_label, ref_types, ref_branch, ref_tors, _a = parsed[0]
    differences = False
    for label, types, branches, torsdof, _atoms in parsed[1:]:
        problems = []
        if types != ref_types:
            keys = set(types) | set(ref_types)
            for key in sorted(keys):
                if types.get(key, 0) != ref_types.get(key, 0):
                    problems.append("%s: %d here vs %d in the crystal file"
                                    % (key, types.get(key, 0), ref_types.get(key, 0)))
        if branches != ref_branch:
            problems.append("BRANCH count %d vs %d" % (branches, ref_branch))
        if torsdof != ref_tors:
            problems.append("TORSDOF %s vs %s" % (torsdof, ref_tors))
        if problems:
            differences = True
            print("  %s DIFFERS from %s:" % (label, ref_label))
            for problem in problems:
                print("      %s" % problem)

    if not differences:
        print("  Every file carries identical atom types and the same torsion")
        print("  tree. Vina sees the same molecule in all of them, so the")
        print("  different global optimum is NOT explained by typing.")
        print("")
        print("  What remains: the internal bond lengths and angles. Crystal")
        print("  geometry comes from experimental refinement, the prepared ones")
        print("  from MMFF94 idealisation. Vina holds those rigid too.")
    else:
        print("")
        print("  Vina scores hydrogen bonds purely from these types. A file with")
        print("  fewer HD or OA atoms cannot form the same hydrogen bonds, so it")
        print("  will find a different optimum however hard it searches.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
