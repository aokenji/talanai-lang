#!/usr/bin/env python3
"""
Which chair is it? Cremer-Pople puckering parameters for the pyranose ring.

WHY THE EARLIER TEST WAS NOT ENOUGH
    Classifying a ring by its torsion sign pattern (+-+-+- versus -+-+-+)
    depends on which direction you traverse the ring and where you start, so it
    cannot distinguish the two chairs. It can only separate "chair-like" from
    "not chair-like". Both 4C1 and 1C4 are chairs.

    That matters here, because 4C1 and 1C4 are related by a ring flip that
    AutoDock Vina cannot perform. If the pipeline built 1C4 while the crystal
    ligand is 4C1, the crystal pose was unreachable no matter the search effort
    and no matter that the pucker amplitude looked correct.

THE MEASURE
    Cremer and Pople, JACS 1975, 97, 1354. For a six-ring, theta near 0 deg is
    4C1 and theta near 180 deg is 1C4, with Q the puckering amplitude. Both are
    independent of traversal direction ONCE the ring is ordered canonically,
    which for a pyranose means starting at the ring oxygen (O5) and walking
    O5, C1, C2, C3, C4, C5.
"""

import math
import os
import sys

from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")
from meeko import PDBQTMolecule, RDKitMolCreate   # noqa: E402

DATA = (r"D:\THESIS_VSC\BALAKATDBV2\talanaidock\dist\TalanaiDock\app"
        r"\docking_assets\docking_data")
HERE = os.path.dirname(os.path.abspath(__file__))

FILES = [
    ("crystal reference", os.path.join(DATA, "validation", "glc_3a4a_601.pdb"), "pdb"),
    ("crystal docking input", os.path.join(DATA, "validation", "glc_ligand.pdbqt"), "pdbqt"),
    ("pipeline, single conf", os.path.join(HERE, "validation-inputs", "glc_reembedded.pdbqt"), "pdbqt"),
    ("pipeline, 50 confs", os.path.join(HERE, "validation-inputs", "glc_multiconf.pdbqt"), "pdbqt"),
    ("docked, single conf", os.path.join(HERE, "validation-run", "orientation", "pipeline_e512_s42.pdbqt"), "pdbqt"),
    ("docked, 50 confs", os.path.join(HERE, "validation-run", "fixtest", "glc_multiconf_e32_s42.pdbqt"), "pdbqt"),
]


def load(path, kind):
    if kind == "pdbqt":
        return Chem.RemoveHs(RDKitMolCreate.from_pdbqt_mol(
            PDBQTMolecule.from_file(path, skip_typing=True))[0])
    return Chem.MolFromPDBFile(path, removeHs=True, sanitize=False)


def ordered_pyranose_ring(mol):
    """Ring atom indices ordered O5, C1, C2, C3, C4, C5."""
    for ring in mol.GetRingInfo().AtomRings():
        if len(ring) != 6:
            continue
        oxygens = [i for i in ring if mol.GetAtomWithIdx(i).GetSymbol() == "O"]
        if len(oxygens) != 1:
            continue
        start = oxygens[0]
        members = set(ring)
        order = [start]
        current, previous = start, None
        for _step in range(5):
            neighbours = [n.GetIdx() for n in mol.GetAtomWithIdx(current).GetNeighbors()
                          if n.GetIdx() in members and n.GetIdx() != previous]
            if not neighbours:
                return None
            previous, current = current, neighbours[0]
            order.append(current)
        # Walk toward C1, the ring carbon bearing an exocyclic oxygen (anomeric).
        anomeric = [i for i in (order[1], order[-1])
                    if any(n.GetSymbol() == "O" and n.GetIdx() not in members
                           for n in mol.GetAtomWithIdx(i).GetNeighbors())]
        if anomeric and order[1] not in anomeric:
            order = [order[0]] + order[1:][::-1]
        return order
    return None


def cremer_pople(coords):
    """Q, theta, phi for a six-membered ring given canonically ordered coords."""
    n = 6
    centre = [sum(c[k] for c in coords) / n for k in range(3)]
    rel = [[c[k] - centre[k] for k in range(3)] for c in coords]

    r1 = [sum(rel[j][k] * math.sin(2 * math.pi * j / n) for j in range(n)) for k in range(3)]
    r2 = [sum(rel[j][k] * math.cos(2 * math.pi * j / n) for j in range(n)) for k in range(3)]
    normal = [r1[1] * r2[2] - r1[2] * r2[1],
              r1[2] * r2[0] - r1[0] * r2[2],
              r1[0] * r2[1] - r1[1] * r2[0]]
    length = math.sqrt(sum(v * v for v in normal))
    unit = [v / length for v in normal]

    z = [sum(rel[j][k] * unit[k] for k in range(3)) for j in range(n)]

    q2cos = math.sqrt(2.0 / n) * sum(z[j] * math.cos(4 * math.pi * j / n) for j in range(n))
    q2sin = -math.sqrt(2.0 / n) * sum(z[j] * math.sin(4 * math.pi * j / n) for j in range(n))
    q3 = (1.0 / math.sqrt(n)) * sum(((-1) ** j) * z[j] for j in range(n))
    q2 = math.sqrt(q2cos ** 2 + q2sin ** 2)

    Q = math.sqrt(q2 ** 2 + q3 ** 2)
    theta = math.degrees(math.atan2(q2, q3)) % 360.0
    phi = math.degrees(math.atan2(q2sin, q2cos)) % 360.0
    return Q, theta, phi


def label_of(theta):
    if theta <= 45:
        return "4C1 chair"
    if theta >= 135:
        return "1C4 chair"
    return "boat / twist-boat"


def main():
    print("")
    print("  Cremer-Pople puckering, ring ordered O5 C1 C2 C3 C4 C5")
    print("  theta near 0 = 4C1 chair, near 180 = 1C4 chair, mid = boat family")
    print("  Vina cannot convert one chair into the other.")
    print("  " + "-" * 70)
    print("  %-24s %7s %8s %8s   %s" % ("file", "Q", "theta", "phi", "conformation"))

    seen = []
    for label, path, kind in FILES:
        if not os.path.isfile(path):
            print("  %-24s MISSING" % label)
            continue
        try:
            mol = load(path, kind)
            order = ordered_pyranose_ring(mol)
            if order is None:
                print("  %-24s no pyranose ring found" % label)
                continue
            conf = mol.GetConformer(0)
            coords = [[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y,
                       conf.GetAtomPosition(i).z] for i in order]
            Q, theta, phi = cremer_pople(coords)
            print("  %-24s %7.3f %8.1f %8.1f   %s"
                  % (label, Q, theta, phi, label_of(theta)))
            seen.append((label, label_of(theta)))
        except Exception as exc:
            print("  %-24s failed: %s: %s" % (label, type(exc).__name__, exc))

    print("  " + "-" * 70)
    if len(seen) >= 2:
        crystal = seen[0][1]
        mismatched = [l for l, k in seen[1:] if k != crystal]
        if mismatched:
            print("  The crystal ring is a %s. These differ:" % crystal)
            for label in mismatched:
                print("    %s" % label)
            print("")
            print("  A ring flip is not something docking can do. Any file above")
            print("  in the wrong chair could never have reached the crystal pose.")
        else:
            print("  Every ring is the same conformation as the crystal (%s)." % crystal)
            print("  Ring pucker is therefore NOT the reason the redock failed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
