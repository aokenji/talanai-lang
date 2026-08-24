#!/usr/bin/env python3
"""
Which prepared ligands carry a genuinely defective ring?

WHY THE FIRST AUDIT WAS NOT GOOD ENOUGH
    audit_ring_conformations.py classified rings by torsion sign pattern and
    mean magnitude. That flags any ring that is not a clean chair, including
    rings that CANNOT be chairs: a ring containing an sp2 carbon is a half
    chair or a sofa by necessity, not by defect. Oleanolic and ursolic acid
    both carry a delta-12,13 double bond, and acarbose carries a cyclohexene,
    so several of those flags were probably false positives.

    Sign patterns also depend on which direction the ring is traversed, so they
    cannot separate the two chairs.

WHAT THIS DOES INSTEAD
    For every non-aromatic six-ring:
      1. Determine whether the ring contains an sp2 atom, from bond orders.
      2. Compute Cremer-Pople Q and theta (JACS 1975, 97, 1354), which are
         geometric and do not depend on traversal direction for the
         chair-versus-boat question.
      3. Judge against what the ring is CAPABLE of:
           saturated ring   -> must be a chair (theta near 0 or 180)
           sp2-containing   -> half chair or sofa expected, not a defect

    Only a saturated ring that is not a chair is called a defect.

    python audit_rings_rigorous.py
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

LIGANDS = [
    ("Rutin",          os.path.join(DATA, "ligands", "rutin.pdbqt")),
    ("Betulinic acid", os.path.join(DATA, "ligands", "betulinic.pdbqt")),
    ("Ursolic acid",   os.path.join(DATA, "ligands", "ursolic.pdbqt")),
    ("Isovitexin",     os.path.join(DATA, "ligands", "isovitexin.pdbqt")),
    ("Vitexin",        os.path.join(DATA, "ligands", "vitexin.pdbqt")),
    ("Oleanolic acid", os.path.join(DATA, "ligands", "oleanolic.pdbqt")),
    ("Luteolin",       os.path.join(DATA, "ligands", "luteolin.pdbqt")),
    ("Quercetin",      os.path.join(DATA, "ligands", "quercetin.pdbqt")),
    ("Kaempferol",     os.path.join(DATA, "ligands", "kaempferol.pdbqt")),
    ("ACARBOSE (ref)", os.path.join(DATA, "reference", "acarbose.pdbqt")),
    ("Spinosin",       os.path.join(HERE, "validation-inputs", "spinosin.pdbqt")),
    ("Isomaltose",     os.path.join(HERE, "validation-inputs", "isomaltose.pdbqt")),
]

CHAIR_MAX_THETA = 35.0        # theta below this, or above 180-this, is a chair
MIN_AMPLITUDE = 0.30          # below this the ring is essentially flat


def cremer_pople(coords):
    n = 6
    centre = [sum(c[k] for c in coords) / n for k in range(3)]
    rel = [[c[k] - centre[k] for k in range(3)] for c in coords]
    r1 = [sum(rel[j][k] * math.sin(2 * math.pi * j / n) for j in range(n)) for k in range(3)]
    r2 = [sum(rel[j][k] * math.cos(2 * math.pi * j / n) for j in range(n)) for k in range(3)]
    normal = [r1[1] * r2[2] - r1[2] * r2[1],
              r1[2] * r2[0] - r1[0] * r2[2],
              r1[0] * r2[1] - r1[1] * r2[0]]
    length = math.sqrt(sum(v * v for v in normal)) or 1e-9
    unit = [v / length for v in normal]
    z = [sum(rel[j][k] * unit[k] for k in range(3)) for j in range(n)]
    q2cos = math.sqrt(2.0 / n) * sum(z[j] * math.cos(4 * math.pi * j / n) for j in range(n))
    q2sin = -math.sqrt(2.0 / n) * sum(z[j] * math.sin(4 * math.pi * j / n) for j in range(n))
    q3 = (1.0 / math.sqrt(n)) * sum(((-1) ** j) * z[j] for j in range(n))
    q2 = math.sqrt(q2cos ** 2 + q2sin ** 2)
    Q = math.sqrt(q2 ** 2 + q3 ** 2)
    theta = math.degrees(math.atan2(q2, q3)) % 360.0
    if theta > 180.0:
        theta = 360.0 - theta
    return Q, theta


def ring_is_unsaturated(mol, ring):
    """True if any bond inside the ring is not single, or any atom is sp2."""
    members = set(ring)
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if a in members and b in members:
            if bond.GetBondType() != Chem.BondType.SINGLE:
                return True
    for index in ring:
        atom = mol.GetAtomWithIdx(index)
        if atom.GetHybridization() == Chem.HybridizationType.SP2:
            return True
        # An exocyclic carbonyl also flattens the ring atom.
        for bond in atom.GetBonds():
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                return True
    return False


def ordered_ring(mol, ring):
    """Walk the ring in connectivity order, starting at a heteroatom if present."""
    members = set(ring)
    hetero = [i for i in ring if mol.GetAtomWithIdx(i).GetSymbol() != "C"]
    start = hetero[0] if hetero else ring[0]
    order, current, previous = [start], start, None
    for _ in range(5):
        neighbours = [n.GetIdx() for n in mol.GetAtomWithIdx(current).GetNeighbors()
                      if n.GetIdx() in members and n.GetIdx() != previous]
        if not neighbours:
            return None
        previous, current = current, neighbours[0]
        order.append(current)
    return order if len(order) == 6 else None


def main():
    print("")
    print("  Ring conformations, judged against what each ring CAN adopt")
    print("  Cremer-Pople theta: near 0 or 180 = chair, near 90 = boat family.")
    print("  A ring with an sp2 atom cannot be a chair, so it is not a defect.")
    print("  " + "-" * 76)
    print("  %-16s %4s %-11s %6s %7s  %s"
          % ("compound", "ring", "saturation", "Q", "theta", "verdict"))

    defects = {}
    for name, path in LIGANDS:
        if not os.path.isfile(path):
            print("  %-16s MISSING" % name)
            continue
        try:
            mol = Chem.RemoveHs(RDKitMolCreate.from_pdbqt_mol(
                PDBQTMolecule.from_file(path, skip_typing=True))[0])
        except Exception as exc:
            print("  %-16s unreadable: %s" % (name, exc))
            continue

        conf = mol.GetConformer(0)
        printed = False
        for index, ring in enumerate(mol.GetRingInfo().AtomRings(), start=1):
            if len(ring) != 6:
                continue
            if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
                continue
            order = ordered_ring(mol, list(ring))
            if order is None:
                continue
            coords = [[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y,
                       conf.GetAtomPosition(i).z] for i in order]
            Q, theta = cremer_pople(coords)
            unsaturated = ring_is_unsaturated(mol, list(ring))
            is_chair = (theta <= CHAIR_MAX_THETA or theta >= 180 - CHAIR_MAX_THETA)

            if unsaturated:
                verdict = "expected, sp2 in ring"
            elif Q < MIN_AMPLITUDE:
                verdict = "flat, unusual"
            elif is_chair:
                verdict = "chair, fine"
            else:
                verdict = "*** DEFECT ***"
                defects.setdefault(name, 0)
                defects[name] += 1

            print("  %-16s %4d %-11s %6.3f %7.1f  %s"
                  % (name if not printed else "", index,
                     "unsaturated" if unsaturated else "saturated",
                     Q, theta, verdict))
            printed = True
        if not printed:
            print("  %-16s    - %-11s %6s %7s  aromatic only, cannot be affected"
                  % (name, "n/a", "-", "-"))

    print("  " + "-" * 76)
    if not defects:
        print("  No saturated ring is misfolded. The first audit's extra flags")
        print("  were sp2-containing rings that cannot be chairs.")
        return 0
    print("  GENUINE DEFECTS, saturated rings that are not chairs:")
    for name, count in defects.items():
        print("    %-18s %d ring(s)" % (name, count))
    print("")
    print("  Vina holds rings rigid, so none of these could have been corrected")
    print("  during docking. Root cause: _prep_assets_local.py:41 embeds a")
    print("  single conformer and minimises locally.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
