"""
Ring-conformation checking for R307.

Kept in its own module because it is the one part of Talanai that needs a
chemistry toolkit. The validator's core stays standard-library only and
`tal check` must keep working with nothing installed, so R307 degrades to
UNVERIFIED when RDKit and Meeko are absent rather than failing.

WHY THIS EXISTS
    AutoDock Vina holds ring bonds rigid. Whatever ring conformation a ligand
    is built with is the one that docks and scores. A saturated six-ring that
    is not a chair was built wrong, and no amount of search fixes it.

THE MEASURE
    Cremer and Pople, JACS 1975, 97, 1354. For a six-ring, theta near 0 or 180
    degrees is a chair, near 90 is the boat family. Unlike a torsion sign
    pattern, theta does not depend on which direction the ring is traversed.

    Rings containing an sp2 atom are skipped: a double bond forces a half chair
    or sofa, which is correct rather than defective.
"""

from __future__ import annotations

import math
import os

EXTENSIONS = (".pdbqt", ".sdf", ".mol2", ".pdb")


def find_ligand_file(directory, compound):
    """Locate a prepared file for a compound, tolerating naming conventions."""
    stems = [compound.filename] if getattr(compound, "filename", None) else []
    name = compound.name
    stems += [name, name.lower(), name.replace("_", ""), name.replace("_", "-"),
              name.split("_")[0], name.split("_")[0].lower()]
    for stem in stems:
        if not stem:
            continue
        if os.path.splitext(stem)[1]:
            candidate = os.path.join(directory, stem)
            if os.path.isfile(candidate):
                return candidate
            continue
        for extension in EXTENSIONS:
            candidate = os.path.join(directory, stem + extension)
            if os.path.isfile(candidate):
                return candidate
    return None


def _cremer_pople_theta(coords):
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
    theta = math.degrees(math.atan2(q2, q3)) % 360.0
    return 360.0 - theta if theta > 180.0 else theta


def saturated_ring_report(path, max_theta):
    """
    {"saturated": n, "defective": [theta, ...]} for one prepared ligand.

    Returns None when the file cannot be read or the toolkit is unavailable,
    so the caller can report UNVERIFIED rather than inventing a pass.
    """
    try:
        from rdkit import Chem, RDLogger
        RDLogger.DisableLog("rdApp.*")
    except ImportError:
        return None

    try:
        if path.lower().endswith(".pdbqt"):
            from meeko import PDBQTMolecule, RDKitMolCreate
            mol = Chem.RemoveHs(RDKitMolCreate.from_pdbqt_mol(
                PDBQTMolecule.from_file(path, skip_typing=True))[0])
        elif path.lower().endswith(".sdf"):
            mol = next(iter(Chem.SDMolSupplier(path, removeHs=True)))
        elif path.lower().endswith(".mol2"):
            mol = Chem.MolFromMol2File(path, removeHs=True)
        else:
            mol = Chem.MolFromPDBFile(path, removeHs=True)
    except Exception:
        return None
    if mol is None or not mol.GetNumConformers():
        return None

    conf = mol.GetConformer(0)
    saturated, defective = 0, []
    for ring in mol.GetRingInfo().AtomRings():
        if len(ring) != 6:
            continue
        if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
            continue
        members = set(ring)
        unsaturated = False
        for bond in mol.GetBonds():
            a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            if a in members and b in members and bond.GetBondType() != Chem.BondType.SINGLE:
                unsaturated = True
                break
        if not unsaturated:
            for index in members:
                atom = mol.GetAtomWithIdx(index)
                if atom.GetHybridization() == Chem.HybridizationType.SP2 or any(
                        bond.GetBondType() == Chem.BondType.DOUBLE
                        for bond in atom.GetBonds()):
                    unsaturated = True
                    break
        if unsaturated:
            continue          # cannot be a chair; not a defect

        order, current, previous = [ring[0]], ring[0], None
        for _step in range(5):
            neighbours = [n.GetIdx() for n in mol.GetAtomWithIdx(current).GetNeighbors()
                          if n.GetIdx() in members and n.GetIdx() != previous]
            if not neighbours:
                break
            previous, current = current, neighbours[0]
            order.append(current)
        if len(order) != 6:
            continue

        saturated += 1
        coords = [[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y,
                   conf.GetAtomPosition(i).z] for i in order]
        theta = _cremer_pople_theta(coords)
        if not (theta <= max_theta or theta >= 180.0 - max_theta):
            defective.append(theta)

    return {"saturated": saturated, "defective": defective}
