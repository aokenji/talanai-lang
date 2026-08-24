#!/usr/bin/env python3
"""
Which prepared ligands carry a non-chair saturated ring?

ESTABLISHED FOR GLUCOSE, 2026-08-03:
  _prep_assets_local.py embeds ONE conformer (AllChem.EmbedMolecule) and then
  runs a local MMFF94 minimisation. MMFF cannot cross a ring-flip barrier, so
  whichever ring conformation ETKDGv3 produced on that single attempt is frozen
  for good. Meeko then writes the ring as rigid, correctly, and AutoDock Vina
  cannot alter ring geometry at all. The docked glucose pose was verified to
  carry torsions identical to its input, to one decimal place.

  Consequence: glucose was prepared as a twist-boat (sign pattern ++-++-, mean
  |torsion| 38.9 deg) rather than the native chair (+-+-+-, 56.4 deg), and its
  crystal pose was therefore unreachable at any exhaustiveness. Confirmed
  empirically up to exhaustiveness 512.

THIS SCRIPT asks how far that reaches. It measures every non-aromatic
six-membered ring in every prepared ligand and flags the ones that are not
chairs. It runs no docking; it is pure geometry on files that already exist.

A chair alternates torsion signs (+-+-+-) with mean |torsion| near 55 deg.
Boats, twist-boats and half-chairs do not.
"""

import os
import sys

from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolTransforms as transforms
RDLogger.DisableLog("rdApp.*")
from meeko import PDBQTMolecule, RDKitMolCreate   # noqa: E402

DATA = (r"D:\THESIS_VSC\BALAKATDBV2\talanaidock\dist\TalanaiDock\app"
        r"\docking_assets\docking_data")
HERE = os.path.dirname(os.path.abspath(__file__))

# The ten screened compounds plus the reference. Spinosin has no original
# prepared file; the one in validation-inputs was built on 2026-08-03.
LIGANDS = [
    ("Rutin",          os.path.join(DATA, "ligands", "rutin.pdbqt")),
    ("Betulinic acid", os.path.join(DATA, "ligands", "betulinic.pdbqt")),
    ("Ursolic acid",   os.path.join(DATA, "ligands", "ursolic.pdbqt")),
    ("Isovitexin",     os.path.join(DATA, "ligands", "isovitexin.pdbqt")),
    ("Luteolin",       os.path.join(DATA, "ligands", "luteolin.pdbqt")),
    ("Quercetin",      os.path.join(DATA, "ligands", "quercetin.pdbqt")),
    ("Kaempferol",     os.path.join(DATA, "ligands", "kaempferol.pdbqt")),
    ("Vitexin",        os.path.join(DATA, "ligands", "vitexin.pdbqt")),
    ("Oleanolic acid", os.path.join(DATA, "ligands", "oleanolic.pdbqt")),
    ("ACARBOSE (ref)", os.path.join(DATA, "reference", "acarbose.pdbqt")),
    ("Spinosin (new)", os.path.join(HERE, "validation-inputs", "spinosin.pdbqt")),
    ("Isomaltose",     os.path.join(HERE, "validation-inputs", "isomaltose.pdbqt")),
]

CHAIR_MIN_MEAN = 45.0     # a real chair sits near 55 deg


def classify(torsions):
    signs = "".join("+" if t > 0 else "-" for t in torsions)
    mean_abs = sum(abs(t) for t in torsions) / 6.0
    alternating = signs in ("+-+-+-", "-+-+-+")
    if alternating and mean_abs >= CHAIR_MIN_MEAN:
        return "chair", signs, mean_abs
    if mean_abs < 20:
        return "near-planar", signs, mean_abs
    return "NOT A CHAIR", signs, mean_abs


def rings_of(path):
    mols = RDKitMolCreate.from_pdbqt_mol(
        PDBQTMolecule.from_file(path, skip_typing=True))
    mol = Chem.RemoveHs(mols[0])
    conf = mol.GetConformer(0)
    out = []
    for ring in mol.GetRingInfo().AtomRings():
        if len(ring) != 6:
            continue
        if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
            continue        # aromatic rings are planar by construction
        ring = list(ring)
        torsions = []
        for i in range(6):
            torsions.append(transforms.GetDihedralDeg(
                conf, ring[i], ring[(i + 1) % 6], ring[(i + 2) % 6],
                ring[(i + 3) % 6]))
        out.append(torsions)
    return out


def main():
    print("")
    print("  Non-aromatic six-ring conformations in the prepared ligands")
    print("  A chair alternates signs (+-+-+-) with mean |torsion| near 55 deg.")
    print("  Vina cannot change ring geometry, so whatever is here is what docked.")
    print("  " + "-" * 74)
    print("  %-16s %6s %-13s %8s  %s" % ("compound", "ring", "class", "mean|t|", "signs"))

    suspect = []
    for name, path in LIGANDS:
        if not os.path.isfile(path):
            print("  %-16s MISSING" % name)
            continue
        try:
            rings = rings_of(path)
        except Exception as exc:
            print("  %-16s could not read: %s" % (name, exc))
            continue
        if not rings:
            print("  %-16s %6s %-13s" % (name, "-", "no saturated 6-ring"))
            continue
        bad = 0
        for index, torsions in enumerate(rings, start=1):
            kind, signs, mean_abs = classify(torsions)
            if kind == "NOT A CHAIR":
                bad += 1
            print("  %-16s %6d %-13s %8.1f  %s"
                  % (name if index == 1 else "", index, kind, mean_abs, signs))
        if bad:
            suspect.append((name, bad, len(rings)))

    print("  " + "-" * 74)
    if not suspect:
        print("  Every saturated six-ring is a chair. Glucose was the exception.")
        return 0
    print("  NON-CHAIR RINGS FOUND, and Vina could not have fixed any of them:")
    for name, bad, total in suspect:
        print("    %-18s %d of %d saturated six-rings" % (name, bad, total))
    print("")
    print("  Root cause: _prep_assets_local.py line 41 embeds a SINGLE conformer")
    print("  (AllChem.EmbedMolecule) and minimises locally. MMFF cannot cross a")
    print("  ring-flip barrier. Fix: EmbedMultipleConfs, minimise all, keep the")
    print("  lowest-energy conformer.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
