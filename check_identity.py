#!/usr/bin/env python3
"""
Is the SMILES-built glucose the same molecule as the crystal glucose?

The chair-conformer preparation still missed the crystal pose (5.90 A) and
scored WORSE (-5.465) than both the twist-boat prep (-5.775) and the crystal
file (-5.905). A correctly puckered molecule that scores worse than a strained
one is a sign that the two are not the same chemical species.

The prime suspect for a pyranose is the ANOMER: alpha-D-glucopyranose has its
anomeric hydroxyl axial, beta has it equatorial. They are different molecules,
they cannot interconvert without breaking a bond, and Vina cannot convert one
into the other any more than it can flip a ring.

This compares canonical SMILES and InChI across every glucose file in play.
"""

import os
import sys

from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")
from meeko import PDBQTMolecule, RDKitMolCreate   # noqa: E402

DATA = (r"D:\THESIS_VSC\BALAKATDBV2\talanaidock\dist\TalanaiDock\app"
        r"\docking_assets\docking_data")
HERE = os.path.dirname(os.path.abspath(__file__))

FILES = [
    ("crystal reference (PDB)", os.path.join(DATA, "validation", "glc_3a4a_601.pdb"), "pdb"),
    ("crystal docking input", os.path.join(DATA, "validation", "glc_ligand.pdbqt"), "pdbqt"),
    ("crystal ref SDF", os.path.join(DATA, "validation", "glc_crystal_ref.sdf"), "sdf"),
    ("pipeline, single conf", os.path.join(HERE, "validation-inputs", "glc_reembedded.pdbqt"), "pdbqt"),
    ("pipeline, 50 confs", os.path.join(HERE, "validation-inputs", "glc_multiconf.pdbqt"), "pdbqt"),
]

# What we believe we asked for, PubChem CID 79025
ALPHA = "OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O"


def load(path, kind):
    if kind == "pdbqt":
        mols = RDKitMolCreate.from_pdbqt_mol(
            PDBQTMolecule.from_file(path, skip_typing=True))
        return Chem.RemoveHs(mols[0])
    if kind == "sdf":
        supplier = Chem.SDMolSupplier(path, removeHs=True, sanitize=False)
        return next(iter(supplier))
    return Chem.MolFromPDBFile(path, removeHs=True, sanitize=False)


def fingerprint(mol):
    """Canonical SMILES and InChI, with stereo assigned from 3D where possible."""
    try:
        work = Chem.Mol(mol)
        Chem.SanitizeMol(work)
        if work.GetNumConformers():
            Chem.AssignStereochemistryFrom3D(work)
        smiles = Chem.MolToSmiles(work)
        try:
            inchi = Chem.MolToInchi(work)
        except Exception:
            inchi = "(inchi failed)"
        return smiles, inchi
    except Exception as exc:
        return "(sanitize failed: %s)" % exc, "(n/a)"


def main():
    print("")
    print("  Chemical identity of every glucose in play")
    print("  Stereochemistry assigned from 3D coordinates where available.")
    print("  " + "-" * 72)

    reference = Chem.MolFromSmiles(ALPHA)
    Chem.AssignStereochemistry(reference, cleanIt=True, force=True)
    ref_smiles = Chem.MolToSmiles(reference)
    ref_inchi = Chem.MolToInchi(reference)
    print("  requested SMILES (CID 79025, alpha-D-glucopyranose)")
    print("    canonical %s" % ref_smiles)
    print("")

    rows = []
    for label, path, kind in FILES:
        if not os.path.isfile(path):
            print("  %-26s MISSING" % label)
            continue
        try:
            mol = load(path, kind)
        except Exception as exc:
            print("  %-26s could not load: %s" % (label, exc))
            continue
        if mol is None:
            print("  %-26s could not parse" % label)
            continue
        smiles, inchi = fingerprint(mol)
        rows.append((label, smiles, inchi))
        print("  %-26s %s" % (label, smiles))

    print("  " + "-" * 72)
    groups = {}
    for label, smiles, _inchi in rows:
        groups.setdefault(smiles, []).append(label)

    if len(groups) == 1:
        print("  All files are the same chemical species.")
        print("  The failure is NOT stereochemistry. Look elsewhere.")
        return 0

    print("  THEY ARE NOT ALL THE SAME MOLECULE. Distinct species found:")
    for index, (smiles, labels) in enumerate(groups.items(), start=1):
        print("")
        print("    species %d: %s" % (index, smiles))
        for label in labels:
            print("       %s" % label)
    print("")
    print("  A pyranose anomer cannot interconvert without breaking a bond, so")
    print("  Vina can no more fix this than it can flip a ring. If the crystal")
    print("  ligand and the prepared ligand are different anomers, the crystal")
    print("  pose was never reachable and no amount of search or conformer")
    print("  generation would have found it.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
