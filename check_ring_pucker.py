#!/usr/bin/env python3
"""
Why can the pipeline-built glucose never reach its crystal pose?

The orientation test showed rank-1 RMSD pinned at 5.35 to 5.36 A across
exhaustiveness 128 and 512 and three seeds, with scores varying by 0.002
kcal/mol. That is a converged search, not an undersampled one, so the crystal
pose is probably not reachable by this ligand at all.

For a pyranose there is one dominant candidate: RING PUCKER. AutoDock Vina
treats ring bonds as rigid, so whatever chair conformation the ligand
preparation generates is the only one the search can ever use. If the
generated pucker differs from the bound pucker, the crystal pose is
unreachable at any exhaustiveness.

This compares the ring geometry of the crystal-derived ligand against the
pipeline-built one directly.
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

TARGETS = [
    ("crystal input  glc_ligand.pdbqt",
     os.path.join(DATA, "validation", "glc_ligand.pdbqt")),
    ("pipeline built glc_reembedded.pdbqt",
     os.path.join(HERE, "validation-inputs", "glc_reembedded.pdbqt")),
    ("pipeline docked pose e512 s42",
     os.path.join(HERE, "validation-run", "orientation", "pipeline_e512_s42.pdbqt")),
]


def load(path):
    mols = RDKitMolCreate.from_pdbqt_mol(
        PDBQTMolecule.from_file(path, skip_typing=True))
    return Chem.RemoveHs(mols[0])


def ring_report(label, mol):
    rings = [r for r in mol.GetRingInfo().AtomRings() if len(r) == 6]
    if not rings:
        print("  %-38s no six-membered ring found" % label)
        return None
    ring = list(rings[0])
    conf = mol.GetConformer(0)
    symbols = "".join(mol.GetAtomWithIdx(i).GetSymbol() for i in ring)

    torsions = []
    for i in range(6):
        a = ring[i]
        b = ring[(i + 1) % 6]
        c = ring[(i + 2) % 6]
        d = ring[(i + 3) % 6]
        torsions.append(transforms.GetDihedralDeg(conf, a, b, c, d))

    mean_abs = sum(abs(t) for t in torsions) / 6.0
    signs = "".join("+" if t > 0 else "-" for t in torsions)
    print("  %-38s ring %s" % (label, symbols))
    print("  %-38s torsions %s" % ("", " ".join("%7.1f" % t for t in torsions)))
    print("  %-38s mean|t| %5.1f   signs %s" % ("", mean_abs, signs))
    return signs, torsions


def main():
    print("")
    print("  Pyranose ring pucker comparison")
    print("  A chair alternates torsion signs (+-+-+-) with mean |torsion|")
    print("  near 55 deg. The two chairs, 4C1 and 1C4, are mirror sign")
    print("  patterns of each other. Vina cannot convert one into the other,")
    print("  because ring bonds are not rotatable.")
    print("  " + "-" * 72)

    seen = []
    for label, path in TARGETS:
        if not os.path.isfile(path):
            print("  %-38s MISSING" % label)
            continue
        result = ring_report(label, load(path))
        if result:
            seen.append((label, result[0]))
        print("")

    print("  " + "-" * 72)
    if len(seen) >= 2:
        reference = seen[0][1]
        for label, signs in seen[1:]:
            same = signs == reference
            print("  %-38s %s the crystal ring conformation"
                  % (label, "MATCHES" if same else "DIFFERS from"))
        if any(s != reference for _l, s in seen[1:]):
            print("")
            print("  If the sign pattern is inverted, the pipeline generated the")
            print("  opposite chair. The crystal pose is then unreachable at any")
            print("  exhaustiveness, and every cyclic ligand in the study is")
            print("  exposed to the same failure mode.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
