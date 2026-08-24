#!/usr/bin/env python3
"""
Ligand preparation with a conformer search, and the experiment that tests it.

WHAT IS BEING FIXED
    The study's preparation (_prep_assets_local.py) embeds ONE conformer with
    AllChem.EmbedMolecule and then runs a LOCAL MMFF94 minimisation. Local
    minimisation cannot cross a ring-flip barrier, Meeko writes rings rigid,
    and AutoDock Vina cannot change ring geometry. So whichever pucker ETKDG
    produced on that single attempt is frozen for the whole study.

    Verified consequence: glucose came out as a twist-boat rather than the
    native chair, and its crystal pose was unreachable at exhaustiveness 512.

THE FIX
    Embed many conformers, minimise all of them, keep the lowest-energy one.
    Everything else is held identical to the original: same ETKDGv3, same seed,
    same MMFF94, same iteration bound, same Meeko call.

    python prep_multiconf.py            prepare and compare glucose
    python prep_multiconf.py --smiles "..." --out foo.pdbqt --name foo
"""

import argparse
import os
import sys

SEED = 0xF00D          # identical to _prep_assets_local.py
NUM_CONFS = 50
MMFF_ITERS = 400

HERE = os.path.dirname(os.path.abspath(__file__))

# alpha-D-glucopyranose, PubChem CID 79025, as used for glc_reembedded.pdbqt
GLUCOSE = "OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O"


def prepare(smiles, num_confs=NUM_CONFS, single=False):
    """
    SMILES to (pdbqt_string, mol, chosen_conf_id, energies).

    single=True reproduces the study's original one-shot behaviour, so the two
    can be compared under otherwise identical settings.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from meeko import MoleculePreparation, PDBQTWriterLegacy

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("RDKit could not parse SMILES")
    mol = Chem.AddHs(mol)

    params = AllChem.ETKDGv3()
    params.randomSeed = SEED
    params.maxIterations = 2000

    if single:
        if AllChem.EmbedMolecule(mol, params) != 0:
            raise RuntimeError("ETKDGv3 embedding failed")
        AllChem.MMFFOptimizeMolecule(mol, maxIters=MMFF_ITERS)
        chosen, energies = 0, None
    else:
        cids = AllChem.EmbedMultipleConfs(mol, numConfs=num_confs, params=params)
        if not len(cids):
            raise RuntimeError("ETKDGv3 embedding failed for every conformer")
        results = AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=MMFF_ITERS)
        # results is [(converged, energy), ...] in conformer order
        energies = [e for _converged, e in results]
        chosen = min(range(len(energies)), key=lambda i: energies[i])
        # Keep only the winner, so the written file is unambiguous.
        keep = mol.GetConformer(cids[chosen])
        winner = Chem.Mol(mol)
        winner.RemoveAllConformers()
        winner.AddConformer(keep, assignId=True)
        mol = winner

    setups = MoleculePreparation().prepare(mol)
    pdbqt, ok, err = PDBQTWriterLegacy.write_string(setups[0])
    if not ok:
        raise RuntimeError("Meeko PDBQT write failed: %s" % err)
    return pdbqt, mol, chosen, energies


def ring_signature(mol):
    """(signs, mean_abs_torsion) for each non-aromatic six-ring."""
    from rdkit.Chem import rdMolTransforms as transforms
    conf = mol.GetConformer(0)
    out = []
    for ring in mol.GetRingInfo().AtomRings():
        if len(ring) != 6:
            continue
        if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
            continue
        ring = list(ring)
        torsions = []
        for i in range(6):
            torsions.append(transforms.GetDihedralDeg(
                conf, ring[i], ring[(i + 1) % 6], ring[(i + 2) % 6],
                ring[(i + 3) % 6]))
        signs = "".join("+" if t > 0 else "-" for t in torsions)
        out.append((signs, sum(abs(t) for t in torsions) / 6.0))
    return out


def describe(label, mol, energies=None, chosen=None):
    print("  %-28s" % label, end="")
    if energies is not None:
        print(" chose conformer %d of %d, energy %.2f (worst %.2f)"
              % (chosen + 1, len(energies), min(energies), max(energies)))
    else:
        print(" single embedding, no alternatives generated")
    for index, (signs, mean_abs) in enumerate(ring_signature(mol), start=1):
        chair = signs in ("+-+-+-", "-+-+-+") and mean_abs >= 45.0
        print("  %-28s ring %d  %-11s mean|t| %5.1f  signs %s"
              % ("", index, "chair" if chair else "NOT A CHAIR", mean_abs, signs))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--smiles", default=GLUCOSE)
    ap.add_argument("--name", default="glc_multiconf")
    ap.add_argument("--out", default=None)
    ap.add_argument("--confs", type=int, default=NUM_CONFS)
    args = ap.parse_args()

    out_dir = os.path.join(HERE, "validation-inputs")
    os.makedirs(out_dir, exist_ok=True)
    out_path = args.out or os.path.join(out_dir, args.name + ".pdbqt")

    print("")
    print("  Ligand preparation: one conformer versus %d" % args.confs)
    print("  ETKDGv3 seed 0x%X, MMFF94 %d iters, Meeko. Only the conformer"
          % (SEED, MMFF_ITERS))
    print("  search differs.")
    print("  " + "-" * 70)

    _pdbqt_one, mol_one, _c, _e = prepare(args.smiles, single=True)
    describe("original, EmbedMolecule", mol_one)
    print("")

    pdbqt, mol_many, chosen, energies = prepare(args.smiles, args.confs)
    describe("fixed, EmbedMultipleConfs", mol_many, energies, chosen)

    with open(out_path, "w", encoding="utf-8") as handle:
        handle.write(pdbqt)
    print("")
    print("  wrote %s" % out_path)

    spread = max(energies) - min(energies)
    print("  conformer energy spread %.2f kcal/mol across %d embeddings"
          % (spread, len(energies)))
    print("")
    return 0


if __name__ == "__main__":
    sys.exit(main())
