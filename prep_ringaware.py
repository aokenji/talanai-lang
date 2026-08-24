#!/usr/bin/env python3
"""
Ring-aware conformer selection, for ligands the energy criterion gets wrong.

THE PROBLEM THIS SOLVES
    prep_multiconf.py embeds 50 conformers and keeps the lowest MMFF energy
    one. That fixed ursolic acid (3 defective rings to 0) and isovitexin
    (1 to 0), but made ACARBOSE worse (1 defective ring to 2).

    Why: for a large flexible molecule, gas-phase MMFF energy is dominated by
    intramolecular contacts, internal hydrogen bonds and folding. A tightly
    folded conformer with a strained ring can score lower than an extended one
    with perfect chairs. Energy alone therefore selects against ring quality
    exactly where ring quality matters most.

    And ring geometry is the one thing docking cannot repair, because Vina
    holds rings rigid.

THE RULE HERE
    Among conformers whose saturated six-rings are ALL chairs, keep the lowest
    energy one. Ring correctness is a hard filter; energy only breaks ties.

    If no conformer qualifies, say so loudly and do not quietly ship the
    lowest-energy one as though it were fine.

    python prep_ringaware.py
"""

import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from rdkit import Chem, RDLogger                       # noqa: E402
from rdkit.Chem import AllChem                          # noqa: E402
RDLogger.DisableLog("rdApp.*")
from meeko import MoleculePreparation, PDBQTWriterLegacy   # noqa: E402

import audit_rings_rigorous as audit                    # noqa: E402

SEED = 0xF00D
MMFF_ITERS = 400
CHAIR_MAX_THETA = 35.0

OUT = os.path.join(HERE, "validation-inputs", "reprepared")

TARGETS = [
    {"name": "Acarbose", "stem": "acarbose", "confs": 300,
     "smiles": ("C[C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@@H]2[C@H](O[C@@H]"
                "([C@@H]([C@H]2O)O)O[C@H]([C@@H](CO)O)[C@@H]([C@H](C=O)O)O)CO)"
                "O)O)N[C@H]3C=C([C@H]([C@@H]([C@H]3O)O)O)CO")},
    {"name": "Ursolic acid", "stem": "ursolic", "confs": 150,
     "smiles": "CC1CCC2(C(=O)O)CCC3(C)C(=CCC4C5(C)CCC(O)C(C)(C)C5CCC43C)C2C1C"},
    {"name": "Isovitexin", "stem": "isovitexin", "confs": 150,
     "smiles": "O=c1cc(-c2ccc(O)cc2)oc2cc(O)c(C3OC(CO)C(O)C(O)C3O)c(O)c12"},
]


def saturated_ring_thetas(mol, conf_id):
    """Cremer-Pople theta for every SATURATED six-ring, in one conformer."""
    conf = mol.GetConformer(conf_id)
    out = []
    for ring in mol.GetRingInfo().AtomRings():
        if len(ring) != 6:
            continue
        if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
            continue
        if audit.ring_is_unsaturated(mol, list(ring)):
            continue          # cannot be a chair, not a defect
        order = audit.ordered_ring(mol, list(ring))
        if order is None:
            continue
        coords = [[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y,
                   conf.GetAtomPosition(i).z] for i in order]
        _Q, theta = audit.cremer_pople(coords)
        out.append(theta)
    return out


def all_chairs(thetas):
    return all(t <= CHAIR_MAX_THETA or t >= 180 - CHAIR_MAX_THETA
               for t in thetas)


def prepare_ring_aware(smiles, num_confs):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("RDKit could not parse SMILES")
    mol = Chem.AddHs(mol)

    params = AllChem.ETKDGv3()
    params.randomSeed = SEED
    params.maxIterations = 2000
    cids = list(AllChem.EmbedMultipleConfs(mol, numConfs=num_confs, params=params))
    if not cids:
        raise RuntimeError("embedding produced no conformers")
    results = AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=MMFF_ITERS)
    energies = [e for _converged, e in results]

    qualifying = []
    for index, cid in enumerate(cids):
        thetas = saturated_ring_thetas(mol, cid)
        if not thetas or all_chairs(thetas):
            qualifying.append(index)

    if qualifying:
        chosen = min(qualifying, key=lambda i: energies[i])
        rule = "lowest energy among all-chair conformers"
    else:
        chosen = min(range(len(energies)), key=lambda i: energies[i])
        rule = "NO ALL-CHAIR CONFORMER FOUND, fell back to lowest energy"

    keep = mol.GetConformer(cids[chosen])
    winner = Chem.Mol(mol)
    winner.RemoveAllConformers()
    winner.AddConformer(keep, assignId=True)

    setups = MoleculePreparation().prepare(winner)
    pdbqt, ok, err = PDBQTWriterLegacy.write_string(setups[0])
    if not ok:
        raise RuntimeError("Meeko PDBQT write failed: %s" % err)

    return {
        "pdbqt": pdbqt, "mol": winner, "chosen": chosen, "rule": rule,
        "n_conformers": len(cids), "n_qualifying": len(qualifying),
        "energy": energies[chosen], "best_energy": min(energies),
        "thetas": saturated_ring_thetas(winner, 0),
    }


def main():
    os.makedirs(OUT, exist_ok=True)
    print("")
    print("  Ring-aware conformer selection")
    print("  Rule: among conformers whose saturated six-rings are ALL chairs,")
    print("  keep the lowest energy. Ring correctness filters, energy breaks ties.")
    print("  " + "-" * 70)

    failures = 0
    for target in TARGETS:
        result = prepare_ring_aware(target["smiles"], target["confs"])
        path = os.path.join(OUT, target["stem"] + ".pdbqt")
        with open(path, "w", encoding="utf-8") as handle:
            handle.write(result["pdbqt"])

        thetas = result["thetas"]
        ok = all_chairs(thetas) if thetas else True
        if not ok:
            failures += 1
        print("  %-15s %3d of %3d conformers had all-chair rings"
              % (target["name"], result["n_qualifying"], result["n_conformers"]))
        print("  %-15s picked #%d, energy %.2f (global best %.2f)"
              % ("", result["chosen"] + 1, result["energy"], result["best_energy"]))
        print("  %-15s saturated ring theta: %s"
              % ("", ", ".join("%.1f" % t for t in thetas) or "none"))
        print("  %-15s %s  ->  %s"
              % ("", result["rule"], "ALL CHAIRS" if ok else "*** STILL DEFECTIVE ***"))
        print("")

    print("  " + "-" * 70)
    if failures:
        print("  %d ligand(s) still defective. More conformers, or a different"
              % failures)
        print("  generator, is needed. Do not dock these as they stand.")
        return 1
    print("  Every saturated ring is a chair. Written to validation-inputs/reprepared/.")
    print("  Note the energy cost: an all-chair conformer sits above the global")
    print("  minimum, which is the point. Vina cannot fix a ring; it can relax")
    print("  a slightly strained pose.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
