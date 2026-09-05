#!/usr/bin/env python3
"""
The deprotonated arm of the enrichment benchmark: betulinic, oleanolic and
ursolic acid as the C-28 carboxylate, which is what the thesis actually
screened.

WHY THIS EXISTS
    The first enrichment run took its actives' SMILES from PubChem, which
    returns the NEUTRAL molecule. So the three triterpene acids were docked as
    COOH, the pre-correction form the 2026-08-13 correction explicitly
    abandoned in favour of COO- "to match the form that predominates (>99%) at
    both physiological pH 7.4 and the ~6.8 used by this project's own
    literature IC50 assays".

    Three of nine actives were therefore not the molecules the thesis
    screened, and they were the three that performed worst: pooled ROC AUC
    0.421, but 0.588 for the six correctly-protonated flavonoids against
    0.082 for these three. The whole negative signal came from an error in
    the benchmark, not from the protocol.

WHY CHARGE MATCHING IS DROPPED HERE, AND WHY THAT IS DEFENSIBLE
    make_decoys.py requires an EXACT formal-charge match. Applied to the
    deprotonated triterpenes it is unusable: a scan of the whole 1.7 million
    molecule ZINC slice found exactly ONE candidate at charge -1 inside the
    property window. Large lipophilic carboxylates barely exist in
    purchasable screening space.

    Dropping it is defensible for THIS engine specifically. AutoDock Vina has
    no electrostatic term: gauss, repulsion, hydrophobic and hydrogen bonding,
    none of which reads a formal charge. What deprotonation actually changes
    in Vina's world is the hydrogen-bond donor count (2 -> 1, since the
    carboxyl no longer donates) and the lipophilicity (logP 7.1 -> 5.8), and
    BOTH of those are matched here, tightly.

    So the criterion being relaxed is one the scoring function cannot see, and
    the criteria that replace it are the ones it can. That is a different act
    from widening a window until enough decoys appear, and the distinction is
    the whole point: state it, do not bury it.

    It is still a relaxation. It is recorded in the manifest, and any report
    of this arm has to carry it.

    python make_decoys_deprotonated.py
"""

import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from rdkit import Chem, RDLogger                          # noqa: E402
from rdkit.Chem import DataStructs                        # noqa: E402
RDLogger.DisableLog("rdApp.*")

import make_decoys                                         # noqa: E402

OUT = os.path.join(HERE, "validation-inputs", "enrichment")
SOURCE = os.path.join(OUT, "decoys.json")
LIBRARY = os.path.join(OUT, "zinc_library.smi")
TARGET = os.path.join(OUT, "decoys_deprotonated.json")

TRITERPENES = ("Betulinic Acid", "Oleanolic Acid", "Ursolic Acid")
PER_ACTIVE = 20
CARBOXYLIC_ACID = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")

# Same windows as the neutral arm, minus charge. See the module docstring.
WINDOWS = {k: v for k, v in make_decoys.WINDOWS.items() if k != "charge"}


def deprotonate(smiles, name):
    """C-28 COOH -> COO-. Exactly one carboxylic acid is expected."""
    mol = Chem.MolFromSmiles(smiles)
    hits = mol.GetSubstructMatches(CARBOXYLIC_ACID)
    if len(hits) != 1:
        raise SystemExit("%s has %d carboxylic acids, expected exactly 1"
                         % (name, len(hits)))
    editable = Chem.RWMol(mol)
    oxygen = hits[0][2]
    editable.GetAtomWithIdx(oxygen).SetFormalCharge(-1)
    editable.GetAtomWithIdx(oxygen).SetNumExplicitHs(0)
    out = editable.GetMol()
    Chem.SanitizeMol(out)
    if Chem.GetFormalCharge(out) != -1:
        raise SystemExit("%s did not come out at charge -1" % name)
    return Chem.MolToSmiles(out)


def within(candidate, target):
    return all(abs(candidate[k] - target[k]) <= tol for k, tol in WINDOWS.items())


def main():
    with open(SOURCE, encoding="utf-8") as handle:
        source = json.load(handle)

    print("")
    print("  STEP 1  deprotonate the C-28 carboxylic acid")
    print("  " + "-" * 74)
    actives = []
    for entry in source["actives"]:
        if entry["name"] not in TRITERPENES:
            continue
        smiles = deprotonate(entry["smiles"], entry["name"])
        mol = Chem.MolFromSmiles(smiles)
        props = make_decoys.props(mol)
        actives.append({**entry, "smiles": smiles, "props": props,
                        "neutral_smiles": entry["smiles"],
                        "protonation": "C-28 carboxylate (COO-), as the thesis screened"})
        print("    %-16s charge %+d  logP %5.2f  HBD %d  (neutral was logP %5.2f, HBD %d)"
              % (entry["name"], props["charge"], props["logp"], props["hbd"],
                 entry["props"]["logp"], entry["props"]["hbd"]))
    if len(actives) != len(TRITERPENES):
        raise SystemExit("expected %d triterpenes, found %d"
                         % (len(TRITERPENES), len(actives)))

    fingerprints = [make_decoys.fingerprint(Chem.MolFromSmiles(a["smiles"]))
                    for a in actives]

    print("")
    print("  STEP 2  matching, charge criterion dropped (see the docstring)")
    print("  " + "-" * 74)
    print("    windows: " + ", ".join("%s +-%s" % (k, v) for k, v in WINDOWS.items()))
    print("    rejecting Tanimoto > %.2f to any active" % make_decoys.MAX_SIMILARITY)

    # Do not reuse a molecule that is already a decoy in the neutral arm: the
    # two arms should not share decoys, or their results are not independent.
    used = {d["smiles"] for group in source["decoys"].values() for d in group}
    picked = {a["name"]: [] for a in actives}
    scanned = 0

    for smiles, ident in make_decoys.read_library(LIBRARY):
        scanned += 1
        if all(len(v) >= PER_ACTIVE for v in picked.values()):
            break
        short = [a for a in actives if len(picked[a["name"]]) < PER_ACTIVE]
        if not short or smiles in used:
            continue
        mol = Chem.MolFromSmiles(smiles)
        if mol is None or abs(mol.GetNumHeavyAtoms() - 33) > WINDOWS["heavy"]:
            continue
        candidate = make_decoys.props(mol)
        target = next((a for a in short if within(candidate, a["props"])), None)
        if target is None:
            continue
        fp = make_decoys.fingerprint(mol)
        if max(DataStructs.BulkTanimotoSimilarity(fp, fingerprints)) > make_decoys.MAX_SIMILARITY:
            continue
        used.add(smiles)
        picked[target["name"]].append({"smiles": smiles, "library_id": ident,
                                       "props": candidate})
        if scanned % 100000 == 0:
            print("    scanned %7d   matched %3d / %d"
                  % (scanned, sum(len(v) for v in picked.values()),
                     PER_ACTIVE * len(actives)))

    print("")
    print("  STEP 3  result")
    print("  " + "-" * 74)
    for a in actives:
        got = len(picked[a["name"]])
        print("    %-16s %3d decoys%s" % (a["name"], got,
                                          "" if got >= PER_ACTIVE else "   <-- SHORT"))

    manifest = {
        "arm": "deprotonated triterpenes",
        "supersedes": ("the three triterpene entries in decoys.json, which were "
                       "docked as neutral COOH and are NOT what the thesis "
                       "screened"),
        "actives": actives,
        "decoys": picked,
        "total_decoys": sum(len(v) for v in picked.values()),
        "matching_windows": WINDOWS,
        "charge_matching": (
            "DROPPED. An exact charge match is unusable here: one candidate at "
            "charge -1 in the whole 1.7M library. Defensible for AutoDock Vina "
            "specifically, which has no electrostatic term, so formal charge is "
            "not a scoring variable. Its physical consequences, hydrogen-bond "
            "donor count and logP, ARE matched. This is still a relaxation and "
            "must be reported with the result."),
        "decoys_disjoint_from_neutral_arm": True,
        "caveat": make_decoys.__doc__.strip().splitlines()[-1] if False else (
            "Decoys are PRESUMED inactive, not measured inactive."),
    }
    with open(TARGET, "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)

    print("")
    print("    %d decoys, %d library molecules scanned" % (manifest["total_decoys"], scanned))
    print("    written to %s" % TARGET)
    print("")
    print("    Next: python run_enrichment.py --decoys %s --out %s"
          % ("validation-inputs/enrichment/decoys_deprotonated.json",
             "validation-run/enrichment-deprotonated"))
    print("")


if __name__ == "__main__":
    main()
