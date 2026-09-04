#!/usr/bin/env python3
"""
Property-matched decoys for an enrichment benchmark.

WHAT AN ENRICHMENT TEST ASKS
    Not "does the score predict potency" - that was already tested
    (run_ic50_correlation.py) and it FAILED: Spearman +0.09 to +0.16 against
    published yeast IC50, indistinguishable from chance.

    This asks something different and easier: can the protocol put KNOWN
    ACTIVES above molecules that merely look like them? That is binary
    discrimination rather than rank ordering, it is what docking is actually
    decent at, and it is the ability the thesis's own claim rests on: a
    shortlist worth following up, not a potency prediction.

WHY THE DECOYS MUST BE PROPERTY MATCHED
    A decoy set drawn at random is trivially easy. Vina's score rises with
    molecular size, so a screen of large actives against small random
    molecules "enriches" beautifully and means nothing: it has rediscovered
    that the actives are big.

    So each decoy is matched to an active on the properties that drive the
    score (heavy atoms, molecular weight, logP, rotatable bonds, H-bond donors
    and acceptors, net charge) and required to be TOPOLOGICALLY DISSIMILAR to
    every active, so it is unlikely to share the real pharmacophore. That is
    the DUD-E construction, and it is what makes the result interpretable.

    Physically decoys are PRESUMED inactive, not measured inactive. With a
    20:1 ratio a few genuine actives may be hiding among them. That biases the
    result DOWNWARD, so it cannot manufacture a pass.

THE LIBRARY
    This project has no compound library on disk, so one has to be supplied.
    Any SMILES file works: one molecule per line, SMILES first, optional id
    second, whitespace separated. ZINC's drug-like subsets and the ChEMBL
    SMILES dump are both suitable and free.

        python make_decoys.py --library zinc_druglike.smi
        python make_decoys.py --library zinc_druglike.smi --per-active 20

    Actives and their SMILES come from PubChem by CID, the CIDs recorded in
    src/data/compounds.js, fetched once and cached in the output directory so
    a rerun is reproducible without the network.

WHAT IS AND IS NOT AN ACTIVE HERE
    An active is a screened compound with a TRACEABLE published IC50 against
    the yeast enzyme, read out of validation-run/ic50/ic50.json. Spinosin is
    screened but has no alpha-glucosidase measurement of any kind, so it is
    docked with the rest and excluded from the active set. Calling it an
    active would be asserting an activity nobody has measured.

    The set is therefore derived from the data, not typed in here. If a
    spinosin IC50 is ever published, drop it into ic50.json and the active
    set grows on its own.
"""

import argparse
import json
import os
import sys
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from rdkit import Chem, RDLogger                            # noqa: E402
from rdkit.Chem import AllChem, Crippen, Descriptors, rdMolDescriptors  # noqa: E402
from rdkit.Chem import DataStructs                          # noqa: E402
RDLogger.DisableLog("rdApp.*")

IC50 = os.path.join(HERE, "validation-run", "ic50", "ic50.json")
OUT = os.path.join(HERE, "validation-inputs", "enrichment")

PUBCHEM = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%d/"
           "property/IsomericSMILES,CanonicalSMILES/JSON")

# The ten screened compounds, by the PubChem CID recorded in compounds.js.
# Every SMILES is fetched from PubChem rather than typed, because a
# hand-copied SMILES is how three ligands lost their stereochemistry once.
SCREENED = [
    ("Quercetin", 5280343), ("Rutin", 5280805), ("Kaempferol", 5280863),
    ("Betulinic Acid", 64971), ("Oleanolic Acid", 10494), ("Luteolin", 5280445),
    ("Vitexin", 5280441), ("Isovitexin", 162350), ("Ursolic Acid", 64945),
    ("Spinosin", 155692),
]

# Matching windows. Deliberately tight on the properties Vina's scoring and
# its rotatable-bond normalisation actually respond to.
WINDOWS = {
    "heavy":  6,      # heavy atom count, the dominant driver of raw score
    "mw":    40.0,    # g/mol
    "logp":   1.0,    # Crippen
    "rotb":   2,      # rotatable bonds, which Vina divides the score by
    "hbd":    1,
    "hba":    2,
    "charge": 0,      # exact match; formal charge changes the interaction set
}

# Above this Tanimoto to ANY active, a candidate is too similar to be a decoy:
# it might genuinely share the pharmacophore, so counting it as a false
# positive would understate the protocol.
MAX_SIMILARITY = 0.35


def props(mol):
    return {
        "heavy":  mol.GetNumHeavyAtoms(),
        "mw":     Descriptors.MolWt(mol),
        "logp":   Crippen.MolLogP(mol),
        "rotb":   rdMolDescriptors.CalcNumRotatableBonds(mol),
        "hbd":    rdMolDescriptors.CalcNumHBD(mol),
        "hba":    rdMolDescriptors.CalcNumHBA(mol),
        "charge": Chem.GetFormalCharge(mol),
    }


def within(candidate, target):
    for key, tol in WINDOWS.items():
        if abs(candidate[key] - target[key]) > tol:
            return False
    return True


def fingerprint(mol):
    return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)


def fetch_smiles(cid):
    with urllib.request.urlopen(PUBCHEM % cid, timeout=30) as response:
        payload = json.load(response)
    row = payload["PropertyTable"]["Properties"][0]
    return row.get("IsomericSMILES") or row["CanonicalSMILES"]


def load_actives():
    """
    The screened compounds, split into actives and not, using ic50.json as the
    only authority on which is which.
    """
    with open(IC50, encoding="utf-8") as handle:
        ic50 = json.load(handle)
    measured = ic50.get("compounds", {})

    cache = os.path.join(OUT, "screened_smiles.json")
    known = {}
    if os.path.isfile(cache):
        with open(cache, encoding="utf-8") as handle:
            known = json.load(handle)

    rows, missing = [], []
    for name, cid in SCREENED:
        if name not in known:
            print("      fetching SMILES for %s (CID %d)" % (name, cid))
            known[name] = fetch_smiles(cid)
        mol = Chem.MolFromSmiles(known[name])
        if mol is None:
            raise SystemExit("could not parse the SMILES for %s" % name)

        entry = measured.get(name) or measured.get(name.replace(" ", "_"))
        yeast = (entry or {}).get("yeast_values") or []
        row = {"name": name, "cid": cid, "smiles": known[name],
               "props": props(mol), "n_yeast_ic50_sources": len(yeast)}
        if yeast:
            rows.append(row)
        else:
            missing.append(row)

    with open(cache, "w", encoding="utf-8") as handle:
        json.dump(known, handle, indent=2, sort_keys=True)
    return rows, missing


def read_library(path):
    """SMILES first token per line. Blank lines and '#' comments skipped."""
    seen = set()
    with open(path, encoding="utf-8", errors="replace") as handle:
        for lineno, line in enumerate(handle, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            smiles = parts[0]
            if smiles in seen:
                continue
            seen.add(smiles)
            yield smiles, (parts[1] if len(parts) > 1 else "lib%d" % lineno)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--library", required=True,
                        help="SMILES file to draw decoys from")
    parser.add_argument("--per-active", type=int, default=20,
                        help="decoys per active (default 20)")
    parser.add_argument("--max-scan", type=int, default=500000,
                        help="stop after this many library molecules")
    args = parser.parse_args()

    os.makedirs(OUT, exist_ok=True)

    print("")
    print("  STEP 1  which screened compounds are known actives")
    print("  " + "-" * 74)
    actives, unmeasured = load_actives()
    for row in actives:
        print("    ACTIVE   %-16s %d published yeast IC50 source(s)"
              % (row["name"], row["n_yeast_ic50_sources"]))
    for row in unmeasured:
        print("    excluded %-16s no yeast alpha-glucosidase measurement exists"
              % row["name"])
    print("")
    print("    %d actives, %d screened compounds excluded from the active set."
          % (len(actives), len(unmeasured)))
    print("    Excluded compounds are still docked; they are simply not")
    print("    evidence of anything in an enrichment calculation.")

    if not actives:
        raise SystemExit("no actives with published IC50; nothing to benchmark")

    active_fps = [fingerprint(Chem.MolFromSmiles(a["smiles"])) for a in actives]
    wanted = args.per_active
    picked = {a["name"]: [] for a in actives}
    taken = set()

    print("")
    print("  STEP 2  property matching against %s" % os.path.basename(args.library))
    print("  " + "-" * 74)
    print("    windows: " + ", ".join("%s +-%s" % (k, v) for k, v in WINDOWS.items()))
    print("    rejecting anything with Tanimoto > %.2f to any active"
          % MAX_SIMILARITY)
    print("")

    scanned = 0
    for smiles, ident in read_library(args.library):
        scanned += 1
        if scanned > args.max_scan:
            print("    reached --max-scan of %d, stopping the scan" % args.max_scan)
            break
        if all(len(v) >= wanted for v in picked.values()):
            break

        short = [a for a in actives if len(picked[a["name"]]) < wanted]
        if not short:
            break

        mol = Chem.MolFromSmiles(smiles)
        if mol is None or mol.GetNumHeavyAtoms() < 8:
            continue
        candidate = props(mol)

        target = next((a for a in short if within(candidate, a["props"])), None)
        if target is None or smiles in taken:
            continue

        fp = fingerprint(mol)
        if max(DataStructs.BulkTanimotoSimilarity(fp, active_fps)) > MAX_SIMILARITY:
            continue

        taken.add(smiles)
        picked[target["name"]].append({"smiles": smiles, "library_id": ident,
                                       "props": candidate})
        if scanned % 25000 == 0:
            have = sum(len(v) for v in picked.values())
            print("    scanned %7d   matched %4d / %d"
                  % (scanned, have, wanted * len(actives)))

    print("")
    print("  STEP 3  result")
    print("  " + "-" * 74)
    short_of_target = []
    for a in actives:
        got = len(picked[a["name"]])
        flag = "" if got >= wanted else "   <-- SHORT"
        if got < wanted:
            short_of_target.append((a["name"], got))
        print("    %-16s %3d decoys%s" % (a["name"], got, flag))

    total = sum(len(v) for v in picked.values())
    manifest = {
        "generated_from": os.path.abspath(args.library),
        "library_molecules_scanned": scanned,
        "decoys_per_active_target": wanted,
        "matching_windows": WINDOWS,
        "max_tanimoto_to_any_active": MAX_SIMILARITY,
        "actives": actives,
        "excluded_from_active_set": [
            {"name": r["name"], "reason":
             "no published yeast alpha-glucosidase IC50 in ic50.json"}
            for r in unmeasured],
        "decoys": picked,
        "total_decoys": total,
        "caveat": (
            "Decoys are PRESUMED inactive, not measured inactive. Any genuine "
            "active hiding in the decoy set counts as a false positive and "
            "pushes the measured enrichment DOWN, so this construction cannot "
            "inflate the result."),
    }
    if short_of_target:
        manifest["incomplete"] = [
            {"active": n, "decoys_found": g, "wanted": wanted}
            for n, g in short_of_target]
        manifest["incomplete_note"] = (
            "Some actives did not reach the requested decoy count. Report the "
            "per-active counts rather than a single ratio, and do not describe "
            "this as a 20:1 benchmark if it is not one.")

    path = os.path.join(OUT, "decoys.json")
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)

    print("")
    print("    %d decoys for %d actives, scanned %d library molecules"
          % (total, len(actives), scanned))
    print("    written to %s" % path)
    if short_of_target:
        print("")
        print("    *** %d active(s) did not reach %d decoys. A bigger or more"
              % (len(short_of_target), wanted))
        print("        chemically diverse library is the fix. Do NOT loosen the")
        print("        matching windows to make up the numbers: that is how a")
        print("        benchmark quietly becomes easy.")
    print("")
    print("    Next: python run_enrichment.py")
    print("")


if __name__ == "__main__":
    main()
