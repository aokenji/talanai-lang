#!/usr/bin/env python3
"""
Build a decoy source library by pulling only the ZINC20 tranches that bracket
the actives.

WHY NOT JUST DOWNLOAD A BIG LIBRARY
    ZINC20's 2D area is organised as 121 tranches, molecular weight crossed
    with logP, each roughly 150 MB. Pulling the lot to find 200 decoys would
    be absurd. The actives here span molecular weight 286 to 610 and logP from
    about -1 (rutin, a glycoside) to about +7 (the triterpenes), so only a
    handful of tranches can possibly contain a property match, and those are
    the only ones fetched. One shard per tranche, a few MB each.

WHY ZINC AND NOT CHEMBL
    ChEMBL is a bioactivity database, so a "decoy" drawn from it is a molecule
    someone already had a reason to test. Some fraction would be real
    alpha-glucosidase actives, and every one of those counts as a false
    positive and pushes the measured enrichment down. ZINC is purchasable
    compound space, mostly never assayed against anything, which is what a
    presumed-inactive decoy is supposed to be. This is the same reason DUD-E
    draws from ZINC.

    Presumed inactive is still not measured inactive. That caveat travels with
    the result and is written into the decoy manifest.

    python fetch_decoy_library.py
    python fetch_decoy_library.py --shards 2      (more molecules per tranche)

Writes validation-inputs/enrichment/zinc_library.smi, ready for:

    python make_decoys.py --library validation-inputs/enrichment/zinc_library.smi
"""

import argparse
import json
import os
import re
import sys
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from rdkit import Chem, RDLogger                     # noqa: E402
from rdkit.Chem import Crippen, Descriptors          # noqa: E402
RDLogger.DisableLog("rdApp.*")

import make_decoys                                    # noqa: E402

BASE = "https://files.docking.org/2D/"
OUT = os.path.join(HERE, "validation-inputs", "enrichment")
LIBRARY = os.path.join(OUT, "zinc_library.smi")
UA = {"User-Agent": "talanai-enrichment/0.1"}

# ZINC20 tranche bin edges. Upper bound of each lettered bin; the last letter
# is everything above. Taken from the tranche browser's own axes.
MW_EDGES = [200, 250, 300, 325, 350, 375, 400, 425, 450, 500]
LOGP_EDGES = [-1, 0, 1, 2, 2.5, 3, 3.5, 4, 4.5, 5]
LETTERS = "ABCDEFGHIJK"


def bin_letter(value, edges):
    for i, edge in enumerate(edges):
        if value <= edge:
            return LETTERS[i]
    return LETTERS[len(edges)]


def tranche_for(mw, logp):
    return bin_letter(mw, MW_EDGES) + bin_letter(logp, LOGP_EDGES)


def neighbours(code):
    """The tranche plus its immediate neighbours, so a property near a bin
    edge does not lose the half of its matches that fell the other side."""
    mw_i, logp_i = LETTERS.index(code[0]), LETTERS.index(code[1])
    out = []
    for dm in (-1, 0, 1):
        for dl in (-1, 0, 1):
            m, l = mw_i + dm, logp_i + dl
            if 0 <= m < len(LETTERS) and 0 <= l < len(LETTERS):
                out.append(LETTERS[m] + LETTERS[l])
    return out


def get(url, limit=None):
    request = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(request, timeout=120) as response:
        return response.read(limit) if limit else response.read()


def shards_in(code):
    try:
        html = get(BASE + code + "/", 200000).decode("utf-8", "replace")
    except Exception as error:
        print("      %s unavailable (%s)" % (code, type(error).__name__))
        return []
    return sorted(set(re.findall(r'href="([^"?][^"/]*\.smi)"', html)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--shards", type=int, default=1,
                        help="shard files per tranche (default 1, a few MB each)")
    parser.add_argument("--only", default="",
                        help="comma separated tranche codes, e.g. KA,KB. Use with "
                             "--shards to top up a sparse corner of property "
                             "space without re-pulling the dense tranches.")
    parser.add_argument("--append", action="store_true",
                        help="add to the existing library instead of replacing it")
    args = parser.parse_args()
    os.makedirs(OUT, exist_ok=True)

    print("")
    print("  STEP 1  where do the actives sit in property space")
    print("  " + "-" * 74)
    actives, unmeasured = make_decoys.load_actives()

    wanted, rows = [], []
    for active in actives:
        mol = Chem.MolFromSmiles(active["smiles"])
        mw, logp = Descriptors.MolWt(mol), Crippen.MolLogP(mol)
        code = tranche_for(mw, logp)
        rows.append((active["name"], mw, logp, code))
        for neighbour in neighbours(code):
            if neighbour not in wanted:
                wanted.append(neighbour)
        print("    %-16s MW %6.1f  logP %5.2f  ->  tranche %s"
              % (active["name"], mw, logp, code))
    for row in unmeasured:
        print("    %-16s not an active (no published IC50), no tranche needed"
              % row["name"])

    if args.only:
        asked = [c.strip().upper() for c in args.only.split(",") if c.strip()]
        wanted = [c for c in asked]
        print("")
        print("    --only given: fetching just %s" % ", ".join(wanted))

    print("")
    print("    %d tranches to fetch, including neighbours of each" % len(wanted))
    print("    " + " ".join(sorted(wanted)))

    print("")
    print("  STEP 2  fetching")
    print("  " + "-" * 74)
    total_bytes, total_lines, fetched = 0, 0, []
    seen = set()
    if args.append and os.path.isfile(LIBRARY):
        with open(LIBRARY, encoding="utf-8") as handle:
            for line in handle:
                parts = line.split()
                if parts:
                    seen.add(parts[0])
        print("    appending to %d molecules already in the library" % len(seen))
    with open(LIBRARY, "a" if args.append else "w", encoding="utf-8") as out:
        for code in sorted(wanted):
            names = shards_in(code)
            if not names:
                continue
            for name in names[:args.shards]:
                url = BASE + code + "/" + name
                try:
                    blob = get(url)
                except Exception as error:
                    print("      %-10s FAILED %s" % (name, type(error).__name__))
                    continue
                text = blob.decode("utf-8", "replace")
                kept = 0
                for line in text.splitlines():
                    line = line.strip()
                    if not line or line.lower().startswith("smiles"):
                        continue
                    parts = line.split()
                    if len(parts) < 2 or parts[0] in seen:
                        continue
                    seen.add(parts[0])
                    out.write("%s %s\n" % (parts[0], parts[1]))
                    kept += 1
                total_bytes += len(blob)
                total_lines += kept
                fetched.append({"tranche": code, "shard": name,
                                "bytes": len(blob), "molecules": kept})
                print("      %-4s %-14s %6.1f MB  %7d molecules"
                      % (code, name, len(blob) / 1e6, kept))

    manifest = {
        "source": "ZINC20 2D tranches, files.docking.org",
        "why_zinc": ("purchasable compound space, largely unassayed, so a "
                     "decoy drawn from it is plausibly inactive. ChEMBL would "
                     "seed real actives into the decoy set."),
        "actives_property_space": [
            {"name": n, "mw": round(m, 2), "logp": round(p, 2), "tranche": c}
            for n, m, p, c in rows],
        "tranches_requested": sorted(wanted),
        "shards_fetched": fetched,
        "total_molecules": total_lines,
        "total_megabytes": round(total_bytes / 1e6, 1),
        "library_file": LIBRARY,
        "caveat": ("Decoys are PRESUMED inactive. Any genuine active in here "
                   "counts as a false positive and lowers the measured "
                   "enrichment, so this cannot inflate the result."),
    }
    manifest_name = ("library_manifest_topup.json" if args.only
                     else "library_manifest.json")
    with open(os.path.join(OUT, manifest_name), "w",
              encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)

    print("")
    print("  %d molecules, %.1f MB, written to" % (total_lines, total_bytes / 1e6))
    print("  %s" % LIBRARY)
    print("")
    print("  Next: python make_decoys.py --library %s"
          % os.path.relpath(LIBRARY, HERE).replace("\\", "/"))
    print("")


if __name__ == "__main__":
    main()
