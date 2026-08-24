#!/usr/bin/env python3
"""
Is the glucose orientation failure a SAMPLING limit or a SCORING limit?

THE QUESTION
    Control 2 showed glucose built through the study's ligand pipeline docks
    5.358 A from its crystal pose, while the same molecule fed as crystal
    coordinates lands at 0.518 A. Both engage the full catalytic triad, and
    they score only 0.13 kcal/mol apart.

    Two explanations, needing two very different fixes:

    SAMPLING  the search never visits the crystal orientation. More effort
              would find it. Fix: raise exhaustiveness and re-run the screen.
    SCORING   the search does visit it but ranks another orientation above it.
              No amount of search fixes that. Fix: pose-level claims must be
              restricted to the ensemble or dropped.

THE DISCRIMINATOR
    Rank-1 RMSD cannot separate these. Best-of-nine can:

      crystal orientation absent from all 9 poses, even at high effort
          -> SAMPLING limit, unfixed at the effort tested
      crystal orientation present among the 9 but never ranked 1st
          -> SCORING limit
      crystal orientation becomes rank 1 as effort rises
          -> SAMPLING limit, and fixed

    So every pose of every ensemble is measured, not just the winner.

    python run_orientation_test.py
"""

import json
import os
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from talanai import control as tal_control   # noqa: E402

ASSETS = (r"D:\THESIS_VSC\BALAKATDBV2\talanaidock\dist\TalanaiDock\app"
          r"\docking_assets")
DATA = os.path.join(ASSETS, "docking_data")
VINA = os.path.join(ASSETS, "vina.exe")
RECEPTOR = os.path.join(DATA, "prepared", "receptor.pdbqt")
REFERENCE = os.path.join(DATA, "validation", "glc_3a4a_601.pdb")

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "validation-run", "orientation")

CENTER = ["21.52", "-7.7", "23.55"]
BOX = ["30", "30", "30"]
THRESHOLD = 2.0

LIGANDS = {
    "pipeline": os.path.join(HERE, "validation-inputs", "glc_reembedded.pdbqt"),
    "crystal":  os.path.join(DATA, "validation", "glc_ligand.pdbqt"),
}

# The pipeline ligand failed, so it gets the escalation. The crystal ligand
# gets one high-effort run as a comparator: if it holds at 0.5 A while the
# pipeline one does not, the difference is the starting conformation, not the
# search budget.
PLAN = (
    [("pipeline", 128, s) for s in (42, 43, 44)]
    + [("pipeline", 512, s) for s in (42, 43, 44)]
    + [("crystal", 512, 42)]
)


def per_pose_rmsd(pose_file, reference):
    """
    Symmetry-corrected, NO superposition, for every pose in the file.

    Returns (values, method). Falls back to Talanai's standard-library lower
    bound if RDKit is unavailable, and says so, because a bound and a true
    value must never be silently interchanged.
    """
    note = ""
    try:
        from rdkit import Chem, RDLogger
        from rdkit.Chem import rdMolAlign
        from meeko import PDBQTMolecule, RDKitMolCreate
        RDLogger.DisableLog("rdApp.*")

        ref = Chem.MolFromPDBFile(reference, removeHs=True, sanitize=False)
        mols = RDKitMolCreate.from_pdbqt_mol(
            PDBQTMolecule.from_file(pose_file, skip_typing=True))
        mol = Chem.RemoveHs(mols[0])
        out = []
        for conf in range(mol.GetNumConformers()):
            out.append(round(rdMolAlign.CalcRMS(mol, ref, prbId=conf, refId=0,
                                                maxMatches=1000000), 3))
        if out:
            return out, "rdkit CalcRMS, symmetry-corrected, no superposition"
        note = "rdkit returned no conformers; "
    except Exception as exc:
        note = "rdkit path failed (%s: %s); " % (type(exc).__name__, exc)

    crystal = tal_control.read_heavy_atoms(reference)
    out, index = [], 1
    while True:
        model = tal_control.read_heavy_atoms(pose_file, model=index)
        if not model:
            break
        value, _m, _n, _e = tal_control.rmsd(model, crystal)
        out.append(round(value, 3) if value is not None else None)
        index += 1
        if index > 30:
            break
    return out, note + "fell back to Talanai LOWER BOUND"


def run(tag, ligand, exhaustiveness, seed):
    os.makedirs(OUT, exist_ok=True)
    stem = "%s_e%d_s%d" % (tag, exhaustiveness, seed)
    pose_out = os.path.join(OUT, stem + ".pdbqt")
    command = [VINA, "--receptor", RECEPTOR, "--ligand", ligand,
               "--center_x", CENTER[0], "--center_y", CENTER[1],
               "--center_z", CENTER[2],
               "--size_x", BOX[0], "--size_y", BOX[1], "--size_z", BOX[2],
               "--exhaustiveness", str(exhaustiveness), "--num_modes", "9",
               "--energy_range", "3", "--seed", str(seed), "--cpu", "4",
               "--out", pose_out]
    started = time.time()
    done = subprocess.run(command, capture_output=True, text=True, timeout=21600)
    elapsed = time.time() - started

    record = {"ligand_kind": tag, "ligand": ligand,
              "exhaustiveness": exhaustiveness, "seed": seed,
              "seconds": round(elapsed, 1), "returncode": done.returncode,
              "command": " ".join(command), "talanai_version": "0.1.0",
              "pose_file": pose_out if os.path.isfile(pose_out) else None}
    if done.returncode != 0 or not os.path.isfile(pose_out):
        record["error"] = (done.stderr or done.stdout)[-1200:]
        return record

    record["rank1_score_kcal_mol"] = tal_control.parse_score(done.stdout)
    values, method = per_pose_rmsd(pose_out, REFERENCE)
    clean = [v for v in values if v is not None]
    record["per_pose_rmsd_A"] = values
    record["rmsd_method"] = method
    record["rank1_rmsd_A"] = values[0] if values else None
    record["best_of_ensemble_rmsd_A"] = min(clean) if clean else None
    record["best_pose_index"] = (values.index(min(clean)) + 1) if clean else None
    record["crystal_orientation_in_ensemble"] = bool(clean and min(clean) < THRESHOLD)
    return record


def main():
    os.makedirs(OUT, exist_ok=True)
    print("")
    print("  Orientation test: sampling limit or scoring limit?")
    print("  receptor %s, box %s A, threshold %.1f A"
          % (os.path.basename(RECEPTOR), BOX[0], THRESHOLD))
    print("  " + "-" * 74)
    print("  %-9s %5s %5s %9s %8s %8s %5s"
          % ("ligand", "exh", "seed", "score", "rank1", "best/9", "at#"))

    results = []
    for tag, exhaustiveness, seed in PLAN:
        ligand = LIGANDS[tag]
        if not os.path.isfile(ligand):
            print("  %-9s MISSING %s" % (tag, ligand))
            continue
        record = run(tag, ligand, exhaustiveness, seed)
        results.append(record)
        with open(os.path.join(OUT, "%s_e%d_s%d.json" % (tag, exhaustiveness, seed)),
                  "w", encoding="utf-8") as handle:
            json.dump(record, handle, indent=2, sort_keys=True)
        if record.get("error"):
            print("  %-9s %5d %5d   FAILED: %s"
                  % (tag, exhaustiveness, seed,
                     record["error"].splitlines()[-1][:38]), flush=True)
            continue
        print("  %-9s %5d %5d %9.3f %8s %8s %5s"
              % (tag, exhaustiveness, seed, record["rank1_score_kcal_mol"],
                 record["rank1_rmsd_A"], record["best_of_ensemble_rmsd_A"],
                 record["best_pose_index"]), flush=True)

    print("  " + "-" * 74)
    pipeline = [r for r in results
                if r.get("ligand_kind") == "pipeline" and not r.get("error")]
    if not pipeline:
        print("  No pipeline runs completed. No verdict.")
        return 1

    rank1_ok = [r for r in pipeline if (r["rank1_rmsd_A"] or 99) < THRESHOLD]
    anywhere = [r for r in pipeline if r["crystal_orientation_in_ensemble"]]

    print("")
    print("  VERDICT")
    if rank1_ok:
        print("  SAMPLING LIMIT, fixable by search effort.")
        print("  The crystal orientation becomes rank 1 in %d of %d pipeline runs."
              % (len(rank1_ok), len(pipeline)))
        print("  Re-run the screen at that exhaustiveness.")
    elif anywhere:
        print("  SCORING LIMIT.")
        print("  The crystal orientation IS found within %.1f A in %d of %d runs,"
              % (THRESHOLD, len(anywhere), len(pipeline)))
        print("  but is never ranked first. More search will not fix this.")
        print("  Pose-level claims must move to the ensemble, or be dropped.")
    else:
        print("  SAMPLING LIMIT, NOT fixed at the effort tested.")
        print("  The crystal orientation appears in no pose of any run, up to")
        print("  exhaustiveness %d. Either far more effort is needed, or the"
              % max(r["exhaustiveness"] for r in pipeline))
        print("  prepared ligand differs from the crystal one in a way that")
        print("  makes the crystal pose unreachable. Check the embedding.")
    print("")
    return 0


if __name__ == "__main__":
    sys.exit(main())
