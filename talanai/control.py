"""
Running the redocking control, which is the thing that unlocks docking.

WHAT A REDOCK IS
    Take the ligand that was crystallised inside the receptor, throw away its
    coordinates, dock it back, and measure how far the predicted pose lands
    from where it actually sits. Under about 2 A means the scoring function
    can recover a known answer in this pocket. It is the positive control.

WHY THIS COMMAND EXISTS SEPARATELY
    Because a control only validates the protocol it was run under. Running it
    from the same .tal that describes the screen is what keeps the two in step,
    and the run record it writes is the citable artefact rule R105 demands.

RMSD, HONESTLY
    Heavy atoms only, matched by atom name, no superposition. That is the right
    measure: the question is whether the pose lands in the same place, not
    whether it has the same shape somewhere else.

    It is NOT symmetry-corrected. TalanaiDock's validate_redock.py uses RDKit
    for that, which is a third-party package this tool cannot depend on. For
    a ligand with topologically equivalent atoms (a phenyl that can flip, a
    carboxylate) the uncorrected value can read high when the pose is right.
    alpha-D-glucopyranose has no such degeneracy, so the two agree there. The
    run record says which method was used so nobody has to guess.
"""

from __future__ import annotations

import json
import math
import os
import subprocess
import time

# Read once, not imported at call time, so a record always names the code
# that produced it. Package initialisation guarantees __version__ already
# exists on the talanai package by the time this module's top level runs, so
# this is not a circular import.
from . import __version__ as TALANAI_VERSION


# AutoDock atom types are not element symbols. Most read correctly from their
# first character, including the hydrogen-bond ACCEPTOR types 'OA', 'NA' and
# 'SA', which mean acceptor oxygen, nitrogen and sulfur and NOT sodium. Three
# do not, and taking the first character of those is wrong:
#
#   A     AROMATIC CARBON. Read as element 'A' it matches nothing in a crystal
#         reference, so rmsd() falls through to case 3, finds no reference atom
#         of that element and returns no value at all. Every aromatic ligand is
#         affected: ten of the fifteen prepared ligands in this study carry
#         'A', and the only ligand ever redocked, alpha-D-glucopyranose, is one
#         of the five that does not. That is why this never surfaced in a run.
#   Cl    read as 'C', which let a chlorine pair against a carbon in case 3.
#   Br    read as 'B'.
_ELEMENT_OVERRIDES = {"A": "C", "CL": "CL", "BR": "BR"}


def _element(name, kind):
    """Element symbol from a PDB atom name or an AutoDock type. 'OA' -> O."""
    kind = (kind or "").strip()
    if kind and kind[0].isalpha():
        return _ELEMENT_OVERRIDES.get(kind.upper(), kind[0].upper())
    name = (name or "").strip()
    return name[:1].upper() if name else "?"


def read_heavy_atoms(path, model=1):
    """
    [(name, x, y, z)] for one model of a PDB or PDBQT file, hydrogens dropped.

    PDBQT keeps PDB column positions and puts the AutoDock atom type where the
    element goes, so 'HD' and 'H' are both hydrogen and both get dropped.
    """
    atoms = []
    current = 1
    with open(path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("MODEL"):
                parts = line.split()
                current = int(parts[1]) if len(parts) > 1 else current
                continue
            if line.startswith("ENDMDL"):
                if current == model and atoms:
                    break
                continue
            if not line.startswith(("ATOM", "HETATM")):
                continue
            if current != model:
                continue
            kind = line[76:79].strip() or line[12:16].strip()
            if kind.upper().startswith("H"):
                continue
            try:
                name = line[12:16].strip()
                kind = line[76:79].strip()
                atoms.append((name,
                              float(line[30:38]),
                              float(line[38:46]),
                              float(line[46:54]),
                              _element(name, kind)))
            except ValueError:
                continue
    return atoms


def _rms(pairs):
    total = 0.0
    for (px, py, pz), (rx, ry, rz) in pairs:
        total += (px - rx) ** 2 + (py - ry) ** 2 + (pz - rz) ** 2
    return math.sqrt(total / len(pairs))


def rmsd(pose, reference):
    """
    Heavy-atom RMSD without superposition.

    Returns (value, method, matched, exact) where `exact` says whether the
    number is a real RMSD or a lower bound.

    Atom correspondence is the hard part and it is where a naive
    implementation quietly lies. Three cases, in order of preference:

    1. Names match one-to-one     -> a true RMSD.
    2. Names do not match, but the two files list the same elements in the
       same order -> a true RMSD by file order.
    3. Neither                    -> a LOWER BOUND, computed as the distance
       from each pose atom to the nearest reference atom of the same element.

    Case 3 matters here: TalanaiDock's docking input names every atom 'C' or
    'O' while the crystal reference uses C1..C6 / O1..O6, and their heavy-atom
    orders differ. Comparing by file order would pair carbons against oxygens
    and return a confident, meaningless number.

    A lower bound is only ever conclusive in one direction. If it exceeds the
    threshold the pose definitively failed, because no atom assignment could
    do better. If it falls under the threshold that is NOT a pass, only an
    absence of proof of failure, and the caller must say so.
    """
    if not pose or not reference:
        return None, "no atoms", 0, False
    if len(pose) != len(reference):
        return None, "atom counts differ (%d pose, %d reference)" % (
            len(pose), len(reference)), 0, False

    # 1. by name
    by_name = {a[0]: (a[1], a[2], a[3]) for a in reference}
    if len(by_name) == len(reference):
        pairs = [((a[1], a[2], a[3]), by_name[a[0]])
                 for a in pose if a[0] in by_name]
        if len(pairs) == len(reference):
            return (_rms(pairs), "heavy-atom, by atom name, no superposition",
                    len(pairs), True)

    # 2. by file order, but only when the element sequences agree
    if [a[4] for a in pose] == [a[4] for a in reference]:
        pairs = [((p[1], p[2], p[3]), (r[1], r[2], r[3]))
                 for p, r in zip(pose, reference)]
        return (_rms(pairs),
                "heavy-atom, by file order with matching element sequence, "
                "no superposition", len(pairs), True)

    # 3. lower bound: nearest reference atom of the same element
    pairs = []
    for atom in pose:
        candidates = [r for r in reference if r[4] == atom[4]]
        if not candidates:
            return None, "no reference atom of element %s" % atom[4], 0, False
        nearest = min(candidates,
                      key=lambda r: (atom[1] - r[1]) ** 2 + (atom[2] - r[2]) ** 2
                      + (atom[3] - r[3]) ** 2)
        pairs.append(((atom[1], atom[2], atom[3]),
                      (nearest[1], nearest[2], nearest[3])))
    return (_rms(pairs),
            "LOWER BOUND, nearest same-element atom, correspondence "
            "unresolved (atom names differ between pose and reference)",
            len(pairs), False)


def build_command(vina, receptor, ligand, center, size, exhaustiveness, seed,
                  out_path, modes=9, energy_range=3, cpu=4):
    return [
        vina,
        "--receptor", receptor,
        "--ligand", ligand,
        "--center_x", str(center[0]),
        "--center_y", str(center[1]),
        "--center_z", str(center[2]),
        "--size_x", str(size[0]),
        "--size_y", str(size[1]),
        "--size_z", str(size[2]),
        "--exhaustiveness", str(int(exhaustiveness)),
        "--num_modes", str(int(modes)),
        "--energy_range", str(int(energy_range)),
        "--seed", str(int(seed)),
        "--cpu", str(int(cpu)),
        "--out", out_path,
    ]


def parse_score(stdout):
    for line in stdout.splitlines():
        parts = line.split()
        if len(parts) >= 2 and parts[0] == "1":
            try:
                return float(parts[1])
            except ValueError:
                continue
    return None


def run(vina, receptor, ligand, reference, center, size, exhaustiveness, seed,
        outdir, modes=9, energy_range=3, cpu=4, timeout=1800, label=""):
    """
    Run one redock and return a run record.

    The record is the artefact R105 asks for: it names every input, every
    setting, the RMSD and the method used to compute it, so that the number in
    a .tal file can be traced to the run that produced it.
    """
    os.makedirs(outdir, exist_ok=True)
    stem = "redock_%s_seed%s" % (label or "control", seed)
    out_path = os.path.join(outdir, stem + ".pdbqt")

    command = build_command(vina, receptor, ligand, center, size,
                            exhaustiveness, seed, out_path, modes,
                            energy_range, cpu)
    started = time.time()
    completed = subprocess.run(command, capture_output=True, text=True,
                               timeout=timeout)
    elapsed = time.time() - started

    record = {
        # Written on every record, success or failure, so that R105's "someone
        # else has to be able to find the run this number came from" also
        # answers "and with what code": talanai_version. rmsd_method defaults
        # to None here and is filled in below only when an RMSD was actually
        # computed, so a failed or errored run never leaves the key absent,
        # only null, which is what let redock_screening_seed42.json's stale
        # 7.046 sit unlabelled next to a method string that no longer exists
        # in this file. See that record's "superseded_by" field.
        "talanai_version": TALANAI_VERSION,
        "kind": "redock_control",
        "label": label,
        "vina": os.path.basename(vina),
        "receptor": receptor,
        "ligand": ligand,
        "reference": reference,
        "box_center": center,
        "box_size": size,
        "exhaustiveness": int(exhaustiveness),
        "num_modes": int(modes),
        "energy_range": int(energy_range),
        "seed": int(seed),
        "cpu": int(cpu),
        "command": " ".join(command),
        "returncode": completed.returncode,
        "seconds": round(elapsed, 1),
        "rank1_score_kcal_mol": parse_score(completed.stdout),
        "pose_file": out_path if os.path.isfile(out_path) else None,
        "rmsd_method": None,
    }

    if completed.returncode != 0:
        record["error"] = completed.stderr[-1500:] or completed.stdout[-1500:]
        return record
    if not os.path.isfile(out_path):
        record["error"] = "vina produced no output pose"
        return record

    pose = read_heavy_atoms(out_path, model=1)
    crystal = read_heavy_atoms(reference)
    value, method, matched, exact = rmsd(pose, crystal)
    record["rank1_rmsd_A"] = round(value, 3) if value is not None else None
    record["rmsd_method"] = method
    record["rmsd_atoms_matched"] = matched
    record["rmsd_is_exact"] = exact
    record["rmsd_is_lower_bound"] = not exact and value is not None
    record["rmsd_symmetry_corrected"] = False

    per_pose = []
    index = 1
    while True:
        model = read_heavy_atoms(out_path, model=index)
        if not model:
            break
        value_i, _method, _n, _exact = rmsd(model, crystal)
        per_pose.append(round(value_i, 3) if value_i is not None else None)
        index += 1
        if index > 20:
            break
    record["per_pose_rmsd_A"] = per_pose
    best = [v for v in per_pose if v is not None]
    record["best_pose_rmsd_A"] = min(best) if best else None
    return record


def store(record, outdir, name):
    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, name)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(record, handle, indent=2, sort_keys=True)
    return path
