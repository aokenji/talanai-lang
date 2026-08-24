"""
Just enough PDB reading to check that a docking box contains the pocket it
claims to contain.

This is not a structure library. It reads ATOM and HETATM records by column,
per the PDB format specification, and returns coordinates. That is all the
box check needs, and it keeps the stdlib-only constraint.
"""

from __future__ import annotations

import os

THREE_LETTER = {
    "ala": "ALA", "arg": "ARG", "asn": "ASN", "asp": "ASP", "cys": "CYS",
    "gln": "GLN", "glu": "GLU", "gly": "GLY", "his": "HIS", "ile": "ILE",
    "leu": "LEU", "lys": "LYS", "met": "MET", "phe": "PHE", "pro": "PRO",
    "ser": "SER", "thr": "THR", "trp": "TRP", "tyr": "TYR", "val": "VAL",
}


def parse_residue_token(token):
    """'Asp215' -> ('ASP', 215). 'A/Asp215' -> ('ASP', 215) with chain A."""
    chain = None
    if "/" in token:
        chain, token = token.split("/", 1)
    letters = "".join(c for c in token if c.isalpha())
    digits = "".join(c for c in token if c.isdigit())
    if not digits:
        return None
    name = THREE_LETTER.get(letters.lower()[:3], letters.upper()[:3])
    return name, int(digits), chain


def read_atoms(path):
    """Yield (record, resname, chain, resseq, x, y, z) for every atom."""
    with open(path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            record = line[:6].strip()
            if record not in ("ATOM", "HETATM"):
                continue
            try:
                resname = line[17:20].strip().upper()
                chain = line[21].strip()
                resseq = int(line[22:26])
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except (ValueError, IndexError):
                continue
            yield record, resname, chain, resseq, x, y, z


def residue_centroid(path, resname, resseq, chain=None):
    """Mean coordinate of one residue, or None when it is not in the file."""
    total = [0.0, 0.0, 0.0]
    count = 0
    for _rec, name, ch, seq, x, y, z in read_atoms(path):
        if seq != resseq or name != resname:
            continue
        if chain and ch != chain:
            continue
        total[0] += x
        total[1] += y
        total[2] += z
        count += 1
    if not count:
        return None
    return [c / count for c in total]


def inside_box(point, center, size):
    """True when point lies within the axis-aligned box."""
    for value, mid, extent in zip(point, center, size):
        if abs(value - mid) > extent / 2.0:
            return False
    return True


def distance_outside(point, center, size):
    """How far outside the box, per axis. Zero when inside."""
    out = []
    for value, mid, extent in zip(point, center, size):
        overshoot = abs(value - mid) - extent / 2.0
        out.append(max(0.0, overshoot))
    return out


def resolve(path_hint, base_dir):
    """Find a receptor file named in a .tal, relative to the file itself."""
    if not path_hint:
        return None
    candidates = [path_hint, os.path.join(base_dir, path_hint)]
    for candidate in candidates:
        if os.path.isfile(candidate):
            return candidate
    return None
