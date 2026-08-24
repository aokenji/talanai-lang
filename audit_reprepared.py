#!/usr/bin/env python3
"""
Did the 50-conformer re-preparation actually fix the rings?

prep_multiconf's own crude check (mean |torsion| < 45) flagged rings in all
three re-prepared ligands, but that check cannot tell a defect from a ring that
contains an sp2 atom and therefore cannot be a chair. This reuses the rigorous
sp2-aware Cremer-Pople classifier on the re-prepared files, so the claim made
to a panel is one that was actually measured.
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import audit_rings_rigorous as audit   # noqa: E402

REPREP = os.path.join(HERE, "validation-inputs", "reprepared")

audit.LIGANDS = [
    ("Ursolic OLD", os.path.join(audit.DATA, "ligands", "ursolic.pdbqt")),
    ("Ursolic NEW", os.path.join(REPREP, "ursolic.pdbqt")),
    ("Isovitexin OLD", os.path.join(audit.DATA, "ligands", "isovitexin.pdbqt")),
    ("Isovitexin NEW", os.path.join(REPREP, "isovitexin.pdbqt")),
    ("Acarbose OLD", os.path.join(audit.DATA, "reference", "acarbose.pdbqt")),
    ("Acarbose NEW", os.path.join(REPREP, "acarbose.pdbqt")),
]

if __name__ == "__main__":
    sys.exit(audit.main())
