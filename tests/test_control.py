#!/usr/bin/env python3
"""
Unit tests for the redocking control: atom typing, structure reading, and the
RMSD itself. Standard library unittest, no pytest.

    python tests/test_control.py

None of this shells out. Running a redock needs Vina; measuring one does not,
and measuring one is where the mistakes have actually happened.

WHY THIS FILE EXISTS
    rmsd() decides whether the protocol is allowed to dock at all, and until
    now nothing tested it. Two real errors came out of that:

    1. A published number was withdrawn. redock_screening_seed42.json records
       rank1_rmsd_A 7.046 superseded by 6.164, because the method then in use
       paired heavy atoms strictly by file order without first checking that
       the two files listed the same elements in the same sequence, so it
       matched pose O5 against reference C6 and reported the result as exact.
    2. Aromatic ligands returned no value at all, because AutoDock's aromatic
       carbon type 'A' was read as an element named 'A'.

    Both are pinned below. Tests named KNOWN GAP assert what the code does
    today rather than what it should do; each says what would be better.
"""

from __future__ import annotations

import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

from talanai import control   # noqa: E402

CATECHIN = os.path.join(ROOT, "validation-inputs", "flavan3ol-pilot",
                        "catechin.pdbqt")
ISOMALTOSE_REF = os.path.join(ROOT, "validation-inputs",
                              "isomaltose_3axh_ref.pdb")
NINE_POSES = os.path.join(ROOT, "examples", "alpha-glucosidase.out", "control",
                          "redock_screening_seed42.pdbqt")


def _atom_line(serial, name, x, y, z, kind, record="ATOM"):
    """
    One fixed-column atom line, filled in only as precisely as control.py
    reads it: record at columns 1-6, name at 13-16, x/y/z at 31-54, and the
    element or AutoDock type at 77-79.
    """
    line = list(" " * 80)
    line[0:len(record)] = list(record)
    line[6:11] = list("%5d" % serial)
    line[12:16] = list("%-4s" % name)
    line[17:20] = list("LIG")
    line[21] = "A"
    line[22:26] = list("%4d" % 1)
    line[30:38] = list("%8.3f" % x)
    line[38:46] = list("%8.3f" % y)
    line[46:54] = list("%8.3f" % z)
    line[76:79] = list("%3s" % kind)
    return "".join(line)


def _write(lines, suffix=".pdb"):
    handle = tempfile.NamedTemporaryFile(
        mode="w", suffix=suffix, delete=False, dir=HERE, encoding="utf-8")
    handle.write("\n".join(lines) + "\n")
    handle.close()
    return handle.name


def _atoms(spec):
    """[(name, x, y, z, element)] the way read_heavy_atoms returns it."""
    return [(n, x, y, z, e) for (n, x, y, z, e) in spec]


class TestElement(unittest.TestCase):

    def test_acceptor_types_are_their_own_element_not_sodium(self):
        # 'NA', 'OA' and 'SA' are AutoDock hydrogen-bond ACCEPTOR nitrogen,
        # oxygen and sulfur. Reading 'NA' as sodium would be a plausible
        # looking "fix" and would be wrong.
        self.assertEqual(control._element("O1", "OA"), "O")
        self.assertEqual(control._element("N1", "NA"), "N")
        self.assertEqual(control._element("S1", "SA"), "S")

    def test_aromatic_carbon_reads_as_carbon(self):
        # REGRESSION. 'A' is AutoDock's aromatic carbon. While it read as an
        # element named 'A', no crystal reference could ever contain a match,
        # so rmsd() returned no value for any aromatic ligand.
        self.assertEqual(control._element("C4", "A"), "C")

    def test_halogens_keep_both_letters(self):
        # REGRESSION. 'Cl' read as 'C' let a chlorine pair against a carbon.
        self.assertEqual(control._element("CL1", "Cl"), "CL")
        self.assertEqual(control._element("BR1", "Br"), "BR")

    def test_single_letter_types_pass_through(self):
        for kind in ("C", "N", "O", "S", "P", "F", "I"):
            self.assertEqual(control._element("X1", kind), kind)

    def test_falls_back_to_the_atom_name_without_an_element_column(self):
        self.assertEqual(control._element("C1", ""), "C")
        self.assertEqual(control._element("O5", "   "), "O")

    def test_nothing_at_all_reads_as_unknown(self):
        self.assertEqual(control._element("", ""), "?")


class TestReadHeavyAtoms(unittest.TestCase):

    def test_hydrogens_are_dropped_by_both_spellings(self):
        path = _write([
            _atom_line(1, "C1", 0.0, 0.0, 0.0, "C"),
            _atom_line(2, "H1", 1.0, 0.0, 0.0, "H"),     # PDB hydrogen
            _atom_line(3, "HO1", 2.0, 0.0, 0.0, "HD"),   # AutoDock polar H
        ])
        try:
            atoms = control.read_heavy_atoms(path)
            self.assertEqual([a[0] for a in atoms], ["C1"])
        finally:
            os.unlink(path)

    def test_reads_a_real_aromatic_pdbqt(self):
        atoms = control.read_heavy_atoms(CATECHIN)
        self.assertEqual(len(atoms), 21)
        # After the 'A' fix this file must contain only real element symbols.
        self.assertEqual(sorted(set(a[4] for a in atoms)), ["C", "O"])

    def test_reads_a_real_crystal_reference_pdb(self):
        atoms = control.read_heavy_atoms(ISOMALTOSE_REF)
        self.assertEqual(len(atoms), 23)
        self.assertEqual([a[0] for a in atoms[:3]], ["C1", "C2", "C3"])

    def test_selects_the_requested_model(self):
        first = control.read_heavy_atoms(NINE_POSES, model=1)
        ninth = control.read_heavy_atoms(NINE_POSES, model=9)
        self.assertEqual(len(first), 12)
        self.assertEqual(len(ninth), 12)
        # Nine docked modes of one ligand: same atoms, different coordinates.
        self.assertNotEqual([a[1:4] for a in first], [a[1:4] for a in ninth])

    def test_malformed_lines_are_dropped_silently(self):
        # KNOWN GAP. A truncated or corrupt file loses atoms with no warning
        # and no count in the run record. The caller then sees only "atom
        # counts differ", which reads as a mismatched reference rather than a
        # damaged file. Better would be to report how many lines were skipped.
        path = _write([
            _atom_line(1, "C1", 0.0, 0.0, 0.0, "C"),
            "ATOM      2  C2  LIG A   1     NOTANUM   0.000   0.000       C",
            "ATOM      3  C3",
            _atom_line(4, "C4", 1.0, 0.0, 0.0, "C"),
        ])
        try:
            self.assertEqual(len(control.read_heavy_atoms(path)), 2)
        finally:
            os.unlink(path)


class TestRmsd(unittest.TestCase):

    def test_identical_coordinates_matched_by_name_are_zero(self):
        atoms = _atoms([("C1", 1.0, 2.0, 3.0, "C"), ("O1", 4.0, 5.0, 6.0, "O")])
        value, method, matched, exact = control.rmsd(atoms, list(atoms))
        self.assertEqual(value, 0.0)
        self.assertIn("by atom name", method)
        self.assertEqual(matched, 2)
        self.assertTrue(exact)

    def test_a_known_offset_measures_that_offset(self):
        pose = _atoms([("C1", 0.0, 0.0, 0.0, "C"), ("O1", 0.0, 0.0, 0.0, "O")])
        ref = _atoms([("C1", 3.0, 0.0, 0.0, "C"), ("O1", 0.0, 4.0, 0.0, "O")])
        value, _method, _matched, exact = control.rmsd(pose, ref)
        # sqrt((9 + 16) / 2)
        self.assertAlmostEqual(value, 3.5355339059, places=7)
        self.assertTrue(exact)

    def test_no_atoms_and_mismatched_counts_return_no_value(self):
        one = _atoms([("C1", 0.0, 0.0, 0.0, "C")])
        for pose, ref in ((one, []), ([], one), (one, one + one)):
            value, method, matched, exact = control.rmsd(pose, ref)
            self.assertIsNone(value)
            self.assertEqual(matched, 0)
            self.assertFalse(exact)
            self.assertTrue(method)

    def test_matching_element_sequence_pairs_by_file_order(self):
        pose = _atoms([("C", 0.0, 0.0, 0.0, "C"), ("O", 1.0, 0.0, 0.0, "O")])
        ref = _atoms([("C1", 0.0, 0.0, 0.0, "C"), ("O1", 1.0, 0.0, 0.0, "O")])
        value, method, _matched, exact = control.rmsd(pose, ref)
        self.assertEqual(value, 0.0)
        self.assertIn("by file order", method)
        self.assertTrue(exact)

    def test_file_order_is_refused_when_element_sequences_differ(self):
        # REGRESSION for the withdrawn 7.046 A. Same atom count, names that do
        # not match, elements in a different order. Pairing by file order here
        # is what matched pose O5 against reference C6. It must not happen,
        # and the answer must not be labelled exact.
        pose = _atoms([("C", 0.0, 0.0, 0.0, "C"),
                       ("C", 1.0, 0.0, 0.0, "C"),
                       ("O", 2.0, 0.0, 0.0, "O")])
        ref = _atoms([("C1", 0.0, 0.0, 0.0, "C"),
                      ("O1", 1.0, 0.0, 0.0, "O"),
                      ("C2", 2.0, 0.0, 0.0, "C")])
        value, method, _matched, exact = control.rmsd(pose, ref)
        self.assertNotIn("by file order", method)
        self.assertIn("LOWER BOUND", method)
        self.assertFalse(exact)
        self.assertIsNotNone(value)

    def test_unresolvable_correspondence_is_reported_as_a_lower_bound(self):
        pose = _atoms([("C", 0.0, 0.0, 0.0, "C"), ("O", 5.0, 0.0, 0.0, "O")])
        ref = _atoms([("O1", 5.5, 0.0, 0.0, "O"), ("C1", 0.5, 0.0, 0.0, "C")])
        value, method, _matched, exact = control.rmsd(pose, ref)
        self.assertAlmostEqual(value, 0.5, places=7)
        self.assertIn("LOWER BOUND", method)
        self.assertFalse(exact)

    def test_a_lower_bound_never_exceeds_a_real_assignment(self):
        # The whole point of the lower bound: it can only ever be too small.
        # Same coordinate set both sides, permuted. Matched by name it is a
        # long way out; matched to the nearest same-element atom it is zero.
        coords = [(0.0, 0.0, 0.0, "C"), (10.0, 0.0, 0.0, "C"),
                  (20.0, 0.0, 0.0, "O")]
        named_pose = _atoms([("A1",) + coords[0], ("A2",) + coords[1],
                             ("A3",) + coords[2]])
        named_ref = _atoms([("A1",) + coords[1], ("A2",) + coords[0],
                            ("A3",) + coords[2]])
        by_name, _m, _n, exact = control.rmsd(named_pose, named_ref)
        self.assertTrue(exact)

        unnamed_pose = _atoms([("C",) + coords[0], ("C",) + coords[1],
                               ("O",) + coords[2]])
        unnamed_ref = _atoms([("X1",) + coords[1], ("X2",) + coords[2],
                              ("X3",) + coords[0]])
        lower, method, _n2, exact2 = control.rmsd(unnamed_pose, unnamed_ref)
        self.assertIn("LOWER BOUND", method)
        self.assertFalse(exact2)
        self.assertLessEqual(lower, by_name)

    def test_an_aromatic_ligand_gets_a_verdict(self):
        # REGRESSION, and the reason this file exists. A real flavan-3-ol from
        # validation-inputs measured against a crystal-style reference built
        # from its own coordinates: the true answer is exactly 0.000 A. While
        # 'A' read as its own element this returned None, and the CLI printed
        # NO VERDICT with the reason "no reference atom of element A".
        pose = control.read_heavy_atoms(CATECHIN)
        carbons = oxygens = 0
        reference = []
        for (_name, x, y, z, element) in pose:
            if element == "C":
                carbons += 1
                reference.append(("C%d" % carbons, x, y, z, "C"))
            else:
                oxygens += 1
                reference.append(("O%d" % oxygens, x, y, z, "O"))
        self.assertGreater(carbons, 0)
        value, _method, matched, exact = control.rmsd(pose, reference)
        self.assertIsNotNone(value)
        self.assertEqual(value, 0.0)
        self.assertEqual(matched, len(pose))
        self.assertTrue(exact)

    def test_duplicate_pose_names_are_accepted_as_exact(self):
        # KNOWN GAP. The guard checks that the REFERENCE has unique names and
        # then counts pairs, not distinct reference atoms, so a pose whose
        # atoms share a name can consume one reference atom twice and leave
        # another never compared. Here reference C2, 9 A away, is ignored and
        # the answer comes back as an exact 0.000. Not reachable from the
        # current production files, because Meeko names every heavy atom 'C'
        # or 'O' and no crystal reference uses those names, so case 1 falls
        # through. It should still refuse rather than rely on that.
        pose = _atoms([("C1", 0.0, 0.0, 0.0, "C"), ("C1", 0.0, 0.0, 0.0, "C")])
        ref = _atoms([("C1", 0.0, 0.0, 0.0, "C"), ("C2", 9.0, 0.0, 0.0, "C")])
        value, _method, _matched, exact = control.rmsd(pose, ref)
        self.assertEqual(value, 0.0)
        self.assertTrue(exact)

    def test_element_sequence_alone_is_blind_to_constitutional_isomers(self):
        # KNOWN GAP, written down at REVIEW-FINDINGS.md:210 and still open.
        # Two different molecules that happen to list the same elements in the
        # same order pair by file order and are reported as exact. Matching
        # element sequences is not proof of correspondence.
        pose = _atoms([("C", 0.0, 0.0, 0.0, "C"), ("O", 1.0, 0.0, 0.0, "O")])
        ref = _atoms([("C1", 5.0, 0.0, 0.0, "C"), ("O1", 6.0, 0.0, 0.0, "O")])
        value, method, _matched, exact = control.rmsd(pose, ref)
        self.assertAlmostEqual(value, 5.0, places=7)
        self.assertIn("by file order", method)
        self.assertTrue(exact)


if __name__ == "__main__":
    unittest.main(verbosity=2)
