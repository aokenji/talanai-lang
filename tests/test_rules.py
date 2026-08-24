#!/usr/bin/env python3
"""
Unit tests for the rules. Standard library unittest, no pytest.

    python tests/test_rules.py

Each test builds the smallest .tal that triggers one rule, so a failure names
the rule that broke rather than "the example file changed".
"""

from __future__ import annotations

import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

from talanai import chem, model, parse, rules   # noqa: E402

BASE = """
study "test"
  claim_scope genus Ziziphus
receptor 3A4A
  prepare raw
  resolution 1.6 A
site
  center 21.52 -7.70 23.55
  size 30 30 30
  must_enclose Asp215
control
  measure redock_rmsd
  require under 2.0 A
  result 0.519 A
  source validation.json
ligands
  prepare meeko
  protonate pH 7.4
  compound Rutin C27H30O16
reference Acarbose
  formula C25H43NO18
  prepare meeko
dock
  engine vina 1.2.7
  exhaustiveness 32
  modes 9
  energy_range 3
  seeds 42 43 44
rank
  by affinity
  by ligand_efficiency
"""


def build(text):
    document = parse.parse_text(text, path="test.tal",
                                known_keys=model.KNOWN_KEYS)
    return model.Experiment(document)


def codes(experiment, level=None):
    found = rules.run_all(experiment)
    return {f.code for f in found if level is None or f.level == level}


def _atom_line(serial, name, resname, chain, resseq, x, y, z, element,
              record="ATOM"):
    """
    One fixed-column PDB atom line, filled in only as precisely as pdb.py and
    control.py actually read: record at columns 1-6, name at 13-16, resName
    at 18-20, chain at 22, resSeq at 23-26, x/y/z at 31-54, element at 77-79.
    Good enough for the R306 fixtures; not a general PDB writer.
    """
    line = list(" " * 80)
    line[0:len(record)] = list(record)
    line[6:11] = list("%5d" % serial)
    line[12:16] = list("%-4s" % name)
    line[17:20] = list("%-3s" % resname)
    line[21] = chain or " "
    line[22:26] = list("%4d" % resseq)
    line[30:38] = list("%8.3f" % x)
    line[38:46] = list("%8.3f" % y)
    line[46:54] = list("%8.3f" % z)
    line[76:79] = list("%3s" % element)
    return "".join(line)


def _write_pdb(lines):
    """A tiny synthetic PDB file under tests/, cleaned up by the caller."""
    handle = tempfile.NamedTemporaryFile(
        mode="w", suffix=".pdb", delete=False, dir=HERE, encoding="utf-8")
    handle.write("\n".join(lines) + "\n")
    handle.close()
    return handle.name


class TestRules(unittest.TestCase):

    def test_baseline_is_accepted(self):
        self.assertEqual(codes(build(BASE), rules.REFUSE), set())

    def test_r101_missing_control_result(self):
        text = BASE.replace("result 0.519 A", "result n.d.")
        self.assertIn("R101", codes(build(text), rules.REFUSE))

    def test_r102_failing_control(self):
        text = BASE.replace("result 0.519 A", "result 3.4 A")
        self.assertIn("R102", codes(build(text), rules.REFUSE))

    def test_r105_control_without_source(self):
        text = BASE.replace("  source validation.json\n", "")
        self.assertIn("R105", codes(build(text), rules.REFUSE))

    def test_r106_control_prep_differs_from_screen(self):
        text = BASE.replace("  result 0.519 A",
                            "  prepare meeko\n  result 0.519 A")
        self.assertIn("R106", codes(build(text), rules.REFUSE))

    def test_r106_matching_prep_passes(self):
        text = BASE.replace("  result 0.519 A",
                            "  prepare raw\n  result 0.519 A")
        self.assertNotIn("R106", codes(build(text), rules.REFUSE))

    def test_r104_control_searched_harder(self):
        text = BASE.replace("  result 0.519 A",
                            "  exhaustiveness 64\n  result 0.519 A")
        self.assertIn("R104", codes(build(text), rules.WARN))

    def test_r201_bare_rmsd(self):
        text = BASE.replace("measure redock_rmsd", "measure rmsd")
        self.assertIn("R201", codes(build(text), rules.REFUSE))

    def test_r301_incomplete_box(self):
        text = BASE.replace("size 30 30 30", "size 30 30")
        self.assertIn("R301", codes(build(text), rules.REFUSE))

    def test_r401_low_exhaustiveness(self):
        text = BASE.replace("exhaustiveness 32", "exhaustiveness 8")
        self.assertIn("R401", codes(build(text), rules.WARN))

    def test_r501_single_seed(self):
        text = BASE.replace("seeds 42 43 44", "seeds 42")
        self.assertIn("R501", codes(build(text), rules.WARN))

    def test_r601_mismatched_preparation(self):
        text = BASE.replace("reference Acarbose\n  formula C25H43NO18\n  prepare meeko",
                            "reference Acarbose\n  formula C25H43NO18\n  prepare raw")
        self.assertIn("R601", codes(build(text), rules.REFUSE))

    def test_r602_receptor_prep_missing(self):
        text = BASE.replace("receptor 3A4A\n  prepare raw", "receptor 3A4A")
        self.assertIn("R602", codes(build(text), rules.REFUSE))

    def test_r701_ranking_without_efficiency(self):
        text = BASE.replace("  by ligand_efficiency\n", "")
        self.assertIn("R701", codes(build(text), rules.WARN))

    def test_r801_species_scope_refused(self):
        text = BASE.replace("claim_scope genus Ziziphus",
                            "claim_scope Ziziphus talanai")
        self.assertIn("R801", codes(build(text), rules.REFUSE))

    def test_r801_missing_scope_refused(self):
        text = BASE.replace("  claim_scope genus Ziziphus\n", "")
        self.assertIn("R801", codes(build(text), rules.REFUSE))

    def test_r803_overclaiming_verb(self):
        text = BASE.replace('study "test"', 'study "Rutin inhibits alpha-glucosidase"')
        self.assertIn("R803", codes(build(text), rules.WARN))

    def test_r804_library_species_scope(self):
        text = "library thesis-five\n  scope Ziziphus talanai\n"
        self.assertIn("R804", codes(build(text), rules.REFUSE))

    def test_r804_library_genus_scope_passes(self):
        text = "library thesis-five\n  scope genus Ziziphus\n"
        self.assertNotIn("R804", codes(build(text), rules.REFUSE))

    def test_r004_duplicate_single_value_key(self):
        text = BASE.replace("exhaustiveness 32", "exhaustiveness 32\n  exhaustiveness 8")
        self.assertIn("R004", codes(build(text), rules.REFUSE))

    def test_unknown_key_warns_not_ignored(self):
        text = BASE.replace("  modes 9", "  modez 9")
        self.assertIn("R003", codes(build(text), rules.WARN))


class TestR306PoseContactsSite(unittest.TestCase):
    """
    R303 checks that the search BOX encloses the catalytic residues. R306
    checks that the rank-1 POSE actually comes near them, using a synthetic
    receptor and pose file built with tempfile so the test does not depend
    on any real structure on disk.
    """

    def setUp(self):
        self.cleanup = []

    def tearDown(self):
        for path in self.cleanup:
            try:
                os.remove(path)
            except OSError:
                pass

    def _receptor(self):
        path = _write_pdb([
            _atom_line(1, "CA", "ASP", "A", 215, 0.0, 0.0, 0.0, "C"),
        ])
        self.cleanup.append(path)
        return path

    def _pose(self, x, y, z):
        path = _write_pdb([
            _atom_line(1, "C1", "LIG", "A", 1, x, y, z, "C", record="HETATM"),
        ])
        self.cleanup.append(path)
        return path

    def _with_receptor_file(self, receptor_path):
        return BASE.replace("receptor 3A4A\n  prepare raw",
                            "receptor 3A4A\n  prepare raw\n  file %s"
                            % receptor_path)

    def test_no_pose_file_is_unverified(self):
        text = self._with_receptor_file(self._receptor())
        found = [f for f in rules.run_all(build(text)) if f.code == "R306"]
        self.assertEqual(len(found), 1)
        self.assertEqual(found[0].level, rules.RECORD)
        self.assertIn("UNVERIFIED", found[0].title)

    def test_pose_contacting_site_passes(self):
        # 1.0 A from the residue atom at the origin: inside the 4.0 A cutoff.
        text = self._with_receptor_file(self._receptor())
        text += "\nresults\n  pose_file %s\n" % self._pose(1.0, 0.0, 0.0)
        found = [f for f in rules.run_all(build(text)) if f.code == "R306"]
        self.assertEqual(len(found), 1)
        self.assertEqual(found[0].level, rules.PASS)

    def test_pose_far_from_site_warns(self):
        # 20 A from the residue atom at the origin: well past the cutoff.
        text = self._with_receptor_file(self._receptor())
        text += "\nresults\n  pose_file %s\n" % self._pose(20.0, 0.0, 0.0)
        found = [f for f in rules.run_all(build(text)) if f.code == "R306"]
        self.assertEqual(len(found), 1)
        self.assertEqual(found[0].level, rules.WARN)


class TestChemistry(unittest.TestCase):

    def test_heavy_atoms(self):
        self.assertEqual(chem.heavy_atoms("C27H30O16"), 43)
        self.assertEqual(chem.heavy_atoms("C15H10O7"), 22)
        self.assertEqual(chem.heavy_atoms("C25H43NO18"), 44)

    def test_unicode_subscripts(self):
        self.assertEqual(chem.heavy_atoms("C₂₇H₃₀O₁₆"), 43)

    def test_ligand_efficiency_inverts_the_ranking(self):
        rutin = chem.ligand_efficiency(-8.857, "C27H30O16")
        kaempferol = chem.ligand_efficiency(-7.479, "C15H10O6")
        self.assertLess(-7.479, -8.857 + 1.5)         # rutin scores stronger
        self.assertGreater(kaempferol, rutin)          # kaempferol is efficient

    def test_species_detection(self):
        self.assertTrue(chem.looks_like_species("Ziziphus talanai"))
        self.assertTrue(chem.looks_like_species("Z. mauritiana"))
        self.assertFalse(chem.looks_like_species("genus Ziziphus"))
        self.assertFalse(chem.looks_like_species("Ziziphus"))


class TestParser(unittest.TestCase):

    def test_repeated_keys_accumulate(self):
        experiment = build(BASE + "\nresults\n  affinity Rutin -8.857\n"
                                  "  affinity Quercetin -7.503\n")
        self.assertEqual(len(experiment.results.all("affinity")), 2)

    def test_comments_and_blank_lines_ignored(self):
        experiment = build("study \"x\"  # trailing comment\n\n  author Ken\n")
        self.assertEqual(experiment.study.one("author"), "Ken")

    def test_wrapped_value_continues(self):
        experiment = build('study "x"\n  claim_scope genus Ziziphus\n'
                           'receptor 3A4A\n  prepare_note first line\n'
                           '               second line\n')
        self.assertIn("second line", experiment.receptor.one("prepare_note"))


class TestR307RingGeometry(unittest.TestCase):
    """
    R307 refuses a ligand whose saturated six-rings are not chairs, because
    Vina holds rings rigid and cannot repair one. It must report UNVERIFIED
    rather than passing when it has nothing to read.
    """

    def test_unverified_without_a_ligand_directory(self):
        # BASE declares no ligands.dir, so there is nothing to check.
        found = rules.run_all(build(BASE))
        r307 = [f for f in found if f.code == "R307"]
        self.assertTrue(r307, "R307 should always report something")
        self.assertEqual(r307[0].level, rules.RECORD)
        self.assertIn("UNVERIFIED", r307[0].title)

    def test_unverified_when_directory_is_missing(self):
        text = BASE.replace("  prepare meeko",
                            "  dir no_such_directory_here\n  prepare meeko")
        found = rules.run_all(build(text))
        r307 = [f for f in found if f.code == "R307"]
        self.assertEqual(r307[0].level, rules.RECORD)
        self.assertIn("UNVERIFIED", r307[0].title)

    def test_never_reports_a_pass_it_did_not_earn(self):
        """The whole convention: no readable ligand means no pass."""
        for text in (BASE,
                     BASE.replace("  prepare meeko",
                                  "  dir nowhere\n  prepare meeko")):
            found = rules.run_all(build(text))
            r307 = [f for f in found if f.code == "R307"]
            self.assertNotEqual(r307[0].level, rules.PASS)

    def test_chair_detection_on_synthetic_geometry(self):
        """
        Cremer-Pople theta must separate a chair from a boat. Built from
        idealised cyclohexane coordinates rather than a real file, so the
        maths is tested independently of any toolkit or parser.
        """
        from talanai.ringcheck import _cremer_pople_theta
        # Ideal chair: alternating +/- 0.25 A out of plane.
        import math as _math
        chair, boat = [], []
        for i in range(6):
            angle = 2 * _math.pi * i / 6
            x, y = _math.cos(angle) * 1.46, _math.sin(angle) * 1.46
            chair.append([x, y, 0.25 if i % 2 == 0 else -0.25])
            # Boat: two atoms up on the same side instead of alternating.
            boat.append([x, y, 0.25 if i in (0, 3) else -0.25 if i in (1, 4) else 0.0])
        theta_chair = _cremer_pople_theta(chair)
        theta_boat = _cremer_pople_theta(boat)
        self.assertTrue(theta_chair <= 35.0 or theta_chair >= 145.0,
                        "ideal chair should read as a chair, got %.1f" % theta_chair)
        self.assertNotAlmostEqual(theta_chair, theta_boat, places=0,
                                  msg="chair and boat must not measure the same")


if __name__ == "__main__":
    unittest.main(verbosity=2)
