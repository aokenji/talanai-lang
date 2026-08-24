"""
The rules. This file is the product.

Everything else in Talanai is plumbing that exists so this file can run. Each
rule below refuses, warns about, or records one way a docking experiment can
fail to mean what it claims to mean. Every one of them has a reason written in
biology, not in software, because the reason is the part a scientist has to
agree with.

SEVERITIES
    REFUSE  the experiment does not run. Not a preference; the result would
            not be interpretable.
    WARN    it runs, and you should know this before you present it.
    RECORD  neither wrong nor optional: a fact that must appear in the run
            record and the methods paragraph.

A check that cannot be performed reports RECORD with the word UNVERIFIED. It
never reports a pass it did not earn.
"""

from __future__ import annotations

import os
import re

from . import chem, control, model, pdb

REFUSE, WARN, RECORD, PASS = "REFUSE", "WARN", "RECORD", "pass"

SEVERITY_ORDER = {REFUSE: 0, WARN: 1, RECORD: 2, PASS: 3}

# --------------------------------------------------------------------------
# Project conventions, not published specifications.
#
# AutoDock Vina documents no scoring-precision figure. Its own accuracy
# against experimental binding affinities is reported at 2 to 3 kcal/mol,
# roughly two orders of magnitude looser than either constant below. Both
# numbers here are judgement calls this project makes for its own purposes,
# recorded so they can be cited accurately and argued with, not mistaken for
# something Vina publishes.
# --------------------------------------------------------------------------

# RA02: how tightly a same-pipeline, same-input RERUN should reproduce its own
# earlier score. This is run-to-run numerical reproducibility (did anything in
# the pipeline change between two runs of itself), which is a different
# quantity from scoring accuracy (how close a Vina score is to a measured
# binding affinity). This project previously called this threshold "scoring
# precision", which implied the latter and is not something Vina publishes.
REPLICATION_REPRODUCIBILITY_LIMIT_KCAL_MOL = 0.1     # project convention

# R307: how far a saturated six-ring's Cremer-Pople theta may sit from a chair
# (0 or 180 degrees) before the ring is called misfolded. 35 degrees is
# generous: a real chair lands within about 15 degrees, and a boat or
# twist-boat sits near 90. The gap between those is wide, so the threshold does
# not need to be delicate.
#
# 2026-08-04: this rule exists because the study's own ligand preparation
# embedded ONE conformer and minimised it locally. Local minimisation cannot
# cross a ring-flip barrier, Meeko writes rings rigid and AutoDock Vina cannot
# alter ring geometry, so whichever pucker that single attempt produced was
# frozen for the whole study. Three ligands were affected, including the
# reference compound, and correcting the reference alone moved another
# compound from "loses to acarbose" to "tied with it".
RING_CHAIR_MAX_THETA_DEG = 35.0                      # project convention

# R306: the heavy-atom distance, in angstrom, under which a docked pose counts
# as "engaging" a named residue rather than merely sharing a box with it. 4.0 A
# is a conventional non-bonded contact cutoff in structural biology, covering
# van der Waals contact and the roughly 2.5 to 3.5 A donor-acceptor spacing of
# a hydrogen bond; it is the same order of magnitude PLIP and similar
# contact-detection tools use for reporting a residue as contacted. It is a
# judgement call, not a PDB or Vina specification, chosen loose enough to
# admit an imperfect but real contact and tight enough to exclude "same
# pocket, not touching."
POSE_CONTACT_DISTANCE_A = 4.0                        # project convention


class Finding:
    def __init__(self, level, code, title, detail="", fix="", line=None):
        self.level = level
        self.code = code
        self.title = title
        self.detail = detail
        self.fix = fix
        self.line = line

    def __repr__(self):
        return "<%s %s %s>" % (self.level, self.code, self.title)


RULES = []


def rule(function=None, scope="experiment"):
    """Register a rule. scope limits it to a kind of file: an experiment, a
    ziziphus/ library, or 'any' for rules about the file format itself."""
    def wrap(target):
        target.scope = scope
        RULES.append(target)
        return target
    return wrap(function) if function is not None else wrap


def run_all(experiment):
    kind = experiment.kind
    findings = []
    for problem_line, message in experiment.document.problems:
        findings.append(Finding(REFUSE, "R000", "the file could not be read",
                                message, "every setting sits under a block",
                                problem_line))
    for function in RULES:
        if getattr(function, "scope", "experiment") not in ("any", kind):
            continue
        try:
            findings.extend(function(experiment) or [])
        except Exception as error:                       # a rule must never
            findings.append(Finding(                     # take the tool down
                WARN, "R999", "rule %s could not run" % function.__name__,
                "%s: %s" % (type(error).__name__, error),
                "this is a bug in Talanai, not in the experiment"))
    return findings


def docking_is_locked(findings):
    return any(f.level == REFUSE for f in findings)


# ==========================================================================
# R0xx  the file itself
# ==========================================================================
@rule
def r001_required_blocks(exp):
    """An experiment that omits any of these is not described, only sketched."""
    out = []
    for name in model.REQUIRED_BLOCKS:
        if exp.block(name) is None:
            out.append(Finding(
                REFUSE, "R001", "missing '%s' block" % name,
                "an experiment is not described without it",
                "add a '%s' block" % name))
    return out


@rule(scope="any")
def r002_unknown_blocks(exp):
    """A misspelt block would otherwise be silently ignored, which is worse."""
    out = []
    for block in exp.document.blocks:
        if block.name not in model.SCHEMA:
            out.append(Finding(
                WARN, "R002", "unknown block '%s'" % block.name,
                "line %d. Talanai ignores it." % block.line,
                "the vocabulary is: %s" % ", ".join(model.EXPERIMENT_BLOCKS),
                block.line))
    return out


@rule(scope="any")
def r003_unknown_keys(exp):
    """Same reasoning one level down: a typo must not read as a default."""
    out = []
    for block in exp.document.blocks:
        known = model.SCHEMA.get(block.name)
        if known is None:
            continue
        for key in block.keys:
            if key not in known:
                out.append(Finding(
                    WARN, "R003", "unknown key '%s' in '%s'" % (key, block.name),
                    "line %d. Talanai ignores it." % block.line_of(key),
                    "known keys here: %s" % ", ".join(sorted(known)) or "none",
                    block.line_of(key)))
    return out


@rule(scope="any")
def r004_repeated_single_value(exp):
    """Two values for one setting means nobody knows which one ran."""
    out = []
    for block in exp.document.blocks:
        for key, values in block.keys.items():
            if key in model.SINGLE_VALUE_KEYS and len(values) > 1:
                out.append(Finding(
                    REFUSE, "R004", "'%s' is set more than once in '%s'"
                    % (key, block.name),
                    "values: %s" % " | ".join(values),
                    "keep one. Which one ran is not a detail.",
                    block.line_of(key)))
    return out


# ==========================================================================
# R1xx  the control, which unlocks everything else
# ==========================================================================
@rule
def r101_control_recorded(exp):
    """
    A protocol that cannot reproduce a known answer has no business
    predicting an unknown one. Redocking the co-crystallised ligand is the
    positive control for docking, exactly as a known-positive lane is the
    control for a blot. Until it is recorded, docking stays locked.
    """
    if exp.control is None:
        return []
    if exp.control_result() is not None:
        return []
    raw = exp.control.one("result") or "empty"
    return [Finding(
        REFUSE, "R101", "the redocking control has no recorded result",
        "control.result is '%s'. Nothing downstream can be trusted yet." % raw,
        "run the redock, store the run record, then put the number here",
        exp.control.line_of("result"))]


@rule
def r102_control_passes(exp):
    """A control that fails is information, not an obstacle. Read it."""
    got, limit = exp.control_result(), exp.control_limit()
    if got is None or limit is None:
        return []
    if got >= limit:
        return [Finding(
            REFUSE, "R102", "the redocking control failed",
            "redock_rmsd %.2f A is not under %.2f A" % (got, limit),
            "the pocket, the preparation or the search is wrong. Fix that "
            "before screening anything.",
            exp.control.line_of("result"))]
    return [Finding(
        PASS, "R102", "control passed",
        "redock_rmsd %.2f A, limit %.2f A" % (got, limit))]


@rule
def r103_control_threshold_conventional(exp):
    """2.0 A is the accepted bar. A looser one has to be argued for, in text."""
    limit = exp.control_limit()
    if limit is None or limit <= model.ACCEPTED_REDOCK_LIMIT:
        return []
    return [Finding(
        WARN, "R103", "the control threshold is looser than convention",
        "require is %.2f A; the field's accepted bar is %.1f A"
        % (limit, model.ACCEPTED_REDOCK_LIMIT),
        "tighten it, or state in the methods why this pocket needs more room")]


@rule
def r106_control_uses_the_screening_receptor(exp):
    """
    THE ONE THAT BITES.

    A control validates the protocol it was run under, not a nearby one. If
    the redock used a differently prepared receptor than the screen, then the
    thing that passed at 2.0 A is not the thing that produced the affinities.

    This is rule R601 applied where it matters most, and it is the same lesson
    this project already paid for in June 2026: a Vina score belongs to a pair
    of prepared files, so changing the receptor preparation changes the
    measurement, including the measurement of whether the protocol works.
    """
    if exp.control is None or exp.receptor is None:
        return []
    control_prep = (exp.control.one("prepare") or "").strip().lower()
    screen_prep = (exp.receptor.one("prepare") or "").strip().lower()
    if not (control_prep and screen_prep):
        return []
    if control_prep == screen_prep:
        return [Finding(PASS, "R106",
                        "control and screen share a receptor preparation",
                        screen_prep)]
    return [Finding(
        REFUSE, "R106", "the control validated a different receptor preparation",
        "control ran on: %s\n            the screen ran on: %s\n"
        "            The redock therefore validates the first protocol, not "
        "the one that produced these affinities."
        % (control_prep, screen_prep),
        "redock on the screening receptor and record that number, or state "
        "plainly that the screening path is unvalidated",
        exp.control.line_of("prepare"))]


@rule
def r104_control_matches_screening_search(exp):
    """
    Redocking in a tight box at high exhaustiveness is standard and sensible:
    it asks whether the scoring function can recover a known pose. But it is
    an easier question than the screen asks. When the two configurations
    differ, the record has to say so, or a reader will reasonably assume the
    reported RMSD validates the screening search.
    """
    if exp.control is None:
        return []
    out = []
    _center, size = exp.box()
    control_size = chem.numbers(exp.control.one("box_size", ""))
    if size and control_size:
        if max(control_size) != max(size):
            out.append(Finding(
                WARN, "R104", "the control used a different box than the screen",
                "control %.0f A, screen %.0f A per side"
                % (max(control_size), max(size)),
                "legitimate and common, but state it in the methods: the "
                "control validates pose recovery, not the screening search",
                exp.control.line_of("box_size")))
    control_exhaustiveness = chem.as_float(exp.control.one("exhaustiveness"))
    screen_exhaustiveness = exp.exhaustiveness()
    if control_exhaustiveness and screen_exhaustiveness:
        if control_exhaustiveness > screen_exhaustiveness:
            out.append(Finding(
                WARN, "R104", "the control searched harder than the screen",
                "control exhaustiveness %.0f, screen %.0f"
                % (control_exhaustiveness, screen_exhaustiveness),
                "the control passed under conditions the screen does not have",
                exp.control.line_of("exhaustiveness")))
    return out


@rule
def r105_control_has_a_source(exp):
    """
    A number with no stored artefact behind it is a claim, not a measurement.
    This is the distinction protocol.js already makes when it holds
    redockRmsd at null despite a live endpoint reporting 0.519.
    """
    if exp.control is None or exp.control_result() is None:
        return []
    if (exp.control.one("source") or "").strip():
        return [Finding(PASS, "R105", "control result is traceable",
                        exp.control.one("source"))]
    return [Finding(
        REFUSE, "R105", "the control result has no source",
        "a live endpoint reading is not a citable artefact. Someone else has "
        "to be able to find the run this number came from.",
        "add 'source' pointing at a stored run record, log or DOI",
        exp.control.line)]


# ==========================================================================
# R2xx  words that mean two things
# ==========================================================================
@rule
def r201_bare_rmsd(exp):
    """
    This project uses RMSD for two unrelated quantities: how close a redocked
    pose lands to the crystal pose, and how far apart the poses in one
    compound's cluster sit (1.61 to 2.00 A). Writing 'rmsd' alone invites
    someone to quote the second where the first was meant, which would be a
    real scientific error rather than a typo.
    """
    if exp.control is None:
        return []
    measure = (exp.control.one("measure") or "").strip()
    if measure.lower() == "rmsd":
        return [Finding(
            REFUSE, "R201", "'rmsd' is ambiguous in this study",
            "it names two different quantities here: redocking accuracy, and "
            "per-compound pose-cluster spread.",
            "write 'redock_rmsd' or 'pose_cluster_rmsd'",
            exp.control.line_of("measure"))]
    if measure:
        return [Finding(PASS, "R202", "control measure named unambiguously",
                        measure)]
    return []


# ==========================================================================
# R3xx  geometry
# ==========================================================================
@rule
def r301_box_specified(exp):
    """A box is three coordinates and three lengths. Anything else is a guess."""
    center, size = exp.box()
    if center and size:
        return []
    return [Finding(
        REFUSE, "R301", "the box is not fully specified",
        "need three numbers for center and three for size",
        "center  21.52 -7.70 23.55   /   size  30 30 30",
        exp.site.line if exp.site else None)]


@rule
def r302_box_declares_residues(exp):
    """
    An unanchored box is the commonest way to produce confident nonsense: the
    search runs perfectly, in the wrong pocket.
    """
    if exp.site is None:
        return []
    if (exp.site.one("must_enclose") or "").strip():
        return []
    return [Finding(
        WARN, "R302", "the box declares no residues to enclose",
        "nothing ties these coordinates to the catalytic site",
        "must_enclose Asp215 Glu277 Asp352", exp.site.line)]


@rule
def r303_box_encloses_site(exp):
    """
    The real geometric check, when the structure is on disk. When it is not,
    this reports UNVERIFIED rather than passing, because an unperformed check
    is not a passed one.
    """
    center, size = exp.box()
    if not (center and size) or exp.site is None:
        return []
    residues = (exp.site.one("must_enclose") or "").split()
    if not residues:
        return []

    base = os.path.dirname(os.path.abspath(exp.path))
    structure = pdb.resolve(exp.receptor.one("file") if exp.receptor else None, base)
    if not structure:
        return [Finding(
            RECORD, "R303", "box encloses the catalytic site: UNVERIFIED",
            "declared residues %s. No receptor file on disk to check against."
            % " ".join(residues),
            "add 'file 3A4A.pdb' to the receptor block to make this a real check")]

    out, checked, missing = [], 0, []
    for token in residues:
        parsed = pdb.parse_residue_token(token)
        if not parsed:
            continue
        resname, resseq, chain = parsed
        centroid = pdb.residue_centroid(structure, resname, resseq, chain)
        if centroid is None:
            missing.append(token)
            continue
        checked += 1
        if not pdb.inside_box(centroid, center, size):
            overshoot = pdb.distance_outside(centroid, center, size)
            out.append(Finding(
                REFUSE, "R303", "the box does not enclose %s" % token,
                "residue centroid lies outside by %.1f, %.1f, %.1f A"
                % tuple(overshoot),
                "move the box centre or enlarge it until the site is inside"))
    if missing:
        out.append(Finding(
            WARN, "R304", "residue not found in the structure",
            "%s not present in %s" % (", ".join(missing),
                                      os.path.basename(structure)),
            "check the numbering and the chain"))
    if not out and checked:
        out.append(Finding(
            PASS, "R303", "box encloses the catalytic site",
            "%d of %d residues verified against %s"
            % (checked, len(residues), os.path.basename(structure))))
    return out


@rule
def r305_blind_docking_declared(exp):
    """
    A very large box is a blind docking run. That is a legitimate method with
    a different interpretation and a much larger search requirement, so it has
    to be stated rather than inferred from the numbers.
    """
    _center, size = exp.box()
    if not size or exp.site is None:
        return []
    if max(size) <= 30.0:
        return []
    if (exp.site.one("blind") or "").strip().lower() in ("yes", "true"):
        return [Finding(RECORD, "R305", "blind docking declared",
                        "box %.0f x %.0f x %.0f A" % tuple(size))]
    return [Finding(
        WARN, "R305", "box is large enough to be blind docking, undeclared",
        "box %.0f x %.0f x %.0f A" % tuple(size),
        "add 'blind yes' if that is the intent, or tighten the box to the pocket")]


@rule
def r306_pose_contacts_site(exp):
    """
    R303 checks that the BOX encloses the catalytic residues. It says
    nothing about whether the POSE that came out of the search actually
    sits near them, which is exactly the failure mode this rule set already
    names in r302's docstring: 'the search runs perfectly, in the wrong
    pocket'. A screen can pass every other rule while every reported pose
    sits in a non-catalytic surface pocket that merely scores well; this
    project's own review measured that a decoy pocket roughly 12 A from the
    catalytic triad scored within 0.02 kcal/mol of the true site for
    glucose, on this exact receptor.

    This is WARN, not REFUSE. A real allosteric result is legitimate, and
    landing in a pocket other than the named one is not automatically
    wrong. But a pose that never comes near the named catalytic residues
    does not support an active-site inhibition claim, and the finding has
    to say so in words rather than let the affinity imply it.

    Reports RECORD with UNVERIFIED when there is no pose file to check,
    exactly as R303 does when there is no receptor structure on disk. It
    never reports a pass it did not earn.
    """
    if exp.site is None:
        return []
    residues = (exp.site.one("must_enclose") or "").split()
    if not residues:
        return []

    base = os.path.dirname(os.path.abspath(exp.path))
    pose_hint = exp.results.one("pose_file") if exp.results else None
    pose_path = pdb.resolve(pose_hint, base)
    if not pose_path:
        return [Finding(
            RECORD, "R306", "pose contacts the catalytic site: UNVERIFIED",
            "no pose file on disk to check against%s."
            % (" ('%s' does not resolve to a file)" % pose_hint if pose_hint
               else " (results.pose_file is not set)"),
            "add 'pose_file' to the results block, pointing at the rank-1 "
            "pose (PDB or PDBQT), to make this a real check")]

    structure = pdb.resolve(exp.receptor.one("file") if exp.receptor else None,
                            base)
    if not structure:
        return [Finding(
            RECORD, "R306", "pose contacts the catalytic site: UNVERIFIED",
            "declared residues %s. No receptor file on disk to check the "
            "pose against." % " ".join(residues),
            "add 'file 3A4A.pdb' to the receptor block to make this a real "
            "check")]

    pose_atoms = control.read_heavy_atoms(pose_path)
    if not pose_atoms:
        return [Finding(
            RECORD, "R306", "pose contacts the catalytic site: UNVERIFIED",
            "'%s' has no readable heavy atoms" % os.path.basename(pose_path),
            "check that the pose file is a valid PDB or PDBQT")]

    site_atoms = list(pdb.read_atoms(structure))
    contacts = []
    for token in residues:
        parsed = pdb.parse_residue_token(token)
        if not parsed:
            continue
        resname, resseq, chain = parsed
        points = [(x, y, z) for _rec, name, ch, seq, x, y, z in site_atoms
                 if seq == resseq and name == resname
                 and (not chain or ch == chain)]
        if not points:
            continue                       # R304 already reports this
        nearest = min(
            ((px - sx) ** 2 + (py - sy) ** 2 + (pz - sz) ** 2) ** 0.5
            for _name, px, py, pz, _element in pose_atoms
            for sx, sy, sz in points)
        contacts.append((token, nearest))

    if not contacts:
        return []

    closest_residue, closest_distance = min(contacts, key=lambda c: c[1])
    detail = ", ".join("%s %.2f A" % (t, d) for t, d in contacts)
    if closest_distance <= POSE_CONTACT_DISTANCE_A:
        return [Finding(
            PASS, "R306", "pose contacts the catalytic site",
            "nearest is %s at %.2f A (contact cutoff %.1f A). %s"
            % (closest_residue, closest_distance, POSE_CONTACT_DISTANCE_A,
               detail))]
    return [Finding(
        WARN, "R306", "the pose does not contact the catalytic site",
        "nearest named residue is %s at %.2f A, past the %.1f A contact "
        "cutoff (%s). A pose that does not reach the catalytic residues "
        "does not support an active-site inhibition claim, whatever the "
        "merit of an allosteric result."
        % (closest_residue, closest_distance, POSE_CONTACT_DISTANCE_A, detail),
        "confirm the pocket is intentional and describe the result as "
        "allosteric, or re-check the box placement and re-dock")]


@rule
def r307_ligand_rings_are_chairs(exp):
    """
    A saturated six-ring that is not a chair was built wrong, and docking
    cannot repair it.

    AutoDock Vina treats ring bonds as rigid. Whatever ring conformation the
    ligand preparation produced is the conformation that docks, scores and
    ends up in the results. A boat or twist-boat is a strained, usually
    non-physiological shape, and the pose it produces is not the pose the real
    molecule would adopt.

    Rings containing an sp2 atom are exempt: a double bond forces a half chair
    or sofa, which is correct rather than defective.

    Checked against the prepared ligand files when `ligands.dir` points at
    them. UNVERIFIED when they are not on disk, never a silent pass.
    """
    if exp.ligands is None:
        return []
    base = os.path.dirname(os.path.abspath(exp.path))
    directory = exp.ligands.one("dir")
    if directory and not os.path.isabs(directory):
        directory = os.path.join(base, directory)
    if not directory or not os.path.isdir(directory):
        return [Finding(
            RECORD, "R307", "ligand ring geometry: UNVERIFIED",
            "no prepared ligand directory on disk to check",
            "add 'dir' to the ligands block so the rings can be verified")]

    try:
        from . import ringcheck
    except ImportError as error:
        return [Finding(RECORD, "R307", "ligand ring geometry: UNVERIFIED",
                        "ring check unavailable: %s" % error)]

    bad, checked, unreadable = [], 0, 0
    for compound in exp.compounds():
        path = ringcheck.find_ligand_file(directory, compound)
        if not path:
            continue
        result = ringcheck.saturated_ring_report(path, RING_CHAIR_MAX_THETA_DEG)
        if result is None:
            unreadable += 1
            continue
        checked += 1
        if result["defective"]:
            bad.append((compound.display_name(), result))

    if not checked:
        return [Finding(
            RECORD, "R307", "ligand ring geometry: UNVERIFIED",
            "no ligand file in %s could be matched and read" % directory,
            "check that prepared files are named after the compounds")]

    if not bad:
        return [Finding(PASS, "R307", "every saturated ligand ring is a chair",
                        "%d ligand(s) checked" % checked)]

    detail = "\n            ".join(
        "%s: %d of %d saturated ring(s), theta %s"
        % (name, len(r["defective"]), r["saturated"],
           ", ".join("%.0f" % t for t in r["defective"]))
        for name, r in bad)
    return [Finding(
        REFUSE, "R307", "ligand ring(s) are not chairs",
        detail,
        "Vina holds rings rigid, so this cannot be fixed during docking. "
        "Re-prepare by generating many conformers and keeping the lowest "
        "energy one WHOSE SATURATED RINGS ARE ALL CHAIRS. Selecting on energy "
        "alone picks a folded, strained conformer for large flexible ligands.",
        exp.ligands.line_of("dir"))]


# ==========================================================================
# R4xx  the search
# ==========================================================================
@rule
def r401_search_budget(exp):
    """
    Vina's default exhaustiveness of 8 was chosen for small boxes. Search
    effort has to scale with the volume being searched or the result is
    whatever the sampler happened to find. Rule of thumb, not a standard.
    """
    volume, exhaustiveness = exp.box_volume(), exp.exhaustiveness()
    if not volume or not exhaustiveness:
        return []
    suggested = max(8, int(round(volume / 1000.0)))
    if exhaustiveness >= suggested:
        return [Finding(PASS, "R401", "search budget suits the box",
                        "%.0f cubic A at exhaustiveness %.0f"
                        % (volume, exhaustiveness))]
    _center, size = exp.box()
    return [Finding(
        WARN, "R401", "search budget is low for this box",
        "box %.0f x %.0f x %.0f A = %.0f cubic A, exhaustiveness %.0f"
        % (size[0], size[1], size[2], volume, exhaustiveness),
        "rule of thumb: near %d for this volume, or shrink the box. Either is "
        "defensible; not knowing which you chose is not." % suggested,
        exp.dock.line_of("exhaustiveness") if exp.dock else None)]


@rule
def r402_modes_recorded(exp):
    """num_modes and energy_range change what you are allowed to say about
    the pose ensemble, so they belong in the record."""
    if exp.dock is None:
        return []
    out = []
    for key, label in (("modes", "num_modes"), ("energy_range", "energy_range")):
        if not exp.dock.has(key):
            out.append(Finding(
                RECORD, "R402", "%s not recorded" % label,
                "the pose ensemble cannot be described without it",
                "add '%s' to the dock block" % key))
    return out


# ==========================================================================
# R5xx  the search is stochastic
# ==========================================================================
@rule
def r501_single_seed(exp):
    """
    Vina searches; where it starts affects where it lands. Repeating one seed
    reproduces the answer exactly and tells you nothing about whether the
    search converged. A fixed seed is a scale that reads the same number every
    time and can still be off.
    """
    seeds = exp.seeds()
    if len(seeds) > 1:
        return [Finding(PASS, "R501", "multiple seeds",
                        "%d seeds: %s" % (len(seeds),
                                          " ".join(str(s) for s in seeds)))]
    return [Finding(
        WARN, "R501", "single seed",
        "seed %s only. That is reproducibility, not convergence."
        % (seeds[0] if seeds else "unset"),
        "add two more seeds. Agreement within about %.1f kcal/mol (this "
        "project's convention for run-to-run reproducibility, not a "
        "published Vina precision) suggests the number is solid; a wide "
        "scatter means exhaustiveness is too low."
        % REPLICATION_REPRODUCIBILITY_LIMIT_KCAL_MOL,
        exp.dock.line_of("seeds") if exp.dock else None)]


# ==========================================================================
# R6xx  preparation is part of the measurement
# ==========================================================================
@rule
def r601_shared_preparation(exp):
    """
    June 2026, this project: TalanaiDock only reproduced the thesis value for
    Quercetin (-7.525) once Meeko preparation was removed from the receptor.
    Same compound, same protein, same engine, different preparation, different
    number. A Vina score belongs to a pair of prepared files, not to a
    molecule, so two preparations are two assays and must not share a column.
    """
    if not (exp.ligands and exp.reference):
        return []
    ligand_prep = (exp.ligands.one("prepare") or "").strip().lower()
    reference_prep = (exp.reference.one("prepare") or "").strip().lower()
    if not (ligand_prep and reference_prep):
        return []
    if ligand_prep != reference_prep:
        return [Finding(
            REFUSE, "R601", "ligands and reference were prepared differently",
            "ligands:   %s\n            reference: %s"
            % (ligand_prep, reference_prep),
            "prepare both the same way, or do not compare them",
            exp.reference.line_of("prepare"))]
    return [Finding(PASS, "R601", "ligands and reference share a preparation",
                    ligand_prep)]


@rule
def r602_receptor_prep_declared(exp):
    """An unrecorded preparation makes the whole run unrepeatable, and it is
    the single most commonly omitted item in published docking methods."""
    if exp.receptor is None:
        return []
    if (exp.receptor.one("prepare") or "").strip():
        return []
    return [Finding(
        REFUSE, "R602", "receptor preparation is not declared",
        "waters, hetatms, hydrogens and charges all change the score",
        "add 'prepare' to the receptor block, even if the answer is 'raw'",
        exp.receptor.line)]


@rule
def r603_ligand_protonation(exp):
    """
    The structure drawn in a paper is not what exists at pH 7.4. Docking the
    neutral form of an acid, or the wrong tautomer, changes the pose and the
    score, and it is almost never stated.
    """
    if exp.ligands is None:
        return []
    if (exp.ligands.one("protonate") or "").strip():
        return []
    return [Finding(
        WARN, "R603", "ligand protonation state not declared",
        "the form that exists in solution is not always the form in the file",
        "add 'protonate pH 7.4'", exp.ligands.line)]


# ==========================================================================
# R7xx  reading the results
# ==========================================================================
@rule
def r701_ranking_includes_efficiency(exp):
    """Vina scores rise with heavy-atom count, so ranking on the raw number
    alone partly ranks by molecular weight."""
    ranking = exp.ranking()
    if not ranking or "affinity" not in ranking:
        return []
    if "ligand_efficiency" in ranking:
        return [Finding(PASS, "R701", "ranked by affinity and by efficiency")]
    return [Finding(
        WARN, "R701", "ranking on raw affinity only",
        "the largest molecule in a set tends to win by default",
        "add 'by ligand_efficiency'",
        exp.rank.line if exp.rank else None)]


@rule
def r702_top_hit_is_heaviest(exp):
    """
    When the best-scoring compound is also the biggest, the size confound is
    live and a reviewer will say so. Recording it is not a criticism of the
    result; it is the answer prepared in advance.
    """
    affinities = exp.affinities()
    compounds = {c.name: c for c in exp.compounds() if c.heavy_atoms}
    scored = {n: a for n, a in affinities.items() if n in compounds}
    if len(scored) < 3:
        return []
    best = min(scored, key=lambda n: scored[n])
    heaviest = max(scored, key=lambda n: compounds[n].heavy_atoms)
    efficiencies = exp.ligand_efficiencies()
    if best != heaviest or not efficiencies:
        return []
    most_efficient = max(efficiencies, key=lambda n: efficiencies[n])
    return [Finding(
        RECORD, "R702", "the strongest binder is also the heaviest compound",
        "%s scores best (%.3f kcal/mol) with %d heavy atoms. Per atom it "
        "ranks %.3f, and %s leads at %.3f."
        % (best.replace("_", " "), scored[best], compounds[best].heavy_atoms,
           efficiencies.get(best, 0.0), most_efficient.replace("_", " "),
           efficiencies[most_efficient]),
        "report both columns so the question is answered before it is asked")]


# ==========================================================================
# R8xx  what the study is allowed to claim
# ==========================================================================
@rule
def r801_claim_scope(exp):
    """
    The compounds here are documented across the genus, in Z. jujuba and
    Z. mauritiana. They were not isolated from Z. talanai material. A
    species-level claim on genus-level evidence is the specific over-claim
    this project exists to avoid.
    """
    if exp.study is None:
        return []
    scope = (exp.study.one("claim_scope") or "").strip()
    title = exp.study.arg or ""
    if not scope:
        return [Finding(
            REFUSE, "R801", "claim scope is not declared",
            "a screen has to say what it is a claim about",
            "claim_scope genus Ziziphus", exp.study.line)]
    if chem.looks_like_species(scope):
        source = (exp.ligands.one("source") if exp.ligands else "") or "the library"
        return [Finding(
            REFUSE, "R801", "claim is species-level, evidence is genus-level",
            "claim_scope is '%s', but the compounds come from %s, reported "
            "across the genus rather than isolated from that species."
            % (scope, source),
            "claim_scope genus Ziziphus", exp.study.line_of("claim_scope"))]
    if re.search(r"talanai", title, re.I) and "genus" not in scope.lower():
        return [Finding(
            WARN, "R801", "the title names a species, the scope does not",
            "title: %s / scope: %s" % (title, scope),
            "make the title say candidate compounds from genus Ziziphus")]
    return [Finding(PASS, "R802", "claim scope declared", scope)]


@rule(scope="library")
def r804_library_is_genus_level(exp):
    """
    The ziziphus/ namespace holds shared, genus-level material. A library that
    declares a species scope would let every experiment importing it inherit a
    species-level claim from genus-level evidence, quietly and at scale.

    This is why the shared namespace is named for the genus rather than the
    species: the boundary is in the name, and this rule keeps it there.
    """
    out = []
    for block in exp.document.find_all("library"):
        scope = (block.one("scope") or "").strip()
        if not scope:
            out.append(Finding(
                REFUSE, "R804", "library '%s' declares no scope" % block.arg,
                "a shared compound set has to say what it is evidence about",
                "scope genus Ziziphus", block.line))
        elif chem.looks_like_species(scope):
            out.append(Finding(
                REFUSE, "R804", "library '%s' declares a species scope" % block.arg,
                "scope is '%s'. Everything importing this library would "
                "inherit a species-level claim." % scope,
                "the ziziphus/ namespace is genus-level. Use 'genus Ziziphus' "
                "and record species in 'reported_in' instead.",
                block.line_of("scope")))
        else:
            out.append(Finding(PASS, "R804",
                               "library '%s' is genus-level" % block.arg, scope))
    return out


@rule
def r803_in_silico_language(exp):
    """
    A docking score is a prediction. Verbs that assert demonstrated activity
    turn a prediction into a finding without any new evidence, which is the
    in-silico over-claim reviewers look for first.
    """
    if exp.study is None:
        return []
    text = " ".join([exp.study.arg or ""] + exp.study.all("claim"))
    banned = re.findall(
        r"\b(inhibits|cures|treats|proves|demonstrates|confirmed|"
        r"more potent than|outperforms)\b", text, re.I)
    if not banned:
        return []
    return [Finding(
        WARN, "R803", "the claim reads as demonstrated, not predicted",
        "found: %s" % ", ".join(sorted(set(w.lower() for w in banned))),
        "docking predicts binding. Say 'is predicted to bind more tightly "
        "than' unless an assay is attached.", exp.study.line)]


# ==========================================================================
# R9xx  the receptor's provenance
# ==========================================================================
@rule
def r901_surrogate_disclosed(exp):
    """
    3A4A is Saccharomyces cerevisiae isomaltase standing in for the human
    intestinal enzyme. That substitution is standard and defensible, and it is
    the first thing a reviewer asks about. Burying it is what causes trouble.
    """
    if exp.receptor is None or not exp.receptor.has("organism"):
        return []
    if exp.receptor.has("surrogate_for"):
        return [Finding(
            PASS, "R901", "surrogate receptor disclosed",
            "%s standing in for %s" % (exp.receptor.one("organism"),
                                       exp.receptor.one("surrogate_for")))]
    return [Finding(
        WARN, "R901", "receptor organism differs from the target, undisclosed",
        "organism is %s" % exp.receptor.one("organism"),
        "add 'surrogate_for' so the methods state it rather than omit it",
        exp.receptor.line)]


@rule
def r902_resolution(exp):
    """Resolution bounds how much can be said about side-chain placement, and
    therefore about the contacts a pose reports."""
    if exp.receptor is None:
        return []
    resolution = chem.as_float(exp.receptor.one("resolution"))
    if resolution is None:
        return [Finding(
            WARN, "R902", "receptor resolution not recorded",
            "it bounds what the contact analysis can claim",
            "add 'resolution 1.6 A'", exp.receptor.line)]
    if resolution > 2.5:
        return [Finding(
            WARN, "R902", "receptor resolution is low for contact analysis",
            "%.2f A. Side-chain positions carry real uncertainty here." % resolution,
            "keep hydrogen-bond claims cautious, or pick a better structure")]
    return [Finding(PASS, "R902", "receptor resolution recorded",
                    "%.2f A" % resolution)]


# ==========================================================================
# RAxx  replication and reporting
# ==========================================================================
@rule
def ra01_replication_kind(exp):
    """
    Re-running one seed through a rebuilt pipeline proves the pipeline is
    faithful. It does not test whether the search converged. Both are worth
    having and neither substitutes for the other, so the record has to say
    which one was done.
    """
    if exp.replication is None:
        return []
    kind = (exp.replication.one("kind") or "").strip().lower()
    if not kind:
        return [Finding(
            WARN, "RA01", "replication kind not stated",
            "same seed and different seed answer different questions",
            "kind same_seed  or  kind different_seed",
            exp.replication.line)]
    if kind == "same_seed":
        return [Finding(
            RECORD, "RA01", "replication is same-seed",
            "this demonstrates the pipeline reproduces the published numbers, "
            "which is a real result.",
            "it does not measure search convergence. Report the two "
            "separately rather than as one claim.")]
    return [Finding(RECORD, "RA01", "replication is %s" % kind)]


@rule
def ra02_replication_deviation(exp):
    """
    This checks run-to-run numerical reproducibility, not scoring precision.
    AutoDock Vina publishes no scoring-precision figure; its own accuracy
    against experimental binding affinities is reported at 2 to 3 kcal/mol.
    What this rule measures is a narrower, different question: whether
    repeating the same pipeline on the same inputs reproduces its own earlier
    number. Agreement inside that project convention is the point of the
    exercise; disagreement beyond it means something in the pipeline changed.
    """
    original, rerun = exp.affinities(), exp.reruns()
    shared = [n for n in rerun if n in original]
    if not shared:
        return []
    worst_name = max(shared, key=lambda n: abs(rerun[n] - original[n]))
    worst = abs(rerun[worst_name] - original[worst_name])
    if worst > REPLICATION_REPRODUCIBILITY_LIMIT_KCAL_MOL:
        return [Finding(
            WARN, "RA02", "replication deviates beyond run-to-run reproducibility",
            "%s differs by %.3f kcal/mol, past the %.1f kcal/mol convention "
            "used here for pipeline reproducibility. That convention is not "
            "a documented Vina precision figure."
            % (worst_name.replace("_", " "), worst,
               REPLICATION_REPRODUCIBILITY_LIMIT_KCAL_MOL),
            "something in the pipeline changed. Find it before publishing both.")]
    return [Finding(
        PASS, "RA02", "replication agrees within run-to-run reproducibility",
        "largest deviation %.3f kcal/mol (%s), across %d compounds, under "
        "the %.1f kcal/mol convention used here. Vina itself documents no "
        "scoring-precision figure; its accuracy against experiment is "
        "reported at 2 to 3 kcal/mol."
        % (worst, worst_name.replace("_", " "), len(shared),
           REPLICATION_REPRODUCIBILITY_LIMIT_KCAL_MOL))]
