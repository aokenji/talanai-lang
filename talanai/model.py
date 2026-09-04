"""
The vocabulary of a .tal file, and typed access to a parsed one.

THE ELEVEN WORDS
    An experiment file is built from eleven block words and nothing else:

        study        what this is, who ran it, what it claims
        receptor     the protein, and how it was prepared
        site         where in the protein we are looking
        control      the redocking control that unlocks docking
        ligands      what is being screened, and how it was prepared
        reference    the comparator drug
        dock         the search settings
        rank         how results are ordered
        report       what to emit
        results      what was obtained
        replication  a re-run, and what kind it was

    Library files in ziziphus/ use one further word, `library`, and are a
    different kind of file: a shared, genus-level compound set.
"""

from __future__ import annotations

from . import chem

# --------------------------------------------------------------------------
# Schema. Keys are listed so the reader can tell a new key from a wrapped
# line, and so rules.py can warn about typos. Unknown keys are a warning,
# never a silent ignore.
# --------------------------------------------------------------------------
SCHEMA = {
    "study": {"author", "date", "doi", "claim_scope", "claim", "note"},
    "receptor": {"name", "organism", "surrogate_for", "resolution", "prepare",
                 "prepare_note", "file", "note", "checksum"},
    "site": {"center", "size", "anchored_on", "anchor_note", "must_enclose",
             "blind", "note"},
    # A control records its OWN configuration, because a control run under
    # different settings validates a different protocol.
    "control": {"redock", "measure", "require", "result", "source",
                "result_note", "note", "prepare", "box_size",
                "exhaustiveness", "receptor_file", "ligand_file",
                "reference_file", "seed", "checksum"},
    "ligands": {"source", "dir", "prepare", "protonate", "compound", "note",
                "checksum"},
    "reference": {"formula", "prepare", "protonate", "role", "note", "file",
                  "checksum"},
    "dock": {"engine", "exhaustiveness", "modes", "energy_range", "seeds",
             "cpu", "vina_path", "note"},
    "rank": {"by", "note"},
    "report": set(),
    "results": {"affinity", "pose_cluster_rmsd", "pose_file", "note"},
    "replication": {"date", "pipeline", "kind", "rerun", "note"},
    # library files
    "library": {"scope", "reported_in", "compound", "note", "source"},
}

EXPERIMENT_BLOCKS = ("study", "receptor", "site", "control", "ligands",
                     "reference", "dock", "rank", "report", "results",
                     "replication")

REQUIRED_BLOCKS = ("study", "receptor", "site", "control", "ligands", "dock")

SINGLE_VALUE_KEYS = {
    "center", "size", "resolution", "prepare", "exhaustiveness", "modes",
    "energy_range", "seeds", "result", "require", "measure", "redock",
    "claim_scope", "organism", "surrogate_for", "kind", "scope", "engine",
    "pose_file",
}

KNOWN_KEYS = {name: keys for name, keys in SCHEMA.items()}

ACCEPTED_REDOCK_LIMIT = 2.0     # angstrom, the field's conventional bar


class Compound:
    def __init__(self, name, formula=None, filename=None, reported_in=None):
        self.name = name
        self.formula = chem.normalise_formula(formula) if formula else None
        self.filename = filename        # optional prepared structure
        self.reported_in = reported_in

    @property
    def heavy_atoms(self):
        return chem.heavy_atoms(self.formula) if self.formula else 0

    def display_name(self):
        return self.name.replace("_", " ")

    def __repr__(self):
        return "<Compound %s %s>" % (self.name, self.formula)


class Experiment:
    """Typed view over a parsed .tal document."""

    def __init__(self, document):
        self.document = document
        self.path = document.path

    # -- blocks ------------------------------------------------------------
    def block(self, name):
        return self.document.find(name)

    @property
    def study(self):
        return self.block("study")

    @property
    def receptor(self):
        return self.block("receptor")

    @property
    def site(self):
        return self.block("site")

    @property
    def control(self):
        return self.block("control")

    @property
    def ligands(self):
        return self.block("ligands")

    @property
    def reference(self):
        return self.block("reference")

    @property
    def dock(self):
        return self.block("dock")

    @property
    def rank(self):
        return self.block("rank")

    @property
    def results(self):
        return self.block("results")

    @property
    def replication(self):
        return self.block("replication")

    # -- derived -----------------------------------------------------------
    @property
    def kind(self):
        """'library' for a ziziphus/ compound set, otherwise 'experiment'.

        A library is a different kind of file and must not be judged against
        the rules for running an experiment.
        """
        if self.block("library") is not None and self.study is None:
            return "library"
        return "experiment"

    @property
    def title(self):
        if self.study:
            return self.study.arg.strip('"')
        library = self.block("library")
        if library:
            return "library %s" % library.arg
        return "untitled"

    def compounds(self):
        """Ligands plus the reference, in file order."""
        out = []
        if self.ligands:
            for entry in self.ligands.all("compound"):
                parts = entry.split()
                if parts:
                    out.append(Compound(
                        parts[0],
                        parts[1] if len(parts) > 1 else None,
                        parts[2] if len(parts) > 2 else None))
        if self.reference and self.reference.arg:
            out.append(Compound(self.reference.arg,
                                self.reference.one("formula"),
                                self.reference.one("file")))
        return out

    def compound_by_name(self, name):
        wanted = name.replace(" ", "_").lower()
        for compound in self.compounds():
            if compound.name.lower() == wanted:
                return compound
        return None

    def affinities(self):
        """name -> float, from the results block."""
        out = {}
        if not self.results:
            return out
        for entry in self.results.all("affinity"):
            parts = entry.split()
            if len(parts) >= 2:
                value = chem.as_float(parts[1])
                if value is not None:
                    out[parts[0]] = value
        return out

    def reruns(self):
        out = {}
        if not self.replication:
            return out
        for entry in self.replication.all("rerun"):
            parts = entry.split()
            if len(parts) >= 2:
                value = chem.as_float(parts[1])
                if value is not None:
                    out[parts[0]] = value
        return out

    def box(self):
        """(center, size) as two lists of three floats, or (None, None)."""
        if not self.site:
            return None, None
        center = chem.numbers(self.site.one("center", ""))
        size = chem.numbers(self.site.one("size", ""))
        return (center if len(center) == 3 else None,
                size if len(size) == 3 else None)

    def box_volume(self):
        _, size = self.box()
        if not size:
            return None
        return size[0] * size[1] * size[2]

    def seeds(self):
        if not self.dock:
            return []
        return [int(s) for s in chem.numbers(self.dock.one("seeds", ""))]

    def exhaustiveness(self):
        return chem.as_float(self.dock.one("exhaustiveness")) if self.dock else None

    def control_result(self):
        """The recorded redock result as a float, or None when not recorded."""
        if not self.control:
            return None
        raw = (self.control.one("result") or "").strip().lower()
        if raw in ("", "n.d.", "nd", "none", "null", "pending", "unknown"):
            return None
        return chem.as_float(raw)

    def control_limit(self):
        if not self.control:
            return None
        return chem.as_float(self.control.one("require", ""))

    def ranking(self):
        return [v.strip().lower() for v in (self.rank.all("by") if self.rank else [])]

    def ligand_efficiencies(self):
        """name -> score per heavy atom, for everything with a formula."""
        out = {}
        affinities = self.affinities()
        for compound in self.compounds():
            score = affinities.get(compound.name)
            if score is None or not compound.formula:
                continue
            value = chem.ligand_efficiency(score, compound.formula)
            if value is not None:
                out[compound.name] = value
        return out
