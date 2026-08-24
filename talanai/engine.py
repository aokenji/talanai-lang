"""
The runner. Offline, local AutoDock Vina, no network.

Design constraints this file exists under:
  - Python standard library only. No pip, no bindings; Vina is invoked as the
    executable it already is.
  - Fully offline. Nothing here contacts the Hugging Face Space or any host.
  - i3 with 8 GB. One ligand at a time by default, and a content-addressed
    cache so an interrupted screen resumes instead of restarting.
  - `tal check` never reaches this module, so validation stays instant and
    needs no docking engine installed.

The cache key is the whole point of the design: it is a hash of everything
that can change the number, which is receptor file and preparation, box,
ligand file and preparation, engine version, exhaustiveness, modes, energy
range and seed. Change any of them and the old result is not reused, because
it is not the same measurement.
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess

CACHE_DIRNAME = ".talanai-cache"

# Where a local Vina tends to live on this machine. TALANAI_VINA overrides all
# of it, and a `vina_path` in the dock block overrides that.
CANDIDATE_PATHS = (
    # The installed TalanaiDock, which is the intended host for this tool.
    r"D:\THESIS_VSC\BALAKATDBV2\talanaidock\dist\TalanaiDock\app\docking_assets\vina.exe",
    r"D:\Local Disk Downloads\TalanaiDock\app\docking_assets\vina.exe",
    r"D:\THESIS_VSC\BALAKATDBV2\talanaidock\odysseus-src\docking_assets\vina.exe",
    r"C:\Program Files\AutoDock Vina\bin\vina.exe",
    "/usr/local/bin/vina",
    "/usr/bin/vina",
)


class EngineMissing(Exception):
    """Raised when no local Vina can be found. Never fall back to a guess."""


class Job:
    """One ligand, one seed, one docking run."""

    def __init__(self, ligand, ligand_file, seed, receptor_file, box, settings):
        self.ligand = ligand
        self.ligand_file = ligand_file
        self.seed = seed
        self.receptor_file = receptor_file
        self.box = box                  # (center, size)
        self.settings = settings        # dict of engine settings

    def fingerprint(self):
        """Everything that can change the number, hashed."""
        center, size = self.box
        payload = {
            "ligand": self.ligand,
            "ligand_file": _file_identity(self.ligand_file),
            "receptor_file": _file_identity(self.receptor_file),
            "center": center,
            "size": size,
            "seed": self.seed,
            "settings": self.settings,
        }
        blob = json.dumps(payload, sort_keys=True, separators=(",", ":"))
        return hashlib.sha256(blob.encode("utf-8")).hexdigest()[:16]

    def config_text(self):
        center, size = self.box
        lines = [
            "receptor = %s" % self.receptor_file,
            "ligand = %s" % self.ligand_file,
            "center_x = %s" % center[0],
            "center_y = %s" % center[1],
            "center_z = %s" % center[2],
            "size_x = %s" % size[0],
            "size_y = %s" % size[1],
            "size_z = %s" % size[2],
            "seed = %s" % self.seed,
        ]
        for key, flag in (("exhaustiveness", "exhaustiveness"),
                          ("modes", "num_modes"),
                          ("energy_range", "energy_range"),
                          ("cpu", "cpu")):
            value = self.settings.get(key)
            if value is not None:
                lines.append("%s = %s" % (flag, int(float(value))))
        return "\n".join(lines) + "\n"

    def command(self, vina, config_path, out_path):
        return [vina, "--config", config_path, "--out", out_path]

    def label(self):
        return "%s seed %s" % (self.ligand, self.seed)


def _file_identity(path):
    """Path plus size plus mtime. Cheap, and catches a swapped input file."""
    if not path or not os.path.isfile(path):
        return {"path": path, "present": False}
    stat = os.stat(path)
    return {"path": os.path.basename(path), "present": True,
            "size": stat.st_size, "mtime": int(stat.st_mtime)}


def find_vina(experiment=None):
    """Locate a local Vina. Returns a path, or raises EngineMissing."""
    declared = None
    if experiment is not None and experiment.dock is not None:
        declared = experiment.dock.one("vina_path")
    for candidate in [declared, os.environ.get("TALANAI_VINA")]:
        if candidate and os.path.isfile(candidate):
            return candidate
    on_path = shutil.which("vina") or shutil.which("vina.exe")
    if on_path:
        return on_path
    for candidate in CANDIDATE_PATHS:
        if os.path.isfile(candidate):
            return candidate
    raise EngineMissing(
        "No local AutoDock Vina found. Talanai runs offline and will not "
        "substitute a remote service.\n"
        "  Set TALANAI_VINA to the executable, or add 'vina_path' to the "
        "dock block, or put vina on PATH.")


def plan(experiment):
    """Build the full job list for an experiment. Does not run anything."""
    base = os.path.dirname(os.path.abspath(experiment.path))
    receptor_file = experiment.receptor.one("file") if experiment.receptor else None
    if receptor_file and not os.path.isabs(receptor_file):
        receptor_file = os.path.join(base, receptor_file)

    ligand_dir = experiment.ligands.one("dir") if experiment.ligands else None
    if ligand_dir and not os.path.isabs(ligand_dir):
        ligand_dir = os.path.join(base, ligand_dir)

    center, size = experiment.box()
    settings = {}
    if experiment.dock:
        for key in ("exhaustiveness", "modes", "energy_range", "cpu"):
            value = experiment.dock.one(key)
            if value is not None:
                settings[key] = value
        settings["engine"] = experiment.dock.one("engine")

    seeds = experiment.seeds() or [42]
    jobs = []
    for compound in experiment.compounds():
        ligand_file = None
        if compound.filename and ligand_dir:
            ligand_file = os.path.join(ligand_dir, compound.filename)
        elif compound.filename:
            ligand_file = compound.filename
        elif ligand_dir:
            stems = (compound.name, compound.name.lower(),
                     compound.name.split("_")[0].lower())
            for stem in stems:
                for extension in (".pdbqt", ".sdf", ".mol2"):
                    candidate = os.path.join(ligand_dir, stem + extension)
                    if os.path.isfile(candidate):
                        ligand_file = candidate
                        break
                if ligand_file:
                    break
            if ligand_file is None:
                ligand_file = os.path.join(ligand_dir, compound.name + ".pdbqt")
        for seed in seeds:
            jobs.append(Job(compound.name, ligand_file, seed, receptor_file,
                            (center, size), settings))
    return jobs


def missing_inputs(jobs):
    """Every input file a plan needs that is not on disk."""
    missing = []
    for job in jobs:
        for path in (job.receptor_file, job.ligand_file):
            if not path:
                missing.append("(not declared)")
            elif not os.path.isfile(path) and path not in missing:
                missing.append(path)
    return sorted(set(missing))


def cache_dir(experiment):
    return os.path.join(os.path.dirname(os.path.abspath(experiment.path)),
                        CACHE_DIRNAME)


def cached(experiment, job):
    path = os.path.join(cache_dir(experiment), job.fingerprint() + ".json")
    if not os.path.isfile(path):
        return None
    try:
        with open(path, encoding="utf-8") as handle:
            return json.load(handle)
    except (OSError, ValueError):
        return None


def store(experiment, job, record):
    directory = cache_dir(experiment)
    os.makedirs(directory, exist_ok=True)
    path = os.path.join(directory, job.fingerprint() + ".json")
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(record, handle, indent=2, sort_keys=True)
    return path


def parse_affinity(stdout):
    """Best-mode affinity from Vina's result table."""
    for line in stdout.splitlines():
        parts = line.split()
        if len(parts) >= 2 and parts[0] == "1":
            try:
                return float(parts[1])
            except ValueError:
                continue
    return None


def execute(experiment, job, vina, workdir, timeout=1800):
    """Run one job. Returns a record dict. Raises on engine failure."""
    os.makedirs(workdir, exist_ok=True)
    stem = "%s_seed%s_%s" % (job.ligand, job.seed, job.fingerprint())
    config_path = os.path.join(workdir, stem + ".conf")
    out_path = os.path.join(workdir, stem + "_out.pdbqt")
    with open(config_path, "w", encoding="utf-8") as handle:
        handle.write(job.config_text())

    completed = subprocess.run(
        job.command(vina, config_path, out_path),
        capture_output=True, text=True, timeout=timeout)

    record = {
        "ligand": job.ligand,
        "seed": job.seed,
        "fingerprint": job.fingerprint(),
        "affinity": parse_affinity(completed.stdout),
        "returncode": completed.returncode,
        "pose_file": out_path if os.path.isfile(out_path) else None,
        "engine": job.settings.get("engine"),
    }
    if completed.returncode != 0 or record["affinity"] is None:
        record["stderr"] = completed.stderr[-2000:]
    return record
