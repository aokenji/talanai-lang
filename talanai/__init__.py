"""
Talanai: a language for describing molecular docking experiments.

A .tal file describes one experiment completely, and the validator refuses it
when the experiment could not mean what it claims to mean. The validator is
the product; docking is the easy part.

    from talanai import load, validate

Python standard library only, by design: this has to run on the interpreter
bundled with TalanaiDock, offline, on a laptop.
"""

from __future__ import annotations

__version__ = "0.1.0"

from . import chem, engine, model, parse, pdb, report, rules   # noqa: F401
from .model import Experiment                                   # noqa: F401
from .rules import REFUSE, WARN, RECORD, PASS                    # noqa: F401


def load(path):
    """Read a .tal file into an Experiment."""
    document = parse.parse_file(path, known_keys=model.KNOWN_KEYS)
    return model.Experiment(document)


def validate(path_or_experiment):
    """Return (experiment, findings, locked)."""
    experiment = (path_or_experiment if isinstance(path_or_experiment, Experiment)
                  else load(path_or_experiment))
    findings = rules.run_all(experiment)
    return experiment, findings, rules.docking_is_locked(findings)
