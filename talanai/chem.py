"""
The small amount of chemistry Talanai needs to do arithmetic on.

Nothing here models a molecule. It counts atoms in a written formula and does
one division, because that division is the difference between "Rutin is the
best binder" and "Rutin is the heaviest molecule in the set".
"""

from __future__ import annotations

import re

# compounds.js writes formulae with unicode subscripts (C27H30O16 as C₂₇H₃₀O₁₆).
SUBSCRIPTS = str.maketrans(
    "₀₁₂₃₄₅₆₇₈₉",
    "0123456789",
)

_TOKEN = re.compile(r"([A-Z][a-z]?)(\d*)")


def normalise_formula(formula):
    """C₂₇H₃₀O₁₆ -> C27H30O16. Idempotent."""
    return (formula or "").translate(SUBSCRIPTS).strip()


def heavy_atoms(formula):
    """Non-hydrogen atom count. Returns 0 for anything unparseable."""
    formula = normalise_formula(formula)
    if not formula or not re.fullmatch(r"(?:[A-Z][a-z]?\d*)+", formula):
        return 0
    total = 0
    for element, count in _TOKEN.findall(formula):
        if not element or element == "H":
            continue
        total += int(count) if count else 1
    return total


def ligand_efficiency(affinity, formula):
    """
    Binding score per heavy atom, the standard size correction.

    AutoDock Vina scores scale with molecular size, so a large molecule tends
    to score well simply by having more atoms in contact. Ligand efficiency
    divides that back out. Vina's score is not a true free energy, so treat
    this as a size correction rather than a thermodynamic quantity.

    Returned positive (0.35 reads better than -0.35 in a table).
    """
    n = heavy_atoms(formula)
    if not n or affinity is None:
        return None
    return abs(float(affinity)) / n


def as_float(text):
    """First number in a string, or None. '-8.857 kcal/mol' -> -8.857"""
    if text is None:
        return None
    match = re.search(r"-?\d+(?:\.\d+)?", str(text))
    return float(match.group(0)) if match else None


def numbers(text):
    """Every number in a string. '30 30 30' -> [30.0, 30.0, 30.0]"""
    return [float(x) for x in re.findall(r"-?\d+(?:\.\d+)?", text or "")]


def looks_like_species(text):
    """
    True when a string names a species rather than a genus.

    'Ziziphus talanai' and 'Z. talanai' are species. 'Ziziphus' and
    'genus Ziziphus' are not. Used to keep the ziziphus/ library honest.
    """
    if not text:
        return False
    cleaned = re.sub(r"\b(genus|from|the)\b", " ", text, flags=re.I).strip()
    return bool(re.search(r"\b([A-Z][a-z]{2,}|[A-Z]\.)\s+[a-z]{3,}\b", cleaned))
