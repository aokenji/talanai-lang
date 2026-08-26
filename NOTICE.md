# Licence scope

The `LICENSE` file is the plain MIT text, deliberately unmodified. GitHub's
licence detector, and Zenodo downstream of it, only recognise a licence when the
file matches the canonical wording, so the scope note lives here instead of
being appended to it.

## What MIT covers

The **software**: the Talanai language, its validator, the `talanai/` package,
the tests, and the `run_*.py` scripts.

## What it does not cover

The **scientific content** those scripts produced. TalanaiHub's data and prose
are CC BY 4.0, and the run records under `validation-run/` and
`validation-inputs/` are research data belonging to the *Ziziphus*
alpha-glucosidase study. They are distributed here because a rule set without
the runs it was built against cannot be checked by anyone, but they carry that
study's CC BY 4.0 terms, not MIT.

## Why the split

CC BY 4.0 is right for data and prose and wrong for software: it is not an
OSI-approved open source licence. The Journal of Open Source Software, the
realistic venue for a tools paper on this after the December defence, requires
one. Licensing the code MIT and leaving the data CC BY 4.0 satisfies both
without relicensing anything that is not ours to relicense.

## Citing

See `CITATION.cff`. The dataset has its own DOI, `10.5281/zenodo.20384660`. A
software DOI for this repository is not minted yet; per the 2026-08-20
publication decision it is cut from a release tag, and a Zenodo record is
permanent.
