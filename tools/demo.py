#!/usr/bin/env python3
"""
Exhibit demo. Four beats, about ninety seconds, no network and no engine.

    python tools/demo.py

The point of the demo is the LANGUAGE, not the verdict. Beat 1 puts a block of
.tal on the screen and lets people read it, because the claim being made is that
a docking protocol can be written so a biologist reads it and a machine checks
it, and the only way to support that claim is to let someone read one.

Beats 2 to 4 then show the checking half: the real file accepts, a one-field
change makes it refuse, and the real file is untouched throughout.

WHAT BEAT 3 IS AND IS NOT
-------------------------
It is a COUNTERFACTUAL, not a historical artifact. It answers "what happens if
the control validated a different receptor preparation than the screen", which
is a property of the configuration and is demonstrable live. It does NOT claim
to be a file that was ever on disk. Say that out loud when presenting: the
difference between a counterfactual and a claim about history is the whole
reason this demo does not ship a reconstructed "before" file.

The real example is never written to. The modified copy goes to a temp file and
is deleted on the way out.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
TAL = os.path.join(ROOT, "tal.py")
EXAMPLE = os.path.join(ROOT, "examples", "alpha-glucosidase.tal")

# The one field the demo changes, and what it becomes. Both halves are asserted
# against the file before anything runs, so a drifted example fails loudly here
# rather than silently demoing nothing.
CONTROL_PREPARE = "  prepare        meeko, polar H and Gasteiger charges"
COUNTERFACTUAL = "  prepare        raw cleaned 3A4A, no hydrogen-bond atom typing"

RULE = "  " + "-" * 70


def show(title, note=""):
    print("")
    print(RULE)
    print("  %s" % title)
    if note:
        print("  %s" % note)
    print(RULE)
    print("")


def wait(pause):
    if pause:
        try:
            input("\n  [enter] ")
        except (EOFError, KeyboardInterrupt):
            print("")
            raise SystemExit(0)


def read_block(path, first, last):
    """Lines first..last inclusive, 1-indexed, as they appear in the file."""
    with open(path, encoding="utf-8") as handle:
        lines = handle.read().split("\n")
    return "\n".join(lines[first - 1:last])


def run_check(path):
    """Run `tal check` and stream it. Exit code 1 means refused, not broken."""
    result = subprocess.run(
        [sys.executable, TAL, "check", path],
        capture_output=True, text=True,
        env=dict(os.environ, PYTHONIOENCODING="utf-8"),
    )
    sys.stdout.write(result.stdout)
    if result.stderr.strip():
        sys.stdout.write(result.stderr)
    return result.returncode


def main():
    parser = argparse.ArgumentParser(description="Talanai exhibit demo.")
    parser.add_argument("--no-pause", action="store_true",
                        help="run all four beats without waiting for a keypress")
    args = parser.parse_args()
    pause = not args.no_pause

    if not os.path.exists(EXAMPLE):
        print("  The example file is missing: %s" % EXAMPLE)
        return 2

    source = open(EXAMPLE, encoding="utf-8").read()
    if source.count(CONTROL_PREPARE) != 1:
        print("  The example has drifted: expected exactly one line reading")
        print("      %s" % CONTROL_PREPARE.strip())
        print("  Found %d. Fix this before demoing, do not present around it."
              % source.count(CONTROL_PREPARE))
        return 2

    # ---- beat 1: the language itself -------------------------------------
    show("1. This is a docking experiment, written in .tal",
         "examples/alpha-glucosidase.tal, lines 39 to 45")
    print(read_block(EXAMPLE, 39, 45))
    print("")
    print("  Read the last line aloud. The search box must enclose Asp215,")
    print("  Glu277 and Asp352, the catalytic residues. Rule R303 opens the")
    print("  receptor file on disk and checks that it does.")
    wait(pause)

    # ---- beat 2: it accepts ----------------------------------------------
    show("2. The validator reads it", "tal check examples/alpha-glucosidase.tal")
    code = run_check(EXAMPLE)
    print("  exit code %d" % code)
    if code != 0:
        print("")
        print("  Expected 0 (accepted). Something changed. Stop and look at it")
        print("  rather than continuing the demo.")
        return 1
    wait(pause)

    # ---- beat 3: change one field, it refuses -----------------------------
    show("3. Now change ONE field",
         "the control's receptor preparation, line 56")
    print("  from:  %s" % CONTROL_PREPARE.strip())
    print("    to:  %s" % COUNTERFACTUAL.strip())
    print("")
    print("  This is a counterfactual, not a file that ever existed. It asks")
    print("  what happens when the control validates a different preparation")
    print("  than the screen actually ran on.")
    wait(pause)

    handle, temp = tempfile.mkstemp(suffix=".tal", prefix="talanai-demo-")
    os.close(handle)
    try:
        with open(temp, "w", encoding="utf-8", newline="") as out:
            out.write(source.replace(CONTROL_PREPARE, COUNTERFACTUAL, 1))
        code = run_check(temp)
        print("  exit code %d" % code)
        if code != 1:
            print("")
            print("  Expected 1 (refused). R106 did not fire. Do not present")
            print("  this until you know why.")
            return 1
    finally:
        try:
            os.remove(temp)
        except OSError:
            pass

    wait(pause)

    # ---- beat 4: nothing was touched -------------------------------------
    show("4. The real file was never written to",
         "tal check examples/alpha-glucosidase.tal")
    code = run_check(EXAMPLE)
    print("  exit code %d" % code)
    print("")
    print("  The modified copy went to a temp file and has been deleted.")
    return 0 if code == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
