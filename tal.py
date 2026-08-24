#!/usr/bin/env python3
"""Entry point. `python tal.py check file.tal`, or use the tal / tal.cmd wrapper."""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from talanai.cli import main  # noqa: E402

if __name__ == "__main__":
    sys.exit(main())
