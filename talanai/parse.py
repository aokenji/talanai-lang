"""
Reading .tal files.

This is a reader, not a parser generator. There is no grammar file, no lexer,
no AST and no dependency. The format is line-oriented on purpose so that the
whole thing fits in one screen and a biologist can predict what it will do:

    #                       starts a comment, to end of line
    receptor 3A4A           a line at column 0 opens a block, with an argument
      resolution  1.6 A     an indented line is a key and its value
      compound    Rutin     a key may repeat; repeats accumulate in order
      note        text      a value may continue on the next indented line if
                  more      that line's first word is not a known key (see
                            Block.fold_continuations)

Everything downstream works on Block objects. Nothing else in Talanai touches
raw text.
"""

from __future__ import annotations


class Block:
    """One block of an experiment file, such as `receptor 3A4A`."""

    def __init__(self, name, arg, line):
        self.name = name
        self.arg = arg
        self.line = line
        self.keys = {}      # key -> list of values, in file order
        self.lines = {}     # key -> line number of first occurrence

    # -- reading -----------------------------------------------------------
    def all(self, key):
        return self.keys.get(key, [])

    def one(self, key, default=None):
        values = self.keys.get(key)
        return values[0] if values else default

    def has(self, key):
        return key in self.keys

    def line_of(self, key):
        return self.lines.get(key, self.line)

    def __repr__(self):
        return "<Block %s %r at line %d>" % (self.name, self.arg, self.line)


class Document:
    """A parsed .tal file."""

    def __init__(self, path, blocks, problems):
        self.path = path
        self.blocks = blocks
        self.problems = problems    # list of (line, message) from the reader

    def find(self, name):
        for block in self.blocks:
            if block.name == name:
                return block
        return None

    def find_all(self, name):
        return [b for b in self.blocks if b.name == name]


def _split_key_value(text):
    parts = text.split(None, 1)
    key = parts[0]
    value = parts[1].strip() if len(parts) > 1 else ""
    return key, value


def parse_text(text, path="<string>", known_keys=None):
    """Turn the contents of a .tal file into a Document.

    known_keys maps a block name to the set of keys it recognises. It is used
    only to decide whether an indented line continues the previous value or
    starts a new key. Validation of key names happens later, in rules.py.
    """
    # Notepad and PowerShell write UTF-8 with a byte-order mark. Left in, it
    # becomes an invisible first character and the first block stops parsing.
    if text.startswith("﻿"):
        text = text[1:]

    blocks, problems = [], []
    current = None
    last_key = None
    last_indent = 0

    for number, raw in enumerate(text.splitlines(), start=1):
        line = raw.split("#", 1)[0].rstrip()
        if not line.strip():
            last_key = None
            continue

        if not line[0].isspace():
            # Column 0: a new block.
            name, arg = _split_key_value(line.strip())
            current = Block(name, arg, number)
            blocks.append(current)
            last_key = None
            continue

        # Indented: a key, or the continuation of the previous value.
        if current is None:
            problems.append((number, "indented line before any block: %s"
                             % line.strip()))
            continue

        key, value = _split_key_value(line.strip())
        indent = len(line) - len(line.lstrip())

        # A continuation is indented DEEPER than the key it continues. Using
        # indentation rather than "is this word a known key" matters: the
        # earlier heuristic quietly folded a misspelt key into the previous
        # value, which is precisely the mistake rule R003 exists to catch.
        is_continuation = last_key is not None and indent > last_indent
        if is_continuation:
            existing = current.keys[last_key]
            existing[-1] = (existing[-1] + " " + line.strip()).strip()
            continue

        current.keys.setdefault(key, []).append(value)
        current.lines.setdefault(key, number)
        last_key = key
        last_indent = indent

    return Document(path, blocks, problems)


def parse_file(path, known_keys=None):
    with open(path, encoding="utf-8") as handle:
        return parse_text(handle.read(), path=path, known_keys=known_keys)
