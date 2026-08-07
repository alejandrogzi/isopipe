#!/usr/bin/env python3
"""Diff published pipeline outputs against committed golden truth.

Usage:
    python3 compare.py GOLD_DIR RESULTS_DIR

Normalizes (gunzip, sort lines, strip trailing whitespace) every text file
under GOLD_DIR, compares it to the same relative path under RESULTS_DIR, and
exits 1 if anything differs or if RESULTS_DIR holds text outputs missing from
gold (i.e. new/renamed outputs).
"""
import gzip
import sys
from pathlib import Path

BINARY_SUFFIXES = (".bam", ".bai", ".pbi", ".bb", ".bw", ".html", ".pdf", ".png", ".mmi", ".sam")
EXCLUDED_DIRS = {"PIPELINE_INFO", "01_PBCCS"}  # ponytail: PBCCS QC embeds runtime ms; not reproducible


def head_diff(path_a, path_b, max_lines=12):
    a, b = set(normalized(path_a)), set(normalized(path_b))
    out = []
    for line in sorted(a - b)[:max_lines]:
        out.append(f"  - {line}")
    for line in sorted(b - a)[:max_lines]:
        out.append(f"  + {line}")
    return out


def normalized(path: Path):
    data = path.read_bytes()
    if path.name.endswith(".gz"):
        data = gzip.decompress(data)
    lines = [l.rstrip() for l in data.decode("utf-8", "replace").splitlines() if l.strip()]
    return sorted(lines)


def collect(root: Path):
    return {
        p.relative_to(root)
        for p in root.rglob("*")
        if p.is_file()
        and not any(part in EXCLUDED_DIRS for part in p.parts)
        and not p.name.startswith(".")
        and not p.name.endswith(BINARY_SUFFIXES)
    }


def main() -> int:
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    if len(args) != 2:
        print(__doc__)
        return 2
    show_diffs = "--show-diffs" in sys.argv
    gold_dir, results_dir = Path(args[0]), Path(args[1])
    if not gold_dir.is_dir() or not results_dir.is_dir():
        print(f"error: {gold_dir} and {results_dir} must exist "
              "(bootstrap: run pipeline on a known-good commit, then copy published outputs into gold/)")
        return 1

    problems = []
    for rel in sorted(collect(gold_dir)):
        got = results_dir / rel
        if not got.is_file():
            problems.append(f"missing: {rel}")
            continue
        if normalized(gold_dir / rel) != normalized(got):
            problems.append(f"changed: {rel}")
            if show_diffs:
                problems.extend("    " + d for d in head_diff(gold_dir / rel, got))
    for rel in sorted(collect(results_dir) - collect(gold_dir)):
        problems.append(f"new output not in gold: {rel}")

    if problems:
        print("\n".join(problems))
        return 1
    print(f"OK: {len(collect(gold_dir))} outputs match gold")
    return 0


if __name__ == "__main__":
    sys.exit(main())
