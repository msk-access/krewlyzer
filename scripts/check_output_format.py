#!/usr/bin/env python3
"""
Static analysis gate: ensure every output-writing call site passes output_format= and compress=.

Run before committing or as part of CI:
    python scripts/check_output_format.py
    # Exit 0 = all OK, Exit 1 = gaps found

Checks:
  - write_table() calls missing output_format=
  - process_region_entropy() calls missing output_format=
  - compute_and_write_gc_factors() calls missing output_format=
  - write_extraction_outputs() calls missing output_format=
  - _core.run_unified_pipeline() calls missing output_format
"""
from __future__ import annotations
import pathlib
import sys

ROOT = pathlib.Path(__file__).parent.parent

# Files to scan (relative to repo root)
SCAN_FILES = [
    "src/krewlyzer/core/unified_processor.py",
    "src/krewlyzer/core/sample_processor.py",
    "src/krewlyzer/core/fsc_processor.py",
    "src/krewlyzer/core/fsr_processor.py",
    "src/krewlyzer/core/region_entropy_processor.py",
    "src/krewlyzer/core/motif_processor.py",
    "src/krewlyzer/core/feature_serializer.py",
    "src/krewlyzer/core/gc_assets.py",
    "src/krewlyzer/wrapper.py",
    "src/krewlyzer/mfsd.py",
    "src/krewlyzer/uxm.py",
    "src/krewlyzer/region_mds.py",
]

# Each entry: (trigger substring in call line, required substring in call block)
# "call block" = the call line + next N lines (to handle multi-line calls)
RULES: list[tuple[str, str, int, list[str]]] = [
    # (trigger, required_param, lookahead_lines, allowed_exception_substrings)
    ("write_table(", "output_format=", 6,
     ["output_utils.py",   # helper definition
      "def write_table",   # function definition line
      "#"]),               # inline comment only lines
    ("process_region_entropy(", "output_format=", 8,
     ["region_entropy_processor.py",  # definition site
      "def process_region_entropy"]),
    ("compute_and_write_gc_factors(", "output_format=", 8,
     ["gc_correction.rs",      # Rust definition (not Python)
      "pon/build.py",           # PON builder: intentionally TSV-only (no user-facing format choice)
      "def compute_and_write",  # definition
      "process_sample(",        # process_sample docstring / pon context
      "write_extraction_outputs() -"]),  # module docstring line
    ("write_extraction_outputs(", "output_format=", 8,
     ["def write_extraction_outputs",          # definition
      "write_extraction_outputs() -"]),        # module-level docstring reference
    ("run_unified_pipeline(", "output_format", 30,
     ["def run_unified_pipeline",  # PyO3 / Python definition
      "rust",                      # any rust path reference
      "run_unified_pipeline() -"]),  # docstring reference
]


def check_file(path: pathlib.Path, rules: list) -> list[tuple[str, int, str]]:
    """Return list of (file, line_num, description) for each violation."""
    violations = []
    lines = path.read_text(encoding="utf-8").split("\n")
    rel = str(path.relative_to(ROOT))

    for trigger, required, lookahead, exceptions in rules:
        for i, line in enumerate(lines):
            stripped = line.strip()
            # Skip blank, comments, definition sites
            if not stripped or stripped.startswith("#"):
                continue
            if any(exc in rel or exc in stripped for exc in exceptions):
                continue
            if trigger not in line:
                continue

            # Check the call block (this line + next N)
            block = "\n".join(lines[i : min(i + lookahead, len(lines))])
            if required not in block:
                violations.append((rel, i + 1, f"'{trigger}' missing '{required}'  →  {stripped[:80]}"))

    return violations


def main() -> int:
    all_violations: list[tuple[str, int, str]] = []
    missing_files: list[str] = []

    for rel in SCAN_FILES:
        p = ROOT / rel
        if not p.exists():
            missing_files.append(rel)
            continue
        all_violations.extend(check_file(p, RULES))

    if missing_files:
        print("⚠️  Missing files (not checked):")
        for f in missing_files:
            print(f"   {f}")
        print()

    if all_violations:
        print(f"❌  {len(all_violations)} output_format gap(s) found:\n")
        for file, lineno, desc in all_violations:
            print(f"  {file}:{lineno}  {desc}")
        print()
        print("Fix: add  output_format=<var>  and  compress=<var>  to each call.")
        return 1

    print(f"✅  All output-writing call sites pass output_format= and compress=  ({len(SCAN_FILES)} files scanned)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
