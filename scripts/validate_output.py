#!/usr/bin/env python3
"""Shim for ``krewlyzer validate-output``.

Present so the gate can be run straight from a checkout, matching
``check_output_format.py``. The implementation lives in the package
(``krewlyzer.validate``) because operators need it on machines that only have
the wheel installed.

    python scripts/validate_output.py RESULTS_DIR [--json-report out.json]

Exit codes: 0 satisfied, 1 contract violation, 2 structural (not comparable).
"""

import sys

import typer

from krewlyzer.validate.cli import validate_output

if __name__ == "__main__":
    app = typer.Typer(add_completion=False)
    app.command()(validate_output)
    sys.exit(app())
