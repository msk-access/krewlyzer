#!/usr/bin/env python3
"""Regenerate the example output report in ``docs/``.

`krewlyzer describe-output` describes one sample's output directory. The docs
need a worked example of what it produces, and that example has to come from
**fabricated** tables.

A real sample's report cannot be committed. The docs site is public, and while
identifier columns are redacted, the report still carries that patient's actual
measurements — every fragment count, every motif frequency. Redaction removes
the name, not the data.

So the example is built from ``tests/invariants/synth_outputs``, the same
fabricated tables the contract-gate tests run against: correct shapes, correct
column names, invented numbers. It exercises the real code path, so the example
is a true rendering of the tool rather than a mock-up of one.

Usage
-----
    python scripts/build_example_report.py            # rewrite the example
    python scripts/build_example_report.py --check    # exit 1 if it is stale
"""

from __future__ import annotations

import argparse
import re
import sys
import tempfile
from pathlib import Path
from typing import List, Optional

_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_ROOT))
sys.path.insert(0, str(_ROOT / "tests"))

OUTPUT = _ROOT / "docs" / "reference" / "example-output-report.md"

#: The sanctioned placeholder. Never a real or realistic identifier -- an
#: invented-but-plausible one is indistinguishable from a real one to everyone
#: downstream of whoever invented it.
SAMPLE = "P-0000000-T01-XS1"

_PREAMBLE = """<!--
GENERATED FILE -- DO NOT EDIT.

Produced by scripts/build_example_report.py from fabricated tables. Edit the
generator or `krewlyzer describe-output`, then run:

    python scripts/build_example_report.py
-->

# Example output report

What [`krewlyzer describe-output`](../cli/index.md) produces for one sample.

```bash
krewlyzer describe-output RESULTS/{sample_id}/ -o report.html
```

!!! info "This example uses fabricated data"
    Every number below is invented. The tables have the right shapes, columns
    and dtypes, but the values carry no biology — this page exists to show the
    *form* of the report, not to be read as a result.

    A report from a real sample is not committed here. Identifier columns are
    redacted, but the measurements themselves are still that patient's data,
    and this site is public.

---

"""


def _write_fixture(directory: Path) -> Path:
    """Fabricate one sample's output tables and return its directory."""
    from invariants.synth_outputs import write_sample

    return write_sample(directory, SAMPLE, sample_idx=0)


def render() -> str:
    from krewlyzer.validate.describe import describe_sample, render_markdown

    with tempfile.TemporaryDirectory() as tmp:
        sample_dir = _write_fixture(Path(tmp))
        body = render_markdown(describe_sample(sample_dir))

    # The temp path is machine- and run-specific; keeping it would make the
    # file churn on every regeneration and leak a local directory layout.
    body = re.sub(r"^`/.*`$", f"`RESULTS/{SAMPLE}`", body, count=1, flags=re.M)
    # The generator's own title is redundant under the page heading above.
    body = re.sub(r"^# Output report.*\n+", "", body, count=1, flags=re.M)
    return _PREAMBLE + body


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument(
        "--check",
        action="store_true",
        help="Do not write; exit 1 if the example is missing or out of date",
    )
    args = ap.parse_args(argv)

    rendered = render()

    # Belt and braces. The fixture uses the placeholder, but this file is
    # destined for a public site and a regression in the redaction path would
    # otherwise publish an identifier before anyone noticed.
    leaked = re.findall(
        r"P-(?!0000000)0\d{6}-[A-Z]\d{2}|C-[0-9A-Z]{6}-[A-Z]\d{3}", rendered
    )
    if leaked:
        print(
            f"error: identifier(s) in the example report: {sorted(set(leaked))}",
            file=sys.stderr,
        )
        return 2

    if args.check:
        if not OUTPUT.exists():
            print(f"{OUTPUT.name} is missing; run scripts/build_example_report.py")
            return 1
        if OUTPUT.read_text(encoding="utf-8") != rendered:
            print(
                f"{OUTPUT.name} is out of date; run scripts/build_example_report.py",
                file=sys.stderr,
            )
            return 1
        print(f"{OUTPUT.name} is up to date")
        return 0

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(rendered, encoding="utf-8")
    print(f"wrote {OUTPUT.relative_to(_ROOT)} ({len(rendered.splitlines())} lines)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
