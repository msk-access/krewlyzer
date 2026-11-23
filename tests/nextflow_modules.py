"""Read the pipeline's Nextflow modules as text: their declared outputs and their script blocks.

Nextflow is installed neither here nor in CI, so the tests that hold the tool-level modules to the
tools (`tests/unit/test_tool_level_modules.py`, `tests/invariants/test_tool_level_outputs.py`)
read the module files. The parsing covers only the forms these modules use, `path("...")` on one
output line with its `emit:` and `optional:`; anything else fails `outputs` loudly rather than
being skipped.
"""

from __future__ import annotations

import fnmatch
import re
from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple

_NEXTFLOW = Path(__file__).resolve().parents[1] / "nextflow"
MODULES = _NEXTFLOW / "modules" / "local" / "krewlyzer"
TOOL_LEVEL = _NEXTFLOW / "subworkflows" / "local" / "tool_level" / "main.nf"

#: The tool-level modules, each with the `krewlyzer` command it runs.
TOOL_MODULES: Tuple[Tuple[str, str], ...] = (
    ("extract", "extract"),
    ("fsc", "fsc"),
    ("fsr", "fsr"),
    ("fsd", "fsd"),
    ("wps", "wps"),
    ("ocf", "ocf"),
    ("motif", "motif"),
    ("region_entropy", "region-entropy"),
    ("region_mds", "region-mds"),
    ("mfsd", "mfsd"),
    ("uxm", "uxm"),
)


@dataclass(frozen=True)
class Output:
    """One declared output file pattern of a module."""

    glob: str
    emit: str
    optional: bool

    def matches(self, filename: str) -> bool:
        return any(fnmatch.fnmatchcase(filename, g) for g in expand_braces(self.glob))


def text(module: str) -> str:
    return (MODULES / module / "main.nf").read_text()


def _section(source: str, start: str, end: str) -> str:
    i = source.index(f"\n    {start}:")
    j = source.index(f"\n    {end}:", i)
    return source[i:j]


def script(module: str) -> str:
    """The `script:` block, up to `stub:`."""
    return _section(text(module), "script", "stub")


def outputs(module: str) -> List[Output]:
    """Every file pattern in the `output:` block (the `versions.yml` line excluded)."""
    found = []
    for line in _section(text(module), "output", "when").splitlines():
        line = line.strip()
        if not line or line.startswith("//") or line == "output:":
            continue
        globs = re.findall(r'path\s*\(?\s*"([^"]+)"', line)
        if len(globs) != 1:
            raise ValueError(f"{module}: cannot read output line {line!r}")
        if globs[0] == "versions.yml":
            continue
        emit = re.search(r"emit:\s*(\w+)", line)
        if emit is None:
            raise ValueError(f"{module}: output without emit: {line!r}")
        optional = re.search(r"optional:\s*true", line) is not None
        found.append(Output(globs[0], emit.group(1), optional))
    return found


def expand_braces(glob: str) -> List[str]:
    """`*.X.{a,b}` as `[*.X.a, *.X.b]`; one level of braces, which is all the modules use."""
    match = re.search(r"\{([^{}]*)\}", glob)
    if match is None:
        return [glob]
    head, tail = glob[: match.start()], glob[match.end() :]
    return [
        g for alt in match.group(1).split(",") for g in expand_braces(head + alt + tail)
    ]


def flags(source: str) -> set:
    """Every `--long-option` written in a block of script."""
    return set(re.findall(r"(?<![\w-])(--[a-z][a-z0-9-]*)", source))


def publish_dir(module: str) -> str:
    """The module's `publishDir` line, or '' when it has none."""
    found = re.findall(r"^\s*publishDir\b.*$", text(module), re.M)
    return found[0].strip() if found else ""
