"""An unknown `params.genome` stops the pipeline at launch.

Every module passes `--genome ${params.genome}` through, and the CLI refuses a genome it has
no assets for. So `--genome hg99` failed every extract task, which the retry-then-ignore error
strategy retried and then dropped: a run that "finished" with no samples. The accepted builds
are listed once, as an enum in `nextflow_schema.json`, and `main.nf` checks `params.genome`
against that list before any task is submitted.

Text analysis, like the other Nextflow tests here: Nextflow is not installed in CI.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

from krewlyzer.assets import Genome, genome_build

_ROOT = Path(__file__).resolve().parents[2]
_NEXTFLOW = _ROOT / "nextflow"
_SCHEMA = _NEXTFLOW / "nextflow_schema.json"
_MAIN = _NEXTFLOW / "main.nf"
_CONFIG = _NEXTFLOW / "nextflow.config"

pytestmark = pytest.mark.skipif(
    not _NEXTFLOW.is_dir(), reason="nextflow/ not present in this checkout"
)


def _genome_property() -> dict:
    return json.loads(_SCHEMA.read_text())["properties"]["genome"]


def test_the_schema_lists_the_genomes_the_cli_accepts():
    enum = _genome_property().get("enum")
    assert enum, "nextflow_schema.json: `genome` has no enum"
    assert sorted(enum) == sorted(g.value for g in Genome)
    for name in enum:
        genome_build(name)  # raises for a build the CLI would refuse


def test_the_default_genome_is_in_the_enum():
    prop = _genome_property()
    default = re.search(r"^\s+genome\s*=\s*'([^']+)'", _CONFIG.read_text(), re.M)
    assert default, "nextflow.config declares no params.genome"
    assert prop["default"] == default.group(1)
    assert default.group(1) in prop["enum"]


def test_main_checks_params_genome_against_the_schema_enum_before_any_task():
    """The check reads the schema's list rather than restating it, so the two cannot drift."""
    text = _MAIN.read_text()
    entry = text[text.index("workflow NFCORE_KREWLYZER") :]
    check = entry[: entry.index("KREWLYZER(")]
    assert (
        "nextflow_schema.json" in check
    ), "main.nf does not read the schema's genome list"
    assert "['properties']['genome']['enum']" in check
    assert "params.genome" in check and "System.exit(1)" in check
