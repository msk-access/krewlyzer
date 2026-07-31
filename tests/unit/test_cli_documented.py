"""Every registered CLI command must appear in the CLI reference.

Four shipped undocumented before this existed — `validate-output`,
`validate-cohort`, `describe-output` and `report`. Each had a CHANGELOG entry
and a `--help` string, and neither is where anyone looks first.

The failure is quiet by construction: a command that works is a command nobody
files a bug about, so the gap only surfaces when someone needed the feature,
could not find it, and did without.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

_ROOT = Path(__file__).resolve().parents[2]
_CLI = _ROOT / "src" / "krewlyzer" / "cli.py"
_REFERENCE = _ROOT / "docs" / "cli" / "index.md"

#: `run-all` has a page of its own rather than an entry in the index table.
_DOCUMENTED_ELSEWHERE = {"run-all"}


def _registered() -> set[str]:
    """Command names as typer sees them.

    Read from the source rather than by importing the app: importing pulls in
    the Rust extension, and a documentation check should not depend on whether
    the extension happens to be built.
    """
    text = _CLI.read_text()
    return set(re.findall(r'app\.command\(name="([^"]+)"\)', text))


def test_the_registry_is_not_empty():
    """Guards the regex above -- a rename would otherwise pass vacuously."""
    found = _registered()
    assert len(found) >= 5, f"only found {found}; has the registration style changed?"


@pytest.mark.parametrize("command", sorted(_registered()))
def test_each_command_appears_in_the_cli_reference(command):
    if command in _DOCUMENTED_ELSEWHERE:
        pytest.skip(f"{command} has a dedicated page")
    text = _REFERENCE.read_text()
    assert f"`{command}`" in text, (
        f"`krewlyzer {command}` is registered but absent from "
        f"docs/cli/index.md. A command nobody can find is a command nobody "
        f"uses."
    )


def test_commands_that_read_output_explain_their_exit_codes():
    """A workflow branches on these, so they are part of the interface.

    `2` means structural and is worth retrying; `1` means the contract was
    violated and is not.
    """
    text = _REFERENCE.read_text()
    section = text[text.index("validate-output") :]
    for code in ("`0`", "`1`", "`2`"):
        assert code in section, f"exit code {code} is not documented"


def test_the_report_command_is_marked_internal():
    """It contains one patient's measurements.

    Someone reading the docs should learn that before they publish one, not
    after.
    """
    text = _REFERENCE.read_text()
    section = text[text.index("### `report`") : text.index("### `validate-output`")]
    assert "Internal use" in section
    assert "not published" in section or "not committed" in section


def test_the_pon_caveat_is_documented_for_the_report():
    """Without a PON most of the interpretable surface is simply absent."""
    text = _REFERENCE.read_text()
    section = text[text.index("### `report`") : text.index("### `validate-output`")]
    assert "panel of normals" in section.lower()
    assert "--pon-model" in section
