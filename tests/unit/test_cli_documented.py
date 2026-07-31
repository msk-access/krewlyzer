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

#: How many commands the CLI has, near enough. The point is not the exact
#: figure but that the parser below sees *most* of them: its first version
#: matched only `app.command(name="x")` and silently skipped the eleven
#: registered as bare `app.command()(fn)`, so it checked 8 of 19 and reported
#: a clean pass. A coverage test that under-collects is worse than none — it
#: converts an unknown into a false assurance.
_MINIMUM_EXPECTED = 15


def _registered() -> set[str]:
    """Command names as typer sees them.

    Read from the source rather than by importing the app: importing pulls in
    the Rust extension, and a documentation check should not depend on whether
    the extension happens to be built.

    Both registration forms count. `app.command()(extract)` takes the function
    name with underscores turned into hyphens, which is how typer derives it.
    """
    text = _CLI.read_text()
    named = re.findall(r'app\.command\(name="([^"]+)"\)\(\w+\)', text)
    plain = re.findall(r"app\.command\(\)\((\w+)\)", text)
    return set(named) | {fn.replace("_", "-") for fn in plain}


def test_the_parser_sees_the_whole_cli():
    """Guards the regexes above; they have already been wrong once.

    Cross-checked against typer itself rather than a hand-written list, so the
    two cannot drift. Skipped when the Rust extension is unbuilt, since that
    makes the import — not the CLI — the thing that failed.
    """
    found = _registered()
    assert len(found) >= _MINIMUM_EXPECTED, (
        f"only found {sorted(found)}; the registration style has changed and "
        "this parser no longer sees the whole CLI"
    )

    typer_app = pytest.importorskip(
        "krewlyzer.cli", reason="Rust extension not built"
    ).app
    actual = {
        c.name or (c.callback.__name__.replace("_", "-"))
        for c in typer_app.registered_commands
    }
    assert found == actual, (
        f"parser and typer disagree — parser only: {sorted(found - actual)}, "
        f"typer only: {sorted(actual - found)}"
    )


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


# ---------------------------------------------------------------------------
# README
# ---------------------------------------------------------------------------
#
# The README is the only page most people read, and it is the one page no check
# covered. Its links are absolute `msk-access.github.io` URLs, so the mkdocs
# build never resolves them and the internal link checker cannot see them —
# six pointed at a docs layout that had been reorganised out of existence.

_README = _ROOT / "README.md"
_DOCS = _ROOT / "docs"
_SITE = "https://msk-access.github.io/krewlyzer/"


def _published_slugs() -> set[str]:
    """Every URL path mkdocs will serve, from the pages that exist."""
    slugs = set()
    for page in _DOCS.rglob("*.md"):
        slug = page.relative_to(_DOCS).with_suffix("").as_posix()
        slugs.add(slug)
        if slug.endswith("/index"):
            slugs.add(slug[: -len("/index")])
        elif slug == "index":
            slugs.add("")
    return slugs


@pytest.mark.parametrize(
    "url",
    sorted(set(re.findall(rf"{re.escape(_SITE)}([^)\s]*)", _README.read_text()))),
)
def test_every_readme_documentation_link_resolves(url):
    slugs = _published_slugs()
    assert url.strip("/") in slugs or url.strip("/") == "", (
        f"README links to {_SITE}{url}, which mkdocs does not publish. "
        f"A 404 on the landing page is the first thing a new reader sees."
    )


@pytest.mark.parametrize("command", sorted(_registered()))
def test_each_command_appears_in_the_readme(command):
    text = _README.read_text()
    assert f"`{command}`" in text, (
        f"`krewlyzer {command}` is registered but absent from README.md. "
        "The README's command table is where people look before the docs site."
    )
