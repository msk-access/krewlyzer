"""Every `params.*` a Nextflow module reads must be declared in the config.

Nextflow returns ``null`` for an undeclared parameter rather than failing, so a
module can read one for years and look fine — the `?: default` swallows it. Two
consequences:

* the parameter is invisible: absent from ``--help``, absent from the docs, and
  not something anyone knows they can set;
* it becomes a hard error the moment ``nextflow.enable.strict`` is switched on.

Three were undeclared when this was written: ``validate_min_samples``,
``silent`` and ``assay``.

These tests are pure text analysis — no Nextflow required, which matters because
it is not installed in CI either.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, Set

import pytest

_ROOT = Path(__file__).resolve().parents[2]
_NEXTFLOW = _ROOT / "nextflow"
_CONFIG = _NEXTFLOW / "nextflow.config"

pytestmark = pytest.mark.skipif(
    not _NEXTFLOW.is_dir(), reason="nextflow/ not present in this checkout"
)


def _declared() -> Set[str]:
    """Names assigned inside the top-level ``params { ... }`` block.

    Scoped to that block deliberately: ``manifest``, ``process`` and
    ``executor`` also use ``name = value``, and counting those would mask a
    genuinely missing parameter.
    """
    text = _CONFIG.read_text()
    block = re.search(r"^params\s*\{(.*?)^\}", text, re.M | re.S)
    assert block, "no top-level params block in nextflow.config"
    return set(re.findall(r"^\s{4}(\w+)\s*=", block.group(1), re.M))


def _referenced() -> Dict[str, Set[str]]:
    """``params.<name>`` reads, keyed by name, valued by the files reading it."""
    out: Dict[str, Set[str]] = {}
    for path in sorted(_NEXTFLOW.rglob("*.nf")):
        rel = path.relative_to(_ROOT).as_posix()
        for name in re.findall(r"params\.(\w+)", path.read_text()):
            out.setdefault(name, set()).add(rel)
    return out


def test_every_referenced_param_is_declared():
    declared = _declared()
    missing = {n: sorted(f) for n, f in _referenced().items() if n not in declared}
    assert not missing, (
        "params read by a module but absent from nextflow.config — they read as "
        "null today and become fatal under nextflow.enable.strict:\n"
        + "\n".join(f"  params.{n}  ← {', '.join(f)}" for n, f in missing.items())
    )


def test_user_facing_params_are_documented():
    """A parameter nobody can discover may as well not exist.

    Excludes nf-core boilerplate and internal plumbing, which are configuration
    rather than knobs a user turns.
    """
    internal = {
        "custom_config_base",
        "custom_config_version",
        "monochrome_logs",
        "test_data_base",
        "help",
    }
    docs = (_ROOT / "docs" / "nextflow" / "parameters.md").read_text()
    undocumented = sorted(p for p in _declared() - internal if p not in docs)
    assert (
        not undocumented
    ), f"declared but undocumented in docs/nextflow/parameters.md: {undocumented}"


def test_the_validation_modules_are_actually_called():
    """They existed, and nothing invoked them.

    `KREWLYZER_VALIDATE_COHORT` was a complete module — process, stub, versions
    — that no workflow included. The gather half of the contract gate had no
    inputs and never ran, which is indistinguishable from passing.
    """
    workflows = "\n".join(
        p.read_text()
        for p in _NEXTFLOW.rglob("*.nf")
        if "modules/local" not in p.as_posix()
    )
    assert (
        "KREWLYZER_VALIDATE_COHORT" in workflows
    ), "the cohort validation module is defined but never included or called"


def test_runall_emits_the_validation_artifacts():
    """`run-all` writes them; the module has to declare them or they are lost.

    They were produced into the work directory and discarded, which removed the
    only affordable route to cohort-level degeneracy checking.
    """
    runall = (_NEXTFLOW / "modules/local/krewlyzer/runall/main.nf").read_text()
    for artifact in ("*.fingerprint.json", "*.validation.json"):
        assert artifact in runall, f"runall does not emit {artifact}"


def test_container_tags_agree_across_modules():
    """One version per release; a stale tag runs the wrong code silently."""
    tags = set()
    for path in _NEXTFLOW.rglob("*.nf"):
        tags |= set(re.findall(r"krewlyzer:([0-9][\w.]*)", path.read_text()))
    assert (
        len(tags) <= 1
    ), f"modules reference more than one container tag: {sorted(tags)}"
