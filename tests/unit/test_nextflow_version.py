"""The pipeline must report the version it actually is.

`manifest.version` in `nextflow.config` sat at 0.8.3 through the whole 0.9.0
release. Nothing functional depended on it, which is exactly why it drifted --
but it is what Nextflow writes into every execution report, trace file and
Tower entry. A 16,000-sample cohort would have been labelled 0.8.3 in its own
provenance, while the containers it ran were 0.9.0.

The Phase 2 version bump is a `sed` over container tags in the modules, and
that pattern never matched this line. Prose in the release guide cannot enforce
it; this can.

Deliberately checks the module container pins too. They are the other half of
the same claim -- what the pipeline says it is, and what it actually runs.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from krewlyzer import __version__

pytestmark = pytest.mark.unit

_ROOT = Path(__file__).resolve().parents[2]
_CONFIG = _ROOT / "nextflow" / "nextflow.config"
_MODULES = _ROOT / "nextflow" / "modules"


def test_manifest_version_matches_the_package():
    text = _CONFIG.read_text()
    match = re.search(r"^\s*version\s*=\s*'([^']+)'", text, re.M)
    assert match, "no manifest version found in nextflow.config"
    assert match.group(1) == __version__, (
        f"nextflow.config manifest.version is {match.group(1)!r} but the package "
        f"is {__version__!r}. Every execution report and trace from this "
        "pipeline would carry the wrong version."
    )


def test_every_module_container_matches_the_package():
    """A stale pin here runs the *previous release's code* under this version.

    Worse than the manifest: that is not a labelling error, it is executing
    defects a release removed. 0.8.3 containers under a 0.9.0 pipeline would
    reproduce every PON defect this release fixed.
    """
    stale: dict[str, set[str]] = {}
    for path in sorted(_MODULES.rglob("main.nf")):
        for tag in re.findall(r"krewlyzer:([0-9]+\.[0-9]+\.[0-9]+)", path.read_text()):
            if tag != __version__:
                stale.setdefault(str(path.relative_to(_ROOT)), set()).add(tag)

    assert not stale, (
        f"module container tags disagree with the package version {__version__}: "
        + "; ".join(f"{p} pins {sorted(v)}" for p, v in stale.items())
    )
