"""Constructs Nextflow 26's config parser refuses.

26.04 rejected `nextflow.config` outright:

    nextflow.config:158:1: Try-catch blocks cannot be mixed with config statements

Not a warning and not a runtime failure -- the file would not parse, so the
pipeline could not start at all. 25.10.3 accepted the same file, which is how it
survived: the only signal was a version nobody had run yet.

This pins the one construct we know 26 forbids. It is not a parser and cannot
claim the file is 26-clean -- the parser stops at the first error, so there may
be more behind any given one. The only way to know is to run
`NXF_VER=26.04.6 nextflow ... -stub-run` on the cluster, and this test exists so
a fix already paid for is not quietly undone between those runs.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

pytestmark = pytest.mark.unit

_CONFIGS = sorted((Path(__file__).resolve().parents[2] / "nextflow").rglob("*.config"))


@pytest.mark.parametrize("path", _CONFIGS, ids=lambda p: p.name)
def test_no_try_catch_in_config(path: Path):
    assert not re.search(r"^\s*try\s*\{", path.read_text(), re.M), (
        f"{path.name} uses a try/catch block. Nextflow 26 refuses to parse the "
        "file -- 'Try-catch blocks cannot be mixed with config statements' -- so "
        "the pipeline cannot start. Use a ternary instead, as nf-core's own "
        "configs now do:\n"
        "    includeConfig !System.getenv('NXF_OFFLINE') && params.x ? \"...\" "
        ': "/dev/null"'
    )


@pytest.mark.parametrize("path", _CONFIGS, ids=lambda p: p.name)
def test_no_bare_environment_variable(path: Path):
    """26 rejects `${HOME}`; it wants `env('HOME')`.

    The second error 26 reported, after the try/catch. Found by grepping the
    pattern rather than waiting for the parser to surface them one run at a
    time -- it stops at the first, so a fix-and-rerun loop costs a cluster
    round trip per instance.

    `env()` is safe on 25: iris.config uses it five times and parses there.
    """
    bare = re.findall(r"\$\{([A-Z][A-Z0-9_]*)\}", path.read_text())
    assert not bare, (
        f"{path.name} interpolates the environment variable(s) {sorted(set(bare))} "
        'directly. Nextflow 26 will not parse the file -- "`HOME` is not '
        "defined (hint: use env('...'))\". Write ${env('HOME')} instead, which "
        "both 25 and 26 accept."
    )


def test_the_institutional_include_survives():
    """The ternary must still actually include the config in the normal case.

    Replacing try/catch with something that never includes anything would make
    this test file pass while removing the queue closure, the retry policy and
    the process_single label -- a far worse outcome than the parse error.
    """
    text = (
        Path(__file__).resolve().parents[2] / "nextflow/nextflow.config"
    ).read_text()
    assert (
        "nfcore_custom.config" in text
    ), "the institutional config is no longer included"
    assert "NXF_OFFLINE" in text, (
        "the include is unconditional. Skipping should require the operator to "
        "say so, not happen on a transient network failure."
    )
