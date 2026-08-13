"""A stub block must not reference a name its script block does not.

SPLIT_MAF's input was renamed from `sample_ids` to `metas`. The script block and
the tag were updated; the **stub** block was not, and still read
`${sample_ids.join(' ')}`. Nothing here caught it:

* no unit test reaches a stub block -- it exists only for `nextflow -stub-run`;
* `test_split_maf_logic.py` exercises the *script* block;
* Groovy resolves `${...}` at task runtime, so the file parses fine and the
  process fails only when it actually executes.

It took a stub run on the cluster to surface, as `No such variable: sample_ids`.

Phrased as a subset check rather than a declaration check on purpose. Parsing
Nextflow's input syntax precisely -- `path(x)`, bare `path x`, multi-element
tuples -- produced false positives on seven modules, and a test that cries wolf
is worse than no test. Comparing the two blocks against each other needs no
parsing and targets the actual failure: a rename applied to one block and not
the other.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

pytestmark = pytest.mark.unit

_MODULES = sorted(
    (Path(__file__).resolve().parents[2] / "nextflow" / "modules").rglob("main.nf")
)

_AMBIENT = {"task", "params", "workflow", "meta", "prefix", "args"}


def _roots(text: str) -> set[str]:
    """Leading identifier of every ${...} interpolation."""
    out: set[str] = set()
    for expr in re.findall(r"\$\{([^}]*)\}", text):
        m = re.match(r"\s*([A-Za-z_][A-Za-z0-9_]*)", expr)
        if m:
            out.add(m.group(1))
    return out


def _blocks(source: str) -> tuple[str, str]:
    """(everything before `stub:`, the stub block)."""
    i = source.find("\n    stub:")
    return (source, "") if i == -1 else (source[:i], source[i:])


@pytest.mark.parametrize("path", _MODULES, ids=lambda p: p.parent.name)
def test_stub_references_nothing_the_script_does_not(path: Path):
    before, stub = _blocks(path.read_text())
    if not stub:
        pytest.skip("no stub block")

    orphans = sorted(_roots(stub) - _roots(before) - _AMBIENT)
    assert not orphans, (
        f"{path.parent.name}'s stub block interpolates {orphans}, which appear "
        "nowhere in its inputs or script. Almost always an input renamed in the "
        "script and missed here -- and because Groovy resolves it at runtime, "
        "only `nextflow -stub-run` would ever find it."
    )
