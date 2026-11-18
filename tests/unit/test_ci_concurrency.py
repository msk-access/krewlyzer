"""A branch push must never displace another run of the same workflow.

GitHub keeps at most one *pending* run per concurrency group (`queue: single`,
the default) and a newly queued run cancels it. With one group per branch, three
quick pushes to develop dropped the middle run. Each push run scans from its own
`github.event.before`, so nothing ever scanned that middle push's commits or
commit messages -- a hole in the PHI backstop, found in review on #86.

So for anything that is not a PR, the group must be unique per run. PR runs may
share a group and cancel each other: a PR run scans `base.sha..HEAD`, the whole
PR, so the surviving run covers what the cancelled one would have.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

pytestmark = pytest.mark.unit

_WORKFLOWS = Path(__file__).resolve().parents[2] / ".github" / "workflows"

#: Workflows whose runs must each complete: the PHI scan lives in the first, and
#: the second records what landed on develop and main.
_EVERY_RUN_MUST_COMPLETE = ("quick-checks.yml", "test.yml")


def _concurrency(name: str) -> dict[str, str]:
    """The workflow-level `concurrency` block, read without a YAML parser.

    Deliberately not `yaml.safe_load`: PyYAML is not a declared dependency, so
    `pytest.importorskip("yaml")` would skip -- green -- in any environment
    that happens not to have it, which is how a guard stops guarding. The
    block is two keys; reading it as text needs nothing.
    """
    text = (_WORKFLOWS / name).read_text(encoding="utf-8")
    m = re.search(r"^concurrency:\n((?:[ \t]+.*\n|[ \t]*\n)+)", text, re.M)
    assert m, f"{name} has no workflow-level concurrency block"
    keys = {}
    for line in m.group(1).splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        k, _, v = s.partition(":")
        keys[k.strip()] = v.strip()
    assert (
        "group" in keys and "cancel-in-progress" in keys
    ), f"{name}: concurrency block lacks group or cancel-in-progress: {keys}"
    return keys


@pytest.mark.parametrize("name", _EVERY_RUN_MUST_COMPLETE)
def test_a_branch_push_gets_a_group_of_its_own(name):
    group = _concurrency(name)["group"]
    assert "github.run_id" in group, (
        f"{name}: concurrency group {group!r} does not fall back to "
        "github.run_id, so branch pushes share a group and a queued run can "
        "be cancelled by the next one"
    )
    for shared in ("github.ref", "github.head_ref", "github.sha"):
        assert shared not in group, (
            f"{name}: concurrency group {group!r} keys on {shared}, which two "
            "pushes to the same branch can share"
        )


@pytest.mark.parametrize("name", _EVERY_RUN_MUST_COMPLETE)
def test_only_pr_runs_cancel_in_progress(name):
    cancel = _concurrency(name)["cancel-in-progress"]
    assert "pull_request" in cancel, (
        f"{name}: cancel-in-progress is {cancel!r}; it must be conditional on "
        "the event being a pull request"
    )
