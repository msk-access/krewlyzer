"""The resource ramps must not become dead code again.

Every `task.attempt` multiplier in `nextflow.config` was inert for the whole of
0.8.x and 0.9.0: `task.attempt` never exceeds 1 without a retry policy, and
there was none. The config read as though a sample that overran or ran out of
memory would get a second, larger attempt. Nothing did, and one failed task out
of a 16,552-sample run terminated the whole thing.

These pin the parts of that fix a reader cannot verify by eye, and that a
future edit could quietly undo -- there is no stub run or integration test in
CI, so nothing else executes this file at all.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

pytestmark = pytest.mark.unit

_CONFIG = Path(__file__).resolve().parents[2] / "nextflow" / "nextflow.config"


@pytest.fixture(scope="module")
def runall_block() -> str:
    text = _CONFIG.read_text()
    start = text.index("withName: 'KREWLYZER_RUNALL'")
    depth, i = 0, text.index("{", start)
    for j in range(i, len(text)):
        if text[j] == "{":
            depth += 1
        elif text[j] == "}":
            depth -= 1
            if depth == 0:
                return text[start : j + 1]
    raise AssertionError("unbalanced KREWLYZER_RUNALL block")


def test_a_retry_policy_exists_at_all(runall_block: str):
    """Without this, every `task.attempt` expression below is dead code."""
    assert "maxRetries" in runall_block, (
        "KREWLYZER_RUNALL has no maxRetries, so task.attempt is pinned at 1 and "
        "the memory/time/queue ramps can never fire."
    )
    assert "errorStrategy" in runall_block


def test_retry_is_limited_to_slurm_kill_signals(runall_block: str):
    """A corrupt BAM must not climb the ladder to fail again on a 7-day queue."""
    codes = re.search(r"task\.exitStatus in \[([0-9,\s]+)\]", runall_block)
    assert codes, "retry is not gated on exit status; any failure would escalate"
    got = {int(c) for c in codes.group(1).replace(" ", "").split(",") if c}
    assert got <= {137, 140, 143}, f"unexpected retry codes: {sorted(got)}"


def test_exhausting_the_ladder_ignores_rather_than_terminates(runall_block: str):
    """`ignore` is only safe because --expect reconciles the cohort afterwards.

    terminate/finish would stop a 16,552-sample run on one bad sample; `ignore`
    without the reconciliation would drop it silently. The pair is the design.
    """
    assert "'ignore'" in runall_block, (
        "KREWLYZER_RUNALL no longer ignores an exhausted sample. If this became "
        "terminate or finish, one bad sample stops the cohort; if the --expect "
        "reconciliation went away instead, ignoring becomes silent loss."
    )


def test_the_first_rung_fits_inside_its_partition(runall_block: str):
    """cmobic_short caps at 3h. Asking for more is rejected at submission."""
    assert "cmobic_short" in runall_block and "cmobic_cpu" in runall_block
    first = re.search(r"task\.attempt == 1 \? (\d+)\.h", runall_block)
    assert first and int(first.group(1)) <= 3, (
        "the first rung requests more than cmobic_short's 3h limit, so SLURM "
        "would refuse the job rather than run it"
    )


def test_cpushort_is_not_a_rung(runall_block: str):
    """8-CPU multi-hour jobs, and --qos=priority does not apply there."""
    assert "cpushort" not in runall_block
