"""RUNALL must not take over queue selection or retry from the cluster config.

`queue`, `errorStrategy` and `maxRetries` come from the institutional config
that `-profile iris` pulls in over HTTPS at parse time. Nothing in this
repository contains them, so `grep` says they do not exist -- which is exactly
the mistake #57 made. It hardcoded a partition ladder and a retry policy in the
`withName` block on the belief that nothing retried.

`withName` outranks a generic `process` block regardless of include order, so
that silently:

* replaced a working 3-attempt policy with 1;
* ignored `--isolated` and `--partition`, which the cluster closure honours;
* pinned attempt 1 to `cmobic_short` alone rather than letting SLURM choose
  between both partitions.

What the pipeline should own is the resource *request*, and only that. The
escalation falls out of it: a growing time request makes SLURM pick the
partition, because attempt 2's 4h no longer fits a 3h queue.

This test pins the absence. A future reader -- human or agent -- who greps for
`maxRetries`, finds nothing, and "fixes" it will fail here and be sent to the
comment that explains why.
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
    depth = 0
    for j in range(text.index("{", start), len(text)):
        if text[j] == "{":
            depth += 1
        elif text[j] == "}":
            depth -= 1
            if depth == 0:
                return text[start : j + 1]
    raise AssertionError("unbalanced KREWLYZER_RUNALL block")


@pytest.mark.parametrize("directive", ["queue", "errorStrategy", "maxRetries"])
def test_cluster_policy_is_left_to_the_cluster_config(runall_block, directive):
    assert not re.search(rf"\b{directive}\s*=", runall_block), (
        f"KREWLYZER_RUNALL sets `{directive}`, which overrides the institutional "
        "config (-profile iris) because withName outranks a generic process "
        "block. iris already sets maxRetries=3, falls back to `ignore`, and "
        "selects the partition from a closure honouring --isolated/--partition. "
        "Setting it here silently replaces all of that. See the comment above "
        "the block; this is the mistake #57 made."
    )


def test_no_partition_name_is_hardcoded(runall_block):
    """Cluster-specific names do not belong in a public pipeline."""
    for name in ("cmobic_short", "cmobic_cpu", "cpushort", "gpushort"):
        assert name not in runall_block, (
            f"'{name}' is hardcoded in KREWLYZER_RUNALL. Partitions belong to "
            "the site config, and pinning one here breaks every other site."
        )


def test_the_resource_ramp_survives(runall_block):
    """The one escalation this file does own.

    A growing time request is what lets SLURM choose the partition: attempt 2
    asks 4h, which no longer fits a 3h queue, so it lands on the long one
    without anybody naming it.
    """
    assert "task.attempt" in runall_block, (
        "the memory/time ramp is gone. With the cluster config's maxRetries=3 "
        "these are live, and they are how a retry gets more room."
    )
    assert re.search(r"memory\s*=\s*\{\s*\d+\.GB \* task\.attempt", runall_block)
    assert re.search(r"time\s*=\s*\{\s*\d+\.h\s*\* task\.attempt", runall_block)
