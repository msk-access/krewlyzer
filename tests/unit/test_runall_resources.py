"""Which resource decisions this pipeline owns, and which belong to the cluster.

Two rules, learned the expensive way.

**Retry policy is the cluster's.** `errorStrategy` and `maxRetries` come from
the institutional config `-profile iris` pulls in over HTTPS at parse time, so
`grep` over this repository says they do not exist. #57 believed that and
hardcoded its own, and because `withName` outranks a generic `process` block it
silently replaced a working 3-attempt policy with 1.

**Resource requests are ours, and must be set with `withName`.** The `withLabel`
blocks in `nextflow.config` are inert under `-profile iris`: the institutional
config is included afterwards and defines the same labels, and at equal
specificity the later definition wins. SPLIT_MAF ran at 12 GB / 2h, never the
4 GB / 1h its label declares.

That stays invisible until a request is illegal. Under `EnforcePartLimits =
ALL`, a partition list enforces the *intersection* of the limits, so a 4h
request against a list containing a 3h partition is refused before it queues.
It cost two things on one 16,552-sample run: every RUNALL retry (fixed with a
per-attempt `queue`), and KREWLYZER_VALIDATE_COHORT on its *first* attempt --
four refusals, then ignored, and the cohort degeneracy report was never
produced with nothing to say so.
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


@pytest.mark.parametrize("directive", ["errorStrategy", "maxRetries"])
def test_retry_policy_is_left_to_the_cluster_config(runall_block, directive):
    """`queue` is deliberately NOT in this list any more.

    It was, when #57 was reverted, on the reasoning that the cluster's own queue
    closure would route a grown request to a partition that could take it. A
    real 16,552-sample run disproved that: with `EnforcePartLimits = ALL` a
    partition list enforces the intersection of the limits, so attempt 2's 4h
    against a list containing a 3h partition was refused outright --
    "Requested time limit is invalid" -- rather than rerouted.

    So the pipeline does now set `queue`, per attempt, which is what #57 got
    right. What it got wrong, and what this still guards, is overriding the
    retry policy: iris gives 3 attempts and a fallback, and a `withName` block
    outranks it silently.
    """
    assert not re.search(rf"\b{directive}\s*=", runall_block), (
        f"KREWLYZER_RUNALL sets `{directive}`, which overrides the institutional "
        "config (-profile iris) because withName outranks a generic process "
        "block. iris already sets maxRetries=3 and falls back to `ignore`; "
        "setting it here silently replaces that. This is the mistake #57 made."
    )


def test_the_queue_ladder_moves_retries_to_the_long_partition(runall_block):
    """The measured half of #57, kept.

    Attempt 1 uses `--partition`; a retry, whose time request has grown past the
    shortest partition's cap, uses `--long_partition` if given. Without this a
    retry is not rerouted but refused, and the sample is dropped.
    """
    assert "task.attempt > 1" in runall_block and "long_partition" in runall_block, (
        "the retry queue ladder is gone. Under EnforcePartLimits=ALL a grown "
        "time request against a partition list is refused, not rerouted, so "
        "every retry would fail to submit and the sample would be dropped."
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


# Processes the run-all path actually executes. A process missing from
# nextflow.config's withName blocks runs on the institutional config's numbers,
# because our withLabel blocks are overridden by it -- see below.
_RUNALL_PATH_PROCESSES = (
    "SPLIT_MAF",
    "FILTER_MAF",
    "KREWLYZER_RUNALL",
    "KREWLYZER_VALIDATE_PON",
    "KREWLYZER_VALIDATE_COHORT",
)


@pytest.mark.parametrize("process", _RUNALL_PATH_PROCESSES)
def test_every_runall_process_pins_its_own_time(process):
    """`withLabel` is inert under `-profile iris`; only `withName` overrides it.

    The institutional config is included after our `process` block and defines
    the same labels, and at equal selector specificity the later definition
    wins. So `withLabel:process_low` here never applied -- SPLIT_MAF ran at
    12 GB / 2h, not the 4 GB / 1h it declares.

    That is invisible until a request is illegal. iris's `process_single` asks
    4h, and under `EnforcePartLimits = ALL` a 4h request against a partition
    list containing a 3h partition is refused before it queues. On a completed
    16,552-sample run, KREWLYZER_VALIDATE_COHORT failed to submit on its FIRST
    attempt, four times, and was ignored -- the cohort degeneracy report was
    never produced and nothing said so.

    `--long_partition` cannot help: it moves retries, and that request was
    already illegal at attempt 1.
    """
    text = _CONFIG.read_text()
    block = re.search(rf"withName:\s*'{process}'\s*\{{(.*?)\n    \}}", text, re.S)
    assert block, (
        f"{process} has no withName block, so it inherits the institutional "
        "config's resources. Ours are overridden and cannot constrain it. If "
        "those numbers exceed the shortest partition in --partition, the "
        "process is refused at submission and -- with the cluster's `ignore` "
        "fallback -- its output silently never appears."
    )
    assert re.search(r"\btime\s*=", block.group(1)), (
        f"{process}'s withName block does not pin `time`. That is the directive "
        "that gets a job refused under EnforcePartLimits=ALL."
    )
