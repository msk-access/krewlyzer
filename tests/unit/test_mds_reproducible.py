"""MDS must not change between two runs on identical counts.

`compute_mds` summed the entropy over a Python `set` of k-mer strings. A set of strings iterates
in the per-process hash-seed order (PYTHONHASHSEED), so the floating-point sum ran in a different
order every run: identical end-motif counts gave MDS values differing in the last digit, and
`mds_z` differed run to run (seen on four real samples, relative 6e-14). Small, but it breaks
every byte-identical comparison between two runs of the same code.
"""

import itertools
import os
import random
import subprocess
import sys
from pathlib import Path

import krewlyzer
from krewlyzer.core.motif_processor import compute_mds

_SRC = str(Path(krewlyzer.__file__).resolve().parents[1])
_SCRIPT = """
import itertools, random
from krewlyzer.core.motif_processor import compute_mds
r = random.Random(7)
counts = {"".join(k): r.randint(1, 10**6) for k in itertools.product("ACGT", repeat=4)}
print(repr(float(compute_mds(counts))))
"""


def _mds_under_hash_seed(seed: int) -> str:
    env = {**os.environ, "PYTHONHASHSEED": str(seed), "PYTHONPATH": _SRC}
    done = subprocess.run(
        [sys.executable, "-c", _SCRIPT],
        env=env,
        capture_output=True,
        text=True,
        check=True,
    )
    return done.stdout.strip()


def test_mds_is_the_same_whatever_the_hash_seed():
    values = {_mds_under_hash_seed(seed) for seed in range(8)}
    assert len(values) == 1, f"MDS depends on the hash seed: {sorted(values)}"


def test_mds_does_not_depend_on_the_order_counts_arrive_in():
    """Extract's counts come from a Rust HashMap, whose order is random per process too."""
    r = random.Random(11)
    counts = {
        "".join(k): r.randint(1, 10**6) for k in itertools.product("ACGT", repeat=4)
    }
    reordered = dict(reversed(list(counts.items())))
    assert compute_mds(counts) == compute_mds(reordered)


def test_different_counts_give_a_different_mds():
    uniform = {"".join(k): 100 for k in itertools.product("ACGT", repeat=4)}
    skewed = {**uniform, "CCCA": 10_000}
    assert compute_mds(uniform) != compute_mds(skewed)
