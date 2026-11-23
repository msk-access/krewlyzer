"""The synthetic run-all keeps every contig and excludes exactly its own blacklist window.

Without an explicit blacklist, run-all takes the bundled hg19 one. Its first entry,
chr1:0-750100, holds the whole 300 kb synthetic chr1, so every run-all and CLI test ran
on half the sample, and each chr1 row compared zeros with zeros. The harness now passes
`AssetBundle.blacklist`, one 10 kb window. This holds both halves of that: no contig is
lost, and the exclusion still bites, on the fragments in that window and no others.
"""

from __future__ import annotations

import gzip

import pytest

from conftest import requires_data

from .runall import synth_run_all
from .synth import BASELINE, BLACKLIST_WINDOW, make_sample

pytestmark = [pytest.mark.invariant, requires_data]


def _intervals(lines):
    return {
        (c, int(s), int(e))
        for c, s, e, *_ in (line.rstrip("\n").split("\t") for line in lines)
    }


def test_run_all_excludes_the_synthetic_window_and_nothing_else(tmp_path):
    sample = make_sample(tmp_path / "in", BASELINE)
    out = tmp_path / "out"
    assert synth_run_all(sample, out) is None

    truth = _intervals(open(sample.truth_bed))
    kept = _intervals(gzip.open(out / "S.bed.gz", "rt"))
    chrom, start, end = BLACKLIST_WINDOW
    in_window = {f for f in truth if f[0] == chrom and f[1] < end and f[2] > start}

    assert in_window, "no fragment in the blacklist window: the exclusion is vacuous"
    assert {c for c, _, _ in kept} == {c.name for c in sample.contigs}
    assert kept == truth - in_window
