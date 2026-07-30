"""The FSC ratio check must fire on the shapes it exists to catch.

The contract-gate tests only ever run ``fsc_gene_ratios_sum_to_one`` against
well-formed tables, so its failure path is never taken there: a mistake in the
comparison or the tolerance would look exactly like a passing check. These
construct the broken tables directly.

Writing them is what showed that an earlier version of the check asserted the
same relation twice. "The six ratios sum to 1" and "the original five sum to
``1 - ultra_long_ratio``" are one equation rearranged, so the second could
never fail on its own; it is now documented on ``KREVIEW_FSC_RATIOS`` instead
of being re-derived at runtime.
"""

from __future__ import annotations

import pandas as pd
import pytest

from krewlyzer.validate.checks import (
    FSC_GENE_RATIOS,
    KREVIEW_FSC_RATIOS,
    fsc_gene_ratios_sum_to_one,
)


def _table(**overrides) -> pd.DataFrame:
    """A table whose six ratios sum to 1, before any override.

    Values carry four decimals deliberately. The tolerance is derived from the
    precision the column appears to have been written at, so a fixture using
    one or two decimals would imply a tolerance wide enough to swallow the very
    deviations these tests inject.
    """
    base = {
        "ultra_short_ratio": [0.1001, 0.2001],
        "core_short_ratio": [0.2000, 0.1500],
        "mono_nucl_ratio": [0.3999, 0.3499],
        "di_nucl_ratio": [0.1500, 0.1500],
        "long_ratio": [0.1000, 0.1000],
        "ultra_long_ratio": [0.0500, 0.0500],
    }
    base.update(overrides)
    return pd.DataFrame(base)


def test_the_two_column_lists_are_consistent():
    """`FSC_GENE_RATIOS` is the kreview five plus ultra_long, in that order."""
    assert FSC_GENE_RATIOS == KREVIEW_FSC_RATIOS + ["ultra_long_ratio"]
    assert len(FSC_GENE_RATIOS) == 6


def test_a_wellformed_table_reports_nothing():
    assert fsc_gene_ratios_sum_to_one(_table()) == []


def test_ratios_that_do_not_sum_to_one_are_reported():
    """A channel dropped from the partition; the six no longer total 1."""
    problems = fsc_gene_ratios_sum_to_one(_table(mono_nucl_ratio=[0.2999, 0.2499]))
    assert problems, "six ratios summing to 0.9 was not reported"
    assert "sum to 1" in problems[0]
    assert "2 row(s)" in problems[0], f"row count not reported: {problems[0]}"


def test_a_single_bad_row_among_good_ones_is_found():
    """The check reports per row, not on the column mean.

    A table where one gene's ratios are broken and the rest are fine averages
    out to something close to 1; only a row-wise comparison sees it.
    """
    df = _table()
    df.loc[0, "long_ratio"] = 0.2000  # this row now sums to 1.1
    problems = fsc_gene_ratios_sum_to_one(df)
    assert problems, "a single broken row was averaged away"
    assert "1 row(s)" in problems[0], problems[0]


def test_the_legacy_five_relation_cannot_break_on_its_own():
    """Documents why there is no second check.

    ``legacy5 - (1 - ultra_long)`` and ``(legacy5 + ultra_long) - 1`` are the
    same quantity, so any table where the five deviate from
    ``1 - ultra_long_ratio`` is a table where the six deviate from 1. If a
    future change makes these genuinely independent -- a seventh channel, or a
    ratio computed against something other than ``total`` -- this test fails
    and the check needs splitting.
    """
    df = _table(mono_nucl_ratio=[0.2999, 0.2499])  # six no longer sum to 1
    six_deviation = (df[FSC_GENE_RATIOS].sum(axis=1) - 1.0).abs()
    legacy_deviation = (
        df[KREVIEW_FSC_RATIOS].sum(axis=1) - (1.0 - df["ultra_long_ratio"])
    ).abs()
    assert (six_deviation - legacy_deviation).abs().max() < 1e-12, (
        "the two relations have become independent; fsc_gene_ratios_sum_to_one "
        "now needs to assert both"
    )


def test_a_missing_sixth_column_is_left_to_the_schema_check():
    """Pre-0.9.0 output has five ratios and no ultra_long_ratio.

    Reporting it here as well would double-count: the schema layer already
    names the absent column, and this check cannot say anything useful about a
    table whose partition it cannot see.
    """
    df = _table().drop(columns=["ultra_long_ratio"])
    assert fsc_gene_ratios_sum_to_one(df) == []


def test_an_empty_table_is_not_an_error():
    """No rows means nothing to contradict, not a violation."""
    assert fsc_gene_ratios_sum_to_one(_table().iloc[0:0]) == []


@pytest.mark.parametrize("decimals", [4, 6])
def test_rounding_at_the_written_precision_is_tolerated(decimals):
    """Ratios are written rounded; the tolerance must absorb that.

    A check that flags ordinary 4dp rounding would be silenced within a week,
    which is worse than not having it.
    """
    df = _table().round(decimals)
    assert fsc_gene_ratios_sum_to_one(df) == []
