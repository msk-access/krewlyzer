"""A z-score of 10^11 is a broken divisor, not a finding about biology.

Every fabricated-sigma defect in 0.9.0 produced *plausible* numbers and
survived every schema check: a hardcoded 167.0/5.0, six sigma floors, a
standard normal substituted for a missing block. This catches the opposite
failure -- a sigma so small the z explodes -- which no schema notices either,
because the column is present, typed and finite.

Measured for scale: real WPS z-scores against the 0.9.0 models reach 6.7, and
34,572 anchors over 6.9M positions produced nothing above 10. The threshold sits
at 100 so it flags arithmetic that has gone wrong by orders of magnitude, not
interesting outliers.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from krewlyzer.validate.checks import IMPLAUSIBLE_Z, plausible_z_scores


def test_ordinary_z_scores_pass():
    frame = pd.DataFrame({"mds_z": [-2.4, 0.1, 1.9, 6.7], "gene": list("abcd")})
    assert plausible_z_scores(frame) == []


def test_an_exploded_z_is_caught_and_named():
    frame = pd.DataFrame({"mds_z": [0.4, 1.2, 3.6e11]})
    problems = plausible_z_scores(frame)
    assert len(problems) == 1
    assert "mds_z" in problems[0]
    assert "3.6e+11" in problems[0] or "3.6e11" in problems[0]


def test_the_threshold_is_where_it_says_it_is():
    assert plausible_z_scores(pd.DataFrame({"a_z": [IMPLAUSIBLE_Z - 0.1]})) == []
    assert plausible_z_scores(pd.DataFrame({"a_z": [IMPLAUSIBLE_Z + 0.1]}))


def test_vector_z_columns_are_checked_element_wise():
    """`wps_nuc_z` is 200 values per row -- a near-zero sigma hides per position.

    Checking only scalars would miss the exact column the WPS vector baseline
    produces, which is the one with a sigma per position.
    """
    good = [np.array([0.1, -0.4, 2.2]), np.array([1.1, 0.0, -3.0])]
    assert plausible_z_scores(pd.DataFrame({"wps_nuc_z": good})) == []

    bad = [np.array([0.1, -0.4, 2.2]), np.array([1.1, 4.6e9, -3.0])]
    problems = plausible_z_scores(pd.DataFrame({"wps_nuc_z": bad}))
    assert problems and "wps_nuc_z" in problems[0]


def test_z_score_the_column_name_used_by_region_entropy():
    """TFBS/ATAC call it `z_score`, not `*_z`."""
    assert plausible_z_scores(pd.DataFrame({"z_score": [1.0, 2.0]})) == []
    assert plausible_z_scores(pd.DataFrame({"z_score": [1.0, 1e7]}))


def test_all_nan_is_not_a_violation():
    """An absent z is the honest outcome this release added everywhere."""
    assert plausible_z_scores(pd.DataFrame({"mds_z": [np.nan, np.nan]})) == []


def test_infinities_do_not_crash_the_check():
    """`x/0.0` gives inf, not a large finite number.

    Left to the non-finite filter deliberately: an inf is already unusable
    downstream and the contract's own NaN checks speak to it. What this check
    exists for is the *finite* absurdity that looks like data.
    """
    assert plausible_z_scores(pd.DataFrame({"mds_z": [np.inf, -np.inf]})) == []


def test_columns_that_merely_end_in_z_are_not_mistaken_for_z_scores():
    """`chrom_z` would be a naming accident, but `alu_count` must not match."""
    frame = pd.DataFrame({"alu_count": [5_000_000], "region_sz": [1e9]})
    assert plausible_z_scores(frame) == []


def test_the_check_runs_on_every_table_not_a_listed_few():
    """A new table must not be able to opt out by omission."""
    from krewlyzer.validate import gate

    source = open((gate.__file__ or "").replace(".pyc", ".py")).read()
    assert 'for name in (*rule.checks, "plausible_z_scores")' in source
    # ...and its columns are read even when the contract does not declare them.
    assert 'c.endswith("_z") or c == "z_score"' in source
