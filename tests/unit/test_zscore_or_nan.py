"""An unmeasurable spread yields no z-score, not a zero one.

Ten baseline classes in `pon/model.py` independently ended with
`if std > 0: return (x - mean) / std` followed by `return 0.0`. Nine copies of
the same three lines, and the same defect nine times over: zero is the single
most confident statement that column can make — "this sample sits exactly at
the healthy baseline" — asserted precisely when the baseline measured no spread
at all, and indistinguishable from a genuine zero.

Same reasoning that removed the σ floors from the builder (`4cd634b`) and
`z_score = 0.0` from region entropy. `zscore_or_nan` is now the one place it
lives, so it can only be wrong once.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from krewlyzer.pon.model import (
    FsdBaseline,
    MdsBaseline,
    OcfBaseline,
    RegionMdsBaseline,
    zscore_or_nan,
)

pytestmark = pytest.mark.unit


# ---------------------------------------------------------------------------
# the helper
# ---------------------------------------------------------------------------


def test_a_real_spread_gives_a_real_z():
    assert zscore_or_nan(12.0, 10.0, 2.0) == pytest.approx(1.0)
    assert zscore_or_nan(8.0, 10.0, 2.0) == pytest.approx(-1.0)


@pytest.mark.parametrize(
    "std,why",
    [
        (0.0, "no spread measured"),
        (-1.0, "negative sigma is not a spread"),
        (float("nan"), "the builder writes NaN where it could not measure"),
        (float("inf"), "infinite sigma is not a measurement either"),
        (None, "absent"),
    ],
)
def test_an_unusable_sigma_gives_nan(std, why):
    assert math.isnan(zscore_or_nan(12.0, 10.0, std)), why


def test_a_tiny_but_real_sigma_is_honoured():
    """The counterpart: no floor sneaks back in.

    `region_mds` sigmas run ~0.007 and `MDS` whole-sample sigmas ~0.0008. A
    guard that rejected "small" sigmas would silence the most precise
    baselines in the model.
    """
    assert zscore_or_nan(0.9530, 0.9524, 0.0008) == pytest.approx(0.75, abs=0.01)


def test_a_missing_mean_or_observation_gives_nan():
    assert math.isnan(zscore_or_nan(float("nan"), 10.0, 2.0))
    assert math.isnan(zscore_or_nan(12.0, float("nan"), 2.0))


# ---------------------------------------------------------------------------
# the classes that use it
# ---------------------------------------------------------------------------


def test_every_baseline_class_agrees_on_an_unusable_sigma():
    """One helper, so they cannot drift apart again."""
    import pandas as pd

    ocf = OcfBaseline(
        regions=pd.DataFrame({"region_id": ["T"], "ocf_mean": [50.0], "ocf_std": [0.0]})
    )
    mds = MdsBaseline(kmer_expected={}, kmer_std={}, mds_mean=0.95, mds_std=0.0)
    gene = RegionMdsBaseline(gene_baseline={"TP53": {"mds_mean": 0.9, "mds_std": 0.0}})

    assert math.isnan(ocf.compute_zscore("T", 50.0))
    assert math.isnan(mds.get_mds_zscore(0.95))
    assert math.isnan(gene.compute_zscore("TP53", 0.9))


def test_an_absent_key_is_still_none_not_nan():
    """The two absences stay distinguishable.

    `None` means "not in this baseline" and `NaN` means "in it, but it measured
    nothing". Both write NaN to the output column, correctly — but a caller
    that collapses them reports "0 regions matched" for regions that matched
    perfectly well, which is the difference between "rebuild the PON" and
    "this region has no variance in the cohort".
    """
    gene = RegionMdsBaseline(gene_baseline={"TP53": {"mds_mean": 0.9, "mds_std": 0.01}})
    assert gene.compute_zscore("NOT_A_GENE", 0.9) is None
    assert gene.compute_zscore("TP53", 0.9) is not None


def test_an_unknown_arm_has_no_expectation():
    """`get_std` returning 0.0 would give a caller infinity, not an absence."""
    fsd = FsdBaseline(size_bins=[65, 70], arms={})
    assert math.isnan(fsd.get_expected("chr99:1-2", 150))
    assert math.isnan(fsd.get_std("chr99:1-2", 150))


def test_a_wps_vector_z_is_nan_where_the_baseline_could_not_measure():
    """The read side must not undo the builder's honesty.

    Since 4cd634b the builder writes NaN for positions with no measurable
    spread. `np.where(std > 0, std, 1.0)` is False for NaN, so a 1.0 default
    turned every one of those into a finite z — silently, at the harder place
    to notice.
    """
    import pandas as pd

    from krewlyzer.pon.model import WpsBaseline

    mean = [10.0, 10.0, 10.0]
    std = [2.0, float("nan"), 0.0]
    baseline = WpsBaseline(
        regions=pd.DataFrame(
            {
                "region_id": ["R"],
                "wps_nuc_mean": [mean],
                "wps_nuc_std": [std],
                "wps_tf_mean": [mean],
                "wps_tf_std": [std],
            }
        ),
        schema_version="2.0",
    )
    z = baseline.compute_z_vector("R", np.array([14.0, 14.0, 14.0]), "wps_nuc")
    assert z is not None
    assert z[0] == pytest.approx(2.0), "a measured position must still score"
    assert np.isnan(z[1]), "NaN sigma must stay NaN, not become 4.0"
    assert np.isnan(z[2]), "zero sigma must be NaN, not 4.0"
