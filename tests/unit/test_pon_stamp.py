"""Stamping a PON asserts it is fit to score against, so it has to be earned.

A PON is built from `develop`, where the version still reads the previous
release — so the model records that, however new the code is. Restamping at
release time is the alternative to bumping the version before a four-hour
build.

The risk that creates is obvious once stated: `stamp-pon --version 0.9.0` on
one of the four models carrying the fabricated `167.0 / 5.0` baseline would
claim exactly the compatibility the guard exists to deny. So stamping runs
`validate-pon` first and refuses on failure.
"""

from __future__ import annotations

import pandas as pd
import pytest

from krewlyzer.validate.pon_stamp import stamp_release

pytestmark = pytest.mark.unit


def _pon(tmp_path, *, nrl_std, version="0.8.3", name="t.pon.parquet"):
    meta = pd.DataFrame(
        [
            {
                "table": "metadata",
                "schema_version": "1.0",
                "assay": "xs2",
                "build_date": "2026-08-03",
                "n_samples": 21.0,
                "reference": "r",
                "panel_mode": True,
                "target_regions_file": "t",
                "krewlyzer_version": version,
                "cohort_digest": "abc123def4567890",
                "cohort_label": "healthy-donors",
            }
        ]
    )
    block = pd.DataFrame(
        {
            "table": "wps_background",
            "group_id": [f"Chr{i}_All" for i in range(1, 6)],
            "nrl_mean": [190.0] * 5,
            "nrl_std": nrl_std,
        }
    )
    path = tmp_path / name
    pd.concat([meta, block], ignore_index=True).to_parquet(path)
    return path


_FITTED = [4.0, 5.1, 6.2, 4.8, 5.5]
_FABRICATED = [5.0] * 5


def test_a_fitted_model_can_be_stamped(tmp_path):
    path = _pon(tmp_path, nrl_std=_FITTED)
    assert stamp_release(path, "0.9.0") == "0.8.3"
    meta = pd.read_parquet(path).query("table == 'metadata'").iloc[0]
    assert meta.krewlyzer_version == "0.9.0"


def test_a_fabricated_model_is_refused(tmp_path):
    """The laundering path this guard closes.

    One σ repeated across every group is the `wps_background` signature. A
    stamp would assert 0.9.0 compatibility for a model that has none.
    """
    path = _pon(tmp_path, nrl_std=_FABRICATED)
    with pytest.raises(ValueError, match="fails validate-pon"):
        stamp_release(path, "0.9.0")
    meta = pd.read_parquet(path).query("table == 'metadata'").iloc[0]
    assert meta.krewlyzer_version == "0.8.3", "the file must be left untouched"


def test_the_refusal_names_what_failed(tmp_path):
    """So the fix is obvious without running validate-pon separately."""
    path = _pon(tmp_path, nrl_std=_FABRICATED)
    with pytest.raises(ValueError) as excinfo:
        stamp_release(path, "0.9.0")
    assert "BLOCK_DEGENERATE" in str(excinfo.value)


def test_force_exists_for_restamping_and_says_so(tmp_path):
    """The one legitimate exception: correcting a version on a model that passed.

    Pinned so that if `force` ever becomes the default path, this fails.
    """
    path = _pon(tmp_path, nrl_std=_FABRICATED)
    assert stamp_release(path, "0.9.0", force=True) == "0.8.3"


def test_a_missing_version_alone_does_not_block(tmp_path):
    """`PON.NO_VERSION` is the condition being fixed, not a reason to refuse.

    Models built before provenance existed have no such column at all.
    """
    path = _pon(tmp_path, nrl_std=_FITTED, version="")
    table = pd.read_parquet(path).drop(columns=["krewlyzer_version"])
    table.to_parquet(path, index=False)
    assert stamp_release(path, "0.9.0") == ""
    meta = pd.read_parquet(path).query("table == 'metadata'").iloc[0]
    assert meta.krewlyzer_version == "0.9.0"


def test_nothing_but_the_version_changes(tmp_path):
    """Every baseline copied through, and the build date left alone.

    `build_date` matters: with the version now meaning "the release this ships
    with", the date is the only remaining record of when it was actually built.
    """
    path = _pon(tmp_path, nrl_std=_FITTED)
    before = pd.read_parquet(path)
    stamp_release(path, "0.9.0")
    after = pd.read_parquet(path)

    b = before[before.table != "metadata"].reset_index(drop=True)
    a = after[after.table != "metadata"].reset_index(drop=True)
    pd.testing.assert_frame_equal(a, b)

    meta = after.query("table == 'metadata'").iloc[0]
    assert meta.build_date == "2026-08-03"
    assert meta.cohort_digest == "abc123def4567890"
    assert meta.n_samples == 21.0


def test_a_dry_run_writes_nothing(tmp_path):
    path = _pon(tmp_path, nrl_std=_FITTED)
    assert stamp_release(path, "0.9.0", dry_run=True) == "0.8.3"
    meta = pd.read_parquet(path).query("table == 'metadata'").iloc[0]
    assert meta.krewlyzer_version == "0.8.3"


def test_something_that_is_not_a_pon_is_refused(tmp_path):
    path = tmp_path / "other.parquet"
    pd.DataFrame({"a": [1]}).to_parquet(path)
    with pytest.raises(ValueError):
        stamp_release(path, "0.9.0")
