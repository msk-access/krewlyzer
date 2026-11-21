"""A PON is identified by what it holds, not by its name (Greptile on #98, finding 3).

The sample record named a PON by file basename and cohort digest. The cohort digest hashes the
donors' names, not the fitted values, so a PON rebuilt in place from the same donors kept both and
passed as the same PON: run-all reused motif scores made under the old baselines, and the two-PON
check let a single tool score against the rebuilt one. `pon_digest` hashes the baseline rows.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from krewlyzer.core.output_provenance import (
    record_names_pon,
    refuse_a_second_pon,
    stamp_metadata,
)
from krewlyzer.core.pon_integration import UnusablePonError, load_pon_model

from . import pon_fixtures

pytestmark = pytest.mark.unit

_FORMAT = dict(output_format="parquet", compress=False)


def _pon(directory, nrl_std):
    directory.mkdir(parents=True, exist_ok=True)
    return pon_fixtures.complete(
        directory, pon_fixtures.background(nrl_std), exclude=("wps_background",)
    )


def _scored_against(tmp_path, pon_path):
    out = tmp_path / "out"
    out.mkdir(exist_ok=True)
    pd.DataFrame([{"sample_id": "s"}]).to_parquet(out / "s.metadata.parquet")
    stamp_metadata(out, "s", pon=load_pon_model(pon_path), pon_path=pon_path, **_FORMAT)
    return out


def test_a_pon_rebuilt_from_the_same_donors_is_not_the_same_pon(tmp_path):
    """Same file name, same cohort digest, different fitted values."""
    old = _pon(tmp_path / "old", [4.0, 5.1, 6.2, 4.8, 5.5])
    rebuilt = _pon(tmp_path / "rebuilt", [3.0, 4.1, 7.2, 5.8, 6.5])
    assert old.name == rebuilt.name
    out = _scored_against(tmp_path, old)
    new_model = load_pon_model(rebuilt)

    assert not record_names_pon(out, "s", new_model, rebuilt, **_FORMAT)
    with pytest.raises(UnusablePonError):
        refuse_a_second_pon(out, "s", new_model, rebuilt, **_FORMAT)


def test_the_same_pon_is_still_the_same_pon(tmp_path):
    pon = _pon(tmp_path / "p", [4.0, 5.1, 6.2, 4.8, 5.5])
    out = _scored_against(tmp_path, pon)
    model = load_pon_model(pon)
    assert record_names_pon(out, "s", model, pon, **_FORMAT)
    refuse_a_second_pon(out, "s", model, pon, **_FORMAT)


def test_restamping_a_pon_keeps_its_digest(tmp_path):
    """`stamp-pon` rewrites the file at every release; a hash of the file would refuse every
    rerun into a directory scored before the stamp."""
    from krewlyzer.validate.pon_stamp import stamp_release

    pon = _pon(tmp_path / "p", [4.0, 5.1, 6.2, 4.8, 5.5])
    before = load_pon_model(pon).pon_digest
    stamp_release(pon, "0.9.9")
    assert load_pon_model(pon).krewlyzer_version == "0.9.9"
    assert load_pon_model(pon).pon_digest == before
    assert len(before) == 64


def test_a_column_empty_in_every_baseline_row_does_not_change_the_digest(tmp_path):
    """`stamp-pon` backfills a missing `krewlyzer_version` column with "" on every row."""
    pon = _pon(tmp_path / "p", [4.0, 5.1, 6.2, 4.8, 5.5])
    before = load_pon_model(pon).pon_digest
    frame = pd.read_parquet(pon)
    frame["unused_column"] = ""
    frame.to_parquet(pon)
    assert load_pon_model(pon).pon_digest == before


def test_a_record_without_a_digest_is_never_reused_by_the_cache(tmp_path):
    """Written before `pon_digest` existed: nothing in it can show the PON is unchanged."""
    pon = _pon(tmp_path / "p", [4.0, 5.1, 6.2, 4.8, 5.5])
    out = _scored_against(tmp_path, pon)
    record = out / "s.metadata.parquet"
    pd.read_parquet(record).drop(columns=["pon_digest"]).to_parquet(record)
    assert not record_names_pon(out, "s", load_pon_model(pon), pon, **_FORMAT)


def test_the_two_pon_check_falls_back_to_the_name_without_a_digest(tmp_path):
    """A single tool into a directory an older run-all scored still runs against the same PON,
    and is still refused another one."""
    pon = _pon(tmp_path / "p", [4.0, 5.1, 6.2, 4.8, 5.5])
    other = _pon(tmp_path / "q", [4.0, 5.1, 6.2, 4.8, 5.5]).rename(
        tmp_path / "q" / "other.pon.parquet"
    )
    out = _scored_against(tmp_path, pon)
    record = out / "s.metadata.parquet"
    pd.read_parquet(record).drop(columns=["pon_digest"]).to_parquet(record)

    refuse_a_second_pon(out, "s", load_pon_model(pon), pon, **_FORMAT)
    with pytest.raises(UnusablePonError):
        refuse_a_second_pon(out, "s", load_pon_model(other), other, **_FORMAT)


def test_the_digest_of_the_fixture_pon_is_pinned():
    """The serialisation is the digest's definition. A change to it (column order, float
    format, NaN marker) changes every digest and makes every recorded PON look foreign.
    """
    from pathlib import Path

    fixture = Path(__file__).parents[1] / "data" / "fixtures" / "test.pon.parquet"
    from krewlyzer.pon.model import PonModel

    assert PonModel.load(fixture).pon_digest == GOLDEN


GOLDEN = "37b6d9eb4e64260ae54fa7ab42bfea914c4bf9a83c113dfa15220582f5ce41d3"


@pytest.mark.parametrize(
    "cell",
    [
        lambda x: np.array([1.0, 2.0, x]),  # a WPS profile
        lambda x: {"AAAA": x, "CCCC": 0.25},  # a k-mer baseline
    ],
    ids=["array", "dict"],
)
def test_a_change_past_the_eighth_digit_changes_the_digest(cell):
    """numpy's `str` keeps eight significant digits: a digest built on it would call two PONs
    whose fitted profiles differ in the ninth the same PON."""
    from krewlyzer.pon.provenance import baseline_digest

    def frame(x):
        return pd.DataFrame(
            {"table": ["metadata", "wps_baseline"], "profile": [None, cell(x)]}
        )

    assert baseline_digest(frame(0.1234567891)) != baseline_digest(frame(0.1234567894))
