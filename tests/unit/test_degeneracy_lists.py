"""List columns with NaN bins must still fingerprint by their values.

The degeneracy gate reduces each array-valued row to one float before comparing
samples. With a plain mean, one NaN bin makes the whole row NaN, and NaN maps to
a single sentinel -- so two samples whose WPS ``prot_frac`` vectors differ
everywhere except in their empty bins would get the same signature, and a
column that varies would be reported as constant across the cohort. Empty bins
are NaN by design (0/0 is not a protection fraction), so every real sample has
some.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from krewlyzer.validate.contract import Kind
from krewlyzer.validate.degeneracy import _nunique, column_signature


def _vectors(scale: float) -> pd.Series:
    rows = []
    for i in range(4):
        v = np.linspace(0.1, 0.9, 200) * scale + 0.01 * i
        v[150:] = np.nan  # empty bins past the last fragment
        rows.append(v.astype(np.float32))
    return pd.Series(rows)


def test_two_samples_differing_outside_their_empty_bins_have_different_signatures():
    a, b = _vectors(1.0), _vectors(0.5)
    assert column_signature(a, Kind.LIST) != column_signature(b, Kind.LIST)


def test_rows_with_empty_bins_are_not_counted_as_one_value():
    assert _nunique(_vectors(1.0), Kind.LIST) == 4


def test_an_all_empty_row_still_reduces_to_nan():
    s = pd.Series([np.full(200, np.nan, dtype=np.float32), np.linspace(0, 1, 200)])
    assert _nunique(s, Kind.LIST) == 2


# --- Columns that are mostly empty by design (`sparse_reason`) -------------
#
# The degeneracy gate fingerprints each sample's first rows. A column that is
# NaN almost everywhere -- WPS `prot_frac_tf` on panel data is NaN in ~99.5% of
# bins, because almost no 35-80 bp fragments reach the anchors -- gives every
# sample the same all-NaN fingerprint, which reads as "identical in every
# sample": a degeneracy error for a column that is behaving exactly as it
# should. `sparse_reason` declares the emptiness, in writing, and the gate then
# judges the column on the samples that have values.

from krewlyzer.validate.contract import CONTRACT, ColumnRule, Vary  # noqa: E402
from krewlyzer.validate.degeneracy import evaluate, observe  # noqa: E402
from krewlyzer.validate.findings import Severity  # noqa: E402


def _empty_rows() -> pd.Series:
    return pd.Series([np.full(200, np.nan, dtype=np.float32) for _ in range(4)])


def _observations(columns):
    return [observe(f"S{i}", s, Kind.LIST) for i, s in enumerate(columns)]


_SPARSE = ColumnRule(
    "prot_frac_tf", Kind.LIST, Vary.BOTH, sparse_reason="rare on panels"
)
_DENSE = ColumnRule("prot_frac_tf", Kind.LIST, Vary.BOTH)


def test_an_ordinary_column_empty_everywhere_is_still_an_error():
    found = evaluate(
        ".WPS.parquet", _DENSE, _observations([_empty_rows(), _empty_rows()]), 2
    )
    assert any(f.severity is Severity.ERROR for f in found), found


def test_a_sparse_column_empty_everywhere_is_a_skip_citing_its_reason():
    found = evaluate(
        ".WPS.parquet", _SPARSE, _observations([_empty_rows(), _empty_rows()]), 2
    )
    assert found and all(f.severity is Severity.SKIP for f in found), found
    assert any("rare on panels" in f.message for f in found), [f.message for f in found]


def test_a_sparse_column_is_judged_on_the_samples_that_have_values():
    # Two samples with distinct values, one empty: varies, so no error.
    found = evaluate(
        ".WPS.parquet",
        _SPARSE,
        _observations([_vectors(1.0), _vectors(0.5), _empty_rows()]),
        2,
    )
    assert not any(f.severity is Severity.ERROR for f in found), found
    # ...and identical values in the samples that have any are still degenerate.
    found = evaluate(
        ".WPS.parquet",
        _SPARSE,
        _observations([_vectors(1.0), _vectors(1.0), _empty_rows()]),
        2,
    )
    assert any(f.severity is Severity.ERROR for f in found), found


def test_the_wps_contract_declares_prot_frac_tf_sparse():
    rules = {
        c.name: c for t in CONTRACT if t.suffix == ".WPS.parquet" for c in t.columns
    }
    assert rules[
        "prot_frac_tf"
    ].sparse_reason, "prot_frac_tf needs a written sparse_reason"
    assert not rules[
        "prot_frac_nuc"
    ].sparse_reason, "prot_frac_nuc is not mostly empty; keep it strict"
