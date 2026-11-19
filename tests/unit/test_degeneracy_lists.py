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
