"""Did the PON actually reach the output, across a real cohort?

The gap this closes: `validate-output` checks that a table is shaped
correctly, and the unit tests check that a scorer *can* run. Neither notices a
scorer that never fired. Measured on the 0.8.3 corpus, four PON blocks were
built, shipped, and produced no column in any of 26 samples.

Every assertion here is over the cohort, because that is the only scale at
which "the code has a scorer" and "the scorer fired" differ.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from .conftest import sample_label

pytestmark = pytest.mark.integration

#: Table suffix -> the column its PON baseline is supposed to produce.
#:
#: Derived from `pon/build.py`: every block it writes should show up somewhere.
#: A table listed here with no such column in any sample means the baseline was
#: computed, stored and shipped for nothing.
EXPECTED_PON_COLUMNS = {
    ".FSR.parquet": ["short_long_log2"],
    ".FSC.parquet": ["ultra_short_log2", "mono_nucl_log2"],
    ".MDS.parquet": ["mds_z"],
    ".MDS.gene.parquet": ["mds_z"],
    ".TFBS.parquet": ["z_score"],
    ".ATAC.parquet": ["z_score"],
    ".correction_factors.parquet": ["expected"],
}


def _present(tables, suffix):
    frames = tables.get(suffix)
    if not frames:
        pytest.skip(f"{suffix} absent from this cohort")
    return frames


@pytest.mark.parametrize(
    "suffix,columns", sorted((s, c) for s, c in EXPECTED_PON_COLUMNS.items())
)
def test_the_pon_column_exists_in_every_sample(corpus_tables, suffix, columns):
    """A scorer that fired in some samples and not others is the harder bug."""
    frames = _present(corpus_tables, suffix)
    for column in columns:
        missing = sum(1 for f in frames if column not in f.columns)
        assert missing == 0, (
            f"{suffix} lacks `{column}` in {missing}/{len(frames)} samples. "
            "The PON was applied to some and not others."
        )


@pytest.mark.parametrize(
    "suffix,columns", sorted((s, c) for s, c in EXPECTED_PON_COLUMNS.items())
)
def test_the_pon_column_is_populated_not_merely_present(corpus_tables, suffix, columns):
    """An all-NaN column is a column that exists and says nothing."""
    frames = _present(corpus_tables, suffix)
    for column in columns:
        empty = [
            i
            for i, f in enumerate(frames)
            if column in f.columns
            and not pd.to_numeric(f[column], errors="coerce").notna().any()
        ]
        assert not empty, (
            f"{suffix}.{column} is entirely NaN in {len(empty)}/{len(frames)} "
            "samples: present, and carrying no information."
        )


@pytest.mark.parametrize(
    "suffix,columns", sorted((s, c) for s, c in EXPECTED_PON_COLUMNS.items())
)
def test_the_pon_column_varies_between_samples(corpus_tables, suffix, columns):
    """Invariant #1, on real data.

    `test_antidegeneracy.py` compares two synthetic profiles. A cohort of real
    samples is the honest test: a metric identical across genuinely different
    people is not measuring them.
    """
    frames = _present(corpus_tables, suffix)
    for column in columns:
        per_sample = [
            float(pd.to_numeric(f[column], errors="coerce").mean(skipna=True))
            for f in frames
            if column in f.columns
        ]
        finite = [v for v in per_sample if np.isfinite(v)]
        if len(finite) < 2:
            pytest.skip(f"{suffix}.{column}: fewer than 2 finite sample means")
        assert len(set(np.round(finite, 12))) > 1, (
            f"{suffix}.{column} has the identical mean in all {len(finite)} "
            "samples. A metric that cannot vary with the data is worse than a "
            "missing one."
        )


def test_no_sample_is_silently_dropped_downstream(corpus):
    """`.metadata.parquet` is the consumer's completion marker (invariant #2).

    A sample without it vanishes from the cohort with no warning, so its
    absence has to be an error somewhere, and this is the only place that sees
    a whole cohort at once.
    """
    missing = [
        sample_label(d)
        for d in corpus
        if not (d / f"{d.name}.metadata.parquet").exists()
    ]
    assert not missing, (
        f"{len(missing)}/{len(corpus)} samples have no .metadata.parquet and "
        f"would be dropped silently: {missing}"
    )
