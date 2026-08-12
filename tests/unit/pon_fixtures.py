"""Minimal-but-complete PON blocks, shared by the gate and stamp tests.

`validate-pon` now checks a packing list: a block absent from the file is a
finding in its own right, because every other check reads the blocks that *are*
there and skips the rest. That was the hole `region_mds` fell through.

Which means a fixture holding only ``metadata`` plus the one block under test
is no longer a valid PON, and its test would fail for a reason that has nothing
to do with what it is testing. `core_blocks` supplies the rest.

Every σ here varies across entries on purpose. A column with one repeated value
is what a fabricated baseline looks like from the outside, and the gate says so
-- correctly. Building fixtures that trip your own degeneracy check wastes an
afternoon; this file exists so nobody repeats it.
"""

from __future__ import annotations

from typing import List

import numpy as np
import pandas as pd

from krewlyzer.validate.pon_gate import CORE_BLOCKS


def _spread(centre: float, n: int) -> List[float]:
    """``n`` distinct values around ``centre``, so no σ column is constant."""
    return [round(centre * (1.0 + 0.11 * i), 6) for i in range(n)]


def core_blocks(exclude: tuple = ()) -> List[pd.DataFrame]:
    """One valid frame per core block, minus anything named in ``exclude``.

    ``exclude`` is how a test says "this block is deliberately missing" without
    hand-building the other seven.
    """
    n = 5
    frames = {
        "gc_bias": pd.DataFrame(
            {
                "table": "gc_bias",
                "gc_bin": [0.30, 0.35, 0.40, 0.45, 0.50],
                "short_expected": _spread(1.0, n),
                "short_std": _spread(0.05, n),
                "intermediate_expected": _spread(1.0, n),
                "intermediate_std": _spread(0.06, n),
                "long_expected": _spread(1.0, n),
                "long_std": _spread(0.07, n),
            }
        ),
        "fsd_baseline": pd.DataFrame(
            {
                "table": "fsd_baseline",
                "arm": ["1p"] * n,
                "size_bin": [65, 70, 75, 80, 85],
                "expected": _spread(0.02, n),
                "std": _spread(0.004, n),
            }
        ),
        "wps_baseline": pd.DataFrame(
            {
                "table": "wps_baseline",
                "region_id": [f"TSS|G{i}|1" for i in range(n)],
                "n_samples": [21] * n,
                "wps_nuc_mean": _spread(0.4, n),
                "wps_nuc_std": _spread(0.09, n),
                "wps_tf_mean": _spread(0.2, n),
                "wps_tf_std": _spread(0.05, n),
            }
        ),
        "wps_shape_baseline": pd.DataFrame(
            {
                "table": "wps_shape_baseline",
                "region_id": [f"TSS|G{i}|1" for i in range(n)],
                "n_samples": [21] * n,
                "log_amplitude_mean": _spread(1.2, n),
                "log_amplitude_std": _spread(0.15, n),
                "shape_corr_fisher_mean": _spread(0.8, n),
                "shape_corr_fisher_std": _spread(0.12, n),
            }
        ),
        "ocf_baseline": pd.DataFrame(
            {
                "table": "ocf_baseline",
                "region_id": ["Lymph", "Liver", "Colon", "Lung", "Breast"],
                "ocf_mean": _spread(10.0, n),
                "ocf_std": _spread(1.5, n),
            }
        ),
        "tfbs_baseline": pd.DataFrame(
            {
                "table": "tfbs_baseline",
                "label": [f"TF{i}" for i in range(n)],
                "entropy_mean": _spread(0.7, n),
                "entropy_std": _spread(0.03, n),
                "n_samples": [21] * n,
            }
        ),
        "atac_baseline": pd.DataFrame(
            {
                "table": "atac_baseline",
                "label": [f"AT{i}" for i in range(n)],
                "entropy_mean": _spread(0.6, n),
                "entropy_std": _spread(0.04, n),
                "n_samples": [21] * n,
            }
        ),
    }
    # `wps_background` is the block most tests are actually about, so it is
    # never supplied here -- the caller builds the version it needs.
    missing = set(CORE_BLOCKS) - set(frames) - {"wps_background"}
    assert not missing, f"core_blocks is short of {sorted(missing)}"

    return [frame for name, frame in frames.items() if name not in exclude]


def background(nrl_std, groups: int = 5) -> pd.DataFrame:
    """A `wps_background` block whose σ is whatever the test needs it to be."""
    ids = [f"Chr{i}_All" for i in range(1, groups + 1)]
    stds = nrl_std if isinstance(nrl_std, list) else [nrl_std] * groups
    return pd.DataFrame(
        {
            "table": "wps_background",
            "group_id": ids,
            "n_samples": [21] * groups,
            "n_at_band_limit": [0] * groups,
            "n_nrl_fitted": [21] * groups,
            "nrl_mean": _spread(190.0, groups),
            "nrl_std": stds,
            "periodicity_mean": _spread(0.44, groups),
            "periodicity_std": _spread(0.02, groups),
        }
    )


def metadata(**overrides) -> pd.DataFrame:
    row = {
        "table": "metadata",
        "schema_version": "1.0",
        "assay": "xs2",
        "n_samples": 21.0,
        "krewlyzer_version": "0.9.0",
        "cohort_digest": "deadbeefcafe0000",
        "cohort_label": "healthy-donors",
        # Fragment BEDs, so the BAM-only blocks are legitimately absent and
        # the gate warns rather than fails. A fixture claiming a BAM cohort
        # would have to carry region_mds and mds_baseline too.
        "input_kind": "bed",
        "panel_mode": False,
    }
    row.update(overrides)
    return pd.DataFrame([row])


def complete(tmp_path, *extra: pd.DataFrame, exclude: tuple = (), **meta):
    """Write a PON that passes the packing-list check, plus ``extra`` blocks."""
    frames = [metadata(**meta)] + core_blocks(exclude=exclude) + list(extra)
    path = tmp_path / "t.pon.parquet"
    pd.concat(frames, ignore_index=True).to_parquet(path)
    return path


__all__ = ["core_blocks", "background", "metadata", "complete", "np"]
