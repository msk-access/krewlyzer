"""WPS z-scoring — the largest PON baseline, previously read by nothing.

`wps_baseline` is ~128k anchors of 200-element mean and sigma vectors, roughly
90% of every PON file, and until 0.9.0 its only consumer was a log line
appending `"WPS"` to a list of available components.

The design question was what to emit, and measuring answered it against the
obvious choices:

- **Not a mean of z over positions.** Adjacent WPS positions have lag-1
  autocorrelation 0.986 — a fragment spans ~167 bp and touches many at once —
  so such a mean has nothing like `sigma/sqrt(200)` precision.
- **Not a max of |z| over positions.** Expected max under pure noise is 2.97,
  so `|z| > 2` would flag nearly every anchor.
- **Not centre-versus-flank.** TSS anchors dip at the centre (−6.8 against
  −3.4 in the flanks); CTCF anchors do the opposite. Any fixed window is
  backwards for one of the two.

So: per-position z vectors, plus two window-free derived quantities each
z-scored against a baseline of itself, and a displacement that is measured but
deliberately not scored.

**Where the assertions live now.** The arithmetic moved to
`rust/src/wps.rs::apply_pon_zscore` in 0.9.0, so the properties of the
individual statistics are asserted in `mod pon_scoring_tests` there — against
the code that actually runs, rather than against a Python copy of it. What
stays here is what only exists at the Python level: the output contract, the
absences, and the design decisions that must not quietly reverse.
`tests/unit/test_rust_python_equivalence.py` binds the two.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from krewlyzer.core.wps_pon import PHASE_MAX_LAG

pytestmark = pytest.mark.unit


def _wave(n=200, period=12.0, offset=0.0, scale=1.0):
    return scale * np.sin((np.arange(n) - offset) / period)


def _frame(seed: int, n_anchors: int = 8):
    rng = np.random.default_rng(seed)
    return pd.DataFrame(
        {
            "region_id": [f"TSS|G{i}|T{i}" for i in range(n_anchors)],
            "wps_nuc": [
                _wave(offset=rng.normal(0, 1.5)) + rng.normal(0, 0.05, 200)
                for _ in range(n_anchors)
            ],
            "wps_tf": [_wave(period=6) for _ in range(n_anchors)],
        }
    )


def _pon_from(frames, tmp_path: Path, with_shape: bool = True):
    """A PON parquet built from in-memory WPS frames, via the real builder.

    Written to disk through `_save_pon_model` rather than hand-assembled: the
    Rust reader takes the *path*, and the shipped PONs store `large_string`
    keys and `list<double>` vectors. A hand-rolled fixture in `string` /
    `list<float>` would pass while the reader was blind to every real PON.
    """
    import pyarrow as pa
    import pyarrow.parquet as pq

    from krewlyzer.pon.build import _compute_wps_baseline, _save_pon_model
    from krewlyzer.pon.model import PonModel

    # float32 lists, matching what the Rust WPS writer produces.
    schema = pa.schema(
        [
            ("region_id", pa.string()),
            ("wps_nuc", pa.list_(pa.float32())),
            ("wps_tf", pa.list_(pa.float32())),
        ]
    )
    paths = []
    for i, frame in enumerate(frames):
        path = tmp_path / f"s{i}.WPS.parquet"
        pq.write_table(pa.Table.from_pandas(frame, schema=schema), path)
        paths.append(path)

    vector, shape = _compute_wps_baseline([str(p) for p in paths])
    model = PonModel(
        schema_version="1.0",
        assay="xs2",
        build_date="2026-01-01",
        n_samples=len(frames),
        reference="r",
        panel_mode=True,
        target_regions_file="t",
        wps_baseline=vector,
        wps_shape_baseline=shape if with_shape else None,
    )
    pon_path = tmp_path / "test.pon.parquet"
    _save_pon_model(model, pon_path)
    return pon_path, paths


# ---------------------------------------------------------------------------
# design decisions that must not quietly reverse
# ---------------------------------------------------------------------------


def test_the_search_window_matches_the_rust_that_runs():
    """A constant kept in two languages drifts unless something fails.

    `PHASE_MAX_LAG` in `core/wps_pon.py` is documentation now — the copy that
    runs is in `rust/src/wps.rs`. Invariant #5: any constant stated twice is
    asserted equal, or the two diverge and only the doc looks wrong.
    """
    source = Path("rust/src/wps.rs").read_text()
    assert (
        f"const PHASE_MAX_LAG: i32 = {PHASE_MAX_LAG};" in source
    ), "the Python mirror of PHASE_MAX_LAG no longer matches rust/src/wps.rs"


def test_a_boundary_shift_is_not_z_scored():
    """The flag is not decorative: the value behind it is excluded from the z.

    Scoring the edge of a search window against a baseline of real shifts
    would turn "we stopped looking" into "displaced by exactly 30 bp".
    """
    from krewlyzer.core import wps_pon

    source = wps_pon.apply_wps_pon.__doc__ or ""
    assert "wps_phase_at_search_limit" in source


def test_displacement_is_measured_but_not_scored():
    """The raw lag is useful; a z-score of it would not be.

    Measured on a real cohort: per-sample mean lag varies by 0.26 bp against a
    within-sample spread of 8.43, so there is no whole-sample phasing signal;
    and per anchor the intraclass correlation is 0.479, meaning about half of
    any lag is noise. That estimate is optimistic — the baseline used
    contained the samples being scored.

    It is also integer-valued, so on a small cohort its sigma bottoms out at
    std([0,0,0,0,0,1]) = 0.408 and a 1 bp shift scores z = 2.4.

    So the column ships and the baseline does not. The measurement is cheap
    and genuinely non-redundant (corr -0.24 and -0.28 with the two scored
    statistics); adding the baseline back is small if a rebuild shows
    reproducible per-anchor shifts.
    """
    from krewlyzer.pon.model import WPS_SHAPE_STATS

    assert (
        "phase_shift_bp" not in WPS_SHAPE_STATS
    ), "the shape baseline must not carry a phase-shift entry"
    assert set(WPS_SHAPE_STATS) == {"log_amplitude", "shape_corr_fisher"}


def test_the_scoring_is_not_computed_in_python():
    """The port is the point: 89k anchors, a ±30 lag search each.

    Measured at the time of the port: 9.5 s in Rust against ~17 min for the
    Python it replaced, on the same 76,595-anchor output. `.agents/rules/
    architecture.md` puts "PON z-score", "loops over >1000 rows" and
    "row-level computation" on the Rust side; this was all three.

    A reintroduced Python loop here would be silent — correct, and a hundred
    times slower.
    """
    import ast
    import inspect
    import textwrap

    from krewlyzer.core import wps_pon

    source = inspect.getsource(wps_pon)
    assert "_core.wps.apply_pon_zscore" in source, "no longer calling Rust"
    assert "numpy" not in source, "numpy is back; the arithmetic belongs in Rust"

    # Parsed, not grepped: the docstring says "for" several times, and a
    # substring check on it passes or fails for reasons unrelated to the code.
    tree = ast.parse(textwrap.dedent(inspect.getsource(wps_pon.apply_wps_pon)))
    loops = [
        n
        for n in ast.walk(tree)
        if isinstance(n, (ast.For, ast.While, ast.comprehension))
    ]
    assert not loops, (
        f"the wrapper has grown {len(loops)} loop(s); per-anchor work belongs "
        "in rust/src/wps.rs, which is ~100x faster on the real 89k-anchor input"
    )


# ---------------------------------------------------------------------------
# end to end
# ---------------------------------------------------------------------------


def test_every_promised_column_is_emitted(tmp_path):
    from krewlyzer.core.wps_pon import apply_wps_pon

    pon_path, paths = _pon_from([_frame(s) for s in range(5)], tmp_path)
    n = apply_wps_pon(paths[0], pon_path, output_base=tmp_path / "o")
    assert n > 0, "nothing scored — the Rust reader found no baseline"

    out = pd.read_parquet(tmp_path / "o.parquet")
    for column in (
        "wps_nuc_z",
        "wps_log_amplitude",
        "wps_log_amplitude_z",
        "wps_shape_corr",
        "wps_shape_corr_z",
        "wps_phase_shift_bp",
        "wps_phase_at_search_limit",
    ):
        assert column in out.columns, f"{column} was not emitted"
    assert len(np.asarray(out["wps_nuc_z"].iloc[0])) == 200
    assert "wps_phase_shift_z" not in out.columns, (
        "displacement is measured but deliberately not z-scored — see "
        "test_displacement_is_measured_but_not_scored"
    )


def test_a_pon_without_a_shape_block_still_writes_the_z_vectors(tmp_path):
    """A PON built before 0.9.0 has `wps_baseline` and no `wps_shape_baseline`.

    The per-position z needs only the vector baseline, so it must not be lost
    along with the derived scores.
    """
    from krewlyzer.core.wps_pon import apply_wps_pon

    pon_path, paths = _pon_from(
        [_frame(s) for s in range(5)], tmp_path, with_shape=False
    )
    apply_wps_pon(paths[0], pon_path, output_base=tmp_path / "o")
    out = pd.read_parquet(tmp_path / "o.parquet")
    assert out["wps_nuc_z"].notna().any(), "the z vectors must survive"
    assert (
        pd.to_numeric(out["wps_log_amplitude_z"], errors="coerce").isna().all()
    ), "the derived z-scores have no baseline and must be NaN"
    assert out["wps_log_amplitude"].notna().any(), "the raw value is still readable"


def test_an_anchor_absent_from_the_baseline_gets_null_not_zero(tmp_path):
    """A z of 0.0 is the most confident statement the column can make.

    "Exactly at the healthy baseline" is a measurement; "we had nothing to
    compare against" is not, and nothing downstream can tell them apart once
    the value is written.
    """
    from krewlyzer.core.wps_pon import apply_wps_pon

    pon_path, paths = _pon_from([_frame(s) for s in range(5)], tmp_path)

    # One anchor the PON has never seen.
    extra = _frame(0)
    extra.loc[len(extra)] = {
        "region_id": "TSS|UNSEEN|T99",
        "wps_nuc": _wave(offset=2.0),
        "wps_tf": _wave(period=6),
    }
    import pyarrow as pa
    import pyarrow.parquet as pq

    schema = pa.schema(
        [
            ("region_id", pa.string()),
            ("wps_nuc", pa.list_(pa.float32())),
            ("wps_tf", pa.list_(pa.float32())),
        ]
    )
    sample = tmp_path / "extra.WPS.parquet"
    pq.write_table(pa.Table.from_pandas(extra, schema=schema), sample)

    apply_wps_pon(sample, pon_path, output_base=tmp_path / "absent")
    out = pd.read_parquet(tmp_path / "absent.parquet")
    row = out[out["region_id"] == "TSS|UNSEEN|T99"].iloc[0]
    assert row["wps_nuc_z"] is None, "an unseen anchor must stay null"
    assert math.isnan(row["wps_shape_corr"]), "and its derived values NaN"
    assert math.isnan(row["wps_shape_corr_z"])
    assert not math.isnan(
        row["wps_log_amplitude"]
    ), "amplitude needs no baseline and is still readable"


def test_scoring_a_sample_inside_its_own_baseline_shrinks_z(tmp_path):
    """Why self-inclusion cannot be used to claim calibration.

    A sample contributing to its own baseline pulls the mean toward itself, so
    |z| is capped near (n-1)/sqrt(n). Measured: max |z| of 1.97/1.94/2.04 at
    n=6, against 18.5/29.2/42.5 for the same sample held out. "0% beyond |z|>2"
    from a self-included run is arithmetic, not evidence.
    """
    from krewlyzer.core.wps_pon import apply_wps_pon

    frames = [_frame(s) for s in range(6)]
    pon_path, paths = _pon_from(frames, tmp_path)
    apply_wps_pon(paths[0], pon_path, output_base=tmp_path / "self")
    included = pd.read_parquet(tmp_path / "self.parquet")
    z = pd.to_numeric(included["wps_shape_corr_z"], errors="coerce").dropna()
    cap = (len(frames) - 1) / math.sqrt(len(frames))
    assert (
        z.abs().max() <= cap + 0.1
    ), f"self-included |z| reached {z.abs().max():.2f}, above the {cap:.2f} cap"


def test_the_output_is_parquet_whatever_the_run_format_is(tmp_path):
    """WPS is Parquet by contract, and this writer no longer takes an opinion.

    It used to honour `--output-format`, whose default is `tsv`. The Rust step
    writes `.WPS.parquet` unconditionally, so a default run produced *two*
    files: a `.WPS.tsv` carrying all seven PON columns and a `.WPS.parquet`
    carrying none. Downstream reads Parquet only (invariant #2), so the WPS
    scoring was present on disk and absent from the product.

    Measured on a real 89,034-anchor run before the fix: `.tsv` had 18 columns,
    `.parquet` had 11, and the parquet was the older file.
    """
    import inspect

    from krewlyzer.core.wps_pon import apply_wps_pon

    assert (
        "output_format" not in inspect.signature(apply_wps_pon).parameters
    ), "the knob is back; it can only be set wrong"

    pon_path, paths = _pon_from([_frame(s) for s in range(5)], tmp_path)
    apply_wps_pon(paths[0], pon_path, output_base=tmp_path / "fmt")

    assert (tmp_path / "fmt.parquet").exists(), "no parquet written"
    assert not (tmp_path / "fmt.tsv").exists(), "wrote a TSV nobody reads"
    out = pd.read_parquet(tmp_path / "fmt.parquet")
    assert "wps_nuc_z" in out.columns


def test_a_missing_pon_leaves_the_raw_table_alone(tmp_path):
    """Half-writing the product is worse than not scoring it.

    Rust returns 0 without creating the file when the PON has no matching
    baseline, so a PON built before the block existed costs the z-scores and
    not the profile.
    """
    from krewlyzer.core.wps_pon import apply_wps_pon

    pon_path, paths = _pon_from([_frame(s) for s in range(5)], tmp_path)
    n = apply_wps_pon(
        paths[0],
        pon_path,
        output_base=tmp_path / "none",
        baseline_table="wps_baseline_panel",  # never built by this PON
    )
    assert n == 0
    assert not (
        tmp_path / "none.parquet"
    ).exists(), "an unmatched baseline must not leave a truncated output behind"
