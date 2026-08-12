"""A `.1`-suffixed column in a shipped file means a frame collision was written.

pandas renames a duplicate column to `name.1` on read and writes that name back
out. So the suffix means some step appended columns that were already there,
and no reader can tell which copy is current -- `read_csv` returns the first,
which is the oldest.

Found in FSD, where `apply_pon_logratio` matched its own `65-69_logR` output as
size bin 65 and every re-run appended another generation: 69 columns to 137 to
273 on a real sample. The check is deliberately generic, so the next step that
appends instead of replacing is caught without anyone thinking to look.
"""

from __future__ import annotations

import pandas as pd

from krewlyzer.validate.checks import no_collided_columns


def test_a_clean_frame_has_nothing_to_report():
    assert no_collided_columns(pd.DataFrame({"region": ["1p"], "65-69": [1.0]})) == []


def test_the_pandas_suffix_is_caught():
    frame = pd.DataFrame([[1.0, 2.0]], columns=["65-69_logR", "65-69_logR.1"])
    problems = no_collided_columns(frame)
    assert problems and "65-69_logR.1" in problems[0]


def test_outright_duplicate_names_are_caught():
    """A Rust writer emits the repeated name directly, with no suffix at all."""
    frame = pd.DataFrame([[1.0, 2.0]], columns=["x", "x"])
    problems = no_collided_columns(frame)
    assert problems and "duplicate column name" in problems[0]


def test_a_legitimate_dotted_name_is_not_flagged():
    """`p.value` has a dot but no numeric suffix and no partner column."""
    assert no_collided_columns(pd.DataFrame({"p.value": [0.1], "chr1.2": [3]})) == []


def test_the_real_fsd_shape_is_caught():
    """Three generations, exactly as measured on a shipped sample."""
    columns = ["region"]
    for gen in ("", ".1", ".2"):
        columns += [f"{s}-{s + 4}_logR{gen}" for s in range(65, 85, 5)]
    frame = pd.DataFrame([[0] * len(columns)], columns=columns)
    assert no_collided_columns(frame)


def test_the_check_runs_on_every_table():
    """A new output must not be able to opt out by omission."""
    from krewlyzer.validate import gate

    source = open((gate.__file__ or "").replace(".pyc", ".py")).read()
    assert '"no_collided_columns"' in source
