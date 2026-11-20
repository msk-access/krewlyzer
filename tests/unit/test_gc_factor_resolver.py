"""One resolver for every GC-factor lookup (redesign item 0).

Production runs write `.correction_factors.tsv.gz` and `.parquet` but never a plain `.tsv`
(0.9.2 IGO: 3,158 folders; 0.9.0 clinical: 33,104 folders). Python consumers that looked for a
plain `.tsv` -- panel WPS, every TFBS and ATAC table, gene FSC -- therefore ran without GC
correction while PON donors (built with TSV forced) were corrected. These tests pin the one
function that finds a factor table, so no call site can guess a filename again.
"""

from pathlib import Path

import pytest

from krewlyzer.core.gc_assets import (
    GcFactorsChoice,
    factor_table_base,
    resolve_gc_factors,
)


def _touch(p: Path) -> Path:
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("x")
    return p


def test_finds_the_tables_production_actually_writes(tmp_path):
    # The defect: only .tsv.gz and .parquet exist, never a plain .tsv.
    bed = _touch(tmp_path / "s.bed.gz")
    _touch(tmp_path / "s.correction_factors.tsv.gz")
    _touch(tmp_path / "s.correction_factors.ontarget.tsv.gz")
    off = resolve_gc_factors(bed, tmp_path, "s", "off_target")
    on = resolve_gc_factors(bed, tmp_path, "s", "on_target")
    assert (off.kind, off.path.name) == ("tsv.gz", "s.correction_factors.tsv.gz")
    assert (on.kind, on.path.name) == ("tsv.gz", "s.correction_factors.ontarget.tsv.gz")


def test_precedence_is_parquet_then_tsv_gz_then_tsv(tmp_path):
    bed = _touch(tmp_path / "s.bed.gz")
    for ext in (".tsv", ".tsv.gz", ".parquet"):
        _touch(tmp_path / f"s.correction_factors{ext}")
    for expected in ("parquet", "tsv.gz", "tsv"):
        c = resolve_gc_factors(bed, tmp_path, "s", "off_target")
        assert c.kind == expected
        c.path.unlink()
    assert resolve_gc_factors(bed, tmp_path, "s", "off_target").kind == "none"


def test_beside_the_bed_wins_over_the_output_dir(tmp_path):
    bed = _touch(tmp_path / "extract" / "s.bed.gz")
    beside = _touch(tmp_path / "extract" / "s.correction_factors.tsv.gz")
    _touch(tmp_path / "out" / "s.correction_factors.parquet")
    c = resolve_gc_factors(bed, tmp_path / "out", "s", "off_target")
    assert c.path == beside


def test_output_dir_is_searched_by_sample_name(tmp_path):
    bed = _touch(tmp_path / "in" / "renamed.bed.gz")
    table = _touch(tmp_path / "out" / "s.correction_factors.ontarget.parquet")
    c = resolve_gc_factors(bed, tmp_path / "out", "s", "on_target")
    assert c.path == table and c.kind == "parquet"


def test_without_an_output_dir_only_the_bed_side_is_searched(tmp_path):
    # The table fed *into* the Rust pass: a table an earlier single-tool run left in the output
    # directory must not be silently reused (unchanged behaviour of resolve_gc_assets).
    bed = _touch(tmp_path / "in" / "s.bed.gz")
    _touch(tmp_path / "out" / "s.correction_factors.parquet")
    assert resolve_gc_factors(bed, None, "s", "off_target").kind == "none"


def test_not_found_says_what_was_searched_without_naming_the_sample(tmp_path):
    bed = _touch(tmp_path / "P-0000000.bed.gz")
    c = resolve_gc_factors(bed, tmp_path, "P-0000000", "on_target")
    assert c.path is None and c.kind == "none" and not c.found
    assert (
        "P-0000000.correction_factors.ontarget.{parquet,tsv.gz,tsv}" in c.searched
    )  # for logs
    assert (
        "P-0000000" not in c.reason and "P-0000000" not in c.record()
    )  # metadata stays constant
    assert c.record() == "none:not found beside the BED or in the output directory"


def test_record_is_the_kind_never_the_sample_name(tmp_path):
    # The metadata column must not embed the sample ID: a constant-by-design column is what lets
    # the degeneracy gate notice a cohort where nothing was corrected.
    bed = _touch(tmp_path / "P-0000000.bed.gz")
    _touch(tmp_path / "P-0000000.correction_factors.parquet")
    assert (
        resolve_gc_factors(bed, tmp_path, "P-0000000", "off_target").record()
        == "parquet"
    )


def test_base_path_is_the_one_place_that_spells_the_name(tmp_path):
    assert (
        factor_table_base(tmp_path, "s", "off_target")
        == tmp_path / "s.correction_factors"
    )
    assert (
        factor_table_base(tmp_path, "s", "on_target")
        == tmp_path / "s.correction_factors.ontarget"
    )


def test_unknown_population_is_refused(tmp_path):
    with pytest.raises(ValueError, match="population"):
        factor_table_base(tmp_path, "s", "ontarget")


def test_choice_is_immutable():
    c = GcFactorsChoice(None, "none", "x")
    with pytest.raises(Exception):
        c.kind = "parquet"  # type: ignore[misc]


def test_a_pon_donor_without_gc_tables_is_refused(tmp_path):
    # A silently uncorrected donor in a baseline is worse than a missing donor.
    from krewlyzer.core.gc_assets import require_gc_for_pon

    found = GcFactorsChoice(tmp_path / "d.correction_factors.parquet", "parquet")
    missing = GcFactorsChoice(
        None, "none", "not found beside the BED or in the output directory"
    )
    require_gc_for_pon({"off_target": found, "on_target": found}, panel_mode=True)
    require_gc_for_pon({"off_target": found, "on_target": missing}, panel_mode=False)
    with pytest.raises(RuntimeError, match="off-target"):
        require_gc_for_pon(
            {"off_target": missing, "on_target": found}, panel_mode=False
        )
    with pytest.raises(RuntimeError, match="on-target"):
        require_gc_for_pon({"off_target": found, "on_target": missing}, panel_mode=True)
