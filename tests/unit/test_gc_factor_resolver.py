"""One resolver for every GC-factor lookup (redesign item 0).

Production runs write `.correction_factors.tsv.gz` and `.parquet` but never a plain `.tsv`
(0.9.2 IGO: 3,158 folders; 0.9.0 clinical: 33,104 folders). Python consumers that looked for a
plain `.tsv` -- panel WPS, every TFBS and ATAC table, gene FSC -- therefore ran without GC
correction while PON donors (built with TSV forced) were corrected. These tests pin the one
function that finds a factor table, so no call site can guess a filename again.
"""

from pathlib import Path

import pytest

from krewlyzer.assets import AssetManager
from krewlyzer.core.gc_assets import (
    GcFactorsChoice,
    factor_table_base,
    resolve_gc_assets,
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


def test_beside_a_linked_bed_is_beside_the_file_it_points_to(tmp_path):
    # Extraction wrote the table with the BED, and the tools judge a linked BED by its file
    # (`refuse_a_bed_the_record_does_not_describe`). Beside the link, another table was read
    # with that BED's fragments (sixth review).
    bed = _touch(tmp_path / "b" / "s.bed.gz")
    table = _touch(tmp_path / "b" / "s.correction_factors.parquet")
    _touch(tmp_path / "c" / "s.correction_factors.parquet")
    (tmp_path / "c" / "s.bed.gz").symlink_to(bed)
    c = resolve_gc_factors(tmp_path / "c" / "s.bed.gz", None, "s", "off_target")
    assert c.path.samefile(table)


@pytest.mark.parametrize(
    "held",
    ["store/SHA256E-s1--0a1b.bed.gz", "store/0a1b2c", "b/s.v1.bed.gz"],
    ids=["annex-object", "dvc-object", "versioned"],
)
def test_beside_a_link_to_a_renamed_file_without_a_record_is_beside_the_link(
    tmp_path, held
):
    # The BED held as a git-annex, DataLad or DVC object, or as a versioned file: following the
    # link looked for the object's table, found none, and the run refitted and deleted the table
    # beside the link (seventh review).
    obj = _touch(tmp_path / held)
    table = _touch(tmp_path / "b" / "s.correction_factors.parquet")
    (tmp_path / "b" / "s.bed.gz").symlink_to(obj)
    c = resolve_gc_factors(tmp_path / "b" / "s.bed.gz", None, "s", "off_target")
    assert c.found and c.path.samefile(table)


def test_a_chain_of_links_is_followed_to_the_last_file_extraction_named(tmp_path):
    # `d/s.bed.gz -> c/s.bed.gz -> b/s.bed.gz -> <annex object>`: b's BED is held in an annex, c
    # links it, and a Nextflow task (d) stages c's link. b is where it was extracted.
    obj = _touch(tmp_path / "store" / "SHA256E-s1--0a1b.bed.gz")
    table = _touch(tmp_path / "b" / "s.correction_factors.parquet")
    (tmp_path / "b" / "s.bed.gz").symlink_to(obj)
    for link, target in (("c", "b"), ("d", "c")):
        _touch(tmp_path / link / "s.correction_factors.parquet")
        (tmp_path / link / "s.bed.gz").symlink_to(tmp_path / target / "s.bed.gz")
    c = resolve_gc_factors(tmp_path / "d" / "s.bed.gz", None, "s", "off_target")
    assert c.found and c.path.samefile(table)


def test_a_loop_of_links_is_not_followed_forever(tmp_path):
    (tmp_path / "x").mkdir()
    (tmp_path / "y").mkdir()
    (tmp_path / "x" / "s.bed.gz").symlink_to(tmp_path / "y" / "s.bed.gz")
    (tmp_path / "y" / "s.bed.gz").symlink_to(tmp_path / "x" / "s.bed.gz")
    c = resolve_gc_factors(tmp_path / "x" / "s.bed.gz", None, "s", "off_target")
    assert c.kind == "none"


def test_a_link_to_another_extractions_bed_is_followed_whatever_its_name(tmp_path):
    # A record beside the target says extraction wrote it there: its table came with it.
    target = _touch(tmp_path / "a" / "t.bed.gz")
    _touch(tmp_path / "a" / "t.metadata.parquet")
    table = _touch(tmp_path / "a" / "t.correction_factors.parquet")
    _touch(tmp_path / "b" / "s.correction_factors.parquet")
    (tmp_path / "b" / "s.bed.gz").symlink_to(target)
    c = resolve_gc_factors(tmp_path / "b" / "s.bed.gz", None, "s", "off_target")
    assert c.found and c.path.samefile(table)


@pytest.mark.parametrize(
    "record",
    ["s.metadata.parquet", "t.metadata.tsv.gz"],
    ids=["its-own", "another-samples"],
)
def test_a_renamed_file_in_a_folder_extraction_wrote_is_not_beside_the_link(
    tmp_path, record
):
    # `b/s.bed.gz -> a/s.v1.bed.gz`, a record in a: the file is a's, not b's BED held under
    # another name, and b's table was read with a's fragments (review of the seventh fix).
    target = _touch(tmp_path / "a" / "s.v1.bed.gz")
    _touch(tmp_path / "a" / record)
    table_b = _touch(tmp_path / "b" / "s.correction_factors.parquet")
    (tmp_path / "b" / "s.bed.gz").symlink_to(target)
    c = resolve_gc_factors(tmp_path / "b" / "s.bed.gz", None, "s", "off_target")
    assert not (c.found and c.path.samefile(table_b))


def test_the_tools_use_the_table_beside_a_link_to_a_store_object_and_do_not_refit(
    tmp_path,
):
    # A refit (`factors_out` set) is the one path to `run_features` removing extract's table from
    # the output folder. The CLI tests cannot see it: the synthetic sample fits no table.
    obj = _touch(tmp_path / "store" / "SHA256E-s1--0a1b.bed.gz")
    table = _touch(tmp_path / "b" / "s.correction_factors.parquet")
    (tmp_path / "b" / "s.bed.gz").symlink_to(obj)
    gc = resolve_gc_assets(
        AssetManager("hg19"), tmp_path / "b", "s", tmp_path / "b" / "s.bed.gz"
    )
    assert gc.factors_input is not None and gc.factors_input.samefile(table)
    assert gc.factors_out is None


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
