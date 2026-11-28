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
    off = resolve_gc_factors(bed, "off_target")
    on = resolve_gc_factors(bed, "on_target")
    assert (off.kind, off.path.name) == ("tsv.gz", "s.correction_factors.tsv.gz")
    assert (on.kind, on.path.name) == ("tsv.gz", "s.correction_factors.ontarget.tsv.gz")


def test_precedence_is_parquet_then_tsv_gz_then_tsv(tmp_path):
    bed = _touch(tmp_path / "s.bed.gz")
    for ext in (".tsv", ".tsv.gz", ".parquet"):
        _touch(tmp_path / f"s.correction_factors{ext}")
    for expected in ("parquet", "tsv.gz", "tsv"):
        c = resolve_gc_factors(bed, "off_target")
        assert c.kind == expected
        c.path.unlink()
    assert resolve_gc_factors(bed, "off_target").kind == "none"


def _record(directory: Path, output_format: str, compress: bool, **columns) -> None:
    """Extraction's record of `s`, in the format a run of this format writes."""
    import pandas as pd

    from krewlyzer.core.output_utils import write_table

    write_table(
        pd.DataFrame([{"sample_id": "s", **columns}]),
        directory / "s.metadata",
        output_format=output_format,
        compress=compress,
    )


@pytest.mark.parametrize(
    "table_ext, output_format, compress",
    [
        (".parquet", "parquet", False),
        (".tsv.gz", "tsv", True),
        (".tsv", "tsv", False),
        (".tsv", "parquet", False),  # a TSV tool run into a Parquet extraction's folder
    ],
)
def test_a_table_extraction_recorded_it_did_not_write_is_not_used(
    tmp_path, table_ext, output_format, compress
):
    # Builds before the refit's removal fitted a table beside the BED after extraction wrote none
    # (run-all's feature step writes into the BED's folder), under extraction's name. A rerun then
    # took it for extraction's table. Extraction's record says it wrote none, so it is not.
    bed = _touch(tmp_path / "s.bed.gz")
    _touch(tmp_path / f"s.correction_factors{table_ext}")
    _record(tmp_path, output_format, compress, gc_correction_computed=False)
    c = resolve_gc_factors(bed, "off_target")
    assert not c.found
    assert c.record() == (
        "none:beside the BED but not extraction's, whose record says it wrote none"
    )


@pytest.mark.parametrize(
    "columns",
    [{"gc_correction_computed": True}, {}, {"gc_correction_computed": None}],
    ids=["wrote-one", "not-recorded", "empty"],
)
def test_a_table_is_extractions_unless_its_record_says_it_wrote_none(tmp_path, columns):
    # The contrast: a record that does not say "none" leaves the table extraction's. Without
    # this, the test above passes on a lookup that never finds a table beside a record.
    bed = _touch(tmp_path / "s.bed.gz")
    _touch(tmp_path / "s.correction_factors.parquet")
    _record(tmp_path, "parquet", False, **columns)
    assert resolve_gc_factors(bed, "off_target").found


def test_a_table_is_judged_by_the_record_of_its_own_format(tmp_path):
    # Two records that disagree, as a TSV run beside a Parquet extraction can leave: each table is
    # judged by the record in its own format, the run that wrote it, not by whichever reads first.
    bed = _touch(tmp_path / "s.bed.gz")
    _record(tmp_path, "tsv", False, gc_correction_computed=False)
    _record(tmp_path, "parquet", False, gc_correction_computed=True)
    tsv = _touch(tmp_path / "s.correction_factors.tsv")
    assert not resolve_gc_factors(bed, "off_target").found
    tsv.unlink()
    _touch(tmp_path / "s.correction_factors.parquet")
    assert resolve_gc_factors(bed, "off_target").kind == "parquet"


def test_a_record_that_cannot_be_read_is_warned_of_and_says_nothing(tmp_path, caplog):
    # Judged as no record, as before the check: the lookup neither fails the run nor quietly
    # drops the table on a damaged record.
    import logging

    bed = _touch(tmp_path / "s.bed.gz")
    _touch(tmp_path / "s.correction_factors.parquet")
    (tmp_path / "s.metadata.parquet").write_bytes(b"not parquet")
    with caplog.at_level(logging.WARNING, logger="core.output_provenance"):
        assert resolve_gc_factors(bed, "off_target").found
    assert any("could not be read" in r.getMessage() for r in caplog.records)


def test_extractions_record_of_off_target_says_nothing_of_the_on_target_table(tmp_path):
    # `gc_correction_computed` is the off-target table's; no feature step ever fitted an
    # on-target one.
    bed = _touch(tmp_path / "s.bed.gz")
    _touch(tmp_path / "s.correction_factors.ontarget.parquet")
    _record(tmp_path, "parquet", False, gc_correction_computed=False)
    assert resolve_gc_factors(bed, "on_target").found


def test_beside_a_linked_bed_is_beside_the_file_it_points_to(tmp_path):
    # Extraction wrote the table with the BED, and the tools judge a linked BED by its file
    # (`refuse_a_bed_the_record_does_not_describe`). Beside the link, another table was read
    # with that BED's fragments (sixth review).
    bed = _touch(tmp_path / "b" / "s.bed.gz")
    table = _touch(tmp_path / "b" / "s.correction_factors.parquet")
    _touch(tmp_path / "c" / "s.correction_factors.parquet")
    (tmp_path / "c" / "s.bed.gz").symlink_to(bed)
    c = resolve_gc_factors(tmp_path / "c" / "s.bed.gz", "off_target")
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
    c = resolve_gc_factors(tmp_path / "b" / "s.bed.gz", "off_target")
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
    c = resolve_gc_factors(tmp_path / "d" / "s.bed.gz", "off_target")
    assert c.found and c.path.samefile(table)


def test_a_loop_of_links_is_not_followed_forever(tmp_path):
    (tmp_path / "x").mkdir()
    (tmp_path / "y").mkdir()
    (tmp_path / "x" / "s.bed.gz").symlink_to(tmp_path / "y" / "s.bed.gz")
    (tmp_path / "y" / "s.bed.gz").symlink_to(tmp_path / "x" / "s.bed.gz")
    c = resolve_gc_factors(tmp_path / "x" / "s.bed.gz", "off_target")
    assert c.kind == "none"


def test_a_link_to_another_extractions_bed_is_followed_whatever_its_name(tmp_path):
    # A record beside the target says extraction wrote it there: its table came with it.
    target = _touch(tmp_path / "a" / "t.bed.gz")
    _touch(tmp_path / "a" / "t.metadata.parquet")
    table = _touch(tmp_path / "a" / "t.correction_factors.parquet")
    _touch(tmp_path / "b" / "s.correction_factors.parquet")
    (tmp_path / "b" / "s.bed.gz").symlink_to(target)
    c = resolve_gc_factors(tmp_path / "b" / "s.bed.gz", "off_target")
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
    c = resolve_gc_factors(tmp_path / "b" / "s.bed.gz", "off_target")
    assert not (c.found and c.path.samefile(table_b))


@pytest.mark.parametrize(
    "middle",
    ["b/s.v1.bed.gz", "beds/s_cfDNA.bed.gz"],
    ids=["renamed-link-beside", "folder-of-renamed-links"],
)
def test_a_chain_through_a_renamed_link_uses_the_table_of_the_bed_it_reaches(
    tmp_path, middle
):
    # `b/s.bed.gz -> <renamed link> -> a/s.bed.gz`: the reader reaches a's BED, so a's table goes
    # with it, not b's (ninth review).
    bed = _touch(tmp_path / "a" / "s.bed.gz")
    _touch(tmp_path / "a" / "s.metadata.parquet")
    table = _touch(tmp_path / "a" / "s.correction_factors.parquet")
    _touch(tmp_path / "b" / "s.correction_factors.parquet")
    link = tmp_path / middle
    link.parent.mkdir(exist_ok=True)
    link.symlink_to(bed)
    (tmp_path / "b" / "s.bed.gz").symlink_to(link)
    c = resolve_gc_factors(tmp_path / "b" / "s.bed.gz", "off_target")
    assert c.found and c.path.samefile(table)


def test_a_versioned_file_beside_the_link_held_in_an_annex_is_the_links_bed(tmp_path):
    # `b/s.bed.gz -> s.v1.bed.gz -> <annex object>`: no record claims anything after the link.
    obj = _touch(tmp_path / "store" / "SHA256E-s1--0a1b.bed.gz")
    table = _touch(tmp_path / "b" / "s.correction_factors.parquet")
    (tmp_path / "b" / "s.v1.bed.gz").symlink_to(obj)
    (tmp_path / "b" / "s.bed.gz").symlink_to(tmp_path / "b" / "s.v1.bed.gz")
    c = resolve_gc_factors(tmp_path / "b" / "s.bed.gz", "off_target")
    assert c.found and c.path.samefile(table)


def test_not_found_says_what_was_searched_without_naming_the_sample(tmp_path):
    bed = _touch(tmp_path / "P-0000000.bed.gz")
    c = resolve_gc_factors(bed, "on_target")
    assert c.path is None and c.kind == "none" and not c.found
    assert (
        "P-0000000.correction_factors.ontarget.{parquet,tsv.gz,tsv}" in c.searched
    )  # for logs
    assert (
        "P-0000000" not in c.reason and "P-0000000" not in c.record()
    )  # metadata stays constant
    assert c.record() == "none:not found beside the BED"


def test_record_is_the_kind_never_the_sample_name(tmp_path):
    # The metadata column must not embed the sample ID: a constant-by-design column is what lets
    # the degeneracy gate notice a cohort where nothing was corrected.
    bed = _touch(tmp_path / "P-0000000.bed.gz")
    _touch(tmp_path / "P-0000000.correction_factors.parquet")
    assert resolve_gc_factors(bed, "off_target").record() == "parquet"


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
    missing = GcFactorsChoice(None, "none", "not found beside the BED")
    require_gc_for_pon({"off_target": found, "on_target": found}, panel_mode=True)
    require_gc_for_pon({"off_target": found, "on_target": missing}, panel_mode=False)
    with pytest.raises(RuntimeError, match="off-target"):
        require_gc_for_pon(
            {"off_target": missing, "on_target": found}, panel_mode=False
        )
    with pytest.raises(RuntimeError, match="on-target"):
        require_gc_for_pon({"off_target": found, "on_target": missing}, panel_mode=True)
