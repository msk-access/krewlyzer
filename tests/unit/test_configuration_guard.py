"""A second configuration may not write into a folder the first one owns (Greptile on #98, F2).

A rerun with different settings writes a different set of tables: a panel folder rerun without
its assay rewrites `TFBS` and not `TFBS.ontarget`, and the old on-target table then reads as this
run's under one `tfbs=ran`. Refusing the second configuration before anything is written means no
such leftover can arise. These pin the comparison itself; the run-level sequences are in
`tests/invariants/test_configuration_guard_runs.py`.
"""

from __future__ import annotations

import pandas as pd
import pytest

from krewlyzer.core.output_provenance import (
    ConfigurationMismatch,
    refuse_a_different_configuration,
)

pytestmark = pytest.mark.unit

_FORMAT = dict(output_format="parquet", compress=False)


def _record(tmp_path, **fields):
    row = {
        "sample_id": "s",
        "genome": "hg19",
        "assay": "xs1",
        "target_regions": "/somewhere/staged/xs1.targets.bed",
    }
    row.update(fields)
    pd.DataFrame([row]).to_parquet(tmp_path / "s.metadata.parquet")
    return tmp_path


def test_a_different_assay_is_refused_and_named(tmp_path):
    out = _record(tmp_path)
    with pytest.raises(ConfigurationMismatch, match="assay"):
        refuse_a_different_configuration(out, "s", {"assay": None}, **_FORMAT)


def test_different_targets_are_refused(tmp_path):
    out = _record(tmp_path)
    with pytest.raises(ConfigurationMismatch, match="target_regions"):
        refuse_a_different_configuration(out, "s", {"target_regions": None}, **_FORMAT)


def test_the_same_settings_are_allowed(tmp_path):
    out = _record(tmp_path)
    refuse_a_different_configuration(
        out,
        "s",
        {"genome": "hg19", "assay": "xs1", "target_regions": "xs1.targets.bed"},
        **_FORMAT,
    )


def test_one_build_under_two_names_is_the_same_genome(tmp_path):
    """Nextflow defaults to hg19; a hand repair may say GRCh37."""
    out = _record(tmp_path)
    refuse_a_different_configuration(out, "s", {"genome": "GRCh37"}, **_FORMAT)
    with pytest.raises(ConfigurationMismatch, match="genome"):
        refuse_a_different_configuration(out, "s", {"genome": "hg38"}, **_FORMAT)


def test_targets_compare_by_file_name_not_path(tmp_path):
    """Nextflow stages the target BED at a new path in every task."""
    out = _record(tmp_path)
    refuse_a_different_configuration(
        out, "s", {"target_regions": "/work/ab/12/xs1.targets.bed"}, **_FORMAT
    )


def test_an_assay_compares_without_case(tmp_path):
    out = _record(tmp_path)
    refuse_a_different_configuration(out, "s", {"assay": "XS1"}, **_FORMAT)


def test_e1_aggregation_is_compared_when_given(tmp_path):
    out = _record(tmp_path, fsc_e1_aggregation=False)
    refuse_a_different_configuration(out, "s", {"fsc_e1_aggregation": False}, **_FORMAT)
    with pytest.raises(ConfigurationMismatch, match="fsc_e1_aggregation"):
        refuse_a_different_configuration(
            out, "s", {"fsc_e1_aggregation": True}, **_FORMAT
        )


def test_a_tsv_record_reads_e1_aggregation_as_a_bool(tmp_path):
    """`bool("False")` is True: read naively, a TSV record's False would refuse its own setting."""
    pd.DataFrame(
        [{"sample_id": "s", "genome": "hg19", "fsc_e1_aggregation": False}]
    ).to_csv(tmp_path / "s.metadata.tsv", sep="\t", index=False)
    refuse_a_different_configuration(
        tmp_path,
        "s",
        {"fsc_e1_aggregation": False},
        output_format="tsv",
        compress=False,
    )


def test_a_field_the_record_does_not_hold_is_not_compared(tmp_path):
    """A record written before `fsc_e1_aggregation` existed."""
    out = _record(tmp_path)
    refuse_a_different_configuration(out, "s", {"fsc_e1_aggregation": True}, **_FORMAT)


def test_a_fresh_folder_has_nothing_to_compare(tmp_path):
    refuse_a_different_configuration(tmp_path, "s", {"assay": "xs2"}, **_FORMAT)


def test_an_unknown_setting_is_a_programming_error(tmp_path):
    out = _record(tmp_path)
    with pytest.raises(ValueError, match="unknown setting"):
        refuse_a_different_configuration(out, "s", {"mapq": 20}, **_FORMAT)


# The output format (the plan "refuse a rerun in a different output format"). A rerun in another
# format rewrites its tables in that format and leaves the first run's beside them, and each
# reader then picks by its own rule: downstream reads the old Parquet, `resolve_table_path` the
# old TSV. `table_formats` is not a `given` setting: the guard takes it from the run's own
# `output_format`/`compress`, so no caller can leave it out.


def _records(tmp_path, *extensions, **fields):
    """`s.metadata` in each extension, with the same row."""
    row = {"sample_id": "s", "genome": "hg19", **fields}
    for ext in extensions:
        path = tmp_path / f"s.metadata{ext}"
        if ext == ".parquet":
            pd.DataFrame([row]).to_parquet(path)
        else:
            pd.DataFrame([row]).to_csv(path, sep="\t", index=False)
    return tmp_path


def test_a_different_table_format_is_refused_and_names_the_flag(tmp_path):
    out = _records(tmp_path, ".parquet", table_formats=".parquet")
    with pytest.raises(ConfigurationMismatch, match="table_formats") as refused:
        refuse_a_different_configuration(
            out, "s", {}, output_format="tsv", compress=False
        )
    assert "--output-format parquet" in str(refused.value)


def test_compression_does_not_change_a_parquet_folder(tmp_path):
    out = _records(tmp_path, ".parquet", table_formats=".parquet")
    refuse_a_different_configuration(
        out, "s", {}, output_format="parquet", compress=True
    )


def test_compressed_tsv_is_another_format(tmp_path):
    out = _records(tmp_path, ".tsv", table_formats=".tsv")
    with pytest.raises(ConfigurationMismatch, match="table_formats"):
        refuse_a_different_configuration(
            out, "s", {}, output_format="tsv", compress=True
        )


def test_without_the_field_the_record_files_say_the_format(tmp_path):
    """A record written before `table_formats` existed: every 0.9.0 folder."""
    out = _records(tmp_path, ".parquet")
    refuse_a_different_configuration(
        out, "s", {}, output_format="parquet", compress=False
    )
    with pytest.raises(ConfigurationMismatch, match="table_formats"):
        refuse_a_different_configuration(
            out, "s", {}, output_format="tsv", compress=False
        )


def test_without_the_field_two_records_are_both(tmp_path):
    """`both` writes the record in both formats, so two records are not by themselves mixed."""
    out = _records(tmp_path, ".parquet", ".tsv")
    refuse_a_different_configuration(out, "s", {}, output_format="both", compress=False)
    with pytest.raises(ConfigurationMismatch, match="table_formats"):
        refuse_a_different_configuration(
            out, "s", {}, output_format="parquet", compress=False
        )


def test_the_table_format_is_not_a_given_setting(tmp_path):
    out = _records(tmp_path, ".parquet")
    with pytest.raises(ValueError, match="unknown setting"):
        refuse_a_different_configuration(
            out, "s", {"table_formats": ".parquet"}, **_FORMAT
        )


def test_a_record_copy_in_a_format_the_field_does_not_name_is_mixed(tmp_path):
    """A TSV run published over a Parquet one leaves two records that disagree; reading only its
    own copy, a `-F tsv` rerun matched and ran (review of the branch)."""
    _records(tmp_path, ".parquet", table_formats=".parquet")
    _records(tmp_path, ".tsv", table_formats=".tsv")
    with pytest.raises(ConfigurationMismatch, match="table_formats") as refused:
        refuse_a_different_configuration(
            tmp_path, "s", {}, output_format="tsv", compress=False
        )
    assert "--output-format" not in str(refused.value)


def test_an_older_two_record_folder_is_not_told_to_rerun_as_both(tmp_path):
    """Two old records may be a `both` run or two runs: `-F both` would keep both runs' stale
    tables, so the only advice is a fresh folder."""
    out = _records(tmp_path, ".parquet", ".tsv")
    with pytest.raises(ConfigurationMismatch) as refused:
        refuse_a_different_configuration(
            out, "s", {}, output_format="parquet", compress=False
        )
    assert "--output-format both" not in str(refused.value)
    assert "fresh output folder" in str(refused.value)
