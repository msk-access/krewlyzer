"""A rerun into a sample's folder never changes the fragment BED under tables made from the old one.

`{sample}.metadata.parquet` is downstream's completion marker (invariant #2): whatever sits beside
it is read as one sample. Two reruns broke that while every stage still said `ran`:

- **run-all reused a BED extracted with other filters.** `extract --mapq 30` then `run-all` (mapq
  20) exited 0: the BED existed, so run-all kept it and extract's record, then read the BAM again
  at mapq 20 for the motif tables. Two filter sets under one record. The filters are now part of
  the configuration a rerun is checked against, and a BED whose record cannot vouch for its
  filters (written before they were recorded truthfully, or missing) is extracted again.
- **`extract` inferred "nothing was made from this BED" from the record's columns.** A run-all
  killed before its final stamp, or a folder written before provenance existed, holds feature
  tables beside a record that looks like extraction's alone, and `extract` with new filters
  replaced the BED and the marker beside them. It now also refuses when any table the output
  contract declares exists for the sample.
"""

from __future__ import annotations

import gzip
import shutil

import pandas as pd
import pytest
from typer.testing import CliRunner

from conftest import requires_data
from krewlyzer.cli import app
from krewlyzer.validate.contract import MOTIF_TABLES, STAGE_TABLES

from .runall import SAMPLE, record, synth_run_all, write_vcf
from .synth import BASELINE, make_sample

pytestmark = [pytest.mark.invariant, requires_data]

#: Each filter: what `extract` is given (not the default), run-all's option at its default, the
#: record column, and a second run-all value.
FILTERS = {
    "mapq": (("--mapq", "30"), ("mapq", 20), "filters_mapq", 30),
    "minlen": (("--minlen", "80"), ("minlen", 65), "filters_min_length", 80),
    "maxlen": (("--maxlen", "600"), ("maxlen", 1000), "filters_max_length", 600),
    "skip_duplicates": (
        ("--no-skip-duplicates",),
        ("skip_duplicates", True),
        "filters_skip_duplicates",
        False,
    ),
    "require_proper_pair": (
        ("--no-require-proper-pair",),
        ("require_proper_pair", True),
        "filters_require_proper_pair",
        False,
    ),
}

FILTER_COLUMNS = [column for _, _, column, _ in FILTERS.values()]


@pytest.fixture(scope="module")
def sample(tmp_path_factory):
    return make_sample(tmp_path_factory.mktemp("rerun_bed_in"), BASELINE)


def _files(out):
    return {p.name: p.stat().st_mtime_ns for p in sorted(out.iterdir()) if p.is_file()}


def _cli(*args):
    return CliRunner().invoke(app, [str(a) for a in args])


def _extract_args(sample, out, blacklist, *extra):
    return (
        "extract",
        "-i",
        sample.bam,
        "-r",
        sample.reference,
        "-o",
        out,
        "-s",
        SAMPLE,
        "-x",
        blacklist,
        "--skip-target-regions",
        "-F",
        "parquet",
        *extra,
    )


def _extract(sample, out, *extra):
    """`extract` with the settings `synth_run_all` uses, so only what a test adds differs."""
    return _cli(*_extract_args(sample, out, sample.assets.blacklist, *extra))


def _motif_args(sample, out, blacklist, *extra):
    """`motif` with the settings `synth_run_all` uses, so only what a test adds differs."""
    return (
        "motif",
        "-i",
        sample.bam,
        "-r",
        sample.reference,
        "-o",
        out,
        "-s",
        SAMPLE,
        "-x",
        blacklist,
        "--skip-target-regions",
        "--skip-pon",
        "-F",
        "parquet",
        *extra,
    )


def _no_blacklist(directory):
    """A blacklist of nothing: another blacklist than the synthetic one, as `-x` may give."""
    path = directory / "nothing.bed"
    path.write_text("")
    return path


def _copied_blacklist(sample, directory):
    """The synthetic blacklist at another path: the same blacklist."""
    copied = directory / "elsewhere" / sample.assets.blacklist.name
    copied.parent.mkdir()
    shutil.copy(sample.assets.blacklist, copied)
    return copied


def _bed_fragments(out):
    with gzip.open(out / f"{SAMPLE}.bed.gz", "rt") as fh:
        return sum(1 for _ in fh)


def _as_an_old_record(out):
    """Rewrite the record as one from before the filters were recorded truthfully: no
    duplicate/proper-pair columns, and 65 whatever the run used."""
    path = out / f"{SAMPLE}.metadata.parquet"
    frame = pd.read_parquet(path).drop(
        columns=["filters_skip_duplicates", "filters_require_proper_pair"]
    )
    frame.assign(filters_min_length=65).to_parquet(path, index=False)


# ── 1. The filters are part of the configuration ──────────────────────────────────────────────


@pytest.mark.parametrize("name", list(FILTERS))
def test_run_all_refuses_a_bed_extracted_with_other_filters(sample, tmp_path, name):
    """The review's case, for each filter: run-all kept extract's BED and its record, and wrote
    motif tables from the BAM under its own filters."""
    extract_args, (field, default), column, _ = FILTERS[name]
    out = tmp_path / "out"
    result = _extract(sample, out, *extract_args)
    assert result.exit_code == 0, result.output
    before = _files(out)

    code = synth_run_all(sample, out, **{field: default})

    assert code == 1, f"run-all accepted a BED made with another {name}"
    assert _files(out) == before, "a refused run-all wrote into the folder"


def test_run_all_with_the_same_filters_reuses_the_bed(sample, tmp_path):
    """The refusal is for other filters: with the same ones, run-all keeps extract's BED."""
    out = tmp_path / "out"
    assert _extract(sample, out, "--minlen", "80").exit_code == 0
    bed = (out / f"{SAMPLE}.bed.gz").stat().st_mtime_ns

    assert synth_run_all(sample, out, minlen=80) is None

    assert (out / f"{SAMPLE}.bed.gz").stat().st_mtime_ns == bed, "BED re-extracted"
    assert int(record(out)["filters_min_length"]) == 80


@pytest.mark.parametrize("name", list(FILTERS))
def test_run_all_records_the_filters_it_used(sample, tmp_path, name):
    """Two run-alls that differ in one filter record two different values (invariant #1)."""
    _, (field, default), column, other = FILTERS[name]
    recorded = []
    for value in (default, other):
        out = tmp_path / str(value)
        assert synth_run_all(sample, out, **{field: value}) is None
        recorded.append(record(out)[column])
    assert [type(default)(v) for v in recorded] == [default, other]


def test_extract_records_the_duplicate_and_proper_pair_filters(sample, tmp_path):
    """The two filters the record held no column for (invariant #1)."""
    recorded = {}
    for name, extra in (
        ("default", ()),
        ("lenient", ("--no-skip-duplicates", "--no-require-proper-pair")),
    ):
        out = tmp_path / name
        assert _extract(sample, out, *extra).exit_code == 0
        row = record(out)
        recorded[name] = (
            bool(row["filters_skip_duplicates"]),
            bool(row["filters_require_proper_pair"]),
        )
    assert recorded == {"default": (True, True), "lenient": (False, False)}


def test_extract_records_which_blacklist_it_applied(sample, tmp_path):
    """The blacklist decides which fragments the BED holds, as the read filters do: two extracts
    with different blacklists record two different values, and the same file at another path the
    same one (invariant #1). By content, not name: two files can share a name, and the bundled one
    sits at another path in every install."""
    recorded = {}
    for name, blacklist in (
        ("synthetic", sample.assets.blacklist),
        ("copy", _copied_blacklist(sample, tmp_path)),
        ("nothing", _no_blacklist(tmp_path)),
    ):
        out = tmp_path / name
        assert _cli(*_extract_args(sample, out, blacklist)).exit_code == 0
        recorded[name] = record(out)["filters_blacklist_sha256"]
    assert recorded["synthetic"] == recorded["copy"]
    assert recorded["synthetic"] != recorded["nothing"]
    assert all(isinstance(v, str) and len(v) == 64 for v in recorded.values())


def test_run_all_refuses_a_bed_extracted_with_another_blacklist(sample, tmp_path):
    """The second review's case: `extract -x nothing.bed`, then run-all with the synthetic
    blacklist, kept a BED holding the fragments run-all would have dropped."""
    out = tmp_path / "out"
    assert _cli(*_extract_args(sample, out, _no_blacklist(tmp_path))).exit_code == 0
    before = _files(out)

    code = synth_run_all(sample, out)

    assert code == 1, "run-all accepted a BED made with another blacklist"
    assert _files(out) == before, "a refused run-all wrote into the folder"


def test_run_all_reuses_a_bed_extracted_with_the_same_blacklist_elsewhere(
    sample, tmp_path
):
    out = tmp_path / "out"
    blacklist = _copied_blacklist(sample, tmp_path)
    assert _cli(*_extract_args(sample, out, blacklist)).exit_code == 0
    bed = (out / f"{SAMPLE}.bed.gz").stat().st_mtime_ns

    assert synth_run_all(sample, out) is None

    assert (out / f"{SAMPLE}.bed.gz").stat().st_mtime_ns == bed, "BED re-extracted"


def test_extract_with_another_blacklist_than_its_record_is_refused(sample, tmp_path):
    out = tmp_path / "out"
    assert _extract(sample, out).exit_code == 0
    before = _files(out)

    result = _cli(*_extract_args(sample, out, _no_blacklist(tmp_path)))

    assert result.exit_code == 1, result.output
    assert "filters_blacklist_sha256" in result.output
    assert _files(out) == before


def test_motif_with_another_blacklist_is_refused(sample, tmp_path):
    """`motif` applies the blacklist itself, reading the BAM: with another one its tables would
    count fragments the sample's BED does not hold."""
    out = tmp_path / "out"
    assert synth_run_all(sample, out) is None
    before = _files(out)

    result = _cli(*_motif_args(sample, out, _no_blacklist(tmp_path)))

    assert result.exit_code != 0, result.output
    assert "filters_blacklist_sha256" in result.output
    assert _files(out) == before


def test_extract_with_other_filters_than_its_record_is_refused(sample, tmp_path):
    """As for the genome: the folder's filters are the folder's, whichever tool wrote them."""
    out = tmp_path / "out"
    assert _extract(sample, out, "--minlen", "80").exit_code == 0
    before = _files(out)

    result = _extract(sample, out)

    assert result.exit_code == 1, result.output
    assert "filters_min_length" in result.output
    assert _files(out) == before


def test_motif_with_other_filters_is_refused(sample, tmp_path):
    """`motif` reads the BAM, not the BED: with other filters its tables would come from
    fragments the sample's BED does not hold."""
    out = tmp_path / "out"
    assert synth_run_all(sample, out) is None
    before = _files(out)

    result = _cli(*_motif_args(sample, out, sample.assets.blacklist, "--minlen", "80"))

    assert result.exit_code != 0, result.output
    assert "filters_min_length" in result.output
    assert _files(out) == before


def _region_mds_args(sample, out, *extra):
    return (
        "region-mds",
        sample.bam,
        sample.reference,
        out,
        "-s",
        SAMPLE,
        "--skip-pon",
        "-F",
        "parquet",
        *extra,
    )


def test_region_mds_with_other_filters_is_refused(sample, tmp_path):
    out = tmp_path / "out"
    assert synth_run_all(sample, out) is None
    before = _files(out)

    result = _cli(*_region_mds_args(sample, out, "--minlen", "80"))

    assert result.exit_code != 0, result.output
    assert "filters_min_length" in result.output
    assert _files(out) == before


def _mfsd(sample, out, vcf, *extra):
    return _cli(
        "mfsd",
        "-i",
        sample.bam,
        "-V",
        vcf,
        "-g",
        sample.reference,
        "-o",
        out,
        "-s",
        SAMPLE,
        "--output-format",
        "parquet",
        *extra,
    )


@pytest.mark.parametrize(
    "option,column",
    [
        (("--mapq", "30"), "filters_mapq"),
        (("--minlen", "80"), "filters_min_length"),
        (("--maxlen", "600"), "filters_max_length"),
    ],
)
def test_mfsd_with_other_filters_is_refused(sample, tmp_path, option, column):
    """`mfsd` reads the BAM, as `motif` does, and run-all gives it the sample's mapq and length
    range (`--mfsd-bam` included). It is not compared on the other two: run-all always runs it
    without the proper-pair filter, and it drops duplicates whatever its flag says."""
    out = tmp_path / "out"
    assert synth_run_all(sample, out) is None
    vcf = write_vcf(tmp_path / "variants.vcf")
    before = _files(out)

    result = _mfsd(sample, out, vcf, *option)

    assert result.exit_code == 1, result.output
    assert column in result.output
    assert _files(out) == before


def test_mfsd_with_the_folders_filters_runs(sample, tmp_path):
    out = tmp_path / "out"
    assert synth_run_all(sample, out) is None
    vcf = write_vcf(tmp_path / "variants.vcf")

    result = _mfsd(sample, out, vcf)

    assert result.exit_code == 0, result.output
    assert (out / f"{SAMPLE}.mFSD.parquet").exists()


@pytest.mark.parametrize("state", ["record_before_filters_were_recorded", "no_record"])
def test_run_all_re_extracts_a_bed_its_record_cannot_vouch_for(sample, tmp_path, state):
    """A record from before the filters were recorded truthfully said 20/65/1000 whatever the run
    used, and has no `filters_skip_duplicates`: its values cannot be compared, so the guard lets
    the run through. The BED it describes is not trusted either: run-all extracts again under its
    own filters, rather than scoring a BED made with unknown ones. A BED with no record beside it
    (a run killed before writing one) is the same case."""
    fresh = tmp_path / "fresh"
    assert synth_run_all(sample, fresh) is None
    expected = _bed_fragments(fresh)

    out = tmp_path / "out"
    assert _extract(sample, out, "--minlen", "150").exit_code == 0
    assert _bed_fragments(out) < expected, "the two filters keep the same fragments"
    if state == "no_record":
        (out / f"{SAMPLE}.metadata.parquet").unlink()
    else:
        _as_an_old_record(out)

    assert synth_run_all(sample, out) is None

    assert (
        _bed_fragments(out) == expected
    ), "run-all scored a BED made with other filters"
    row = record(out)
    assert int(row["total_fragments"]) == expected
    assert all(c in row for c in FILTER_COLUMNS)


def test_an_old_record_is_not_refused_for_its_filters(sample, tmp_path):
    """Its 65 may not be what the run used, so a run with another value is let through (and
    re-extracts, above) rather than refused for a value nobody can trust."""
    out = tmp_path / "out"
    assert _extract(sample, out, "--minlen", "80").exit_code == 0
    _as_an_old_record(out)

    assert synth_run_all(sample, out, minlen=80) is None

    assert int(record(out)["filters_min_length"]) == 80


@pytest.mark.parametrize("bed", ["kept", "deleted"])
def test_run_all_does_not_re_extract_a_bed_tables_were_made_from(sample, tmp_path, bed):
    """A run-all folder from before the filters were recorded truthfully: its BED cannot be
    trusted, and cannot be replaced either, since every table beside it was made from it. A run
    killed after extracting would leave a new BED and a fresh record over the old tables. Refused,
    as `extract` is, before anything is written. The same when the BED was deleted to save space:
    the tables were still made from it, with filters nobody knows."""
    out = tmp_path / "out"
    assert synth_run_all(sample, out) is None
    _as_an_old_record(out)
    if bed == "deleted":
        (out / f"{SAMPLE}.bed.gz").unlink()
    before = _files(out)

    code = synth_run_all(sample, out)

    assert code == 1, "run-all replaced a BED its tables were made from"
    assert _files(out) == before


def test_run_all_re_extracts_a_deleted_bed_its_record_vouches_for(sample, tmp_path):
    """The refusal above is for a record that cannot say which filters made the tables. A current
    record can: the configuration guard has refused any other filters, so a deleted BED is
    extracted again with the same ones, beside the tables made from the first."""
    out = tmp_path / "out"
    assert synth_run_all(sample, out, minlen=80) is None
    expected = _bed_fragments(out)
    (out / f"{SAMPLE}.bed.gz").unlink()

    assert synth_run_all(sample, out, minlen=80) is None

    assert _bed_fragments(out) == expected
    assert int(record(out)["filters_min_length"]) == 80


def test_a_tool_reading_the_bam_warns_when_the_filters_cannot_be_compared(
    sample, tmp_path, caplog
):
    """`motif` cannot extract the BED again, as run-all does: over a record whose filters cannot
    be trusted it is let through, and says so."""
    import logging

    out = tmp_path / "out"
    assert _extract(sample, out).exit_code == 0
    _as_an_old_record(out)
    with caplog.at_level(logging.WARNING):
        result = _cli(*_motif_args(sample, out, sample.assets.blacklist))
    assert result.exit_code == 0, result.output
    assert any(
        "filters" in r.getMessage() and r.levelno >= logging.WARNING
        for r in caplog.records
    ), "a filter check that could not be made went unsaid"


# ── 1b. Tables made with filters nobody recorded ─────────────────────────────────────────────
#
# Over a record that cannot say which filters made the sample's fragments, run-all refuses a
# folder holding tables (above); a tool reading the BAM refused nothing, and wrote its own beside
# them. One check now serves all four: `output_provenance.refuse_tables_of_unrecorded_filters`.


def _an_old_run_all_folder(sample, out):
    assert synth_run_all(sample, out) is None
    _as_an_old_record(out)
    return _files(out)


def _refused_for_unrecorded_filters(result, own):
    """Refused by the shared check, naming tables another step made and none of the tool's own."""
    from rich.text import Text

    # As text, and on one line: Rich styles the message under GITHUB_ACTIONS, and wraps it (in a
    # box, for a usage error).
    message = " ".join(Text.from_ansi(result.output).plain.replace("│", " ").split())
    assert result.exit_code != 0, message
    assert "made with read filters its record does not name" in message, message
    assert f"among them {SAMPLE}." in message, message
    assert not any(f"{SAMPLE}{suffix}" in message for suffix in own), message


def test_motif_refuses_a_folder_of_tables_made_with_unrecorded_filters(
    sample, tmp_path
):
    """The second review's case: over an old record, `motif` warned and wrote motif tables from
    the BAM beside tables made from a BED nobody knows the filters of."""
    out = tmp_path / "out"
    before = _an_old_run_all_folder(sample, out)

    result = _cli(*_motif_args(sample, out, sample.assets.blacklist))

    _refused_for_unrecorded_filters(result, MOTIF_TABLES)
    assert _files(out) == before


def test_region_mds_refuses_a_folder_of_tables_made_with_unrecorded_filters(
    sample, tmp_path
):
    out = tmp_path / "out"
    before = _an_old_run_all_folder(sample, out)

    result = _cli(*_region_mds_args(sample, out))

    _refused_for_unrecorded_filters(result, STAGE_TABLES["region_mds"])
    assert _files(out) == before


def test_mfsd_refuses_a_folder_of_tables_made_with_unrecorded_filters(sample, tmp_path):
    out = tmp_path / "out"
    before = _an_old_run_all_folder(sample, out)
    vcf = write_vcf(tmp_path / "variants.vcf")

    result = _mfsd(sample, out, vcf)

    _refused_for_unrecorded_filters(result, STAGE_TABLES["mfsd"])
    assert _files(out) == before


def test_a_tools_own_earlier_tables_do_not_refuse_it(sample, tmp_path):
    """`motif` replaces its own tables, so those alone are no reason to refuse: an extract-and-
    motif folder whose record cannot vouch for its filters is still rerun by `motif`."""
    out = tmp_path / "out"
    assert _extract(sample, out).exit_code == 0
    assert _cli(*_motif_args(sample, out, sample.assets.blacklist)).exit_code == 0
    _as_an_old_record(out)

    result = _cli(*_motif_args(sample, out, sample.assets.blacklist))

    assert result.exit_code == 0, result.output


def test_a_tool_into_a_folder_without_a_record_is_not_refused(sample, tmp_path):
    """Without a record there is nothing to compare, and a single tool writes none
    (`record_stage_outcomes`), so it leaves no completion marker over the tables beside it: tools
    run by hand into one folder still run. run-all writes one, and refuses (above)."""
    out = tmp_path / "out"
    assert synth_run_all(sample, out) is None
    for record_file in out.glob(f"{SAMPLE}.metadata.*"):
        record_file.unlink()

    result = _cli(*_motif_args(sample, out, sample.assets.blacklist))

    assert result.exit_code == 0, result.output


def test_run_all_into_a_folder_of_tables_without_a_record_is_refused(sample, tmp_path):
    """run-all writes the completion marker: over tables nobody recorded the filters of, it would
    vouch for them as one sample's."""
    out = tmp_path / "out"
    assert synth_run_all(sample, out) is None
    for record_file in out.glob(f"{SAMPLE}.metadata.*"):
        record_file.unlink()
    before = _files(out)

    assert synth_run_all(sample, out) == 1
    assert _files(out) == before


# ── 2. extract refuses a folder with tables made from its BED ────────────────────────────────


def _feature_suffixes():
    from krewlyzer.validate.contract import feature_table_suffixes

    return sorted(feature_table_suffixes())


@pytest.fixture(scope="module")
def extracted(sample, tmp_path_factory):
    out = tmp_path_factory.mktemp("extracted") / "out"
    result = _extract(sample, out)
    assert result.exit_code == 0, result.output
    return out


@pytest.mark.parametrize("suffix", _feature_suffixes())
def test_extract_refuses_a_folder_holding_a_feature_table(
    sample, extracted, tmp_path, suffix
):
    """Extraction's own record and one table something read this BED to make. The record alone
    looks like a rerun of `extract`."""
    out = tmp_path / "out"
    shutil.copytree(extracted, out)
    (out / f"{SAMPLE}{suffix}").write_bytes(b"")
    before = _files(out)

    result = _extract(sample, out)

    assert result.exit_code == 1, result.output
    assert f"{SAMPLE}{suffix}" in result.output
    assert _files(out) == before


@pytest.mark.parametrize("extension", [".tsv", ".tsv.gz"])
def test_extract_refuses_a_feature_table_in_any_format(
    sample, extracted, tmp_path, extension
):
    out = tmp_path / "out"
    shutil.copytree(extracted, out)
    (out / f"{SAMPLE}.FSD{extension}").write_bytes(b"")

    result = _extract(sample, out)

    assert result.exit_code == 1, result.output
    assert f"{SAMPLE}.FSD{extension}" in result.output


def test_extract_refuses_a_run_all_folder_written_before_provenance(
    sample, extracted, tmp_path
):
    """A real run-all folder whose record holds extraction's columns only: a run-all killed
    between writing the record and stamping it, or a folder from before provenance existed.
    """
    out = tmp_path / "out"
    assert synth_run_all(sample, out) is None
    keep = set(pd.read_parquet(extracted / f"{SAMPLE}.metadata.parquet").columns)
    path = out / f"{SAMPLE}.metadata.parquet"
    frame = pd.read_parquet(path)
    frame[[c for c in frame.columns if c in keep]].to_parquet(path, index=False)
    before = _files(out)

    result = _extract(sample, out)

    assert result.exit_code == 1, result.output
    assert _files(out) == before


def test_feature_table_suffixes_cover_the_contract_but_not_extractions_own():
    from krewlyzer.validate.contract import (
        CONTRACT,
        NOT_CONSUMED,
        OCF_PREFERENCE,
        STAGE_TABLES,
        feature_table_suffixes,
    )

    declared = (
        {t.suffix for t in CONTRACT}
        | set(NOT_CONSUMED)
        | set(OCF_PREFERENCE)
        | {s for tables in STAGE_TABLES.values() for s in tables}
    )
    assert set(feature_table_suffixes()) == declared - {
        ".metadata.parquet",
        ".correction_factors.parquet",
    }
