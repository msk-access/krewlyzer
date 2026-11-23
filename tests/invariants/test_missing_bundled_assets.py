"""A named assay or the default blacklist whose bundled file is missing is an error.

A genome the AssetManager knows, with a data folder that lacks the files: a pip install
without the data clone (the wheel has no `data/`), a folder laid out for another release, or
a clone that was never `git lfs pull`ed, whose files are pointers rather than data.
`resolve_target_regions` took a missing targets file for a named panel assay as WGS, at
debug level, and `resolve_exclude_regions` only warned that the blacklist was missing. So
`extract -A xs1` and `run-all -A xs1` exited 0, with ``assay=xs1`` and ``panel_mode=False``
in the record and no blacklist applied: a panel sample measured as WGS, with nothing to say
so (AGENTS.md invariant #1).

``KREWLYZER_DATA_DIR`` points here at a folder that has every category and genome directory
but none of the files, so the AssetManager accepts it, the genome is known, and only the
files are missing. An explicit `--target-regions` that does not exist is the same defect by
another route, and is refused the same way.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from typer.testing import CliRunner

from krewlyzer.cli import app

from .runall import SAMPLE, lfs_pointer, record, synth_run_all
from .synth import BASELINE, make_sample

pytestmark = pytest.mark.invariant

#: `extract` resolves targets and the blacklist whatever `--gc-correct` says.
GC_FLAGS = pytest.mark.parametrize("gc_flag", ["--gc-correct", "--no-gc-correct"])

#: Every category the AssetManager reads, laid out as `data/{category}/{genome}/`.
CATEGORIES = (
    "ATAC",
    "ChromosomeArms",
    "ChromosomeBins",
    "MethMark",
    "OpenChromatinRegion",
    "TFBS",
    "WpsAnchors",
    "WpsBackground",
    "exclude-regions",
    "gc",
    "genes",
    "pon",
    "targets",
)


@pytest.fixture(scope="module")
def sample(tmp_path_factory):
    return make_sample(tmp_path_factory.mktemp("missing_assets_in"), BASELINE)


@pytest.fixture
def empty_data_dir(tmp_path, monkeypatch) -> Path:
    data = tmp_path / "data"
    for category in CATEGORIES:
        for genome in ("GRCh37", "GRCh38"):
            (data / category / genome).mkdir(parents=True)
    monkeypatch.setenv("KREWLYZER_DATA_DIR", str(data))
    return data


@pytest.fixture
def pointer_data_dir(empty_data_dir) -> Path:
    """A clone that skipped `git lfs pull`: the files are there, as LFS pointers."""
    for path in (
        empty_data_dir / "targets" / "GRCh37" / "xs1.targets.bed.gz",
        empty_data_dir / "exclude-regions" / "GRCh37" / "hg19-blacklist.v2.bed.gz",
    ):
        lfs_pointer(path)
    return empty_data_dir


@pytest.fixture
def dataless_package(tmp_path, monkeypatch) -> Path:
    """A pip install without the data clone: the package has no `data/` folder at all."""
    import krewlyzer.assets as assets_module

    package = tmp_path / "site-packages" / "krewlyzer"
    package.mkdir(parents=True)
    monkeypatch.setattr(assets_module, "__file__", str(package / "assets.py"))
    monkeypatch.delenv("KREWLYZER_DATA_DIR", raising=False)
    return package


def _said(text: str, *phrases: str) -> bool:
    """Every phrase is in the console text. Rich wraps at the terminal width and folds a word
    longer than a line (a path), so compare with all whitespace removed."""
    flat = "".join(text.split())
    return all("".join(p.split()) in flat for p in phrases)


def _refused(result) -> None:
    """Exit 1 through `typer.Exit`. An uncaught exception is exit 1 under CliRunner too, with
    no traceback in the output, so the exit code alone cannot tell the two apart."""
    assert result.exit_code == 1, result.output
    assert isinstance(result.exception, SystemExit), repr(result.exception)


def _extract(sample, out: Path, *extra):
    args = [
        "extract",
        "-i",
        sample.bam,
        "--reference",
        sample.reference,
        "--output",
        out,
        "--sample-name",
        SAMPLE,
        "--output-format",
        "parquet",
        "--threads",
        1,
        *extra,
    ]
    return CliRunner().invoke(app, [str(a) for a in args])


# --------------------------------------------------------------------------- assay targets


@GC_FLAGS
def test_extract_refuses_a_named_assay_without_its_targets(
    sample, tmp_path, empty_data_dir, gc_flag
):
    out = tmp_path / "out"
    result = _extract(sample, out, "-A", "xs1", "-x", sample.assets.blacklist, gc_flag)
    _refused(result)
    assert _said(result.output, "xs1", "--target-regions")
    assert not (out / f"{SAMPLE}.bed.gz").exists()


def test_run_all_refuses_a_named_assay_without_its_targets(
    sample, tmp_path, empty_data_dir, capsys
):
    out = tmp_path / "out"
    code = synth_run_all(sample, out, assay="xs1", skip_target_regions=False)
    printed = capsys.readouterr().err  # the console writes to stderr
    assert code == 1, printed
    assert _said(printed, "xs1", "--target-regions")
    assert not (out / f"{SAMPLE}.metadata.parquet").exists()


@GC_FLAGS
def test_wgs_needs_no_targets(sample, tmp_path, empty_data_dir, gc_flag):
    """The control: `wgs` has no bundled targets by design, and is still WGS."""
    out = tmp_path / "out"
    result = _extract(sample, out, "-A", "wgs", "-x", sample.assets.blacklist, gc_flag)
    assert result.exit_code == 0, result.output
    assert bool(record(out)["panel_mode"]) is False

    code = synth_run_all(
        sample, tmp_path / "run_all", assay="wgs", skip_target_regions=False
    )
    assert code is None
    row = record(tmp_path / "run_all")
    assert row["assay"] == "wgs" and bool(row["panel_mode"]) is False


def test_skip_target_regions_still_runs_a_named_assay_as_wgs(
    sample, tmp_path, empty_data_dir
):
    """What the error tells the user to do, if WGS is what they meant."""
    code = synth_run_all(
        sample, tmp_path / "out", assay="xs1", skip_target_regions=True
    )
    assert code is None
    assert bool(record(tmp_path / "out")["panel_mode"]) is False


# --------------------------------------------------------------------------- blacklist


@GC_FLAGS
def test_extract_refuses_to_run_without_the_genomes_blacklist(
    sample, tmp_path, empty_data_dir, gc_flag
):
    out = tmp_path / "out"
    result = _extract(sample, out, gc_flag)
    _refused(result)
    assert _said(result.output, "blacklist", "--exclude-regions")
    assert not (out / f"{SAMPLE}.bed.gz").exists()


def test_run_all_refuses_to_run_without_the_genomes_blacklist(
    sample, tmp_path, empty_data_dir, capsys
):
    out = tmp_path / "out"
    code = synth_run_all(sample, out, exclude_regions=None)
    printed = capsys.readouterr().err  # the console writes to stderr
    assert code == 1, printed
    assert _said(printed, "blacklist", "--exclude-regions")
    assert not (out / f"{SAMPLE}.metadata.parquet").exists()


@GC_FLAGS
def test_an_explicit_exclude_file_needs_no_bundled_blacklist(
    sample, tmp_path, empty_data_dir, gc_flag
):
    """`-x` is the way out the error names, an empty file included (it excludes nothing)."""
    empty = tmp_path / "nothing.bed"
    empty.write_text("")
    result = _extract(sample, tmp_path / "out", "-x", empty, gc_flag)
    assert result.exit_code == 0, result.output
    assert synth_run_all(sample, tmp_path / "run_all", exclude_regions=empty) is None


# --------------------------------------------------------------------------- LFS pointers


@GC_FLAGS
def test_extract_refuses_targets_that_are_lfs_pointers(
    sample, tmp_path, pointer_data_dir, gc_flag
):
    """A pointer exists, so an existence check passed it: `panel_mode=True` on a text file."""
    out = tmp_path / "out"
    result = _extract(sample, out, "-A", "xs1", "-x", sample.assets.blacklist, gc_flag)
    _refused(result)
    assert _said(result.output, "git lfs pull")
    assert not (out / f"{SAMPLE}.bed.gz").exists()


@GC_FLAGS
def test_extract_refuses_a_blacklist_that_is_an_lfs_pointer(
    sample, tmp_path, pointer_data_dir, gc_flag
):
    """It was applied as given and excluded nothing: every blacklisted fragment kept, exit 0."""
    out = tmp_path / "out"
    result = _extract(sample, out, gc_flag)
    _refused(result)
    assert _said(result.output, "git lfs pull")
    assert not (out / f"{SAMPLE}.bed.gz").exists()


def test_run_all_refuses_lfs_pointers(sample, tmp_path, pointer_data_dir, capsys):
    code = synth_run_all(sample, tmp_path / "a", assay="xs1", skip_target_regions=False)
    assert code == 1 and _said(capsys.readouterr().err, "git lfs pull")
    code = synth_run_all(sample, tmp_path / "b", exclude_regions=None)
    assert code == 1 and _said(capsys.readouterr().err, "git lfs pull")


# --------------------------------------------------------------------------- explicit targets


@GC_FLAGS
def test_extract_refuses_an_explicit_targets_file_that_does_not_exist(
    sample, tmp_path, gc_flag
):
    """It was warned about and then run as WGS: `assay=xs1`, `panel_mode=False`, exit 0."""
    out = tmp_path / "out"
    missing = tmp_path / "typo.bed"
    result = _extract(
        sample, out, "-A", "xs1", "-T", missing, "-x", sample.assets.blacklist, gc_flag
    )
    _refused(result)
    assert _said(result.output, "typo.bed")
    assert not (out / f"{SAMPLE}.bed.gz").exists()


def test_run_all_refuses_an_explicit_targets_file_that_does_not_exist(
    sample, tmp_path, capsys
):
    out = tmp_path / "out"
    code = synth_run_all(
        sample,
        out,
        assay="xs1",
        skip_target_regions=False,
        target_regions=tmp_path / "typo.bed",
    )
    assert code == 1 and _said(capsys.readouterr().err, "typo.bed")
    assert not (out / f"{SAMPLE}.metadata.parquet").exists()


# --------------------------------------------------------------------------- no data folder


def test_a_dataless_install_refuses_a_named_assay_and_says_why(
    sample, tmp_path, dataless_package
):
    out = tmp_path / "out"
    result = _extract(sample, out, "-A", "xs1", "-x", sample.assets.blacklist)
    _refused(result)
    assert _said(result.output, f"no data folder at {dataless_package / 'data'}")
    assert _said(result.output, "pip install")


def test_a_dataless_install_runs_when_nothing_bundled_is_needed(
    sample, tmp_path, dataless_package
):
    """Explicit inputs and no bundled lookup: the folder's absence is no reason to stop. A
    check in the AssetManager's constructor stopped this, and every command with it."""
    result = _extract(
        sample,
        tmp_path / "out",
        "-A",
        "wgs",
        "-x",
        sample.assets.blacklist,
        "--no-gc-correct",
    )
    assert result.exit_code == 0, result.output
