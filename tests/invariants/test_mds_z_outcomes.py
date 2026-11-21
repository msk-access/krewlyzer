"""The MDS z-scores carry on past their own failure, so they record what happened.

`mds_z` is a column, not a table, which is why the contract never required it: a PON-scored
sample whose z-score failed passed every check. Both whole-sample z-scores (off- and on-target)
and Region-MDS's gene z-scores now come back as an outcome alongside the tables.
"""

from __future__ import annotations

from types import SimpleNamespace

import pandas as pd
import pytest

from krewlyzer.pon.model import MdsBaseline

from .synth import BASELINE, make_sample

pytestmark = pytest.mark.invariant


@pytest.fixture(scope="module")
def extracted(tmp_path_factory):
    from krewlyzer.core.sample_processor import extract_sample

    sample = make_sample(tmp_path_factory.mktemp("mds_in"), BASELINE)
    result = extract_sample(
        input_path=sample.bam,
        reference=sample.reference,
        threads=1,
        target_regions=sample.assets.targets,
        write_bed=False,
        output_dir=None,
        sample_name="S",
        genome="hg19",
    )
    assert result.is_panel_mode and sum(result.em_counts_ontarget.values()) > 0
    return result


def _pon(mean=0.9, std=0.05, ontarget=None):
    return SimpleNamespace(
        mds_baseline=MdsBaseline(mds_mean=mean, mds_std=std),
        mds_baseline_ontarget=ontarget,
        breakpoint_motif_baseline=None,
        breakpoint_motif_baseline_ontarget=None,
    )


def _write(extracted, out, pon, fmt="parquet"):
    from krewlyzer.core.sample_processor import write_motif_tables

    _, outcomes = write_motif_tables(
        extracted, out, pon=pon, include_ontarget=True, output_format=fmt
    )
    return outcomes


def test_a_z_score_that_is_computed_records_ran(extracted, tmp_path):
    outcomes = _write(extracted, tmp_path, _pon())
    assert outcomes["mds_z"].status == "ran"
    assert (
        outcomes["mds_z_ontarget"].status == "ran"
    ), "falls back to the genome-wide baseline"
    assert "mds_z" in pd.read_parquet(tmp_path / "S.MDS.parquet").columns


def test_a_malformed_baseline_records_failed_and_the_tables_are_still_written(
    extracted, tmp_path
):
    """Two inputs, two outputs: the same sample against a baseline that cannot be used."""
    outcomes = _write(extracted, tmp_path, _pon(mean="not-a-number"))
    assert outcomes["mds_z"].status == "failed"
    assert outcomes["mds_z"].error_class.endswith("TypeError")  # numpy's UFuncTypeError
    assert outcomes["mds_z_ontarget"].status == "failed"
    assert "mds_z" not in pd.read_parquet(tmp_path / "S.MDS.parquet").columns
    assert (tmp_path / "S.MDS.ontarget.parquet").exists()


def test_a_pon_without_an_mds_baseline_records_unavailable(extracted, tmp_path):
    pon = _pon()
    pon.mds_baseline = None
    outcomes = _write(extracted, tmp_path, pon)
    assert outcomes["mds_z"].status == "unavailable"
    assert outcomes["mds_z_ontarget"].status == "unavailable"


def test_no_pon_records_not_requested(extracted, tmp_path):
    outcomes = _write(extracted, tmp_path, None)
    assert outcomes["mds_z"].status == "not_requested"
    assert outcomes["mds_z_ontarget"].status == "not_requested"


def test_the_off_target_z_score_reads_the_table_this_run_wrote(extracted, tmp_path):
    """The append read Parquet-first: a TSV run beside an earlier run's `.MDS.parquet` read that,
    found its `mds_z`, and skipped the append -- the fresh TSV shipped with no z-score.
    """
    pd.DataFrame({"sample": ["S"], "MDS": [0.5], "mds_z": [9.9]}).to_parquet(
        tmp_path / "S.MDS.parquet"
    )
    _write(extracted, tmp_path, _pon(), fmt="tsv")
    fresh = pd.read_csv(tmp_path / "S.MDS.tsv", sep="\t", comment="#")
    assert "mds_z" in fresh.columns


# ---------------------------------------------------------------------------
# Region-MDS: its gene z-scoring is an outcome of its own
# ---------------------------------------------------------------------------


def _region_pon(path, with_region_mds: bool):
    rows = [
        {
            "table": "metadata",
            "schema_version": "1.0",
            "assay": "xs1",
            "n_samples": 21.0,
            "krewlyzer_version": "0.9.0",
        }
    ]
    if with_region_mds:
        rows += [
            {
                "table": "region_mds",
                "gene": f"GENE{i}",
                "mds_mean": 0.9,
                "mds_std": 0.05 + 0.01 * i,
                "mds_e1_mean": 0.9,
                "mds_e1_std": 0.06,
                "n_samples": 21,
            }
            for i in range(2)
        ]
    pd.DataFrame(rows).to_parquet(path)
    return path


@pytest.mark.parametrize(
    "with_baseline, expected", [(True, "ran"), (False, "unavailable")]
)
def test_region_mds_returns_what_its_z_scoring_did(tmp_path, with_baseline, expected):
    from krewlyzer.region_mds import region_mds_stage

    sample = make_sample(tmp_path / "in", BASELINE)
    pon = _region_pon(tmp_path / "r.pon.parquet", with_baseline)
    out = tmp_path / "out"
    *_, z = region_mds_stage(
        sample.bam,
        sample.reference,
        out,
        gene_bed=sample.assets.genes,
        sample_name="S",
        pon_model=pon,
        output_format="parquet",
        silent=True,
    )
    assert z.status == expected
    columns = pd.read_parquet(out / "S.MDS.gene.parquet").columns
    assert ("mds_z" in columns) is (expected == "ran")


def test_region_mds_without_a_pon_did_not_request_z_scores(tmp_path):
    from krewlyzer.region_mds import region_mds_stage

    sample = make_sample(tmp_path / "in", BASELINE)
    *_, z = region_mds_stage(
        sample.bam,
        sample.reference,
        tmp_path / "out",
        gene_bed=sample.assets.genes,
        sample_name="S",
        output_format="parquet",
        silent=True,
    )
    assert z.status == "not_requested"


def test_region_mds_z_scoring_that_fails_is_recorded_not_swallowed(
    tmp_path, monkeypatch
):
    """The one patched failure in this suite. Every input that reaches the gene z-score has been
    through the PON loader, which refuses a malformed baseline; what is left is the write
    failing. It used to log a warning and count as scored."""
    from krewlyzer.pon.model import RegionMdsBaseline
    from krewlyzer.region_mds import region_mds_stage

    def broken(self, gene, value):
        raise ZeroDivisionError("a baseline no loader would accept")

    monkeypatch.setattr(RegionMdsBaseline, "compute_zscore", broken)
    sample = make_sample(tmp_path / "in", BASELINE)
    *_, z = region_mds_stage(
        sample.bam,
        sample.reference,
        tmp_path / "out",
        gene_bed=sample.assets.genes,
        sample_name="S",
        pon_model=_region_pon(tmp_path / "r.pon.parquet", True),
        output_format="parquet",
        silent=True,
    )
    assert (z.status, z.error_class) == ("failed", "ZeroDivisionError")
