"""run-all's optional stages call functions that raise, not typer commands.

run-all called the `uxm`, `mfsd` and `region_mds` typer commands as plain functions. Every
option it did not pass stayed a typer `OptionInfo`, so UXM's `threads > 0` raised `TypeError`
on every run-all with a bisulfite BAM, and run-all's catch-all logged it as "UXM failed" and
carried on: UXM has not run under run-all since 0.5.0. Region-MDS was called without
`pon_variant`, so an `OptionInfo` reached the PON lookup in its place. And the commands turn
every failure into `typer.Exit`, which would hide the real error from anything recording it.
"""

from __future__ import annotations

import pytest

from conftest import requires_data

from .runall import synth_run_all, write_vcf
from .synth import BASELINE, make_sample

pytestmark = pytest.mark.invariant


@requires_data
def test_run_all_runs_uxm(tmp_path):
    """Pins the TypeError: before the split no run-all since 0.5.0 produced a UXM table."""
    sample = make_sample(tmp_path / "in", BASELINE)
    out = tmp_path / "out"
    assert synth_run_all(sample, out, bisulfite_bam=sample.bam) is None
    assert (out / "S.UXM.parquet").exists(), "run-all's UXM step did not run"


@requires_data
def test_run_all_passes_its_pon_variant_to_region_mds(tmp_path, monkeypatch):
    """run-all resolved its PON for `--pon-variant`; Region-MDS resolved again for the default."""
    import krewlyzer.region_mds as region_mds_module

    seen = []
    real = region_mds_module.resolve_pon_model

    def spy(*args, **kwargs):
        seen.append(kwargs.get("variant"))
        return real(*args, **kwargs)

    monkeypatch.setattr(region_mds_module, "resolve_pon_model", spy)
    sample = make_sample(tmp_path / "in", BASELINE)
    synth_run_all(sample, tmp_path / "out", assay="xs1", pon_variant="duplex")
    assert seen == ["duplex"]


def test_the_stage_bodies_raise_the_real_error(tmp_path):
    """A recorded failure needs the real class; the commands convert everything to `Exit`."""
    from krewlyzer.mfsd import mfsd_stage
    from krewlyzer.region_mds import region_mds_stage
    from krewlyzer.uxm import uxm_stage

    sample = make_sample(tmp_path / "in", BASELINE)
    missing = tmp_path / "absent.bam"
    with pytest.raises(FileNotFoundError):
        uxm_stage(
            missing,
            tmp_path / "out",
            "S",
            None,
            "hg19",
            output_format="tsv",
            compress=False,
        )
    with pytest.raises(FileNotFoundError):
        mfsd_stage(
            missing,
            write_vcf(tmp_path / "v.vcf"),
            tmp_path / "out",
            "S",
            output_format="tsv",
            compress=False,
        )
    with pytest.raises(FileNotFoundError):
        mfsd_stage(
            sample.bam,
            tmp_path / "absent.vcf",
            tmp_path / "out",
            "S",
            output_format="tsv",
            compress=False,
        )
    with pytest.raises(FileNotFoundError):
        region_mds_stage(
            missing,
            sample.reference,
            tmp_path / "out",
            sample_name="S",
            output_format="tsv",
            compress=False,
        )


def test_the_commands_still_exit_on_failure(tmp_path):
    """The CLI behaviour is unchanged: a failing command exits non-zero, without a traceback."""
    from typer.testing import CliRunner

    from krewlyzer.cli import app

    runner = CliRunner()
    missing = str(tmp_path / "absent.bam")
    for args in (
        ["uxm", "-i", missing, "-o", str(tmp_path / "o"), "-s", "S"],
        ["mfsd", "-i", missing, "-V", missing, "-o", str(tmp_path / "o"), "-s", "S"],
    ):
        result = runner.invoke(app, args)
        assert result.exit_code == 1, (args[0], result.output)
        assert not isinstance(result.exception, FileNotFoundError), args[0]


def test_region_mds_keeps_its_exit_codes(tmp_path):
    """Usage errors (a missing argument path) exit 2, as `BadParameter` did before; anything
    else in the stage -- a malformed gene BED, an unknown genome -- exits 1, as before.
    """
    from typer.testing import CliRunner

    from krewlyzer.cli import app

    sample = make_sample(tmp_path / "in", BASELINE)
    malformed = tmp_path / "genes.bed"
    malformed.write_text("chr1\tnot-a-number\t100\tGENE\n")
    base = [str(sample.bam), str(sample.reference), str(tmp_path / "o"), "-s", "S"]
    runner = CliRunner()
    missing_bam = [
        "region-mds",
        str(tmp_path / "absent.bam"),
        *base[1:],
        "--gene-bed",
        str(sample.assets.genes),
    ]
    assert runner.invoke(app, missing_bam).exit_code == 2
    assert (
        runner.invoke(
            app, ["region-mds", *base, "--gene-bed", str(malformed)]
        ).exit_code
        == 1
    )
    unknown_genome = [
        "region-mds",
        *base,
        "--gene-bed",
        str(sample.assets.genes),
        "-G",
        "hg99",
    ]
    assert runner.invoke(app, unknown_genome).exit_code == 1
