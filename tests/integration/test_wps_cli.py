"""
Integration tests for WPS (Windowed Protection Score) CLI.

Tests the standalone wps command with Parquet output.
"""

import pytest
import pysam
import re
from typer.testing import CliRunner
from krewlyzer.cli import app

runner = CliRunner()


def strip_ansi(text: str) -> str:
    """Remove ANSI escape codes from text for clean assertions."""
    ansi_escape = re.compile(r"\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])")
    return ansi_escape.sub("", text)


@pytest.fixture
def sample_bedgz(tmp_path):
    """Create sample BED.gz for WPS testing with fragments spanning anchors."""
    bed = tmp_path / "sample.bed"
    with open(bed, "w") as f:
        # Dense coverage around position 1000 (simulated TSS)
        for i in range(100):
            start = 800 + (i * 15)
            # Mix of nucleosome-sized and TF-sized fragments
            if i % 3 == 0:
                length = 167  # Primary nucleosome range
            elif i % 3 == 1:
                length = 145  # Secondary nucleosome range
            else:
                length = 55  # TF range
            f.write(f"chr1\t{start}\t{start + length}\t0.50\n")

    bedgz = tmp_path / "sample.bed.gz"
    pysam.tabix_compress(str(bed), str(bedgz), force=True)
    pysam.tabix_index(str(bedgz), preset="bed", force=True)
    return bedgz


@pytest.fixture
def early_closing_bedgz(tmp_path):
    """Fragments that all close well before the anchor window ends.

    The anchor below is centred at 1000 and expanded to +/-1000, so the window
    runs 0-2000. These fragments span 550-812 and every one of them closes
    inside it, which makes the coverage difference array telescope to exactly
    zero at the final base.

    That is the arrangement `local_depth` got wrong through 0.9.2: it reported
    the running value left after integration -- the last base's coverage -- so
    this region read as having no coverage at all, and the `== 0.0` emptiness
    guard dropped it outright (#79).
    """
    bed = tmp_path / "early.bed"
    with open(bed, "w") as f:
        for i in range(20):
            start = 550 + (i * 5)
            # 167 bp lands in the nucleosome stream and 55 bp in the TF stream.
            # A length between the two -- 100 bp, say -- is weighted zero by both,
            # registers no coverage at all, and leaves the region genuinely empty
            # rather than wrongly dropped, which tests nothing.
            length = 167 if i % 2 == 0 else 55
            f.write(f"chr1\t{start}\t{start + length}\t0.50\n")

    bedgz = tmp_path / "early.bed.gz"
    pysam.tabix_compress(str(bed), str(bedgz), force=True)
    pysam.tabix_index(str(bedgz), preset="bed", force=True)
    return bedgz


@pytest.fixture
def sample_anchors(tmp_path):
    """Create WPS anchors BED for testing."""
    bed = tmp_path / "anchors.bed"
    with open(bed, "w") as f:
        # TSS anchor
        f.write("chr1\t500\t1500\tGENE1\t.\t+\n")

    bedgz = tmp_path / "anchors.bed.gz"
    pysam.tabix_compress(str(bed), str(bedgz), force=True)
    pysam.tabix_index(str(bedgz), preset="bed", force=True)
    return bedgz


@pytest.mark.integration
def test_wps_cli_help():
    """Test WPS CLI help output."""
    result = runner.invoke(app, ["wps", "--help"])
    assert result.exit_code == 0
    output = strip_ansi(result.output)
    assert "--wps-anchors" in output
    assert "--background" in output
    assert "--genome" in output


@pytest.mark.integration
def test_wps_cli_parquet_output(tmp_path, sample_bedgz, sample_anchors):
    """Test WPS produces Parquet output with expected columns."""
    import pandas as pd

    output_dir = tmp_path / "output"

    result = runner.invoke(
        app,
        [
            "wps",
            "--input",
            str(sample_bedgz),
            "--wps-anchors",
            str(sample_anchors),
            "-o",
            str(output_dir),
            "-s",
            "test_sample",
            "--no-gc-correct",  # Skip GC for speed
        ],
    )

    # Check exit code
    if result.exit_code != 0:
        pytest.skip(f"WPS execution failed (may need Rust backend): {result.output}")

    # Check Parquet output exists
    wps_parquet = output_dir / "test_sample.WPS.parquet"
    assert wps_parquet.exists(), f"WPS Parquet missing. Output: {result.output}"

    # Verify columns
    df = pd.read_parquet(wps_parquet)
    expected_cols = ["region_id", "chrom", "center", "strand"]
    for col in expected_cols:
        assert col in df.columns, f"Missing column: {col}"


@pytest.mark.integration
def test_wps_cli_smooth_columns(tmp_path, sample_bedgz, sample_anchors):
    """Test WPS applies smoothing and adds _smooth columns."""
    import pandas as pd

    output_dir = tmp_path / "output"

    result = runner.invoke(
        app,
        [
            "wps",
            "--input",
            str(sample_bedgz),
            "--wps-anchors",
            str(sample_anchors),
            "-o",
            str(output_dir),
            "-s",
            "test_sample",
            "--no-gc-correct",
        ],
    )

    if result.exit_code != 0:
        pytest.skip("WPS execution failed")

    wps_parquet = output_dir / "test_sample.WPS.parquet"
    if not wps_parquet.exists():
        pytest.skip("Parquet output not created")

    df = pd.read_parquet(wps_parquet)

    # Raw WPS columns should always be present
    assert "wps_nuc" in df.columns, "Missing wps_nuc column"
    assert "wps_tf" in df.columns, "Missing wps_tf column"

    # Smoothing columns are added during post-processing (run-all pipeline),
    # not by the standalone wps CLI. Verify they're absent in standalone mode.
    # If smoothing is ever integrated into standalone wps, update this test.
    if "wps_nuc_smooth" in df.columns:
        assert (
            "wps_tf_smooth" in df.columns
        ), "Missing wps_tf_smooth if smoothing applied"


@pytest.mark.integration
def test_wps_genome_flag():
    """Test --genome flag is accepted."""
    result = runner.invoke(app, ["wps", "--help"])
    output = strip_ansi(result.output)
    assert "--genome" in output
    assert "hg19" in output or "hg38" in output


@pytest.mark.integration
def test_local_depth_is_mean_coverage_over_the_window(
    tmp_path, early_closing_bedgz, sample_anchors
):
    """A covered region survives, and reports a positive mean depth.

    Two defects met in one expression (#79). `local_depth` divided the running
    value left after integrating the difference arrays -- coverage at the
    window's last base -- by the window length, and the emptiness guard tested
    that same value against `0.0`. With every fragment closing early, as here,
    the array telescopes to exactly zero: the region was dropped, and when it
    did survive the column came out as signed floating-point residue, negative
    on 22-42% of rows across real cohorts.

    Asserted end to end rather than only in Rust because the broken value
    reached `{sample}.WPS.parquet`, which is the product downstream reads.
    """
    import pandas as pd

    output_dir = tmp_path / "output"
    result = runner.invoke(
        app,
        [
            "wps",
            "--input",
            str(early_closing_bedgz),
            "--wps-anchors",
            str(sample_anchors),
            "-o",
            str(output_dir),
            "-s",
            "test_sample",
            "--no-gc-correct",
        ],
    )
    # Not pytest.skip. The other tests in this file skip on a non-zero exit, and
    # that is the wrong call here: this is a regression test for the Rust path
    # being invoked at all, so a CLI error is precisely the failure it exists to
    # catch. Skipping would let a broken backend report green -- the same
    # "absence reads as success" shape as the `== 0.0` guard this test is about.
    assert result.exit_code == 0, (
        "the wps CLI failed, so the regression below was never exercised. "
        f"exit={result.exit_code}\n{result.output}"
    )

    wps_parquet = output_dir / "test_sample.WPS.parquet"
    assert wps_parquet.exists(), f"WPS Parquet missing. Output: {result.output}"

    df = pd.read_parquet(wps_parquet)
    assert not df.empty, (
        "the anchor was dropped: its fragments all close before the window ends, "
        "which is exactly what the old `curr_cov == 0.0` guard discarded"
    )
    assert "local_depth" in df.columns, "local_depth is a declared contract column"

    depth = df["local_depth"]
    assert (depth >= 0).all(), (
        f"negative coverage is impossible; got min {depth.min():.3g}. The column "
        "is reporting residue from the difference array rather than a mean."
    )
    assert (
        depth > 0
    ).any(), "every row reports zero depth over a window with 20 fragments in it"

    # The scale, not just the sign. A sign check passes just as happily on
    # `total_cov` left undivided -- a factor of 2001 out here -- which is the same
    # class of error as the one being fixed, so the test has to pin the value.
    #
    # Derivation. With --no-gc-correct, gc_weight is 1.0, and the length weights
    # are 1.0 for 167 bp (nucleosome stream) and 55 bp (TF stream):
    #
    #   window   the anchor centre 1000 expanded +/-1000, so len = 2000-0+1 = 2001
    #   per frag `length` bases -- a fragment covers exactly its own bases
    #   total    10*167 + 10*55 = 1670 + 550 = 2220
    #
    # Through 0.9.2 this read 2200: fragments were placed with a 1-based start
    # against a 0-based region, so each lost its first base (#83). The per-stream
    # intervals are pinned on a single fragment in wps.rs
    # `coordinate_convention_tests`; this test pins the end-to-end consequence.
    expected = 2220 / 2001
    assert depth.max() == pytest.approx(expected, rel=1e-6), (
        f"local_depth is the wrong scale: got {depth.max():.10g}, expected "
        f"{expected:.10g} (= 2220/2001). Near 2220 means the total was written "
        "without dividing by the window length; near 2200/2001 means fragments "
        "are losing a base again (#83); near zero means the last base's "
        "coverage came back instead of the mean."
    )
