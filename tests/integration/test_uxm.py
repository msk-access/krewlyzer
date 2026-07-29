import pysam
from krewlyzer.cli import app
from typer.testing import CliRunner


def create_mock_bam_uxm(path):
    header = {"HD": {"VN": "1.0"}, "SQ": [{"LN": 5000, "SN": "chr1"}]}

    with pysam.AlignmentFile(str(path), "wb", header=header) as outf:
        # Read 1: Methylated (M) - XM: ZZZ (all methylated)
        a = pysam.AlignedSegment()
        a.query_name = "read_M"
        a.query_sequence = "CGATA"
        a.flag = 0
        a.reference_id = 0
        a.reference_start = 1000
        a.mapping_quality = 60
        a.cigar = ((0, 5),)
        a.set_tag("XM", "ZZZ")  # 3 CpGs methylated
        outf.write(a)

        # Read 2: Unmethylated (U) - XM: zzz
        b = pysam.AlignedSegment()
        b.query_name = "read_U"
        b.query_sequence = "CGATA"
        b.flag = 0
        b.reference_id = 0
        b.reference_start = 1000
        b.mapping_quality = 60
        b.cigar = ((0, 5),)
        b.set_tag("XM", "zzz")  # 3 CpGs unmethylated
        outf.write(b)

    pysam.sort("-o", str(path), str(path))
    pysam.index(str(path))


def test_uxm_integration(tmp_path):
    # Single BAM file (not directory)
    bam_file = tmp_path / "sample.bam"
    create_mock_bam_uxm(bam_file)

    mark_file = tmp_path / "markers.bed"
    with open(mark_file, "w") as f:
        f.write("chr1\t1000\t1005\n")

    output_dir = tmp_path / "output"

    runner = CliRunner()
    # New CLI: uxm -i <bam> -o <output_dir> -s <sample_name> -m <markers>
    result = runner.invoke(
        app,
        [
            "uxm",
            "-i",
            str(bam_file),
            "-o",
            str(output_dir),
            "-s",
            "sample",
            "-m",
            str(mark_file),
        ],
    )

    assert result.exit_code == 0, f"CLI failed: {result.output}"

    # Check output
    output_file = output_dir / "sample.UXM.tsv"
    assert output_file.exists()


def create_mock_bam_uxm_mixed(path):
    """BAM with one M, one U and one genuinely MIXED fragment.

    The mixed read has 2/4 CpGs methylated (ratio 0.50), which must be
    classified X under the Loyfer thresholds (U <= 0.25, M >= 0.75).
    """
    header = {"HD": {"VN": "1.0"}, "SQ": [{"LN": 5000, "SN": "chr1"}]}

    with pysam.AlignmentFile(str(path), "wb", header=header) as outf:
        for name, xm in [("read_M", "ZZZZ"), ("read_U", "zzzz"), ("read_X", "ZZzz")]:
            a = pysam.AlignedSegment()
            a.query_name = name
            a.query_sequence = "CGATA"
            a.flag = 0
            a.reference_id = 0
            a.reference_start = 1000
            a.mapping_quality = 60
            a.cigar = ((0, 5),)
            a.set_tag("XM", xm)
            outf.write(a)

    pysam.sort("-o", str(path), str(path))
    pysam.index(str(path))


def test_uxm_mixed_fragments_are_classified_x(tmp_path):
    """Regression: the X (mixed-methylation) class must be reachable.

    Historically both thresholds were passed as 0.5, and because the backend
    tests `ratio >= methy_threshold` first, every fragment collapsed into
    M or U -- the published X column was identically zero for every sample.
    """
    bam_file = tmp_path / "sample.bam"
    create_mock_bam_uxm_mixed(bam_file)

    mark_file = tmp_path / "markers.bed"
    with open(mark_file, "w") as f:
        f.write("chr1\t1000\t1005\n")

    output_dir = tmp_path / "output"

    runner = CliRunner()
    result = runner.invoke(
        app,
        [
            "uxm",
            "-i",
            str(bam_file),
            "-o",
            str(output_dir),
            "-s",
            "sample",
            "-m",
            str(mark_file),
        ],
    )
    assert result.exit_code == 0, f"CLI failed: {result.output}"

    output_file = output_dir / "sample.UXM.tsv"
    assert output_file.exists()

    rows = [ln.strip().split("\t") for ln in output_file.read_text().splitlines()]
    header, data = rows[0], rows[1]
    record = dict(zip(header, data))

    u, x, m = float(record["U"]), float(record["X"]), float(record["M"])

    # One fragment in each class -> each fraction is 1/3.
    assert x > 0.0, f"X class is unreachable (U={u}, X={x}, M={m})"
    assert u > 0.0 and m > 0.0, f"expected all three classes populated: {record}"
    # Output is rounded to 6 dp, so allow a little slack on the sum.
    assert abs(u + x + m - 1.0) < 1e-5, f"fractions must sum to 1: {record}"
