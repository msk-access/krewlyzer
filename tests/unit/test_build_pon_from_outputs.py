"""``build-pon --from-outputs`` aggregates a directory nobody has to re-extract.

The in-process builder re-runs feature extraction itself, 55-97 minutes per
sample and roughly 15 hours for a 47-sample cohort, to produce inputs
``run-all`` already writes. Every defect found in the 0.9.0 models so far has
been in aggregation rather than extraction, so each one cost a full rebuild to
re-check. This route costs minutes.

The fixture writes ``.parquet`` and ``.tsv.gz`` and *no* plain ``.tsv``, because
that is what a real ``run-all`` directory looks like. An earlier version of the
collector hard-coded ``.tsv`` and found nothing in the one layout it exists to
read; the two Rust readers that parse plain TSV silently returned nothing for
the same reason.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from krewlyzer.pon import from_outputs

# Enough samples to clear MIN_SAMPLES_PER_KEY with one to spare.
N_SAMPLES = 4
GENES = ("TP53", "EGFR", "KRAS")
KMERS = ["".join(p) for p in __import__("itertools").product("ACTG", repeat=4)]


def _write_pair(frame: pd.DataFrame, base: Path) -> None:
    """As ``run-all`` writes: parquet plus gzipped TSV, never a bare ``.tsv``."""
    frame.to_parquet(Path(f"{base}.parquet"))
    with gzip.open(Path(f"{base}.tsv.gz"), "wt") as handle:
        frame.to_csv(handle, sep="\t", index=False)


def _sample_dir(root: Path, index: int, *, band_limited: bool = False) -> Path:
    """One run-all-shaped sample directory whose values vary with ``index``.

    Each entry is drawn independently rather than shifted by a per-sample
    constant. A constant offset gives every gene and exon the *identical*
    sigma, which `validate-pon` correctly reports as `PON.BLOCK_DEGENERATE` --
    the signature of a fabricated baseline. The first version of this fixture
    did exactly that and the gate caught it, which is the gate working.
    """
    rng = np.random.default_rng(index)
    stem = f"DONOR{index}-T"
    directory = root / stem
    directory.mkdir(parents=True)
    base = directory / stem

    def spread(centre: float, scale: float, n: int) -> np.ndarray:
        """Values whose cross-sample spread differs per entry."""
        return centre + rng.normal(0.0, scale, n) * np.linspace(0.5, 2.0, n)

    bins = [f"{s}-{s + 4}" for s in range(65, 105, 5)]
    _write_pair(
        pd.DataFrame(
            {
                "mean_gc": [0.38, 0.42, 0.46, 0.50],
                "ultra_short": rng.integers(50, 90, 4),
                "core_short": rng.integers(200, 400, 4),
                "mono_nucl": rng.integers(800, 1200, 4),
                "di_nucl": rng.integers(100, 200, 4),
                "long": rng.integers(30, 70, 4),
            }
        ),
        Path(f"{base}.fsc_counts"),
    )
    _write_pair(
        pd.DataFrame(
            [
                {"arm": arm, **{b: float(rng.integers(10, 500)) for b in bins}}
                for arm in ("1p", "1q", "2p")
            ]
        ),
        Path(f"{base}.FSD"),
    )

    # WPS goes to a Rust *parquet* reader, so it is parquet-only in run-all too.
    #
    # The list type must be float32. Written as float64 the reader finds no
    # usable column, skips every anchor, and reports "backed by fewer than 3
    # samples" -- which reads as a cohort that is too small rather than a
    # schema that does not match. Met while writing this fixture.
    import pyarrow as pa
    import pyarrow.parquet as pq

    pq.write_table(
        pa.Table.from_pandas(
            pd.DataFrame(
                {
                    "region_id": [f"TSS|{g}|1" for g in GENES],
                    # Per anchor as well as per sample. One shared waveform
                    # gives every anchor the identical amplitude and shape
                    # correlation, which the gate reports as
                    # PON.BLOCK_DEGENERATE on wps_shape_baseline -- correctly.
                    "wps_nuc": [
                        list(
                            (1.0 + 0.3 * g + 0.05 * index)
                            * np.sin(np.arange(200) / (10.0 + g) + index + g)
                        )
                        for g in range(len(GENES))
                    ],
                    "wps_tf": [
                        list(
                            (0.8 + 0.2 * g)
                            * np.sin(np.arange(200) / (5.0 + g) + index * 0.7)
                        )
                        for g in range(len(GENES))
                    ],
                }
            ),
            schema=pa.schema(
                [
                    ("region_id", pa.string()),
                    ("wps_nuc", pa.list_(pa.float32())),
                    ("wps_tf", pa.list_(pa.float32())),
                ]
            ),
        ),
        Path(f"{base}.WPS.parquet"),
    )

    pd.DataFrame(
        {
            "group_id": ["Global_All", "Family_AluY", "Chr8_All"],
            "nrl_bp": [181.0 + index * 1.7, 170.0 + index * 0.9, 250.0],
            "periodicity_score": list(spread(0.45, 0.02, 3)),
            # Chr8_All sits at the search-band edge for every sample, which is
            # what the real xs2 cohort does for four of its groups.
            "nrl_at_band_limit": [False, False, True if band_limited else False],
        }
    ).to_parquet(Path(f"{base}.WPS_background.parquet"))

    counts = rng.integers(500, 1500, len(KMERS)).astype(float)
    freq = counts / counts.sum()
    _write_pair(
        pd.DataFrame({"Motif": KMERS, "Frequency": np.round(freq, 6)}),
        Path(f"{base}.EndMotif"),
    )
    entropy = -float(np.sum(freq * np.log2(freq + 1e-12))) / np.log2(len(KMERS))
    _write_pair(
        pd.DataFrame([{"Sample": stem, "MDS": round(entropy, 6)}]), Path(f"{base}.MDS")
    )

    _write_pair(
        pd.DataFrame(
            [
                {"gene": g, "mds_mean": m, "mds_e1": e}
                for g, m, e in zip(
                    GENES,
                    spread(0.84, 0.01, len(GENES)),
                    spread(0.80, 0.01, len(GENES)),
                )
            ]
        ),
        Path(f"{base}.MDS.gene"),
    )
    _write_pair(
        pd.DataFrame(
            [
                {"gene": g, "name": f"{g}_exon{e}", "mds": v}
                for (g, e), v in zip(
                    [(g, e) for g in GENES for e in (1, 2)],
                    spread(0.83, 0.01, len(GENES) * 2),
                )
            ]
        ),
        Path(f"{base}.MDS.exon"),
    )
    _write_pair(
        pd.DataFrame(
            {"region_id": ["Lymph", "Liver"], "ocf": list(spread(10.0, 1.5, 2))}
        ),
        Path(f"{base}.OCF"),
    )
    for feature in ("TFBS", "ATAC"):
        _write_pair(
            pd.DataFrame(
                {"label": ["CTCF", "AR"], "entropy": list(spread(0.68, 0.03, 2))}
            ),
            Path(f"{base}.{feature}"),
        )
    _write_pair(
        pd.DataFrame(
            [
                {"gene": g, "normalized_depth": d}
                for g, d in zip(GENES, spread(1.0, 0.05, len(GENES)))
            ]
        ),
        Path(f"{base}.FSC.gene"),
    )
    # On-target counterparts.
    #
    # An assay like xs1 resolves bundled target regions, so a run-all directory
    # for it is panel mode and carries every one of these. Omitting them left
    # the built model short of nine on-target blocks -- which the gate's
    # packing-list check reported, correctly, as an incomplete PON.
    _write_pair(
        pd.DataFrame(
            {
                "mean_gc": [0.38, 0.42, 0.46, 0.50],
                "ultra_short": rng.integers(40, 80, 4),
                "core_short": rng.integers(150, 350, 4),
                "mono_nucl": rng.integers(600, 1000, 4),
                "di_nucl": rng.integers(80, 180, 4),
                "long": rng.integers(20, 60, 4),
            }
        ),
        Path(f"{base}.fsc_counts.ontarget"),
    )
    _write_pair(
        pd.DataFrame(
            [
                {"arm": arm, **{b: float(rng.integers(10, 400)) for b in bins}}
                for arm in ("1p", "1q", "2p")
            ]
        ),
        Path(f"{base}.FSD.ontarget"),
    )
    pq.write_table(
        pa.Table.from_pandas(
            pd.DataFrame(
                {
                    "region_id": [f"PANEL|{g}|1" for g in GENES],
                    "wps_nuc": [
                        list(
                            (1.1 + 0.25 * g + 0.04 * index)
                            * np.sin(np.arange(200) / (11.0 + g) + index + g)
                        )
                        for g in range(len(GENES))
                    ],
                    "wps_tf": [
                        list(
                            (0.7 + 0.15 * g)
                            * np.sin(np.arange(200) / (4.0 + g) + index * 0.6)
                        )
                        for g in range(len(GENES))
                    ],
                }
            ),
            schema=pa.schema(
                [
                    ("region_id", pa.string()),
                    ("wps_nuc", pa.list_(pa.float32())),
                    ("wps_tf", pa.list_(pa.float32())),
                ]
            ),
        ),
        Path(f"{base}.WPS.panel.parquet"),
    )
    for suffix, centre in ((".OCF.ontarget", 11.0), (".OCF.offtarget", 9.0)):
        _write_pair(
            pd.DataFrame(
                {
                    "region_id": ["Lymph", "Liver"],
                    "ocf": list(spread(centre, 1.5, 2)),
                }
            ),
            Path(f"{base}{suffix}"),
        )
    for feature in ("TFBS", "ATAC"):
        _write_pair(
            pd.DataFrame(
                {"label": ["CTCF", "AR"], "entropy": list(spread(0.66, 0.03, 2))}
            ),
            Path(f"{base}.{feature}.ontarget"),
        )
    on_counts = rng.integers(400, 1200, len(KMERS)).astype(float)
    on_freq = on_counts / on_counts.sum()
    _write_pair(
        pd.DataFrame({"Motif": KMERS, "Frequency": np.round(on_freq, 6)}),
        Path(f"{base}.EndMotif.ontarget"),
    )
    on_entropy = -float(np.sum(on_freq * np.log2(on_freq + 1e-12))) / np.log2(
        len(KMERS)
    )
    _write_pair(
        pd.DataFrame([{"Sample": stem, "MDS": round(on_entropy, 6)}]),
        Path(f"{base}.MDS.ontarget"),
    )

    _write_pair(
        pd.DataFrame(
            [
                {
                    "chrom": "17",
                    "start": 7_000_000 + e * 1000,
                    "end": 7_000_500 + e * 1000,
                    "normalized_depth": d,
                }
                for e, d in enumerate(spread(1.0, 0.05, 3))
            ]
        ),
        Path(f"{base}.FSC.regions"),
    )
    return directory


@pytest.fixture
def corpus(tmp_path):
    root = tmp_path / "runall"
    for i in range(N_SAMPLES):
        _sample_dir(root, i, band_limited=True)
    return root


def test_every_sample_is_complete_and_collected(corpus):
    dirs = from_outputs.discover_samples(corpus)
    assert len(dirs) == N_SAMPLES
    assert all(from_outputs.incomplete_reason(d) is None for d in dirs)

    collected, skipped = from_outputs.collect(dirs)
    assert not skipped
    assert len(collected.sample_ids) == N_SAMPLES
    # Every input, for every sample. A short one means a suffix went stale.
    for name in (
        "gc",
        "fsd",
        "ocf",
        "mds",
        "tfbs",
        "atac",
        "fsc_gene",
        "fsc_region",
        "fsd_paths",
        "wps_paths",
        "wps_background_paths",
        "mds_gene_paths",
        "mds_exon_paths",
    ):
        assert len(getattr(collected, name)) == N_SAMPLES, f"{name} is short"


def test_no_bare_tsv_exists_in_the_fixture(corpus):
    """The premise. If this fails the fixture stopped resembling run-all."""
    assert not list(corpus.glob("*/*.tsv")), "fixture has plain TSV; run-all does not"
    assert list(corpus.glob("*/*.tsv.gz")) and list(corpus.glob("*/*.parquet"))


def test_a_half_written_directory_is_refused_by_name(corpus):
    """Silently aggregating it makes the cohort smaller than its metadata."""
    victim = sorted(corpus.iterdir())[0]
    for leftover in victim.glob("*.EndMotif.*"):
        leftover.unlink()

    reason = from_outputs.incomplete_reason(victim)
    assert reason is not None and "EndMotif" in reason

    collected, skipped = from_outputs.collect(from_outputs.discover_samples(corpus))
    assert len(skipped) == 1 and skipped[0][0] == victim
    assert len(collected.sample_ids) == N_SAMPLES - 1


def test_a_run_that_never_finished_is_refused(corpus):
    """The marker is the *last* file written, so its absence means unfinished."""
    victim = sorted(corpus.iterdir())[1]
    for leftover in victim.glob("*.MDS.gene.*"):
        leftover.unlink()
    assert "did not finish" in (from_outputs.incomplete_reason(victim) or "")


def test_the_motif_record_comes_back_off_disk(corpus):
    """The one input the in-process path keeps only in memory."""
    collected, _ = from_outputs.collect(from_outputs.discover_samples(corpus))
    record = collected.mds[0]
    assert len(record["kmers"]) == 256
    assert 0.0 < record["mds"] < 1.0
    assert sum(record["kmers"].values()) == pytest.approx(1.0, abs=1e-3)


def test_the_mds_score_is_recomputed_when_its_table_is_absent(corpus):
    """``compute_mds`` normalises, so frequencies and counts give one answer.

    Exact, not an approximation -- which is why a missing MDS table is worth
    recovering from rather than refusing.
    """
    victim = sorted(corpus.iterdir())[0]
    stem = from_outputs.sample_stem(victim)
    written = from_outputs._motif_record(victim, stem)
    for leftover in victim.glob(f"{stem}.MDS.parquet"):
        leftover.unlink()
    for leftover in victim.glob(f"{stem}.MDS.tsv.gz"):
        leftover.unlink()
    recomputed = from_outputs._motif_record(victim, stem)

    assert recomputed is not None
    assert recomputed["mds"] == pytest.approx(written["mds"], abs=1e-5)


# ---------------------------------------------------------------------------
# the whole command
# ---------------------------------------------------------------------------


def _build(corpus: Path, output: Path):
    """Always panel mode, on an explicit target BED.

    Without one, panel mode depends on whether the assay's bundled targets
    happen to be present -- which they are locally and are not in CI, where
    the payload is LFS-backed. The same test then meant two different things
    in two places, and the difference only surfaced as a CI-only failure.

    The BED is written beside the corpus and contains one region. The panel
    baselines are aggregated from the fixture's on-target tables, so nothing
    reads its coordinates; it exists to pin the mode.
    """
    from typer.testing import CliRunner

    from krewlyzer.cli import app

    targets = corpus.parent / "targets.bed"
    if not targets.exists():
        targets.write_text("17\t7000000\t7010000\tTP53\n")

    return CliRunner().invoke(
        app,
        [
            "build-pon",
            "--from-outputs",
            str(corpus),
            "--assay",
            "xs1",
            "-r",
            "/dev/null",
            "-o",
            str(output),
            "--target-regions",
            str(targets),
            "--cohort-label",
            "synthetic",
        ],
    )


def test_a_model_is_built_and_passes_its_own_gate(corpus, tmp_path):
    """End to end, and then through `validate-pon`.

    A build that produces a model the gate rejects has not succeeded, whatever
    its exit code said.
    """
    from krewlyzer.validate.pon_gate import check_pon, exit_code

    output = tmp_path / "synthetic.pon.parquet"
    result = _build(corpus, output)
    assert result.exit_code == 0, result.output
    assert output.exists()

    findings = check_pon(output)
    assert exit_code(findings) == 0, [f"{f.id} {f.table}.{f.column}" for f in findings]


def test_the_band_limited_group_survives_with_no_nrl(corpus, tmp_path):
    """Part A, reached through the whole command rather than the builder."""
    output = tmp_path / "synthetic.pon.parquet"
    assert _build(corpus, output).exit_code == 0

    table = pd.read_parquet(output)
    block = table[table["table"] == "wps_background"]
    row = block[block["group_id"] == "Chr8_All"]
    assert len(row) == 1, "the row was dropped rather than emptied"
    assert np.isnan(row.iloc[0]["nrl_mean"]), "the band edge was reported as a mean"
    assert row.iloc[0]["n_at_band_limit"] == N_SAMPLES
    # ...while the groups that measured one keep theirs.
    measured = block[block["group_id"] == "Global_All"].iloc[0]
    assert np.isfinite(measured["nrl_mean"]) and measured["nrl_mean"] < 250.0


def test_giving_both_inputs_is_refused(corpus, tmp_path):
    """Which cohort would the digest then describe?"""
    from typer.testing import CliRunner

    from krewlyzer.cli import app

    listing = tmp_path / "samples.txt"
    listing.write_text("/nonexistent.bam\n")
    result = CliRunner().invoke(
        app,
        [
            "build-pon",
            str(listing),
            "--from-outputs",
            str(corpus),
            "--assay",
            "xs1",
            "-r",
            "/dev/null",
            "-o",
            str(tmp_path / "x.pon.parquet"),
        ],
    )
    assert result.exit_code == 2


def test_an_incomplete_cohort_stops_the_build_unless_allowed(corpus, tmp_path):
    """A dropped sample must be a decision, not a side effect."""
    victim = sorted(corpus.iterdir())[0]
    for leftover in victim.glob("*.EndMotif.*"):
        leftover.unlink()

    assert _build(corpus, tmp_path / "a.pon.parquet").exit_code == 1


def test_aggregation_needs_no_extraction_assets(corpus, tmp_path):
    """No BAM is opened, so no reference, bin file or anchor set is needed.

    This failed in CI while passing locally: the bin file is LFS-backed, the
    runner has no LFS payload, and `build-pon` refused before reaching the
    branch that would never have read it. An aggregation route that depends on
    a data package it does not touch is not much of an aggregation route.
    """
    from typer.testing import CliRunner

    from krewlyzer.cli import app

    output = tmp_path / "no-assets.pon.parquet"
    result = CliRunner().invoke(
        app,
        [
            "build-pon",
            "--from-outputs",
            str(corpus),
            "--assay",
            "xs1",
            "-r",
            "/dev/null",
            "-o",
            str(output),
            "--bin-file",
            str(tmp_path / "definitely-not-here.bed.gz"),
        ],
    )
    assert result.exit_code == 0, result.output
    assert output.exists()
