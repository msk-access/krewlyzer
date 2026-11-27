"""A panel donor, extracted the way `extract` and the PON build extract it.

Shared by the tests of what each extraction path writes beside the BED: the GC factor tables
(`test_factor_table_writes.py`) and the record a kept `build-pon` cache carries
(`test_pon_cache_record.py`).
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional, Tuple

from .runall import SAMPLE
from .synth import Contig, SynthProfile, SynthSample, make_sample


def dense_panel(directory: Path) -> Tuple[SynthSample, Path]:
    """A panel whose targets hold all but a sliver of the fragments: extraction's off-target tally
    is too small to fit, while the whole BED (on-target fragments included) fits a table.
    """
    profile = SynthProfile(
        "dense", n_fragments=60_000, length_mixture=((167, 2, 1.0),), seed=7
    )
    contigs = (Contig("chr1", 300_000), Contig("chr2", 300_000))
    dense = make_sample(directory / "in", profile, contigs=contigs)
    targets = directory / "targets.bed"
    targets.write_text("chr1\t0\t300000\nchr2\t0\t298500\n")
    return dense, targets


def factor_tables(out: Path, *populations: str) -> List[str]:
    """The names of the factor tables in ``out`` for ``populations``, every format."""
    from krewlyzer.core.gc_assets import FACTOR_KINDS, factor_table_base

    names = {
        f"{factor_table_base(out, SAMPLE, population).name}{ext}"
        for population in populations
        for _, ext in FACTOR_KINDS
    }
    return sorted(p.name for p in out.iterdir() if p.name in names)


def extract(
    sample: SynthSample,
    out: Path,
    targets: Path,
    *extra: str,
    output_format: str = "parquet",
):
    """`krewlyzer extract` on the sample's BAM. Returns the CLI result; ``extra`` are more options."""
    from typer.testing import CliRunner

    from krewlyzer.cli import app

    return CliRunner().invoke(
        app,
        [
            "extract",
            "-i",
            str(sample.bam),
            "-r",
            str(sample.reference),
            "-o",
            str(out),
            "-s",
            SAMPLE,
            "-x",
            str(sample.assets.blacklist),
            "-T",
            str(targets),
            "-F",
            output_format,
            *extra,
        ],
    )


def build_donor(
    sample: SynthSample,
    out: Path,
    targets: Path,
    *,
    input_path: Optional[Path] = None,
    arms: Optional[Path] = None,
) -> None:
    """`build-pon`'s per-donor path (`process_sample`): extraction (from the BAM unless
    `input_path` is a BED) and its factor tables, every feature off unless `arms` turns FSD on (and
    with it the PON's GC gate). With the sample's synthetic blacklist, as run-all and `extract` are
    given it. ``out`` may already exist and hold an earlier run's files.
    """
    from krewlyzer.core.sample_processor import SampleParams, process_sample

    out.mkdir(parents=True, exist_ok=True)
    process_sample(
        input_path=input_path or sample.bam,
        output_dir=out,
        sample_name=SAMPLE,
        reference=sample.reference,
        params=SampleParams(
            threads=1, exclude_regions=sample.assets.blacklist, arms_file=arms
        ),
        target_regions=targets,
        enable_fsc=False,
        enable_fsr=False,
        enable_fsd=arms is not None,
        enable_wps=False,
        enable_ocf=False,
        enable_tfbs=False,
        enable_atac=False,
        pon_mode=True,
        output_format="tsv",
        compress=False,
    )
