"""Drive the real `run-all` over a synthetic sample.

`run_all` is a typer command, so calling it as a function leaves every option
it is not given as a typer `OptionInfo`. That is exactly how run-all's own UXM
call broke: it passed no `threads`, and `OptionInfo > 0` raised `TypeError` on
every run. So this fills in every option from the command's own declared
defaults, and a test overrides only what it is about.
"""

from __future__ import annotations

import inspect
from pathlib import Path
from typing import Dict, Optional

import pandas as pd
from typer.models import OptionInfo

from krewlyzer.wrapper import run_all

from .synth import SynthSample

SAMPLE = "S"


def _declared_defaults() -> Dict[str, object]:
    defaults = {}
    for name, param in inspect.signature(run_all).parameters.items():
        if isinstance(param.default, OptionInfo) and param.default.default is not ...:
            defaults[name] = param.default.default
    return defaults


def write_arms(path: Path, sample: SynthSample) -> Path:
    """Arms named the way the validator insists on (`1p`, not `chr1p`)."""
    with open(path, "w") as fh:
        for contig in sample.contigs:
            mid = contig.length // 2
            arm = contig.name.removeprefix("chr")
            fh.write(f"{contig.name}\t0\t{mid}\t{arm}p\n")
            fh.write(f"{contig.name}\t{mid}\t{contig.length}\t{arm}q\n")
    return path


def write_vcf(path: Path, n_variants: int = 1) -> Path:
    lines = ["##fileformat=VCFv4.2", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]
    for i in range(n_variants):
        lines.append(f"chr1\t{1100 + 997 * i}\t.\tA\tT\t.\t.\t.")
    path.write_text("\n".join(lines) + "\n")
    return path


def synth_run_all(sample: SynthSample, output: Path, **overrides) -> Optional[int]:
    """Run run-all into ``output``. Returns the exit code if it exited, else None.

    Parquet, no PON, no on/off-target split and no TFBS/ATAC unless a test asks:
    the bundled panel assets are on hg19 coordinates and the synthetic contigs are
    300 kb, so those stages would only add time.
    """
    import typer

    arms = write_arms(sample.directory / "arms.validated.bed", sample)
    options = _declared_defaults()
    options.update(
        bam_input=sample.bam,
        reference=sample.reference,
        output=output,
        sample_name=SAMPLE,
        bin_input=sample.assets.bins,
        arms_file=arms,
        ocr_file=sample.assets.regions,
        wps_anchors=sample.assets.anchors,
        wps_background=sample.assets.alu,
        output_format="parquet",
        skip_pon=True,
        skip_target_regions=True,
        no_tfbs=True,
        no_atac=True,
        threads=1,
    )
    options.update(overrides)
    try:
        run_all(**options)
    except (SystemExit, typer.Exit) as exc:
        return getattr(exc, "code", None) or getattr(exc, "exit_code", None)
    return None


def record(output: Path, sample: str = SAMPLE) -> Dict[str, object]:
    """The sample's metadata row."""
    return pd.read_parquet(output / f"{sample}.metadata.parquet").iloc[0].to_dict()


def outcomes(output: Path, sample: str = SAMPLE) -> Dict[str, tuple]:
    """``{stage: (status, error_class)}`` as the sample's record has them."""
    from krewlyzer.validate.contract import STAGES

    row = record(output, sample)
    found = {}
    for stage in STAGES:
        status = row.get(f"{stage}_status")
        if isinstance(status, str) and status:
            error = row.get(f"{stage}_error")
            found[stage] = (status, error if isinstance(error, str) else "")
    return found


def lfs_pointer(path: Path) -> Path:
    """What a bundled asset is in a clone that skipped `git lfs pull`: text, not gzip.

    The realistic way for an asset-reading stage to fail, and a real failing input rather than a
    patched exception.
    """
    path.write_text(
        "version https://git-lfs.github.com/spec/v1\n"
        "oid sha256:" + "0" * 64 + "\nsize 1024\n"
    )
    return path
