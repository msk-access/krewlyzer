"""Run the real pipeline over a synthetic sample.

Shared by the invariant tests so each one states its assertion and nothing
else. Deliberately drives the same ``run_features`` the CLI and ``run-all``
both use -- a test against a reimplementation would prove nothing about what
ships.
"""

from __future__ import annotations

import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import pysam

from krewlyzer import _core
from krewlyzer.core.unified_processor import FeatureOutputs, run_features

from .synth import SynthProfile, SynthSample, make_sample


@dataclass
class PipelineRun:
    sample: SynthSample
    bed_gz: Path
    outputs: FeatureOutputs
    output_dir: Path


def extract_to_bedgz(sample: SynthSample, directory: Path) -> Path:
    """Extract fragments and leave a tabix-indexed BED where features can read it."""
    directory.mkdir(parents=True, exist_ok=True)
    raw = directory / "fragments.bed"

    _core.extract_motif.process_bam_parallel(
        str(sample.bam),
        str(sample.reference),
        20,  # mapq
        65,  # min_len
        1000,  # max_len
        4,  # kmer
        1,  # threads
        str(raw),
        None,  # motif prefix
        None,  # exclude
        None,  # target regions
        True,  # skip duplicates
        True,  # require proper pair
        True,  # silent
    )

    # The Rust writer emits BGZF whatever the extension says, so compressing
    # again would produce a double-wrapped file that tabix rejects.
    bed_gz = directory / "fragments.bed.gz"
    with open(raw, "rb") as fh:
        already_compressed = fh.read(2) == b"\x1f\x8b"
    if already_compressed:
        shutil.move(str(raw), str(bed_gz))
    else:
        pysam.tabix_compress(str(raw), str(bed_gz), force=True)
    pysam.tabix_index(str(bed_gz), preset="bed", force=True)
    return bed_gz


def run_pipeline(
    directory: Path,
    profile: SynthProfile,
    sample: Optional[SynthSample] = None,
    **feature_flags,
) -> PipelineRun:
    """Synthesize a sample, extract it, and run the requested features."""
    directory.mkdir(parents=True, exist_ok=True)
    sample = sample or make_sample(directory / "input", profile)
    bed_gz = extract_to_bedgz(sample, directory / "extract")
    out_dir = directory / "features"

    flags = {"enable_fsc": True, "enable_fsr": True, "enable_fsd": True}
    flags.update(feature_flags)

    outputs = run_features(
        bed_path=bed_gz,
        output_dir=out_dir,
        sample_name=profile.name,
        fsc_bins=sample.assets.bins,
        fsd_arms=sample.assets.arms,
        gc_correct=False,  # no GC model for a synthetic reference
        threads=1,
        **flags,
    )
    return PipelineRun(
        sample=sample, bed_gz=bed_gz, outputs=outputs, output_dir=out_dir
    )
