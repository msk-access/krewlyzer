"""Build a PON by aggregating ``run-all`` output directories.

``build-pon SAMPLE_LIST`` re-runs feature extraction itself, four samples at a
time on one node: 55-97 minutes per sample, roughly 15 hours for a 47-sample
cohort. Meanwhile ``run-all`` computes every one of those features already, and
a ``run-all`` output directory is a strict superset of what the PON builder
consumes -- verified against the 0.8.3 corpus, all 18 inputs present, 83 files
against the PON path's 31.

So the expensive half is redundant. This module is the cheap half on its own: a
pure aggregation over per-sample directories, reading only files. Nothing here
opens a BAM.

That buys three things the previous plan needed separate machinery for:

* **Re-aggregating after a fix costs minutes.** Every defect found in the
  0.9.0 models so far -- the fabricated ``wps_background``, the sigma floors,
  the band-limited NRLs -- was in aggregation, not extraction. Each one
  previously meant a 15-hour rebuild to re-check.
* **Any existing ``run-all`` output can seed a PON**, including the 0.8.3
  corpus already on the share.
* **Per-sample cluster parallelism** comes free, because generating the
  directories is now someone else's job -- an sbatch array, Nextflow, or the
  in-process path itself.

The in-process path stays. It is the only way to build a PON without a
cluster, and `scripts/build_pon.sh` still uses it.

One difference is unavoidable and worth stating plainly. ``.EndMotif.tsv``
stores ``Frequency`` rounded to six decimals, while the in-process path passes
unrounded ``kmer_frequencies`` through memory. Any file-based route inherits
that rounding. At a typical k-mer frequency of 0.0039 and a cross-donor sigma
of 2.7e-4, one rounding step is 0.37% of sigma: the two routes agree on sigma
to within 0.08% and on z to within about 0.002. Pinned by
``test_the_two_routes_agree``. Since Parquet is the product and the product is
what rounds, this is the precision that actually exists.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd

from krewlyzer.core.output_utils import TABLE_EXTENSIONS, read_table, resolve_table_path

logger = logging.getLogger("build-pon")

#: The output whose presence means a sample finished.
#:
#: Measured from the write order: ``bed.gz`` -> WPS/FSC/FSD/OCF -> TFBS ->
#: ATAC -> MDS last. Checking the *last* output rather than counting files is
#: what stops a half-written directory being aggregated as though it were
#: complete -- the failure that makes a cohort quietly smaller than it claims.
#:
#: Held without an extension, and resolved through `resolve_table_path`.
#: Writers honour ``--output-format`` and ``--compress``, so the same logical
#: output lands as ``.tsv``, ``.tsv.gz`` or ``.parquet``. A real ``run-all``
#: directory carries ``.parquet`` and ``.tsv.gz`` and *no* bare ``.tsv``, so a
#: hard-coded ``.MDS.gene.tsv`` would have found nothing in the one layout this
#: module exists to read.
COMPLETION_MARKER = ".MDS.gene"

#: Files an aggregation cannot proceed without, beyond the marker.
#:
#: Deliberately short. A sample missing an optional on-target table is a
#: panel-mode question, not a broken sample, and refusing it would make this
#: route unusable on WGS output.
REQUIRED_SUFFIXES = (
    ".fsc_counts",
    ".FSD",
    ".WPS",
    ".WPS_background",
    ".EndMotif",
)


@dataclass
class CollectedInputs:
    """Everything the baseline builders read, and nothing else.

    One field per collector in ``build_pon``'s aggregation stage, named to
    match so the two routes can be compared by eye. Both routes fill this;
    only the filling differs.
    """

    # Per-sample tables passed to Rust by path
    fsd_paths: List[str] = field(default_factory=list)
    fsd_ontarget_paths: List[str] = field(default_factory=list)
    wps_paths: List[str] = field(default_factory=list)
    wps_panel_paths: List[str] = field(default_factory=list)
    wps_background_paths: List[str] = field(default_factory=list)
    mds_gene_paths: List[str] = field(default_factory=list)
    mds_exon_paths: List[str] = field(default_factory=list)

    # Per-sample frames and dicts aggregated in Python
    gc: List[dict] = field(default_factory=list)
    gc_ontarget: List[dict] = field(default_factory=list)
    fsd: List[pd.DataFrame] = field(default_factory=list)
    fsd_ontarget: List[pd.DataFrame] = field(default_factory=list)
    ocf: List[pd.DataFrame] = field(default_factory=list)
    ocf_ontarget: List[pd.DataFrame] = field(default_factory=list)
    ocf_offtarget: List[pd.DataFrame] = field(default_factory=list)
    mds: List[dict] = field(default_factory=list)
    mds_ontarget: List[dict] = field(default_factory=list)
    #: Breakpoint 4-mers, which are a different distribution from end 4-mers
    #: and need their own baseline -- see `core/motif_pon.py`.
    breakpoint_motif: List[dict] = field(default_factory=list)
    breakpoint_motif_ontarget: List[dict] = field(default_factory=list)
    tfbs: List[dict] = field(default_factory=list)
    tfbs_ontarget: List[dict] = field(default_factory=list)
    atac: List[dict] = field(default_factory=list)
    atac_ontarget: List[dict] = field(default_factory=list)
    fsc_gene: List[dict] = field(default_factory=list)
    fsc_region: List[dict] = field(default_factory=list)

    #: Sample stems, in directory order. Feeds the cohort digest, so it must
    #: hold what was actually aggregated rather than what was asked for.
    sample_ids: List[str] = field(default_factory=list)


def discover_samples(directory: Path) -> List[Path]:
    """Per-sample subdirectories of ``directory``, sorted.

    A ``run-all`` cohort is one directory per sample; so is the
    ``--keep-sample-outputs`` cache the in-process path writes. Sorted so a
    rebuild aggregates in the same order and the cohort digest is stable.
    """
    if not directory.is_dir():
        raise NotADirectoryError(f"not a directory: {directory}")
    return sorted(p for p in directory.iterdir() if p.is_dir())


def sample_stem(sample_dir: Path) -> Optional[str]:
    """The ``{stem}`` that this directory's files are named for.

    Taken from the marker file rather than the directory name: the two agree
    for ``run-all`` output but need not in general, and the files are what the
    readers below will open.
    """
    for extension in TABLE_EXTENSIONS:
        hits = sorted(sample_dir.glob(f"*{COMPLETION_MARKER}{extension}"))
        if hits:
            return hits[0].name[: -len(f"{COMPLETION_MARKER}{extension}")]
    return None


def incomplete_reason(sample_dir: Path) -> Optional[str]:
    """Why this directory cannot be aggregated, or None when it can.

    Returns a reason rather than a bool so the refusal can name the file. A
    build that drops a sample silently is how a cohort ends up smaller than
    its own metadata claims -- the same class as the downstream reader that
    swallows exceptions and yields an empty feature dict.
    """
    stem = sample_stem(sample_dir)
    if stem is None:
        return f"no {COMPLETION_MARKER} -- the run did not finish"
    missing = [
        s
        for s in REQUIRED_SUFFIXES
        if resolve_table_path(sample_dir / f"{stem}{s}") is None
    ]
    if missing:
        return f"missing {', '.join(missing)}"
    return None


def _read(base: Path) -> Optional[pd.DataFrame]:
    """Read whatever was written for ``base``, whichever extension it took.

    ``base`` carries no extension. `read_table` resolves it and is
    parquet-first, which is the right preference here: this module reads
    directories someone else finished writing, so "whatever was produced for
    this output" is exactly the question. (Immediately *after* a write the
    answer differs -- that is what `read_exact_table` is for.)

    Returns None with a warning rather than raising. A warning rather than a
    debug line: the output resolved, so failing to parse it is a fact about
    the cohort, not noise.
    """
    path = resolve_table_path(base)
    if path is None:
        return None
    try:
        return read_table(path)
    except Exception as exc:
        logger.warning(f"  could not read {path.name}: {type(exc).__name__}: {exc}")
        return None


def _gc_record(frame: pd.DataFrame) -> Optional[dict]:
    """The GC-bias input, in the same shape the in-process path builds.

    Kept identical to `build_pon`'s inline version deliberately -- the two
    routes must produce the same model, and a second spelling of the same
    grouping is how they would drift.
    """
    needed = ["ultra_short", "core_short", "mono_nucl", "di_nucl", "long", "mean_gc"]
    if not all(c in frame.columns for c in needed):
        return None
    return {
        "gc": frame["mean_gc"].values,
        "short": (frame["ultra_short"] + frame["core_short"]).values,
        "intermediate": frame["mono_nucl"].values,
        "long": (frame["di_nucl"] + frame["long"]).values,
    }


def _ocf_record(frame: pd.DataFrame) -> Optional[pd.DataFrame]:
    """OCF as ``region_id``/``ocf``, accepting the older ``tissue``/``OCF``."""
    if "tissue" in frame.columns and "OCF" in frame.columns:
        frame = frame.rename(columns={"tissue": "region_id", "OCF": "ocf"})
    if "region_id" in frame.columns and "ocf" in frame.columns:
        return frame[["region_id", "ocf"]]
    return None


def _motif_record(
    sample_dir: Path,
    stem: str,
    suffix: str = "",
    table: str = "EndMotif",
) -> Optional[dict]:
    """K-mer frequencies and the MDS score, from the files rather than memory.

    This is the one input the in-process path never wrote: ``process_sample``
    does not call ``write_motif_outputs``, so a ``--keep-sample-outputs`` cache
    has no ``EndMotif`` table even though the counts reached the model through
    memory. A ``run-all`` directory has it.

    The MDS score is read when present and recomputed otherwise. Recomputing is
    exact, not an approximation: ``compute_mds`` normalises its input, so
    frequencies and counts give the same entropy.
    """
    from krewlyzer.core.motif_processor import compute_mds

    frame = _read(sample_dir / f"{stem}.{table}{suffix}")
    if frame is None or not {"Motif", "Frequency"}.issubset(frame.columns):
        return None
    kmers = dict(zip(frame["Motif"], frame["Frequency"]))

    # MDS is defined on *end* motifs, so it is not computed for breakpoints:
    # a number there would be a different statistic wearing the same name.
    if table != "EndMotif":
        return {"kmers": kmers, "mds": None}

    score: Optional[float] = None
    mds_frame = _read(sample_dir / f"{stem}.MDS{suffix}")
    if mds_frame is not None and "MDS" in mds_frame.columns and len(mds_frame):
        score = float(mds_frame["MDS"].iloc[0])
    if score is None:
        score = float(compute_mds(kmers))
        logger.debug(f"  {stem}: no MDS{suffix} table, recomputed from EndMotif")

    return {"kmers": kmers, "mds": score}


def _entropy_record(base: Path) -> Optional[dict]:
    """``{label: entropy}``, via the same two functions ``run-all`` uses.

    ``process_sample`` builds ``tfbs_data`` by reading the file it just wrote,
    so calling the same pair here is not a reimplementation -- it is the same
    code on the same input, and cannot drift from it.
    """
    from krewlyzer.core.region_entropy_processor import (
        extract_entropy_data,
        load_entropy_tsv,
    )

    path = resolve_table_path(base)
    if path is None:
        return None
    try:
        # `load_entropy_tsv` handles the parquet case itself; going through it
        # rather than `_read` keeps this identical to what `process_sample`
        # does with the very same file.
        return extract_entropy_data(load_entropy_tsv(path))
    except Exception as exc:
        logger.warning(f"  could not read {path.name}: {type(exc).__name__}: {exc}")
        return None


def _depth_map(base: Path, key: str) -> Optional[Dict[str, float]]:
    """``{key: normalized_depth}`` for the FSC gene and region tables.

    ``key`` is ``"gene"`` for the gene table; the region table has no single
    id column, so its key is built from coordinates exactly as the in-process
    path builds it.
    """
    frame = _read(base)
    if frame is None or "normalized_depth" not in frame.columns:
        return None
    if key == "gene":
        if "gene" not in frame.columns:
            return None
        ids = frame["gene"]
    else:
        if not {"chrom", "start", "end"}.issubset(frame.columns):
            return None
        ids = (
            frame["chrom"].astype(str)
            + ":"
            + frame["start"].astype(str)
            + "-"
            + frame["end"].astype(str)
        )
    return dict(zip(ids, frame["normalized_depth"]))


def collect(sample_dirs: List[Path]) -> Tuple[CollectedInputs, List[Tuple[Path, str]]]:
    """Read every baseline input from a list of per-sample directories.

    Returns the collected inputs and the directories that were refused, each
    with the reason. Callers decide whether a refusal is fatal; nothing is
    dropped without being reported.
    """
    collected = CollectedInputs()
    skipped: List[Tuple[Path, str]] = []

    for sample_dir in sample_dirs:
        reason = incomplete_reason(sample_dir)
        if reason:
            skipped.append((sample_dir, reason))
            continue

        stem = sample_stem(sample_dir)
        assert stem is not None  # incomplete_reason would have caught it
        collected.sample_ids.append(stem)

        def path(suffix: str) -> Path:
            return sample_dir / f"{stem}{suffix}"

        # --- tables Rust reads by path -----------------------------------
        path_targets: Tuple[Tuple[str, List[str]], ...] = (
            (".FSD", collected.fsd_paths),
            (".FSD.ontarget", collected.fsd_ontarget_paths),
            (".WPS", collected.wps_paths),
            (".WPS.panel", collected.wps_panel_paths),
            (".WPS_background", collected.wps_background_paths),
            (".MDS.gene", collected.mds_gene_paths),
            (".MDS.exon", collected.mds_exon_paths),
        )
        # Paths handed on verbatim. WPS goes to a Rust parquet reader and FSD
        # to a Rust plain-TSV one; `_compute_fsd_baseline` normalises whatever
        # arrives, because a `run-all` directory holds `.FSD.parquet` and
        # `.FSD.tsv.gz` and no plain `.tsv` at all.
        for suffix, path_bucket in path_targets:
            resolved = resolve_table_path(path(suffix))
            if resolved is not None:
                path_bucket.append(str(resolved))

        # --- GC bias ------------------------------------------------------
        gc_targets: Tuple[Tuple[str, List[dict]], ...] = (
            (".fsc_counts", collected.gc),
            (".fsc_counts.ontarget", collected.gc_ontarget),
        )
        for suffix, gc_bucket in gc_targets:
            frame = _read(path(suffix))
            gc_record = _gc_record(frame) if frame is not None else None
            if gc_record is not None:
                gc_bucket.append(gc_record)

        # --- FSD frames (the Python-side aggregate) -----------------------
        fsd_targets: Tuple[Tuple[str, List[pd.DataFrame]], ...] = (
            (".FSD", collected.fsd),
            (".FSD.ontarget", collected.fsd_ontarget),
        )
        for suffix, frames in fsd_targets:
            frame = _read(path(suffix))
            if frame is not None:
                frames.append(frame)

        # --- OCF ----------------------------------------------------------
        ocf_targets: Tuple[Tuple[str, List[pd.DataFrame]], ...] = (
            (".OCF", collected.ocf),
            (".OCF.ontarget", collected.ocf_ontarget),
            (".OCF.offtarget", collected.ocf_offtarget),
        )
        for suffix, frames in ocf_targets:
            frame = _read(path(suffix))
            ocf_frame = _ocf_record(frame) if frame is not None else None
            if ocf_frame is not None:
                frames.append(ocf_frame)

        # --- motifs -------------------------------------------------------
        motif_targets: Tuple[Tuple[str, str, List[dict]], ...] = (
            ("EndMotif", "", collected.mds),
            ("EndMotif", ".ontarget", collected.mds_ontarget),
            ("BreakPointMotif", "", collected.breakpoint_motif),
            ("BreakPointMotif", ".ontarget", collected.breakpoint_motif_ontarget),
        )
        for table, suffix, motif_bucket in motif_targets:
            motif = _motif_record(sample_dir, stem, suffix, table=table)
            if motif is not None:
                motif_bucket.append(motif)

        # --- region entropy -------------------------------------------------
        entropy_targets: Tuple[Tuple[str, List[dict]], ...] = (
            (".TFBS", collected.tfbs),
            (".TFBS.ontarget", collected.tfbs_ontarget),
            (".ATAC", collected.atac),
            (".ATAC.ontarget", collected.atac_ontarget),
        )
        for suffix, entropy_bucket in entropy_targets:
            entropy = _entropy_record(path(suffix))
            if entropy:
                entropy_bucket.append(entropy)

        # --- FSC depth maps -------------------------------------------------
        gene_map = _depth_map(path(".FSC.gene"), "gene")
        if gene_map:
            collected.fsc_gene.append(gene_map)
        region_map = _depth_map(path(".FSC.regions"), "region")
        if region_map:
            collected.fsc_region.append(region_map)

    return collected, skipped


def describe(collected: CollectedInputs) -> str:
    """A one-line-per-input summary, for the build log.

    Printed before aggregation so a cohort that is short of one input is
    visible up front rather than inferred afterwards from a missing block.
    """
    counts: List[Any] = [
        ("gc", len(collected.gc)),
        ("fsd", len(collected.fsd_paths)),
        ("wps", len(collected.wps_paths)),
        ("wps_background", len(collected.wps_background_paths)),
        ("ocf", len(collected.ocf)),
        ("mds", len(collected.mds)),
        ("mds_gene", len(collected.mds_gene_paths)),
        ("mds_exon", len(collected.mds_exon_paths)),
        ("tfbs", len(collected.tfbs)),
        ("atac", len(collected.atac)),
        ("fsc_gene", len(collected.fsc_gene)),
        ("fsc_region", len(collected.fsc_region)),
    ]
    n = len(collected.sample_ids)
    short = [f"{name} {got}/{n}" for name, got in counts if got != n]
    full = ", ".join(f"{name} {got}" for name, got in counts)
    line = f"collected from {n} samples: {full}"
    if short:
        line += f"\n  short of the cohort: {', '.join(short)}"
    return line
