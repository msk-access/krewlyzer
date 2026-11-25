"""
GC Correction Asset Resolution Helper.

Provides a centralized function for resolving GC correction assets
to avoid code duplication across standalone tools.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, TYPE_CHECKING
import logging

from .output_provenance import bed_as_extracted

if TYPE_CHECKING:
    from ..assets import AssetManager

logger = logging.getLogger("krewlyzer.core.gc_assets")


#: Share of fragments on contigs the reference lacks above which it is warned about. Extract itself
#: refuses a wrong reference (a primary contig absent, or a length mismatch) before reading; what
#: remains is decoys and alts, recorded in the metadata whatever the share.
CONTIG_MISMATCH_WARN_SHARE = 0.01


def check_reference_contigs(n_missing: int, n_fetch_failed: int, n_total: int) -> float:
    """Report extract's contig counts; return the share of fragments on contigs the reference lacks.

    Those fragments, and any whose reference fetch failed, have no GC and enter no GC model. A
    wrong reference never reaches here: extract refuses it from the BAM header.
    """
    share = n_missing / n_total if n_total else 0.0
    if share > CONTIG_MISMATCH_WARN_SHARE:
        logger.warning(
            f"{share:.1%} of fragments are on contigs missing from the reference "
            "(decoys, alts, unplaced); their GC is unresolved and they enter no GC model"
        )
    if n_fetch_failed:
        logger.warning(
            f"{n_fetch_failed} fragments' reference sequence could not be fetched; "
            "their GC is unresolved and they enter no GC model"
        )
    return share


#: The populations a GC factor table is fitted on. Panel runs fit both.
POPULATIONS = ("off_target", "on_target")
#: Search order for an existing table: full-precision Parquet first, then compressed TSV, then TSV.
FACTOR_KINDS = (("parquet", ".parquet"), ("tsv.gz", ".tsv.gz"), ("tsv", ".tsv"))


def factor_table_base(directory: Path, sample: str, population: str) -> Path:
    """`{directory}/{sample}.correction_factors[.ontarget]`, without an extension.

    The only place a factor-table filename is spelled: writers add the extension their
    `--output-format`/`--compress` call for, and readers go through `resolve_gc_factors`.
    """
    if population not in POPULATIONS:
        raise ValueError(
            f"unknown GC factor population {population!r}; expected one of {POPULATIONS}"
        )
    suffix = ".correction_factors" + ("" if population == "off_target" else ".ontarget")
    return Path(directory) / f"{sample}{suffix}"


def remove_factor_tables(
    directory: Path, sample: str, populations: tuple = POPULATIONS
) -> list:
    """Remove a sample's factor tables, in every format, before a run writes its own.

    A table belongs to the BED it was fitted on. The resolver finds one in any format, Parquet
    first, so a table an earlier run left here was reused whenever this run wrote none (no
    fittable bin, no observations, a failed fit) or wrote another format, and the metadata then
    recorded the correction as applied. Returns the paths removed.
    """
    removed = []
    for population in populations:
        base = factor_table_base(directory, sample, population)
        for _, ext in FACTOR_KINDS:
            path = Path(f"{base}{ext}")
            if path.exists():
                path.unlink()
                removed.append(path)
    if removed:
        logger.info(
            f"Removed {len(removed)} GC factor table(s) left by an earlier run: "
            + ", ".join(p.name for p in removed)
        )
    return removed


def gc_populations_needed(panel_mode: bool) -> tuple:
    """The factor tables a sample must have to count as GC-corrected: panel runs fit both."""
    return POPULATIONS if panel_mode else ("off_target",)


def factor_table_path(
    directory: Path, sample: str, population: str, ext: str = ".tsv"
) -> Path:
    """Path handed to a factor *writer*. The Rust writer replaces the extension with the one the
    run's `--output-format`/`--compress` call for, so `.tsv` here is a base, not a format choice.
    """
    base = factor_table_base(directory, sample, population)
    return base.with_name(base.name + ext)


@dataclass(frozen=True)
class GcFactorsChoice:
    """Which GC factor table a run uses for one population, or why there is none."""

    path: Optional[Path]
    kind: str  # "parquet" | "tsv.gz" | "tsv" | "recomputed" | "none"
    reason: str = (
        ""  # why there is no table; never a filename (filenames carry the sample ID)
    )
    searched: tuple = ()  # the names looked for, for log messages only

    @property
    def found(self) -> bool:
        return self.path is not None

    def record(self) -> str:
        """The metadata value: the kind, never a filename (filenames carry the sample ID)."""
        return self.kind if self.found else f"none:{self.reason}"


def resolve_gc_factors(
    bed_path: Path, output_dir: Optional[Path], sample: str, population: str
) -> GcFactorsChoice:
    """Find the GC factor table for `population`, in any format run-all writes.

    Search order: beside the BED, under the BED's stem (the table extract wrote with it), then
    `output_dir` under the sample name. Within each, Parquet, then `.tsv.gz`, then `.tsv`.
    Pass `output_dir=None` to search beside the BED only.

    Every reader goes through this. Before it existed, three call sites looked only for a plain
    `.tsv`, which production (`--output-format both --compress`) never writes, and silently ran
    uncorrected.

    Beside the BED means where extraction wrote it, as the BED check places it
    (`output_provenance.bed_as_extracted`): a link is followed to its end, and the table is the
    one beside the last file along it that an extraction claims, so a file of another name no
    extraction claims (an annex or DVC object, or a versioned file beside the link) is read with
    the table beside the link. In a
    Nextflow task, a BED the pipeline extracted is linked in from that extraction's folder, with
    its tables; a samplesheet BED staged as a link is read with the table beside its source.
    """
    bed_path = bed_as_extracted(bed_path)
    places = [(bed_path.parent, bed_path.name.removesuffix(".gz").removesuffix(".bed"))]
    if output_dir is not None:
        places.append((Path(output_dir), sample))
    searched = []
    for directory, stem in places:
        base = factor_table_base(directory, stem, population)
        for kind, ext in FACTOR_KINDS:
            candidate = base.with_name(base.name + ext)
            if candidate.exists():
                return GcFactorsChoice(candidate, kind)
        searched.append(base.name + ".{parquet,tsv.gz,tsv}")
    where = (
        "beside the BED or in the output directory"
        if output_dir is not None
        else "beside the BED"
    )
    return GcFactorsChoice(
        None, "none", f"not found {where}", tuple(dict.fromkeys(searched))
    )


def resolve_feature_gc_tables(
    bed_path: Path,
    output_dir: Path,
    sample: str,
    gc_correct: bool,
    panel_mode: bool,
    refitted: bool,
) -> "tuple[GcFactorsChoice, GcFactorsChoice]":
    """The (off-target, on-target) factor tables run_features' Python-side consumers use.

    `refitted` is True when no table sat beside the BED, so the Rust pass refitted factors from
    the BED itself (Phase 1, a different tally from extract's); that choice is recorded as
    `recomputed` and warned about, because its values differ from run-all's.
    """
    if not gc_correct:
        off = on = GcFactorsChoice(None, "none", "gc correction disabled")
        return off, on
    off = resolve_gc_factors(bed_path, output_dir, sample, "off_target")
    on = (
        resolve_gc_factors(bed_path, output_dir, sample, "on_target")
        if panel_mode
        else GcFactorsChoice(None, "none", "not panel mode")
    )
    # Found must mean usable. The consumers of these tables (gene FSC, panel WPS, TFBS, ATAC) each
    # catch their own failure and carry on, so an unreadable table that merely existed made every
    # one of them fail quietly while the metadata recorded the correction as applied.
    from krewlyzer import _core

    for choice in (off, on):
        if choice.path is not None:
            try:
                _core.gc.validate_factor_table(str(choice.path))
            except Exception as e:
                label = "off-target" if choice is off else "on-target"
                raise RuntimeError(
                    f"{label} GC factor table {choice.path.name} cannot be read ({e}); "
                    "refusing to run rather than record a correction that did not happen"
                ) from e
    if off.found and refitted:
        off = GcFactorsChoice(off.path, "recomputed")
        logger.warning(
            "GC factors were refitted from this BED because extract's table was not beside it; "
            "values differ from run-all. Run the tool on the extract output directory to match."
        )
    return off, on


def require_gc_for_pon(choices: dict, panel_mode: bool) -> None:
    """Refuse to build a PON baseline from a donor whose GC factors were not found.

    The donor's features would enter the baseline uncorrected and corrected samples would be
    scored against them -- the mismatch item 0 removes. A missing donor is better than that.
    """
    for population in gc_populations_needed(panel_mode):
        choice = choices.get(population)
        if choice is None or not choice.found:
            reason = choice.reason if choice is not None else "not resolved"
            raise RuntimeError(
                f"PON donor has no {population.replace('_', '-')} GC factor table ({reason}); "
                "refusing to build a baseline from uncorrected features"
            )


@dataclass
class GcAssets:
    """Resolved GC correction assets."""

    gc_ref: Optional[Path] = None
    valid_regions: Optional[Path] = None
    factors_out: Optional[Path] = None
    factors_input: Optional[Path] = None
    gc_correct_enabled: bool = False


def resolve_gc_assets(
    assets: "AssetManager",
    output_dir: Path,
    sample_name: str,
    bedgz_input: Path,
    gc_correct: bool = True,
    genome: str = "hg19",
) -> GcAssets:
    """
    Resolve GC correction assets with automatic pre-computed factors detection.

    This centralizes the GC asset resolution logic used by all standalone tools
    (fsc, fsd, fsr, wps, ocf, etc.)

    Args:
        assets: AssetManager instance
        output_dir: Output directory for correction_factors.tsv
        sample_name: Sample name for output file naming
        bedgz_input: Input BED.gz file (checked for adjacent correction_factors.tsv)
        gc_correct: Whether GC correction is requested
        genome: Genome build for logging

    Returns:
        GcAssets dataclass with resolved paths and enabled status
    """
    result = GcAssets()

    if not gc_correct:
        return result

    # Pre-computed factors beside the input BED (extract writes them there). Only beside the BED:
    # a table an earlier single-tool run left in output_dir must not be silently reused.
    pre = resolve_gc_factors(bedgz_input, None, sample_name, "off_target")
    if pre.path is not None:
        result.factors_input = pre.path
        result.gc_correct_enabled = True
        logger.info(f"Using pre-computed GC factors: {pre.path.name}")
        return result

    # Resolve bundled GC assets
    try:
        result.gc_ref = assets.resolve("gc_reference")
        result.valid_regions = assets.resolve("valid_regions")
        result.factors_out = factor_table_path(output_dir, sample_name, "off_target")
        result.gc_correct_enabled = True
        logger.info(f"GC correction enabled ({genome})")
    except FileNotFoundError as e:
        logger.warning(f"GC correction assets not found: {e}")
        logger.warning(
            "Proceeding without GC correction. Use --no-gc-correct to suppress."
        )
        result.gc_correct_enabled = False

    return result
