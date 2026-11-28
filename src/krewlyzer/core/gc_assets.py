"""
GC factor tables: where a sample's table is, or why there is none.

Only extraction fits a table, and writes it beside the BED it fitted it from. A feature step
reads that table through `resolve_gc_factors`, or, with none, runs without GC correction.
run_features (run-all and the unified-pipeline tools) and `region-entropy` find and check their
tables with `resolve_feature_gc_tables` and say which they use and why one is missing with
`report_gc_tables`; run-all's mFSD takes run_features' choice. Extraction records the choice
(`GcFactorsChoice.record`) by the same call, before any feature is computed, and a single tool
whose choice differs stops (`output_provenance.refuse_a_different_gc_choice`). The single `mfsd`
tool takes its table by name (`--correction-factors`).
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Optional
import logging

from .output_provenance import bed_as_extracted, extraction_wrote_none_beside

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
#: Why a table beside the BED is not used: its extraction recorded writing none.
NOT_EXTRACTIONS = "beside the BED but not extraction's, whose record says it wrote none"


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


def remove_factor_tables(directory: Path, sample: str) -> list:
    """Remove a sample's factor tables, in every format, before a run writes its own.

    A table belongs to the BED it was fitted on. The resolver finds one in any format, Parquet
    first, so a table an earlier run left here was reused whenever this run wrote none (no
    fittable bin, no observations, a failed fit) or wrote another format, and the metadata then
    recorded the correction as applied. Returns the paths removed.
    """
    removed = []
    for population in POPULATIONS:
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
    kind: str  # "parquet" | "tsv.gz" | "tsv" | "none"
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


def resolve_gc_factors(bed_path: Path, population: str) -> GcFactorsChoice:
    """Find the GC factor table for `population`, in any format run-all writes.

    Only beside the BED, under the BED's stem: the table extraction fitted from those fragments
    and wrote with them. Parquet, then `.tsv.gz`, then `.tsv`. A table elsewhere, such as in the
    folder a tool writes to, may be another extraction's, fitted from other fragments, so it is
    never used.

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

    An off-target table beside a record in which extraction says it wrote none is not
    extraction's (`output_provenance.extraction_wrote_none_beside`): builds before the refit's
    removal fitted one there, from the whole BED, when run-all's feature step found none, and a
    rerun reusing that BED took it for extraction's. No step fits an on-target table but
    extraction, and the record says nothing of it.
    """
    bed_path = bed_as_extracted(bed_path)
    stem = bed_path.name.removesuffix(".gz").removesuffix(".bed")
    base = factor_table_base(bed_path.parent, stem, population)
    searched = (base.name + ".{parquet,tsv.gz,tsv}",)
    for kind, ext in FACTOR_KINDS:
        candidate = base.with_name(base.name + ext)
        if candidate.exists():
            if population == "off_target" and extraction_wrote_none_beside(
                bed_path.parent, stem, ext
            ):
                return GcFactorsChoice(None, "none", NOT_EXTRACTIONS, searched)
            return GcFactorsChoice(candidate, kind)
    return GcFactorsChoice(None, "none", "not found beside the BED", searched)


def resolve_feature_gc_tables(
    bed_path: Path, gc_correct: bool, panel_mode: bool
) -> "tuple[GcFactorsChoice, GcFactorsChoice]":
    """The (off-target, on-target) factor tables a sample's features correct with: run_features'
    Rust pass and its Python-side consumers (gene FSC, panel WPS, TFBS, ATAC), the single tool
    `region-entropy`, and run-all's mFSD, which takes run_features' choice.
    """
    if not gc_correct:
        off = on = GcFactorsChoice(None, "none", "gc correction disabled")
        return off, on
    off = resolve_gc_factors(bed_path, "off_target")
    on = (
        resolve_gc_factors(bed_path, "on_target")
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
    return off, on


def report_gc_tables(
    off: GcFactorsChoice,
    on: GcFactorsChoice,
    *,
    gc_correct: bool,
    panel_mode: bool,
    off_effect: str,
    on_effect: str,
    log: logging.Logger,
) -> None:
    """Log the tables `resolve_feature_gc_tables` chose, and warn of a missing one: why it is
    missing, what runs without it (``off_effect``, ``on_effect``: the caller's features), and
    whether re-extracting can give one.

    One text for every consumer. region-entropy had its own, which said only that TFBS and ATAC
    ran uncorrected, nothing of when re-extracting helps, and nothing at all of a missing
    on-target table.
    """
    if gc_correct and not off.found:
        log.warning(
            f"GC factors (off-target) {off.reason} (looked for {', '.join(off.searched)}). "
            f"{off_effect} Only extraction fits one, and its log says why it wrote none. With "
            "too few off-target fragments to fit (a panel whose targets hold almost every "
            "fragment) there is nothing to fit, and re-extracting changes nothing. If it ran "
            "with --no-gc-correct, or its table was not kept beside the BED, re-extract with GC "
            "correction and keep the table beside the BED."
        )
    if gc_correct and panel_mode and not on.found:
        log.warning(
            f"GC factors (on-target) {on.reason} (looked for {', '.join(on.searched)}); "
            f"{on_effect} With too few on-target fragments to fit there is nothing to fit. If "
            "extraction ran with --no-gc-correct or without these target regions, or its "
            "tables were not kept beside the BED, re-extract with GC correction and the target "
            "regions, and keep its tables beside the BED."
        )
    for label, choice in (("off-target", off), ("on-target", on)):
        if choice.path is not None:
            log.info(f"  GC factors ({label}): {choice.path.name} [{choice.kind}]")


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
