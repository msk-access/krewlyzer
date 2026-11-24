"""Record, in the product itself, what produced a sample and what it was scored against.

`{sample}.metadata.parquet` is the downstream consumer's **completion marker**:
a sample without it is dropped from the cohort silently (invariant #2). It
recorded the genome, the assay, the filters and whether GC correction ran — and
nothing at all about the code version or the PON.

The version *was* recorded, in `{sample}.features.json`. Downstream reads
Parquet only and never that file, so the one fact needed to tell two runs apart
lived precisely where the product ignores it.

That became load-bearing in 0.9.0. This release changed what several columns
*mean* — every FSD `_logR` was previously scored against the wrong size bin,
`pon_stability` was wrong by a median of 4709%, and WPS gained seven columns.
Every one of those defects produced output that was present, finite and
plausible. Given a results directory, there was no way to answer *was this
produced before or after the fix?* except by re-running it.

Five fields answer it:

``krewlyzer_version``      which code wrote this sample
``pon_applied``            whether any feature was z-scored at all
``pon_model``              which model, by **basename** — never the full path,
                           which would pin the operator's home directory into a
                           shipped table
``pon_cohort_digest``      which healthy cohort that model was fitted to; the
                           salted digest from `pon/provenance.py`, so it answers
                           "the same cohort?" without naming anyone
``pon_krewlyzer_version``  the release the model was stamped with

`pon_applied` is the one a reader acts on. A `--skip-pon` run and a run with no
PON are both legitimately unscored, and look identical to a reader counting
columns. (A PON that was asked for and cannot be used stops the run instead:
`pon_integration.UnusablePonError`.) Recording
the distinction is also what lets `validate-output` require the PON-derived
columns when — and only when — a PON was actually used.

The same problem one level down is the optional stages run-all carries on past (UXM, mFSD,
Region-MDS, gene FSC, panel WPS, TFBS, ATAC, and the MDS z-scores): "not requested" and
"requested and failed" both leave a table absent. `record_stage_outcomes` writes which, as
`{stage}_status` and `{stage}_error`.
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Mapping, Optional, Sequence, Tuple

from .output_utils import (
    read_exact_table,
    TABLE_EXTENSIONS,
    resolve_table_path,
    table_extensions,
    table_formats,
    table_paths,
    write_table,
)

logger = logging.getLogger("core.output_provenance")

#: The columns this module adds. Named once so the gate and the tests can both
#: refer to them without restating the list (invariant #5, one level down).
PROVENANCE_COLUMNS = (
    "krewlyzer_version",
    "pon_applied",
    "pon_model",
    "pon_cohort_digest",
    "pon_digest",
    "pon_krewlyzer_version",
)


def provenance_fields(
    pon: Any = None,
    pon_path: Optional[Path] = None,
) -> Dict[str, Any]:
    """The provenance to record for one sample.

    ``pon`` is the loaded model and ``pon_path`` the file it came from; pass
    neither for an unscored run. Both are required for ``pon_applied`` to be
    True: a path without a loaded model is not a scored run.
    """
    from krewlyzer import __version__

    applied = pon is not None and pon_path is not None
    return {
        "krewlyzer_version": __version__,
        "pon_applied": applied,
        # Basename only. The full path would put the operator's home directory
        # into a table that ships, which `scripts/check_phi_guard.sh` exists to
        # keep out of this project.
        "pon_model": Path(pon_path).name if (applied and pon_path) else "",
        "pon_cohort_digest": (
            str(getattr(pon, "cohort_digest", "") or "") if applied else ""
        ),
        # What the PON holds, not who it was fitted to: a PON rebuilt from the same donors keeps
        # its name and cohort digest and changes this.
        "pon_digest": str(getattr(pon, "pon_digest", "") or "") if applied else "",
        "pon_krewlyzer_version": (
            str(getattr(pon, "krewlyzer_version", "") or "") if applied else ""
        ),
    }


def record_path(
    output_dir: Path, sample: str, *, output_format: str, compress: bool
) -> Optional[Path]:
    """The metadata file a run of this format reads its record from: its own, else whatever exists.

    One rule for the stamp and the second-PON check, so the two never read different copies (a
    TSV run beside an earlier run's `.parquet`, or the reverse).
    """
    base = output_dir / f"{sample}.metadata"
    own = table_paths(base, output_format=output_format, compress=compress)[0]
    return own if own.exists() else resolve_table_path(base)


def _flag(value: Any) -> bool:
    """A recorded boolean, with an empty cell false: `bool(NaN)` is True, so an empty
    `pon_applied` read back from TSV counted as applied."""
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return False
    return bool(value)


def _as_text(value: Any) -> str:
    """A recorded field as text, with None and NaN both empty: an empty field reads back from
    TSV as NaN, which `str()` turns into "nan" and made the same PON look like another.
    """
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return ""
    return str(value)


def refuse_a_second_pon(
    output_dir: Path,
    sample: str,
    pon: Any,
    pon_path: Path,
    *,
    output_format: str,
    compress: bool,
) -> None:
    """Refuse to score a sample against a PON other than the one its record names.

    `pon_applied` is one record per sample. A single tool scoring against another PON would
    leave tables from two PONs under a record that names one, and nothing would flag it.
    """
    from .pon_integration import UnusablePonError

    existing = record_path(
        output_dir, sample, output_format=output_format, compress=compress
    )
    frame = _read_record(existing)
    if frame is None or frame.empty or "pon_applied" not in frame.columns:
        return
    if not _flag(frame["pon_applied"].iloc[0]):
        return
    ours = provenance_fields(pon, pon_path)
    recorded = {
        key: frame[key].iloc[0] if key in frame.columns else None
        for key in ("pon_model", "pon_cohort_digest", "pon_digest")
    }
    # The digest when both sides have one; a record written before it existed falls back to the
    # name and cohort digest, so a single tool still runs in a directory an older run-all scored.
    if _as_text(recorded["pon_digest"]) and ours["pon_digest"]:
        compared: tuple = ("pon_digest",)
    else:
        compared = ("pon_model", "pon_cohort_digest")
    if any(_as_text(recorded[k]) != _as_text(ours[k]) for k in compared):
        raise UnusablePonError(
            f"{sample} was scored against {recorded['pon_model']} "
            f"(cohort {recorded['pon_cohort_digest']}); refusing to score it against "
            f"{pon_path.name} as well. Re-run run-all to rescore the whole sample, or use "
            "the same PON."
        )


class ConfigurationMismatch(ValueError):
    """A run whose settings differ from the ones the folder's record holds. Raised before
    anything is written."""


#: The read filters extraction applied, as the record names them. They decide which fragments
#: the BED holds, and every table is made from those fragments or, for the motif and Region-MDS
#: tables, from the BAM read again under the filters that tool takes. A rerun with other filters
#: therefore mixes two fragment sets under one record: run-all kept a BED `extract --mapq 30` wrote
#: and made its motif tables at mapq 20. (Region-MDS takes no duplicate or proper-pair filter: it
#: always drops duplicates, whatever the record says. Neither it nor mFSD applies the blacklist.)
#:
#: The blacklist is one of them: extraction drops the fragments it covers, as it drops those under
#: the mapq. Recorded by content (`blacklist_digest`), not by name: two files can share a name, and
#: the bundled one sits at another path in every install.
#:
#: Recorded truthfully only since the duplicate and proper-pair columns were added: before, every
#: record said 20/65/1000 whatever the run used, and held neither column, nor any blacklist. So a
#: record holding all six can be compared, and one that does not holds values nobody can trust
#: (`filters_vouched`).
FILTER_FIELDS = (
    "filters_mapq",
    "filters_min_length",
    "filters_max_length",
    "filters_skip_duplicates",
    "filters_require_proper_pair",
    "filters_blacklist_sha256",
)

#: The settings that decide which tables a run writes, or what they are made from, as the record
#: names them. A rerun with other values writes another set of tables, and whatever the first run
#: wrote that the second does not is left looking current: a panel folder rerun without its assay
#: rewrites `TFBS` but not `TFBS.ontarget`, under one `tfbs=ran` (Greptile on #98). `genome`,
#: `assay`, `target_regions` and the filters are written at extraction; `fsc_e1_aggregation` by
#: run-all's own stamp.
CONFIGURATION_FIELDS = (
    "genome",
    "assay",
    "target_regions",
    "fsc_e1_aggregation",
    *FILTER_FIELDS,
)


def blacklist_digest(path: Optional[Path]) -> str:
    """The `filters_blacklist_sha256` of a blacklist: the SHA-256 of its bytes, or "" for none
    applied. The same file at any path gives the same value, and an edited one another.
    """
    if path is None:
        return ""
    import hashlib

    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def filter_settings(
    *,
    mapq: int,
    minlen: int,
    maxlen: int,
    skip_duplicates: Optional[bool] = None,
    require_proper_pair: Optional[bool] = None,
    blacklist: Optional[Path] = None,
) -> Dict[str, Any]:
    """The `FILTER_FIELDS` a run takes, for `refuse_a_different_configuration`. A filter the entry
    point does not take (None) is left out, so it is not compared. ``blacklist`` is the file the
    run applies (`resolve_exclude_regions`), compared by `blacklist_digest`."""
    values = {
        "filters_mapq": mapq,
        "filters_min_length": minlen,
        "filters_max_length": maxlen,
        "filters_skip_duplicates": skip_duplicates,
        "filters_require_proper_pair": require_proper_pair,
        "filters_blacklist_sha256": (
            blacklist_digest(blacklist) if blacklist is not None else None
        ),
    }
    return {k: v for k, v in values.items() if v is not None}


def filters_vouched(frame: Any) -> bool:
    """Whether a record's filters are what its run applied: it holds all six `FILTER_FIELDS`,
    which only a record written since they were recorded truthfully does."""
    return frame is not None and all(f in frame.columns for f in FILTER_FIELDS)


def record_vouches_for_filters(
    output_dir: Path, sample: str, *, output_format: str, compress: bool
) -> bool:
    """Whether the sample's record says which filters its BED was extracted with.

    For a cache: run-all reuses a BED only when this holds and the configuration guard has passed,
    so the BED's filters are this run's. Without a record, or with one from before the filters
    were recorded truthfully, a BED beside it was made with filters nobody knows.
    """
    frame = _read_record(
        record_path(output_dir, sample, output_format=output_format, compress=compress)
    )
    return frame is not None and not frame.empty and filters_vouched(frame)


#: The record column naming the extensions the folder's tables are written in (`table_formats`).
#: Not one of `CONFIGURATION_FIELDS`: every run decides its own format, so the guard takes it from
#: the run's `output_format`/`compress` rather than from ``given``, and every writer of the record
#: records it. A rerun in another format leaves the first run's tables beside its own, and each
#: reader then picks by its own rule: downstream reads the old Parquet.
TABLE_FORMATS = "table_formats"


def _as_table_formats(extensions: Sequence[str]) -> str:
    """The `table_formats` a run writing exactly these extensions records, else the extensions
    joined: a set no single run writes, which equals no run's value."""
    for output_format in ("parquet", "tsv", "both"):
        for compress in (False, True):
            if set(
                table_extensions(output_format=output_format, compress=compress)
            ) == set(extensions):
                return table_formats(output_format=output_format, compress=compress)
    return ",".join(sorted(extensions))


def _recorded_table_formats(
    output_dir: Path, sample: str, frame: Any
) -> Tuple[str, bool]:
    """The folder's `table_formats`, and whether it is certain enough to advise a rerun with it.

    As recorded, unless a copy of the record exists in an extension the value does not name: a
    run in another format wrote that copy (a TSV run published over a Parquet one), so the folder
    is mixed, and the value says so and equals no run's. A record written before the field
    existed (every 0.9.0 folder) is judged by its own extensions, since `write_table` writes the
    record in every format a run writes: one ``.parquet`` record is ``.parquet``, a ``.parquet``
    and a ``.tsv`` are ``both``. The residual: an older Parquet run followed by an older TSV run
    also leaves two records, so a ``both`` rerun is let through, and no refusal advises one.
    """
    base = Path(output_dir) / f"{sample}.metadata"
    present = [e for e in TABLE_EXTENSIONS if (base.parent / (base.name + e)).exists()]
    value = (
        _as_text(frame[TABLE_FORMATS].iloc[0]) if TABLE_FORMATS in frame.columns else ""
    )
    if value:
        named = value.split(",")
        stray = [e for e in present if e not in named]
        if stray:
            return f"{_as_table_formats(named + stray)} (mixed)", False
        return value, True
    return _as_table_formats(present), len(present) <= 1


def _format_flags(value: str) -> Optional[str]:
    """The options a run needs to write `table_formats` ``value``, or None for a mixed folder."""
    for output_format in ("parquet", "tsv", "both"):
        for compress in (False, True):
            if table_formats(output_format=output_format, compress=compress) == value:
                return f"--output-format {output_format}" + (
                    " --compress" if compress else ""
                )
    return None


def _setting(field: str, value: Any) -> str:
    """A setting in the form two runs can compare."""
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return ""
    text = str(value).strip()
    if field == "genome":
        from krewlyzer.assets import genome_build

        try:
            return genome_build(text) if text else ""
        except ValueError:
            return text.lower()
    if field == "assay":
        return text.lower()
    if field == "target_regions":
        # By file name: Nextflow stages the BED at a new path in every task, and the record holds
        # whatever path extraction was given.
        return Path(text).name if text else ""
    if field in ("filters_mapq", "filters_min_length", "filters_max_length"):
        return text
    if field == "filters_blacklist_sha256":
        return text.lower()
    # fsc_e1_aggregation and the two filter flags. Not `bool(text)`: "False" read back from a TSV
    # record is truthy.
    return "true" if text.lower() in ("true", "1", "1.0") else "false"


def recorded_settings(
    output_dir: Path, sample: str, *, output_format: str, compress: bool
) -> Dict[str, Any]:
    """The `CONFIGURATION_FIELDS` the folder's record holds, for a writer that replaces the record
    but does not decide all of them. `extract` no longer reaches here over run-all's record (it
    refuses it); run-all does, when it re-extracts a deleted BED, and stamps its own settings only
    after its feature pass. Carrying them over means a run-all killed in between still leaves the
    folder's E1 setting in the record, so a rerun with the other one is refused."""
    frame = _read_record(
        record_path(output_dir, sample, output_format=output_format, compress=compress)
    )
    if frame is None or frame.empty:
        return {}
    return {f: frame[f].iloc[0] for f in CONFIGURATION_FIELDS if f in frame.columns}


def refuse_a_different_configuration(
    output_dir: Path,
    sample: str,
    given: Mapping[str, Any],
    *,
    output_format: str,
    compress: bool,
) -> None:
    """Refuse a run whose settings differ from the ones the folder's record holds.

    Called by every entry point that writes a settings-dependent set of tables, before it writes
    anything: run-all, `extract`, the `run_features` tools, `motif`, `region-entropy`,
    `region-mds`, `uxm` (the genome picks its markers) and `mfsd` (its mapq and lengths); not
    run-all's own calls or the PON build, which own the record they write.
    ``given`` holds the settings this entry point takes, by `CONFIGURATION_FIELDS` name; a setting
    it does not take is not compared. Nothing to compare in a fresh folder, nor for a field the
    record does not hold (written before the field existed), nor for the filters of a record that
    does not hold all six (`filters_vouched`): their values may not be what that run applied.
    run-all then extracts its BED again rather than reuse it (`record_vouches_for_filters`).

    The output format is always compared, from ``output_format``/``compress``, and is not a
    ``given`` setting (`TABLE_FORMATS`): a record without the field is judged by the extensions of
    its record files (`_recorded_table_formats`).

    Raises `ConfigurationMismatch`, naming each differing setting with both values, and saying to
    use a fresh output folder: a legitimate change of settings changes every table anyway.
    """
    unknown = sorted(set(given) - set(CONFIGURATION_FIELDS))
    if unknown:
        raise ValueError(f"unknown setting(s) {unknown}; one of {CONFIGURATION_FIELDS}")
    frame = _read_record(
        record_path(output_dir, sample, output_format=output_format, compress=compress)
    )
    if frame is None or frame.empty:
        return
    differing = []
    vouched = filters_vouched(frame)
    if not vouched and any(f in FILTER_FIELDS for f in given):
        # run-all extracts the BED again in this case; a tool reading the BAM cannot.
        logger.warning(
            f"{sample}'s record in {Path(output_dir).name}/ does not say which read filters its "
            "BED was extracted with (written before the filters were recorded), so this run's "
            "filters cannot be compared with them"
        )
    for field, value in given.items():
        if field not in frame.columns or (field in FILTER_FIELDS and not vouched):
            continue
        recorded, ours = _setting(field, frame[field].iloc[0]), _setting(field, value)
        if recorded != ours:
            differing.append(
                f"{field}={recorded or '(none)'} there, {ours or '(none)'} here"
            )
    rerun_with = ""
    recorded_formats, certain = _recorded_table_formats(output_dir, sample, frame)
    our_formats = table_formats(output_format=output_format, compress=compress)
    if recorded_formats != our_formats:
        differing.append(
            f"{TABLE_FORMATS}={recorded_formats} there, {our_formats} here"
        )
        flags = _format_flags(recorded_formats) if certain else None
        rerun_with = (
            f" This folder's tables are {recorded_formats}: rerun with {flags}."
            if flags
            else f" This folder already holds tables in {recorded_formats}."
        )
    if differing:
        raise ConfigurationMismatch(
            f"{sample} in {Path(output_dir).name}/ was written with other settings: "
            + "; ".join(differing)
            + ". A rerun with other settings would leave the first run's extra tables looking "
            "current." + rerun_with + " Use a fresh output folder, or remove this one."
        )


class RecordInUse(ValueError):
    """A record later stages wrote into, which an extraction would replace. Raised before
    anything is written."""


def feature_tables(output_dir: Path, sample: str) -> list:
    """The names of the sample's declared tables made from its fragments: every table the output
    contract declares but extraction's own (`contract.feature_table_suffixes`), in any format.
    Undeclared ones (`.TFBS.raw`, `features.json`) are not looked for: a folder holding them holds
    declared tables too. One of them beside a BED means the BED cannot be replaced without leaving
    it stale."""
    from krewlyzer.validate.contract import feature_table_suffixes

    folder = Path(output_dir)
    return [
        name
        for suffix in feature_table_suffixes()
        for ext in TABLE_EXTENSIONS
        if (
            folder / (name := f"{sample}{suffix.removesuffix('.parquet')}{ext}")
        ).exists()
    ]


def refuse_to_replace_a_later_record(output_dir: Path, sample: str) -> None:
    """Refuse an extraction into a folder where a later stage read the BED or wrote the record.

    `extract` writes the record from scratch, with its own fields only. Over run-all's record it
    dropped every provenance column and every `{stage}_status` (about thirty columns), and the
    sample still looked complete to a reader of the completion marker alone. Updating only its own
    fields, as single tools do, would not do either: a new BED (other filters) would leave every
    table the later stages made from the old one under their recorded `ran`.

    Two signs of a later stage, either of which refuses:

    - **A feature table** for the sample, any table the output contract declares but extraction's
      own (`contract.feature_table_suffixes`), in any format. The record cannot be the only sign:
      a run-all killed between writing the record and stamping it, or a folder written before
      provenance existed, holds tables beside a record that looks like extraction's alone.
    - **A record holding `PROVENANCE_COLUMNS` or a stage status**, which run-all and the single
      tools add. Every copy is read, not only this run's format: a single tool updates only the
      copies in its own (`record_stage_outcomes`).

    A folder extraction wrote alone has neither, so a rerun of `extract` is allowed (with the
    folder's settings: `refuse_a_different_configuration` runs first).
    """
    from krewlyzer.validate.contract import STAGES

    folder = Path(output_dir)
    tables = feature_tables(folder, sample)
    if tables:
        raise RecordInUse(
            f"{folder.name}/ holds tables made from this sample's BED ({len(tables)}, among them "
            f"{', '.join(tables[:3])}). Extracting again would replace the BED and the record "
            "under them, and leave them looking current. Use a fresh output folder, or remove "
            "this one."
        )

    later = set(PROVENANCE_COLUMNS) | {f"{stage}_status" for stage in STAGES}
    base = folder / f"{sample}.metadata"
    for ext in TABLE_EXTENSIONS:
        path = base.parent / (base.name + ext)
        if not path.exists():
            continue
        frame = _read_record(path)
        found = sorted(later & set(frame.columns)) if frame is not None else []
        if found:
            raise RecordInUse(
                f"{path.name} in {Path(output_dir).name}/ holds what later stages recorded "
                f"({len(found)} columns, among them {', '.join(found[:3])}). Extracting again "
                "would replace that record with extraction's alone, and leave the tables those "
                "stages made from this folder's BED looking current. Use a fresh output folder, "
                "or remove this one."
            )


def record_names_pon(
    output_dir: Path,
    sample: str,
    pon: Any,
    pon_path: Optional[Path],
    *,
    output_format: str,
    compress: bool,
) -> bool:
    """Whether the sample's record says it was scored against exactly this PON (or, for
    ``pon=None``, against none).

    For a cache: a table made under one PON is not reusable under another. run-all reused its
    motif tables whatever the PON, so a run without one and then a run with one shipped
    ``pon_applied=True`` beside an MDS with no z-score. No record, or one without provenance,
    names nothing, so nothing cached can be trusted against it.
    """
    existing = record_path(
        output_dir, sample, output_format=output_format, compress=compress
    )
    frame = _read_record(existing)
    if frame is None or frame.empty or "pon_applied" not in frame.columns:
        return False
    ours = provenance_fields(pon, pon_path)
    if _flag(frame["pon_applied"].iloc[0]) != ours["pon_applied"]:
        return False
    if not ours["pon_applied"]:
        return True
    # By content only. A record without a digest (written before it existed) cannot show the PON
    # is unchanged, so nothing cached under it is reused.
    recorded = _as_text(
        frame["pon_digest"].iloc[0] if "pon_digest" in frame.columns else None
    )
    return bool(recorded) and recorded == ours["pon_digest"]


def recorded_status(
    output_dir: Path,
    sample: str,
    stage: str,
    *,
    output_format: str,
    compress: bool,
) -> Optional[str]:
    """The status the sample's record holds for ``stage``, or None when it records none."""
    existing = record_path(
        output_dir, sample, output_format=output_format, compress=compress
    )
    column = f"{stage}_status"
    frame = _read_record(existing)
    if frame is None or frame.empty or column not in frame.columns:
        return None
    return _as_text(frame[column].iloc[0]) or None


def record_pon_used(
    output_dir: Path,
    sample: str,
    pon: Any,
    pon_path: Optional[Path],
    *,
    output_format: str,
    compress: bool,
) -> bool:
    """A standalone scorer's entry in the sample's record, after it has scored.

    `motif`, `region-entropy` and `region-mds` do not go through `run_features`, so without this
    they checked the record (`refuse_a_second_pon`) but never wrote it: after `extract`, one could
    score with PON A and a later tool with PON B, unnoticed. Raise-only, like every single tool.
    """
    return stamp_metadata(
        output_dir,
        sample,
        pon=pon,
        pon_path=pon_path if pon is not None else None,
        output_format=output_format,
        compress=compress,
        owns_record=False,
    )


def stamp_metadata(
    output_dir: Path,
    sample: str,
    pon: Any = None,
    pon_path: Optional[Path] = None,
    *,
    output_format: str,
    compress: bool,
    gc_factors: Optional[dict] = None,
    owns_record: bool = True,
    settings: Optional[Mapping[str, Any]] = None,
) -> bool:
    """Add the provenance columns to `{sample}.metadata.*`. True if written.

    `gc_factors` maps "off_target"/"on_target" to the `GcFactorsChoice` the run used; it records
    which kind of factor table each population was corrected with, or why none was (redesign
    item 0: production ran panel WPS, TFBS and ATAC uncorrected and nothing said so).

    A no-op when there is no metadata table, which is the normal case for a
    single-feature CLI run: those write one feature table and no completion
    marker. Logged at debug rather than warned, because it is not a fault.

    `owns_record` is True for run-all and the PON build, which rewrite the whole record. A single
    tool (False) fills it in where none exists, and otherwise may only raise ``pon_applied``, so
    that flag does not depend on the order single tools ran in. Which PON is recorded is kept
    consistent by `refuse_a_second_pon`, which stops a single tool before any work.

    Read-modify-write rather than appended at creation: the metadata table is
    written during extraction, before any PON is loaded, so the PON fields are
    not knowable yet at that point.

    ``settings`` are `CONFIGURATION_FIELDS` this run decided that extraction could not record
    (run-all's E1 aggregation), written only when the run owns the record: a single tool must not
    rewrite the configuration a later rerun is checked against.
    """
    base = output_dir / f"{sample}.metadata"
    # Read the file this run's format writes (`record_path`), then write all of them. Resolving
    # `.tsv` first and then reading Parquet-first meant a TSV run beside an earlier run's
    # `.parquet` read the stale row and wrote it back over this run's.
    existing = record_path(
        output_dir, sample, output_format=output_format, compress=compress
    )
    if existing is None:
        logger.debug(
            f"No metadata table for {sample}; skipping provenance stamp. "
            "Expected for a single-feature run, which writes no completion marker."
        )
        return False

    frame = _read_record(existing)
    if frame is None or frame.empty:
        logger.warning(
            f"{existing.name} could not be read; provenance not recorded. "
            "Downstream cannot tell which build produced this sample."
        )
        return False

    fields = provenance_fields(pon, pon_path)
    recorded = "pon_applied" in frame.columns
    if not owns_record and recorded:
        # A single tool keeps the sample's record (run-all's, or an earlier tool's): rewriting
        # pon_applied and gc_factors_* from one tool's run made validate-output stop checking
        # every other table. It may only *raise* pon_applied, so any order of single tools ends
        # with the same record: a table one of them left unscored is then flagged, not hidden.
        if not (fields["pon_applied"] and not _flag(frame["pon_applied"].iloc[0])):
            logger.info(
                f"{existing.name} already records how this sample was run; a single tool keeps it"
            )
            return False
        gc_factors = (
            None  # raising the PON record says nothing about this tool's GC tables
        )

    for name, value in fields.items():
        frame[name] = value
    if owns_record:
        for name, value in (settings or {}).items():
            if name not in CONFIGURATION_FIELDS:
                raise ValueError(
                    f"unknown setting {name!r}; one of {CONFIGURATION_FIELDS}"
                )
            frame[name] = value
    for population, choice in (gc_factors or {}).items():
        tag = population.replace("_", "")  # off_target -> offtarget
        frame[f"gc_factors_{tag}"] = choice.record()
        frame[f"gc_correction_applied_{tag}"] = bool(choice.found)
    frame[TABLE_FORMATS] = table_formats(output_format=output_format, compress=compress)

    write_table(frame, base, output_format=output_format, compress=compress)

    if fields["pon_applied"]:
        logger.info(
            f"Provenance: krewlyzer {fields['krewlyzer_version']}, scored against "
            f"{fields['pon_model']} (cohort {fields['pon_cohort_digest'] or 'unrecorded'})"
        )
    else:
        logger.info(
            f"Provenance: krewlyzer {fields['krewlyzer_version']}, no PON applied — "
            "every z-score column is absent by design, not by failure."
        )
    return True


#: Columns read back from a TSV record as text: a digit-only cohort digest such as '0123…' read
#: back as an integer, and a status column must never become a number either.
_RECORD_TEXT = (
    "pon_model",
    "pon_cohort_digest",
    "pon_digest",
    "pon_krewlyzer_version",
    "krewlyzer_version",
)


def _read_record(path: Optional[Path]):
    """A sample record, its text columns read as text, every stage column included.

    As text: a digit-only digest such as '0123…' reads back from TSV as an integer, and an empty
    one as NaN. One dtype map for every reader, so a column added to the record is read the same
    way everywhere.
    """
    if path is None:
        return None
    from krewlyzer.validate.contract import STAGES

    text = {name: str for name in _RECORD_TEXT}
    text.update({f"{s}_{f}": str for s in STAGES for f in ("status", "error")})
    return read_exact_table(path, dtype=text)


@dataclass(frozen=True)
class Outcome:
    """What one optional stage did: a status from `contract.STAGE_STATUSES`, and for a failure
    the exception's **class name**.

    Never the message. Messages carry BAM, reference and output paths, and sample directories
    are named for the patient (invariant #4); the full message stays in the log.
    """

    status: str
    error_class: str = ""

    def __post_init__(self) -> None:
        from krewlyzer.validate.contract import STAGE_STATUSES

        if self.status not in STAGE_STATUSES:
            raise ValueError(
                f"unknown stage status {self.status!r}; one of {STAGE_STATUSES}"
            )
        if (self.status == "failed") != bool(self.error_class):
            raise ValueError(
                "a failed outcome carries the exception's class name, and only a failed one"
            )
        if self.error_class and not self.error_class.isidentifier():
            raise ValueError(
                "error_class is an exception's class name, never its message: "
                f"{self.error_class[:40]!r}"
            )

    @classmethod
    def failed(cls, exc: BaseException) -> "Outcome":
        return cls("failed", type(exc).__name__)


def record_stage_outcomes(
    output_dir: Path,
    sample: str,
    outcomes: Mapping[str, Outcome],
    *,
    output_format: str,
    compress: bool,
    owns_record: bool = True,
) -> bool:
    """Write `{stage}_status` and `{stage}_error` into `{sample}.metadata.*`. True if written.

    The one writer for stage outcomes.

    ``owns_record`` is True for run-all, which records every stage on every run. It reads the
    record its own format writes (`record_path`), exactly, and writes all of its formats, as
    `stamp_metadata` does.

    A single tool (False) passes only the stages whose outputs it rewrote, and records each as it
    passes it, whatever the status. `krewlyzer fsc` without `--assay` did not run gene FSC, so it
    passes no `gene_fsc` (run-all may have run it). `motif --skip-pon` rewrote `MDS` without
    `mds_z`, so it passes `mds_z=not_requested`: the z-score stage lives in a table the tool just
    rewrote, and run-all's `ran` there would describe a column that is gone. (A filter here once
    dropped every `not_requested` from a single tool, which kept that `ran` -- Greptile on #98.)
    The last writer wins, because it rewrote that stage's outputs -- in the format it wrote them. So it updates
    only the copies of the record in its own format(s) (`tsv` the `.tsv`, or the `.tsv.gz`
    with `--compress`; `parquet` the `.parquet`; `both` both), and creates none. A tool at
    `-F tsv` after a Parquet run-all rewrote `.TFBS.tsv`, not the `.TFBS.parquet` downstream
    reads; recording `ran` in `.metadata.parquet` would have hidden that stale table from the
    gate. And a `.metadata.parquet` created in a TSV directory is a completion marker.

    A no-op when there is no metadata table: a single tool run into an empty directory writes no
    completion marker.
    """
    from krewlyzer.validate.contract import STAGES

    unknown = sorted(set(outcomes) - set(STAGES))
    if unknown:
        raise ValueError(
            f"unknown stage(s) {unknown}; declare them in validate/contract.py first"
        )
    if not outcomes:
        return False

    base = output_dir / f"{sample}.metadata"
    if owns_record:
        existing = record_path(
            output_dir, sample, output_format=output_format, compress=compress
        )
        targets = [(existing, output_format, compress)] if existing is not None else []
    else:
        written_by_this_run = [
            (
                (path, "parquet", False)
                if path.suffix == ".parquet"
                else (path, "tsv", compress)
            )
            for path in table_paths(
                base, output_format=output_format, compress=compress
            )
        ]
        targets = [t for t in written_by_this_run if t[0].exists()]
    if not targets:
        logger.debug(f"No metadata table for {sample}; stage outcomes not recorded.")
        return False

    written = False
    for path, fmt, gz in targets:
        frame = _read_record(path)
        if frame is None or frame.empty:
            logger.warning(
                f"{path.name} could not be read; stage outcomes not recorded in it. "
                "Downstream cannot tell a failed stage from one that was not requested."
            )
            continue
        for stage, outcome in outcomes.items():
            frame[f"{stage}_status"] = outcome.status
            frame[f"{stage}_error"] = outcome.error_class
        # The run's formats, not this copy's: one copy of a `both` record is `.parquet`.
        frame[TABLE_FORMATS] = table_formats(
            output_format=output_format, compress=compress
        )
        write_table(frame, base, output_format=fmt, compress=gz)
        written = True
    if not written:
        return False
    failed = [
        f"{k} ({o.error_class})" for k, o in outcomes.items() if o.status == "failed"
    ]
    logger.info(
        "Stage outcomes: "
        + ", ".join(f"{k}={o.status}" for k, o in outcomes.items())
        + (f"; failed: {', '.join(failed)}" if failed else "")
    )
    return True
