"""
Shared output utilities for Krewlyzer Python modules.

Provides two boilerplate-free I/O helpers used by all Python feature writers:

- ``write_table``:  Write a DataFrame to TSV and/or Parquet based on
  ``--output-format`` and ``--compress`` flags. Replaces all ``df.to_csv()``
  call sites in FSR, FSC, region entropy, and metadata writers.

- ``read_table``:  Parquet-first reader that falls back to plain TSV or
  gzip-compressed TSV. Replaces ``pd.read_csv(path, sep="\\t")`` everywhere a
  file may now be in any of the three formats.

Design principles
-----------------
* Single responsibility — each helper does one thing well.
* No hidden state — all behaviour is controlled by explicit arguments.
* Logging at DEBUG level for every write/read so operators can trace I/O.
* Graceful degradation — ``read_table`` returns ``None`` (not an exception)
  when no candidate exists; callers can then decide what to do.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import pandas as pd

logger = logging.getLogger("krewlyzer.core.output_utils")


# ---------------------------------------------------------------------------
# Writer
# ---------------------------------------------------------------------------


def write_table(
    df: "pd.DataFrame",
    output_path: Path,
    output_format: str = "tsv",
    compress: bool = False,
    float_format: str | None = None,
) -> None:
    """Write *df* to TSV and/or Parquet based on *output_format*.

    Parameters
    ----------
    df:
        DataFrame to persist.
    output_path:
        Base path **without** a trailing extension, e.g.
        ``Path("out/sample.FSR")``.  The correct extension(s) are appended
        automatically:

        * ``tsv``     → ``{output_path}.tsv``  (or ``.tsv.gz`` with *compress*)
        * ``parquet`` → ``{output_path}.parquet``
        * ``both``    → both of the above

    output_format:
        One of ``"tsv"``, ``"parquet"``, or ``"both"``.  Any other value
        is treated as ``"tsv"`` with a warning.
    compress:
        When *True* and TSV is being written, compress with gzip
        (``.tsv.gz``).  Has no effect on Parquet output.
    float_format:
        Optional ``printf``-style format string for floating-point columns
        (e.g. ``"%.6f"``).  Forwarded to ``DataFrame.to_csv``.

    Raises
    ------
    ValueError
        If *output_format* is not one of the recognised values (after warning
        it will fall back to ``"tsv"``, so this is only raised when the caller
        explicitly wants strict validation — currently unused).
    IOError
        If the underlying file-system write fails.

    Examples
    --------
    >>> write_table(df, Path("out/sample.FSR"), output_format="both", compress=False)
    # → writes out/sample.FSR.tsv  and  out/sample.FSR.parquet
    >>> write_table(df, Path("out/sample.FSR"), output_format="tsv", compress=True)
    # → writes out/sample.FSR.tsv.gz
    """

    if output_format not in ("tsv", "parquet", "both"):
        logger.warning(
            "write_table: unknown output_format %r — defaulting to 'tsv'", output_format
        )
        output_format = "tsv"
    # Normalize base path: strip any existing known extension so callers that pass
    # 'sample.FSC.tsv' (legacy) and callers that pass 'sample.EndMotif' (no extension)
    # both work correctly.
    # IMPORTANT: do NOT use .with_suffix('') — Python's with_suffix() replaces only the
    # last dot-segment, so 'sample.EndMotif'.with_suffix('') gives 'sample' (treats
    # 'EndMotif' as an extension). Instead, strip only well-known format extensions.
    _name = output_path.name
    for _known in (".tsv.gz", ".tsv", ".parquet"):
        if _name.endswith(_known):
            _name = _name[: -len(_known)]
            break
    _base = output_path.parent / _name

    if output_format in ("tsv", "both"):
        ext = ".tsv.gz" if compress else ".tsv"
        tsv_path = _base.parent / (_base.name + ext)
        logger.debug(
            "write_table: writing TSV (%d rows \u00d7 %d cols) \u2192 %s",
            len(df),
            len(df.columns),
            tsv_path,
        )
        df.to_csv(
            tsv_path,
            sep="\t",
            index=False,
            float_format=float_format,
            compression="gzip" if compress else None,
        )
        logger.info("write_table: wrote %s", tsv_path.name)

    if output_format in ("parquet", "both"):
        parquet_path = _base.parent / (_base.name + ".parquet")
        logger.debug(
            "write_table: writing Parquet (%d rows \u00d7 %d cols) \u2192 %s",
            len(df),
            len(df.columns),
            parquet_path,
        )
        df.to_parquet(parquet_path, index=False, engine="pyarrow")
        logger.info("write_table: wrote %s", parquet_path.name)


def cleanup_intermediate_tsv(
    tsv_path: Path,
    output_format: str,
    compress: bool,
) -> None:
    """Remove the raw ``.tsv`` left by Rust after ``write_table()`` produced the target format.

    Call this **after** ``write_table()`` whenever the source was a Rust-produced
    ``.tsv`` that has been re-written as ``.tsv.gz`` and/or ``.parquet``.

    Deletes *tsv_path* when:

    - ``output_format == "parquet"`` — only ``.parquet`` was written.
    - ``compress is True`` — ``.tsv.gz`` replaces the original ``.tsv``.

    Does **not** delete when ``output_format`` is ``"tsv"`` or ``"both"`` and
    ``compress`` is ``False``, because the ``.tsv`` *is* the desired output.

    Args:
        tsv_path: Path to the intermediate ``.tsv`` file.
        output_format: ``"tsv"``, ``"parquet"``, or ``"both"``.
        compress: Whether gzip compression was applied.
    """
    if not tsv_path.exists():
        return
    if output_format == "parquet" or compress:
        tsv_path.unlink()
        logger.debug("Cleaned up intermediate TSV: %s", tsv_path.name)


# ---------------------------------------------------------------------------
# Reader
# ---------------------------------------------------------------------------


def read_table(path: Path, **csv_kwargs) -> "pd.DataFrame | None":
    """Parquet-first reader with TSV and ``.tsv.gz`` fallback.

    Tries candidate paths in the following priority order:

    1. ``{path.stem}.parquet``      (format-free stem + .parquet)
    2. ``{path}``                   (as given, typically ``.tsv``)
    3. ``{path}.gz``                (gzip-compressed)

    Parameters
    ----------
    path:
        The *expected* TSV path (e.g. ``sample.FSR.tsv``).  The Parquet
        candidate is derived by replacing the ``.tsv`` suffix with
        ``.parquet``.
    **csv_kwargs:
        Extra keyword arguments forwarded to ``pd.read_csv`` when a TSV or
        ``.tsv.gz`` candidate is selected (e.g. ``dtype``, ``usecols``).

    Returns
    -------
    DataFrame | None
        The loaded DataFrame, or ``None`` if no candidate was found on disk.

    Examples
    --------
    >>> df = read_table(Path("out/sample.FSR.tsv"))
    >>> if df is None:
    ...     raise FileNotFoundError("FSR output missing")
    """
    import pandas as pd  # local import

    # NOTE: Do NOT use path.with_suffix('.parquet') — same compound-extension bug.
    # Use name+ext pattern to find the Parquet sibling of a TSV base path.
    parquet_candidate = (
        path.parent / (path.stem + ".parquet")
        if path.suffix == ".tsv"
        else path.parent / (path.name + ".parquet")
    )
    candidates = [
        parquet_candidate,
        path,
        Path(str(path) + ".gz"),
    ]

    for candidate in candidates:
        if candidate.exists():
            logger.debug("read_table: found %s", candidate)
            if candidate.suffix == ".parquet":
                df = pd.read_parquet(candidate)
                logger.debug(
                    "read_table: loaded Parquet %s (%d rows × %d cols)",
                    candidate.name,
                    len(df),
                    len(df.columns),
                )
                return df
            # Plain TSV or gzip TSV
            from typing import Literal, cast

            compression_arg = cast(
                "Literal['gzip'] | None",
                "gzip" if str(candidate).endswith(".gz") else None,
            )
            df = pd.read_csv(  # type: ignore[assignment]
                candidate, sep="\t", compression=compression_arg, **csv_kwargs
            )
            logger.debug(
                "read_table: loaded TSV %s (%d rows × %d cols)",
                candidate.name,
                len(df),
                len(df.columns),
            )
            return df

    logger.debug(
        "read_table: no file found for %s (tried %d candidates)", path, len(candidates)
    )
    return None
