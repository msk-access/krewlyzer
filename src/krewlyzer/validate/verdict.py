"""Cross-axis agreement: does more than one lens point the same way?

The strongest claim a single-sample fragmentomics report can make is that
**independent measurements converge**. Fragment size, nuclease cutting, tissue
shedding and chromatin accessibility are four different physical observations of
the same event — tumour cells dying and shedding DNA. One elevated metric is a
hypothesis; four agreeing is evidence.

So the summary is *"N of M axes agree"*, with every axis's value shown, never a
single score. A composite number would hide exactly the disagreement a reader
needs to see.

Three things this deliberately does not do
------------------------------------------

**It does not invent a threshold.** ``z_threshold`` defaults to 2.0 because that
is the conventional flag, and the report says so. It is an argument, not a
constant, and the rendered output labels it as a convention rather than a
cut-off. Every numeric band in this project that was treated as a cut-off turned
out to be a display default or refuted outright.

**It does not treat absence as disagreement.** An axis with no PON z-score is
*not assessable* — mFSD needs a variant list, OCF is hg19-only, z-scores need a
panel of normals. Scoring those as "does not agree" would make a thinner run
look like a healthier one.

**It does not decide anything clinical.** It reports convergence between
measurements. What that means for a patient is not krewlyzer's call.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import List, Optional

import pandas as pd

from krewlyzer.core.output_utils import read_table
from krewlyzer.validate.contract import OCF_PREFERENCE

#: Conventional |z| flag. An argument everywhere it is used, never a constant
#: baked into a comparison -- see the module docstring.
DEFAULT_Z_THRESHOLD = 2.0


class Support(str, Enum):
    """What one axis says."""

    #: Moves in the direction associated with tumour cfDNA, past the threshold.
    AGREES = "agrees"
    #: Measured, but not past the threshold.
    QUIET = "quiet"
    #: Measured, and moving the *other* way. Worth seeing, not hiding.
    OPPOSES = "opposes"
    #: The input this axis needs was not present.
    NOT_ASSESSABLE = "not assessable"


@dataclass
class Axis:
    """One independent line of evidence."""

    name: str
    #: Column read, qualified by table -- so a reader can check it themselves.
    source: str
    #: Which way tumour cfDNA moves this metric, in words.
    direction: str
    support: Support
    value: Optional[float] = None
    #: What was examined, e.g. "strongest of 23 cancer types".
    detail: Optional[str] = None
    #: Why it could not be assessed, when that is the case.
    reason: Optional[str] = None


@dataclass
class Verdict:
    axes: List[Axis]
    z_threshold: float

    @property
    def assessable(self) -> List[Axis]:
        return [a for a in self.axes if a.support is not Support.NOT_ASSESSABLE]

    @property
    def agreeing(self) -> List[Axis]:
        return [a for a in self.axes if a.support is Support.AGREES]

    @property
    def opposing(self) -> List[Axis]:
        return [a for a in self.axes if a.support is Support.OPPOSES]

    @property
    def headline(self) -> str:
        n, m = len(self.agreeing), len(self.assessable)
        if m == 0:
            return "No axis could be assessed"
        return f"{n} of {m} assessable axes agree"

    @property
    def summary(self) -> str:
        """One sentence a reader can act on, without overclaiming."""
        n, m = len(self.agreeing), len(self.assessable)
        if m == 0:
            return (
                "Nothing here can be read as evidence either way: no axis had "
                "the inputs it needs."
            )
        if n == 0:
            return (
                "No axis exceeds the flag threshold. Consistent with a sample "
                "carrying little or no tumour cfDNA, and equally consistent "
                "with a burden below what these measurements resolve."
            )
        if n == 1:
            return (
                "One axis is elevated. A single axis is a hypothesis, not "
                "evidence — convergence is what distinguishes signal from a "
                "measurement artefact."
            )
        if self.opposing:
            return (
                f"{n} axes agree while {len(self.opposing)} move the other way. "
                "Worth resolving before treating this as convergent evidence."
            )
        return (
            f"{n} independent measurements agree. They observe different "
            "physical properties, so a shared artefact is an unlikely "
            "explanation."
        )


def _first_value(df: Optional[pd.DataFrame], column: str) -> Optional[float]:
    if df is None or column not in df.columns or df.empty:
        return None
    series = pd.to_numeric(df[column], errors="coerce").dropna()
    return float(series.iloc[0]) if not series.empty else None


def _extreme(
    df: Optional[pd.DataFrame], column: str, label_col: str, largest: bool
) -> tuple[Optional[float], Optional[str]]:
    """The most extreme row, and which label it belongs to.

    Reported rather than a mean: a tissue-specific or factor-specific signal is
    the point, and averaging 808 transcription factors would dilute exactly the
    thing being measured.
    """
    if df is None or column not in df.columns or df.empty:
        return None, None
    values = pd.to_numeric(df[column], errors="coerce")
    valid = values.dropna()
    if valid.empty:
        return None, None
    idx = valid.idxmax() if largest else valid.idxmin()
    label = str(df.loc[idx, label_col]) if label_col in df.columns else None
    return float(values.loc[idx]), label


def _classify(
    value: Optional[float], threshold: float, higher_is_tumour: bool
) -> Support:
    if value is None:
        return Support.NOT_ASSESSABLE
    signed = value if higher_is_tumour else -value
    if signed >= threshold:
        return Support.AGREES
    if signed <= -threshold:
        return Support.OPPOSES
    return Support.QUIET


def compute_verdict(
    directory: Path,
    sample_id: Optional[str] = None,
    z_threshold: float = DEFAULT_Z_THRESHOLD,
) -> Verdict:
    """Read the four axes from a finished sample directory."""
    directory = Path(directory)
    sample_id = sample_id or directory.name

    def table(suffix: str) -> Optional[pd.DataFrame]:
        return read_table(directory / f"{sample_id}{suffix}")

    axes: List[Axis] = []

    # -- 1. fragment size ---------------------------------------------------
    fsr = table(".FSR.parquet")
    log2 = None
    if fsr is not None and "short_long_log2" in fsr.columns:
        values = pd.to_numeric(fsr["short_long_log2"], errors="coerce").dropna()
        log2 = float(values.median()) if not values.empty else None
    axes.append(
        Axis(
            name="Fragment size",
            source="FSR.short_long_log2 (median over windows)",
            direction="higher — tumour cfDNA is shorter",
            # Already PON-normalised, so ~0 is the healthy expectation and the
            # value is directly comparable to a z. Not a z-score, so the
            # threshold is a looser analogy here; stated in the report.
            support=_classify(log2, z_threshold, higher_is_tumour=True),
            value=log2,
            detail="genome-wide median; PON-normalised, so ~0 is healthy",
            reason=None if log2 is not None else "FSR.parquet absent or unnormalised",
        )
    )

    # -- 2. nuclease cutting ------------------------------------------------
    mds_z = _first_value(table(".MDS.parquet"), "mds_z")
    axes.append(
        Axis(
            name="Nuclease cutting",
            source="MDS.mds_z",
            direction="lower — cutting becomes stereotyped",
            support=_classify(mds_z, z_threshold, higher_is_tumour=False),
            value=mds_z,
            detail="motif diversity across the whole sample",
            reason=None if mds_z is not None else "no PON z-score (needs --pon-model)",
        )
    )

    # -- 3. tissue of origin ------------------------------------------------
    ocf = next(
        (df for df in (table(sfx) for sfx in OCF_PREFERENCE) if df is not None), None
    )
    ocf_z, tissue = _extreme(ocf, "ocf_z", "tissue", largest=True)
    axes.append(
        Axis(
            name="Tissue shedding",
            source="OCF.ocf_z (most elevated tissue)",
            direction="higher — that tissue is shedding DNA",
            support=_classify(ocf_z, z_threshold, higher_is_tumour=True),
            value=ocf_z,
            detail=f"strongest tissue: {tissue}" if tissue else None,
            reason=(
                None
                if ocf_z is not None
                else "OCF absent (hg19 assets only) or carries no PON z-score"
            ),
        )
    )

    # -- 4. chromatin accessibility -----------------------------------------
    atac = table(".ATAC.parquet")
    atac_z, cancer_type = _extreme(atac, "z_score", "label", largest=True)
    axes.append(
        Axis(
            name="Chromatin accessibility",
            source="ATAC.z_score (most elevated cancer type)",
            direction="higher — size entropy shifts at that type's peaks",
            support=_classify(atac_z, z_threshold, higher_is_tumour=True),
            value=atac_z,
            detail=f"strongest type: {cancer_type}" if cancer_type else None,
            reason=None if atac_z is not None else "ATAC absent or no PON z-score",
        )
    )

    return Verdict(axes=axes, z_threshold=z_threshold)
