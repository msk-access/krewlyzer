"""Charts for the sample report, as Plotly figures.

Why Plotly rather than static images
------------------------------------
The data is high-cardinality in a way that defeats a static chart: 808
transcription factors cannot be labelled, and 28,823 FSR windows cannot be read
off an axis. Hover and zoom are not decoration here, they are the only way to
get from "there is a tail" to "the tail is these regions".

Self-contained either way. ``include_plotlyjs="inline"`` embeds the runtime, so
the report opens from a filesystem, inside a container, on a cluster with no
network. That costs ~4.6 MB once; each additional figure is ~7 KB.

Design constraints, in order of how much trouble ignoring them causes
---------------------------------------------------------------------

**A wrong chart is worse than a missing one, because it persuades.** Everything
here refuses rather than approximates: a constant column is reported as
constant, not drawn as a flat line that looks like a measurement; an absent
table produces a stated reason, not an empty axis.

**No threshold is drawn as a cut-off.** Reference lines appear only where the
value is a literature anchor, and they are labelled as anchors. Every numeric
band in this project that got treated as a cut-off turned out to be a display
default or refuted outright — the documented ATAC/TFBS "normal" entropy range
flags a healthy sample as abnormal.

**Charts follow the page theme.** Traces carry explicit colours; everything
structural is left to the layout template, which the page swaps on toggle.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Sequence

import pandas as pd

from krewlyzer.validate.contract import OCF_PREFERENCE

#: The accent, matching the docs site. "Krew" is blood in Polish.
ACCENT = "#ef5552"
ACCENT_DARK = "#ff7875"
#: Semantic, kept distinct from the accent so "the other way" never reads as
#: "the brand colour".
OPPOSE = "#3b7ea1"
MUTED = "#9aa0a6"


@dataclass
class Chart:
    """A rendered figure, or a stated reason there isn't one."""

    #: Which output table this belongs beside. Charts live with their data.
    suffix: str
    title: str
    caption: str
    #: Plotly figure, ready for `to_html`. None when not drawn.
    figure: Any = None
    reason: Optional[str] = None

    @property
    def drawn(self) -> bool:
        return self.figure is not None


def _plotly():
    try:
        import plotly.graph_objects as go

        return go
    except ImportError:
        return None


_NO_PLOTLY = "plotly not installed — `pip install krewlyzer[report]`"


def _no(suffix: str, title: str, caption: str, reason: str) -> Chart:
    return Chart(suffix=suffix, title=title, caption=caption, reason=reason)


def _is_constant(values: "pd.Series | Sequence[float]") -> bool:
    series = pd.Series(list(values)).dropna()
    return len(series) > 1 and series.nunique() == 1


def _layout(go, height: int = 280, **kwargs) -> Dict:
    """Shared layout. Colours come from the page template, not from here."""
    base = dict(
        height=height,
        margin=dict(l=52, r=16, t=8, b=40),
        showlegend=False,
        hovermode="closest",
        dragmode="pan",
        font=dict(size=12),
        xaxis=dict(showgrid=False, zeroline=False),
        yaxis=dict(gridwidth=1, zeroline=False),
    )
    base.update(kwargs)
    return base


# ---------------------------------------------------------------------------
# Charts, each keyed to the table it explains
# ---------------------------------------------------------------------------


def fragment_size_density(fsd: Optional[pd.DataFrame]) -> Chart:
    suffix, title = ".FSD.parquet", "Fragment size distribution"
    caption = (
        "Counts per 5 bp bin, summed across regions. The mononucleosomal peak "
        "near 166 bp dominates healthy plasma; tumour cfDNA shifts mass toward "
        "~145 bp. The dashed lines are literature anchors, not thresholds — "
        "drag to pan, scroll to zoom."
    )
    go = _plotly()
    if go is None:
        return _no(suffix, title, caption, _NO_PLOTLY)
    if fsd is None or fsd.empty:
        return _no(suffix, title, caption, "FSD table absent")

    bins = sorted(
        (int(m.group(1)), str(c))
        for c in fsd.columns
        if (m := re.fullmatch(r"(\d+)-(\d+)", str(c)))
    )
    if not bins:
        return _no(suffix, title, caption, "no size-bin columns in FSD")

    x = [b for b, _ in bins]
    y = [float(pd.to_numeric(fsd[c], errors="coerce").sum()) for _, c in bins]
    if sum(y) == 0:
        return _no(suffix, title, caption, "every size bin is empty")

    fig = go.Figure(
        go.Scatter(
            x=x,
            y=y,
            mode="lines",
            fill="tozeroy",
            line=dict(color=ACCENT, width=2),
            fillcolor="rgba(239,85,82,0.16)",
            hovertemplate="%{x} bp<br>%{y:,.0f} fragments<extra></extra>",
        )
    )
    for anchor, label in ((166, "166 bp healthy mode"), (145, "145 bp tumour mode")):
        fig.add_vline(
            x=anchor,
            line=dict(color=MUTED, width=1, dash="dash"),
            annotation_text=label,
            annotation_position="top",
            annotation_font_size=10,
        )
    fig.update_layout(
        **_layout(go, 320, xaxis_title="fragment length (bp)", yaxis_title="fragments")
    )
    fig.update_xaxes(range=[min(x), min(max(x), 600)])
    return Chart(suffix, title, caption, fig)


def size_channel_composition(fsc: Optional[pd.DataFrame]) -> Chart:
    suffix, title = ".FSC.parquet", "Size-channel composition"
    caption = (
        "The six non-overlapping channels as a share of all counted fragments. "
        "They partition the total, so the bar sums to 100%."
    )
    go = _plotly()
    if go is None:
        return _no(suffix, title, caption, _NO_PLOTLY)
    if fsc is None or fsc.empty:
        return _no(suffix, title, caption, "FSC table absent")

    channels = [
        "ultra_short",
        "core_short",
        "mono_nucl",
        "di_nucl",
        "long",
        "ultra_long",
    ]
    missing = [c for c in channels if c not in fsc.columns]
    if missing:
        return _no(suffix, title, caption, f"channels absent: {', '.join(missing)}")

    totals = [float(pd.to_numeric(fsc[c], errors="coerce").sum()) for c in channels]
    grand = sum(totals)
    if grand == 0:
        return _no(suffix, title, caption, "no fragments counted in any channel")

    fig = go.Figure()
    # One hue stepped in opacity: six distinct colours would imply six unrelated
    # categories rather than one ordered size axis.
    for i, (name, total) in enumerate(zip(channels, totals)):
        share = total / grand * 100
        fig.add_bar(
            x=[share],
            y=[""],
            orientation="h",
            name=name,
            marker=dict(color=ACCENT, opacity=0.35 + 0.13 * i),
            hovertemplate=f"{name}<br>%{{x:.1f}}%<br>{total:,.0f} fragments<extra></extra>",
            text=[name if share > 7 else ""],
            textposition="inside",
            insidetextanchor="middle",
            textfont=dict(size=10),
        )
    fig.update_layout(
        **_layout(go, 130, barmode="stack", xaxis_title="% of counted fragments"),
    )
    fig.update_xaxes(range=[0, 100])
    return Chart(suffix, title, caption, fig)


def short_long_spread(fsr: Optional[pd.DataFrame]) -> Chart:
    suffix, title = ".FSR.parquet", "Short:long ratio across windows"
    caption = (
        "`short_long_log2` per window, already PON-normalised — so ~0 is the "
        "healthy expectation. The spread and the right tail carry the signal: "
        "focal high-burden regions appear there even when the median is "
        "unremarkable."
    )
    go = _plotly()
    if go is None:
        return _no(suffix, title, caption, _NO_PLOTLY)
    if fsr is None or "short_long_log2" not in getattr(fsr, "columns", []):
        return _no(suffix, title, caption, "short_long_log2 absent — needs --pon-model")

    values = (
        pd.to_numeric(fsr["short_long_log2"], errors="coerce")
        .replace([float("inf"), float("-inf")], pd.NA)
        .dropna()
    )
    if values.empty:
        return _no(suffix, title, caption, "no finite values")
    if _is_constant(values):
        return _no(
            suffix, title, caption, f"constant at {values.iloc[0]:.4g} — not plotted"
        )

    fig = go.Figure(
        go.Histogram(
            x=values,
            nbinsx=80,
            marker=dict(color=ACCENT, opacity=0.8),
            hovertemplate="log2 %{x:.3f}<br>%{y:,} windows<extra></extra>",
        )
    )
    fig.add_vline(x=0, line=dict(color=MUTED, width=1))
    median = float(values.median())
    fig.add_vline(
        x=median,
        line=dict(color=OPPOSE, width=2),
        annotation_text=f"median {median:+.3f}",
        annotation_position="top",
        annotation_font_size=10,
    )
    fig.update_layout(
        **_layout(
            go,
            280,
            xaxis_title="short:long log2 (PON-normalised)",
            yaxis_title="windows",
        )
    )
    return Chart(suffix, title, caption, fig)


def end_motif_spectrum(motif: Optional[pd.DataFrame]) -> Chart:
    suffix, title = ".EndMotif.parquet", "End-motif spectrum"
    caption = (
        "All 256 4-mers at fragment 5′ ends, ranked by frequency. A healthy "
        "sample is broad; stereotyped cutting concentrates mass on the left, "
        "which is what lowers MDS. Hover to identify a motif."
    )
    go = _plotly()
    if go is None:
        return _no(suffix, title, caption, _NO_PLOTLY)
    if motif is None or "Frequency" not in getattr(motif, "columns", []):
        return _no(suffix, title, caption, "EndMotif table absent")

    work = motif.copy()
    work["Frequency"] = pd.to_numeric(work["Frequency"], errors="coerce")
    work = work.dropna(subset=["Frequency"]).sort_values("Frequency", ascending=False)
    if work.empty:
        return _no(suffix, title, caption, "no finite frequencies")
    if _is_constant(work["Frequency"]):
        return _no(
            suffix, title, caption, "every motif has the same frequency — not plotted"
        )

    labels = (
        work["Motif"].astype(str) if "Motif" in work.columns else work.index.astype(str)
    )
    fig = go.Figure(
        go.Scatter(
            x=list(range(1, len(work) + 1)),
            y=work["Frequency"],
            mode="lines",
            fill="tozeroy",
            customdata=labels,
            line=dict(color=ACCENT, width=1.6),
            fillcolor="rgba(239,85,82,0.16)",
            hovertemplate="%{customdata}<br>rank %{x}<br>%{y:.5f}<extra></extra>",
        )
    )
    uniform = 1.0 / len(work)
    fig.add_hline(
        y=uniform,
        line=dict(color=MUTED, width=1, dash="dash"),
        annotation_text="uniform",
        annotation_position="right",
        annotation_font_size=10,
    )
    fig.update_layout(
        **_layout(go, 260, xaxis_title="4-mers, ranked", yaxis_title="frequency")
    )
    return Chart(suffix, title, caption, fig)


def _ranked_bars(
    suffix: str,
    df: Optional[pd.DataFrame],
    label_col: str,
    value_col: str,
    title: str,
    caption: str,
    top: Optional[int] = None,
    absent: str = "table absent",
) -> Chart:
    go = _plotly()
    if go is None:
        return _no(suffix, title, caption, _NO_PLOTLY)
    if df is None or df.empty or value_col not in df.columns:
        return _no(suffix, title, caption, absent)

    work = df[[label_col, value_col]].copy()
    work[value_col] = pd.to_numeric(work[value_col], errors="coerce")
    work = work.dropna()
    if work.empty:
        return _no(suffix, title, caption, f"no finite {value_col} values")
    if _is_constant(work[value_col]):
        return _no(suffix, title, caption, f"{value_col} is constant — not plotted")

    total = len(work)
    if top is not None and total > top:
        work = work.reindex(work[value_col].abs().sort_values(ascending=False).index)
        work = work.head(top)
        caption += f" Showing the {top} largest |{value_col}| of {total:,}."
    work = work.sort_values(value_col)

    fig = go.Figure(
        go.Bar(
            x=work[value_col],
            y=work[label_col].astype(str),
            orientation="h",
            marker=dict(color=[ACCENT if v >= 0 else OPPOSE for v in work[value_col]]),
            hovertemplate="%{y}<br>" + value_col + " %{x:+.3f}<extra></extra>",
        )
    )
    fig.add_vline(x=0, line=dict(color=MUTED, width=1))
    fig.update_layout(
        **_layout(go, max(180, 22 * len(work) + 60), xaxis_title=value_col)
    )
    fig.update_yaxes(automargin=True)
    return Chart(suffix, title, caption, fig)


def tissue_shedding(suffix: str, ocf: Optional[pd.DataFrame]) -> Chart:
    return _ranked_bars(
        suffix,
        ocf,
        "tissue",
        "ocf_z",
        "Tissue shedding",
        "Orientation-aware fragmentation per tissue, as a PON z-score. A rise "
        "means that tissue is shedding DNA; a fall in T-cell is equally "
        "informative, since it reflects dilution of the normal haematopoietic "
        "background by tumour DNA.",
        absent="OCF absent (hg19 only) or carries no z-score",
    )


def accessibility_by_cancer_type(atac: Optional[pd.DataFrame]) -> Chart:
    return _ranked_bars(
        ".ATAC.parquet",
        atac,
        "label",
        "z_score",
        "Accessibility by cancer type",
        "Fragment-size entropy at each TCGA type's accessible peaks, as a "
        "z-score. Absolute entropy bands are refuted; only the z-scores are "
        "interpretable.",
        absent="ATAC absent or carries no z-score",
    )


def top_transcription_factors(tfbs: Optional[pd.DataFrame]) -> Chart:
    return _ranked_bars(
        ".TFBS.parquet",
        tfbs,
        "label",
        "z_score",
        "Transcription factors, strongest shifts",
        "Fragment-size entropy at each factor's binding sites, as a z-score.",
        top=25,
        absent="TFBS absent or carries no z-score",
    )


def fsc_short_long_by_bin(fsc: Optional[pd.DataFrame]) -> Chart:
    """Act 1.3 — the regional view FSR's genome median cannot show."""
    suffix, title = ".FSC.parquet", "Short:long across genomic bins"
    caption = (
        "Per-bin log2 of short (65–149 bp) over long (261–1000 bp) counts. The "
        "genome median hides heterogeneity; a tumour with focal high-burden "
        "regions shows a right-skewed tail here even when the median is "
        "unremarkable. Raw counts, not PON-normalised — read the *spread*, not "
        "the offset."
    )
    go = _plotly()
    if go is None:
        return _no(suffix, title, caption, _NO_PLOTLY)
    needed = ["ultra_short", "core_short", "long", "ultra_long"]
    if fsc is None or any(c not in getattr(fsc, "columns", []) for c in needed):
        return _no(suffix, title, caption, "FSC size channels absent")

    import numpy as np

    short = pd.to_numeric(fsc["ultra_short"], errors="coerce") + pd.to_numeric(
        fsc["core_short"], errors="coerce"
    )
    long_ = pd.to_numeric(fsc["long"], errors="coerce") + pd.to_numeric(
        fsc["ultra_long"], errors="coerce"
    )
    keep = (short > 0) & (long_ > 0)
    if not keep.any():
        return _no(suffix, title, caption, "no bin has both short and long fragments")
    ratio = np.log2(short[keep] / long_[keep])
    if _is_constant(ratio):
        return _no(suffix, title, caption, "constant across bins — not plotted")

    fig = go.Figure(
        go.Histogram(
            x=ratio,
            nbinsx=90,
            marker=dict(color=ACCENT, opacity=0.8),
            hovertemplate="log2 %{x:.2f}<br>%{y:,} bins<extra></extra>",
        )
    )
    median = float(ratio.median())
    fig.add_vline(
        x=median,
        line=dict(color=OPPOSE, width=2),
        annotation_text=f"median {median:+.2f}",
        annotation_position="top",
        annotation_font_size=10,
    )
    fig.update_layout(
        **_layout(go, 280, xaxis_title="log2(short / long) per bin", yaxis_title="bins")
    )
    return Chart(suffix, title, caption, fig)


def mutant_fragment_sizes(mfsd: Optional[pd.DataFrame]) -> Chart:
    """Act 1.4 — the cleanest test there is, at single-variant resolution.

    A scatter rather than one row per variant, because the variant count is not
    bounded: a panel run may carry one somatic variant or several hundred, and
    a per-variant row would be a legible chart at n=3 and a 12,000-pixel one at
    n=300. One point per variant scales across that whole range and keeps the
    question — is ALT shorter than REF — readable as a single axis.
    """
    suffix, title = ".mFSD.parquet", "Mutant vs reference fragment size"
    caption = (
        "At a known somatic variant, ALT-carrying fragments are tumour DNA by "
        "definition and REF-carrying fragments are mostly normal, so comparing "
        "their sizes tests the size thesis without any genome-wide average. "
        "Each point is one variant; left of the line means its ALT fragments "
        "are shorter, which is the tumour direction. Hover for the locus."
    )
    go = _plotly()
    if go is None:
        return _no(suffix, title, caption, _NO_PLOTLY)
    if mfsd is None or mfsd.empty:
        return _no(
            suffix, title, caption, "mFSD absent — needs a variant list (--variants)"
        )

    need = {"REF_MeanSize", "ALT_MeanSize", "ALT_Count"}
    if not need <= set(mfsd.columns):
        return _no(
            suffix,
            title,
            caption,
            f"mFSD lacks {', '.join(sorted(need - set(mfsd.columns)))}",
        )

    work = mfsd.copy()
    for col in ("REF_MeanSize", "ALT_MeanSize", "ALT_Count"):
        work[col] = pd.to_numeric(work[col], errors="coerce")

    # A variant with no ALT fragment has ALT_MeanSize == 0. That is an absence
    # of measurement, not a measurement of zero, and plotting it would put a
    # fabricated point at the shortest possible size -- the most abnormal end
    # of the very axis being read.
    total = len(work)
    measurable = (work["ALT_Count"] > 0) & (work["REF_MeanSize"] > 0)
    unobserved = work[~measurable].copy()
    work = work[measurable].dropna(subset=["REF_MeanSize", "ALT_MeanSize"])
    if work.empty and unobserved.empty:
        return _no(suffix, title, caption, "no variant has usable size data")

    delta = work["ALT_MeanSize"] - work["REF_MeanSize"]

    locus: "pd.Series | pd.Index" = (
        work["Chrom"].astype(str) + ":" + work["Pos"].astype(str)
        if {"Chrom", "Pos"} <= set(work.columns)
        else work.index.astype(str)
    )

    # KS_Valid marks variants with enough fragments for the test to mean
    # anything. Shown as opacity rather than filtered out: a low-support
    # variant is still a data point, it just carries less weight.
    valid = (
        work["KS_Valid"].astype(bool)
        if "KS_Valid" in work.columns
        else pd.Series(True, index=work.index)
    )
    custom = list(
        zip(
            locus,
            work["ALT_MeanSize"].round(1),
            work["REF_MeanSize"].round(1),
            work["ALT_Count"].astype(int),
            ["yes" if v else "no" for v in valid],
        )
    )

    fig = go.Figure(
        go.Scatter(
            x=delta,
            y=work["ALT_Count"],
            mode="markers",
            marker=dict(
                color=[ACCENT if d < 0 else OPPOSE for d in delta],
                size=10,
                opacity=[0.9 if v else 0.35 for v in valid],
                line=dict(width=0),
            ),
            customdata=custom,
            hovertemplate=(
                "%{customdata[0]}<br>ALT %{customdata[1]} bp · "
                "REF %{customdata[2]} bp<br>difference %{x:+.1f} bp"
                "<br>%{customdata[3]:,} ALT fragments"
                "<br>KS valid: %{customdata[4]}<extra></extra>"
            ),
        )
    )
    # Variants with no ALT-supporting fragment are shown, not dropped -- their
    # absence is a finding. But they are placed at the origin with their own
    # symbol, never on the measurement scale: ALT_MeanSize is 0 for these, and
    # plotting 0 - REF would put a fabricated point at the far left, the most
    # tumour-like end of the very axis being read.
    if not unobserved.empty:
        un_locus: "pd.Series | pd.Index" = (
            unobserved["Chrom"].astype(str) + ":" + unobserved["Pos"].astype(str)
            if {"Chrom", "Pos"} <= set(unobserved.columns)
            else unobserved.index.astype(str)
        )
        fig.add_trace(
            go.Scatter(
                x=[0] * len(unobserved),
                y=[0] * len(unobserved),
                mode="markers",
                marker=dict(
                    color=MUTED,
                    size=11,
                    symbol="x-thin",
                    line=dict(width=1.6, color=MUTED),
                ),
                customdata=list(un_locus),
                hovertemplate=(
                    "%{customdata}<br><b>no ALT-supporting fragment</b>"
                    "<br>no size to compare — position is not a measurement"
                    "<extra></extra>"
                ),
            )
        )

    fig.add_vline(
        x=0,
        line=dict(color=MUTED, width=1),
        annotation_text="no difference",
        annotation_position="top",
        annotation_font_size=10,
    )
    fig.update_layout(
        **_layout(
            go,
            320,
            xaxis_title="ALT − REF mean fragment length (bp)",
            yaxis_title="ALT-supporting fragments",
        )
    )
    # Log scale would drop the zero-count markers entirely, so it is only used
    # when every variant was observed.
    if (
        unobserved.empty
        and len(work) > 1
        and work["ALT_Count"].max() > 10 * max(work["ALT_Count"].min(), 1)
    ):
        fig.update_yaxes(type="log")

    shorter = int((delta < 0).sum())
    caption += (
        f" {shorter} of {len(work)} measured variants have shorter ALT fragments."
    )
    if not unobserved.empty:
        caption += (
            f" The {len(unobserved)} grey cross{'es' if len(unobserved) > 1 else ''} "
            f"at the origin had no ALT-supporting fragment out of {total} total — "
            "shown because that is itself a finding, positioned off the "
            "measurement scale because there is no size to place them at."
        )
    if not valid.all():
        caption += (
            f" Faded points ({int((~valid).sum())}) lack the fragment support "
            "for a valid KS test."
        )
    return Chart(suffix, title, caption, fig)


def jagged_end_composition(one_mer: Optional[pd.DataFrame]) -> Chart:
    """Act 2.2 — cut chemistry, orthogonal to motif diversity."""
    suffix, title = ".EndMotif1mer.parquet", "Terminal base composition"
    caption = (
        "The single 5′ base of each fragment. DNASE1L3 leaves single-stranded "
        "overhangs; filling them in during library prep writes templated bases "
        "at the end, so this reads the chemistry of the cut rather than its "
        "diversity. The C fraction is the quantity of interest."
    )
    go = _plotly()
    if go is None:
        return _no(suffix, title, caption, _NO_PLOTLY)
    if one_mer is None or one_mer.empty:
        return _no(suffix, title, caption, "EndMotif1mer absent")

    cols = {c.lower(): c for c in one_mer.columns}
    base, frac = cols.get("base"), cols.get("fraction")
    if not (base and frac):
        return _no(suffix, title, caption, "no base/fraction columns")

    work = one_mer[[base, frac]].copy()
    work[frac] = pd.to_numeric(work[frac], errors="coerce")
    work = work.dropna()
    if work.empty:
        return _no(suffix, title, caption, "no finite fractions")

    fig = go.Figure(
        go.Bar(
            x=work[base].astype(str),
            y=work[frac],
            marker=dict(
                color=[ACCENT if str(b).upper() == "C" else MUTED for b in work[base]]
            ),
            hovertemplate="%{x}<br>%{y:.4f}<extra></extra>",
        )
    )
    fig.add_hline(
        y=0.25,
        line=dict(color=MUTED, width=1, dash="dash"),
        annotation_text="even",
        annotation_position="right",
        annotation_font_size=10,
    )
    fig.update_layout(
        **_layout(go, 220, xaxis_title="terminal base", yaxis_title="fraction")
    )
    return Chart(suffix, title, caption, fig)


def gene_level_mds(mds_gene: Optional[pd.DataFrame]) -> Chart:
    """Act 2.3 — localise aberrant cutting to a driver."""
    suffix, title = ".MDS.gene.parquet", "Motif diversity per gene"
    caption = (
        "Where global MDS says cutting is abnormal, this says where. `mds_e1` "
        "is the first exon specifically — the promoter's nucleosome-depleted "
        "region, where accessibility differences are largest. Lower is the "
        "abnormal direction."
    )
    go = _plotly()
    if go is None:
        return _no(suffix, title, caption, _NO_PLOTLY)
    if mds_gene is None or "mds_mean" not in getattr(mds_gene, "columns", []):
        return _no(suffix, title, caption, "MDS.gene absent")

    work = mds_gene.copy()
    work["mds_mean"] = pd.to_numeric(work["mds_mean"], errors="coerce")
    work = work.dropna(subset=["mds_mean"])
    if work.empty:
        return _no(suffix, title, caption, "no finite mds_mean values")
    if _is_constant(work["mds_mean"]):
        return _no(suffix, title, caption, "mds_mean is constant — not plotted")

    work = (
        work.sort_values("mds_mean").head(25).sort_values("mds_mean", ascending=False)
    )
    fig = go.Figure()
    fig.add_bar(
        x=work["mds_mean"],
        y=work["gene"].astype(str),
        orientation="h",
        marker=dict(color=ACCENT, opacity=0.85),
        name="gene mean",
        hovertemplate="%{y}<br>gene mean %{x:.4f}<extra></extra>",
    )
    if "mds_e1" in work.columns:
        e1 = pd.to_numeric(work["mds_e1"], errors="coerce")
        drawn = e1.notna()
        if drawn.any():
            fig.add_trace(
                go.Scatter(
                    x=e1[drawn],
                    y=work.loc[drawn, "gene"].astype(str),
                    mode="markers",
                    name="E1",
                    marker=dict(color=OPPOSE, size=8, symbol="diamond"),
                    hovertemplate="%{y}<br>E1 %{x:.4f}<extra></extra>",
                )
            )
            caption += (
                f" Diamonds mark `mds_e1`, present for {int(drawn.sum())} of "
                f"{len(work)} shown — NaN elsewhere means the gene has no "
                "captured first exon, which is not the same as zero."
            )
    fig.update_layout(
        **_layout(
            go, max(220, 22 * len(work) + 70), xaxis_title="motif diversity score"
        )
    )
    fig.update_yaxes(automargin=True)
    return Chart(suffix, title, caption + " Lowest 25 genes by mean.", fig)


#: Base pairs per stored WPS bin.
#:
#: Every WPS profile in the output -- foreground anchors, panel anchors and the
#: Alu background -- is 200 bins over a 2000 bp window (`rust/src/wps.rs`:
#: `NUM_BINS`, `WINDOW_SIZE`, and "2000bp profile / 10bp bins = 200 bins" for
#: the background). The plots below used the bin index directly as the x value
#: while labelling the axis "bp", which understated every distance tenfold: a
#: TSS dip at bin +15 read as +15 bp when it is +150 bp, and a +-1000 bp window
#: read as +-100 bp.
#:
#: Pinned in `validate/claims.py` against the Rust constants, so a change to
#: the window or the bin count fails a test rather than silently rescaling
#: every WPS figure.
WPS_BIN_BP = 10


def nucleosome_profile(bg: Optional[pd.DataFrame]) -> Chart:
    """Act 3 — reinstated. Every derived WPS metric was constant before 0.9.0."""
    suffix, title = ".WPS_background.parquet", "Nucleosome protection over Alu"
    caption = (
        "WPS stacked over Alu elements: a genome-wide chromatin readout. The "
        "periodicity of this wave gives the nucleosome repeat length. This "
        "panel was withdrawn from earlier reporting because every derived "
        "metric was structurally constant; the underlying stacking was always "
        "sound, and the derived values are data-dependent as of 0.9.0."
    )
    go = _plotly()
    if go is None:
        return _no(suffix, title, caption, _NO_PLOTLY)
    if bg is None or "stacked_wps_nuc" not in getattr(bg, "columns", []):
        return _no(suffix, title, caption, "WPS_background absent")

    row = None
    for _, candidate in bg.iterrows():
        value = candidate["stacked_wps_nuc"]
        if value is not None and len(value) > 1:
            row = candidate
            break
    if row is None:
        return _no(suffix, title, caption, "no stacked profile has any positions")

    profile = list(row["stacked_wps_nuc"])
    centre = len(profile) // 2
    x = [(i - centre) * WPS_BIN_BP for i in range(len(profile))]
    group = str(row.get("group_id", "Alu"))

    fig = go.Figure(
        go.Scatter(
            x=x,
            y=profile,
            mode="lines",
            line=dict(color=ACCENT, width=1.4),
            hovertemplate="%{x} bp from centre<br>WPS %{y:.2f}<extra></extra>",
        )
    )
    fig.add_hline(y=0, line=dict(color=MUTED, width=1))
    nrl = row.get("nrl_bp")
    censored = bool(row.get("nrl_at_band_limit", False))
    if nrl is not None and not censored:
        caption += f" Repeat length for `{group}`: {float(nrl):.0f} bp."
    elif censored:
        caption += (
            f" `{group}` is flagged `nrl_at_band_limit` — the spectral peak sat "
            "on the edge of the search band, so no repeat length was measured."
        )
    fig.update_layout(
        **_layout(
            go,
            260,
            xaxis_title="position relative to Alu centre (bp)",
            yaxis_title="stacked WPS",
        )
    )
    return Chart(suffix, title, caption, fig)


def methylation_composition(uxm: Optional[pd.DataFrame]) -> Chart:
    """Act 4.4 — cell-type fingerprint, when a bisulfite BAM was supplied."""
    suffix, title = ".UXM.parquet", "Methylation composition"
    caption = (
        "Every fragment spanning a marker region is classified Unmethylated, "
        "miXed or Methylated. Methylation is the most cell-type-specific mark "
        "in the genome, so the U/M pattern is a tissue fingerprint."
    )
    go = _plotly()
    if go is None:
        return _no(suffix, title, caption, _NO_PLOTLY)
    if uxm is None or uxm.empty:
        return _no(
            suffix,
            title,
            caption,
            "UXM absent — needs a methylation-aware BAM (--bisulfite-bam)",
        )

    fractions = [c for c in ("U", "X", "M") if c in uxm.columns]
    if not fractions:
        return _no(suffix, title, caption, "no U/X/M columns")

    means = {c: float(pd.to_numeric(uxm[c], errors="coerce").mean()) for c in fractions}
    if "X" in means and means["X"] == 0:
        caption += (
            " The X class is structurally zero as shipped — the classifier is "
            "threshold-parameterised and correct, but the call site collapses "
            "it. Read this as a binary U/M split."
        )
    fig = go.Figure(
        go.Bar(
            x=list(means.values()),
            y=list(means.keys()),
            orientation="h",
            marker=dict(color=[ACCENT, MUTED, OPPOSE][: len(means)]),
            hovertemplate="%{y}<br>mean fraction %{x:.3f}<extra></extra>",
        )
    )
    fig.update_layout(**_layout(go, 190, xaxis_title="mean fraction across markers"))
    return Chart(suffix, title, caption, fig)


def fsd_region_heatmap(fsd: Optional[pd.DataFrame]) -> Chart:
    """Per-arm size distributions, which the summed density curve hides.

    Taken from the existing output EDA notebook, whose per-region heatmap is
    the one thing a genome-wide density trace cannot show: whether the size
    shift is uniform or concentrated on particular arms.
    """
    suffix, title = ".FSD.parquet", "Fragment size by region"
    caption = (
        "Each row is a chromosome arm, each column a 5 bp size bin, normalised "
        "within the row so arms of different depth are comparable. The summed "
        "curve above answers *what* the size distribution is; this answers "
        "whether it is the same everywhere."
    )
    go = _plotly()
    if go is None:
        return _no(suffix, title, caption, _NO_PLOTLY)
    if fsd is None or fsd.empty or "region" not in fsd.columns:
        return _no(suffix, title, caption, "FSD table absent or has no region column")

    bins = sorted(
        (int(m.group(1)), str(c))
        for c in fsd.columns
        if (m := re.fullmatch(r"(\d+)-(\d+)", str(c)))
    )
    if len(bins) < 2:
        return _no(suffix, title, caption, "fewer than two size-bin columns")
    if len(fsd) < 2:
        return _no(suffix, title, caption, "only one region — nothing to compare")

    labels = [b for _, b in bins]
    matrix = fsd[labels].apply(pd.to_numeric, errors="coerce")
    totals = matrix.sum(axis=1)
    usable = totals > 0
    if not usable.any():
        return _no(suffix, title, caption, "every region is empty")
    # Row-normalise: without it the heatmap shows sequencing depth per arm,
    # not the size distribution, and the two look identical at a glance.
    fractions = matrix[usable].div(totals[usable], axis=0)

    fig = go.Figure(
        go.Heatmap(
            z=fractions.values,
            x=[b for b, _ in bins],
            y=fsd.loc[usable, "region"].astype(str),
            colorscale="Reds",
            hovertemplate="%{y}<br>%{x} bp<br>%{z:.4f} of that region<extra></extra>",
            colorbar=dict(title=dict(text="fraction", side="right"), thickness=10),
        )
    )
    fig.update_layout(
        **_layout(
            go,
            max(260, 13 * int(usable.sum()) + 90),
            xaxis_title="fragment length (bp)",
            showlegend=False,
        )
    )
    fig.update_xaxes(range=[min(b for b, _ in bins), 600])
    fig.update_yaxes(automargin=True, tickfont=dict(size=9))
    dropped = int((~usable).sum())
    if dropped:
        caption += f" {dropped} region(s) with no fragments are omitted."
    return Chart(suffix, title, caption, fig)


def _translucent(hex_colour: str, alpha: float) -> str:
    """`#ef5552` -> `rgba(239,85,82,0.15)`.

    Plotly rejects 8-digit hex outright, so the alpha has to be carried in
    rgba() rather than appended to the colour constant.
    """
    h = hex_colour.lstrip("#")
    r, g, b = (int(h[i : i + 2], 16) for i in (0, 2, 4))
    return f"rgba({r},{g},{b},{alpha})"


def _capture_note(frame: "pd.DataFrame") -> str:
    """State how much of the profile sits on capture bait, and per anchor.

    `capture_mask` is 1 only for bins inside a bait and away from its edges
    (`rust/src/wps.rs`). Nothing in the Python package read it before this.

    Reported per *anchor* as well as per bin, because the bin fraction alone
    misleads in both directions. A WPS window is 2000 bp and a bait is of order
    100 bp, so even a perfectly bait-centred anchor can only ever have a small
    minority of its bins captured -- on one XS2 sample the best anchor reached
    103 bins of 200 and the median among capturing anchors was 7. A low bin
    percentage is therefore the expected geometry, not evidence that the
    anchors were badly chosen, and an earlier draft of this note said "not the
    assay's intended capture", which was simply wrong for the panel set: those
    anchors are exactly what the assay targets.

    What the anchor count does distinguish is real. On that sample 5 of 166
    panel anchors touched a bait against 33 of 58,522 genome-wide -- the panel
    set is some thirty times more bait-proximal, and both curves still rest
    mostly on off-target coverage.

    Returns "" when the column is absent, and when every position is captured:
    a WGS run has no baits, so the whole notion does not apply and a note would
    be noise.
    """
    if "capture_mask" not in getattr(frame, "columns", []):
        return ""
    import numpy as np

    masks = [np.asarray(m, dtype=float) for m in frame["capture_mask"] if m is not None]
    if not masks:
        return ""
    stacked = np.vstack(masks)
    n_bins, n_total = int(stacked.sum()), int(stacked.size)
    if n_bins == n_total:
        return ""
    n_anchors = int((stacked.sum(axis=1) > 0).sum())

    # Counts, not a bare percentage: 506 of 11,704,400 renders as "0.00%" under
    # every fixed-decimal format and then reads as rounding rather than a fact.
    pct = 100.0 * n_bins / n_total if n_total else 0.0
    # States what was measured and the geometry that produced it, and stops
    # there. An earlier draft ended "read this as a fragmentomic profile, not a
    # targeted measurement", which is a verdict on usability -- not this
    # toolkit's call to make. Measure it, report what it rests on, and let the
    # reader and the downstream consumer decide what it supports.
    return (
        f" Bait coverage: {n_anchors:,} of {len(stacked):,} anchors have any "
        f"bait-covered position, {n_bins:,} of {n_total:,} bins ({pct:.3g}%). "
        "A WPS window spans 2000 bp and a capture bait is typically far "
        "narrower, so most of the window lies outside bait regardless."
    )


def _anchor_profile_chart(
    wps: "Optional[pd.DataFrame]", suffix: str, title: str, caption: str
) -> Chart:
    """Mean WPS around anchors, grouped by `region_type`.

    Shared by the genome-wide and panel anchor sets: same columns, same
    binning, same degeneracy rule, so a fix to one is a fix to both.
    """
    go = _plotly()
    if go is None:
        return _no(suffix, title, caption, _NO_PLOTLY)
    if wps is None or "wps_nuc" not in getattr(wps, "columns", []):
        return _no(suffix, title, caption, "WPS absent")

    import numpy as np

    caption += _capture_note(wps)

    groups = (
        wps.groupby("region_type") if "region_type" in wps.columns else [("all", wps)]
    )
    fig = go.Figure()
    drawn_any = False
    for i, (name, frame) in enumerate(groups):
        stacked = np.array(
            [
                np.asarray(v, dtype=float)
                for v in frame["wps_nuc"]
                if v is not None and len(v)
            ]
        )
        if stacked.size == 0:
            continue
        mean = np.nanmean(stacked, axis=0)
        if np.allclose(mean, mean[0]):
            continue
        centre = len(mean) // 2
        x = [(j - centre) * WPS_BIN_BP for j in range(len(mean))]

        # A +-1 SEM band, so the reader can see how much the mean is worth
        # instead of being told. It is the difference between a curve over
        # 34,851 anchors and one over 105, and describing that in prose invites
        # the caption to draw the conclusion for them.
        if len(stacked) > 1:
            sem = np.nanstd(stacked, axis=0, ddof=1) / np.sqrt(len(stacked))
            fig.add_trace(
                go.Scatter(
                    x=x + x[::-1],
                    y=list(mean + sem) + list((mean - sem)[::-1]),
                    fill="toself",
                    fillcolor=_translucent(ACCENT if i == 0 else OPPOSE, 0.15),
                    line=dict(width=0),
                    hoverinfo="skip",
                    showlegend=False,
                    name=f"{str(name)} +-1 SEM",
                )
            )
        fig.add_trace(
            go.Scatter(
                x=x,
                y=mean,
                mode="lines",
                # str() first, here as in the tooltip below: the groupby key
                # is Any, and interpolating bytes renders b'TSS'. The tooltip
                # already did this; the legend was added later and repeated
                # the mistake the comment was written about.
                name=f"{str(name)} (n={len(stacked)})",
                line=dict(color=ACCENT if i == 0 else OPPOSE, width=1.8),
                hovertemplate=(
                    f"{str(name)}<br>%{{x}} bp<br>WPS %{{y:.2f}}<extra></extra>"
                ),
            )
        )
        drawn_any = True
    if not drawn_any:
        return _no(suffix, title, caption, "no anchor group has a varying profile")

    fig.add_hline(y=0, line=dict(color=MUTED, width=1))
    fig.add_vline(
        x=0,
        line=dict(color=MUTED, width=1, dash="dash"),
        annotation_text="anchor",
        annotation_position="top",
        annotation_font_size=10,
    )
    fig.update_layout(
        **_layout(
            go,
            300,
            xaxis_title="position relative to anchor (bp)",
            yaxis_title="mean WPS",
            showlegend=True,
            legend=dict(orientation="h", y=1.12, x=0),
        ),
    )
    return Chart(suffix, title, caption, fig)


def wps_anchor_profile(wps: Optional[pd.DataFrame]) -> Chart:
    """Act 3.1 — the foreground WPS profile, averaged over anchors.

    The audit withdrew every *derived* WPS metric because each was constant.
    The stacking itself was always sound, and this is the panel that was
    explicitly still worth plotting.
    """
    return _anchor_profile_chart(
        wps,
        ".WPS.parquet",
        "Protection around TSS and CTCF",
        "Windowed protection score averaged across genome-wide anchors, by "
        "anchor type. A positioned nucleosome shields the DNA beneath it, so "
        "the dip at the anchor traces the nucleosome-depleted region at an "
        "active promoter, and the flanking ridges are the nucleosomes either "
        "side of it.",
    )


def wps_panel_profile(wps: Optional[pd.DataFrame]) -> Chart:
    """The same profile over the panel's own anchors.

    Worth its own figure rather than a footnote. The genome-wide set is
    dominated by off-bait positions on a targeted assay, so its curve is an
    off-target readout; these anchors are the ones the panel was designed
    around, and on the sample this was written against they carry far more
    amplitude (262 and 199, against 1.6 and 14.7) from a far smaller n
    (105 and 61). Two different measurements, and previously only one was
    plotted while both were written.
    """
    return _anchor_profile_chart(
        wps,
        ".WPS.panel.parquet",
        "Protection around the panel's anchors",
        "As the genome-wide profile, but over the anchors the panel targets. "
        "Far fewer anchors, so the mean is noisier, but they are the positions "
        "the assay was designed to cover.",
    )


def gc_correction_curve(factors: Optional[pd.DataFrame]) -> Chart:
    """Is the GC model sane? Every corrected count depends on the answer."""
    suffix, title = ".correction_factors.parquet", "GC correction factors"
    caption = (
        "The weight applied to each fragment before it is counted, by GC "
        "content. Everything GC-corrected downstream — FSC, gene-level "
        "coverage — is multiplied by this, so a wild curve propagates "
        "everywhere. Values near 1 mean little bias to correct."
    )
    go = _plotly()
    if go is None:
        return _no(suffix, title, caption, _NO_PLOTLY)
    if factors is None or "correction_factor" not in getattr(factors, "columns", []):
        return _no(suffix, title, caption, "correction factors absent")
    if "gc_percent" not in factors.columns:
        return _no(suffix, title, caption, "no gc_percent column")

    work = (
        factors[["gc_percent", "correction_factor"]]
        .apply(pd.to_numeric, errors="coerce")
        .dropna()
    )
    if work.empty:
        return _no(suffix, title, caption, "no finite correction factors")
    if _is_constant(work["correction_factor"]):
        value = work["correction_factor"].iloc[0]
        return _no(
            suffix,
            title,
            caption,
            f"every factor is {value:.4g} — no GC correction was fitted",
        )

    summary = work.groupby("gc_percent")["correction_factor"].agg(
        ["median", "min", "max"]
    )
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=list(summary.index) + list(summary.index[::-1]),
            y=list(summary["max"]) + list(summary["min"][::-1]),
            fill="toself",
            fillcolor="rgba(239,85,82,0.13)",
            line=dict(width=0),
            hoverinfo="skip",
            showlegend=False,
        )
    )
    fig.add_trace(
        go.Scatter(
            x=summary.index,
            y=summary["median"],
            mode="lines",
            line=dict(color=ACCENT, width=2),
            hovertemplate="GC %{x}%<br>median factor %{y:.3f}<extra></extra>",
        )
    )
    fig.add_hline(
        y=1.0,
        line=dict(color=MUTED, width=1, dash="dash"),
        annotation_text="no correction",
        annotation_position="right",
        annotation_font_size=10,
    )
    fig.update_layout(
        **_layout(
            go, 260, xaxis_title="GC content (%)", yaxis_title="correction factor"
        )
    )
    caption += " The band spans the range across fragment-length bins."
    return Chart(suffix, title, caption, fig)


def _first_present(tables: Dict, *suffixes: str):
    """The first suffix that resolved to a frame, and which one it was.

    Explicitly not ``a or b``: a DataFrame has no truth value, so that raises
    rather than falling through — and an *empty* frame is a real answer here,
    distinct from an absent one.
    """
    for suffix in suffixes:
        df = tables.get(suffix)
        if df is not None:
            return suffix, df
    return suffixes[0], None


def build_charts(tables: Dict) -> List[Chart]:
    """Every chart, each tagged with the table it belongs beside."""
    g = tables.get
    ocf_suffix, ocf = _first_present(tables, *OCF_PREFERENCE)
    # Ordered as the five acts of the audit notebook's report structure:
    # size, cutting, nucleosome positioning, tissue and accessibility.
    return [
        # Act 1 -- fragment size
        fragment_size_density(g(".FSD.parquet")),
        fsd_region_heatmap(g(".FSD.parquet")),
        size_channel_composition(g(".FSC.parquet")),
        fsc_short_long_by_bin(g(".FSC.parquet")),
        short_long_spread(g(".FSR.parquet")),
        mutant_fragment_sizes(g(".mFSD.parquet")),
        # Act 2 -- nuclease cutting
        end_motif_spectrum(g(".EndMotif.parquet")),
        jagged_end_composition(g(".EndMotif1mer.parquet")),
        gene_level_mds(g(".MDS.gene.parquet")),
        # Act 3 -- nucleosome positioning. Withdrawn from earlier reporting
        # because every derived metric was constant; reinstated in 0.9.0.
        wps_anchor_profile(g(".WPS.parquet")),
        wps_panel_profile(g(".WPS.panel.parquet")),
        nucleosome_profile(g(".WPS_background.parquet")),
        # Act 4 -- tissue and accessibility
        tissue_shedding(ocf_suffix, ocf),
        accessibility_by_cancer_type(g(".ATAC.parquet")),
        top_transcription_factors(g(".TFBS.parquet")),
        methylation_composition(g(".UXM.parquet")),
        # Diagnostic: everything GC-corrected depends on this curve.
        gc_correction_curve(g(".correction_factors.parquet")),
    ]


def charts_by_suffix(charts: Sequence[Chart]) -> Dict[str, List[Chart]]:
    """Group charts under the table they explain, so they render together."""
    out: Dict[str, List[Chart]] = {}
    for chart in charts:
        out.setdefault(chart.suffix, []).append(chart)
    return out
