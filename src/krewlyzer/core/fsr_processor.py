"""
FSR (Fragment Size Ratio) processor.

Computes short/long fragment ratios for cancer biomarker analysis.
Uses PoN normalization BEFORE ratio calculation for accurate comparison.
Used by both standalone fsr.py and run-all wrapper.py.
"""

from pathlib import Path
import pandas as pd
import numpy as np
import logging

from .output_utils import write_table  # TSV/Parquet unified writer

logger = logging.getLogger("core.fsr_processor")

# Using FSC channels for consistent terminology
CHANNELS = ["ultra_short", "core_short", "mono_nucl", "di_nucl", "long"]


def process_fsr(
    counts_df: pd.DataFrame,
    output_path: Path,
    windows: int = 100000,
    continue_n: int = 50,
    pon=None,
    output_format: str = "tsv",
    compress: bool = False,
) -> Path:
    """
    Process fragment counts into FSR ratios.

    Takes raw bin-level counts and computes short/long ratios with proper
    normalization order:
    1. Aggregate bins into windows (continue_n bins per window)
    2. Normalize counts to PoN FIRST (if provided)
    3. THEN compute ratios from normalized values

    This order is critical for accurate cross-sample comparison.

    Region labels use actual genomic coordinates from the input DataFrame
    (matching the FSC approach), not a synthetic formula.

    Args:
        counts_df: DataFrame with columns [chrom, start, end, ultra_short,
                   core_short, mono_nucl, di_nucl, long, total, mean_gc]
        output_path: Path to write FSR output
        windows: Retained for API compatibility; coordinates are now read
            from the DataFrame directly.
        continue_n: Number of bins to aggregate per window (default 50).
            Panel mode typically uses 1 for per-bin resolution.
        pon: Optional PON model for count normalization
        output_format: Output format ('tsv' or 'parquet')
        compress: Whether to gzip-compress TSV output

    Returns:
        Path to the output file
    """
    logger.info(f"Processing FSR: {len(counts_df)} bins → {output_path}")
    logger.debug(
        f"FSR params: windows={windows}, continue_n={continue_n}, "
        f"pon={'yes' if pon else 'no'}, format={output_format}"
    )

    # Validate required columns (start/end needed for region labels)
    required_cols = ["chrom", "start", "end"] + CHANNELS + ["total"]
    missing = [c for c in required_cols if c not in counts_df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    # Aggregate bins into windows
    results = []
    for chrom, group in counts_df.groupby("chrom", sort=False):
        n_bins = len(group)
        n_windows = n_bins // continue_n

        if n_windows == 0:
            continue

        trunc_len = n_windows * continue_n

        # Use actual genomic coordinates from the input bins.
        # Takes the start of the first bin and end of the last bin
        # in each aggregated window (matches FSC approach at
        # fsc_processor.py:L228-229).
        starts = (
            group["start"].to_numpy()[:trunc_len].reshape(n_windows, continue_n)[:, 0]
        )
        ends = group["end"].to_numpy()[:trunc_len].reshape(n_windows, continue_n)[:, -1]

        # Aggregate all channels
        channel_sums = {}
        for ch in CHANNELS:
            mat = group[ch].to_numpy()[:trunc_len].reshape(n_windows, continue_n)
            channel_sums[ch] = mat.sum(axis=1).astype(float)

        # For FSR, "short" = ultra_short + core_short (65-149bp)
        # and "long" = di_nucl + long (221-400bp)
        short_counts = channel_sums["ultra_short"] + channel_sums["core_short"]
        long_counts = channel_sums["di_nucl"] + channel_sums["long"]

        total_mat = group["total"].to_numpy()[:trunc_len].reshape(n_windows, continue_n)
        total_counts = total_mat.sum(axis=1).astype(float)

        for i in range(n_windows):
            region = f"{chrom}:{starts[i]}-{ends[i]}"

            s = short_counts[i]
            long_val = long_counts[i]
            t = total_counts[i]

            # CORRECT ORDER: Normalize to PoN FIRST
            if pon is not None:
                short_mean = pon.get_mean("short") if hasattr(pon, "get_mean") else None
                long_mean = pon.get_mean("long") if hasattr(pon, "get_mean") else None

                if short_mean and short_mean > 0:
                    s_norm = s / short_mean
                else:
                    s_norm = s

                if long_mean and long_mean > 0:
                    l_norm = long_val / long_mean
                else:
                    l_norm = long_val
            else:
                s_norm = s
                l_norm = long_val

            # THEN compute ratio from normalized values
            if l_norm > 0:
                short_long_ratio = s_norm / l_norm
            else:
                short_long_ratio = s_norm if s_norm > 0 else 0.0

            # Log2 ratio for ML
            if short_long_ratio > 0:
                short_long_log2 = np.log2(short_long_ratio)
            else:
                short_long_log2 = 0.0

            # Fractions of total
            short_frac = s / t if t > 0 else 0.0
            long_frac = long_val / t if t > 0 else 0.0

            results.append(
                {
                    "region": region,
                    "short_count": int(s),
                    "long_count": int(long_val),
                    "total_count": int(t),
                    "short_norm": s_norm,
                    "long_norm": l_norm,
                    "short_long_ratio": short_long_ratio,
                    "short_long_log2": short_long_log2,
                    "short_frac": short_frac,
                    "long_frac": long_frac,
                }
            )

    if not results:
        logger.warning("No valid windows found for FSR")
        empty_df = pd.DataFrame(
            columns=[
                "region",
                "short_count",
                "long_count",
                "total_count",
                "short_norm",
                "long_norm",
                "short_long_ratio",
                "short_long_log2",
                "short_frac",
                "long_frac",
            ]
        )
        write_table(
            empty_df, output_path, output_format=output_format, compress=compress
        )
        return output_path

    results_df = pd.DataFrame(results)
    write_table(
        results_df,
        output_path,
        output_format=output_format,
        compress=compress,
        float_format="%.6f",
    )

    logger.info(f"FSR complete: {len(results_df)} windows → {output_path}")
    return output_path
