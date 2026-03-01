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

logger = logging.getLogger("core.fsr_processor")

# Using FSC channels for consistent terminology
CHANNELS = ["ultra_short", "core_short", "mono_nucl", "di_nucl", "long"]


def process_fsr(
    counts_df: pd.DataFrame,
    output_path: Path,
    windows: int = 100000,
    continue_n: int = 50,
    pon=None,
) -> Path:
    """
    Process fragment counts into FSR ratios.

    Takes raw bin-level counts and computes short/long ratios with proper
    normalization order:
    1. Aggregate into windows
    2. Normalize counts to PoN FIRST (if provided)
    3. THEN compute ratios from normalized values

    This order is critical for accurate cross-sample comparison.

    Args:
        counts_df: DataFrame with columns [chrom, start, end, ultra_short,
                   core_short, mono_nucl, di_nucl, long, total, mean_gc]
        output_path: Path to write FSR.tsv output
        windows: Window size in bp (default 100000)
        continue_n: Number of bins to aggregate (default 50)
        pon: Optional PON model for count normalization

    Returns:
        Path to the output file
    """
    logger.info(f"Processing FSR: {len(counts_df)} bins → {output_path}")

    # Validate required columns
    required_cols = ["chrom"] + CHANNELS + ["total"]
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
        window_starts = np.arange(n_windows) * continue_n * windows
        window_ends = (np.arange(n_windows) + 1) * continue_n * windows - 1

        # Aggregate all channels
        channel_sums = {}
        for ch in CHANNELS:
            mat = group[ch].values[:trunc_len].reshape(n_windows, continue_n)
            channel_sums[ch] = mat.sum(axis=1).astype(float)

        # For FSR, "short" = ultra_short + core_short (65-149bp)
        # and "long" = di_nucl + long (221-400bp)
        short_counts = channel_sums["ultra_short"] + channel_sums["core_short"]
        long_counts = channel_sums["di_nucl"] + channel_sums["long"]

        total_mat = group["total"].values[:trunc_len].reshape(n_windows, continue_n)
        total_counts = total_mat.sum(axis=1).astype(float)

        for i in range(n_windows):
            region = f"{chrom}:{window_starts[i]}-{window_ends[i]}"

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
        pd.DataFrame(
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
        ).to_csv(output_path, sep="\t", index=False)
        return output_path

    results_df = pd.DataFrame(results)
    results_df.to_csv(output_path, sep="\t", index=False, float_format="%.6f")

    logger.info(f"FSR complete: {len(results_df)} windows → {output_path}")
    return output_path
