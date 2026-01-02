"""
WPS (Windowed Protection Score) processor.

Provides WPS-specific processing including:
- PON z-score overlay
- Savitzky-Golay smoothing for noise reduction
Used by both standalone wps.py and run-all wrapper.py.
"""

from pathlib import Path
from typing import Optional
import pandas as pd
import numpy as np
import logging

logger = logging.getLogger("core.wps_processor")

# Savitzky-Golay defaults
SAVGOL_WINDOW = 11  # Window length (must be odd)
SAVGOL_POLYORDER = 3  # Polynomial order


# =============================================================================
# Unified Post-Processing (shared by standalone wps and run-all)
# =============================================================================

def post_process_wps(
    wps_parquet: Path,
    wps_background_parquet: Optional[Path] = None,
    pon_baseline_parquet: Optional[Path] = None,
    smooth: bool = True,
    extract_periodicity: bool = True
) -> dict:
    """
    Unified WPS post-processing pipeline.
    
    Called by both standalone `krewlyzer wps` and `run-all` for consistent output.
    
    Steps:
    1. Savitzky-Golay smoothing (adds *_smooth columns)
    2. PoN subtraction if baseline provided (adds *_delta, *_z columns)
    3. FFT periodicity for background (adds nrl_* columns)
    
    Args:
        wps_parquet: Path to foreground WPS Parquet (e.g., sample.WPS.parquet)
        wps_background_parquet: Path to background WPS Parquet (e.g., sample.WPS_background.parquet)
        pon_baseline_parquet: Path to PoN baseline Parquet (optional)
        smooth: Apply Savitzky-Golay smoothing
        extract_periodicity: Extract FFT periodicity from background
        
    Returns:
        dict with processing summary and metrics
    """
    result = {
        "wps_parquet": str(wps_parquet),
        "smoothed": False,
        "pon_subtracted": False,
        "periodicity_extracted": False,
        "periodicity_score": None
    }
    
    # 1. Smooth foreground WPS
    if smooth and wps_parquet.exists():
        try:
            logger.info(f"Smoothing foreground WPS: {wps_parquet}")
            df = pd.read_parquet(wps_parquet)
            
            for col in ["wps_nuc", "wps_tf", "prot_frac_nuc", "prot_frac_tf"]:
                if col in df.columns:
                    smoothed_vectors = []
                    for vec in df[col]:
                        if vec is not None and len(vec) > 0:
                            arr = np.array(vec, dtype=np.float32)
                            smoothed = savgol_smooth(arr)
                            smoothed_vectors.append(smoothed.tolist())
                        else:
                            smoothed_vectors.append(vec)
                    df[f"{col}_smooth"] = smoothed_vectors
            
            df.to_parquet(wps_parquet, index=False)
            result["smoothed"] = True
            logger.info(f"Added smooth columns to {wps_parquet}")
        except Exception as e:
            logger.warning(f"Failed to smooth foreground WPS: {e}")
    
    # 2. Smooth and add periodicity to background WPS
    if wps_background_parquet and wps_background_parquet.exists():
        try:
            logger.info(f"Processing background WPS: {wps_background_parquet}")
            df = pd.read_parquet(wps_background_parquet)
            
            # Smooth stacked profiles
            if smooth:
                for col in ["stacked_wps_nuc", "stacked_wps_tf"]:
                    if col in df.columns:
                        smoothed_vectors = []
                        for vec in df[col]:
                            if vec is not None and len(vec) > 0:
                                arr = np.array(vec, dtype=np.float32)
                                smoothed = savgol_smooth(arr)
                                smoothed_vectors.append(smoothed.tolist())
                            else:
                                smoothed_vectors.append(vec)
                        df[f"{col}_smooth"] = smoothed_vectors
            
            # Extract periodicity
            if extract_periodicity and "stacked_wps_nuc" in df.columns:
                periods, amplitudes, snrs, quality_scores = [], [], [], []
                for vec in df["stacked_wps_nuc"]:
                    if vec is not None and len(vec) > 0:
                        arr = np.array(vec, dtype=np.float64)
                        metrics = extract_periodicity_fft(arr, bin_size_bp=10.0)
                        periods.append(metrics["period_bp"])
                        amplitudes.append(metrics["amplitude"])
                        snrs.append(metrics["snr"])
                        quality_scores.append(metrics["quality_score"])
                    else:
                        periods.append(0.0)
                        amplitudes.append(0.0)
                        snrs.append(0.0)
                        quality_scores.append(0.0)
                
                df["nrl_period_bp"] = periods
                df["nrl_amplitude"] = amplitudes
                df["nrl_snr"] = snrs
                df["nrl_quality"] = quality_scores
                result["periodicity_extracted"] = True
                
                # Get Global_All score
                global_row = df[df.get("group_id", "") == "Global_All"]
                if len(global_row) > 0:
                    result["periodicity_score"] = global_row.iloc[0].get("nrl_quality", 0.0)
            
            df.to_parquet(wps_background_parquet, index=False)
            logger.info(f"Processed background WPS: {wps_background_parquet}")
        except Exception as e:
            logger.warning(f"Failed to process background WPS: {e}")
    
    # 3. PoN subtraction
    if pon_baseline_parquet and Path(pon_baseline_parquet).exists() and wps_parquet.exists():
        try:
            logger.info(f"Subtracting PoN baseline: {pon_baseline_parquet}")
            subtract_pon_baseline(wps_parquet, pon_baseline_parquet, wps_parquet)
            result["pon_subtracted"] = True
        except Exception as e:
            logger.warning(f"Failed to subtract PoN baseline: {e}")
    
    return result


# =============================================================================
# Savitzky-Golay Smoothing for WPS Parquet
# =============================================================================

def savgol_smooth(
    profile: np.ndarray,
    window_length: int = SAVGOL_WINDOW,
    polyorder: int = SAVGOL_POLYORDER
) -> np.ndarray:
    """
    Apply Savitzky-Golay filter to a 1D profile.
    
    Args:
        profile: 1D numpy array (e.g., 200 bins of WPS values)
        window_length: Filter window length (must be odd, > polyorder)
        polyorder: Polynomial order for local fitting
        
    Returns:
        Smoothed profile of same length
    """
    from scipy.signal import savgol_filter
    
    if len(profile) < window_length:
        logger.warning(f"Profile length ({len(profile)}) < window ({window_length}), skipping")
        return profile
    
    return savgol_filter(profile, window_length, polyorder)


def smooth_wps_parquet(
    input_path: Path,
    output_path: Optional[Path] = None,
    window_length: int = SAVGOL_WINDOW,
    polyorder: int = SAVGOL_POLYORDER,
    columns: list = None
) -> pd.DataFrame:
    """
    Smooth WPS vector columns in a Parquet file.
    
    Args:
        input_path: Path to WPS Parquet (foreground or background)
        output_path: Optional output path (defaults to input with .smoothed suffix)
        window_length: Savitzky-Golay window
        polyorder: Polynomial order
        columns: List of vector columns to smooth (auto-detected if None)
        
    Returns:
        DataFrame with smoothed profiles
    """
    if columns is None:
        columns = ["wps_nuc", "wps_tf", "prot_frac_nuc", "prot_frac_tf"]
    
    logger.info(f"Loading WPS Parquet: {input_path}")
    df = pd.read_parquet(input_path)
    
    n_smoothed = 0
    for col in columns:
        if col not in df.columns:
            continue
            
        logger.info(f"Smoothing column: {col}")
        
        # Apply smoothing to each row's vector
        smoothed_vectors = []
        for vec in df[col]:
            if vec is not None and len(vec) > 0:
                arr = np.array(vec, dtype=np.float32)
                smoothed = savgol_smooth(arr, window_length, polyorder)
                smoothed_vectors.append(smoothed.tolist())
            else:
                smoothed_vectors.append(vec)
        
        df[f"{col}_smooth"] = smoothed_vectors
        n_smoothed += 1
    
    if output_path is None:
        output_path = input_path.with_suffix('.smoothed.parquet')
    
    logger.info(f"Writing smoothed Parquet: {output_path}")
    df.to_parquet(output_path, index=False)
    
    logger.info(f"Smoothed {n_smoothed} columns, {len(df)} rows")
    return df


def smooth_background_parquet(
    input_path: Path,
    output_path: Optional[Path] = None,
    window_length: int = SAVGOL_WINDOW,
    polyorder: int = SAVGOL_POLYORDER
) -> pd.DataFrame:
    """
    Smooth background WPS Parquet (hierarchical Alu stacking).
    """
    return smooth_wps_parquet(
        input_path,
        output_path,
        window_length,
        polyorder,
        columns=["stacked_wps_nuc", "stacked_wps_tf"]
    )


# =============================================================================
# PoN Subtraction for WPS Parquet
# =============================================================================

def subtract_pon_baseline(
    sample_parquet: Path,
    pon_baseline_parquet: Path,
    output_path: Optional[Path] = None,
    columns: list = None
) -> pd.DataFrame:
    """
    Subtract PoN (Panel of Normals) baseline from sample WPS profiles.
    
    For each region, computes:
    - delta_* = sample - pon_mean (raw difference)
    - z_* = (sample - pon_mean) / pon_std (z-score)
    
    Args:
        sample_parquet: Path to sample WPS Parquet
        pon_baseline_parquet: Path to PoN baseline Parquet (aggregated healthy samples)
        output_path: Output path (defaults to .pon_subtracted.parquet)
        columns: Vector columns to process (default: wps_nuc, wps_tf)
        
    Returns:
        DataFrame with PoN-subtracted profiles
    """
    if columns is None:
        columns = ["wps_nuc", "wps_tf"]
    
    logger.info(f"Loading sample: {sample_parquet}")
    sample_df = pd.read_parquet(sample_parquet)
    
    logger.info(f"Loading PoN baseline: {pon_baseline_parquet}")
    pon_df = pd.read_parquet(pon_baseline_parquet)
    
    # Build PoN lookup by region_id
    pon_lookup = {}
    region_id_col = "region_id" if "region_id" in pon_df.columns else "group_id"
    
    for _, row in pon_df.iterrows():
        region_id = row[region_id_col]
        pon_lookup[region_id] = {
            col: {
                "mean": np.array(row.get(f"{col}_mean", row.get(col, [])), dtype=np.float32),
                "std": np.array(row.get(f"{col}_std", [1.0] * len(row.get(col, []))), dtype=np.float32)
            }
            for col in columns if col in pon_df.columns or f"{col}_mean" in pon_df.columns
        }
    
    logger.info(f"PoN baseline contains {len(pon_lookup)} regions")
    
    # Process each sample row
    for col in columns:
        if col not in sample_df.columns:
            continue
        
        delta_col = []
        z_col = []
        
        for idx, row in sample_df.iterrows():
            region_id = row.get("region_id", row.get("group_id", str(idx)))
            sample_vec = np.array(row[col], dtype=np.float32)
            
            if region_id in pon_lookup and col in pon_lookup[region_id]:
                pon_mean = pon_lookup[region_id][col]["mean"]
                pon_std = pon_lookup[region_id][col]["std"]
                
                # Ensure same length
                if len(pon_mean) == len(sample_vec):
                    delta = sample_vec - pon_mean
                    # Avoid division by zero
                    pon_std_safe = np.where(pon_std > 1e-6, pon_std, 1.0)
                    z = delta / pon_std_safe
                else:
                    delta = sample_vec
                    z = np.zeros_like(sample_vec)
            else:
                delta = sample_vec
                z = np.zeros_like(sample_vec)
            
            delta_col.append(delta.tolist())
            z_col.append(z.tolist())
        
        sample_df[f"{col}_delta"] = delta_col
        sample_df[f"{col}_z"] = z_col
        logger.info(f"Added {col}_delta and {col}_z columns")
    
    if output_path is None:
        output_path = sample_parquet.with_suffix('.pon_subtracted.parquet')
    
    logger.info(f"Writing PoN-subtracted Parquet: {output_path}")
    sample_df.to_parquet(output_path, index=False)
    
    return sample_df



# =============================================================================
# FFT Periodicity Extraction
# =============================================================================

def extract_periodicity_fft(
    profile: np.ndarray,
    bin_size_bp: float = 10.0,
    min_period_bp: float = 150.0,
    max_period_bp: float = 250.0
) -> dict:
    """
    Extract nucleosome repeat length using FFT.
    
    Nucleosomes are spaced ~190bp apart. This function finds the dominant
    periodicity in the WPS signal using Fast Fourier Transform.
    
    Args:
        profile: 1D numpy array of WPS values (e.g., 200 bins)
        bin_size_bp: Size of each bin in base pairs (default: 10bp)
        min_period_bp: Minimum period to search (default: 150bp)
        max_period_bp: Maximum period to search (default: 250bp)
        
    Returns:
        dict with:
            - period_bp: Dominant period in base pairs
            - amplitude: FFT amplitude at dominant frequency
            - snr: Signal-to-noise ratio (peak vs background)
            - quality_score: 0-1 score for periodicity quality
    """
    from scipy.fft import rfft, rfftfreq
    
    n = len(profile)
    if n < 10:
        return {"period_bp": 0.0, "amplitude": 0.0, "snr": 0.0, "quality_score": 0.0}
    
    # Detrend (remove DC offset and linear trend)
    x = np.arange(n)
    coeffs = np.polyfit(x, profile, 1)
    detrended = profile - np.polyval(coeffs, x)
    
    # Z-score normalize to ensure consistent FFT amplitude regardless of Alu count
    # This makes quality_score stable whether we stack 50K or 1M Alus
    std = np.std(detrended)
    if std > 1e-9:
        normalized = detrended / std
    else:
        # Flat profile (no signal)
        return {"period_bp": 0.0, "amplitude": 0.0, "snr": 0.0, "quality_score": 0.0}
    
    # Apply Hann window to reduce spectral leakage
    window = np.hanning(n)
    windowed = normalized * window
    
    # Compute FFT
    fft_vals = rfft(windowed)
    freqs = rfftfreq(n, d=bin_size_bp)  # Frequencies in cycles/bp
    
    # Convert to periods
    with np.errstate(divide='ignore', invalid='ignore'):
        periods = np.where(freqs > 0, 1.0 / freqs, np.inf)
    
    # Find peaks in the target range
    amplitudes = np.abs(fft_vals)
    
    # Mask for target period range
    mask = (periods >= min_period_bp) & (periods <= max_period_bp) & (freqs > 0)
    
    if not np.any(mask):
        return {"period_bp": 0.0, "amplitude": 0.0, "snr": 0.0, "quality_score": 0.0}
    
    # Find maximum amplitude in range
    masked_amps = np.where(mask, amplitudes, 0)
    peak_idx = np.argmax(masked_amps)
    peak_amplitude = amplitudes[peak_idx]
    peak_period = periods[peak_idx]
    
    # Compute SNR (peak vs median background)
    background = np.median(amplitudes[mask])
    snr = peak_amplitude / background if background > 0 else 0.0
    
    # Quality score: SNR-based, capped at 1.0
    # SNR > 3 indicates clear periodicity (empirical threshold for nucleosomes)
    quality_score = min(1.0, snr / 3.0)
    
    return {
        "period_bp": float(peak_period),
        "amplitude": float(peak_amplitude),
        "snr": float(snr),
        "quality_score": float(quality_score)
    }


def add_periodicity_to_parquet(
    input_path: Path,
    output_path: Optional[Path] = None,
    column: str = "stacked_wps_nuc",
    bin_size_bp: float = 10.0
) -> pd.DataFrame:
    """
    Add FFT periodicity metrics to a WPS Parquet file.
    
    Args:
        input_path: Path to WPS Parquet (typically background with stacked profiles)
        output_path: Output path (defaults to input with .periodicity suffix)
        column: Column name containing WPS vectors
        bin_size_bp: Bin size in base pairs
        
    Returns:
        DataFrame with added periodicity columns
    """
    logger.info(f"Extracting FFT periodicity from: {input_path}")
    df = pd.read_parquet(input_path)
    
    if column not in df.columns:
        logger.warning(f"Column {column} not found, skipping periodicity extraction")
        return df
    
    periods = []
    amplitudes = []
    snrs = []
    quality_scores = []
    
    for vec in df[column]:
        if vec is not None and len(vec) > 0:
            arr = np.array(vec, dtype=np.float64)
            result = extract_periodicity_fft(arr, bin_size_bp)
            periods.append(result["period_bp"])
            amplitudes.append(result["amplitude"])
            snrs.append(result["snr"])
            quality_scores.append(result["quality_score"])
        else:
            periods.append(0.0)
            amplitudes.append(0.0)
            snrs.append(0.0)
            quality_scores.append(0.0)
    
    df["nrl_period_bp"] = periods  # NRL = Nucleosome Repeat Length
    df["nrl_amplitude"] = amplitudes
    df["nrl_snr"] = snrs
    df["nrl_quality"] = quality_scores
    
    if output_path is None:
        output_path = input_path.with_suffix('.periodicity.parquet')
    
    logger.info(f"Writing periodicity Parquet: {output_path}")
    df.to_parquet(output_path, index=False)
    
    # Log summary stats
    mean_period = np.mean([p for p in periods if p > 0])
    mean_quality = np.mean(quality_scores)
    logger.info(f"Mean NRL: {mean_period:.1f}bp, Mean quality: {mean_quality:.3f}")
    
    return df


def compute_sample_periodicity_score(
    background_parquet: Path,
    expected_nrl: float = 190.0,
    tolerance_bp: float = 20.0
) -> dict:
    """
    Compute overall periodicity score for a sample from background Parquet.
    
    Args:
        background_parquet: Path to WPS_background.parquet
        expected_nrl: Expected nucleosome repeat length (default: 190bp)
        tolerance_bp: Allowed deviation from expected (default: 20bp)
        
    Returns:
        dict with sample-level periodicity metrics
    """
    df = pd.read_parquet(background_parquet)
    
    # Look for Global_All group
    global_row = df[df.get("group_id", df.get("region_group", "")) == "Global_All"]
    
    if len(global_row) == 0:
        logger.warning("Global_All group not found in background Parquet")
        return {"periodicity_score": 0.0, "nrl_bp": 0.0, "deviation_bp": 0.0}
    
    row = global_row.iloc[0]
    
    # Get stacked profile
    col = "stacked_wps_nuc" if "stacked_wps_nuc" in df.columns else None
    if col is None:
        return {"periodicity_score": 0.0, "nrl_bp": 0.0, "deviation_bp": 0.0}
    
    profile = np.array(row[col], dtype=np.float64)
    result = extract_periodicity_fft(profile, bin_size_bp=10.0)
    
    nrl = result["period_bp"]
    deviation = abs(nrl - expected_nrl)
    
    # Score based on quality and deviation from expected
    deviation_penalty = max(0.0, 1.0 - deviation / tolerance_bp)
    periodicity_score = result["quality_score"] * deviation_penalty
    
    return {
        "periodicity_score": periodicity_score,
        "nrl_bp": nrl,
        "deviation_bp": deviation,
        "snr": result["snr"],
        "amplitude": result["amplitude"]
    }
