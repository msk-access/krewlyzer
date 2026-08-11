"""Build a synthetic output cohort that satisfies the contract.

Driven from ``CONTRACT`` itself rather than a hand-written file list, so a new
TableRule is covered the moment it is declared. Per-suffix builders exist only
where a domain check constrains the contents (motif alphabets, chr prefixes,
sums that must land on 1).

Values are varied by ``sample_idx`` so a well-formed cohort is *not* degenerate;
tests that want degeneracy break that on purpose.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd

from krewlyzer.validate.contract import CONTRACT, Kind, TableRule

CHROMS = [f"chr{i}" for i in range(1, 6)]


def _acgt_4mers() -> List[str]:
    bases = "ACGT"
    return [a + b + c + d for a in bases for b in bases for c in bases for d in bases]


def _normalised(n: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    raw = rng.random(n) + 0.05
    return np.round(raw / raw.sum(), 6)


def _rows_for(rule: TableRule) -> int:
    if rule.rows.exactly is not None:
        return rule.rows.exactly
    return max(rule.rows.at_least or 1, 4)


def build_frame(
    rule: TableRule, sample_idx: int, pon_applied: bool = True
) -> pd.DataFrame:
    n = _rows_for(rule)
    rng = np.random.default_rng(1000 + sample_idx)
    suffix = rule.suffix
    offset = sample_idx * 0.017

    if "EndMotif1mer" in suffix:
        return pd.DataFrame(
            {
                "base": list("ACGT"),
                "count": [10, 20, 30, 40],
                "fraction": _normalised(4, sample_idx),
            }
        )
    if "EndMotif" in suffix or "BreakPointMotif" in suffix:
        return pd.DataFrame(
            {"Motif": _acgt_4mers(), "Frequency": _normalised(256, sample_idx)}
        )
    if suffix.startswith(".FSD"):
        frame = {"region": [f"{c}:0-100000" for c in CHROMS[:n]]}
        bins = {f"{65 + 5 * i}-{69 + 5 * i}": rng.random(n) * 10 for i in range(6)}
        frame.update(bins)
        frame["total"] = np.sum(list(bins.values()), axis=0)
        if pon_applied:
            # One log-ratio per bin, plus the arm's stability. The per-bin names
            # are dynamic so the contract does not declare them; the
            # `fsd_only_size_bins` check covers their shape.
            for name, counts in bins.items():
                frame[f"{name}_logR"] = np.log2((counts + 1) / (rng.random(n) * 10 + 1))
            frame["pon_stability"] = rng.random(n) * 5 + offset
        return pd.DataFrame(frame)
    if suffix.startswith(".FSR"):
        short = rng.random(n) * 100 + offset
        long_ = rng.random(n) * 50 + offset
        total = short + long_ + 10
        return pd.DataFrame(
            {
                "region": [f"{c}:0-100000" for c in CHROMS[:n]],
                "total_count": total,
                "short_long_ratio": short / long_,
                "short_long_log2": np.log2(short / long_),
                "short_frac": short / total,
                "long_frac": long_ / total,
            }
        )
    if suffix.startswith(".FSC.gene") or suffix.startswith(".FSC.regions"):
        # Six channels, matching the genome-bin table. The gene path emitted
        # five until the bands were aligned in 0.9.0; five ratios summing to 1
        # is exactly what a well-formed cohort must no longer look like.
        names = [
            "ultra_short",
            "core_short",
            "mono_nucl",
            "di_nucl",
            "long",
            "ultra_long",
        ]
        ratios = np.array(
            [_normalised(len(names), sample_idx * 7 + i) for i in range(n)]
        )
        frame = {"gene": [f"GENE{i}" for i in range(n)], "total": rng.random(n) * 500}
        if suffix.startswith(".FSC.regions"):
            # Annotations copied from the gene BED asset. Constant across
            # samples by construction, which is why the contract declares them
            # vary=NEVER; a region may carry both flags or neither.
            frame["strand"] = ["+" if i % 2 == 0 else "-" for i in range(n)]
            # int 0/1, as `gene_bed.py` produces them: it parses
            # `fields[8] == "1"` into a Python bool, which lands as int64.
            # This fixture wrote strings and agreed with a contract that also
            # said string -- the two matched each other and neither matched the
            # code, which is how the mismatch survived.
            frame["is_e1"] = [1 if i % 3 == 0 else 0 for i in range(n)]
            frame["is_alt_e1"] = [1 if i % 4 == 0 else 0 for i in range(n)]
        for i, name in enumerate(names):
            frame[f"{name}_ratio"] = ratios[:, i]
        frame["normalized_depth"] = rng.random(n) * 3
        return pd.DataFrame(frame)
    if suffix.startswith(".FSC"):
        channels = [
            "ultra_short",
            "core_short",
            "mono_nucl",
            "di_nucl",
            "long",
            "ultra_long",
        ]
        data = {c: np.round(rng.random(n) * 100 + offset, 2) for c in channels}
        frame = {
            "chrom": CHROMS[:n],
            "start": np.arange(n) * 1000,
            "end": np.arange(n) * 1000 + 1000,
        }
        frame.update(data)
        frame["total"] = np.round(np.sum(list(data.values()), axis=0), 2)
        for c in channels:
            frame[f"{c}_log2"] = np.log2(data[c] + 1)
        return pd.DataFrame(frame)
    if suffix.startswith(".WPS_background"):
        return pd.DataFrame(
            {
                "group_id": [f"AluY_{i}" for i in range(n)],
                "nrl_bp": 180.0 + rng.random(n) * 20 + offset,
                "nrl_deviation_bp": rng.random(n) * 10 + offset,
                "periodicity_score": rng.random(n) + offset,
                "adjusted_score": rng.random(n) + offset,
                "fragment_ratio": rng.random(n) * 0.01 + offset,
            }
        )
    if suffix.startswith(".WPS"):
        vectors = [rng.random(8) + offset for _ in range(n)]
        frame = {
            "region_type": (["TSS", "CTCF"] * n)[:n],
            "wps_nuc": vectors,
            "wps_tf": [v * 0.5 for v in vectors],
            "prot_frac_nuc": [v * 0.3 for v in vectors],
            "prot_frac_tf": [v * 0.2 for v in vectors],
        }
        if "panel" in suffix:
            frame["local_depth"] = rng.random(n) * 100
        if pon_applied:
            # What PON scoring produces. Written here so the default cohort is a
            # *scored* one -- the realistic case, and the one where the gate
            # requires these columns.
            frame["wps_nuc_z"] = [rng.standard_normal(8) + offset for _ in range(n)]
            frame["wps_log_amplitude"] = rng.random(n) * 2 + offset
            frame["wps_log_amplitude_z"] = rng.standard_normal(n) + offset
            frame["wps_shape_corr"] = rng.random(n) * 0.5 + 0.4 + offset
            frame["wps_shape_corr_z"] = rng.standard_normal(n) + offset
            frame["wps_phase_shift_bp"] = rng.integers(-5, 6, n).astype(float)
            # All-False: no anchor's search hit its own edge, the healthy case
            # the contract declares as a legitimate constant.
            frame["wps_phase_at_search_limit"] = [False] * n
        return pd.DataFrame(frame)
    if suffix.startswith(".OCF") and "sync" in suffix:
        # Per-position orientation profile. `position` is the offset grid and
        # is identical in every sample by construction, so it must not be
        # randomised or the vary=NEVER declaration fails for the wrong reason.
        tissues = [f"tissue_{i}" for i in range(n)]
        positions = list(range(-n, n))
        rows = [(t, p_) for t in tissues for p_ in positions]
        return pd.DataFrame(
            {
                "tissue": [t for t, _ in rows],
                "position": [p_ for _, p_ in rows],
                "left_count": rng.random(len(rows)) * 100 + offset,
                "left_norm": rng.random(len(rows)) + offset,
                "right_count": rng.random(len(rows)) * 100 + offset,
                "right_norm": rng.random(len(rows)) + offset,
            }
        )
    if suffix.startswith(".OCF"):
        frame = {
            "tissue": [f"tissue_{i}" for i in range(n)],
            "OCF": rng.random(n) * 20 - 10 + offset,
        }
        if pon_applied:
            frame["ocf_z"] = rng.standard_normal(n) + offset
        return pd.DataFrame(frame)
    if suffix.startswith(".metadata"):
        row = {
            "sample_id": f"S{sample_idx}",
            "total_fragments": 1_000_000 + sample_idx * 137,
            "genome": "hg19",
            # Provenance, as `run_features` stamps it. Without these the gate
            # cannot tell which build produced the directory, and cannot decide
            # whether the PON-derived columns ought to be present.
            "krewlyzer_version": "0.9.0",
            "pon_applied": pon_applied,
            "pon_model": "xs2.duplex.pon.parquet" if pon_applied else "",
            "pon_cohort_digest": "0123456789abcdef" if pon_applied else "",
            "pon_krewlyzer_version": "0.9.0" if pon_applied else "",
        }
        return pd.DataFrame([row])

    # Generic: satisfy the declared columns and nothing more.
    frame: Dict[str, object] = {}
    for col in rule.columns:
        if col.kind is Kind.STRING:
            frame[col.name] = [f"{col.name}_{i}" for i in range(n)]
        elif col.kind is Kind.LIST:
            frame[col.name] = [rng.random(8) + offset for _ in range(n)]
        else:
            frame[col.name] = rng.random(n) * 10 + offset
    return pd.DataFrame(frame)


def write_sample(
    root: Path, sample: str, sample_idx: int, pon_applied: bool = True
) -> Path:
    sample_dir = root / sample
    sample_dir.mkdir(parents=True, exist_ok=True)
    for rule in CONTRACT:
        build_frame(rule, sample_idx, pon_applied=pon_applied).to_parquet(
            sample_dir / f"{sample}{rule.suffix}", index=False
        )
    return sample_dir


def write_cohort(root: Path, n_samples: int = 3, pon_applied: bool = True) -> List[str]:
    """A well-formed cohort. Scored by default, which is the realistic case.

    Pass ``pon_applied=False`` for a `--skip-pon` cohort: the PON-derived
    columns are then legitimately absent, and the gate must not ask for them.
    """
    names = [f"S{i}" for i in range(n_samples)]
    for idx, name in enumerate(names):
        write_sample(root, name, idx, pon_applied=pon_applied)
    return names
