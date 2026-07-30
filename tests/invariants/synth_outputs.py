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


def build_frame(rule: TableRule, sample_idx: int) -> pd.DataFrame:
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
        ratios = np.array([_normalised(5, sample_idx * 7 + i) for i in range(n)])
        frame = {"gene": [f"GENE{i}" for i in range(n)], "total": rng.random(n) * 500}
        for i, name in enumerate(
            ["ultra_short", "core_short", "mono_nucl", "di_nucl", "long"]
        ):
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
        return pd.DataFrame(frame)
    if suffix.startswith(".metadata"):
        return pd.DataFrame(
            [
                {
                    "sample_id": f"S{sample_idx}",
                    "total_fragments": 1_000_000 + sample_idx * 137,
                    "genome": "hg19",
                }
            ]
        )

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


def write_sample(root: Path, sample: str, sample_idx: int) -> Path:
    sample_dir = root / sample
    sample_dir.mkdir(parents=True, exist_ok=True)
    for rule in CONTRACT:
        build_frame(rule, sample_idx).to_parquet(
            sample_dir / f"{sample}{rule.suffix}", index=False
        )
    return sample_dir


def write_cohort(root: Path, n_samples: int = 3) -> List[str]:
    names = [f"S{i}" for i in range(n_samples)]
    for idx, name in enumerate(names):
        write_sample(root, name, idx)
    return names
