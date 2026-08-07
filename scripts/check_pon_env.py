#!/usr/bin/env python
"""Does this install actually have the PON fixes? Ask the code, not the version.

`krewlyzer --version` reports the previous release until the bump, and a CLI
flag only proves the *Python* side is current. Neither can see a stale compiled
extension -- `_core*.so` is gitignored and rebuilt by the build backend, so a
`git pull` alone leaves the old binary in place and every Rust fix silently
absent.

That is not hypothetical. A cohort was re-aggregated after two Rust-adjacent
PRs merged and came back byte-identical, and telling "the build is stale" from
"the fix does not touch this path" took a column-by-column diff against the
previous model. This script answers it in a second.

Each probe drives the real backend on input whose correct answer is known, and
fails loudly with the remedy. Run it before any PON build:

    python scripts/check_pon_env.py && sbatch scripts/build_pon_array.sh ...
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

REMEDY = "    pip install -e '.[all]'   (or: maturin develop --release)"

WPS_SCHEMA = pa.schema(
    [
        ("region_id", pa.string()),
        ("wps_nuc", pa.list_(pa.float32())),
        ("wps_tf", pa.list_(pa.float32())),
    ]
)


def _wps_cohort(directory: Path, n: int, value: float) -> list[str]:
    """``n`` donors whose profiles are identical, so every sigma must be NaN."""
    paths = []
    for i in range(n):
        path = directory / f"s{i}.WPS.parquet"
        pq.write_table(
            pa.Table.from_pandas(
                pd.DataFrame(
                    {
                        "region_id": ["TSS|PROBE|1"],
                        "wps_nuc": [[value] * 200],
                        "wps_tf": [[value] * 200],
                    }
                ),
                schema=WPS_SCHEMA,
            ),
            path,
        )
        paths.append(str(path))
    return paths


def probe_vector_sigma() -> tuple[bool, str]:
    """`element_wise_std`: identical donor vectors must give NaN, not residue.

    21 copies of 0.1 -- the xs2 cohort size -- used to yield 7.6e-9, which
    passes every `std > 0` guard and becomes a divisor. A one-unit deviation
    against it scores z ~ 1e11.
    """
    from krewlyzer import _core

    with tempfile.TemporaryDirectory() as directory:
        paths = _wps_cohort(Path(directory), 21, 0.1)
        result = _core.pon_builder.compute_wps_baseline(paths)

    sigma = np.asarray(result["TSS|PROBE|1"]["wps_nuc_std"], dtype=float)
    finite = sigma[np.isfinite(sigma)]
    if finite.size:
        return (
            False,
            f"wps_nuc_std has {finite.size} finite values, max {finite.max():.3e}",
        )
    return True, "identical donor vectors give NaN"


def probe_entropy_zscore() -> tuple[bool, str]:
    """`region_entropy`: an unmeasurable sigma must give NaN, not 0.0.

    A z-score of zero says "exactly at the healthy baseline" -- the most
    confident statement the column can make, from a comparison that could not
    be made.
    """
    from krewlyzer import _core

    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        raw = root / "s.TFBS.tsv"
        pd.DataFrame({"label": ["PROBE"], "entropy": [0.85]}).to_csv(
            raw, sep="\t", index=False
        )
        pon = root / "probe.pon.parquet"
        pd.concat(
            [
                pd.DataFrame([{"table": "metadata", "schema_version": "1.0"}]),
                pd.DataFrame(
                    {
                        "table": "tfbs_baseline",
                        "label": ["PROBE"],
                        "entropy_mean": [0.71],
                        "entropy_std": [0.0],
                    }
                ),
            ],
            ignore_index=True,
        ).to_parquet(pon)

        out = root / "s.TFBS.z.tsv"
        _core.region_entropy.apply_pon_zscore(
            str(raw), str(pon), str(out), "tfbs_baseline"
        )
        if not out.exists():
            return False, "apply_pon_zscore wrote nothing"
        z = pd.read_csv(out, sep="\t")["z_score"].iloc[0]

    if not np.isnan(z):
        return False, f"unmeasurable sigma scored z = {z}, expected NaN"
    return True, "unmeasurable sigma gives no z-score"


def probe_fsd_accepts_any_format() -> tuple[bool, str]:
    """`_compute_fsd_baseline` must read parquet, not silently skip it.

    The Rust reader takes plain TSV; a run-all directory has none. Handing it
    parquet used to yield zero arms.
    """
    from krewlyzer.pon.build import _compute_fsd_baseline

    bins = [f"{s}-{s + 4}" for s in range(65, 105, 5)]
    with tempfile.TemporaryDirectory() as directory:
        paths = []
        for i in range(3):
            frame = pd.DataFrame(
                [
                    {
                        "arm": arm,
                        **{b: float((j + 1) * (i + 1)) for j, b in enumerate(bins)},
                    }
                    for arm in ("1p", "1q")
                ]
            )
            path = Path(directory) / f"s{i}.FSD.parquet"
            frame.to_parquet(path)
            paths.append(str(path))
        baseline = _compute_fsd_baseline([], paths)

    if not baseline or len(baseline.arms) != 2:
        return False, "parquet FSD tables produced no baseline"
    return True, "FSD reads parquet as well as TSV"


PROBES = (
    ("vector sigma (element_wise_std)", probe_vector_sigma),
    ("entropy z-score (region_entropy.rs)", probe_entropy_zscore),
    ("FSD format normalisation", probe_fsd_accepts_any_format),
)


def main() -> int:
    try:
        import krewlyzer

        print(
            f"krewlyzer {getattr(krewlyzer, '__version__', '?')} "
            f"from {Path(krewlyzer.__file__).parent}"
        )
    except Exception as exc:  # pragma: no cover - import failure is the answer
        print(f"FATAL: cannot import krewlyzer: {exc}")
        return 3

    failed = []
    for name, probe in PROBES:
        try:
            ok, detail = probe()
        except Exception as exc:
            ok, detail = False, f"{type(exc).__name__}: {exc}"
        print(f"  [{'ok  ' if ok else 'FAIL'}] {name}: {detail}")
        if not ok:
            failed.append(name)

    if failed:
        print(f"\n{len(failed)} probe(s) failed. This install predates the PON")
        print("fixes, most likely because the compiled extension was not rebuilt.")
        print("A `git pull` does not rebuild it:\n")
        print(REMEDY)
        return 1

    print("\nAll probes pass -- safe to build a PON with this install.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
