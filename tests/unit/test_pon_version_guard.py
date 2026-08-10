"""0.9.0 changed what the features mean, so an older PON measures something else.

Recording `krewlyzer_version` was half a guard. `stamp-pon` writes it and
`validate-pon` flags it missing, but nothing stopped a run scoring against a
model built before the meaning changed -- a fabricated `wps_background`, six
floored sigmas, a region-MDS fitted over 65-400bp while samples are measured
over 65-1000bp. Every z-score against such a model is wrong in a way no schema
check can see, because the columns are all present and finite.

Refusing, not warning: a warning in a log nobody reads is not a guard, and a
plausible wrong number is the exact failure this release exists to remove.
"""

from __future__ import annotations

import pandas as pd
import pytest

from krewlyzer.pon.provenance import (
    ALLOW_OLD_PON_ENV,
    MIN_PON_VERSION,
    check_pon_version,
    format_version,
    parse_version,
)

# ---------------------------------------------------------------------------
# the comparison
# ---------------------------------------------------------------------------


def test_versions_compare_as_numbers_not_strings():
    """`"0.10.0" < "0.9.0"` is true as strings and false as releases.

    That ordering error would silently reverse the guard on the first
    double-digit minor -- accepting exactly the models it exists to refuse.
    """
    assert parse_version("0.10.0") > parse_version("0.9.0")
    assert "0.10.0" < "0.9.0", "the string comparison this avoids"


@pytest.mark.parametrize(
    "text,expected",
    [
        ("0.9.0", (0, 9, 0)),
        ("1.2.3", (1, 2, 3)),
        ("0.9.0+local", (0, 9, 0)),
        ("0.9.0-rc1", (0, 9, 0)),
        ("", None),
        ("   ", None),
        ("unknown", None),
        ("v0.9.0", None),
    ],
)
def test_version_parsing(text, expected):
    assert parse_version(text) == expected


def test_the_floor_is_a_tuple_not_a_version_string():
    """It is a compatibility floor, not the package version.

    It stays at 0.9.0 when krewlyzer is 1.5.0, because what changed is what the
    features mean and that happened once. A string would also trip
    `test_no_module_restates_the_version`, correctly.
    """
    assert isinstance(MIN_PON_VERSION, tuple)
    assert format_version(MIN_PON_VERSION) == "0.9.0"


# ---------------------------------------------------------------------------
# the decision
# ---------------------------------------------------------------------------


def test_a_current_pon_is_accepted():
    assert check_pon_version("0.9.0") is None
    assert check_pon_version("1.4.2") is None


def test_an_older_pon_is_refused_and_told_why():
    complaint = check_pon_version("0.8.3", "xs1.duplex.pon.parquet")
    assert complaint
    assert "0.8.3" in complaint and "0.9.0" in complaint
    assert "xs1.duplex.pon.parquet" in complaint, "does not name the file"
    assert "build-pon" in complaint, "does not say what to do"


def test_a_pon_with_no_version_is_refused():
    """Absent is not old, but it is equally unusable.

    Every PON built before provenance existed carries the defects; without a
    version there is no way to tell one from a good model.
    """
    for recorded in ("", "   ", "unknown"):
        assert check_pon_version(recorded), f"{recorded!r} was accepted"


def test_the_escape_hatch_works_and_must_be_deliberate(monkeypatch):
    """Without a documented way out, someone edits the parquet instead.

    Then the version stops meaning anything at all, which is worse than an old
    model knowingly used.
    """
    assert check_pon_version("0.8.3")
    monkeypatch.setenv(ALLOW_OLD_PON_ENV, "1")
    assert check_pon_version("0.8.3") is None


def test_an_empty_escape_hatch_does_not_count(monkeypatch):
    monkeypatch.setenv(ALLOW_OLD_PON_ENV, "")
    assert check_pon_version("0.8.3")


# ---------------------------------------------------------------------------
# the loader, and everything that reaches a PON
# ---------------------------------------------------------------------------


def _pon(tmp_path, version):
    frame = pd.DataFrame(
        [
            {
                "table": "metadata",
                "schema_version": "1.0",
                "assay": "xs1",
                "n_samples": 21.0,
                "krewlyzer_version": version,
            }
        ]
    )
    path = tmp_path / "t.pon.parquet"
    frame.to_parquet(path)
    return path


def test_the_loader_refuses_an_old_model(tmp_path, caplog):
    from krewlyzer.core.pon_integration import load_pon_model

    assert load_pon_model(_pon(tmp_path, "0.8.3")) is None
    assert any("0.8.3" in r.message for r in caplog.records)


def test_the_loader_accepts_a_current_model(tmp_path):
    from krewlyzer.core.pon_integration import load_pon_model

    assert load_pon_model(_pon(tmp_path, "0.9.0")) is not None


def test_every_subcommand_goes_through_the_guarded_loader():
    """`PonModel.load` bypasses the guard, so nothing outside it may call it.

    `motif` and `region-mds` both did, which meant those subcommands would
    score against a model `run-all` refuses -- invariant #6, single-tool and
    run-all behaving differently.
    """
    import subprocess

    hits = subprocess.run(
        ["grep", "-rn", "PonModel.load", "src/krewlyzer"],
        capture_output=True,
        text=True,
    ).stdout.splitlines()
    offenders = []
    for hit in hits:
        parts = hit.split(":", 2)
        if len(parts) < 3:
            continue  # grep printed something without a body
        path, _lineno, code = parts
        if path.endswith("core/pon_integration.py"):
            continue  # the guard itself
        if code.strip().startswith("#"):
            continue  # a comment naming it, not a call
        offenders.append(hit)
    assert not offenders, "these bypass the version guard:\n" + "\n".join(offenders)


def test_a_refused_pon_does_not_reach_the_rust_scorers():
    """The Rust scorers take a parquet *path* and read it themselves.

    Setting that path regardless of whether the Python load succeeded meant a
    refused model was still normalised against by FSD, WPS, OCF and region
    entropy: the Python half stopped and the Rust half carried on.
    """
    from krewlyzer.core import unified_processor

    source = open((unified_processor.__file__ or "").replace(".pyc", ".py")).read()
    assert "pon_parquet = resolved_pon_model if pon else None" in source
    assert "pon_parquet = resolved_pon_model\n" not in source
