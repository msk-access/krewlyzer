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
import re

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


def test_the_loader_refuses_an_old_model(tmp_path):
    """Refused means the run stops. It used to return None, and the run finished unscored with
    pon_applied=False -- indistinguishable from a deliberate --skip-pon."""
    from krewlyzer.core.pon_integration import UnusablePonError, load_pon_model

    with pytest.raises(UnusablePonError, match="0.8.3"):
        load_pon_model(_pon(tmp_path, "0.8.3"))


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
    import ast
    import inspect

    from krewlyzer.core import unified_processor

    tree = ast.parse(inspect.getsource(unified_processor))

    # Parsed, not grepped. The invariant is "`pon_parquet` is never bound to a
    # path unless the guarded load returned a model" -- a substring check pins
    # one spelling of that and passes the moment someone writes an equivalent
    # one. Every binding must be either None or conditional on `pon`.
    bindings = [
        node.value
        for node in ast.walk(tree)
        if isinstance(node, ast.Assign)
        for target in node.targets
        if isinstance(target, ast.Name) and target.id == "pon_parquet"
    ]
    assert bindings, "pon_parquet is gone; this test needs rewriting"
    for value in bindings:
        is_none = isinstance(value, ast.Constant) and value.value is None
        is_guarded = isinstance(value, ast.IfExp) and any(
            isinstance(n, ast.Name) and n.id == "pon" for n in ast.walk(value.test)
        )
        assert is_none or is_guarded, (
            "pon_parquet is bound unconditionally: a refused PON would reach "
            f"the Rust scorers again ({ast.dump(value)[:90]})"
        )

    # And the scorers that take the path must be handed *that* variable, not a
    # loaded model. WPS joined this set in 0.9.0 when its scoring moved to
    # Rust, so it is named explicitly rather than assumed.
    source = inspect.getsource(unified_processor)
    for call in ("apply_wps_pon(",):
        assert call in source, f"{call} is gone; the guard's coverage changed"
    calls = re.findall(r"apply_wps_pon\((.*?)\n\s*\)", source, re.S)
    assert calls, "no apply_wps_pon call sites found"
    for c in calls:
        assert "pon_parquet" in c, (
            "apply_wps_pon is not being handed the guarded path; a refused PON "
            f"would still be scored against:\n{c}"
        )


def test_region_entropy_goes_through_the_guarded_loader_too():
    """It handed the PON path straight to Rust, which never checks the version, so a refused
    model still z-scored TFBS and ATAC (review of #88 PR A)."""
    import inspect

    from krewlyzer import region_entropy

    assert "load_requested_pon(" in inspect.getsource(region_entropy)


def test_an_unusable_pon_stops_a_single_tool_before_it_writes_anything(tmp_path):
    """Loading after the pass left unscored tables under their final names, beside whatever
    marker the directory already had."""
    from krewlyzer.core.pon_integration import UnusablePonError
    from krewlyzer.core.sample_processor import extract_sample
    from krewlyzer.core.unified_processor import run_features
    from tests.invariants.synth import make_sample

    sample = make_sample(tmp_path / "s")
    extract_sample(
        sample.bam,
        sample.reference,
        output_dir=tmp_path / "bed",
        sample_name="s",
        threads=1,
    )
    junk = tmp_path / "junk.pon.parquet"
    junk.write_bytes(b"not a parquet file")
    out = tmp_path / "out"
    with pytest.raises(UnusablePonError):
        run_features(
            bed_path=tmp_path / "bed" / "s.bed.gz",
            output_dir=out,
            sample_name="s",
            fsc_bins=sample.assets.bins,
            enable_fsc=True,
            enable_fsd=True,
            fsd_arms=sample.assets.arms,
            gc_correct=False,
            pon_model=junk,
            threads=1,
        )
    assert not out.exists() or not any(out.iterdir()), sorted(
        p.name for p in out.iterdir()
    )


def test_a_single_tool_refuses_to_score_a_sample_against_a_second_pon(tmp_path):
    """The sample's record names the PON it was scored against. A single tool scoring it against
    another one would leave tables from two PONs under one record, with nothing to flag it
    (review of #88 PR A). Refused before any work; the same PON again is fine."""
    import pandas as pd

    from krewlyzer.core.output_provenance import stamp_metadata
    from krewlyzer.core.pon_integration import UnusablePonError, load_requested_pon

    x_dir, y_dir = tmp_path / "x", tmp_path / "y"
    x_dir.mkdir(), y_dir.mkdir()
    pon_x, pon_y = _pon(x_dir, "0.9.0"), _pon(y_dir, "0.9.0")
    pon_y = pon_y.rename(y_dir / "other.pon.parquet")
    out = tmp_path / "out"
    out.mkdir()
    pd.DataFrame([{"sample_id": "s"}]).to_parquet(out / "s.metadata.parquet")
    model_x = load_requested_pon(pon_x)
    stamp_metadata(out, "s", pon=model_x, pon_path=pon_x, output_format="parquet")

    assert load_requested_pon(pon_x, record_dir=out, sample="s") is not None
    with pytest.raises(
        UnusablePonError,
        match="t.pon.parquet.*other.pon.parquet|other.pon.parquet.*t.pon.parquet",
    ):
        load_requested_pon(pon_y, record_dir=out, sample="s")


def test_the_refusal_reads_the_record_this_run_writes_not_a_stale_copy(tmp_path):
    """A stale `.tsv` record from an earlier run (PON X) beside this run's `.parquet` (PON Y):
    checking Y against the `.tsv` refused the run's own PON."""
    import pandas as pd

    from krewlyzer.core.output_provenance import refuse_a_second_pon, stamp_metadata
    from krewlyzer.core.pon_integration import load_requested_pon

    x_dir, y_dir = tmp_path / "x", tmp_path / "y"
    x_dir.mkdir(), y_dir.mkdir()
    pon_x = _pon(x_dir, "0.9.0")
    pon_y = _pon(y_dir, "0.9.0").rename(y_dir / "other.pon.parquet")
    out = tmp_path / "out"
    out.mkdir()
    pd.DataFrame([{"sample_id": "s"}]).to_csv(
        out / "s.metadata.tsv", sep="\t", index=False
    )
    stamp_metadata(
        out, "s", pon=load_requested_pon(pon_x), pon_path=pon_x, output_format="tsv"
    )
    pd.DataFrame([{"sample_id": "s"}]).to_parquet(out / "s.metadata.parquet")
    stamp_metadata(
        out, "s", pon=load_requested_pon(pon_y), pon_path=pon_y, output_format="parquet"
    )

    refuse_a_second_pon(
        out, "s", load_requested_pon(pon_y), pon_y, output_format="parquet"
    )


def test_run_all_owns_the_record_for_region_mds():
    """run-all rewrites the whole record, so its region-MDS call must not run the second-PON
    check (a stale copy could refuse run-all's own PON, and run-all's catch would hide it).
    """
    import inspect

    from krewlyzer import wrapper

    source = inspect.getsource(wrapper)
    assert "owns_sample_record=True" in source.split("_region_mds(", 1)[1].split(
        ")", 1
    )[0] or (
        "owns_sample_record=True"
        in source[source.index("_region_mds(") : source.index("_region_mds(") + 1500]
    )
