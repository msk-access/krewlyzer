"""Numeric constants must match what the registry says they are.

Pure regex over the Rust sources plus imports of the Python package, so this
runs without building the extension -- which matters, because a stale build is
exactly when you most want to know that a constant moved.
"""

from __future__ import annotations

import pytest

from krewlyzer.validate.claims import (
    CLAIMS,
    KNOWN_DIVERGENCES,
    ConstantNotFound,
    RustConst,
    RustPattern,
    _RUST_SRC,
)

pytestmark = pytest.mark.unit

_HAVE_RUST_SRC = _RUST_SRC.is_dir()
_needs_rust_src = pytest.mark.skipif(
    not _HAVE_RUST_SRC, reason="Rust sources absent (installed wheel)"
)


def _needs_sources(ref) -> bool:
    return isinstance(ref, (RustConst, RustPattern))


@pytest.mark.parametrize("claim", CLAIMS, ids=lambda c: c.id)
def test_constant_matches_the_registry(claim):
    """The implementation still defines what the registry says it does."""
    if _needs_sources(claim.impl) and not _HAVE_RUST_SRC:
        pytest.skip("Rust sources absent")

    actual = claim.impl.resolve()

    assert actual == claim.value, (
        f"{claim.id}: registry says {claim.value!r}, implementation says "
        f"{actual!r}.\n"
        f"Why this constant matters: {claim.why}\n"
        "If the change is intentional, update claims.py and any documentation "
        "quoting the old value in the same commit."
    )


@pytest.mark.parametrize("div", KNOWN_DIVERGENCES, ids=lambda d: d.id)
def test_known_divergence_is_still_present(div):
    """A recorded disagreement must still be there -- or be removed from the list.

    Asserting that something is still broken looks odd, but it is the only way
    to stop the list rotting. Fixing the underlying issue fails this test, which
    is the prompt to delete the entry; without it a resolved divergence would
    sit in the registry indefinitely, teaching readers to distrust it.
    """
    if (_needs_sources(div.left) or _needs_sources(div.right)) and not _HAVE_RUST_SRC:
        pytest.skip("Rust sources absent")

    left = div.left.resolve()
    right = div.right.resolve()

    if left == div.left_value and right == div.right_value:
        return  # still diverging exactly as recorded

    pytest.fail(
        f"{div.id}: the recorded divergence has changed.\n"
        f"  expected left={div.left_value!r} right={div.right_value!r}\n"
        f"  found    left={left!r} right={right!r}\n"
        f"{div.why}\n"
        "If this is now fixed, delete the entry from KNOWN_DIVERGENCES. If it "
        "moved rather than resolved, update the recorded values."
    )


@_needs_rust_src
def test_a_renamed_constant_fails_rather_than_skipping():
    """Drift must be loud.

    A registry that quietly passes when it cannot find a constant is worse than
    no registry -- renaming is one of the ways a documented value goes stale.
    """
    with pytest.raises(ConstantNotFound):
        RustConst("wps.rs", "NO_SUCH_CONSTANT_EXISTS").resolve()

    with pytest.raises(ConstantNotFound):
        RustConst("no_such_file.rs", "ANYTHING").resolve()


def test_registry_ids_are_unique():
    ids = [c.id for c in CLAIMS] + [d.id for d in KNOWN_DIVERGENCES]
    duplicates = {i for i in ids if ids.count(i) > 1}
    assert not duplicates, f"duplicate registry ids: {sorted(duplicates)}"


def test_every_entry_explains_itself():
    """`why` is the whole value of the registry: a bare number tells the next
    reader nothing about whether they may change it."""
    for entry in (*CLAIMS, *KNOWN_DIVERGENCES):
        assert len(entry.why) > 30, f"{entry.id}: 'why' is too thin to be useful"


def test_the_python_sigma_floor_matches_the_rust_one():
    """One rule, two languages, and nothing else keeps them equal.

    The builder refuses to *write* a residue sigma; the read side refuses to
    *divide* by one, which is what makes the four PONs shipped before 0.9.0
    safe without a rebuild. If the two constants drift, one half of that
    protection silently stops applying.
    """
    import re
    from pathlib import Path

    from krewlyzer.pon.model import SIGMA_FLOOR

    rust = (Path(__file__).resolve().parents[2] / "rust/src/pon_builder.rs").read_text()
    match = re.search(r"pub const SIGMA_FLOOR: f32 = ([0-9.e-]+);", rust)
    assert match, "SIGMA_FLOOR is gone from rust/src/pon_builder.rs"
    assert (
        float(match.group(1)) == SIGMA_FLOOR
    ), f"Rust says {match.group(1)}, Python says {SIGMA_FLOOR}"
