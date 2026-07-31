"""Every contracted table needs an interpretation, and it must not drift.

`meaning.py` is a second structure keyed by output suffix, alongside
`contract.py`. Two structures keyed by the same thing is exactly how the FSC
bands, the version string and the MDS scale all drifted in this repository, so
it only survives with a test that makes divergence impossible.
"""

from __future__ import annotations

import pytest

from krewlyzer.validate.contract import CONTRACT
from krewlyzer.validate.meaning import MEANINGS

_SUFFIXES = {rule.suffix for rule in CONTRACT}


def test_every_contracted_table_has_a_meaning():
    """A new output must not ship without a reader being told what it is."""
    missing = sorted(_SUFFIXES - set(MEANINGS))
    assert not missing, (
        f"contracted tables with no entry in meaning.py: {missing}. "
        "A table nobody can interpret is one nobody will use correctly."
    )


def test_no_meaning_describes_a_table_that_does_not_exist():
    """The other direction: a renamed suffix leaves an orphan behind."""
    orphans = sorted(set(MEANINGS) - _SUFFIXES)
    assert (
        not orphans
    ), f"meaning.py describes tables absent from the contract: {orphans}"


@pytest.mark.parametrize("suffix", sorted(MEANINGS))
def test_each_meaning_is_substantive(suffix):
    """Guards against a placeholder added only to satisfy the test above.

    An entry may be short *if* it is an explicit cross-reference — the
    `.ontarget` variants are the same measurement restricted to captured
    regions, and "As FSC, over captured regions." says that better than a
    padded paragraph would. Anything else has to stand on its own.
    """
    m = MEANINGS[suffix]
    text = m.measures.strip()
    assert text.endswith("."), f"{suffix}: `measures` is not a sentence"
    if not text.startswith("As "):
        assert len(text) > 40, (
            f"{suffix}: `measures` is too short to be useful and is not a "
            "cross-reference"
        )


def test_directions_are_stated_wherever_they_apply():
    """The notebook's central point: a wrong direction is the commonest misread.

    Only `metadata` is legitimately directionless — it is provenance, not a
    measurement. Anything else lacking a direction is an omission.
    """
    silent = {s for s, m in MEANINGS.items() if not m.cancer_direction}
    assert silent == {".metadata.parquet"}, (
        f"tables with no stated cancer direction: {sorted(silent)}. "
        "Every measurement has one, and leaving it unsaid is how MDS came to be "
        "documented backwards."
    )


def test_mds_direction_is_recorded_as_lower():
    """The specific error this registry exists because of.

    `region-mds.md` asserted the opposite of its own threshold table for a year.
    Pinned explicitly rather than left to prose review.
    """
    for suffix in (".MDS.parquet", ".MDS.ontarget.parquet", ".MDS.exon.parquet"):
        direction = MEANINGS[suffix].cancer_direction or ""
        assert (
            "lower" in direction.lower()
        ), f"{suffix} does not record lower MDS as the tumour direction"


def test_no_thresholds_are_recorded():
    """Directions are robust; magnitudes are cohort-specific.

    Every numeric band examined in the source notebook was a display default or
    refuted outright — the documented ATAC/TFBS entropy range flags a healthy
    N(167,30) distribution as abnormal. A number here would acquire an authority
    it has not earned, so the registry carries none.
    """
    import re

    # A bare decimal or a comparison against one. Sizes in bp and citation years
    # are descriptive, not thresholds, so they are allowed.
    threshold = re.compile(r"[<>]=?\s*\d|z\s*[<>]|\b\d+\.\d+\b")
    offenders = []
    for suffix, m in MEANINGS.items():
        for field in (m.measures, m.cancer_direction or "", m.caveat or ""):
            text = re.sub(
                r"\d+\s*bp|\b(19|20)\d{2}\b|log2\(256\)|N\(\d+,\s*\d+\)", "", field
            )
            if threshold.search(text):
                offenders.append((suffix, field))
    assert not offenders, f"threshold-like values in meaning.py: {offenders}"
