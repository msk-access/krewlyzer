"""Output-contract validation.

``krewlyzer validate-output`` answers a question the pipeline's own tests
cannot: is a *finished* output directory something a downstream consumer can
actually rely on? See :mod:`krewlyzer.validate.contract` for what is asserted
and why.
"""

from .contract import CONTRACT, COMPLETION_MARKER
from .findings import Category, Finding, Severity
from .gate import (
    EXIT_PASS,
    EXIT_STRUCTURAL,
    EXIT_VIOLATION,
    Fingerprint,
    Result,
    check_sample,
    evaluate_cohort,
    run,
)

__all__ = [
    "CONTRACT",
    "COMPLETION_MARKER",
    "Category",
    "EXIT_PASS",
    "EXIT_STRUCTURAL",
    "EXIT_VIOLATION",
    "Finding",
    "Fingerprint",
    "Result",
    "Severity",
    "check_sample",
    "evaluate_cohort",
    "run",
]
