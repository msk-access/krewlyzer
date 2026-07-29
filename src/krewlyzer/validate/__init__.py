"""Output-contract validation.

``krewlyzer validate-output`` answers a question the pipeline's own tests
cannot: is a *finished* output directory something a downstream consumer can
actually rely on? See :mod:`krewlyzer.validate.contract` for what is asserted
and why.
"""

from .contract import CONTRACT, COMPLETION_MARKER, kreview_suffixes
from .findings import Category, Finding, Severity
from .gate import EXIT_PASS, EXIT_STRUCTURAL, EXIT_VIOLATION, Result, run

__all__ = [
    "CONTRACT",
    "COMPLETION_MARKER",
    "Category",
    "EXIT_PASS",
    "EXIT_STRUCTURAL",
    "EXIT_VIOLATION",
    "Finding",
    "Result",
    "Severity",
    "kreview_suffixes",
    "run",
]
