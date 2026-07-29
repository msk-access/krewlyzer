"""Findings produced by the contract gate."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, List, Optional


class Severity(str, Enum):
    ERROR = "error"  # contract violation -> exit 1
    WARN = "warn"  # worth seeing, does not gate
    SKIP = "skip"  # could not be evaluated (e.g. too few samples)


class Category(str, Enum):
    COMPLETION = "completion"
    MISSING = "missing"
    SCHEMA = "schema"
    DOMAIN = "domain"
    DEGENERACY = "degeneracy"
    STRUCTURAL = "structural"  # -> exit 2


@dataclass
class Finding:
    id: str
    severity: Severity
    category: Category
    message: str
    table: Optional[str] = None
    column: Optional[str] = None
    samples: List[str] = field(default_factory=list)
    evidence: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "id": self.id,
            "severity": self.severity.value,
            "category": self.category.value,
            "message": self.message,
            "table": self.table,
            "column": self.column,
            "samples_affected": self.samples,
            "evidence": self.evidence,
        }
