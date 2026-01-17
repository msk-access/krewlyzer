"""
Output format handling for krewlyzer.

Provides unified output format control with smart defaults:
- TSV: Human-readable tables (FSD, FSR, FSC, Motif, OCF, UXM)
- Parquet: Vector data (WPS), models (PON)
- JSON: ML pipelines, APIs (unified sample JSON)

Usage:
    writer = OutputWriter(OutputFormat.AUTO)
    writer.write(df, Path("sample.FSD"), tool="fsd")  # -> sample.FSD.tsv
"""

from enum import Enum
from pathlib import Path
from typing import Optional
import pandas as pd
import logging

logger = logging.getLogger(__name__)


class OutputFormat(Enum):
    """Supported output formats."""
    TSV = "tsv"
    PARQUET = "parquet"
    JSON = "json"
    AUTO = "auto"


# Smart defaults per tool
TOOL_DEFAULTS = {
    "fsd": OutputFormat.TSV,
    "fsr": OutputFormat.TSV,
    "fsc": OutputFormat.TSV,
    "wps": OutputFormat.PARQUET,  # Vector data - keep Parquet
    "motif": OutputFormat.TSV,
    "ocf": OutputFormat.TSV,
    "uxm": OutputFormat.TSV,
    "mfsd": OutputFormat.TSV,
    "gc_factors": OutputFormat.TSV,  # Was CSV, now TSV
}


class OutputWriter:
    """
    Unified output writer supporting multiple formats.
    
    Attributes:
        format: Default output format (AUTO uses smart per-tool defaults)
    """
    
    def __init__(self, format: OutputFormat = OutputFormat.AUTO):
        self.format = format
    
    def write(
        self,
        data: pd.DataFrame,
        path: Path,
        tool: Optional[str] = None,
        index: bool = False
    ) -> Path:
        """
        Write DataFrame in specified format.
        
        Args:
            data: DataFrame to write
            path: Output path (suffix will be replaced based on format)
            tool: Tool name for smart default resolution
            index: Whether to include index in output
        
        Returns:
            Actual path written (with correct suffix)
        """
        fmt = self._resolve_format(tool)
        
        if fmt == OutputFormat.TSV:
            out_path = path.with_suffix(".tsv")
            data.to_csv(out_path, sep="\t", index=index)
            logger.debug(f"Wrote TSV: {out_path}")
            
        elif fmt == OutputFormat.PARQUET:
            out_path = path.with_suffix(".parquet")
            data.to_parquet(out_path, index=index)
            logger.debug(f"Wrote Parquet: {out_path}")
            
        elif fmt == OutputFormat.JSON:
            out_path = path.with_suffix(".json")
            data.to_json(out_path, orient="records", indent=2)
            logger.debug(f"Wrote JSON: {out_path}")
        
        else:
            raise ValueError(f"Unknown format: {fmt}")
        
        return out_path
    
    def _resolve_format(self, tool: Optional[str]) -> OutputFormat:
        """Resolve AUTO format to concrete format."""
        if self.format == OutputFormat.AUTO:
            return TOOL_DEFAULTS.get(tool, OutputFormat.TSV)
        return self.format
    
    @classmethod
    def from_string(cls, format_str: str) -> "OutputWriter":
        """Create OutputWriter from string format name."""
        try:
            fmt = OutputFormat(format_str.lower())
        except ValueError:
            logger.warning(f"Unknown format '{format_str}', using AUTO")
            fmt = OutputFormat.AUTO
        return cls(fmt)


def parse_output_format(format_str: Optional[str]) -> OutputFormat:
    """Parse format string to OutputFormat enum."""
    if format_str is None:
        return OutputFormat.AUTO
    try:
        return OutputFormat(format_str.lower())
    except ValueError:
        logger.warning(f"Unknown format '{format_str}', using AUTO")
        return OutputFormat.AUTO
