"""
End-to-end tests for run-all with PON model.

Tests full pipeline integration with PON normalization.
"""

import pytest
import pandas as pd
from pathlib import Path

# Skip entire module if test data not available
pytestmark = [
    pytest.mark.skipif(
        not Path("/Users/shahr2/Documents/Github/krewlyzer/tests/data").exists(),
        reason="Test data not available",
    ),
    pytest.mark.slow,
]


class TestOcfProcessor:
    """Unit tests for OCF processor."""

    def test_process_ocf_with_pon(self, tmp_path):
        """Test OCF z-score computation using Rust implementation."""
        from krewlyzer.core.ocf_processor import process_ocf_with_pon
        import pyarrow as pa
        import pyarrow.parquet as pq

        # Create test OCF file
        ocf_data = pd.DataFrame(
            {"tissue": ["region1", "region2", "region3"], "OCF": [0.6, 0.25, 0.5]}
        )
        ocf_path = tmp_path / "test.OCF.tsv"
        ocf_data.to_csv(ocf_path, sep="\t", index=False)

        # Create PON parquet with ocf_baseline table
        baseline_df = pd.DataFrame(
            {
                "table": ["ocf_baseline", "ocf_baseline"],
                "region_id": ["region1", "region2"],
                "ocf_mean": [0.5, 0.3],
                "ocf_std": [0.1, 0.05],
            }
        )
        pon_path = tmp_path / "test_pon.parquet"
        table = pa.Table.from_pandas(baseline_df)
        pq.write_table(table, pon_path)

        # Process
        output_path = tmp_path / "test.OCF.zscore.tsv"
        n_matched = process_ocf_with_pon(ocf_path, pon_path, output_path)

        # Verify
        assert n_matched == 2  # 2 regions matched
        result = pd.read_csv(output_path, sep="\t")
        assert "ocf_z" in result.columns
        # region1: (0.6 - 0.5) / 0.1 = 1.0
        assert pytest.approx(result.iloc[0]["ocf_z"], rel=0.01) == 1.0
        # region2: (0.25 - 0.3) / 0.05 = -1.0
        assert pytest.approx(result.iloc[1]["ocf_z"], rel=0.01) == -1.0


class TestPonModelOnTarget:
    """Tests for on-target baseline fields."""

    def test_ontarget_fields_exist(self):
        """Test PonModel has ontarget fields."""
        from krewlyzer.pon.model import PonModel

        pon = PonModel()
        assert hasattr(pon, "gc_bias_ontarget")
        assert hasattr(pon, "fsd_baseline_ontarget")
        assert pon.gc_bias_ontarget is None
        assert pon.fsd_baseline_ontarget is None
