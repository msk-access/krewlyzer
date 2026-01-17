"""
Unit tests for output format classes.

Tests OutputWriter and FeatureSerializer functionality.
"""

import pytest
import pandas as pd
import numpy as np
import json
from pathlib import Path

from krewlyzer.core.output_format import (
    OutputFormat, OutputWriter, TOOL_DEFAULTS, parse_output_format
)
from krewlyzer.core.feature_serializer import FeatureSerializer, NumpyEncoder


class TestOutputFormat:
    """Tests for OutputFormat enum and TOOL_DEFAULTS."""
    
    def test_output_format_values(self):
        """Test OutputFormat enum has expected values."""
        assert OutputFormat.TSV.value == "tsv"
        assert OutputFormat.PARQUET.value == "parquet"
        assert OutputFormat.JSON.value == "json"
        assert OutputFormat.AUTO.value == "auto"
    
    def test_tool_defaults(self):
        """Test default formats for each tool."""
        assert TOOL_DEFAULTS["fsd"] == OutputFormat.TSV
        assert TOOL_DEFAULTS["fsr"] == OutputFormat.TSV
        assert TOOL_DEFAULTS["wps"] == OutputFormat.PARQUET  # Vector data
        assert TOOL_DEFAULTS["motif"] == OutputFormat.TSV
    
    def test_parse_output_format(self):
        """Test format string parsing."""
        assert parse_output_format("tsv") == OutputFormat.TSV
        assert parse_output_format("TSV") == OutputFormat.TSV  # Case insensitive
        assert parse_output_format("parquet") == OutputFormat.PARQUET
        assert parse_output_format("json") == OutputFormat.JSON
        assert parse_output_format("auto") == OutputFormat.AUTO
        assert parse_output_format(None) == OutputFormat.AUTO
        assert parse_output_format("invalid") == OutputFormat.AUTO  # Fallback


class TestOutputWriter:
    """Tests for OutputWriter class."""
    
    def test_write_tsv(self, tmp_path):
        """Test TSV output."""
        writer = OutputWriter(OutputFormat.TSV)
        df = pd.DataFrame({"a": [1, 2], "b": [3, 4]})
        
        path = writer.write(df, tmp_path / "test.FSD", tool="fsd")
        
        assert path.suffix == ".tsv"
        assert path.exists()
        loaded = pd.read_csv(path, sep="\t")
        assert len(loaded) == 2
        assert list(loaded.columns) == ["a", "b"]
    
    def test_write_parquet(self, tmp_path):
        """Test Parquet output."""
        writer = OutputWriter(OutputFormat.PARQUET)
        df = pd.DataFrame({"x": [1.0, 2.0], "y": [3.0, 4.0]})
        
        path = writer.write(df, tmp_path / "test.WPS", tool="wps")
        
        assert path.suffix == ".parquet"
        assert path.exists()
        loaded = pd.read_parquet(path)
        assert len(loaded) == 2
    
    def test_write_json(self, tmp_path):
        """Test JSON output."""
        writer = OutputWriter(OutputFormat.JSON)
        df = pd.DataFrame({"id": ["a", "b"], "value": [1, 2]})
        
        path = writer.write(df, tmp_path / "test", tool="fsd")
        
        assert path.suffix == ".json"
        assert path.exists()
        with open(path) as f:
            data = json.load(f)
        assert len(data) == 2
    
    def test_auto_format_resolution(self, tmp_path):
        """Test AUTO format resolves to tool defaults."""
        writer = OutputWriter(OutputFormat.AUTO)
        df = pd.DataFrame({"a": [1]})
        
        # FSD should get TSV
        path = writer.write(df, tmp_path / "fsd_test", tool="fsd")
        assert path.suffix == ".tsv"
        
        # WPS should get Parquet
        path = writer.write(df, tmp_path / "wps_test", tool="wps")
        assert path.suffix == ".parquet"
    
    def test_from_string(self):
        """Test creating writer from string."""
        writer = OutputWriter.from_string("parquet")
        assert writer.format == OutputFormat.PARQUET
        
        writer = OutputWriter.from_string("invalid")
        assert writer.format == OutputFormat.AUTO


class TestNumpyEncoder:
    """Tests for NumpyEncoder JSON encoder."""
    
    def test_numpy_array(self):
        """Test numpy array encoding."""
        data = {"arr": np.array([1, 2, 3])}
        result = json.dumps(data, cls=NumpyEncoder)
        assert "[1, 2, 3]" in result
    
    def test_numpy_int64(self):
        """Test numpy int64 encoding."""
        data = {"val": np.int64(42)}
        result = json.dumps(data, cls=NumpyEncoder)
        assert "42" in result
    
    def test_numpy_float64(self):
        """Test numpy float64 encoding."""
        data = {"val": np.float64(3.14)}
        result = json.dumps(data, cls=NumpyEncoder)
        assert "3.14" in result


class TestFeatureSerializer:
    """Tests for FeatureSerializer class."""
    
    def test_init(self):
        """Test initialization."""
        serializer = FeatureSerializer("sample_001", "0.3.2")
        assert serializer.sample_id == "sample_001"
        assert serializer.version == "0.3.2"
        assert serializer.features == {}
    
    def test_add_metadata(self):
        """Test adding metadata."""
        serializer = FeatureSerializer("test")
        serializer.add_metadata("total_fragments", 1000000)
        serializer.add_metadata("gc_corrected", True)
        
        assert serializer.metadata["total_fragments"] == 1000000
        assert serializer.metadata["gc_corrected"] == True
    
    def test_add_fsd(self):
        """Test adding FSD data."""
        serializer = FeatureSerializer("test")
        fsd_df = pd.DataFrame({
            "region": ["1p", "1q", "2p"],
            "65-69": [100, 110, 95],
            "70-74": [150, 160, 140],
        })
        
        serializer.add_fsd(fsd_df)
        
        assert "fsd" in serializer.features
        assert serializer.features["fsd"]["arms"] == ["1p", "1q", "2p"]
        assert "65-69" in serializer.features["fsd"]["size_bins"]
        assert len(serializer.features["fsd"]["counts"]) == 3
    
    def test_add_fsr(self):
        """Test adding FSR data."""
        serializer = FeatureSerializer("test")
        fsr_df = pd.DataFrame({
            "region": ["1p", "1q"],
            "short": [0.34, 0.32],
            "long": [0.24, 0.26],
        })
        
        serializer.add_fsr(fsr_df)
        
        assert "fsr" in serializer.features
        assert len(serializer.features["fsr"]) == 2
        assert serializer.features["fsr"][0]["region"] == "1p"
    
    def test_add_motif(self):
        """Test adding motif data."""
        serializer = FeatureSerializer("test")
        edm_df = pd.DataFrame([{"AAAA": 0.02, "AAAC": 0.015}])
        bpm_df = pd.DataFrame([{"AAAA": 0.018}])
        
        serializer.add_motif(edm_df, bpm_df, mds=0.87, mds_z=-0.5)
        
        assert serializer.features["motif"]["mds"] == 0.87
        assert serializer.features["motif"]["mds_z"] == -0.5
        assert "AAAA" in serializer.features["motif"]["edm"]
    
    def test_to_dict(self):
        """Test conversion to dictionary."""
        serializer = FeatureSerializer("sample", "0.3.2")
        serializer.add_metadata("total", 100)
        
        result = serializer.to_dict()
        
        assert result["schema_version"] == "1.0"
        assert result["sample_id"] == "sample"
        assert result["krewlyzer_version"] == "0.3.2"
        assert "timestamp" in result
        assert result["metadata"]["total"] == 100
    
    def test_save_and_load(self, tmp_path):
        """Test saving and loading JSON."""
        serializer = FeatureSerializer("test_sample")
        serializer.add_metadata("fragments", 5000000)
        
        fsd_df = pd.DataFrame({
            "region": ["1p"],
            "65-69": [100],
        })
        serializer.add_fsd(fsd_df)
        
        # Save
        path = serializer.save(tmp_path / "test_sample")
        
        assert path.exists()
        assert path.suffix == ".json"
        assert "features" in path.name
        
        # Load and verify
        with open(path) as f:
            data = json.load(f)
        
        assert data["sample_id"] == "test_sample"
        assert data["metadata"]["fragments"] == 5000000
        assert data["features"]["fsd"]["arms"] == ["1p"]
    
    def test_from_outputs(self, tmp_path):
        """Test creating serializer from existing outputs."""
        sample_id = "test_sample"
        
        # Create mock output files
        fsd_df = pd.DataFrame({"region": ["1p"], "65-69": [10]})
        fsd_df.to_csv(tmp_path / f"{sample_id}.FSD.tsv", sep="\t", index=False)
        
        mds_df = pd.DataFrame({"Sample": [sample_id], "MDS": [0.85]})
        mds_df.to_csv(tmp_path / f"{sample_id}.MDS.tsv", sep="\t", index=False)
        
        # Load
        serializer = FeatureSerializer.from_outputs(sample_id, tmp_path)
        
        assert "fsd" in serializer.features
        assert serializer.features["motif"]["mds"] == 0.85


class TestFeatureSerializerEdgeCases:
    """Edge case tests for FeatureSerializer."""
    
    def test_empty_dataframe(self):
        """Test handling empty DataFrames."""
        serializer = FeatureSerializer("test")
        serializer.add_fsd(pd.DataFrame())
        serializer.add_fsr(None)
        
        assert "fsd" not in serializer.features
        assert "fsr" not in serializer.features
    
    def test_wps_parquet(self, tmp_path):
        """Test WPS parquet loading."""
        sample_id = "wps_test"
        
        wps_df = pd.DataFrame({
            "region_id": ["TSS_1", "TSS_2"],
            "wps_nuc_mean": [0.5, 0.6],
            "wps_tf_mean": [0.3, 0.4],
        })
        wps_df.to_parquet(tmp_path / f"{sample_id}.WPS.parquet")
        
        serializer = FeatureSerializer.from_outputs(sample_id, tmp_path)
        
        assert "wps" in serializer.features
        assert serializer.features["wps"]["regions"] == ["TSS_1", "TSS_2"]
