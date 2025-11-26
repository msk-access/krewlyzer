import pytest
from pathlib import Path
from krewlyzer.wrapper import run_all
from unittest.mock import patch, MagicMock

@patch("krewlyzer.wrapper.motif")
@patch("krewlyzer.wrapper.fsc")
@patch("krewlyzer.wrapper.fsr")
@patch("krewlyzer.wrapper.fsd")
@patch("krewlyzer.wrapper.wps")
@patch("krewlyzer.wrapper.ocf")
@patch("krewlyzer.wrapper.uxm")
@patch("krewlyzer.wrapper.mfsd")
def test_run_all_logic(mock_mfsd, mock_uxm, mock_ocf, mock_wps, mock_fsd, mock_fsr, mock_fsc, mock_motif, tmp_path):
    bam_file = tmp_path / "test.bam"
    bam_file.touch()
    reference = tmp_path / "ref.fa"
    reference.touch()
    output = tmp_path / "output"
    variant_input = tmp_path / "variants.vcf"
    variant_input.touch()
    
    # Run with variant input
    run_all(
        bam_file=bam_file,
        reference=reference,
        output=output,
        variant_input=variant_input,
        threads=2,
        pe_type="PE"
    )
    
    # Check all called
    mock_motif.assert_called_once()
    mock_fsc.assert_called_once()
    mock_fsr.assert_called_once()
    mock_fsd.assert_called_once()
    mock_wps.assert_called_once()
    mock_ocf.assert_called_once()
    mock_uxm.assert_called_once()
    mock_mfsd.assert_called_once()
    
    # Check mfsd args
    args, kwargs = mock_mfsd.call_args
    assert kwargs['bam_path'] == bam_file
    assert kwargs['input_file'] == variant_input
    
    # Run without variant input
    mock_mfsd.reset_mock()
    run_all(
        bam_file=bam_file,
        reference=reference,
        output=output,
        variant_input=None,
        threads=2
    )
    mock_mfsd.assert_not_called()
