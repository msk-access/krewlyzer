
import os
import tempfile
import gzip
import shutil
import pysam
import pytest
from pathlib import Path
from click.testing import CliRunner
import typer

from krewlyzer.wrapper import run_all

def create_dummy_data(temp_dir):
    """
    Creates dummy input files for run-all.
    """
    temp_path = Path(temp_dir)
    
    # 1. Dummy BAM & Ref
    bam = temp_path / "test_sample.bam"
    
    # Create valid BAM with header
    header = { 'HD': {'VN': '1.0'}, 'SQ': [{'LN': 100000, 'SN': 'chr1'}] }
    with pysam.AlignmentFile(str(bam), "wb", header=header) as outf:
        pass
    pysam.index(str(bam)) # Create index (requires valid BAM)
    
    ref = temp_path / "genome.fa"
    with open(ref, "w") as f:
        f.write(">chr1\n" + "N" * 100000 + "\n")
    pysam.faidx(str(ref))
    
    # 2. Dummy Output Dir
    out_dir = temp_path / "output"
    out_dir.mkdir()
    
    # 3. Pre-create BED file to skip 'extract' step
    bed_gz = out_dir / "test_sample.bed.gz"
    
    # Write some valid BED data
    fragments = [
        ("chr1", 1000, 1167, 0.5), # 167bp
        ("chr1", 2000, 2100, 0.4), # 100bp
    ]
    
    with gzip.open(bed_gz, 'wt') as f:
        for chrom, start, end, gc in fragments:
            f.write(f"{chrom}\t{start}\t{end}\t{gc:.4f}\n")
            
    # 4. Dummy Resources (Bins, Arms, etc.)
    # We need valid resource files because the pipeline will read them.
    
    # FSC Bins
    bins = temp_path / "bins.bed"
    with open(bins, "w") as f:
        f.write("chr1\t0\t5000\tBin1\n")
        
    # FSD Arms
    arms = temp_path / "arms.bed"
    with open(arms, "w") as f:
        f.write("chr1\t0\t10000\tArm1\n")

    # OCF Regions
    ocr = temp_path / "ocr.bed"
    with open(ocr, "w") as f:
        f.write("chr1\t100\t200\tTissueA\n")

    # WPS Genes (TranscriptAnno) - TSV
    wps = temp_path / "wps.tsv"
    with open(wps, "w") as f:
        f.write("Gene1\tchr1\t1000\t2000\t+\n")
        
    return bam, ref, out_dir, bins, arms, ocr, wps
    # The wrapper uses default path resolved relative to package.
    # We can't easily inject custom WPS resource via CLI arguments for `run-all` 
    # because `run-all` doesn't expose `tsv_input` for WPS!
    # Wait, `run-all` does NOT expose WPS resource override.
    # It hardcodes defaults or uses what's in `wps.py`.
    # My updated wrapper uses `pkg_dir / ...`.
    # If I run this test, it will look for data in `krewlyzer/data/...`.
    # Those files SHOULD exist in the repo.
    # If they are missing in the environment, the test will fail.
    # Validating existence of default resources in the environment might be tricky if not installed.
    # However, since I am running in the repo, `krewlyzer/data` should be accessible relative to `krewlyzer/wrapper.py`.
    
    return bam, ref, out_dir, bins, arms, ocr, wps

def test_run_all_unified():
    with tempfile.TemporaryDirectory() as temp_dir:
        bam, ref, out_dir, bins, arms, ocr, wps = create_dummy_data(temp_dir)
        
        print(f"Testing in {temp_dir}")
        
        # Call run_all directly
        try:
            run_all(
                bam_input=bam,
                reference=ref,
                output=out_dir,
                bin_input=bins,
                arms_file=arms,
                ocr_file=ocr,
                wps_file=wps,
                sample_name="test_sample",
                mapq=20,
                minlen=65,
                maxlen=400,
                skip_duplicates=True,
                require_proper_pair=True,
                threads=0,
                exclude_regions=None,
                bisulfite_bam=None,
                variants=None,
                chromosomes=None
            )
        except SystemExit as e:
            if e.code != 0:
                print(f"run_all failed with exit code {e.code}")
                raise e
        except Exception as e:
            print(f"run_all failed with exception: {e}")
            raise e
            
        # Verify Outputs
        out_dir = Path(out_dir)
        # Note: FSC output standardization to .tsv renamed these:
        fsc_out = out_dir / "test_sample.FSC.tsv"
        assert fsc_out.exists(), "FSC output missing"
        
        fsr_out = out_dir / "test_sample.FSR.tsv"
        assert fsr_out.exists(), "FSR output missing"
        
        fsd_out = out_dir / "test_sample.FSD.tsv"
        assert fsd_out.exists(), "FSD output missing"
        
if __name__ == "__main__":
    test_run_all_unified()
