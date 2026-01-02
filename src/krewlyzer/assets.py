from pathlib import Path
from enum import Enum
from typing import Optional
import logging

logger = logging.getLogger("krewlyzer.assets")

class Genome(str, Enum):
    HG19 = "hg19"
    GRCH37 = "GRCh37"
    HG38 = "hg38"
    GRCH38 = "GRCh38"

class AssetManager:
    """
    Manages access to bundled data assets for specific genomes.
    Resolves paths based on reorganization: data/{Category}/{Genome}/{Filename}
    """
    def __init__(self, genome: str):
        self.raw_genome = genome
        if genome.lower() in ("hg19", "grch37"):
            self.genome_dir = "GRCh37"
            self.file_prefix = "hg19"
        elif genome.lower() in ("hg38", "grch38"):
            self.genome_dir = "GRCh38"
            self.file_prefix = "hg38"
        else:
            raise ValueError(f"Unsupported genome: {genome}. Must be hg19/GRCh37 or hg38/GRCh38.")
            
        # Assets are located in krewlyzer/data/...
        self.base_path = Path(__file__).parent / "data"

    def _get_path(self, category_dir: str, filename: str) -> Path:
        path = self.base_path / category_dir / self.genome_dir / filename
        return path

    @property
    def arms(self) -> Path:
        """Chromosome arms BED"""
        return self._get_path("ChromosomeArms", f"{self.file_prefix}.arms.bed.gz")

    @property
    def valid_regions(self) -> Path:
        """Valid regions BED for GC correction (generated)"""
        return self._get_path("gc", f"valid_regions_{self.file_prefix}.bed.gz")
    
    @property
    def gc_reference(self) -> Path:
        """GC reference Parquet for GC correction"""
        return self._get_path("gc", f"ref_genome_GC_{self.file_prefix}.parquet")

    @property
    def bins_100kb(self) -> Path:
        """100kb fixed windows BED"""
        return self._get_path("ChromosomeBins", f"{self.file_prefix}_window_100kb.bed.gz")
        
    @property
    def exclude_regions(self) -> Path:
        """Blacklist/exclude regions BED"""
        return self._get_path("exclude-regions", f"{self.file_prefix}-blacklist.v2.bed.gz")
        
    @property
    def transcript_anno(self) -> Path:
        """Transcript annotation (default 1kb)"""
        return self._get_path("TranscriptAnno", f"transcriptAnno-{self.file_prefix}-1kb.tsv")
    
    @property
    def ocf_regions(self) -> Path:
        """Open Chromatin Regions BED (currently hg19 specific filename used for both if missing hg38 specific?)"""
        # TODO: Add hg38 specific OCF file if available
        return self._get_path("OpenChromatinRegion", "7specificTissue.all.OC.bed.gz")

    @property
    def wps_anchors(self) -> Path:
        """WPS anchors BED.GZ (merged TSS+CTCF) for dual-stream chromatin profiling"""
        return self._get_path("WpsAnchors", f"{self.file_prefix}.wps_anchors.bed.gz")

    @property
    def wps_background(self) -> Path:
        """WPS background BED.GZ (full-length Alu elements) for global stacking"""
        return self._get_path("WpsBackground", f"{self.file_prefix}.alu_consensus.bed.gz")

    @property
    def methylation_markers(self) -> Path:
        """Methylation markers BED for UXM analysis"""
        # Note: Currently hg19 only - hg38 markers need to be generated
        return self.base_path / "methylation-markers" / f"uxm_markers_{self.file_prefix}.bed"
        
    def resolve(self, asset_name: str) -> Path:
        """
        Resolve an asset by name and validate it exists.
        
        Args:
            asset_name: Name of the asset (e.g., 'bins_100kb', 'gc_reference', 'valid_regions')
            
        Returns:
            Path to the asset
            
        Raises:
            FileNotFoundError: If the asset doesn't exist
            ValueError: If the asset name is unknown
        """
        asset_map = {
            "arms": self.arms,
            "valid_regions": self.valid_regions,
            "gc_reference": self.gc_reference,
            "gc_correction": self.gc_reference,  # Alias
            "bins_100kb": self.bins_100kb,
            "exclude_regions": self.exclude_regions,
            "transcript_anno": self.transcript_anno,
            "ocf_regions": self.ocf_regions,
            "methylation_markers": self.methylation_markers,
            "wps_anchors": self.wps_anchors,
            "wps_background": self.wps_background,
        }
        
        if asset_name not in asset_map:
            raise ValueError(f"Unknown asset: {asset_name}. Available: {list(asset_map.keys())}")
        
        path = asset_map[asset_name]
        if not path.exists():
            raise FileNotFoundError(f"Asset '{asset_name}' not found at: {path}")
        
        return path
        
    def check_exists(self, path: Path, description: str) -> bool:
        if not path.exists():
            logger.warning(f"Missing bundled asset for {description}: {path}")
            return False
        return True
