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
    def ocf_regions(self) -> Path:
        """Open Chromatin Regions BED (currently hg19/GRCh37 only)"""
        return self._get_path("OpenChromatinRegion", "7specificTissue.all.OC.bed.gz")
    
    @property
    def ocf_available(self) -> bool:
        """Check if OCF regions are available for this genome"""
        return self.ocf_regions.exists()
    
    @property
    def tfbs_regions(self) -> Path:
        """GTRD Top 5k TFBS metaclusters for size entropy analysis"""
        filename = f"Homo_sapiens_meta_clusters_{self.file_prefix}_midpoint_top5k_sorted.bed.gz"
        return self._get_path("TFBS", filename)
    
    @property
    def tfbs_available(self) -> bool:
        """Check if TFBS regions are available for this genome"""
        return self.tfbs_regions.exists()
    
    @property
    def atac_regions(self) -> Path:
        """TCGA ATAC-seq cancer peak atlas for size entropy analysis"""
        filename = f"TCGA_ATAC_peak.{self.file_prefix}.bed.gz"
        return self._get_path("ATAC", filename)
    
    @property
    def atac_available(self) -> bool:
        """Check if ATAC regions are available for this genome"""
        return self.atac_regions.exists()

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
        """Methylation markers BED for UXM analysis (Atlas U250 markers)"""
        return self._get_path("MethMark", f"Markers.U250.{self.file_prefix}.bed.gz")
    
    # =========================================================================
    # Assay-aware asset resolution (for MSK-ACCESS panels)
    # =========================================================================
    
    def get_gene_bed(self, assay: str) -> Path:
        """
        Get bundled gene BED file for a specific assay.
        
        Args:
            assay: Assay code (xs1, xs2)
            
        Returns:
            Path to gene BED file
            
        Raises:
            FileNotFoundError: If the asset doesn't exist
        """
        path = self._get_path("genes", f"{assay}.genes.bed.gz")
        if not path.exists():
            raise FileNotFoundError(f"Gene BED not found for assay '{assay}': {path}")
        return path
    
    def get_target_bed(self, assay: str) -> Path:
        """
        Get bundled target regions BED for a specific assay.
        
        Args:
            assay: Assay code (xs1, xs2)
            
        Returns:
            Path to target BED file
        """
        path = self._get_path("targets", f"{assay}.targets.bed")
        if not path.exists():
            raise FileNotFoundError(f"Target BED not found for assay '{assay}': {path}")
        return path
    
    def get_wps_anchors(self, assay: str = None) -> Path:
        """
        Get WPS anchors BED, optionally filtered for a specific assay.
        
        Args:
            assay: Optional assay code (xs1, xs2). If None, returns genome-wide anchors.
            
        Returns:
            Path to WPS anchors file
        """
        if assay:
            path = self._get_path("WpsAnchors", f"{assay}.wps_anchors.bed.gz")
            if path.exists():
                return path
            logger.warning(f"Panel-specific WPS anchors not found for '{assay}', using genome-wide")
        return self.wps_anchors
    
    def get_pon(self, assay: str) -> Path:
        """
        Get bundled PON model for a specific assay.
        
        Args:
            assay: Assay code (xs1, xs2)
            
        Returns:
            Path to PON parquet file
        """
        # New naming convention: pon/GRCh37/xs1.pon.parquet
        path = self._get_path("pon", f"{assay}.pon.parquet")
        if path.exists():
            return path
        raise FileNotFoundError(f"PON not found for assay '{assay}': {path}")
    
    def list_available_assays(self) -> list:
        """List all available assays based on bundled gene BED files."""
        genes_dir = self.base_path / "genes" / self.genome_dir
        if not genes_dir.exists():
            return []
        return [
            p.stem.replace(".genes.bed", "").replace(".genes", "")
            for p in genes_dir.glob("*.genes.bed.gz")
        ]
        
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
            "bins_100kb": self.bins_100kb,
            "exclude_regions": self.exclude_regions,

            "ocf_regions": self.ocf_regions,
            "methylation_markers": self.methylation_markers,
            "wps_anchors": self.wps_anchors,
            "wps_background": self.wps_background,
            "tfbs_regions": self.tfbs_regions,
            "atac_regions": self.atac_regions,
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
    
    def validate(self, assay: Optional[str] = None) -> dict:
        """
        Validate bundled asset file formats.
        
        Validates that bundled data files match expected schemas.
        Use this to catch issues if someone modifies the data files.
        
        Args:
            assay: Optional assay code to also validate assay-specific assets
            
        Returns:
            Dict mapping asset name to validation result (True/False)
            
        Raises:
            ValueError: If any bundled file has invalid format (with details)
        """
        from .core.asset_validation import validate_file, FileSchema
        
        results = {}
        
        # Define asset-to-schema mapping for bundled files
        genome_assets = [
            ("arms", self.arms, FileSchema.ARMS_BED),
            ("bins_100kb", self.bins_100kb, FileSchema.BED3),
            ("wps_anchors", self.wps_anchors, FileSchema.WPS_ANCHORS),
            ("wps_background", self.wps_background, FileSchema.WPS_ANCHORS),
            ("ocf_regions", self.ocf_regions, FileSchema.REGION_BED),
            ("tfbs_regions", self.tfbs_regions, FileSchema.REGION_BED),
            ("atac_regions", self.atac_regions, FileSchema.REGION_BED),
        ]
        
        logger.info(f"Validating bundled assets for genome: {self.genome_dir}")
        
        for name, path, schema in genome_assets:
            if path.exists():
                try:
                    validate_file(path, schema)
                    results[name] = True
                except ValueError as e:
                    logger.error(f"Bundled asset validation failed: {name}")
                    raise ValueError(f"Bundled asset '{name}' has invalid format: {e}")
            else:
                results[name] = None  # Not available
        
        # Validate assay-specific assets if requested
        if assay:
            assay_assets = [
                (f"gene_bed_{assay}", self.get_gene_bed(assay), FileSchema.GENE_BED),
                (f"target_bed_{assay}", self.get_target_bed(assay), FileSchema.TARGETS_BED),
            ]
            
            for name, path, schema in assay_assets:
                try:
                    if path and path.exists():
                        validate_file(path, schema)
                        results[name] = True
                    else:
                        results[name] = None
                except Exception:
                    results[name] = None
        
        passed = sum(1 for v in results.values() if v is True)
        skipped = sum(1 for v in results.values() if v is None)
        logger.info(f"Validation complete: {passed} passed, {skipped} not available")
        
        return results

