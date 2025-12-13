"""
Legacy unit tests for mFSD Python implementation.
These tests are for the original Python implementation which has been
replaced by the Rust backend. The integration test (test_mfsd_integration.py)
now covers the end-to-end functionality.

These tests are kept for reference but marked as skipped.
"""
import pytest

# Skip all tests in this module - replaced by Rust implementation
pytestmark = pytest.mark.skip(reason="mFSD logic moved to Rust backend. Use test_mfsd_integration.py instead.")


# Original tests below for reference
class MockRead:
    def __init__(self, query_sequence, reference_start, reference_end, template_length, mapping_quality=60, is_duplicate=False, is_unmapped=False, is_secondary=False):
        self.query_sequence = query_sequence
        self.reference_start = reference_start
        self.reference_end = reference_end
        self.template_length = template_length
        self.mapping_quality = mapping_quality
        self.is_duplicate = is_duplicate
        self.is_unmapped = is_unmapped
        self.is_secondary = is_secondary
        self.query_length = len(query_sequence)
        
    def get_aligned_pairs(self, matches_only=True):
        pairs = []
        for i in range(len(self.query_sequence)):
            pairs.append((i, self.reference_start + i))
        return pairs


def test_classify_read():
    """Obsolete: classify_read is now in Rust"""
    pass


def test_parse_input_file():
    """Obsolete: parse_input_file is now in Rust"""
    pass


def test_calc_mfsd_integration():
    """Obsolete: Use test_mfsd_integration.py instead"""
    pass
