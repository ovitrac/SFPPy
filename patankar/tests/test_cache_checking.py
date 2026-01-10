#!/usr/bin/env python3
"""
Tests for cache checking methods in loadpubchem.py (v1.52)

Tests the efficient cache status checking methods used for
prioritization workflows.
"""

import pytest
import sys
import os

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from patankar.loadpubchem import (
    CompoundIndex, get_default_index,
    is_substance_cached, batch_check_cached,
    is_nonpure_substance
)


class TestCompoundIndexIsCached:
    """Tests for CompoundIndex.is_cached() method."""

    @pytest.fixture
    def db(self):
        """Get the default CompoundIndex instance."""
        return get_default_index()

    def test_is_cached_returns_dict(self, db):
        """is_cached() should return a dict with expected keys."""
        result = db.is_cached("anything")
        assert isinstance(result, dict)
        assert 'cached' in result
        assert 'source' in result
        assert 'cids' in result
        assert 'has_simple' in result
        assert 'has_full' in result

    def test_is_cached_by_cid_numeric(self, db):
        """Numeric CID queries should check cache file directly."""
        # Get a known cached CID
        cached_cids = db.list_cached_cids()
        if not cached_cids:
            pytest.skip("No cached compounds available")

        cid = cached_cids[0]
        result = db.is_cached(str(cid))

        assert result['cached'] is True
        assert result['source'] == 'cid_file'
        assert cid in result['cids']
        assert result['has_simple'] is True

    def test_is_cached_by_name(self, db):
        """Name queries should check the synonym index."""
        # Try common substances that should be in cache
        test_names = ["limonene", "BHT", "Irganox 1076", "formaldehyde"]

        found_one = False
        for name in test_names:
            result = db.is_cached(name)
            if result['cached']:
                found_one = True
                assert result['source'] == 'index'
                assert len(result['cids']) > 0
                break

        if not found_one:
            pytest.skip("None of the test substances are cached")

    def test_is_cached_unknown_substance(self, db):
        """Unknown substances should return cached=False."""
        result = db.is_cached("unknown_substance_xyz_12345")

        assert result['cached'] is False
        assert result['source'] is None
        assert result['cids'] == []
        assert result['has_simple'] is False
        assert result['has_full'] is False

    def test_is_cached_empty_query(self, db):
        """Empty queries should return cached=False."""
        result = db.is_cached("")
        assert result['cached'] is False

        result = db.is_cached(None)
        assert result['cached'] is False

    def test_is_cached_whitespace_handling(self, db):
        """Queries with whitespace should be handled correctly."""
        cached_cids = db.list_cached_cids()
        if not cached_cids:
            pytest.skip("No cached compounds available")

        cid = cached_cids[0]
        result_clean = db.is_cached(str(cid))
        result_space = db.is_cached(f"  {cid}  ")

        assert result_clean['cached'] == result_space['cached']


class TestCompoundIndexListCachedCids:
    """Tests for CompoundIndex.list_cached_cids() method."""

    @pytest.fixture
    def db(self):
        """Get the default CompoundIndex instance."""
        return get_default_index()

    def test_list_cached_cids_returns_list(self, db):
        """list_cached_cids() should return a list."""
        result = db.list_cached_cids()
        assert isinstance(result, list)

    def test_list_cached_cids_contains_integers(self, db):
        """list_cached_cids() should return integers."""
        result = db.list_cached_cids()
        if result:
            assert all(isinstance(cid, int) for cid in result)

    def test_list_cached_cids_is_sorted(self, db):
        """list_cached_cids() should return sorted CIDs."""
        result = db.list_cached_cids()
        if len(result) > 1:
            assert result == sorted(result)

    def test_list_cached_cids_no_duplicates(self, db):
        """list_cached_cids() should not contain duplicates."""
        result = db.list_cached_cids()
        assert len(result) == len(set(result))


class TestCompoundIndexCacheStats:
    """Tests for CompoundIndex.cache_stats() method."""

    @pytest.fixture
    def db(self):
        """Get the default CompoundIndex instance."""
        return get_default_index()

    def test_cache_stats_returns_dict(self, db):
        """cache_stats() should return a dict with expected keys."""
        result = db.cache_stats()

        assert isinstance(result, dict)
        assert 'total_compounds' in result
        assert 'index_entries' in result
        assert 'simple_files' in result
        assert 'full_files' in result

    def test_cache_stats_values_are_integers(self, db):
        """cache_stats() should return integer values."""
        result = db.cache_stats()

        assert isinstance(result['total_compounds'], int)
        assert isinstance(result['index_entries'], int)
        assert isinstance(result['simple_files'], int)
        assert isinstance(result['full_files'], int)

    def test_cache_stats_values_non_negative(self, db):
        """cache_stats() values should be non-negative."""
        result = db.cache_stats()

        assert result['total_compounds'] >= 0
        assert result['index_entries'] >= 0
        assert result['simple_files'] >= 0
        assert result['full_files'] >= 0

    def test_cache_stats_consistency(self, db):
        """total_compounds should equal simple_files."""
        result = db.cache_stats()
        assert result['total_compounds'] == result['simple_files']


class TestIsSubstanceCached:
    """Tests for is_substance_cached() standalone function."""

    def test_is_substance_cached_returns_dict(self):
        """is_substance_cached() should return a dict with expected keys."""
        result = is_substance_cached("anything")

        assert isinstance(result, dict)
        assert 'cached' in result
        assert 'pure' in result
        assert 'nonpure' in result
        assert 'source' in result
        assert 'cids' in result
        assert 'nonpure_type' in result

    def test_is_substance_cached_pure_compound(self):
        """Pure compounds should be detected correctly."""
        db = get_default_index()
        cached_cids = db.list_cached_cids()

        if not cached_cids:
            pytest.skip("No cached compounds available")

        cid = cached_cids[0]
        result = is_substance_cached(str(cid))

        assert result['cached'] is True
        assert result['pure'] is True
        assert result['nonpure'] is False
        assert 'pure_' in result['source']
        assert cid in result['cids']
        assert result['nonpure_type'] is None

    def test_is_substance_cached_nonpure_compound(self):
        """Non-pure compounds should be detected correctly."""
        # Castor oil (8001-79-4) is in the NonPure cache
        if not is_nonpure_substance("8001-79-4"):
            pytest.skip("Castor oil not in NonPure cache")

        result = is_substance_cached("8001-79-4")

        assert result['cached'] is True
        assert result['pure'] is False
        assert result['nonpure'] is True
        assert result['source'] == 'nonpure'
        assert result['cids'] == []
        assert result['nonpure_type'] is not None

    def test_is_substance_cached_unknown(self):
        """Unknown substances should return cached=False."""
        result = is_substance_cached("unknown_substance_xyz_99999")

        assert result['cached'] is False
        assert result['pure'] is False
        assert result['nonpure'] is False
        assert result['source'] is None
        assert result['cids'] == []
        assert result['nonpure_type'] is None

    def test_is_substance_cached_with_preloaded_index(self):
        """Should accept a pre-loaded CompoundIndex for efficiency."""
        db = get_default_index()
        cached_cids = db.list_cached_cids()

        if not cached_cids:
            pytest.skip("No cached compounds available")

        cid = cached_cids[0]
        result = is_substance_cached(str(cid), pure_index=db)

        assert result['cached'] is True


class TestBatchCheckCached:
    """Tests for batch_check_cached() function."""

    def test_batch_check_cached_returns_dict(self):
        """batch_check_cached() should return a dict with expected keys."""
        result = batch_check_cached(["test1", "test2"])

        assert isinstance(result, dict)
        assert 'cached' in result
        assert 'not_cached' in result
        assert 'pure' in result
        assert 'nonpure' in result
        assert 'stats' in result

    def test_batch_check_cached_stats_structure(self):
        """batch_check_cached() stats should have expected keys."""
        result = batch_check_cached(["test1"])

        stats = result['stats']
        assert 'total' in stats
        assert 'cached_count' in stats
        assert 'not_cached_count' in stats
        assert 'pure_count' in stats
        assert 'nonpure_count' in stats

    def test_batch_check_cached_total_consistency(self):
        """cached + not_cached should equal total queries."""
        queries = ["unknown1", "unknown2", "unknown3"]
        result = batch_check_cached(queries)

        assert result['stats']['total'] == len(queries)
        assert (result['stats']['cached_count'] +
                result['stats']['not_cached_count']) == len(queries)

    def test_batch_check_cached_with_mixed_queries(self):
        """Should correctly categorize a mix of cached and uncached."""
        db = get_default_index()
        cached_cids = db.list_cached_cids()

        if not cached_cids:
            pytest.skip("No cached compounds available")

        # Mix of known cached CID and unknown substance
        queries = [str(cached_cids[0]), "unknown_xyz_99999"]
        result = batch_check_cached(queries)

        assert str(cached_cids[0]) in result['cached']
        assert "unknown_xyz_99999" in result['not_cached']
        assert result['stats']['cached_count'] >= 1
        assert result['stats']['not_cached_count'] >= 1

    def test_batch_check_cached_empty_list(self):
        """Empty query list should return empty results."""
        result = batch_check_cached([])

        assert result['cached'] == []
        assert result['not_cached'] == []
        assert result['stats']['total'] == 0

    def test_batch_check_cached_with_preloaded_index(self):
        """Should accept a pre-loaded CompoundIndex for efficiency."""
        db = get_default_index()

        # This should work without loading the index again
        result = batch_check_cached(["test1", "test2"], pure_index=db)

        assert result['stats']['total'] == 2

    def test_batch_check_cached_lists_are_lists(self):
        """All return lists should be actual lists."""
        result = batch_check_cached(["test"])

        assert isinstance(result['cached'], list)
        assert isinstance(result['not_cached'], list)
        assert isinstance(result['pure'], list)
        assert isinstance(result['nonpure'], list)


class TestCacheCheckingIntegration:
    """Integration tests for cache checking workflow."""

    def test_prioritization_workflow(self):
        """Simulate a typical prioritization workflow."""
        # Get some known cached CIDs
        db = get_default_index()
        cached_cids = db.list_cached_cids()[:3]  # Take up to 3

        if not cached_cids:
            pytest.skip("No cached compounds available")

        # Create a list with cached and uncached substances
        all_substances = [str(cid) for cid in cached_cids] + [
            "unknown_compound_1",
            "unknown_compound_2"
        ]

        # Run batch check
        result = batch_check_cached(all_substances, pure_index=db)

        # Verify we can identify what needs fetching
        to_fetch = result['not_cached']
        already_have = result['cached']

        assert len(to_fetch) == 2  # Our unknown compounds
        assert len(already_have) == len(cached_cids)

        # The unknown compounds should be in not_cached
        assert "unknown_compound_1" in to_fetch
        assert "unknown_compound_2" in to_fetch

    def test_cache_stats_matches_list_count(self):
        """cache_stats total should match list_cached_cids length."""
        db = get_default_index()

        stats = db.cache_stats()
        cids = db.list_cached_cids()

        assert stats['total_compounds'] == len(cids)


class TestEdgeCases:
    """Tests for edge cases and error handling."""

    def test_is_cached_with_special_characters(self):
        """Queries with special characters should not crash."""
        db = get_default_index()

        # These should not raise exceptions
        result = db.is_cached("test-compound")
        assert isinstance(result, dict)

        result = db.is_cached("test/compound")
        assert isinstance(result, dict)

        result = db.is_cached("test compound with spaces")
        assert isinstance(result, dict)

    def test_is_substance_cached_type_safety(self):
        """Non-string queries should be handled gracefully."""
        result = is_substance_cached(123)  # Integer instead of string
        assert result['cached'] is False

        result = is_substance_cached(None)
        assert result['cached'] is False

        result = is_substance_cached(["list"])  # List instead of string
        assert result['cached'] is False

    def test_batch_check_with_duplicates(self):
        """Duplicate queries should be handled correctly."""
        db = get_default_index()
        cached_cids = db.list_cached_cids()

        if not cached_cids:
            pytest.skip("No cached compounds available")

        cid_str = str(cached_cids[0])
        queries = [cid_str, cid_str, cid_str]  # Same query 3 times

        result = batch_check_cached(queries, pure_index=db)

        # Should process all 3 (even though duplicates)
        assert result['stats']['total'] == 3
        assert result['stats']['cached_count'] == 3


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
