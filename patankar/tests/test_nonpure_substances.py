#!/usr/bin/env python3
"""
Tests for non-pure substance handling in loadpubchem.py (v1.51)

Tests both pure compounds (from PubChem) and non-pure substances
(mixtures, polymers, UVCB) from the cache.NonPure database.
"""

import pytest
import sys
import os

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from patankar import loadpubchem
from patankar.loadpubchem import (
    NonPureIndex, NonPureSubstanceError,
    is_nonpure_substance, get_nonpure_index,
    migrant, migrantToxtree
)


class TestNonPureIndex:
    """Tests for NonPureIndex class."""

    def test_index_loads(self):
        """NonPureIndex should load from nonpure_index.json."""
        db = NonPureIndex()
        assert len(db) > 0, "NonPureIndex should contain entries"

    def test_is_nonpure_castor_oil(self):
        """Castor oil (8001-79-4) should be recognized as non-pure."""
        db = NonPureIndex()
        assert db.is_nonpure("8001-79-4"), "Castor oil should be non-pure"

    def test_is_nonpure_formaldehyde(self):
        """Formaldehyde (50-00-0) should NOT be recognized as non-pure."""
        db = NonPureIndex()
        assert not db.is_nonpure("50-00-0"), "Formaldehyde should be pure"

    def test_get_metadata(self):
        """get() should return metadata for non-pure substances."""
        db = NonPureIndex()
        info = db.get("8001-79-4")
        assert info is not None, "Should return metadata"
        assert info.get("name") == "Castor oil"
        assert info.get("type") == "mixture"

    def test_contains_operator(self):
        """'in' operator should work for non-pure substances."""
        db = NonPureIndex()
        assert "8001-79-4" in db
        assert "50-00-0" not in db


class TestHelperFunctions:
    """Tests for helper functions."""

    def test_is_nonpure_substance_helper(self):
        """is_nonpure_substance() helper should work correctly."""
        assert is_nonpure_substance("8001-79-4") is True
        assert is_nonpure_substance("50-00-0") is False

    def test_get_nonpure_index_singleton(self):
        """get_nonpure_index() should return same instance."""
        db1 = get_nonpure_index()
        db2 = get_nonpure_index()
        assert db1 is db2, "Should return same singleton instance"


class TestMigrantPureCompounds:
    """Tests for migrant class with pure compounds."""

    def test_pure_compound_lookup(self):
        """Pure compounds should be looked up from PubChem cache."""
        m = migrant("50-00-0")  # Formaldehyde
        assert m.compound == "50-00-0"
        assert m.M is not None
        assert m.M == pytest.approx(30.026, rel=0.01)

    def test_pure_compound_has_structure(self):
        """Pure compounds should have molecular structure."""
        m = migrant("50-00-0")
        assert m.smiles is not None or m.formula is not None


class TestMigrantNonPureCompounds:
    """Tests for migrant class with non-pure compounds."""

    def test_nonpure_rejected_by_default(self):
        """Non-pure compounds should be rejected with pure_only=True (default)."""
        with pytest.raises(NonPureSubstanceError):
            migrant("8001-79-4")  # Castor oil

    def test_nonpure_allowed_with_flag(self):
        """Non-pure compounds should be allowed with pure_only=False."""
        m = migrant("8001-79-4", pure_only=False, verbose=False)
        assert m.compound == "8001-79-4"
        assert getattr(m, '_is_nonpure', False) is True

    def test_nonpure_has_limited_data(self):
        """Non-pure compounds should have limited data (no M, no SMILES)."""
        m = migrant("8001-79-4", pure_only=False, verbose=False)
        assert m.M is None
        assert m.smiles is None

    def test_nonpure_has_metadata(self):
        """Non-pure compounds should have _nonpure_info attribute."""
        m = migrant("8001-79-4", pure_only=False, verbose=False)
        assert hasattr(m, '_nonpure_info')
        assert m._nonpure_info.get("name") == "Castor oil"


class TestMigrantToxtreePure:
    """Tests for migrantToxtree class with pure compounds."""

    def test_pure_compound_toxtree(self):
        """Pure compounds should have ToxTree analysis."""
        m = migrantToxtree("50-00-0", verbose=False)
        assert m.ToxTree is not None
        assert m.CramerClass is not None
        assert m.TTC is not None


class TestMigrantToxtreeNonPure:
    """Tests for migrantToxtree class with non-pure compounds."""

    def test_nonpure_rejected_by_default(self):
        """Non-pure compounds should be rejected with pure_only=True (default)."""
        with pytest.raises(NonPureSubstanceError):
            migrantToxtree("8001-79-4")

    def test_nonpure_allowed_with_flag(self):
        """Non-pure compounds should be allowed with pure_only=False."""
        m = migrantToxtree("8001-79-4", pure_only=False, verbose=False)
        assert m.compound == "8001-79-4"

    def test_nonpure_no_toxtree(self):
        """Non-pure compounds should have no ToxTree analysis."""
        m = migrantToxtree("8001-79-4", pure_only=False, verbose=False)
        assert m.ToxTree is None
        assert m.CramerClass is None
        assert m.TTC is None


class TestSubstanceCategories:
    """Tests for different non-pure substance categories."""

    @pytest.mark.parametrize("cas,expected_type", [
        ("8001-79-4", "mixture"),      # Castor oil
        ("9002-88-4", "polymer"),      # Polyethylene
        ("65997-06-0", "mixture"),     # Rosins
    ])
    def test_substance_type(self, cas, expected_type):
        """Different non-pure substances should have correct type."""
        db = NonPureIndex()
        info = db.get(cas)
        if info:  # Only test if substance is in our database
            assert info.get("type") == expected_type


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
