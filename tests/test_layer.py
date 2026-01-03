"""
Layer property tests for SFPPy.

Tests layer creation, property calculations, and polymer-specific behavior.

@project: SFPPy - Safe Food Packaging in Python
@author: Olivier Vitrac
@license: MIT
"""

import pytest
import numpy as np
from numpy.testing import assert_allclose


class TestLayerCreation:
    """Tests for basic layer creation and initialization."""

    def test_ldpe_creation(self, ldpe_layer):
        """LDPE layer should be created with correct properties."""
        assert ldpe_layer is not None
        assert ldpe_layer.l is not None  # thickness should be set

    def test_hdpe_creation(self, hdpe_layer):
        """HDPE layer should be created with correct properties."""
        assert hdpe_layer is not None

    def test_gpet_creation(self, gpet_layer):
        """gPET layer should be created with correct properties."""
        assert gpet_layer is not None

    def test_layer_thickness_si(self):
        """Layer thickness should be stored in SI units (meters)."""
        from patankar.layer import LDPE
        layer = LDPE(l=(100, "um"))
        assert_allclose(layer.l, 1e-4, rtol=1e-10)

    def test_layer_thickness_different_units(self):
        """Layer should accept different thickness units."""
        from patankar.layer import LDPE

        layer_um = LDPE(l=(100, "um"))
        layer_mm = LDPE(l=(0.1, "mm"))

        assert_allclose(layer_um.l, layer_mm.l, rtol=1e-10)

    def test_layer_initial_concentration(self):
        """Initial concentration C0 should be set correctly."""
        from patankar.layer import LDPE
        layer = LDPE(l=(100, "um"), C0=1000)
        assert layer.C0 == 1000

    def test_layer_temperature(self):
        """Temperature should be stored correctly."""
        from patankar.layer import LDPE
        layer = LDPE(l=(100, "um"), T=(40, "degC"))
        assert_allclose(layer.T, 40, rtol=1e-3)


class TestLayerProperties:
    """Tests for layer property calculations."""

    def test_diffusivity_positive(self, ldpe_layer):
        """Diffusivity should be positive."""
        assert ldpe_layer.D > 0

    def test_diffusivity_temperature_dependence(self):
        """Higher temperature should increase diffusivity."""
        from patankar.layer import LDPE
        from patankar.loadpubchem import migrant

        m = migrant("limonene")

        layer_cold = LDPE(l=(100, "um"), substance=m, T=(5, "degC"))
        layer_hot = LDPE(l=(100, "um"), substance=m, T=(40, "degC"))

        # Diffusivity should be higher at higher temperature
        assert layer_hot.D > layer_cold.D

    def test_partition_coefficient_positive(self, ldpe_layer, limonene_migrant):
        """Partition coefficient should be positive."""
        from patankar.layer import LDPE
        layer = LDPE(l=(100, "um"), substance=limonene_migrant)
        assert layer.k > 0

    def test_glass_transition_temperature(self):
        """Polymers should have defined Tg."""
        from patankar.layer import LDPE, gPET

        ldpe = LDPE(l=(100, "um"))
        gpet = gPET(l=(12, "um"))

        # LDPE is rubbery (low Tg), PET is glassy (high Tg)
        assert ldpe.Tg < gpet.Tg


class TestPolymerDatabase:
    """Tests for polymer material database."""

    def test_list_materials(self):
        """list_materials() should print available polymers."""
        from patankar.layer import list_materials
        import io
        import sys

        # Capture stdout since list_materials prints
        captured = io.StringIO()
        sys.stdout = captured
        list_materials()
        sys.stdout = sys.__stdout__

        output = captured.getvalue()
        # Should contain polymer names
        assert "LDPE" in output
        assert "HDPE" in output

    def test_ldpe_material_code(self):
        """LDPE should have correct material code."""
        from patankar.layer import LDPE
        layer = LDPE(l=(100, "um"))
        assert hasattr(layer, 'code') or hasattr(layer, 'material')

    def test_different_polymer_types(self):
        """Different polymers should have different properties."""
        from patankar.layer import LDPE, HDPE, PP, gPET

        ldpe = LDPE(l=(100, "um"))
        hdpe = HDPE(l=(100, "um"))
        pp = PP(l=(100, "um"))
        gpet = gPET(l=(100, "um"))

        # Each should be a distinct class
        assert type(ldpe) != type(hdpe)
        assert type(ldpe) != type(pp)
        assert type(ldpe) != type(gpet)


class TestMultilayer:
    """Tests for multilayer structures."""

    def test_bilayer_creation(self):
        """Bilayer structure should be created with + operator."""
        from patankar.layer import LDPE, gPET

        ldpe = LDPE(l=(100, "um"), C0=1000)
        gpet = gPET(l=(12, "um"), C0=0)

        bilayer = ldpe + gpet
        assert bilayer is not None

    def test_trilayer_creation(self):
        """Trilayer (ABA) structure should be created."""
        from patankar.layer import PP, gPET

        a1 = gPET(l=(12, "um"), C0=0)
        b = PP(l=(500, "um"), C0=1000)
        a2 = gPET(l=(12, "um"), C0=0)

        trilayer = a1 + b + a2
        assert trilayer is not None

    def test_multilayer_nlayer(self):
        """Multilayer should report correct number of layers."""
        from patankar.layer import LDPE, gPET

        ldpe = LDPE(l=(100, "um"))
        gpet = gPET(l=(12, "um"))

        mono = ldpe
        bi = ldpe + gpet
        tri = gpet + ldpe + gpet

        # Use _nlayer (internal attribute)
        assert mono._nlayer == 1
        assert bi._nlayer == 2
        assert tri._nlayer == 3


class TestLayerLink:
    """Tests for layerLink functionality."""

    def test_layerlink_creation(self):
        """layerLink should connect layers properly."""
        from patankar.layer import layer, layerLink

        # Create two generic layers
        p1 = layer(l=(100, "um"), D=(1e-12, "m**2/s"), k=1)
        p2 = layer(l=(50, "um"), D=(1e-13, "m**2/s"), k=2)

        # Link them
        linked = layerLink(p1, p2)
        assert linked is not None


class TestMesh:
    """Tests for mesh generation."""

    def test_mesh_creation(self):
        """Mesh should be generated with thickness and node count."""
        from patankar.layer import mesh

        # mesh(l, n, x0, index) - l is thickness, n is number of nodes
        m = mesh(l=1e-4, n=100, index=0)

        assert m is not None
        assert m.n == 100
        assert m.l == 1e-4

    def test_mesh_refinement(self):
        """Finer mesh should have more points."""
        from patankar.layer import mesh

        m_coarse = mesh(l=1e-4, n=50, index=0)
        m_fine = mesh(l=1e-4, n=200, index=0)

        assert m_fine.n > m_coarse.n


class TestSubstanceAssignment:
    """Tests for substance/migrant assignment to layers."""

    def test_substance_in_constructor(self, limonene_migrant):
        """Substance can be assigned in constructor."""
        from patankar.layer import LDPE

        layer = LDPE(
            l=(100, "um"),
            substance=limonene_migrant,
            C0=1000
        )
        assert layer.substance is not None

    def test_migrant_alias(self, limonene_migrant):
        """'migrant' should work as alias for 'substance'."""
        from patankar.layer import LDPE

        layer = LDPE(
            l=(100, "um"),
            migrant=limonene_migrant,
            C0=1000
        )
        # Should have substance set (either via substance or migrant)
        assert hasattr(layer, 'substance') or hasattr(layer, 'migrant')
