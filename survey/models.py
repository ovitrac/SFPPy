"""
survey/models.py — Data Models for Survey Computations
======================================================

Defines the core data structures for survey-scale migration estimation:
- LayerSpec: Single layer specification
- PackagingSpec: Monolayer or multilayer packaging
- SubstanceSpec: Substance with inferred parameters
- PriorSpec: Triangular prior specification
- SurveyConfig: Full survey configuration

@project: SFPPy/INSERM — Survey-scale exposure estimation
@author: Olivier Vitrac, PhD, HDR
@email: olivier.vitrac@gmail.com
@license: MIT
"""

from dataclasses import dataclass, field, asdict
from typing import List, Optional, Dict, Any, Tuple
import hashlib
import json


# =============================================================================
# Layer and Packaging Specifications
# =============================================================================

@dataclass(frozen=True)
class LayerSpec:
    """
    Single layer specification.

    Attributes
    ----------
    polymer : str
        Polymer identifier (e.g., 'PP', 'LDPE', 'gPET', 'wPET').
    thickness_m : float
        Layer thickness in meters.
    C0 : float
        Initial concentration in the layer (mg/kg or arbitrary units).
    temperature_degC : float
        Temperature for D inference (°C).
    """
    polymer: str
    thickness_m: float
    C0: float = 0.0
    temperature_degC: float = 40.0

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class PackagingSpec:
    """
    Packaging specification (monolayer or multilayer).

    Attributes
    ----------
    layers : List[LayerSpec]
        List of layers from food contact side to outer side.
    h_m_s : float
        Mass transfer coefficient at food/polymer interface (m/s).
    surface_area_m2 : float
        Contact surface area (m²).
    food_volume_m3 : float
        Food volume (m³).
    contact_temperature_degC : float
        Contact temperature (°C).
    cf0 : float
        Initial food concentration.
    food_simulant : str
        Food simulant name for k inference.
    """
    layers: List[LayerSpec]
    h_m_s: float = 1e-7
    surface_area_m2: float = 0.06
    food_volume_m3: float = 0.001
    contact_temperature_degC: float = 40.0
    cf0: float = 0.0
    food_simulant: str = "oliveoil"

    @property
    def is_monolayer(self) -> bool:
        return len(self.layers) == 1

    @property
    def n_layers(self) -> int:
        return len(self.layers)

    @property
    def total_thickness_m(self) -> float:
        return sum(layer.thickness_m for layer in self.layers)

    @property
    def food_contact_layer(self) -> LayerSpec:
        """Layer in contact with food (index 0)."""
        return self.layers[0]

    def to_dict(self) -> Dict[str, Any]:
        return {
            'layers': [l.to_dict() for l in self.layers],
            'h_m_s': self.h_m_s,
            'surface_area_m2': self.surface_area_m2,
            'food_volume_m3': self.food_volume_m3,
            'contact_temperature_degC': self.contact_temperature_degC,
            'cf0': self.cf0,
            'food_simulant': self.food_simulant,
        }


# =============================================================================
# Substance Specification
# =============================================================================

@dataclass
class SubstanceSpec:
    """
    Substance specification with inferred parameters.

    Attributes
    ----------
    id : str
        Unique identifier (name, CAS, or generated).
    mass_g_mol : float
        Molar mass (g/mol).
    name : str, optional
        Human-readable name.
    cas : str, optional
        CAS registry number.
    D : float, optional
        Inferred diffusivity (m²/s). Set during inference.
    k : float, optional
        Inferred partition coefficient. Set during inference.
    k0 : float, optional
        Inferred k0 parameter. Set during inference.
    """
    id: str
    mass_g_mol: float
    name: Optional[str] = None
    cas: Optional[str] = None
    D: Optional[float] = None
    k: Optional[float] = None
    k0: Optional[float] = None

    @classmethod
    def from_mass(cls, mass_g_mol: float, idx: int = 0) -> "SubstanceSpec":
        """Create SubstanceSpec from molar mass only."""
        return cls(
            id=f"M{mass_g_mol:.1f}_{idx}",
            mass_g_mol=mass_g_mol,
            name=f"Substance (M={mass_g_mol:.1f})",
        )

    @classmethod
    def from_name(cls, name: str, mass_g_mol: float, cas: str = None) -> "SubstanceSpec":
        """Create SubstanceSpec from name and mass."""
        return cls(
            id=cas or name.lower().replace(" ", "_"),
            mass_g_mol=mass_g_mol,
            name=name,
            cas=cas,
        )

    def canonical_id(self) -> str:
        """Return canonical identifier for uniqueness checking."""
        if self.cas:
            return f"CAS:{self.cas}"
        return f"M:{self.mass_g_mol:.2f}"

    def to_dict(self) -> Dict[str, Any]:
        return {
            'id': self.id,
            'mass_g_mol': self.mass_g_mol,
            'name': self.name,
            'cas': self.cas,
            'D': self.D,
            'k': self.k,
            'k0': self.k0,
        }


# =============================================================================
# Prior Specification
# =============================================================================

@dataclass
class PriorSpec:
    """
    Triangular prior specification.

    Represents Triangular(min=0, mode, max) distribution
    discretized into a finite grid with cell-integrated weights.

    Attributes
    ----------
    mode : float
        Mode of the triangular distribution.
    max_val : float
        Maximum value of the triangular distribution.
    n_low : int
        Number of grid points below mode.
    n_high : int
        Number of grid points above mode.
    name : str
        Name for display/logging.
    """
    mode: float
    max_val: float
    n_low: int = 15
    n_high: int = 15
    name: str = "prior"

    def to_dict(self) -> Dict[str, Any]:
        return {
            'mode': self.mode,
            'max_val': self.max_val,
            'n_low': self.n_low,
            'n_high': self.n_high,
            'name': self.name,
        }


# =============================================================================
# Survey Configuration
# =============================================================================

@dataclass
class SurveyConfig:
    """
    Complete survey configuration.

    Attributes
    ----------
    name : str
        Survey name.
    packaging : PackagingSpec
        Packaging specification.
    time_prior : PriorSpec
        Contact time prior (seconds).
    conc_prior : PriorSpec
        Initial concentration prior (mg/kg).
    pdf_bins : int
        Number of bins for PDF output.
    n_fo : int
        Number of Fourier number grid points.
    fo_max_factor : float
        Factor to extend Fo_max beyond computed value.
    fo_min_floor : float
        Minimum Fo value (floor for log grid).
    cache_dir : str
        Directory for master curve cache.
    """
    name: str
    packaging: PackagingSpec
    time_prior: PriorSpec
    conc_prior: PriorSpec
    pdf_bins: int = 250
    n_fo: int = 200
    fo_max_factor: float = 1.5
    fo_min_floor: float = 1e-15
    cache_dir: str = ".survey_cache"

    def to_dict(self) -> Dict[str, Any]:
        return {
            'name': self.name,
            'packaging': self.packaging.to_dict(),
            'time_prior': self.time_prior.to_dict(),
            'conc_prior': self.conc_prior.to_dict(),
            'pdf_bins': self.pdf_bins,
            'n_fo': self.n_fo,
            'fo_max_factor': self.fo_max_factor,
            'fo_min_floor': self.fo_min_floor,
            'cache_dir': self.cache_dir,
        }


# =============================================================================
# Inference Hooks
# =============================================================================

class SubstanceModel:
    """
    Individualized inference hooks for D, k, k0.

    Replace these methods to branch/debranch assumptions:
    - infer_D: typically Piringer (polymer, M, T)
    - infer_k: partition coefficient (polymer/food)
    - infer_k0: food-side reference parameter

    The default implementation uses:
    - D: Piringer model
    - k: 1.0 (neutral)
    - k0: 1.0 (neutral)
    """

    def __init__(self, food_simulant: str = "oliveoil"):
        self.food_simulant = food_simulant
        self._migrant_cache = {}

    def infer_D(self, polymer: str, mass_g_mol: float, temperature_degC: float) -> float:
        """Infer diffusivity using Piringer model."""
        from patankar.property import Dpiringer
        return float(Dpiringer.evaluate(polymer=polymer, M=mass_g_mol, T=temperature_degC))

    def infer_k(self, polymer: str, mass_g_mol: float, temperature_degC: float,
                substance_name: str = None) -> float:
        """
        Infer partition coefficient k.

        If substance_name is provided, attempts to compute k from
        substance/food simulant pair. Otherwise returns 1.0.
        """
        if substance_name and self.food_simulant:
            try:
                return self._compute_k_from_food(substance_name, temperature_degC)
            except Exception:
                pass
        return 1.0

    def infer_k0(self, polymer: str, mass_g_mol: float, temperature_degC: float,
                 substance_name: str = None) -> float:
        """Infer k0 parameter. Default: 1.0."""
        return 1.0

    def _compute_k_from_food(self, substance_name: str, T: float) -> float:
        """Compute k from substance/food pair using SFPPy food class."""
        from patankar.loadpubchem import migrant
        from patankar import food

        # Cache migrant lookup
        if substance_name not in self._migrant_cache:
            self._migrant_cache[substance_name] = migrant(substance_name)
        m = self._migrant_cache[substance_name]

        # Get food class
        simulant_map = {
            "oliveoil": food.oliveoil,
            "ethanol50": food.ethanol50,
            "water": food.water,
            "isooctane": food.isooctane,
            "ethanol": food.ethanol,
            "tenax": food.tenax,
        }
        name_lower = self.food_simulant.lower().replace(" ", "").replace("-", "")
        if name_lower not in simulant_map:
            return 1.0

        food_class = simulant_map[name_lower]
        F = food_class(substance=m, contacttemperature=(T, "degC"))

        k0_val = F.k0
        if hasattr(k0_val, '__len__'):
            k0_val = float(k0_val[0])
        return float(k0_val)

    def infer_all(self, substance: SubstanceSpec, layer: LayerSpec) -> SubstanceSpec:
        """
        Infer all parameters (D, k, k0) for a substance.

        Returns a new SubstanceSpec with inferred values.
        """
        D = self.infer_D(layer.polymer, substance.mass_g_mol, layer.temperature_degC)
        k = self.infer_k(layer.polymer, substance.mass_g_mol, layer.temperature_degC,
                        substance.name)
        k0 = self.infer_k0(layer.polymer, substance.mass_g_mol, layer.temperature_degC,
                          substance.name)

        return SubstanceSpec(
            id=substance.id,
            mass_g_mol=substance.mass_g_mol,
            name=substance.name,
            cas=substance.cas,
            D=D,
            k=k,
            k0=k0,
        )


# =============================================================================
# Reference Layer Selection (Multilayer)
# =============================================================================

def select_reference_layer(
    layers: List[LayerSpec],
    substance_model: SubstanceModel,
    surrogate_mass: float = 136.23,  # Limonene
) -> Tuple[int, List[Dict[str, Any]]]:
    """
    Select reference layer for multilayer packaging.

    The reference layer is the one with minimum permeability proxy:
        i_ref = argmin( D_i / (k_i * l_i) )

    This identifies the rate-limiting barrier.

    Parameters
    ----------
    layers : List[LayerSpec]
        Layer specifications.
    substance_model : SubstanceModel
        Model for parameter inference.
    surrogate_mass : float
        Surrogate substance mass for D inference (default: limonene).

    Returns
    -------
    Tuple[int, List[Dict]]
        (i_ref, layer_details) where layer_details contains
        D, k, permeability_proxy for each layer.
    """
    import numpy as np

    layer_details = []
    permeability_proxies = []

    for i, layer in enumerate(layers):
        # Infer D using surrogate
        D = substance_model.infer_D(layer.polymer, surrogate_mass, layer.temperature_degC)
        k = 1.0  # Assume k=1 for same substance in all layers
        l_m = layer.thickness_m

        # Permeability proxy
        if k > 0 and l_m > 0:
            perm_proxy = D / (k * l_m)
        else:
            perm_proxy = float('inf')

        permeability_proxies.append(perm_proxy)
        layer_details.append({
            'index': i,
            'polymer': layer.polymer,
            'thickness_m': l_m,
            'D': D,
            'k': k,
            'permeability_proxy': perm_proxy,
        })

    i_ref = int(np.argmin(permeability_proxies))
    return i_ref, layer_details
