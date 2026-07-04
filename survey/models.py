"""
survey/models.py — Data Models for Survey Computations
======================================================

Defines the core data structures for survey-scale migration estimation:
- LayerSpec: Single layer specification
- PackagingSpec: Monolayer or multilayer packaging
- SubstanceSpec: Substance with inferred parameters
- PriorSpec: Triangular prior specification
- SurveyConfig: Full survey configuration
- ComponentSpec: Single packaging component (Phase 4)
- MultiComponentJob: Multi-component packaging job (Phase 4)

@project: SFPPy — Survey-scale exposure estimation
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
    weight : float
        Occurrence-based weight for family CDF computation.
        Default 1.0 (equiprobable). Use occurrence values (0.5, 1, 2)
        to weight substances by likelihood of presence.
    exchangeable : bool
        Whether substance is exchangeable within the family.
        - True (default): Alternatives - concentration is fractionated among
          exchangeable substances (represents uncertainty about which is used).
        - False: Always present - uses full family concentration.
        Monomers are non-exchangeable by definition.
    """
    id: str
    mass_g_mol: float
    name: Optional[str] = None
    cas: Optional[str] = None
    D: Optional[float] = None
    k: Optional[float] = None
    k0: Optional[float] = None
    weight: float = 1.0
    exchangeable: bool = True

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
            'weight': self.weight,
            'exchangeable': self.exchangeable,
        }


# =============================================================================
# Prior Specification
# =============================================================================

@dataclass
class PriorSpec:
    """
    Triangular prior specification.

    Represents Triangular(min_val, mode, max_val) distribution
    discretized into a finite grid with cell-integrated weights.

    Attributes
    ----------
    mode : float
        Mode of the triangular distribution.
    max_val : float
        Maximum value of the triangular distribution.
    min_val : float
        Minimum value (lower bound). Default 0.0 for backward compatibility
        with pre-D0 scenarios. Added 2026-04-08 (D0): the solver now honours
        this field. YAML scenarios with min > 0 are discretized from min,
        not from 0.
    n_low : int
        Number of grid points between min_val and mode.
    n_high : int
        Number of grid points between mode and max_val.
    spacing_low : str
        'linear' (default) or 'log'. Log-spaced low segment for wide-span
        priors (D6 of singleton cp0 revision).
    spacing_high : str
        'linear' (default) or 'log'.
    name : str
        Name for display/logging.
    """
    mode: float
    max_val: float
    min_val: float = 0.0
    n_low: int = 15
    n_high: int = 15
    spacing_low: str = 'linear'
    spacing_high: str = 'linear'
    name: str = "prior"

    def to_dict(self) -> Dict[str, Any]:
        return {
            'min_val': self.min_val,
            'mode': self.mode,
            'max_val': self.max_val,
            'n_low': self.n_low,
            'n_high': self.n_high,
            'spacing_low': self.spacing_low,
            'spacing_high': self.spacing_high,
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
    quantity : int
        Number of food units per pack (dbs `Quantity`). >1 marks a shared
        OVERPACK; used by the physical engine to remove an overpack paper/board
        layer at the oven step (WP-physical, OV 2026-06-25). Default 1.
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
    quantity: int = 1

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

    # Mapping from polymer code to layer class
    _polymer_class_map = None

    def __init__(self, food_simulant: str = "oliveoil"):
        self.food_simulant = food_simulant
        self._migrant_cache = {}
        self._init_polymer_map()

    def _init_polymer_map(self):
        """Initialize polymer code to layer class mapping."""
        if SubstanceModel._polymer_class_map is not None:
            return
        from patankar import layer as layer_module
        SubstanceModel._polymer_class_map = {
            # Polyolefins
            "LDPE": layer_module.LDPE,
            "HDPE": layer_module.HDPE,
            "LLDPE": layer_module.LLDPE,
            "PP": layer_module.PP,
            "oPP": getattr(layer_module, 'oPP', layer_module.PP),
            # Polyesters
            "PET": layer_module.gPET,
            "gPET": layer_module.gPET,
            "wPET": layer_module.wPET,
            "rPET": layer_module.rPET,
            "PBT": layer_module.PBT,
            "PEN": layer_module.PEN,
            # Styrenes
            "PS": layer_module.PS,
            "HIPS": layer_module.HIPS,
            "SBS": layer_module.SBS,
            # Acrylics
            "PMMA": layer_module.PMMA,
            "PVAc": layer_module.PVAc,
            # Polyamides
            "PA": layer_module.PA6,
            "PA6": layer_module.PA6,
            "PA66": layer_module.PA66,
        }

    def _get_migrant(self, substance_name: str = None, cas: str = None):
        """
        Get or create cached migrant object.

        Prioritizes CAS lookup over name lookup since CAS numbers are
        more reliable for PubChem queries.

        Parameters
        ----------
        substance_name : str, optional
            Substance name for fallback lookup.
        cas : str, optional
            CAS registry number (preferred).

        Returns
        -------
        migrant
            Migrant object from PubChem.
        """
        from patankar.loadpubchem import migrant

        # Use CAS as primary key if available
        cache_key = cas or substance_name
        if not cache_key:
            raise ValueError("Either cas or substance_name must be provided")

        if cache_key not in self._migrant_cache:
            # Try CAS first, then name
            if cas:
                try:
                    self._migrant_cache[cache_key] = migrant(cas)
                except Exception:
                    if substance_name:
                        self._migrant_cache[cache_key] = migrant(substance_name)
                    else:
                        raise
            else:
                self._migrant_cache[cache_key] = migrant(substance_name)

        return self._migrant_cache[cache_key]

    def infer_D(self, polymer: str, mass_g_mol: float, temperature_degC: float,
                substance_name: str = None, cas: str = None) -> float:
        """Infer diffusivity. Piringer for plastics; for paper/board (NOT Piringer polymers) use the
        SFPPy Cardboard/Paper layer class D; fall back to a conservative 1e-14 if neither resolves.
        Previously Piringer-only and uncaught — it crashed reference-layer selection on paper/board."""
        from patankar.property import Dpiringer
        p = str(polymer)
        if any(t in p.upper() for t in ("PAPER", "BOARD", "CARTON")):
            try:
                from patankar import layer as _lm
                import numpy as _np
                cls = _lm.Cardboard if "BOARD" in p.upper() else _lm.Paper
                m = self._get_migrant(substance_name=substance_name, cas=cas) if (substance_name or cas) else None
                lay = cls(l=1e-4, substance=m, T=temperature_degC) if m is not None else cls(l=1e-4, T=temperature_degC)
                return float(_np.asarray(lay.D).ravel()[0])
            except Exception:
                return 1e-14
        try:
            return float(Dpiringer.evaluate(polymer=polymer, M=mass_g_mol, T=temperature_degC))
        except Exception:
            return 1e-14

    def infer_k(self, polymer: str, mass_g_mol: float, temperature_degC: float,
                substance_name: str = None, cas: str = None) -> float:
        """
        Infer polymer-side Henry-like coefficient k.

        Uses Flory-Huggins theory via the layer class to compute k
        based on polarity indices and molar volumes of the migrant
        and polymer monomer.

        Parameters
        ----------
        polymer : str
            Polymer code (e.g., 'PET', 'LDPE', 'PP').
        mass_g_mol : float
            Molar mass (g/mol).
        temperature_degC : float
            Temperature (°C).
        substance_name : str, optional
            Substance name for migrant lookup (fallback).
        cas : str, optional
            CAS registry number (preferred for lookup).

        Returns
        -------
        float
            Polymer-side Henry-like coefficient k.
        """
        if (substance_name or cas) and polymer:
            try:
                return self._compute_k_from_layer(polymer, substance_name, temperature_degC, cas)
            except Exception:
                pass
        return 1.0

    def infer_k0(self, polymer: str, mass_g_mol: float, temperature_degC: float,
                 substance_name: str = None, cas: str = None) -> float:
        """
        Infer food-side Henry-like coefficient k0.

        Uses Flory-Huggins theory via the food class to compute k0
        based on polarity indices and molar volumes of the migrant
        and food simulant.

        At equilibrium: k * Cp_eq = k0 * CF_eq
        Partition coefficient K = k / k0

        Parameters
        ----------
        polymer : str
            Polymer code (unused, for signature consistency).
        mass_g_mol : float
            Molar mass (g/mol).
        temperature_degC : float
            Temperature (°C).
        substance_name : str, optional
            Substance name for migrant lookup (fallback).
        cas : str, optional
            CAS registry number (preferred for lookup).

        Returns
        -------
        float
            Food-side Henry-like coefficient k0.
        """
        if (substance_name or cas) and self.food_simulant:
            try:
                return self._compute_k0_from_food(substance_name, temperature_degC, cas)
            except Exception:
                pass
        return 1.0

    def _compute_k_from_layer(self, polymer: str, substance_name: str, T: float,
                               cas: str = None) -> float:
        """
        Compute polymer-side k from substance/polymer pair using SFPPy layer.

        Uses Flory-Huggins theory (kFHP) via layer class to compute the
        Henry-like coefficient in the polymer phase.

        Parameters
        ----------
        polymer : str
            Polymer code (e.g., 'PET', 'LDPE').
        substance_name : str
            Substance name for migrant lookup (fallback).
        T : float
            Contact temperature (°C).
        cas : str, optional
            CAS registry number (preferred for lookup).

        Returns
        -------
        float
            Polymer-side Henry-like coefficient k.
        """
        # Get layer class for polymer
        polymer_upper = polymer.upper().strip()
        if polymer_upper not in SubstanceModel._polymer_class_map:
            # Try without case sensitivity
            for key in SubstanceModel._polymer_class_map:
                if key.upper() == polymer_upper:
                    polymer_upper = key
                    break
            else:
                return 1.0

        layer_class = SubstanceModel._polymer_class_map[polymer_upper]

        # Get migrant (CAS preferred)
        m = self._get_migrant(substance_name=substance_name, cas=cas)

        # Create layer with substance - k is computed automatically via FH
        lay = layer_class(substance=m, T=T)

        k_val = lay.k
        if hasattr(k_val, '__len__'):
            k_val = float(k_val[0])
        return float(k_val)

    def _compute_k0_from_food(self, substance_name: str, T: float, cas: str = None) -> float:
        """
        Compute food-side k0 from substance/food simulant pair using SFPPy.

        Uses Flory-Huggins theory (kFHP) via food class to compute the
        Henry-like coefficient in the food phase.

        Parameters
        ----------
        substance_name : str
            Substance name for migrant lookup (fallback).
        T : float
            Contact temperature (°C).
        cas : str, optional
            CAS registry number (preferred for lookup).

        Returns
        -------
        float
            Food-side Henry-like coefficient k0.
        """
        from patankar import food

        # Get migrant (CAS preferred)
        m = self._get_migrant(substance_name=substance_name, cas=cas)

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
        Uses CAS for migrant lookup (more reliable than names).
        Preserves weight and exchangeable attributes from input.
        """
        D = self.infer_D(layer.polymer, substance.mass_g_mol, layer.temperature_degC)
        k = self.infer_k(layer.polymer, substance.mass_g_mol, layer.temperature_degC,
                        substance_name=substance.name, cas=substance.cas)
        k0 = self.infer_k0(layer.polymer, substance.mass_g_mol, layer.temperature_degC,
                          substance_name=substance.name, cas=substance.cas)

        return SubstanceSpec(
            id=substance.id,
            mass_g_mol=substance.mass_g_mol,
            name=substance.name,
            cas=substance.cas,
            D=D,
            k=k,
            k0=k0,
            weight=substance.weight,
            exchangeable=substance.exchangeable,
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


# =============================================================================
# Multi-Component Packaging (Phase 4)
# =============================================================================

@dataclass
class ComponentSpec:
    """
    Single packaging component specification for multi-component jobs.

    A component represents one physical part of a packaging system
    (e.g., bottle body, cap, lid) that contacts the same food.

    Attributes
    ----------
    pack_code : str
        Unique component identifier (e.g., 'E01_S1', 'E01_S2').
    polymer : str
        Polymer identifier (e.g., 'gPET', 'HDPE').
    thickness_m : float
        Component thickness in meters.
    surface_area_m2 : float
        Contact surface area in m².
    temperature_degC : float
        Contact temperature (°C).
    substances : List[SubstanceSpec]
        Substances present in this component.
    """
    pack_code: str
    polymer: str
    thickness_m: float
    surface_area_m2: float
    temperature_degC: float
    substances: List[SubstanceSpec] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        return {
            'pack_code': self.pack_code,
            'polymer': self.polymer,
            'thickness_m': self.thickness_m,
            'surface_area_m2': self.surface_area_m2,
            'temperature_degC': self.temperature_degC,
            'substances': [s.to_dict() for s in self.substances],
        }

    @property
    def n_substances(self) -> int:
        return len(self.substances)


@dataclass
class MultiComponentJob:
    """
    Multi-component packaging job specification.

    Groups multiple packaging components that contact the same food product.
    Used for aggregating migration contributions via tensor product.

    Attributes
    ----------
    name : str
        Job name (typically '{food_code}_combined').
    food_code : str
        Food product identifier (e.g., 'E01').
    components : List[ComponentSpec]
        List of packaging components.
    food_volume_m3 : float
        Total food volume in m³.
    food_simulant : str
        Food simulant name for k inference.
    contact_temperature_degC : float
        Contact temperature (°C).
    time_prior : PriorSpec
        Contact time prior (seconds).
    conc_prior : PriorSpec
        Initial concentration prior (mg/kg).
    h_m_s : float
        Mass transfer coefficient at food/polymer interface (m/s).
    cf0 : float
        Initial food concentration.
    """
    name: str
    food_code: str
    components: List[ComponentSpec]
    food_volume_m3: float
    food_simulant: str
    contact_temperature_degC: float
    time_prior: PriorSpec
    conc_prior: PriorSpec
    h_m_s: float = 1e-7
    cf0: float = 0.0

    @property
    def is_multicomponent(self) -> bool:
        """True if job has multiple components."""
        return len(self.components) > 1

    @property
    def n_components(self) -> int:
        return len(self.components)

    @property
    def pack_codes(self) -> List[str]:
        return [c.pack_code for c in self.components]

    @property
    def polymers(self) -> List[str]:
        return [c.polymer for c in self.components]

    @property
    def total_surface_area_m2(self) -> float:
        return sum(c.surface_area_m2 for c in self.components)

    def get_all_substances(self) -> List[SubstanceSpec]:
        """Get all unique substances across all components."""
        seen = set()
        result = []
        for comp in self.components:
            for sub in comp.substances:
                if sub.id not in seen:
                    seen.add(sub.id)
                    result.append(sub)
        return result

    def get_substance_components(self, substance_id: str) -> List[ComponentSpec]:
        """Get components containing a specific substance."""
        return [c for c in self.components
                if any(s.id == substance_id for s in c.substances)]

    def to_dict(self) -> Dict[str, Any]:
        return {
            'name': self.name,
            'food_code': self.food_code,
            'components': [c.to_dict() for c in self.components],
            'food_volume_m3': self.food_volume_m3,
            'food_simulant': self.food_simulant,
            'contact_temperature_degC': self.contact_temperature_degC,
            'time_prior': self.time_prior.to_dict(),
            'conc_prior': self.conc_prior.to_dict(),
            'h_m_s': self.h_m_s,
            'cf0': self.cf0,
        }

    def to_survey_config(self, component: ComponentSpec) -> SurveyConfig:
        """
        Create a SurveyConfig for a single component.

        Used to run per-component simulations before aggregation.
        """
        layer = LayerSpec(
            polymer=component.polymer,
            thickness_m=component.thickness_m,
            temperature_degC=component.temperature_degC,
        )
        packaging = PackagingSpec(
            layers=[layer],
            h_m_s=self.h_m_s,
            surface_area_m2=component.surface_area_m2,
            food_volume_m3=self.food_volume_m3,
            contact_temperature_degC=self.contact_temperature_degC,
            cf0=self.cf0,
            food_simulant=self.food_simulant,
        )
        return SurveyConfig(
            name=f"{self.name}_{component.pack_code}",
            packaging=packaging,
            time_prior=self.time_prior,
            conc_prior=self.conc_prior,
        )
