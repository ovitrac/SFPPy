"""
survey/cache.py — Persistent Cache for Master Curves and PDFs
=============================================================

Provides persistent caching with:
- MasterCurveCache: Stores computed master curves (NPZ + JSON metadata)
- PdfCache: Stores probability outputs (PDF/CDF)
- Lock mechanism to prevent duplicate computation
- Resume support via manifest

@project: SFPPy/INSERM — Survey-scale exposure estimation
@author: Olivier Vitrac, PhD, HDR
@email: olivier.vitrac@gmail.com
@license: MIT
"""

import os
import json
import time
from pathlib import Path
from typing import Tuple, Optional, Dict, Any, List
from dataclasses import dataclass, asdict

import numpy as np

from survey.fingerprints import fingerprint_physics


# =============================================================================
# Master Curve Key
# =============================================================================

@dataclass(frozen=True)
class MasterCurveKey:
    """
    Key that uniquely identifies one master curve g(Fo) = CF/CP0.

    All parameters that affect the solver must be included.
    """
    polymer: str
    mass_g_mol: float
    D: float
    k: float
    k0: float
    lP_m: float
    Fo_max: float
    h: float
    surface_area: float
    food_volume: float
    contact_temperature_degC: float
    CF0: float
    n_fo: int

    def stable_hash(self) -> str:
        """Compute stable hash for cache lookup."""
        return fingerprint_physics(
            polymer=self.polymer,
            mass_g_mol=self.mass_g_mol,
            D=self.D,
            k=self.k,
            k0=self.k0,
            lP_m=self.lP_m,
            Fo_max=self.Fo_max,
            n_fo=self.n_fo,
            h=self.h,
            surface_area=self.surface_area,
            food_volume=self.food_volume,
            contact_temperature_degC=self.contact_temperature_degC,
            CF0=self.CF0,
        )

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


# =============================================================================
# Master Curve Cache
# =============================================================================

class MasterCurveCache:
    """
    Persistent cache of master curves stored as NPZ + JSON metadata.

    Layout:
        cache_dir/
            master_curves/
                <hash>.npz   (fo_grid, cf_over_cp0)
                <hash>.json  (metadata / human-readable)
                <hash>.lock  (temporary lock during computation)
    """

    def __init__(self, cache_dir: str):
        self.root = Path(cache_dir).expanduser().resolve() / "master_curves"
        self.root.mkdir(parents=True, exist_ok=True)
        self._stats = {"hits": 0, "misses": 0}

    def _paths(self, key: MasterCurveKey) -> Tuple[Path, Path, Path]:
        """Get paths for NPZ, JSON, and lock files."""
        h = key.stable_hash()
        npz = self.root / f"{h}.npz"
        js = self.root / f"{h}.json"
        lock = self.root / f"{h}.lock"
        return npz, js, lock

    def exists(self, key: MasterCurveKey) -> bool:
        """Check if master curve exists in cache."""
        npz, js, _ = self._paths(key)
        return npz.exists() and js.exists()

    def load(self, key: MasterCurveKey) -> Tuple[np.ndarray, np.ndarray]:
        """
        Load master curve from cache.

        Returns
        -------
        Tuple[np.ndarray, np.ndarray]
            (fo_grid, cf_over_cp0)

        Raises
        ------
        FileNotFoundError
            If cache entry doesn't exist.
        """
        npz, _, _ = self._paths(key)
        if not npz.exists():
            raise FileNotFoundError(f"Cache entry not found: {key.stable_hash()}")
        data = np.load(npz)
        self._stats["hits"] += 1
        return data["fo_grid"], data["cf_over_cp0"]

    def save(self, key: MasterCurveKey, fo_grid: np.ndarray, cf_over_cp0: np.ndarray) -> None:
        """
        Save master curve to cache.

        Parameters
        ----------
        key : MasterCurveKey
            Cache key.
        fo_grid : np.ndarray
            Fourier number grid.
        cf_over_cp0 : np.ndarray
            Normalized concentration values.
        """
        npz, js, _ = self._paths(key)
        np.savez_compressed(npz, fo_grid=fo_grid, cf_over_cp0=cf_over_cp0)
        js.write_text(
            json.dumps(key.to_dict(), indent=2, sort_keys=True),
            encoding="utf-8"
        )
        self._stats["misses"] += 1

    def acquire_lock(self, key: MasterCurveKey) -> bool:
        """
        Try to acquire exclusive lock for computation.

        Returns
        -------
        bool
            True if lock acquired, False if already locked.
        """
        _, _, lock = self._paths(key)
        try:
            fd = os.open(str(lock), os.O_CREAT | os.O_EXCL | os.O_WRONLY)
            os.close(fd)
            return True
        except FileExistsError:
            return False

    def release_lock(self, key: MasterCurveKey) -> None:
        """Release computation lock."""
        _, _, lock = self._paths(key)
        try:
            lock.unlink()
        except FileNotFoundError:
            pass

    def wait_for_result(self, key: MasterCurveKey, timeout_s: float = 60.0) -> bool:
        """
        Wait for another process to complete computation.

        Parameters
        ----------
        key : MasterCurveKey
            Cache key.
        timeout_s : float
            Maximum wait time in seconds.

        Returns
        -------
        bool
            True if result appeared, False if timeout.
        """
        start = time.time()
        while time.time() - start < timeout_s:
            if self.exists(key):
                return True
            time.sleep(0.1)
        return False

    @property
    def stats(self) -> Dict[str, int]:
        """Return cache statistics."""
        return dict(self._stats)

    def reset_stats(self) -> None:
        """Reset cache statistics."""
        self._stats = {"hits": 0, "misses": 0}

    def list_entries(self) -> List[str]:
        """List all cache entry hashes."""
        return [p.stem for p in self.root.glob("*.npz")]

    def clear(self) -> int:
        """
        Clear all cache entries.

        Returns
        -------
        int
            Number of entries removed.
        """
        count = 0
        for npz in self.root.glob("*.npz"):
            npz.unlink()
            count += 1
        for js in self.root.glob("*.json"):
            js.unlink()
        for lock in self.root.glob("*.lock"):
            lock.unlink()
        return count


# =============================================================================
# PDF Cache
# =============================================================================

class PdfCache:
    """
    Cache for probability outputs (PDF, CDF, quantiles).

    Layout:
        cache_dir/
            pdf_cache/
                <fingerprint>.npz
                <fingerprint>.json
    """

    def __init__(self, cache_dir: str):
        self.root = Path(cache_dir).expanduser().resolve() / "pdf_cache"
        self.root.mkdir(parents=True, exist_ok=True)

    def _paths(self, fingerprint: str) -> Tuple[Path, Path]:
        npz = self.root / f"{fingerprint}.npz"
        js = self.root / f"{fingerprint}.json"
        return npz, js

    def exists(self, fingerprint: str) -> bool:
        npz, js = self._paths(fingerprint)
        return npz.exists()

    def load(self, fingerprint: str) -> Dict[str, np.ndarray]:
        """Load PDF results from cache."""
        npz, _ = self._paths(fingerprint)
        if not npz.exists():
            raise FileNotFoundError(f"PDF cache entry not found: {fingerprint}")
        return dict(np.load(npz))

    def save(
        self,
        fingerprint: str,
        pdf_bin_centers: np.ndarray,
        pdf: np.ndarray,
        cdf: np.ndarray,
        metadata: Dict[str, Any] = None,
    ) -> None:
        """Save PDF results to cache."""
        npz, js = self._paths(fingerprint)
        np.savez_compressed(
            npz,
            pdf_bin_centers=pdf_bin_centers,
            pdf=pdf,
            cdf=cdf,
        )
        if metadata:
            js.write_text(json.dumps(metadata, indent=2), encoding="utf-8")


# =============================================================================
# Survey State (for Resume)
# =============================================================================

@dataclass
class SurveyState:
    """
    Survey computation state for resume capability.

    Tracks which substances have been computed and current progress.
    """
    survey_id: str
    config_fingerprint: str
    total_substances: int
    completed_substances: List[str]  # List of substance IDs
    pending_substances: List[str]
    timestamp: float = 0.0

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "SurveyState":
        return cls(**d)


class SurveyStateManager:
    """
    Manages survey state for resume capability.

    Layout:
        cache_dir/
            survey_state/
                <survey_id>.json
    """

    def __init__(self, cache_dir: str):
        self.root = Path(cache_dir).expanduser().resolve() / "survey_state"
        self.root.mkdir(parents=True, exist_ok=True)

    def _path(self, survey_id: str) -> Path:
        return self.root / f"{survey_id}.json"

    def save(self, state: SurveyState) -> None:
        """Save survey state."""
        state.timestamp = time.time()
        path = self._path(state.survey_id)
        path.write_text(json.dumps(state.to_dict(), indent=2), encoding="utf-8")

    def load(self, survey_id: str) -> Optional[SurveyState]:
        """Load survey state if exists."""
        path = self._path(survey_id)
        if not path.exists():
            return None
        data = json.loads(path.read_text(encoding="utf-8"))
        return SurveyState.from_dict(data)

    def delete(self, survey_id: str) -> None:
        """Delete survey state."""
        path = self._path(survey_id)
        if path.exists():
            path.unlink()

    def list_states(self) -> List[str]:
        """List all saved survey states."""
        return [p.stem for p in self.root.glob("*.json")]
