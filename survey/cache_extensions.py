"""
survey.cache_extensions
=======================

Backward-compatible extensions to the bilayer master-surface cache
schema (`.survey_cache/bilayer_curves/<hash>.npz`).

Purpose
-------
After the v3 survey-engine consolidation, the bilayer cache
gains four redundant-but-explicit fields that make downstream
disentangling and validation cheap and easy to reason about:

- ``surface_col0`` : ``surface[:, 0]`` — the step-1 marginal slice
  (CF/Cp0 at fo2 = 0).  Pre-sliced so disentangling reads one 1-D
  array instead of loading the full 2-D ``surface``.
- ``D_ref_T1``, ``D_ref_T2`` : the diffusivity of the min-permeability
  reference layer at T1 (and T2 for the two-step variant) — the same
  D used to size ``Fo1_max`` and ``Fo2_max``.  Stored explicitly so
  downstream consumers do *not* have to re-derive it from the cache
  key dict.
- ``l_ref``, ``i_ref`` : the reference layer's thickness and index.
- ``cf_sat_bound`` : the mass-balance ceiling ``V_P · Cp0_max / V_F``
  used by the CP-4 control-point check in
  ``survey/step_disentangling.py``.
- ``schema_version`` : ``"cache_ext_v1"`` marker so legacy entries
  are detected on read and upgraded transparently.

Author:  Olivier Vitrac, PhD, HDR
Contact: olivier.vitrac@gmail.com
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

import numpy as np


SCHEMA_VERSION = "cache_ext_v1"


@dataclass(frozen=True)
class ExtendedCacheEntry:
    """
    Strongly-typed view over a (possibly legacy) bilayer cache file.

    All numerical fields are surfaced unconditionally; legacy entries
    have ``schema_version is None`` and any fields not present in the
    raw NPZ are derived on the fly (see ``load_bilayer_cache``).
    """

    # Core grid + surface (always present)
    fo1_grid: np.ndarray
    fo2_grid: np.ndarray
    surface: np.ndarray

    # Extension fields (always populated by the reader; derived for
    # legacy entries)
    surface_col0: np.ndarray
    D_ref_T1: float
    D_ref_T2: Optional[float]
    l_ref: float
    i_ref: int
    cf_sat_bound: Optional[float]
    schema_version: Optional[str]

    # Raw key dict from the JSON sidecar (for callers that want it)
    key_meta: Dict[str, Any]


# --------------------------------------------------------------------- #
# Reader                                                                #
# --------------------------------------------------------------------- #

def load_bilayer_cache(
    cache_npz: Path,
    *,
    upgrade_on_read: bool = False,
) -> ExtendedCacheEntry:
    """
    Read a bilayer cache entry and surface the extended schema.

    Legacy entries (no ``schema_version``) have missing fields derived
    on the fly from the JSON sidecar's key dict.  When
    ``upgrade_on_read=True``, the derived fields are written back into
    the NPZ atomically, so subsequent reads skip the derivation.

    Parameters
    ----------
    cache_npz : Path
        Path to ``<hash>.npz`` in ``bilayer_curves/``.
    upgrade_on_read : bool, default False
        Re-save the NPZ with the full schema after reading.  Off by
        default so reads stay side-effect-free; the bugfix workflow
        runs a separate ``upgrade_cache_dir(...)`` pass.

    Returns
    -------
    ExtendedCacheEntry
    """
    cache_npz = Path(cache_npz)
    if not cache_npz.exists():
        raise FileNotFoundError(cache_npz)

    sidecar = cache_npz.with_suffix(".json")
    key_meta: Dict[str, Any] = {}
    if sidecar.exists():
        try:
            key_meta = json.loads(sidecar.read_text())
        except Exception:
            key_meta = {}

    raw = np.load(cache_npz, allow_pickle=False)
    keys = set(raw.files)

    if "fo1_grid" not in keys or "fo2_grid" not in keys or "surface" not in keys:
        raise ValueError(
            f"{cache_npz} is not a two-step bilayer cache "
            "(missing fo1_grid / fo2_grid / surface)."
        )

    fo1_grid = np.asarray(raw["fo1_grid"], dtype=float)
    fo2_grid = np.asarray(raw["fo2_grid"], dtype=float)
    surface = np.asarray(raw["surface"], dtype=float)

    # Extension fields — read if present, else derive
    if "schema_version" in keys:
        schema_version = str(raw["schema_version"])
    else:
        schema_version = None

    if "surface_col0" in keys:
        surface_col0 = np.asarray(raw["surface_col0"], dtype=float)
    else:
        surface_col0 = surface[:, 0].copy()

    if "D_ref_T1" in keys:
        D_ref_T1 = float(raw["D_ref_T1"])
        l_ref = float(raw["l_ref"]) if "l_ref" in keys else _l_ref_from_key(key_meta)
        i_ref = int(raw["i_ref"]) if "i_ref" in keys else _i_ref_from_key(key_meta)
    else:
        D_ref_T1, l_ref, i_ref = _ref_layer_from_key(key_meta)

    if "D_ref_T2" in keys:
        D_ref_T2: Optional[float] = float(raw["D_ref_T2"])
    else:
        D_ref_T2 = _D_ref_T2_from_key(key_meta, i_ref)

    cf_sat_bound = float(raw["cf_sat_bound"]) if "cf_sat_bound" in keys else None

    entry = ExtendedCacheEntry(
        fo1_grid=fo1_grid,
        fo2_grid=fo2_grid,
        surface=surface,
        surface_col0=surface_col0,
        D_ref_T1=D_ref_T1,
        D_ref_T2=D_ref_T2,
        l_ref=l_ref,
        i_ref=i_ref,
        cf_sat_bound=cf_sat_bound,
        schema_version=schema_version,
        key_meta=key_meta,
    )

    if upgrade_on_read and schema_version != SCHEMA_VERSION:
        save_extended_cache(cache_npz, entry)

    return entry


# --------------------------------------------------------------------- #
# Writer (used by upgrade pass + by workers_v3 after the bugfix)        #
# --------------------------------------------------------------------- #

def save_extended_cache(cache_npz: Path, entry: ExtendedCacheEntry) -> None:
    """
    Atomically re-save a cache entry with the full extended schema.

    Preserves the original ``surface`` array bit-for-bit; only adds
    the redundant ``surface_col0``, ``D_ref_T*``, ``l_ref``, ``i_ref``,
    ``cf_sat_bound``, ``schema_version`` fields.
    """
    payload: Dict[str, Any] = {
        "fo1_grid":       entry.fo1_grid,
        "fo2_grid":       entry.fo2_grid,
        "surface":        entry.surface,
        "surface_col0":   entry.surface_col0,
        "D_ref_T1":       np.float64(entry.D_ref_T1),
        "l_ref":          np.float64(entry.l_ref),
        "i_ref":          np.int64(entry.i_ref),
        "schema_version": np.array(SCHEMA_VERSION, dtype="<U16"),
    }
    if entry.D_ref_T2 is not None:
        payload["D_ref_T2"] = np.float64(entry.D_ref_T2)
    if entry.cf_sat_bound is not None:
        payload["cf_sat_bound"] = np.float64(entry.cf_sat_bound)

    tmp = cache_npz.with_suffix(".npz.tmp")
    np.savez_compressed(tmp, **payload)
    tmp.replace(cache_npz)


# --------------------------------------------------------------------- #
# Batch upgrade pass                                                    #
# --------------------------------------------------------------------- #

def upgrade_cache_dir(
    bilayer_dir: Path,
    *,
    cf_sat_bound_provider=None,
) -> Tuple[int, int, int]:
    """
    Upgrade every entry in ``<.survey_cache>/bilayer_curves/`` to the
    extended schema.

    Parameters
    ----------
    bilayer_dir : Path
        ``bilayer_curves/`` directory.
    cf_sat_bound_provider : callable(key_meta) -> Optional[float], optional
        Optional hook to compute the mass-balance ceiling per entry
        from the JSON sidecar's key dict.  When omitted, ``cf_sat_bound``
        is left None on legacy entries and recomputable later.

    Returns
    -------
    (n_total, n_upgraded, n_already)
    """
    n_total = 0
    n_upgraded = 0
    n_already = 0
    for npz in sorted(Path(bilayer_dir).glob("*.npz")):
        n_total += 1
        try:
            entry = load_bilayer_cache(npz, upgrade_on_read=False)
        except Exception:
            continue  # not a bilayer surface entry (skip)
        if entry.schema_version == SCHEMA_VERSION:
            n_already += 1
            continue
        # Derive cf_sat_bound if a provider is supplied
        if cf_sat_bound_provider is not None and entry.cf_sat_bound is None:
            try:
                cf_sat = cf_sat_bound_provider(entry.key_meta)
            except Exception:
                cf_sat = None
            entry = ExtendedCacheEntry(**{**entry.__dict__, "cf_sat_bound": cf_sat})
        save_extended_cache(npz, entry)
        n_upgraded += 1
    return n_total, n_upgraded, n_already


# --------------------------------------------------------------------- #
# Internal helpers — derive missing fields from the JSON sidecar        #
# --------------------------------------------------------------------- #

def _ref_layer_from_key(key_meta: Dict[str, Any]) -> Tuple[float, float, int]:
    """
    Apply the min-permeability rule from ``_build_bilayer_tasks``
    (survey_v2.py:262–271) to the cache key dict.

    Returns ``(D_ref_T1, l_ref, i_ref)``.
    """
    if not key_meta:
        raise ValueError(
            "Cannot derive reference layer: cache sidecar JSON missing "
            "or empty.  Re-run with `upgrade_on_read=True` after the "
            "bugfix run so the cache carries D_ref_T1, l_ref, i_ref "
            "explicitly."
        )

    # Two-step keys have D_*_T1; 1-step keys have D_*.
    if "D_1_T1" in key_meta and "D_2_T1" in key_meta:
        D_1 = float(key_meta["D_1_T1"])
        D_2 = float(key_meta["D_2_T1"])
    else:
        D_1 = float(key_meta["D_1"])
        D_2 = float(key_meta["D_2"])

    k_1 = float(key_meta["k_1"])
    k_2 = float(key_meta["k_2"])
    l_1 = float(key_meta["l_1_m"])
    l_2 = float(key_meta["l_2_m"])

    perm_0 = D_1 / (k_1 * l_1) if k_1 * l_1 > 0 else float("inf")
    perm_1 = D_2 / (k_2 * l_2) if k_2 * l_2 > 0 else float("inf")

    if perm_1 < perm_0:
        return D_2, l_2, 1
    else:
        return D_1, l_1, 0


def _l_ref_from_key(key_meta: Dict[str, Any]) -> float:
    return _ref_layer_from_key(key_meta)[1]


def _i_ref_from_key(key_meta: Dict[str, Any]) -> int:
    return _ref_layer_from_key(key_meta)[2]


def _D_ref_T2_from_key(key_meta: Dict[str, Any], i_ref: int) -> Optional[float]:
    """For two-step keys, D_ref_T2 is the same layer's D at T2."""
    if not key_meta:
        return None
    if i_ref == 0 and "D_1_T2" in key_meta:
        return float(key_meta["D_1_T2"])
    if i_ref == 1 and "D_2_T2" in key_meta:
        return float(key_meta["D_2_T2"])
    return None
