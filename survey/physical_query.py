"""Physical-primary survey query (WP-physical P2).

Computes the family CF tensors by evaluating the real-units chain at the scenario's ACTUAL discretised
time-prior points (`chain_engine.chain_at_times`) — NO Fo-surface, NO interpolation over time. Reuses
the survey's own family-combine (`_combine_substance_curves` / `_combine_substance_surfaces`), the
per-substance contribution factor (`_per_substance_factors`) and the CP0 prior, so the ONLY difference
from the legacy Fo-surface query is that `g` is exact at the grid (no interp) and real-units (no Fo-axis
normalisation — removes the T5 bug class).

Three drop-in functions, flag-gated by `BilayerSurvey._physical_engine`:
  physical_cf_tensor(sv)            -> family CF tensor (== _compute_cf_tensor, full key parity)
  physical_per_substance_cf_tensors -> per-CAS decomposition (== per_substance_cf_tensors)
  physical_step_resolved_cf_tensors -> CF1/CF2/CF12 step-resolved (== step_resolved_cf_tensors)

Reheat is NOT capped here (the query must reproduce the pipeline for the scenario as given; the food
≤100 °C cap is a separate scenario-design correction). @author: Olivier Vitrac — MIT.
"""
from __future__ import annotations
from typing import Dict, Any, List
import os, json, hashlib
from pathlib import Path
import numpy as np

_T2_FLOOR_S = 1e-3   # tiny storage time standing in for "oven alone" (no storage)

# PET conservatism (OV 2026-06-25). For the Piringer D, model PET as rPET (RUBBERY PET) — the conservative
# Piringer parameter set: D ~27x glassy gPET (1.18e-17 vs 4.36e-19 m^2/s @20degC, M~226). gPET (glassy) is
# the LEAST conservative (lowest D) and is therefore NOT the default. This stays WITHIN the "D = Piringer"
# rule — it is NOT the native free-volume (hFV) PET D, which is valid for TOLUENE ONLY (the retracted
# "1e10" alarm). wPET is the hFV/toluene variant and is deliberately NOT used here.
#   k: the FH partition is Tg-independent (gPET.k == wPET.k), so k is resolved via the GLASSY class gPET —
#   rPET's native .k errors (an accented `_chemicalsubstance` string fails the PubChem lookup), so k is
#   NOT routed through rPET. Hence D and k normalise to DIFFERENT PET classes.
# Scenario codes "PET"/"gPET" (converted raw "PET"); explicit "rPET"/"wPET" pass through. Physical path only.
_PET_FOR_D = {"PET": "rPET", "GPET": "rPET"}     # conservative rubbery Piringer D
_PET_FOR_K = {"PET": "gPET", "RPET": "gPET"}     # FH k via the glassy class whose name resolves


def _norm_polymer_D(polymer: str) -> str:
    """PET-family → rPET for the conservative rubbery Piringer D. Non-PET codes pass through."""
    return _PET_FOR_D.get(str(polymer).upper().strip(), polymer)


def _norm_polymer_k(polymer: str) -> str:
    """PET-family k via gPET (same FH k as any PET variant; rPET's representative-substance name is
    broken). Non-PET codes pass through."""
    return _PET_FOR_K.get(str(polymer).upper().strip(), polymer)


# Overpack removal (WP-6, OV 2026-06-25). A paper/board layer in a pack holding MORE THAN ONE food unit
# (dbs Quantity > 1) is a shared OVERPACK: it is present during storage but NOT carried into the oven
# (it is left behind with the carton). At the oven (step-2) the overpack layer(s) are dropped; the
# retained layers + the food carry forward (the overpack's absorbed mass leaves with it — discarded, not
# redistributed). This generalises the converter's manual `t2_mode=0` trick (which can only blank the
# WHOLE oven step, i.e. the monolayer case) to the multilayer case (drop the board, keep the primary
# plastic migrating in the oven). Physical-engine path only (flag-gated).
_OVERPACK_POLYMERS = ("PAPER", "PAPERBOARD", "CARDBOARD", "BOARD")


def _is_overpack_polymer(polymer: str) -> bool:
    s = str(polymer).upper().strip()
    return any(tok in s for tok in ("PAPER", "BOARD", "CARTON"))

# -------------------------------------------------------------------------------------------------
# P3 cache: g (CF/CP0) keyed by (physical params + time grid). g is CP0-independent, so the CP0 prior
# is applied post-hoc and is NOT part of the key — reuse over CP0 is automatic. In-memory (per-process)
# layer over an atomic .npz disk cache under <cache_dir>/physical_g/. Concurrent workers writing the
# same key are idempotent (same data); atomic os.replace avoids torn reads.
# -------------------------------------------------------------------------------------------------
_GVER = "phys-g-v3"   # v3: cache also carries the genuine storage-step C0eq + default flag
_MEM: Dict[str, Any] = {}
# Parallel cache (keyed identically to _MEM): the genuine storage-step equilibrium
# scale C0eq and whether the chain defaulted C0eq (resume/open step). g itself is
# C0eq-independent (g = CF/CP0), so this is metadata for the guardrail, persisted
# to the NPZ — NOT part of g.
_C0EQ: Dict[str, Any] = {}
_STATS = {"mem_hit": 0, "disk_hit": 0, "miss": 0}


def cache_stats() -> Dict[str, int]:
    return dict(_STATS)


def reset_cache_stats() -> None:
    for k in _STATS:
        _STATS[k] = 0


def _spec_sig(specs):
    return [[float(L.l_m), float(L.D_T1), (None if L.D_T2 is None else float(L.D_T2)),
             float(L.k), float(L.C0)] for L in specs]


def _chain_key(ctx, specs, t1_grid) -> str:
    s2 = _step2_specs(ctx, specs)                       # overpack-removed oven layer set (or None)
    obj = {
        "layers": _spec_sig(specs),
        "s2": (None if s2 is None else _spec_sig(s2)),  # distinguishes overpack-removed scenarios
        "k0": ctx["k0"], "h": ctx["h"], "S": ctx["S"], "V": ctx["V"], "cf0": ctx["cf0"],
        "T1": ctx["T1"], "T2": ctx["T2"], "focal": ctx["focal"], "two_step": ctx["two_step"],
        "t1": [round(float(x), 6) for x in np.asarray(t1_grid, float)],
        "t2": ([round(float(x), 6) for x in np.asarray(ctx["time2_vals"], float)]
               if ctx["two_step"] else None),
        "ver": _GVER,
    }
    return hashlib.sha256(json.dumps(obj, sort_keys=True, default=str).encode()).hexdigest()[:24]


def _cached_chain(sv, ctx, sub, t1_grid):
    """Chain output for one substance at the given t1 grid, memoised + disk-cached.
    Returns (g1, g12) for two-step, g1 for single-step (== chain_at_times)."""
    from survey.chain_engine import chain_at_times
    specs = _layer_specs(sv, ctx, sub)
    key = _chain_key(ctx, specs, t1_grid)
    if key in _MEM:
        _STATS["mem_hit"] += 1
        return _MEM[key]
    cache_dir = getattr(sv.config, "cache_dir", None)
    path = Path(cache_dir) / "physical_g" / f"{key}.npz" if cache_dir else None
    if path is not None and path.exists():
        try:
            d = np.load(path)
            val = (d["g1"], d["g12"]) if "g12" in d.files else d["g1"]
            _MEM[key] = val; _STATS["disk_hit"] += 1
            if "c0eq" in d.files:
                _C0EQ[key] = (float(d["c0eq"]), bool(d["c0eq_default"]))
            return val
        except Exception:
            pass  # corrupt/partial → recompute
    _STATS["miss"] += 1
    common = dict(layers=specs, k0=ctx["k0"], h=ctx["h"], surface_area=ctx["S"],
                  food_volume=ctx["V"], T1_degC=ctx["T1"], CF0=ctx["cf0"], focal_layer=ctx["focal"])
    c0eq_out = {}
    if ctx["two_step"]:
        g1, g12 = chain_at_times(t1_grid=t1_grid, T2_degC=ctx["T2"], t2_grid=ctx["time2_vals"],
                                 layers_step2=_step2_specs(ctx, specs), c0eq_out=c0eq_out, **common)
        val = (np.asarray(g1, float), np.asarray(g12, float))
    else:
        val = np.asarray(chain_at_times(t1_grid=t1_grid, c0eq_out=c0eq_out, **common), float)
    c0eq = float(c0eq_out.get("C0eq", 1.0))
    c0eq_def = bool(c0eq_out.get("C0eq_is_default", False))
    _C0EQ[key] = (c0eq, c0eq_def)
    if path is not None:
        try:
            path.parent.mkdir(parents=True, exist_ok=True)
            tmp = path.parent / f"{path.name}.tmp.{os.getpid()}.npz"   # ends .npz → savez won't append
            if ctx["two_step"]:
                np.savez(tmp, g1=val[0], g12=val[1], c0eq=c0eq, c0eq_default=c0eq_def)
            else:
                np.savez(tmp, g1=val, c0eq=c0eq, c0eq_default=c0eq_def)
            os.replace(tmp, path)
        except Exception:
            pass  # cache write is best-effort
    _MEM[key] = val
    return val


def _resolve(sv):
    """Common context: priors, temperatures, food params, focal — from the survey."""
    from survey.priors import discretize_prior
    pkg = sv.config.packaging
    ctx = {
        "pkg": pkg, "layers_geom": pkg.layers, "focal": int(getattr(sv, "_focal_layer", 0)),
        "subs": sv._substances, "k0": float(sv._substances[0].k0 or 1.0),
        "h": float(pkg.h_m_s), "S": float(pkg.surface_area_m2),
        "V": float(pkg.food_volume_m3), "cf0": float(pkg.cf0),
        "quantity": int(getattr(sv.config, "quantity", 1) or 1),
        # Source layer = where C_P0/C_O0 starts at t=0 (default = focal/food-contact). For an
        # OVERPACK-ORIGIN scenario the source sits on the OUTER (overpack) layer while the food still
        # contacts `focal`; the converter sets `_source_layer` = last and `_overpack_removed` = True
        # (NOTE_overpack_geometry §23.0). Negative indices allowed (Python convention).
        "source_layer": getattr(sv, "_source_layer", None),
        "overpack_removed": bool(getattr(sv, "_overpack_removed", False)),
    }
    # Per-polymer Henry coefficient k via kFHP (assign substance → layer.k per polymer), not a single
    # substance.k for every layer (the master-curve / monolayer holdover). D stays Piringer (the
    # conservative default; SFPPy's hFV for PET is toluene-only, so we do NOT use the gPET-class D).
    from survey.models import SubstanceModel
    ctx["sm"] = SubstanceModel(food_simulant=getattr(pkg, "food_simulant", "water"))
    ctx["time_vals"], ctx["time_w"] = discretize_prior(sv.config.time_prior)
    ctx["conc_vals"], ctx["conc_w"] = discretize_prior(sv.config.conc_prior)
    ctx["two_step"] = bool(sv.is_twostep and sv._step2_time_prior)
    if ctx["two_step"]:
        ctx["T1"] = float(sv._steps[0]["temperature_degC"])
        ctx["T2"] = float(sv._steps[1]["temperature_degC"])
        ctx["time2_vals"], ctx["time2_w"] = discretize_prior(sv._step2_time_prior)
    else:
        ctx["T1"] = float(pkg.contact_temperature_degC); ctx["T2"] = None
    return ctx


def _layer_k(ctx, sub, polymer):
    """Per-polymer Henry coefficient via kFHP (= layer(substance=migrant).k); falls back to the
    single substance.k if the migrant can't be resolved (offline / no cache). PET-family k via the
    glassy class gPET (FH k is Tg-independent; rPET's name is broken)."""
    polymer = _norm_polymer_k(polymer)
    cas = getattr(sub, "cas", None); name = getattr(sub, "name", None)
    mass = float(getattr(sub, "mass_g_mol", 0) or 0)
    try:
        return float(ctx["sm"].infer_k(polymer, mass, ctx["T1"], substance_name=name, cas=cas))
    except Exception:
        return float(sub.k or 1.0)


def _layer_D(sub, polymer, T):
    """Diffusivity for one layer. PET-family → rPET Piringer (conservative, §D4). Paper/board are NOT
    Piringer polymers, so their D comes from the SFPPy Cardboard/Paper layer class (FH/hFV) — needed for
    the overpack source. Everything else uses Piringer. Falls back to substance.D|1e-14."""
    if _is_overpack_polymer(polymer):                       # paper / board: use the SFPPy class D
        try:
            from patankar import layer as lm
            from patankar.loadpubchem import migrant
            cls = lm.Cardboard if "BOARD" in str(polymer).upper() else lm.Paper
            mig = migrant(getattr(sub, "cas", None) or getattr(sub, "name", None))
            return float(np.asarray(cls(l=1e-4, substance=mig, T=T).D).ravel()[0])
        except Exception:
            return float(getattr(sub, "D", None) or 1e-14)
    from patankar.property import Dpiringer
    try:
        return float(Dpiringer.evaluate(_norm_polymer_D(polymer), sub.mass_g_mol, T))
    except Exception:
        return float(getattr(sub, "D", None) or 1e-14)


def _layer_specs(sv, ctx, sub):
    from survey.chain_engine import LayerSpec
    n = len(ctx["layers_geom"])
    src = ctx["focal"] if ctx.get("source_layer") is None else (ctx["source_layer"] % n)  # C0 placement
    out = []
    for i, L in enumerate(ctx["layers_geom"]):
        out.append(LayerSpec(
            l_m=float(L.thickness_m), D_T1=_layer_D(sub, L.polymer, ctx["T1"]),
            k=_layer_k(ctx, sub, L.polymer), C0=(1.0 if i == src else 0.0),
            D_T2=(_layer_D(sub, L.polymer, ctx["T2"]) if ctx["two_step"] else None)))
    return out


def _overpack_indices(ctx):
    """Indices of overpack layers to drop at the oven step.
    Preferred (role-based, NOTE §23): an OVERPACK-ORIGIN scenario sets `_overpack_removed=True` and
    `_source_layer` = the overpack (outer) layer ⇒ drop that layer. Legacy heuristic (paper/board when
    Quantity>1) kept as a fallback when no role is set."""
    n = len(ctx["layers_geom"])
    if ctx.get("overpack_removed") and ctx.get("source_layer") is not None:
        return [ctx["source_layer"] % n]                      # drop the overpack (source) layer
    if not ctx["two_step"] or ctx["quantity"] <= 1:
        return []
    return [i for i, L in enumerate(ctx["layers_geom"]) if _is_overpack_polymer(L.polymer)]


def _step2_specs(ctx, full_specs):
    """Retained LayerSpec list for the OVEN step after removing overpack layer(s); None if no removal.
    Returns [] (empty) for a monolayer overpack (nothing left in the oven)."""
    drop = _overpack_indices(ctx)
    if not drop:
        return None
    drop = set(drop)
    return [s for i, s in enumerate(full_specs) if i not in drop]


def reference_layer_record(sv, sub):
    """The most-impervious layer (lowest permeability = D/(l·k)) for a substance, using per-polymer k —
    stored for review (OV: the layer object exposes `permeability`). Returns {index, polymer, permeability}."""
    ctx = _resolve(sv)
    specs = _layer_specs(sv, ctx, sub)
    perms = [(s.D_T1 / (s.l_m * s.k)) if (s.l_m * s.k) > 0 else float("inf") for s in specs]
    i_ref = int(np.argmin(perms))
    return {"index": i_ref, "polymer": ctx["layers_geom"][i_ref].polymer,
            "permeability": float(perms[i_ref]), "k_per_layer": [float(s.k) for s in specs]}


def _g_per_substance(sv, ctx, sub):
    """Physical g for one substance: g1[n1] (1-step) or g12[n1,n2] (2-step) — exact at the prior,
    via the (physical key + grid) cache."""
    r = _cached_chain(sv, ctx, sub, ctx["time_vals"])
    return r[1] if ctx["two_step"] else r


def _c0eq_for_substance(sv, ctx, sub):
    """Genuine storage-step C0eq + default flag for one substance (read after its
    chain has run; see _C0EQ). Falls back to (nan, two_step) if absent."""
    key = _chain_key(ctx, _layer_specs(sv, ctx, sub), ctx["time_vals"])
    return _C0EQ.get(key, (float("nan"), bool(ctx["two_step"])))


def _diag_keys(ctx):
    subs = ctx["subs"]
    sw = np.array([float(getattr(s, "weight", 1.0)) for s in subs], dtype=float)
    return {
        "D_vals": np.array([float(s.D or 1e-14) for s in subs], dtype=float),
        "conc_vals": ctx["conc_vals"], "conc_weights": ctx["conc_w"],
        "time_vals": ctx["time_vals"], "time_weights": ctx["time_w"],
        "substance_weights": (sw / sw.sum() if sw.sum() > 0 else sw),
        "n_exchangeable": sum(1 for s in subs if getattr(s, "exchangeable", True)),
        "n_non_exchangeable": sum(1 for s in subs if not getattr(s, "exchangeable", True)),
    }


def physical_cf_tensor(sv) -> Dict[str, Any]:
    """Family CF tensor via the physical chain at the discretised prior (full key parity)."""
    ctx = _resolve(sv); subs = ctx["subs"]; conc = ctx["conc_vals"]; diag = _diag_keys(ctx)
    if ctx["two_step"]:
        surfaces = [_g_per_substance(sv, ctx, s) for s in subs]
        G = sv._combine_substance_surfaces(surfaces)
        CF = np.einsum("ij,k->ijk", G, conc)
        w = np.einsum("i,j,k->ijk", ctx["time_w"], ctx["time2_w"], ctx["conc_w"])
        w = (w / w.sum()).ravel()
        c0eq = np.array([_c0eq_for_substance(sv, ctx, s) for s in subs], float)
        return {"CF_tensor": CF, "CF_samples": CF.ravel(), "weights": w,
                "time2_vals": ctx["time2_vals"], "time2_weights": ctx["time2_w"],
                "C0eq": c0eq[:, 0], "C0eq_is_default": bool(c0eq[:, 1].any()), **diag}
    g_tj = np.empty((len(ctx["time_vals"]), len(subs)), float)
    for j, s in enumerate(subs):
        g_tj[:, j] = _g_per_substance(sv, ctx, s)
    G = sv._combine_substance_curves(g_tj)
    CF = np.outer(G, conc)
    w = np.outer(ctx["time_w"], ctx["conc_w"]).ravel(); w = w / w.sum() if w.sum() > 0 else w
    l_focal2 = float(ctx["layers_geom"][ctx["focal"]].thickness_m) ** 2
    Fo = np.outer(ctx["time_vals"], diag["D_vals"] / l_focal2) if l_focal2 > 0 else None
    c0eq = np.array([_c0eq_for_substance(sv, ctx, s) for s in subs], float)
    return {"CF_tensor": CF, "CF_samples": CF.ravel(), "weights": w, "Fo": Fo,
            "C0eq": c0eq[:, 0], "C0eq_is_default": bool(c0eq[:, 1].any()), **diag}


def physical_per_substance_cf_tensors(sv) -> Dict[str, Dict[str, np.ndarray]]:
    """Per-CAS decomposition (Σ_i CF_i == family) via the physical chain."""
    ctx = _resolve(sv); subs = ctx["subs"]; conc = ctx["conc_vals"]
    factor = sv._per_substance_factors()
    out: Dict[str, Dict[str, np.ndarray]] = {}

    def _store(cas, cf, extra):
        if cas in out:
            out[cas]["CF_tensor"] = out[cas]["CF_tensor"] + cf
            out[cas]["CF_samples"] = out[cas]["CF_tensor"].ravel()
        else:
            out[cas] = {"CF_tensor": cf, "CF_samples": cf.ravel(), **extra}

    if ctx["two_step"]:
        w3 = np.einsum("i,j,k->ijk", ctx["time_w"], ctx["time2_w"], ctx["conc_w"]); w3 /= w3.sum()
        # step2_active=True: the per-substance CF_tensor is 3-D (n_t1, n_t2, n_cp0);
        # without this flag the aggregator (combine_sources_shared_t) treats the
        # source as step-2-inactive and iterates t1 only, then raises "t2_idx
        # required for step-2-active source with 3-D CF_tensor". With it set, the
        # combiner marginalises over the oven-time (t2) prior — the same path the
        # family R5b aggregation uses — so per-substance 2-step groups aggregate.
        extra0 = {"weights": w3.ravel(), "time_vals": ctx["time_vals"], "time_weights": ctx["time_w"],
                  "time2_vals": ctx["time2_vals"], "time2_weights": ctx["time2_w"],
                  "conc_vals": conc, "conc_weights": ctx["conc_w"], "step2_active": True}
        for j, s in enumerate(subs):
            cf = np.einsum("ij,k->ijk", factor[j] * _g_per_substance(sv, ctx, s), conc)
            cas = getattr(s, "cas", None) or getattr(s, "id", f"idx{j}")
            _store(cas, cf, dict(extra0, exchangeable=bool(getattr(s, "exchangeable", True)),
                                 weight=float(getattr(s, "weight", 1.0))))
    else:
        w = np.outer(ctx["time_w"], ctx["conc_w"]).ravel()
        extra0 = {"weights": w, "time_vals": ctx["time_vals"], "time_weights": ctx["time_w"],
                  "conc_vals": conc, "conc_weights": ctx["conc_w"]}
        for j, s in enumerate(subs):
            cf = np.outer(factor[j] * _g_per_substance(sv, ctx, s), conc)
            cas = getattr(s, "cas", None) or getattr(s, "id", f"idx{j}")
            _store(cas, cf, dict(extra0, exchangeable=bool(getattr(s, "exchangeable", True)),
                                 weight=float(getattr(s, "weight", 1.0))))
    return out


def physical_step_resolved_cf_tensors(sv) -> Dict[str, np.ndarray]:
    """Step-resolved CF1 (storage alone), CF2 (oven alone), CF12 (full) via the physical chain."""
    ctx = _resolve(sv)
    if not ctx["two_step"]:
        raise ValueError("step_resolved requires a two-step scenario")
    subs = ctx["subs"]; conc = ctx["conc_vals"]
    t1v, t2v = ctx["time_vals"], ctx["time2_vals"]
    surf_full, g1_list, g2_list = [], [], []
    for s in subs:
        g1, g12 = _cached_chain(sv, ctx, s, t1v)                     # main chain (cached)
        surf_full.append(np.asarray(g12, float))
        g1_list.append(np.asarray(g1, float))                       # storage alone, per t1
        # oven alone (no storage): t1 -> ~0, then oven over t2 (distinct key → cached separately)
        _, g12_oven = _cached_chain(sv, ctx, s, np.array([_T2_FLOOR_S]))
        g2_list.append(np.asarray(g12_oven, float)[0])              # oven alone, per t2
    g12c = sv._combine_substance_surfaces(surf_full)
    g1c = sv._combine_substance_curves(np.asarray(g1_list).T)
    g2c = sv._combine_substance_curves(np.asarray(g2_list).T)
    w12 = np.einsum("i,j,k->ijk", ctx["time_w"], ctx["time2_w"], ctx["conc_w"]); w12 /= w12.sum()
    w1 = np.outer(ctx["time_w"], ctx["conc_w"]).ravel()
    w2 = np.outer(ctx["time2_w"], ctx["conc_w"]).ravel()
    return {"CF1": np.outer(g1c, conc), "CF2": np.outer(g2c, conc),
            "CF12": np.einsum("ij,k->ijk", g12c, conc),
            "w1": w1, "w2": w2, "w12": w12.ravel(),
            "t1_vals": t1v, "t2_vals": t2v, "conc_vals": conc}


def _wquantile(x, w, q):
    """Weighted quantile of flattened samples x with weights w (finite, normalised internally)."""
    x = np.asarray(x, float).ravel(); w = np.asarray(w, float).ravel()
    m = np.isfinite(x); x, w = x[m], w[m]
    if x.size == 0 or w.sum() <= 0:
        return float("nan")
    i = np.argsort(x); xs = x[i]; cw = np.cumsum(w[i]); cw /= cw[-1]
    return float(np.interp(q, cw, xs))


def worst_exchangeable_member(sv, qs=(0.5, 0.95, 0.99)):
    """Worst-exchangeable-member exposure — the conservative companion to the family AVERAGE (T2/C2).

    The family combine reports exchangeable additives as an occurrence-weighted AVERAGE (Σ p_i g_i): the
    EXPECTED exposure over which additive is present. At a quantile this UNDER-reports the worst plausible
    member (exchangeability = uncertainty on occurrence → exactly one is present, unknown which). This
    function reports, for each exchangeable additive, its CONDITIONAL-PRESENT exposure (g_i ⊗ conc, factor
    = 1, i.e. "if THIS additive is the one present"), takes its own weighted quantile, and returns the MAX
    across exchangeable members — "exposure if the worst plausible additive is present". Monomers
    (non-exchangeable) are excluded: they are genuinely co-present (summed), no occurrence ambiguity.

    Read-only; does NOT change the family combine. Returns the family quantile alongside for comparison
    (ratio worst/family is the T2 dilution factor). Physical-engine path (flag-gated).
    """
    ctx = _resolve(sv); subs = ctx["subs"]; conc = ctx["conc_vals"]
    exch = [(j, s) for j, s in enumerate(subs) if getattr(s, "exchangeable", True)]
    if not exch:
        return {"has_exchangeable": False, "n_exchangeable": 0,
                "note": "no exchangeable additives (all co-present monomers); worst-member not defined"}
    if ctx["two_step"]:
        w = np.einsum("i,j,k->ijk", ctx["time_w"], ctx["time2_w"], ctx["conc_w"]).ravel()
        def cf_cond(s):
            return np.einsum("ij,k->ijk", _g_per_substance(sv, ctx, s), conc).ravel()
    else:
        w = np.outer(ctx["time_w"], ctx["conc_w"]).ravel()
        def cf_cond(s):
            return np.outer(_g_per_substance(sv, ctx, s), conc).ravel()
    per_member = {}
    for j, s in exch:
        cf = cf_cond(s)
        cas = getattr(s, "cas", None) or getattr(s, "id", f"idx{j}")
        per_member[cas] = {"name": getattr(s, "name", None),
                           **{f"q{int(q*100)}": _wquantile(cf, w, q) for q in qs}}
    fam = physical_cf_tensor(sv)
    fam_q = {f"q{int(q*100)}": _wquantile(fam["CF_samples"], fam["weights"], q) for q in qs}
    qref = f"q{int(qs[1]*100)}" if len(qs) > 1 else f"q{int(qs[0]*100)}"     # rank by the middle q (q95)
    worst_cas = max(per_member, key=lambda c: (per_member[c][qref]
                                               if np.isfinite(per_member[c][qref]) else -np.inf))
    worst = {"cas": worst_cas, "name": per_member[worst_cas]["name"],
             **{k: per_member[worst_cas][k] for k in per_member[worst_cas] if k != "name"}}
    ratio = (worst[qref] / fam_q[qref]) if (np.isfinite(fam_q[qref]) and fam_q[qref] > 0) else float("nan")
    return {"has_exchangeable": True, "n_exchangeable": len(exch),
            "family_quantile": fam_q, "worst_member": worst, "per_member": per_member,
            "ratio_worst_over_family": ratio, "ranked_by": qref}
