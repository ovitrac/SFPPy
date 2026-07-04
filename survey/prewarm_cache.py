#!/usr/bin/env python3
"""survey.prewarm_cache — robust, reusable OFFLINE-FIRST substance-cache prewarmer.

Per-scenario substance resolution makes a **live PubChem PUG-REST** call (~1-3 s cold,
and its ServerBusy/timeout retry loop balloons when PubChem throttles a large run) and
spawns a **local ToxTree Java** subprocess (~1.6 s cold, JVM-per-call, can hang). Across
a ~10^5-scenario survey those costs and throttling stall the worker pool.

This module resolves every *unique* substance ONCE via :class:`patankar.loadpubchem.migrantToxtree`,
populating both ``cache.PubChem`` and ``cache.ToxTree`` so the subsequent run is fully
local — no live PubChem, no per-scenario JVM startup, no throttling, no hangs.

Robust + reusable:
  • RESUMABLE — each CAS outcome is journalled to a JSON state file; a re-run skips
    terminally-resolved CAS (ok / not_found) and RETRIES only transient ones
    (timeout / server-busy / network). Safe to interrupt and relaunch.
  • PER-CAS TIMEOUT — SIGALRM so a hung PubChem retry or ToxTree JVM never stalls the pass.
  • RETRY w/ BACKOFF for transient failures.
  • SURROGATE-AWARE — applies :data:`survey.io.SURROGATE_CAS` so the production-relevant
    target (the surrogate) is what gets cached.
  • Programmatic API (:func:`prewarm`) + CLI; CAS from a dbs sheet, scenario YAMLs, or a list.

CLI:
  python -m survey.prewarm_cache --xlsx path/to/dbs.xlsx --state dbs.prewarm.json --timeout 30
  python -m survey.prewarm_cache --scenarios-dir path/to/scenarios
"""
from __future__ import annotations

import argparse
import json
import re
import signal
import sys
import time
from pathlib import Path
from typing import Dict, Iterable, List, Optional

from patankar.loadpubchem import migrantToxtree

try:
    from survey.io import SURROGATE_CAS
except Exception:  # pragma: no cover - survey.io always present in the package
    SURROGATE_CAS = {}

__all__ = ["prewarm", "collect_cas", "expand_surrogates",
           "cas_from_xlsx", "cas_from_scenarios", "cas_from_file"]

CAS_RE = re.compile(r"\b(\d{2,7}-\d{2}-\d)\b")
_TRANSIENT = ("timeout", "serverbusy", "server busy", "temporarily",
              "connection", "network", "rate limit", "ratelimit", "throttl")
_TERMINAL = ("ok", "not_found")


class _Timeout(Exception):
    pass


# --------------------------------------------------------------------------- inputs
def cas_from_xlsx(xlsx: str, sheet: str = "FCC") -> set:
    import pandas as pd
    df = pd.read_excel(xlsx, sheet_name=sheet)
    df.columns = [str(c).strip() for c in df.columns]
    col = next((c for c in df.columns if "cas" in c.lower()), None)
    if col is None:
        raise ValueError(f"no CAS column in {xlsx}[{sheet}]")
    return {str(x).strip() for x in df[col].dropna()}


def cas_from_scenarios(scenarios_dir: str) -> set:
    out = set()
    for p in Path(scenarios_dir).glob("*.yml"):
        out |= set(CAS_RE.findall(p.read_text()))
    return out


def cas_from_file(path: str) -> set:
    return {ln.strip() for ln in Path(path).read_text().splitlines() if ln.strip()}


def expand_surrogates(cas: Iterable[str]) -> List[str]:
    """Map each well-formed CAS to what production will actually resolve: a surrogated
    CAS is REPLACED by its target(s) (the ladder always redirects it, so the raw CAS is
    never resolved at run time and need not be cached); all others pass through."""
    targets = set()
    for c in cas:
        if not CAS_RE.fullmatch(str(c)):
            continue
        s = SURROGATE_CAS.get(c)
        if s:
            targets |= set(s) if isinstance(s, (tuple, list)) else {s}
        else:
            targets.add(c)
    return sorted(targets)


def collect_cas(xlsx: Optional[str] = None, sheet: str = "FCC",
                scenarios_dir: Optional[str] = None, cas_file: Optional[str] = None) -> List[str]:
    cas: set = set()
    if xlsx:
        cas |= cas_from_xlsx(xlsx, sheet)
    if scenarios_dir:
        cas |= cas_from_scenarios(scenarios_dir)
    if cas_file:
        cas |= cas_from_file(cas_file)
    return expand_surrogates(cas)


# --------------------------------------------------------------------------- resolve
def resolve_one(cas: str, timeout: int = 30) -> tuple:
    """Resolve one CAS via migrantToxtree with an optional hard SIGALRM timeout.

    Returns ``(status, detail)`` with status ∈ {ok, not_found, timeout, transient, error}.
    SIGALRM is only armed on the main thread (the normal prewarm context).
    """
    armed = False
    if timeout > 0:
        try:
            def _alarm(_s, _f):
                raise _Timeout()
            old = signal.signal(signal.SIGALRM, _alarm)
            signal.setitimer(signal.ITIMER_REAL, timeout)
            armed = True
        except ValueError:
            armed = False  # not main thread → no hard timeout (migrant has its own retry cap)
    try:
        m = migrantToxtree(cas, verbose=False, raiseerror=False)
        cid = getattr(m, "cid", None)
        return ("ok", str(cid)) if cid else ("not_found", "")
    except _Timeout:
        return ("timeout", f">{timeout}s")
    except Exception as e:
        low = str(e).lower()
        # not-found is terminal (don't retry). Includes the loadpubchem not-found
        # signatures: a clean "not found"/"does not exist", an empty-object return,
        # and the `all_data` UnboundLocalError raised on a missing compound.
        if any(k in low for k in ("not found", "does not exist", "empty object", "all_data")):
            return ("not_found", str(e)[:160])
        if any(t in low for t in _TRANSIENT):
            return ("transient", str(e)[:160])
        return ("error", str(e)[:160])
    finally:
        if armed:
            signal.setitimer(signal.ITIMER_REAL, 0)
            signal.signal(signal.SIGALRM, old)


# --------------------------------------------------------------------------- driver
def prewarm(cas: Iterable[str], state_path: Optional[str] = None, *,
            timeout: int = 30, retries: int = 2, force: bool = False,
            verbose: bool = True) -> Dict[str, dict]:
    """Resolve each CAS once (surrogate targets assumed already expanded), journalling
    outcomes to ``state_path`` for resumability. Returns the full journal dict."""
    cas = list(dict.fromkeys(cas))  # dedupe, keep order
    state: Dict[str, dict] = {}
    sp = Path(state_path) if state_path else None
    if sp and sp.exists() and not force:
        try:
            state = json.loads(sp.read_text())
        except Exception:
            state = {}

    todo = [c for c in cas if force or state.get(c, {}).get("status") not in _TERMINAL]
    if verbose:
        print(f"prewarm: {len(cas)} targets; {len(cas)-len(todo)} done, {len(todo)} to resolve"
              + (f"; journal={sp}" if sp else ""), flush=True)

    counts: Dict[str, int] = {}
    t0 = time.time()
    for i, c in enumerate(todo, 1):
        status = detail = None
        attempt = 0
        for attempt in range(1, retries + 1):
            status, detail = resolve_one(c, timeout)
            if status in ("ok", "not_found", "error"):
                break
            time.sleep(min(5.0, 0.5 * 2 ** attempt))  # backoff for transient/timeout
        state[c] = {"status": status, "detail": detail, "attempts": attempt}
        counts[status] = counts.get(status, 0) + 1
        if sp and (i % 25 == 0 or i == len(todo)):
            sp.write_text(json.dumps(state, indent=0))  # checkpoint
        if verbose and (i % 25 == 0 or i == len(todo)):
            el = time.time() - t0
            print(f"  {i}/{len(todo)}  {dict(sorted(counts.items()))}  ({el/max(i,1):.2f} s/cas)",
                  flush=True)

    if sp:
        sp.write_text(json.dumps(state, indent=0))
    if verbose:
        tot: Dict[str, int] = {}
        for v in state.values():
            tot[v["status"]] = tot.get(v["status"], 0) + 1
        n_retry = sum(1 for v in state.values() if v["status"] in ("timeout", "transient"))
        print(f"done in {time.time()-t0:.0f}s. totals: {dict(sorted(tot.items()))}"
              + (f" — {n_retry} transient/timeout left (re-run to retry)" if n_retry else ""),
              flush=True)
    return state


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--xlsx", help="dbs workbook (CAS from --sheet)")
    ap.add_argument("--sheet", default="FCC")
    ap.add_argument("--scenarios-dir", help="scrape CAS from scenario YAMLs")
    ap.add_argument("--cas-file", help="text file, one CAS per line")
    ap.add_argument("--state", default=None, help="resume journal (JSON)")
    ap.add_argument("--timeout", type=int, default=30, help="per-CAS hard timeout s (0=off)")
    ap.add_argument("--retries", type=int, default=2)
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args(argv)
    if not (a.xlsx or a.scenarios_dir or a.cas_file):
        ap.error("provide at least one of --xlsx / --scenarios-dir / --cas-file")
    cas = collect_cas(a.xlsx, a.sheet, a.scenarios_dir, a.cas_file)
    state = a.state or ((a.xlsx or a.scenarios_dir or a.cas_file) + ".prewarm.json")
    prewarm(cas, state, timeout=a.timeout, retries=a.retries, force=a.force)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
