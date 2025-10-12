# ── patankar.private.message.py ─────────────────────────────────────────────────
"""
Repetition-limited, shared, print-compatible messaging with history and stats.

This module provides a drop-in, `print`-compatible function that prevents
log flooding by limiting how many times identical messages are emitted while
still recording *all attempts* (printed or suppressed) in a shared history.
It is designed for multi-module scientific/engineering codebases.

All modules `import`-ing this file share the same global logger instance
(singleton). You can:
  • Set a global repetition limit and override it *per-caller* (module).
  • Maintain a rolling history with timestamps, caller, file:line, tags, etc.
  • Keep counters that increment even when messages are suppressed.
  • Compute stats grouped by message *key* or by *caller*.
  • Export the full history as a TSV file.
  • Temporarily mute output while still tracking (suppress/disable).
  • Reset counters after a configurable “silent” time window (global,
    per-caller, or per-call), so recurring messages can re-print after cooling.

Quick start
-----------
>>> import patankar.private.message as message
>>> message.initialize(maxrepetitions=2)
>>> message.print("Substance not found")      # prints
>>> message.print("Substance not found")      # prints
>>> message.print("Substance not found")      # suppressed (still counted & recorded)
>>> message.history(3)
2025-08-31 21:50:01 mymod.py:42 | Substance not found
2025-08-31 21:50:01 mymod.py:43 | Substance not found
2025-08-31 21:50:02 mymod.py:44 | Substance not found (su

>>> message.stats(by="key")
[('Substance not found', 2, 1, 3)]

>>> message.stats(by="caller")
[('mymodule.py', 2, 1, 3)]

>>> message.save("messages.tsv")  # write history to disk

Per-caller limits
-----------------
>>> message.set_limit(1)                        # global default
>>> message.set_limit(5, caller="toxicheck.py") # only for that module file

Silent delay reset (cooldown)
-----------------------------
If `tsilent` is set (globally, per-caller, or per-call), and more than
`tsilent` seconds have elapsed since a key was *last printed*, the per-key
counters are reset *before* the next attempt, allowing it to print again.

# prints twice, suppresses; after 10s of no prints for this *key*, the counter resets
>>> message.set_silence(10.0)               # global cooldown of 10s
>>> message.print("inventory miss X")       # prints
>>> message.print("inventory miss X")       # suppressed
# ... wait 10s without printing this key ...
>>> message.print("inventory miss X")       # prints again (counters reset for this key)

Silencing machinery (advanced)
-------------------
# Per-caller silent reset (module-specific):
>> message.set_silence_for("toxpredict.py", 5.0)   # only in toxpredict.py
# Per-call silent reset (strongest precedence):
>> message.print("inventory miss for CAS 123-45-6", tsilent=3.0)
Manually reset one key (e.g., after handling a burst):
>> message.reset_key("inventory miss for CAS 123-45-6")
# Optionally purge stale timestamps (keeps dict small in long runs):
>> message.purge_idle_keys(3600)  # forget keys idle for > 1 hour

Public API
----------
- initialize, print, clear, history, get_history
- stats, stats_for_caller, callers
- set_limit, reset_counts, set_history_size
- enable, disable, suppress, unsuppress
- save
- set_silence, set_silence_for, get_silence, reset_key, purge_idle_keys
- limit (context manager for temporary per-caller/global limits)

Notes
-----
- “Caller” is the file name of the immediate caller (e.g., `module.py`).
- History size is bounded with a deque (default 2000).
- Thread-safe via a re-entrant lock.

Typical use in SFPPy
---------------------
# in module A
import patankar.private.message as message
message.set_limit(2, caller="moduleA.py")  # only affects module A
message.print("tox alert: X missed threshold")

# in module B
import patankar.private.message as message
message.set_limit(10, caller="moduleB.py")
message.print("tox alert: X missed threshold")


Key inspection helpers
-----------------------
message.stats(by="key")         # messages sorted by total (printed+suppressed)
message.stats(by="caller")      # callers sorted by total
message.stats_for_caller("moduleA.py", top=20)
message.callers()               # [(caller, printed, suppressed, total), ...]
message.history(20)             # last 20, all callers
message.history(filter_caller="moduleA.py")


---
@version: 1.50
@project: SFPPy - SafeFoodPackaging Portal in Python initiative
@author: INRAE\\olivier.vitrac@agroparistech.fr
@licence: MIT
@Date: 2024-02-10
@rev: 2025-08-31
"""

from __future__ import annotations
import sys
import time
import threading
import traceback
import builtins
from collections import defaultdict, deque

__all__ = ['callers', 'clear', 'disable', 'enable', 'get_history', 'get_silence', 'history', 'initialize', 'limit', 'print', 'purge_idle_keys', 'reset_counts', 'reset_key', 'save', 'set_history_size', 'set_limit', 'set_silence', 'set_silence_for', 'stats', 'stats_for_caller', 'suppress', 'unsuppress']

__project__ = "SFPPy"
__author__ = "Olivier Vitrac"
__copyright__ = "Copyright 2022"
__credits__ = ["Olivier Vitrac"]
__license__ = "MIT"
__maintainer__ = "Olivier Vitrac"
__email__ = "olivier.vitrac@agroparistech.fr"
__version__ = "1.50"


class _MessageLogger:
    """
    Core singleton implementing repetition-limited printing, history, and stats.

    This class is used through the module-level façade functions (e.g.,
    `message.print`, `message.stats`) that delegate to a single global instance.

    Data recorded per attempt (history record dict):
        "ts"            float UNIX timestamp of the attempt
        "text"          message text (tagged if tag provided)
        "key"           canonical message key (possibly transformed by `key=`)
        "src"           "filename:lineno" of the call site
        "caller"        caller identifier (filename only)
        "tag"           optional short category string
        "suppressed"    bool, True if not emitted to stream
        "count_at_emit" int, count of this key *after* increment at this attempt
        "stream"        "stdout", "stderr", or "file"

    Counters:
        • Per-key totals (attempted) and suppressed counts.
        • Per-caller totals (attempted) and suppressed counts.
      Counters increment even when messages are suppressed.

    Limits:
        • Global repetition limit for all keys.
        • Per-caller limit (overrides global for that caller).
        • Per-call `once=True` forces a limit of 1 for that specific call.

    Silent-delay reset (cooldown):
        • Global `tsilent` seconds, optional per-caller override, and a
          per-call `tsilent=...` that takes precedence.
        • If more than `tsilent` seconds have elapsed since a key was last
          *printed*, the counters for that key are reset just before counting
          the new attempt. This allows it to print again after a quiet period.
    """

    def __init__(self):
        """Initialize all shared state with sensible defaults."""
        self._lock = threading.RLock()
        self._counts = defaultdict(int)               # key -> total attempts
        self._suppressed_counts = defaultdict(int)    # key -> suppressed attempts
        self._counts_by_caller = defaultdict(int)     # caller -> total attempts
        self._suppressed_by_caller = defaultdict(int) # caller -> suppressed attempts

        self._hist = deque(maxlen=2000)               # history of records (dict)
        self._hist_size = 2000

        self._maxrepetitions = 1                      # global per-key limit
        self._per_caller_limit = {}                   # caller -> limit

        self._enabled = True                          # if False, track but never emit
        self._globally_suppressed = False             # temporary mute (track anyway)

        # Silent-delay (cooldown) machinery:
        self._tsilent_default = None                  # None => disabled
        self._tsilent_by_caller = {}                  # caller -> seconds
        self._last_print_ts = {}                      # key -> last *printed* timestamp

    # ── configuration ────────────────────────────────────────────────────────
    def initialize(self, maxrepetitions=1, history_size=2000, enabled=True):
        """
        Initialize or reconfigure the logger.

        Parameters
        ----------
        maxrepetitions : int
            Default maximum emitted times per unique key. Attempts beyond
            this limit are suppressed but still tracked.
        history_size : int
            Bounded history length (oldest entries drop when exceeded).
        enabled : bool
            If False, suppress actual emissions (history and counters still recorded).
        """
        with self._lock:
            self._maxrepetitions = max(0, int(maxrepetitions))
            self._hist_size = max(0, int(history_size))
            self._enabled = bool(enabled)
            # keep existing content but resize deque
            self._hist = deque(self._hist, maxlen=self._hist_size)

    def set_limit(self, n, caller=None):
        """
        Set repetition limit.

        Parameters
        ----------
        n : int
            Maximum emitted times per key.
        caller : str or None
            If None, set global default. If a string (e.g., "module.py"),
            set/override the limit for that caller only.
        """
        with self._lock:
            n = max(0, int(n))
            if caller is None:
                self._maxrepetitions = n
            else:
                self._per_caller_limit[caller] = n

    def set_history_size(self, n):
        """
        Change the maximum number of records retained in history.

        Parameters
        ----------
        n : int
            New bound for the history deque.
        """
        with self._lock:
            n = max(0, int(n))
            old = list(self._hist)
            self._hist_size = n
            self._hist = deque(old[-n:], maxlen=n)

    def enable(self):
        """
        Enable actual emissions to the chosen stream (default stdout).

        Tracking and history are *always* on; this only controls emission.
        """
        with self._lock:
            self._enabled = True

    def disable(self):
        """
        Disable actual emissions (no printing), but keep tracking and history.
        """
        with self._lock:
            self._enabled = False

    def suppress(self):
        """
        Temporarily mute emissions (like a global 'quiet' mode).
        Tracking and history continue.
        """
        with self._lock:
            self._globally_suppressed = True

    def unsuppress(self):
        """Re-enable emissions after `suppress()`."""
        with self._lock:
            self._globally_suppressed = False

    # ── silent delay (cooldown) configuration ────────────────────────────────
    def set_silence(self, seconds):
        """
        Set a global silent delay (cooldown) in seconds.

        If not None, and if a key was last *printed* more than `seconds` ago,
        its per-key counters are reset just before the next attempt, allowing
        it to print again.

        Parameters
        ----------
        seconds : float or None
            Global cooldown. Use None to disable globally.
        """
        with self._lock:
            self._tsilent_default = None if seconds is None else float(seconds)

    def set_silence_for(self, caller, seconds):
        """
        Set a per-caller silent delay (cooldown), overriding the global value.

        Parameters
        ----------
        caller : str
            Caller identifier (usually a filename, e.g., "module.py").
        seconds : float or None
            Cooldown for this caller. Use None to remove the override.
        """
        with self._lock:
            if seconds is None:
                self._tsilent_by_caller.pop(caller, None)
            else:
                self._tsilent_by_caller[caller] = float(seconds)

    def get_silence(self, caller=None):
        """
        Get the effective silent delay for a caller (or the global default).

        Parameters
        ----------
        caller : str or None
            If provided and present in per-caller overrides, return it;
            otherwise return the global default.

        Returns
        -------
        float or None
            Cooldown in seconds, or None if disabled.
        """
        with self._lock:
            if caller is not None and caller in self._tsilent_by_caller:
                return self._tsilent_by_caller[caller]
            return self._tsilent_default

    def reset_key(self, key):
        """
        Reset counters and last-printed timestamp for a specific message key.

        Parameters
        ----------
        key : str
            The canonical key used for counting (after optional key-normalization).
        """
        with self._lock:
            self._counts.pop(key, None)
            self._suppressed_counts.pop(key, None)
            self._last_print_ts.pop(key, None)

    def purge_idle_keys(self, idle_seconds):
        """
        Remove last-printed timestamps for keys idle longer than `idle_seconds`.

        Does not affect history; only forgets _last_print_ts entries to keep
        memory small in very long runs.

        Parameters
        ----------
        idle_seconds : float
            Entries older than now - idle_seconds are purged.
        """
        with self._lock:
            tnow = time.time()
            stale = [k for k, ts in self._last_print_ts.items()
                     if (tnow - ts) > float(idle_seconds)]
            for k in stale:
                self._last_print_ts.pop(k, None)

    # ── utils ────────────────────────────────────────────────────────────────
    def _caller_id(self, explicit=None):
        """
        Determine the caller identifier (module filename).

        Parameters
        ----------
        explicit : str or None
            If provided, use as caller; else infer from stack.

        Returns
        -------
        str
            Caller identifier (filename).
        """
        if explicit:
            return explicit
        # stack: [ ... this fn, print(), external caller ]
        frm = traceback.extract_stack(limit=4)[0]
        filename = frm.filename.rsplit("/", 1)[-1]
        return filename

    # ── core print ───────────────────────────────────────────────────────────
    def print(self, *args, **kwargs):
        """
        Print with repetition limiting, history, and optional cooldown.

        Signature-compatible with built-in print:
            print(*args, sep=' ', end='\\n', file=None, flush=False)

        Extra keyword-only arguments:
        --------------------------------
        key : callable or None
            A function str -> str to normalize the message for de-duplication
            (e.g., strip numbers or IDs). If None, the exact text is used.
        tag : str or None
            Optional short category prefix (visible in history).
        once : bool
            If True, enforces a one-time emission limit for this call.
        caller : str or None
            Override/force the caller identifier for this emission.
        tsilent : float or None
            Per-call cooldown (seconds). If set, it takes precedence over any
            per-caller or global cooldown for this call only.
        """
        sep = kwargs.pop("sep", " ")
        end = kwargs.pop("end", "\n")
        file = kwargs.pop("file", None)
        flush = kwargs.pop("flush", False)
        key_fn = kwargs.pop("key", None)
        tag = kwargs.pop("tag", None)
        once = kwargs.pop("once", False)
        caller_override = kwargs.pop("caller", None)
        tsilent_call = kwargs.pop("tsilent", None)
        if kwargs:
            raise TypeError("Unexpected kwargs: %r" % list(kwargs.keys()))

        text = sep.join(map(str, args))
        k = key_fn(text) if key_fn else text

        # source with line number for history; per-caller limit uses filename only
        frm = traceback.extract_stack(limit=3)[0]
        src = "%s:%d" % (frm.filename.rsplit("/", 1)[-1], frm.lineno)
        caller = self._caller_id(caller_override)

        with self._lock:
            # Resolve effective cooldown (per-call > per-caller > global)
            if tsilent_call is not None:
                eff_tsilent = float(tsilent_call)
            else:
                eff_tsilent = self._tsilent_by_caller.get(caller, self._tsilent_default)

            tnow = time.time()

            # If cooldown elapsed since last *printed* emission of this key, reset counts
            if eff_tsilent is not None:
                last_ts = self._last_print_ts.get(k)
                if (last_ts is not None) and ((tnow - last_ts) > eff_tsilent):
                    self._counts.pop(k, None)
                    self._suppressed_counts.pop(k, None)

            limit_n = 1 if once else self._per_caller_limit.get(caller, self._maxrepetitions)

            # Increment counters (attempts always counted)
            self._counts[k] += 1
            self._counts_by_caller[caller] += 1
            n_seen = self._counts[k]

            suppressed = (n_seen > limit_n) or self._globally_suppressed or (not self._enabled)

            rec = {
                "ts": tnow,
                "text": text if tag is None else "[%s] %s" % (tag, text),
                "key": k,
                "src": src,
                "caller": caller,
                "tag": tag,
                "suppressed": suppressed,
                "count_at_emit": n_seen,
                "stream": "stderr" if (file is sys.stderr) else (
                    "stdout" if (file in (None, sys.stdout)) else "file")
            }
            self._hist.append(rec)

            if suppressed:
                self._suppressed_counts[k] += 1
                self._suppressed_by_caller[caller] += 1
            else:
                target = sys.stdout if file is None else file
                builtins.print(text, sep="", end=end, file=target, flush=flush)
                # record last printed time for cooldown checks
                self._last_print_ts[k] = tnow

    # ── history & stats ──────────────────────────────────────────────────────
    def clear(self):
        """
        Clear all internal state:
          - history, all counters, per-caller limits, and last-printed times.
        """
        with self._lock:
            self._counts.clear()
            self._suppressed_counts.clear()
            self._counts_by_caller.clear()
            self._suppressed_by_caller.clear()
            self._hist.clear()
            self._per_caller_limit.clear()
            self._last_print_ts.clear()

    def reset_counts(self):
        """
        Reset all counters (printed/suppressed by key and by caller).

        History is preserved; only counts are cleared.
        """
        with self._lock:
            self._counts.clear()
            self._suppressed_counts.clear()
            self._counts_by_caller.clear()
            self._suppressed_by_caller.clear()

    def get_history(self, n=None, include_suppressed=True, filter_caller=None):
        """
        Return recorded history as a list of dict records.

        Parameters
        ----------
        n : int or None
            If None, return full history; else return only the last n records.
        include_suppressed : bool
            If False, filter out suppressed entries.
        filter_caller : str or None
            If provided, return only records for this caller.

        Returns
        -------
        list of dict
            Each dict contains the fields described in the class docstring.
        """
        with self._lock:
            data = list(self._hist) if n is None else list(self._hist)[-int(n):]
            if not include_suppressed:
                data = [r for r in data if not r["suppressed"]]
            if filter_caller:
                data = [r for r in data if r["caller"] == filter_caller]
            return data

    def history(self, n=None, include_suppressed=True, timestamps=True, filter_caller=None):
        """
        Pretty-print history to stdout.

        Parameters
        ----------
        n : int or None
            If None, print full history; else last n entries.
        include_suppressed : bool
            If False, skip suppressed entries.
        timestamps : bool
            If True, prefix each line with a human-readable timestamp.
        filter_caller : str or None
            If provided, only show entries for this caller.
        """
        recs = self.get_history(n=n, include_suppressed=include_suppressed,
                                filter_caller=filter_caller)
        for r in recs:
            ts = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(r["ts"])) if timestamps else ""
            prefix = ("%s " % ts) if timestamps else ""
            sup = " (suppressed)" if r["suppressed"] else ""
            sys.stdout.write("%s%s | %s%s\n" % (prefix, r["src"], r["text"], sup))

    def callers(self):
        """
        Return aggregated counts per caller.

        Returns
        -------
        list of tuples
            (caller, printed, suppressed, total), sorted by total desc.
        """
        with self._lock:
            out = []
            for c, tot in self._counts_by_caller.items():
                sup = self._suppressed_by_caller.get(c, 0)
                printed = tot - sup
                out.append((c, printed, sup, tot))
            out.sort(key=lambda t: t[3], reverse=True)
            return out

    def stats(self, top=None, min_count=1, by="key"):
        """
        Return aggregated counts sorted by total.

        Parameters
        ----------
        top : int or None
            If provided, return only the first `top` entries.
        min_count : int
            Minimum total (printed+suppressed) required to appear.
        by : {"key", "caller"}
            Group by message key or by caller.

        Returns
        -------
        list of tuples
            If by="key":    (key, printed, suppressed, total)
            If by="caller": (caller, printed, suppressed, total)
        """
        with self._lock:
            out = []
            if by == "key":
                for k, tot in self._counts.items():
                    sup = self._suppressed_counts.get(k, 0)
                    printed = tot - sup
                    if (tot >= min_count) or (sup >= min_count):
                        out.append((k, printed, sup, tot))
            elif by == "caller":
                for c, tot in self._counts_by_caller.items():
                    sup = self._suppressed_by_caller.get(c, 0)
                    printed = tot - sup
                    if (tot >= min_count) or (sup >= min_count):
                        out.append((c, printed, sup, tot))
            else:
                raise ValueError('by must be "key" or "caller"')
            out.sort(key=lambda t: t[3], reverse=True)
            return out[:top] if top else out

    def stats_for_caller(self, caller, top=None, min_count=1):
        """
        Return per-key stats restricted to a given caller.

        Parameters
        ----------
        caller : str
            Caller identifier (filename).
        top : int or None
            If provided, only the first `top` entries are returned.
        min_count : int
            Minimum total attempts to include.

        Returns
        -------
        list of tuples
            (key, printed, suppressed, total), sorted by total desc.
        """
        with self._lock:
            agg = defaultdict(lambda: [0, 0])  # key -> [tot, sup]
            for r in self._hist:
                if r["caller"] != caller:
                    continue
                agg[r["key"]][0] += 1
                if r["suppressed"]:
                    agg[r["key"]][1] += 1
            out = []
            for k, pair in agg.items():
                tot, sup = pair
                printed = tot - sup
                if (tot >= min_count) or (sup >= min_count):
                    out.append((k, printed, sup, tot))
            out.sort(key=lambda t: t[3], reverse=True)
            return out[:top] if top else out

    def save(self, path, include_suppressed=True, filter_caller=None):
        """
        Save history to a TSV file on disk.

        Parameters
        ----------
        path : str
            Output file path.
        include_suppressed : bool
            If False, skip suppressed entries.
        filter_caller : str or None
            If provided, export only records for this caller.
        """
        recs = self.get_history(include_suppressed=include_suppressed,
                                filter_caller=filter_caller)
        with open(path, "w", encoding="utf-8") as f:
            f.write("ts\tstream\tsrc\tcaller\ttag\tcount\tSuppressed\ttext\n")
            for r in recs:
                f.write("%d\t%s\t%s\t%s\t%s\t%d\t%d\t%s\n" % (
                    int(r["ts"]), r["stream"], r["src"], r["caller"], r["tag"] or "",
                    r["count_at_emit"], 1 if r["suppressed"] else 0, r["text"]))

    # ── context manager for temporary limits ─────────────────────────────────
    def limit(self, n, caller=None):
        """
        Context manager to temporarily override repetition limits.

        Parameters
        ----------
        n : int
            Temporary maximum emitted times per key.
        caller : str or None
            If None, override the global default; else override only for
            the specified caller (filename). Previous values are restored
            on exit.

        Examples
        --------
        >>> with message.limit(0):             # globally mute emissions (still tracked)
        ...     message.print("no output now")
        >>> with message.limit(10, caller="tox.py"):
        ...     message.print("verbose for tox only")
        """
        logger = self
        class _LimitCtx:
            def __enter__(self_):
                logger._lock.acquire()
                self_.caller = caller or None
                if self_.caller is None:
                    self_.prev = logger._maxrepetitions
                    logger._maxrepetitions = max(0, int(n))
                else:
                    self_.prev = logger._per_caller_limit.get(self_.caller, None)
                    logger._per_caller_limit[self_.caller] = max(0, int(n))
                return logger
            def __exit__(self_, exc_type, exc, tb):
                if self_.caller is None:
                    logger._maxrepetitions = self_.prev
                else:
                    if self_.prev is None:
                        logger._per_caller_limit.pop(self_.caller, None)
                    else:
                        logger._per_caller_limit[self_.caller] = self_.prev
                logger._lock.release()
        return _LimitCtx()


# ── module-level façade (singleton) ──────────────────────────────────────────
_LOGGER = _MessageLogger()


def initialize(maxrepetitions=1, history_size=2000, enabled=True):
    """
    Initialize/reconfigure the shared logger. See `_MessageLogger.initialize`.
    """
    _LOGGER.initialize(maxrepetitions=maxrepetitions,
                       history_size=history_size,
                       enabled=enabled)


def print(*args, **kwargs):
    """
    Print with repetition limiting, history, and optional cooldown.

    Same signature as built-in `print` + extras: key, tag, once, caller, tsilent.
    See `_MessageLogger.print` for details.
    """
    _LOGGER.print(*args, **kwargs)


def clear():
    """
    Clear all state: history, counters, per-caller limits, last-printed times.
    """
    _LOGGER.clear()


def reset_counts():
    """
    Reset per-key and per-caller counters (history preserved).
    """
    _LOGGER.reset_counts()


def set_limit(n, caller=None):
    """
    Set repetition limit globally or for a specific caller.
    """
    _LOGGER.set_limit(n, caller=caller)


def set_history_size(n):
    """
    Change the maximum number of records retained in history.
    """
    _LOGGER.set_history_size(n)


def enable():
    """
    Enable actual emissions (tracking is always on).
    """
    _LOGGER.enable()


def disable():
    """
    Disable actual emissions (tracking continues).
    """
    _LOGGER.disable()


def suppress():
    """
    Temporarily mute emissions (tracking continues).
    """
    _LOGGER.suppress()


def unsuppress():
    """
    Re-enable emissions after `suppress()`.
    """
    _LOGGER.unsuppress()


def history(n=None, include_suppressed=True, timestamps=True, filter_caller=None):
    """
    Pretty-print history. See `_MessageLogger.history` for parameters.
    """
    _LOGGER.history(n=n, include_suppressed=include_suppressed,
                    timestamps=timestamps, filter_caller=filter_caller)


def get_history(n=None, include_suppressed=True, filter_caller=None):
    """
    Return history records as a list of dicts. See `_MessageLogger.get_history`.
    """
    return _LOGGER.get_history(n=n, include_suppressed=include_suppressed,
                               filter_caller=filter_caller)


def stats(top=None, min_count=1, by="key"):
    """
    Aggregated counts sorted by total. See `_MessageLogger.stats`.
    """
    return _LOGGER.stats(top=top, min_count=min_count, by=by)


def stats_for_caller(caller, top=None, min_count=1):
    """
    Per-key stats restricted to a given caller. See `_MessageLogger.stats_for_caller`.
    """
    return _LOGGER.stats_for_caller(caller, top=top, min_count=min_count)


def callers():
    """
    Aggregated counts per caller. See `_MessageLogger.callers`.
    """
    return _LOGGER.callers()


def save(path, include_suppressed=True, filter_caller=None):
    """
    Save history to a TSV file on disk. See `_MessageLogger.save`.
    """
    _LOGGER.save(path=path, include_suppressed=include_suppressed,
                 filter_caller=filter_caller)


def set_silence(seconds):
    """
    Set global silent delay (cooldown) in seconds, or None to disable.
    """
    _LOGGER.set_silence(seconds)


def set_silence_for(caller, seconds):
    """
    Set per-caller silent delay (cooldown), or None to remove override.
    """
    _LOGGER.set_silence_for(caller, seconds)


def get_silence(caller=None):
    """
    Get the effective cooldown for a caller or the global default.
    """
    return _LOGGER.get_silence(caller=caller)


def reset_key(key):
    """
    Reset counters and last-printed timestamp for a specific message key.
    """
    _LOGGER.reset_key(key)


def purge_idle_keys(idle_seconds):
    """
    Forget last-printed timestamps for keys idle longer than `idle_seconds`.
    """
    _LOGGER.purge_idle_keys(idle_seconds)


def limit(n, caller=None):
    """
    Context manager to temporarily override repetition limits (global or per-caller).
    """
    return _LOGGER.limit(n, caller=caller)
# ─────────────────────────────────────────────────────────────────────────────
