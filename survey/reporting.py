"""
survey/reporting.py — Progress Reporting and Display
=====================================================

Provides progress bars and summary displays for survey computations.

@project: SFPPy/INSERM — Survey-scale exposure estimation
@author: Olivier Vitrac, PhD, HDR
@email: olivier.vitrac@gmail.com
@license: MIT
"""

import sys
import time
from typing import Optional


class ProgressReporter:
    """
    Text-based progress reporter with ETA and cache statistics.

    Provides friendly feedback during long computations.
    """

    def __init__(self, total_steps: int, desc: str = "Processing"):
        """
        Initialize progress reporter.

        Parameters
        ----------
        total_steps : int
            Total number of steps to complete.
        desc : str
            Description text for progress bar.
        """
        self.total = int(total_steps)
        self.desc = desc
        self.current = 0
        self.start_time = time.time()
        self.last_print = 0.0
        self.cache_hits = 0
        self.cache_misses = 0
        self._print_bar()

    def hit(self, n: int = 1) -> None:
        """Record cache hit(s)."""
        self.cache_hits += n
        self._maybe_print()

    def miss(self, n: int = 1) -> None:
        """Record cache miss(es)."""
        self.cache_misses += n
        self._maybe_print()

    def update(self, step: int = 1) -> None:
        """Update progress by step count."""
        self.current += int(step)
        self._maybe_print()

    def _maybe_print(self) -> None:
        """Print if enough time has elapsed."""
        now = time.time()
        if now - self.last_print > 0.15 or self.current >= self.total:
            self._print_bar()
            self.last_print = now

    def _print_bar(self) -> None:
        """Print progress bar to stdout."""
        elapsed = time.time() - self.start_time
        pct = min(1.0, self.current / self.total) if self.total > 0 else 0.0

        # Build progress bar
        bar_len = 30
        filled = int(bar_len * pct)
        bar = "\u2588" * filled + "\u2591" * (bar_len - filled)

        # Calculate ETA
        if elapsed > 1e-9 and self.current > 0:
            rate = self.current / elapsed
            eta = (self.total - self.current) / rate if rate > 0 else 0.0
        else:
            eta = 0.0

        eta_str = f"{int(eta // 60)}m {int(eta % 60):02d}s"

        msg = (
            f"\r{self.desc:<18} |{bar}| {pct*100:5.1f}% "
            f"[{self.current}/{self.total}] "
            f"ETA: {eta_str} | cache hit/miss: {self.cache_hits}/{self.cache_misses} "
        )
        sys.stdout.write(msg)
        sys.stdout.flush()

    def done(self) -> None:
        """Mark progress as complete."""
        self.current = self.total
        self._print_bar()
        sys.stdout.write("\n")
        sys.stdout.flush()


def format_survey_preview(
    name: str,
    n_substances: int,
    n_layers: int,
    i_ref: int,
    polymer: str,
    thickness_m: float,
    fo_range: tuple,
    cache_dir: str,
    cache_entries: int,
) -> str:
    """
    Format survey preview for display.

    Parameters
    ----------
    name : str
        Survey name.
    n_substances : int
        Number of substances.
    n_layers : int
        Number of layers.
    i_ref : int
        Reference layer index.
    polymer : str
        Polymer identifier.
    thickness_m : float
        Total thickness (m).
    fo_range : tuple
        (Fo_min, Fo_max) range.
    cache_dir : str
        Cache directory.
    cache_entries : int
        Number of existing cache entries.

    Returns
    -------
    str
        Formatted preview string.
    """
    lines = [
        f"Survey: {name}",
        "=" * 50,
        f"  Substances:      {n_substances}",
        f"  Layers:          {n_layers} ({'monolayer' if n_layers == 1 else 'multilayer'})",
    ]

    if n_layers > 1:
        lines.append(f"  Reference layer: {i_ref}")

    lines.extend([
        f"  Polymer:         {polymer}",
        f"  Thickness:       {thickness_m * 1e6:.1f} um",
        f"  Fo range:        [{fo_range[0]:.2e}, {fo_range[1]:.2e}]",
        "",
        f"  Cache dir:       {cache_dir}",
        f"  Cache entries:   {cache_entries}",
        "=" * 50,
    ])

    return "\n".join(lines)


def format_pdf_summary(
    pdf_bins: int,
    q50: float,
    q95: float,
    q99: float,
    mean: float,
) -> str:
    """
    Format PDF summary for display.

    Parameters
    ----------
    pdf_bins : int
        Number of PDF bins.
    q50, q95, q99 : float
        Quantiles.
    mean : float
        Mean value.

    Returns
    -------
    str
        Formatted summary string.
    """
    return (
        f"PDF Summary ({pdf_bins} bins):\n"
        f"  Mean:  {mean:.4f}\n"
        f"  q50:   {q50:.4f}\n"
        f"  q95:   {q95:.4f}\n"
        f"  q99:   {q99:.4f}"
    )
