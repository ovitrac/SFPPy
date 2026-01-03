"""
survey/visualization.py — Visualization Functions for Survey Results
====================================================================

Generates PDF, PNG, and SVG plots of survey PDF/CDF distributions.

@project: SFPPy/INSERM — Survey-scale exposure estimation
@author: Olivier Vitrac, PhD, HDR
@email: olivier.vitrac@gmail.com
@license: MIT
"""

from pathlib import Path
from typing import Optional, Tuple, Dict, Any, List
import numpy as np

try:
    import matplotlib
    matplotlib.use('Agg')  # Non-interactive backend
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


def check_matplotlib():
    """Check if matplotlib is available."""
    if not HAS_MATPLOTLIB:
        raise ImportError(
            "matplotlib is required for visualization. "
            "Install with: pip install matplotlib"
        )


def plot_pdf_cdf(
    centers: np.ndarray,
    pdf: np.ndarray,
    cdf: np.ndarray,
    title: str = "Migration Survey Results",
    xlabel: str = "Concentration CF (mg/kg food)",
    quantiles: Optional[Dict[float, float]] = None,
    figsize: Tuple[int, int] = (12, 5),
) -> "matplotlib.figure.Figure":
    """
    Create dual-panel plot with PDF and CDF.

    Parameters
    ----------
    centers : np.ndarray
        PDF bin centers.
    pdf : np.ndarray
        PDF values.
    cdf : np.ndarray
        CDF values.
    title : str
        Plot title.
    xlabel : str
        X-axis label.
    quantiles : dict, optional
        Quantile values to mark: {0.50: value, 0.95: value, ...}
    figsize : tuple
        Figure size (width, height).

    Returns
    -------
    matplotlib.figure.Figure
        Created figure.
    """
    check_matplotlib()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)

    # PDF plot
    ax1.fill_between(centers, pdf, alpha=0.3, color='steelblue', label='PDF')
    ax1.plot(centers, pdf, color='steelblue', linewidth=1.5)
    ax1.set_xlabel(xlabel, fontsize=11)
    ax1.set_ylabel('Probability Density', fontsize=11)
    ax1.set_title('Probability Density Function (PDF)', fontsize=12)
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(left=0)

    # Add quantile markers to PDF
    if quantiles:
        colors = {0.50: 'green', 0.95: 'orange', 0.99: 'red'}
        for q, val in quantiles.items():
            color = colors.get(q, 'gray')
            ax1.axvline(val, color=color, linestyle='--', linewidth=1.5, alpha=0.7,
                       label=f'q{int(q*100)} = {val:.3f}')
        ax1.legend(loc='upper right', fontsize=9)

    # CDF plot
    ax2.plot(centers, cdf, color='darkblue', linewidth=2, label='CDF')
    ax2.set_xlabel(xlabel, fontsize=11)
    ax2.set_ylabel('Cumulative Probability', fontsize=11)
    ax2.set_title('Cumulative Distribution Function (CDF)', fontsize=12)
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(left=0)
    ax2.set_ylim(0, 1.05)

    # Add quantile markers to CDF
    if quantiles:
        for q, val in quantiles.items():
            color = colors.get(q, 'gray')
            ax2.axhline(q, color=color, linestyle=':', linewidth=1, alpha=0.5)
            ax2.axvline(val, color=color, linestyle='--', linewidth=1.5, alpha=0.7)
            ax2.plot(val, q, 'o', color=color, markersize=8)

    fig.suptitle(title, fontsize=14, fontweight='bold')
    fig.tight_layout()

    return fig


def save_figure(
    fig: "matplotlib.figure.Figure",
    output_dir: Path,
    basename: str,
    formats: List[str] = ['pdf', 'png', 'svg'],
    dpi: int = 150,
) -> Dict[str, Path]:
    """
    Save figure in multiple formats.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Figure to save.
    output_dir : Path
        Output directory.
    basename : str
        Base filename (without extension).
    formats : list
        Output formats (pdf, png, svg).
    dpi : int
        Resolution for raster formats.

    Returns
    -------
    dict
        Mapping format -> output path.
    """
    check_matplotlib()

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    paths = {}
    for fmt in formats:
        path = output_dir / f"{basename}.{fmt}"
        fig.savefig(path, format=fmt, dpi=dpi, bbox_inches='tight')
        paths[fmt] = path

    plt.close(fig)
    return paths


def generate_survey_plots(
    survey: "Survey",
    output_dir: Path,
    basename: Optional[str] = None,
    formats: List[str] = ['pdf', 'png', 'svg'],
) -> Dict[str, Path]:
    """
    Generate and save plots for a completed survey.

    Parameters
    ----------
    survey : Survey
        Completed survey instance (compute() must have been called).
    output_dir : Path
        Output directory.
    basename : str, optional
        Base filename. If None, uses survey name.
    formats : list
        Output formats.

    Returns
    -------
    dict
        Mapping format -> output path.
    """
    check_matplotlib()

    if basename is None:
        basename = survey.config.name.replace(' ', '_')

    results = survey.results
    centers = results['pdf_bin_centers']
    pdf = results['pdf']
    cdf = results['cdf']

    quantiles = {
        0.50: survey.quantile(0.50),
        0.95: survey.quantile(0.95),
        0.99: survey.quantile(0.99),
    }

    title = f"Survey: {survey.config.name}"
    fig = plot_pdf_cdf(centers, pdf, cdf, title=title, quantiles=quantiles)

    return save_figure(fig, output_dir, basename, formats=formats)
