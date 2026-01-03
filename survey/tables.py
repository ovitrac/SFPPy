"""
survey/tables.py — Markdown Table Generation for Survey Results
===============================================================

Generates markdown tables summarizing survey inputs and outputs.

@project: SFPPy/INSERM — Survey-scale exposure estimation
@author: Olivier Vitrac, PhD, HDR
@email: olivier.vitrac@gmail.com
@license: MIT
"""

from pathlib import Path
from typing import List, Dict, Any, Optional
from dataclasses import dataclass
import numpy as np
from datetime import datetime


@dataclass
class SurveyResult:
    """Container for survey results with inputs."""
    name: str
    # Packaging inputs
    polymer: str
    thickness_um: float
    surface_area_dm2: float
    food_volume_L: float
    temperature_C: float
    food_simulant: str
    n_layers: int
    i_ref: int
    # Prior inputs
    time_mode_days: float
    time_max_days: float
    conc_mode_mg_kg: float
    conc_max_mg_kg: float
    # Substances
    n_substances: int
    substances: List[str]
    # Results
    mean: float
    q50: float
    q75: float
    q90: float
    q95: float
    q99: float
    cf_max: float


def extract_survey_result(survey: "Survey") -> SurveyResult:
    """
    Extract inputs and results from a Survey object.

    Parameters
    ----------
    survey : Survey
        Completed survey instance.

    Returns
    -------
    SurveyResult
        Extracted data container.
    """
    config = survey.config
    pkg = config.packaging
    ref_layer = survey.ref_layer
    results = survey.results

    # Compute statistics
    weights = results['weights']
    cf_samples = results['CF_samples']
    mean = float(np.average(cf_samples, weights=weights))
    cf_max = float(np.max(cf_samples))

    return SurveyResult(
        name=config.name,
        polymer=ref_layer.polymer,
        thickness_um=pkg.total_thickness_m * 1e6,
        surface_area_dm2=pkg.surface_area_m2 * 100,  # m² → dm²
        food_volume_L=pkg.food_volume_m3 * 1000,     # m³ → L
        temperature_C=pkg.contact_temperature_degC,
        food_simulant=pkg.food_simulant,
        n_layers=pkg.n_layers,
        i_ref=survey._i_ref,
        time_mode_days=config.time_prior.mode / 86400,
        time_max_days=config.time_prior.max_val / 86400,
        conc_mode_mg_kg=config.conc_prior.mode,
        conc_max_mg_kg=config.conc_prior.max_val,
        n_substances=survey.n_substances,
        substances=[s.id for s in survey.substances],
        mean=mean,
        q50=survey.quantile(0.50),
        q75=survey.quantile(0.75),
        q90=survey.quantile(0.90),
        q95=survey.quantile(0.95),
        q99=survey.quantile(0.99),
        cf_max=cf_max,
    )


def format_inputs_table(results: List[SurveyResult]) -> str:
    """
    Generate markdown table for survey inputs.

    Parameters
    ----------
    results : List[SurveyResult]
        List of survey results.

    Returns
    -------
    str
        Markdown table string.
    """
    lines = [
        "## Survey Inputs",
        "",
        "| Scenario | Polymer | Thickness (µm) | Area (dm²) | Volume (L) | Temp (°C) | Simulant | Layers | t_mode (d) | t_max (d) | C0_mode | C0_max | N_sub |",
        "|----------|---------|----------------|------------|------------|-----------|----------|--------|------------|-----------|---------|--------|-------|",
    ]

    for r in results:
        lines.append(
            f"| {r.name} | {r.polymer} | {r.thickness_um:.1f} | "
            f"{r.surface_area_dm2:.2f} | {r.food_volume_L:.3f} | "
            f"{r.temperature_C:.0f} | {r.food_simulant} | {r.n_layers} | "
            f"{r.time_mode_days:.1f} | {r.time_max_days:.1f} | "
            f"{r.conc_mode_mg_kg:.1f} | {r.conc_max_mg_kg:.1f} | {r.n_substances} |"
        )

    return "\n".join(lines)


def format_outputs_table(results: List[SurveyResult]) -> str:
    """
    Generate markdown table for survey outputs (quantiles).

    Parameters
    ----------
    results : List[SurveyResult]
        List of survey results.

    Returns
    -------
    str
        Markdown table string.
    """
    lines = [
        "## Survey Results — Migration Quantiles (mg/kg food)",
        "",
        "| Scenario | Mean | q50 | q75 | q90 | q95 | q99 | Max |",
        "|----------|------|-----|-----|-----|-----|-----|-----|",
    ]

    for r in results:
        lines.append(
            f"| {r.name} | {r.mean:.4f} | {r.q50:.4f} | {r.q75:.4f} | "
            f"{r.q90:.4f} | {r.q95:.4f} | {r.q99:.4f} | {r.cf_max:.4f} |"
        )

    return "\n".join(lines)


def format_substance_table(results: List[SurveyResult]) -> str:
    """
    Generate markdown table listing substances per scenario.

    Parameters
    ----------
    results : List[SurveyResult]
        List of survey results.

    Returns
    -------
    str
        Markdown table string.
    """
    lines = [
        "## Substances per Scenario",
        "",
        "| Scenario | N | Substances |",
        "|----------|---|------------|",
    ]

    for r in results:
        subs = ", ".join(r.substances[:5])
        if len(r.substances) > 5:
            subs += f" ... (+{len(r.substances) - 5} more)"
        lines.append(f"| {r.name} | {r.n_substances} | {subs} |")

    return "\n".join(lines)


def generate_summary_report(
    results: List[SurveyResult],
    title: str = "Survey Batch Results",
    output_dir: Optional[Path] = None,
) -> str:
    """
    Generate complete markdown summary report.

    Parameters
    ----------
    results : List[SurveyResult]
        List of survey results.
    title : str
        Report title.
    output_dir : Path, optional
        Output directory (for relative image paths).

    Returns
    -------
    str
        Complete markdown report.
    """
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    sections = [
        f"# {title}",
        "",
        f"Generated: {timestamp}",
        "",
        "---",
        "",
        format_inputs_table(results),
        "",
        "---",
        "",
        format_outputs_table(results),
        "",
        "---",
        "",
        format_substance_table(results),
        "",
    ]

    # Add image references if output_dir provided
    if output_dir:
        sections.extend([
            "---",
            "",
            "## Visualizations",
            "",
        ])
        for r in results:
            # Use scenario name for subdirectory, which matches how run_all_examples saves files
            # Files are in: output/<scenario_stem>/<scenario_stem>.png
            # The scenario name from YAML may differ from filename stem
            basename = r.name.replace(' ', '_')
            sections.append(f"### {r.name}")
            sections.append("")
            sections.append(f"![{r.name}]({basename}/{basename}.png)")
            sections.append("")

    sections.extend([
        "---",
        "",
        "*Generated by SFPPy Survey Module*",
    ])

    return "\n".join(sections)


def save_summary_report(
    results: List[SurveyResult],
    output_path: Path,
    title: str = "Survey Batch Results",
) -> None:
    """
    Save markdown summary report to file.

    Parameters
    ----------
    results : List[SurveyResult]
        List of survey results.
    output_path : Path
        Output file path.
    title : str
        Report title.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    report = generate_summary_report(
        results, title=title, output_dir=output_path.parent
    )
    output_path.write_text(report, encoding='utf-8')
