#!/usr/bin/env python3
"""
run_batch.py — Batch Simulation Runner for SFPPy Survey Module
==============================================================

Run migration simulations from YAML scenarios or XLSX spreadsheets.
Supports parallel execution, visualization, and comprehensive reporting.

Usage Examples:
---------------
    # Run all YAML files in a directory
    python survey/run_batch.py scenarios/

    # Run a single YAML scenario
    python survey/run_batch.py scenario.yml

    # Run from XLSX spreadsheet (Family Editor format)
    python survey/run_batch.py families.xlsx

    # Run from PF Jobs JSON (Survey Simulator export)
    python survey/run_batch.py pf_jobs.json

    # Custom output directory and parallel workers
    python survey/run_batch.py scenarios/ -o results/ -w 4

    # Quiet mode (no progress output)
    python survey/run_batch.py scenarios/ -q

Supported Input Formats:
------------------------
1. YAML scenarios (.yml, .yaml) — Direct Survey.from_scenario() format
2. XLSX spreadsheets (.xlsx) — Family Editor export format
3. JSON PF jobs (.json) — Survey Simulator export format

Output Structure:
-----------------
    output_dir/
    ├── SUMMARY.md              # Markdown report with all results
    ├── SUMMARY.json            # Machine-readable summary
    ├── scenario_name_1/
    │   ├── scenario_name_1.npz        # Raw results (PDF, CDF, samples)
    │   ├── scenario_name_1_manifest.json  # Metadata and fingerprints
    │   ├── scenario_name_1.pdf        # PDF plot (if matplotlib available)
    │   ├── scenario_name_1.png        # PNG plot (if matplotlib available)
    │   └── scenario_name_1.svg        # SVG plot (if matplotlib available)
    ├── scenario_name_2/
    │   └── ...
    └── ...

@project: SFPPy — Survey-scale exposure estimation
@author: Olivier Vitrac, PhD, HDR
@email: olivier.vitrac@gmail.com
@license: MIT
"""

import sys
import os
import argparse
import json
import yaml
import time
import traceback
from pathlib import Path
from datetime import datetime
from typing import List, Dict, Any, Optional, Tuple
from dataclasses import dataclass, field, asdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing as mp

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))


# ============================================================================
# Data Models
# ============================================================================

@dataclass
class BatchResult:
    """Result from a single scenario simulation."""
    name: str
    success: bool
    error: Optional[str] = None
    duration_s: float = 0.0
    quantiles: Dict[str, float] = field(default_factory=dict)
    n_substances: int = 0
    polymer: str = ""
    simulant: str = ""
    thickness_um: float = 0.0
    output_dir: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class BatchSummary:
    """Summary of batch execution."""
    total: int = 0
    success: int = 0
    failed: int = 0
    duration_s: float = 0.0
    results: List[BatchResult] = field(default_factory=list)
    timestamp: str = ""
    input_path: str = ""
    output_dir: str = ""

    def to_dict(self) -> Dict[str, Any]:
        return {
            "total": self.total,
            "success": self.success,
            "failed": self.failed,
            "duration_s": self.duration_s,
            "timestamp": self.timestamp,
            "input_path": self.input_path,
            "output_dir": self.output_dir,
            "results": [r.to_dict() for r in self.results],
        }


# ============================================================================
# Input Loaders
# ============================================================================

def load_yaml_scenarios(path: Path) -> List[Tuple[str, Path]]:
    """
    Load YAML scenario files from path.

    Parameters
    ----------
    path : Path
        Path to YAML file or directory containing YAML files.

    Returns
    -------
    List[Tuple[str, Path]]
        List of (scenario_name, scenario_path) tuples.
    """
    scenarios = []

    if path.is_file():
        if path.suffix.lower() in ('.yml', '.yaml'):
            name = path.stem
            scenarios.append((name, path))
    elif path.is_dir():
        for yaml_file in sorted(path.glob('*.yml')) + sorted(path.glob('*.yaml')):
            name = yaml_file.stem
            scenarios.append((name, yaml_file))

    return scenarios


def load_xlsx_scenarios(xlsx_path: Path, temp_dir: Path) -> List[Tuple[str, Path]]:
    """
    Load scenarios from XLSX spreadsheet (Family Editor format).

    Converts each family in the spreadsheet to a YAML scenario file.

    Parameters
    ----------
    xlsx_path : Path
        Path to XLSX file.
    temp_dir : Path
        Directory to write temporary YAML files.

    Returns
    -------
    List[Tuple[str, Path]]
        List of (scenario_name, scenario_path) tuples.
    """
    from survey.spreadsheet import read_spreadsheet, SpreadsheetData

    data = read_spreadsheet(xlsx_path)
    scenarios = []

    for family in data.families:
        # Create scenario YAML for each family
        scenario = {
            "name": family.name,
            "physics": {
                "monolayer": {
                    "polymer": family.polymer,
                    "thickness_m": family.thickness_um * 1e-6,
                    "temperature_degC": family.temperature_C,
                },
                "interface": {
                    "h_m_s": 1e-7,
                    "surface_area_m2": 0.06,
                    "food_volume_m3": family.food_volume_ml * 1e-6,
                    "contact_temperature_degC": family.temperature_C,
                    "cf0": 0.0,
                    "food_simulant": family.food_simulant or "ethanol50",
                },
            },
            "priors": {
                "time_s": {
                    "triangular": {
                        "min": 0.0,
                        "mode": family.contact_days * 86400,
                        "max": family.contact_days * 2 * 86400,
                    },
                    "grid": {"nlow": 15, "nhigh": 15},
                },
                "cp0_av": {
                    "triangular": {
                        "min": 0.0,
                        "mode": 50.0,
                        "max": 200.0,
                    },
                    "grid": {"nlow": 15, "nhigh": 15},
                },
            },
            "family": {
                "substances": [
                    {"cas": sub.cas} if sub.cas else {"name": sub.name}
                    for sub in family.substances
                ]
            },
            "solver": {
                "pdf_bins": 250,
                "fo_grid": {"n_fo": 200, "fo_max_factor": 1.5},
            },
        }

        # Write YAML
        yaml_path = temp_dir / f"{family.name.replace(' ', '_')}.yml"
        with open(yaml_path, 'w') as f:
            yaml.dump(scenario, f, default_flow_style=False, allow_unicode=True)

        scenarios.append((family.name, yaml_path))

    return scenarios


def load_json_pf_jobs(json_path: Path, temp_dir: Path) -> List[Tuple[str, Path]]:
    """
    Load scenarios from PF Jobs JSON (Survey Simulator export format).

    Parameters
    ----------
    json_path : Path
        Path to JSON file.
    temp_dir : Path
        Directory to write temporary YAML files.

    Returns
    -------
    List[Tuple[str, Path]]
        List of (scenario_name, scenario_path) tuples.
    """
    with open(json_path, 'r') as f:
        pf_jobs = json.load(f)

    # Ensure it's a list
    if isinstance(pf_jobs, dict) and 'jobs' in pf_jobs:
        pf_jobs = pf_jobs['jobs']
    elif not isinstance(pf_jobs, list):
        pf_jobs = [pf_jobs]

    scenarios = []

    for job in pf_jobs:
        name = job.get('name', f"job_{len(scenarios)}")

        # Parse units (handle both µm and um)
        thickness_unit = job.get('thickness_unit', 'µm')
        thickness_val = job.get('thickness_value', 100)
        if 'µm' in thickness_unit or 'um' in thickness_unit:
            thickness_m = thickness_val * 1e-6
        elif 'mm' in thickness_unit:
            thickness_m = thickness_val * 1e-3
        else:
            thickness_m = thickness_val * 1e-6  # assume µm

        # Parse food volume
        volume_unit = job.get('food_volume_unit', 'mL')
        volume_val = job.get('food_volume_value', 1000)
        if 'mL' in volume_unit or 'ml' in volume_unit:
            volume_m3 = volume_val * 1e-6
        elif 'L' in volume_unit:
            volume_m3 = volume_val * 1e-3
        else:
            volume_m3 = volume_val * 1e-6  # assume mL

        # Parse surface area
        area_unit = job.get('surface_area_unit', 'cm²')
        area_val = job.get('surface_area_value', 600)
        if 'cm²' in area_unit or 'cm2' in area_unit:
            area_m2 = area_val * 1e-4
        elif 'dm²' in area_unit or 'dm2' in area_unit:
            area_m2 = area_val * 1e-2
        elif 'm²' in area_unit or 'm2' in area_unit:
            area_m2 = area_val
        else:
            area_m2 = area_val * 1e-4  # assume cm²

        # Parse time
        time_unit = job.get('time_unit', 'days')
        time_mode = job.get('time_mode', 30)
        time_max = job.get('time_max', time_mode * 2)
        if 'day' in time_unit.lower():
            time_mode_s = time_mode * 86400
            time_max_s = time_max * 86400
        elif 'hour' in time_unit.lower():
            time_mode_s = time_mode * 3600
            time_max_s = time_max * 3600
        else:
            time_mode_s = time_mode * 86400  # assume days
            time_max_s = time_max * 86400

        # Build substances with concentration priors
        substances_list = job.get('substances', [])
        c0_values = []
        family_substances = []

        for sub in substances_list:
            identifier = sub.get('identifier') or sub.get('cas') or sub.get('name')
            c0_mode = sub.get('c0_mode', 50)
            c0_values.append(c0_mode)

            # Try CAS first, then name
            if identifier and '-' in str(identifier):
                family_substances.append({"cas": identifier})
            else:
                family_substances.append({"name": identifier})

        # Compute concentration prior from substances
        import numpy as np
        c0_mode_avg = float(np.mean(c0_values)) if c0_values else 50.0
        c0_max_avg = float(np.max([sub.get('c0_max', 200) for sub in substances_list])) if substances_list else 200.0

        scenario = {
            "name": name,
            "physics": {
                "monolayer": {
                    "polymer": job.get('polymer', 'LDPE'),
                    "thickness_m": thickness_m,
                    "temperature_degC": job.get('layer_temperature_C', 25),
                },
                "interface": {
                    "h_m_s": job.get('h_m_s', 1e-7),
                    "surface_area_m2": area_m2,
                    "food_volume_m3": volume_m3,
                    "contact_temperature_degC": job.get('contact_temperature_C', 25),
                    "cf0": 0.0,
                    "food_simulant": job.get('simulant', 'ethanol50'),
                },
            },
            "priors": {
                "time_s": {
                    "triangular": {
                        "min": 0.0,
                        "mode": time_mode_s,
                        "max": time_max_s,
                    },
                    "grid": {"nlow": 15, "nhigh": 15},
                },
                "cp0_av": {
                    "triangular": {
                        "min": 0.0,
                        "mode": c0_mode_avg,
                        "max": c0_max_avg,
                    },
                    "grid": {"nlow": 15, "nhigh": 15},
                },
            },
            "family": {
                "substances": family_substances if family_substances else [{"mass_g_mol": 200}]
            },
            "solver": {
                "pdf_bins": 250,
                "fo_grid": {"n_fo": 200, "fo_max_factor": 1.5},
            },
        }

        # Write YAML
        safe_name = name.replace(' ', '_').replace('/', '_')
        yaml_path = temp_dir / f"{safe_name}.yml"
        with open(yaml_path, 'w') as f:
            yaml.dump(scenario, f, default_flow_style=False, allow_unicode=True)

        scenarios.append((name, yaml_path))

    return scenarios


# ============================================================================
# Simulation Runner
# ============================================================================

def run_single_scenario(args: Tuple[str, Path, Path, bool]) -> BatchResult:
    """
    Run a single scenario simulation (worker function).

    Parameters
    ----------
    args : tuple
        (scenario_name, scenario_path, output_dir, generate_plots)

    Returns
    -------
    BatchResult
        Simulation result.
    """
    scenario_name, scenario_path, output_dir, generate_plots = args
    start_time = time.time()

    try:
        # Import here to avoid multiprocessing issues
        from survey import Survey
        from survey.tables import extract_survey_result

        # Load and compute survey
        survey = Survey.from_scenario(scenario_path)
        survey.compute(parallel=False)  # No nested parallelism

        # Extract configuration info
        config = survey.config
        polymer = config.packaging.layers[0].polymer if config.packaging.layers else "unknown"
        simulant = config.packaging.food_simulant
        thickness_um = config.packaging.layers[0].thickness_m * 1e6 if config.packaging.layers else 0
        n_substances = len(survey.substances)

        # Extract quantiles
        quantiles = {
            "q50": float(survey.quantile(0.50)),
            "q75": float(survey.quantile(0.75)),
            "q90": float(survey.quantile(0.90)),
            "q95": float(survey.quantile(0.95)),
            "q99": float(survey.quantile(0.99)),
        }

        # Create output directory
        scenario_dir = output_dir / scenario_name.replace(' ', '_')
        scenario_dir.mkdir(parents=True, exist_ok=True)

        # Save NPZ results
        basename = scenario_name.replace(' ', '_')
        npz_path = scenario_dir / f"{basename}.npz"
        survey.save(npz_path)

        # Save manifest
        manifest_path = scenario_dir / f"{basename}_manifest.json"
        survey.save_manifest(manifest_path)

        # Generate visualizations if requested
        if generate_plots:
            try:
                from survey.visualization import generate_survey_plots, HAS_MATPLOTLIB
                if HAS_MATPLOTLIB:
                    generate_survey_plots(
                        survey,
                        scenario_dir,
                        basename=basename,
                        formats=['pdf', 'png', 'svg'],
                    )
            except ImportError:
                pass  # matplotlib not available

        duration = time.time() - start_time

        return BatchResult(
            name=scenario_name,
            success=True,
            duration_s=duration,
            quantiles=quantiles,
            n_substances=n_substances,
            polymer=polymer,
            simulant=simulant,
            thickness_um=thickness_um,
            output_dir=str(scenario_dir),
        )

    except Exception as e:
        duration = time.time() - start_time
        return BatchResult(
            name=scenario_name,
            success=False,
            error=str(e),
            duration_s=duration,
        )


def run_batch(
    scenarios: List[Tuple[str, Path]],
    output_dir: Path,
    n_workers: int = 1,
    generate_plots: bool = True,
    verbose: bool = True,
) -> BatchSummary:
    """
    Run batch simulations.

    Parameters
    ----------
    scenarios : List[Tuple[str, Path]]
        List of (scenario_name, scenario_path) tuples.
    output_dir : Path
        Output directory.
    n_workers : int
        Number of parallel workers.
    generate_plots : bool
        Generate PDF/PNG/SVG plots.
    verbose : bool
        Print progress.

    Returns
    -------
    BatchSummary
        Summary of batch execution.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    summary = BatchSummary(
        total=len(scenarios),
        timestamp=datetime.now().isoformat(),
        output_dir=str(output_dir),
    )

    start_time = time.time()

    # Prepare worker arguments
    worker_args = [
        (name, path, output_dir, generate_plots)
        for name, path in scenarios
    ]

    if n_workers > 1 and len(scenarios) > 1:
        # Parallel execution
        if verbose:
            print(f"Running {len(scenarios)} scenarios with {n_workers} workers...")

        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            futures = {
                executor.submit(run_single_scenario, args): args[0]
                for args in worker_args
            }

            for i, future in enumerate(as_completed(futures), 1):
                name = futures[future]
                try:
                    result = future.result()
                    summary.results.append(result)

                    if result.success:
                        summary.success += 1
                        status = "OK"
                    else:
                        summary.failed += 1
                        status = f"FAIL: {result.error}"

                    if verbose:
                        print(f"  [{i}/{len(scenarios)}] {name}: {status} ({result.duration_s:.1f}s)")

                except Exception as e:
                    summary.failed += 1
                    summary.results.append(BatchResult(
                        name=name,
                        success=False,
                        error=str(e),
                    ))
                    if verbose:
                        print(f"  [{i}/{len(scenarios)}] {name}: FAIL: {e}")
    else:
        # Sequential execution
        if verbose:
            print(f"Running {len(scenarios)} scenarios sequentially...")

        for i, args in enumerate(worker_args, 1):
            name = args[0]
            result = run_single_scenario(args)
            summary.results.append(result)

            if result.success:
                summary.success += 1
                status = "OK"
            else:
                summary.failed += 1
                status = f"FAIL: {result.error}"

            if verbose:
                print(f"  [{i}/{len(scenarios)}] {name}: {status} ({result.duration_s:.1f}s)")

    summary.duration_s = time.time() - start_time

    return summary


# ============================================================================
# Report Generation
# ============================================================================

def generate_markdown_report(summary: BatchSummary) -> str:
    """Generate Markdown summary report."""
    lines = [
        "# SFPPy Batch Simulation Report",
        "",
        f"**Generated:** {summary.timestamp}",
        f"**Duration:** {summary.duration_s:.1f}s",
        "",
        "## Summary",
        "",
        f"| Metric | Value |",
        f"|--------|-------|",
        f"| Total scenarios | {summary.total} |",
        f"| Successful | {summary.success} |",
        f"| Failed | {summary.failed} |",
        "",
        "## Results",
        "",
        "| Scenario | Polymer | Simulant | Thickness (µm) | Substances | Q95 (mg/kg) | Status |",
        "|----------|---------|----------|----------------|------------|-------------|--------|",
    ]

    for r in sorted(summary.results, key=lambda x: x.name):
        if r.success:
            q95 = r.quantiles.get('q95', 0)
            status = "✅"
            lines.append(
                f"| {r.name} | {r.polymer} | {r.simulant} | {r.thickness_um:.0f} | "
                f"{r.n_substances} | {q95:.4f} | {status} |"
            )
        else:
            lines.append(
                f"| {r.name} | - | - | - | - | - | ❌ {r.error[:30]}... |"
            )

    # Detailed quantiles for successful runs
    lines.extend([
        "",
        "## Quantile Details",
        "",
        "| Scenario | Q50 | Q75 | Q90 | Q95 | Q99 |",
        "|----------|-----|-----|-----|-----|-----|",
    ])

    for r in sorted(summary.results, key=lambda x: x.name):
        if r.success:
            lines.append(
                f"| {r.name} | {r.quantiles.get('q50', 0):.4f} | "
                f"{r.quantiles.get('q75', 0):.4f} | {r.quantiles.get('q90', 0):.4f} | "
                f"{r.quantiles.get('q95', 0):.4f} | {r.quantiles.get('q99', 0):.4f} |"
            )

    return "\n".join(lines)


# ============================================================================
# Main Entry Point
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Run batch migration simulations from YAML/XLSX/JSON",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s scenarios/                  # Run all YAML in directory
  %(prog)s scenario.yml                # Run single YAML scenario
  %(prog)s families.xlsx               # Run from XLSX spreadsheet
  %(prog)s pf_jobs.json                # Run from PF Jobs JSON
  %(prog)s scenarios/ -o results/ -w 4 # Custom output and workers
        """,
    )

    parser.add_argument(
        "input",
        type=Path,
        help="Input path: YAML file, YAML directory, XLSX file, or JSON file",
    )
    parser.add_argument(
        "-o", "--output",
        type=Path,
        default=None,
        help="Output directory (default: input_name_results/)",
    )
    parser.add_argument(
        "-w", "--workers",
        type=int,
        default=None,
        help="Number of parallel workers (default: CPU_count/2)",
    )
    parser.add_argument(
        "--no-plots",
        action="store_true",
        help="Skip plot generation",
    )
    parser.add_argument(
        "-q", "--quiet",
        action="store_true",
        help="Suppress progress output",
    )

    args = parser.parse_args()

    input_path = args.input.resolve()
    verbose = not args.quiet
    generate_plots = not args.no_plots

    # Determine output directory
    if args.output:
        output_dir = args.output.resolve()
    else:
        output_dir = input_path.parent / f"{input_path.stem}_results"

    # Determine number of workers
    if args.workers:
        n_workers = args.workers
    else:
        n_workers = max(1, (os.cpu_count() or 2) // 2)

    # Print header
    if verbose:
        print(f"""
╔══════════════════════════════════════════════════════════════╗
║           SFPPy Batch Simulation Runner                      ║
╚══════════════════════════════════════════════════════════════╝

Input:   {input_path}
Output:  {output_dir}
Workers: {n_workers}
Plots:   {'Yes' if generate_plots else 'No'}
""")

    # Load scenarios based on input type
    scenarios = []
    temp_dir = None

    try:
        if input_path.suffix.lower() in ('.yml', '.yaml'):
            # Single YAML file
            scenarios = load_yaml_scenarios(input_path)
        elif input_path.suffix.lower() == '.xlsx':
            # XLSX spreadsheet
            import tempfile
            temp_dir = Path(tempfile.mkdtemp())
            scenarios = load_xlsx_scenarios(input_path, temp_dir)
        elif input_path.suffix.lower() == '.json':
            # JSON PF jobs
            import tempfile
            temp_dir = Path(tempfile.mkdtemp())
            scenarios = load_json_pf_jobs(input_path, temp_dir)
        elif input_path.is_dir():
            # Directory of YAML files
            scenarios = load_yaml_scenarios(input_path)
        else:
            print(f"[ERROR] Unsupported input format: {input_path.suffix}")
            sys.exit(1)

        if not scenarios:
            print(f"[ERROR] No scenarios found in {input_path}")
            sys.exit(1)

        if verbose:
            print(f"Found {len(scenarios)} scenarios to process\n")

        # Run batch
        summary = run_batch(
            scenarios,
            output_dir,
            n_workers=n_workers,
            generate_plots=generate_plots,
            verbose=verbose,
        )
        summary.input_path = str(input_path)

        # Generate reports
        markdown_report = generate_markdown_report(summary)

        # Save reports
        report_md_path = output_dir / "SUMMARY.md"
        with open(report_md_path, 'w') as f:
            f.write(markdown_report)

        report_json_path = output_dir / "SUMMARY.json"
        with open(report_json_path, 'w') as f:
            json.dump(summary.to_dict(), f, indent=2)

        # Print summary
        if verbose:
            print(f"""
╔══════════════════════════════════════════════════════════════╗
║                     Batch Complete                           ║
╠══════════════════════════════════════════════════════════════╣
║  Total:     {summary.total:4d}                                          ║
║  Success:   {summary.success:4d}                                          ║
║  Failed:    {summary.failed:4d}                                          ║
║  Duration:  {summary.duration_s:6.1f}s                                      ║
╠══════════════════════════════════════════════════════════════╣
║  Reports:                                                    ║
║    {str(report_md_path)[:54]:<54s} ║
║    {str(report_json_path)[:54]:<54s} ║
╚══════════════════════════════════════════════════════════════╝
""")

        # Exit with error code if any failures
        sys.exit(0 if summary.failed == 0 else 1)

    finally:
        # Cleanup temp directory
        if temp_dir and temp_dir.exists():
            import shutil
            shutil.rmtree(temp_dir, ignore_errors=True)


if __name__ == "__main__":
    main()
