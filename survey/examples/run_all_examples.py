#!/usr/bin/env python3
"""
run_all_examples.py — Test Runner for Survey Examples
=====================================================

Processes all example scenarios and generates:
- PDF, PNG, SVG visualizations
- Markdown summary tables with inputs and quantiles

Usage:
    python run_all_examples.py [--output-dir OUTPUT]

@project: SFPPy/INSERM — Survey-scale exposure estimation
@author: Olivier Vitrac, PhD, HDR
@email: olivier.vitrac@gmail.com
@license: MIT
"""

import sys
import argparse
from pathlib import Path
from datetime import datetime

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from survey import Survey
from survey.visualization import generate_survey_plots, HAS_MATPLOTLIB
from survey.tables import extract_survey_result, generate_summary_report, save_summary_report, SurveyResult
from typing import List, Dict, Any


# Example scenarios to process
EXAMPLES = [
    "basic_monolayer.yml",
    "multilayer_barrier.yml",
]


def find_examples_dir() -> Path:
    """Find the examples directory."""
    script_dir = Path(__file__).resolve().parent
    return script_dir


def process_example(
    example_path: Path,
    output_dir: Path,
    verbose: bool = True,
) -> SurveyResult:
    """
    Process a single example scenario.

    Parameters
    ----------
    example_path : Path
        Path to example YAML file.
    output_dir : Path
        Output directory for artifacts.
    verbose : bool
        Print progress information.

    Returns
    -------
    SurveyResult
        Extracted results.
    """
    if verbose:
        print(f"\n{'='*60}")
        print(f"Processing: {example_path.name}")
        print(f"{'='*60}")

    # Load and compute survey
    survey = Survey.from_scenario(example_path)
    if verbose:
        print(survey.preview())

    survey.compute(parallel=True)

    if verbose:
        print("\n" + survey.summary())

    # Generate visualizations
    # Use scenario name from YAML for consistency with markdown tables
    basename = survey.config.name.replace(' ', '_')
    scenario_dir = output_dir / basename

    if HAS_MATPLOTLIB:
        if verbose:
            print(f"\nGenerating visualizations...")
        paths = generate_survey_plots(
            survey,
            scenario_dir,
            basename=basename,
            formats=['pdf', 'png', 'svg'],
        )
        if verbose:
            for fmt, path in paths.items():
                print(f"  {fmt.upper()}: {path}")
    else:
        if verbose:
            print("\n[WARNING] matplotlib not available - skipping visualizations")

    # Save NPZ results
    npz_path = scenario_dir / f"{basename}.npz"
    survey.save(npz_path)
    if verbose:
        print(f"  NPZ: {npz_path}")

    # Save manifest
    manifest_path = scenario_dir / f"{basename}_manifest.json"
    survey.save_manifest(manifest_path)
    if verbose:
        print(f"  Manifest: {manifest_path}")

    # Extract results for summary
    return extract_survey_result(survey)


def run_all_examples(
    output_dir: Path,
    examples: List[str] = EXAMPLES,
    verbose: bool = True,
) -> List[SurveyResult]:
    """
    Run all example scenarios.

    Parameters
    ----------
    output_dir : Path
        Output directory for artifacts.
    examples : list
        List of example filenames.
    verbose : bool
        Print progress information.

    Returns
    -------
    List[SurveyResult]
        Results from all examples.
    """
    examples_dir = find_examples_dir()
    results = []

    for example in examples:
        example_path = examples_dir / example
        if not example_path.exists():
            if verbose:
                print(f"[WARNING] Example not found: {example_path}")
            continue

        try:
            result = process_example(example_path, output_dir, verbose=verbose)
            results.append(result)
        except Exception as e:
            if verbose:
                print(f"[ERROR] Failed to process {example}: {e}")
            import traceback
            traceback.print_exc()

    return results


def main():
    parser = argparse.ArgumentParser(
        description="Run all survey examples and generate outputs"
    )
    parser.add_argument(
        "--output-dir", "-o",
        type=Path,
        default=Path(__file__).parent / "output",
        help="Output directory for artifacts (default: examples/output)",
    )
    parser.add_argument(
        "--quiet", "-q",
        action="store_true",
        help="Suppress verbose output",
    )
    args = parser.parse_args()

    output_dir = args.output_dir.resolve()
    verbose = not args.quiet

    print(f"""
╔══════════════════════════════════════════════════════════════╗
║         SFPPy Survey Examples — Test Runner                  ║
╚══════════════════════════════════════════════════════════════╝

Output directory: {output_dir}
Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
""")

    # Run all examples
    results = run_all_examples(output_dir, verbose=verbose)

    if not results:
        print("\n[ERROR] No examples were processed successfully.")
        sys.exit(1)

    # Generate summary report
    print(f"\n{'='*60}")
    print("Generating Summary Report")
    print(f"{'='*60}")

    report_path = output_dir / "SUMMARY.md"
    save_summary_report(results, report_path, title="Survey Examples — Batch Results")
    print(f"\nSummary report: {report_path}")

    # Print summary tables to console
    report = generate_summary_report(results)
    print("\n" + report)

    print(f"""
╔══════════════════════════════════════════════════════════════╗
║                     Test Run Complete                        ║
╠══════════════════════════════════════════════════════════════╣
║  Examples processed: {len(results):3d}                                     ║
║  Output directory:   {str(output_dir)[:40]:<40s} ║
╚══════════════════════════════════════════════════════════════╝
""")


if __name__ == "__main__":
    main()
