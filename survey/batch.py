#!/usr/bin/env python3
"""
Batch Simulation Runner for SFPPy Survey Module
================================================

Simplified batch runner for the web simulator app.
Uses the existing Survey class and io module.

Author: Olivier Vitrac
"""

import os
import sys
import time
import traceback
import tempfile
import multiprocessing as mp
from pathlib import Path
from typing import Dict, List, Optional, Any, Callable, Union
from dataclasses import dataclass, field
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
import json
import yaml

import numpy as np

# Import survey modules
from survey.survey import Survey
from survey.io import load_scenario, load_substances_from_scenario
from survey.models import (
    LayerSpec,
    PackagingSpec,
    SubstanceSpec,
    PriorSpec,
    SurveyConfig,
)
from survey.spreadsheet import SpreadsheetData, read_spreadsheet


@dataclass
class SimulationTask:
    """A single simulation task."""
    task_id: str
    family_name: str
    config: Dict[str, Any]
    output_dir: Optional[Path] = None

    def to_dict(self) -> Dict[str, Any]:
        return {
            "task_id": self.task_id,
            "family_name": self.family_name,
            "config": self.config,
            "output_dir": str(self.output_dir) if self.output_dir else None,
        }


@dataclass
class SimulationResult:
    """Result of a single simulation."""
    task_id: str
    family_name: str
    success: bool
    error: Optional[str] = None
    duration_s: float = 0.0
    n_samples: int = 0
    quantiles: Dict[str, float] = field(default_factory=dict)
    output_files: List[str] = field(default_factory=list)
    # Store centers and distributions for visualization
    centers: Optional[List[float]] = None
    pdf: Optional[List[float]] = None
    cdf: Optional[List[float]] = None

    def to_dict(self) -> Dict[str, Any]:
        return {
            "task_id": self.task_id,
            "family_name": self.family_name,
            "success": self.success,
            "error": self.error,
            "duration_s": self.duration_s,
            "n_samples": self.n_samples,
            "quantiles": self.quantiles,
            "output_files": self.output_files,
            "centers": self.centers,
            "pdf": self.pdf,
            "cdf": self.cdf,
        }


@dataclass
class BatchProgress:
    """Progress tracking for batch execution."""
    total: int = 0
    completed: int = 0
    failed: int = 0
    running: int = 0
    start_time: Optional[float] = None
    results: List[SimulationResult] = field(default_factory=list)

    @property
    def percent(self) -> float:
        if self.total == 0:
            return 0.0
        return 100.0 * self.completed / self.total

    @property
    def elapsed_s(self) -> float:
        if self.start_time is None:
            return 0.0
        return time.time() - self.start_time

    def to_dict(self) -> Dict[str, Any]:
        return {
            "total": self.total,
            "completed": self.completed,
            "failed": self.failed,
            "running": self.running,
            "percent": self.percent,
            "elapsed_s": self.elapsed_s,
            "results": [r.to_dict() for r in self.results],
        }


def get_cpu_count() -> int:
    """Get number of available CPUs."""
    return os.cpu_count() or 1


def get_default_workers() -> int:
    """Get default number of workers (50% of CPUs)."""
    return max(1, get_cpu_count() // 2)


def _create_scenario_yaml(config: Dict[str, Any], temp_dir: Path) -> Path:
    """
    Create a temporary scenario YAML file from config dictionary.

    The scenario format matches what Survey.from_scenario() expects.
    This uses the exact structure from io.py:
    - physics.monolayer or physics.multilayer
    - physics.interface
    - priors.time_s and priors.cp0_av
    - solver options
    - family with substances or masses_g_mol
    """
    # Extract parameters
    name = config.get("name", "survey")
    polymer = config.get("polymer", "LDPE")
    thickness_um = config.get("thickness_um", 100)
    temperature_C = config.get("temperature_C", 25)
    contact_days = config.get("contact_days", 10)
    food_volume_ml = config.get("food_volume_ml", 1000)
    food_density = config.get("food_density", 1.0)
    food_simulant = config.get("food_simulant") or config.get("simulant", "ethanol50")
    substances = config.get("substances", {})

    # Convert units
    thickness_m = thickness_um * 1e-6
    food_volume_m3 = food_volume_ml * 1e-6
    contact_s = contact_days * 86400  # days to seconds

    # Compute concentration prior from substances
    c0_values = []
    for sub_config in substances.values():
        c0_values.append(sub_config.get("C0_likely", 100))

    # Convert to native Python floats to avoid numpy serialization issues
    c0_mode = float(np.mean(c0_values)) if c0_values else 100.0
    c0_max = float(max(sub_config.get("C0_max", 1000) for sub_config in substances.values())) if substances else 1000.0

    # Build scenario structure matching io.py format
    scenario = {
        "name": name,

        # Physics section
        "physics": {
            "monolayer": {
                "polymer": polymer,
                "thickness_m": thickness_m,
                "C0": 0.0,  # Actual C0 from priors
                "temperature_degC": temperature_C,
            },
            "interface": {
                "h_m_s": 1e-7,
                "surface_area_m2": 0.06,  # ~600 cm² typical
                "food_volume_m3": food_volume_m3,
                "contact_temperature_degC": temperature_C,
                "cf0": 0.0,
                "food_simulant": food_simulant,
            },
        },

        # Priors section
        "priors": {
            "time_s": {
                "triangular": {
                    "mode": contact_s,
                    "max": contact_s * 2,
                },
                "grid": {
                    "nlow": 15,
                    "nhigh": 15,
                },
            },
            "cp0_av": {
                "triangular": {
                    "mode": c0_mode,
                    "max": c0_max,
                },
                "grid": {
                    "nlow": 15,
                    "nhigh": 15,
                },
            },
        },

        # Solver settings
        "solver": {
            "pdf_bins": 250,
            "fo_grid": {
                "n_fo": 200,
                "fo_max_factor": 1.5,
                "fo_min_floor": 1e-15,
            },
            "cache_dir": str(temp_dir / ".cache"),
        },

        # Family with substances
        "family": {
            "substances": [
                {"cas": identifier}
                for identifier in substances.keys()
            ]
        },
    }

    # Write YAML
    yaml_path = temp_dir / f"{name}.yml"
    with open(yaml_path, "w", encoding="utf-8") as f:
        yaml.dump(scenario, f, default_flow_style=False, allow_unicode=True)

    return yaml_path


def _run_single_simulation(task: SimulationTask) -> SimulationResult:
    """
    Run a single simulation task (worker function).

    This function runs in a separate process.
    """
    start_time = time.time()

    try:
        config = task.config
        n_samples = config.get("n_samples", 1000)

        # Create temporary directory for scenario
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)

            # Create scenario YAML
            yaml_path = _create_scenario_yaml(config, temp_path)

            # Load and run survey
            survey = Survey.from_scenario(yaml_path)
            survey.compute(parallel=False)  # No nested parallelism

            # Extract results (using .results property, not method)
            results = survey.results
            centers = results['pdf_bin_centers'].tolist()
            pdf = results['pdf'].tolist()
            cdf = results['cdf'].tolist()

            # Compute quantiles
            quantiles = {
                "q50": float(survey.quantile(0.50)),
                "q75": float(survey.quantile(0.75)),
                "q90": float(survey.quantile(0.90)),
                "q95": float(survey.quantile(0.95)),
                "q99": float(survey.quantile(0.99)),
                "mean": float(np.average(results['CF_samples'], weights=results['weights'])),
                "std": float(np.sqrt(np.average(
                    (results['CF_samples'] - np.average(results['CF_samples'], weights=results['weights']))**2,
                    weights=results['weights']
                ))),
                "min": float(np.min(results['CF_samples'])),
                "max": float(np.max(results['CF_samples'])),
            }

        # Generate output files if output_dir specified
        output_files = []
        if task.output_dir:
            task.output_dir.mkdir(parents=True, exist_ok=True)

            # Save plots using visualization module
            try:
                from survey.visualization import generate_survey_plots
                # Create a simple wrapper object for the plot function
                class SurveyWrapper:
                    def __init__(self, centers, pdf, cdf, config):
                        self.config = type('Config', (), {'name': config.get('name', 'survey')})()
                        self.total_migration = None  # Not used directly
                        self._results = {
                            'pdf_bin_centers': np.array(centers),
                            'pdf': np.array(pdf),
                            'cdf': np.array(cdf),
                        }

                    def get_results(self):
                        return self._results

                wrapper = SurveyWrapper(centers, pdf, cdf, config)
                plot_files = generate_survey_plots(
                    wrapper,
                    task.output_dir,
                    task.family_name,
                    formats=["png", "svg"]
                )
                output_files.extend([str(f) for f in plot_files])
            except Exception as e:
                pass  # Plots are optional

        duration = time.time() - start_time

        return SimulationResult(
            task_id=task.task_id,
            family_name=task.family_name,
            success=True,
            duration_s=duration,
            n_samples=n_samples,
            quantiles=quantiles,
            output_files=output_files,
            centers=centers,
            pdf=pdf,
            cdf=cdf,
        )

    except Exception as e:
        duration = time.time() - start_time
        return SimulationResult(
            task_id=task.task_id,
            family_name=task.family_name,
            success=False,
            error=f"{type(e).__name__}: {str(e)}\n{traceback.format_exc()}",
            duration_s=duration,
        )


class BatchRunner:
    """
    Batch simulation runner with multiprocessing support.

    Usage:
        runner = BatchRunner(n_workers=4)
        runner.add_task(task)
        results = runner.run()
    """

    def __init__(
        self,
        n_workers: Optional[int] = None,
        output_dir: Optional[Path] = None,
    ):
        """
        Initialize batch runner.

        Args:
            n_workers: Number of parallel workers. Default: 50% of CPUs.
            output_dir: Base directory for output files.
        """
        self.n_workers = n_workers or get_default_workers()
        self.output_dir = Path(output_dir) if output_dir else None
        self.tasks: List[SimulationTask] = []
        self.progress = BatchProgress()
        self._progress_callback: Optional[Callable[[BatchProgress], None]] = None

    def add_task(self, task: SimulationTask) -> None:
        """Add a simulation task."""
        if self.output_dir and task.output_dir is None:
            task.output_dir = self.output_dir / task.family_name
        self.tasks.append(task)

    def add_from_config(self, config: Dict[str, Any], task_id: Optional[str] = None) -> None:
        """Add a task from a configuration dictionary."""
        family_name = config.get("name", f"family_{len(self.tasks)}")
        task = SimulationTask(
            task_id=task_id or f"task_{len(self.tasks):04d}",
            family_name=family_name,
            config=config,
        )
        self.add_task(task)

    def add_from_spreadsheet(self, spreadsheet_path: Union[str, Path]) -> int:
        """
        Add tasks from a spreadsheet file.

        Returns:
            Number of tasks added.
        """
        data = read_spreadsheet(spreadsheet_path)
        configs = data.to_survey_configs()

        for i, config in enumerate(configs):
            self.add_from_config(config, task_id=f"task_{len(self.tasks):04d}")

        return len(configs)

    def add_from_yaml_folder(self, folder_path: Union[str, Path]) -> int:
        """
        Add tasks from all YAML files in a folder.

        Returns:
            Number of tasks added.
        """
        folder = Path(folder_path)
        count = 0

        for yaml_file in sorted(folder.glob("*.yml")) + sorted(folder.glob("*.yaml")):
            try:
                with open(yaml_file, 'r', encoding='utf-8') as f:
                    data = yaml.safe_load(f)

                # Convert family format to survey config format
                config = {
                    "name": data.get("name", yaml_file.stem),
                    "description": data.get("description", ""),
                    "polymer": data.get("polymer", "generic"),
                    "thickness_um": data.get("thickness_um", 100),
                    "temperature_C": data.get("temperature_C", 25),
                    "contact_days": data.get("contact_days", 10),
                    "food_volume_ml": data.get("food_volume_ml", 1000),
                    "food_density": data.get("food_density", 1.0),
                    "substances": {},
                }

                # Add substances
                substances = data.get("substances", {})
                if isinstance(substances, dict):
                    config["substances"] = substances
                elif isinstance(substances, list):
                    for sub in substances:
                        if isinstance(sub, dict) and "id" in sub:
                            config["substances"][sub["id"]] = sub

                self.add_from_config(config, task_id=f"task_{len(self.tasks):04d}")
                count += 1
            except Exception as e:
                print(f"Warning: Failed to load {yaml_file}: {e}", file=sys.stderr)

        return count

    def set_progress_callback(self, callback: Callable[[BatchProgress], None]) -> None:
        """Set a callback function for progress updates."""
        self._progress_callback = callback

    def _update_progress(self) -> None:
        """Update progress and call callback if set."""
        if self._progress_callback:
            self._progress_callback(self.progress)

    def run(self, n_samples: int = 1000) -> List[SimulationResult]:
        """
        Run all tasks in parallel.

        Args:
            n_samples: Number of Monte Carlo samples per simulation.

        Returns:
            List of simulation results.
        """
        if not self.tasks:
            return []

        # Update configs with n_samples
        for task in self.tasks:
            task.config["n_samples"] = n_samples

        # Initialize progress
        self.progress = BatchProgress(
            total=len(self.tasks),
            start_time=time.time(),
        )
        self._update_progress()

        results = []

        # Use ProcessPoolExecutor for parallel execution
        with ProcessPoolExecutor(max_workers=self.n_workers) as executor:
            # Submit all tasks
            future_to_task = {
                executor.submit(_run_single_simulation, task): task
                for task in self.tasks
            }

            self.progress.running = len(future_to_task)
            self._update_progress()

            # Collect results as they complete
            for future in as_completed(future_to_task):
                task = future_to_task[future]

                try:
                    result = future.result()
                except Exception as e:
                    result = SimulationResult(
                        task_id=task.task_id,
                        family_name=task.family_name,
                        success=False,
                        error=f"Process error: {str(e)}",
                    )

                results.append(result)
                self.progress.results.append(result)
                self.progress.completed += 1
                self.progress.running -= 1

                if not result.success:
                    self.progress.failed += 1

                self._update_progress()

        return results

    def run_sequential(self, n_samples: int = 1000) -> List[SimulationResult]:
        """
        Run all tasks sequentially (for debugging).

        Args:
            n_samples: Number of Monte Carlo samples per simulation.

        Returns:
            List of simulation results.
        """
        if not self.tasks:
            return []

        # Update configs with n_samples
        for task in self.tasks:
            task.config["n_samples"] = n_samples

        # Initialize progress
        self.progress = BatchProgress(
            total=len(self.tasks),
            start_time=time.time(),
        )
        self._update_progress()

        results = []

        for task in self.tasks:
            self.progress.running = 1
            self._update_progress()

            result = _run_single_simulation(task)
            results.append(result)
            self.progress.results.append(result)
            self.progress.completed += 1
            self.progress.running = 0

            if not result.success:
                self.progress.failed += 1

            self._update_progress()

        return results


def run_batch_from_spreadsheet(
    spreadsheet_path: Union[str, Path],
    output_dir: Optional[Union[str, Path]] = None,
    n_workers: Optional[int] = None,
    n_samples: int = 1000,
    progress_callback: Optional[Callable[[BatchProgress], None]] = None,
) -> List[SimulationResult]:
    """
    Convenience function to run batch simulations from a spreadsheet.

    Args:
        spreadsheet_path: Path to XLSX or ODS file.
        output_dir: Directory for output files.
        n_workers: Number of parallel workers (default: 50% of CPUs).
        n_samples: Number of Monte Carlo samples.
        progress_callback: Optional callback for progress updates.

    Returns:
        List of simulation results.
    """
    runner = BatchRunner(
        n_workers=n_workers,
        output_dir=Path(output_dir) if output_dir else None,
    )

    runner.add_from_spreadsheet(spreadsheet_path)

    if progress_callback:
        runner.set_progress_callback(progress_callback)

    return runner.run(n_samples=n_samples)


def run_batch_from_yaml_folder(
    folder_path: Union[str, Path],
    output_dir: Optional[Union[str, Path]] = None,
    n_workers: Optional[int] = None,
    n_samples: int = 1000,
    progress_callback: Optional[Callable[[BatchProgress], None]] = None,
) -> List[SimulationResult]:
    """
    Convenience function to run batch simulations from YAML files.

    Args:
        folder_path: Directory containing YAML family files.
        output_dir: Directory for output files.
        n_workers: Number of parallel workers (default: 50% of CPUs).
        n_samples: Number of Monte Carlo samples.
        progress_callback: Optional callback for progress updates.

    Returns:
        List of simulation results.
    """
    runner = BatchRunner(
        n_workers=n_workers,
        output_dir=Path(output_dir) if output_dir else None,
    )

    runner.add_from_yaml_folder(folder_path)

    if progress_callback:
        runner.set_progress_callback(progress_callback)

    return runner.run(n_samples=n_samples)
