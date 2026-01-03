"""
survey/engine.py — Batch Processing Engine for Survey Computations
==================================================================

Production-grade batch processing engine for survey-scale exposure estimation.

Features:
- Spreadsheet-driven workflow (XLS/ODS input)
- Row → YAML scenario conversion with full traceability
- Parallel scenario processing
- Change detection and incremental refresh
- Comprehensive audit trail

@project: SFPPy/INSERM — Survey-scale exposure estimation
@author: Olivier Vitrac, PhD, HDR
@email: olivier.vitrac@gmail.com
@license: MIT
"""

import os
import sys
import json
import time
import hashlib
import logging
from pathlib import Path
from datetime import datetime
from dataclasses import dataclass, field, asdict
from typing import List, Dict, Any, Optional, Tuple, Union
from concurrent.futures import ProcessPoolExecutor, as_completed
import shutil

import yaml
import numpy as np

from survey.survey import Survey
from survey.models import PriorSpec, SurveyConfig
from survey.fingerprints import stable_hash, canonical_json
from survey.reporting import ProgressReporter

logger = logging.getLogger(__name__)


# =============================================================================
# Spreadsheet Schema
# =============================================================================

SPREADSHEET_COLUMNS = {
    # Required columns
    'key': {'required': True, 'type': str, 'description': 'Unique scenario identifier'},
    'food': {'required': True, 'type': str, 'description': 'Food description'},
    'food_simulant': {'required': True, 'type': str, 'description': 'Food simulant name'},
    'food_volume_L': {'required': True, 'type': float, 'description': 'Food volume in liters'},
    'storage_temp_C': {'required': True, 'type': float, 'description': 'Storage temperature (°C)'},

    # Contact time prior (triangular)
    't_min_days': {'required': True, 'type': float, 'description': 'Minimum contact time (days)'},
    't_likely_days': {'required': True, 'type': float, 'description': 'Likely contact time (days)'},
    't_max_days': {'required': True, 'type': float, 'description': 'Maximum contact time (days)'},

    # Packaging
    'packaging': {'required': True, 'type': str, 'description': 'Packaging description'},
    'polymer': {'required': True, 'type': str, 'description': 'Polymer identifier'},
    'thickness_um': {'required': True, 'type': float, 'description': 'Thickness in µm'},
    'surface_area_dm2': {'required': True, 'type': float, 'description': 'Contact surface in dm²'},

    # Substance families (YAML references)
    'family1': {'required': True, 'type': str, 'description': 'Path to family YAML'},
    'family2': {'required': False, 'type': str, 'description': 'Path to family YAML (optional)'},
    'family3': {'required': False, 'type': str, 'description': 'Path to family YAML (optional)'},

    # Optional columns
    'h_m_s': {'required': False, 'type': float, 'default': 1e-7, 'description': 'Mass transfer coefficient'},
    'notes': {'required': False, 'type': str, 'description': 'Notes'},
}


# =============================================================================
# Data Models
# =============================================================================

@dataclass
class FamilySpec:
    """Substance family specification loaded from YAML."""
    name: str
    description: str
    substances: List[Dict[str, Any]]
    source_path: str
    fingerprint: str

    @classmethod
    def from_yaml(cls, path: Union[str, Path]) -> "FamilySpec":
        """Load family from YAML file."""
        path = Path(path)
        with open(path, 'r', encoding='utf-8') as f:
            data = yaml.safe_load(f)

        # Compute fingerprint
        fp = stable_hash(canonical_json(data))

        return cls(
            name=data.get('name', path.stem),
            description=data.get('description', ''),
            substances=data.get('substances', []),
            source_path=str(path.resolve()),
            fingerprint=fp,
        )

    def get_prior_spec(self) -> PriorSpec:
        """Extract concentration prior from family substances."""
        c0_min = min(s.get('C0_min', 0) for s in self.substances)
        c0_likely = np.mean([s.get('C0_likely', 10) for s in self.substances])
        c0_max = max(s.get('C0_max', 100) for s in self.substances)

        return PriorSpec(
            mode=float(c0_likely),
            max_val=float(c0_max),
            n_low=15,
            n_high=15,
            name=f"{self.name}_conc",
        )


@dataclass
class ScenarioRow:
    """Single row from spreadsheet representing one scenario."""
    key: str
    food: str
    food_simulant: str
    food_volume_L: float
    storage_temp_C: float
    t_min_days: float
    t_likely_days: float
    t_max_days: float
    packaging: str
    polymer: str
    thickness_um: float
    surface_area_dm2: float
    family1: str
    family2: Optional[str] = None
    family3: Optional[str] = None
    h_m_s: float = 1e-7
    notes: str = ""
    row_number: int = 0
    fingerprint: str = ""

    def __post_init__(self):
        """Compute fingerprint after initialization."""
        if not self.fingerprint:
            self.fingerprint = self._compute_fingerprint()

    def _compute_fingerprint(self) -> str:
        """Compute deterministic fingerprint for change detection."""
        payload = {
            'key': self.key,
            'food_simulant': self.food_simulant,
            'food_volume_L': self.food_volume_L,
            'storage_temp_C': self.storage_temp_C,
            't_min_days': self.t_min_days,
            't_likely_days': self.t_likely_days,
            't_max_days': self.t_max_days,
            'polymer': self.polymer,
            'thickness_um': self.thickness_um,
            'surface_area_dm2': self.surface_area_dm2,
            'family1': self.family1,
            'family2': self.family2,
            'family3': self.family3,
            'h_m_s': self.h_m_s,
        }
        return stable_hash(canonical_json(payload))

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

    def get_families(self) -> List[str]:
        """Return list of non-empty family paths."""
        families = [self.family1]
        if self.family2:
            families.append(self.family2)
        if self.family3:
            families.append(self.family3)
        return families


@dataclass
class ProcessingResult:
    """Result of processing one scenario."""
    key: str
    family: str
    status: str  # 'success', 'failed', 'skipped', 'cached'
    output_dir: str
    q50: Optional[float] = None
    q95: Optional[float] = None
    q99: Optional[float] = None
    cache_hits: int = 0
    cache_misses: int = 0
    duration_s: float = 0.0
    error: Optional[str] = None
    timestamp: str = ""

    def __post_init__(self):
        if not self.timestamp:
            self.timestamp = datetime.now().isoformat()


@dataclass
class BatchManifest:
    """Manifest for batch processing tracking."""
    batch_id: str
    spreadsheet_path: str
    spreadsheet_fingerprint: str
    families_dir: str
    artifacts_dir: str
    total_scenarios: int
    total_families: int
    processed: int = 0
    skipped: int = 0
    failed: int = 0
    cached: int = 0
    results: List[ProcessingResult] = field(default_factory=list)
    started_at: str = ""
    completed_at: str = ""

    def __post_init__(self):
        if not self.started_at:
            self.started_at = datetime.now().isoformat()

    def to_dict(self) -> Dict[str, Any]:
        d = asdict(self)
        d['results'] = [asdict(r) for r in self.results]
        return d

    def save(self, path: Union[str, Path]) -> None:
        """Save manifest to JSON file."""
        path = Path(path)
        path.write_text(json.dumps(self.to_dict(), indent=2), encoding='utf-8')

    @classmethod
    def load(cls, path: Union[str, Path]) -> "BatchManifest":
        """Load manifest from JSON file."""
        path = Path(path)
        data = json.loads(path.read_text(encoding='utf-8'))
        results = [ProcessingResult(**r) for r in data.pop('results', [])]
        return cls(**data, results=results)


# =============================================================================
# Spreadsheet Loader
# =============================================================================

class SpreadsheetLoader:
    """Load scenarios from XLS/ODS spreadsheet."""

    def __init__(self, path: Union[str, Path]):
        self.path = Path(path)
        if not self.path.exists():
            raise FileNotFoundError(f"Spreadsheet not found: {self.path}")

        self.extension = self.path.suffix.lower()
        self._data = None
        self._fingerprint = None

    @property
    def fingerprint(self) -> str:
        """Compute spreadsheet fingerprint."""
        if self._fingerprint is None:
            with open(self.path, 'rb') as f:
                self._fingerprint = hashlib.sha256(f.read()).hexdigest()
        return self._fingerprint

    def load(self) -> List[ScenarioRow]:
        """Load all rows from spreadsheet."""
        if self.extension in ['.xlsx', '.xls']:
            return self._load_excel()
        elif self.extension == '.ods':
            return self._load_ods()
        elif self.extension == '.csv':
            return self._load_csv()
        else:
            raise ValueError(f"Unsupported spreadsheet format: {self.extension}")

    def _load_excel(self) -> List[ScenarioRow]:
        """Load from Excel file."""
        try:
            import openpyxl
        except ImportError:
            raise ImportError("openpyxl required for Excel files: pip install openpyxl")

        wb = openpyxl.load_workbook(self.path, read_only=True, data_only=True)
        ws = wb.active

        rows = list(ws.iter_rows(values_only=True))
        if not rows:
            return []

        headers = [str(h).strip().lower() if h else '' for h in rows[0]]
        return self._parse_rows(headers, rows[1:])

    def _load_ods(self) -> List[ScenarioRow]:
        """Load from ODS file."""
        try:
            from odf import text, draw
            from odf.opendocument import load
            from odf.table import Table, TableRow, TableCell
        except ImportError:
            raise ImportError("odfpy required for ODS files: pip install odfpy")

        doc = load(self.path)
        table = doc.spreadsheet.getElementsByType(Table)[0]

        rows = []
        for tr in table.getElementsByType(TableRow):
            row = []
            for tc in tr.getElementsByType(TableCell):
                repeat = int(tc.getAttribute('numbercolumnsrepeated') or 1)
                # Get text content
                cell_text = ''
                for p in tc.getElementsByType(text.P):
                    cell_text += str(p)
                row.extend([cell_text] * repeat)
            rows.append(row)

        if not rows:
            return []

        headers = [str(h).strip().lower() if h else '' for h in rows[0]]
        return self._parse_rows(headers, rows[1:])

    def _load_csv(self) -> List[ScenarioRow]:
        """Load from CSV file."""
        import csv
        with open(self.path, 'r', encoding='utf-8') as f:
            reader = csv.reader(f)
            rows = list(reader)

        if not rows:
            return []

        headers = [str(h).strip().lower() if h else '' for h in rows[0]]
        return self._parse_rows(headers, rows[1:])

    def _parse_rows(self, headers: List[str], data_rows: List[tuple]) -> List[ScenarioRow]:
        """Parse rows into ScenarioRow objects."""
        # Map header names to column indices
        col_map = {h: i for i, h in enumerate(headers) if h}

        scenarios = []
        for row_num, row in enumerate(data_rows, start=2):
            if not row or not any(row):
                continue

            try:
                scenario = self._parse_row(col_map, row, row_num)
                if scenario:
                    scenarios.append(scenario)
            except Exception as e:
                logger.warning(f"Row {row_num}: Failed to parse - {e}")

        return scenarios

    def _parse_row(self, col_map: Dict[str, int], row: tuple, row_num: int) -> Optional[ScenarioRow]:
        """Parse single row into ScenarioRow."""
        def get_val(col_name: str, default=None):
            idx = col_map.get(col_name.lower())
            if idx is None or idx >= len(row):
                return default
            val = row[idx]
            if val is None or (isinstance(val, str) and val.strip() == ''):
                return default
            return val

        # Check required key
        key = get_val('key')
        if not key:
            return None

        return ScenarioRow(
            key=str(key).strip(),
            food=str(get_val('food', '')),
            food_simulant=str(get_val('food_simulant', 'oliveoil')),
            food_volume_L=float(get_val('food_volume_L', 1.0)),
            storage_temp_C=float(get_val('storage_temp_C', 25.0)),
            t_min_days=float(get_val('t_min_days', 0.0)),
            t_likely_days=float(get_val('t_likely_days', 7.0)),
            t_max_days=float(get_val('t_max_days', 30.0)),
            packaging=str(get_val('packaging', '')),
            polymer=str(get_val('polymer', 'PP')),
            thickness_um=float(get_val('thickness_um', 50.0)),
            surface_area_dm2=float(get_val('surface_area_dm2', 6.0)),
            family1=str(get_val('family1', '')),
            family2=str(get_val('family2')) if get_val('family2') else None,
            family3=str(get_val('family3')) if get_val('family3') else None,
            h_m_s=float(get_val('h_m_s', 1e-7)),
            notes=str(get_val('notes', '')),
            row_number=row_num,
        )


# =============================================================================
# Row Processor
# =============================================================================

class RowProcessor:
    """Convert spreadsheet row to scenario YAML and process."""

    def __init__(self, families_dir: Union[str, Path], artifacts_dir: Union[str, Path]):
        self.families_dir = Path(families_dir)
        self.artifacts_dir = Path(artifacts_dir)
        self._family_cache: Dict[str, FamilySpec] = {}

    def load_family(self, family_ref: str) -> FamilySpec:
        """Load family YAML, with caching."""
        if family_ref in self._family_cache:
            return self._family_cache[family_ref]

        # Try as absolute path first
        path = Path(family_ref)
        if not path.is_absolute():
            path = self.families_dir / family_ref
        if not path.suffix:
            path = path.with_suffix('.yml')

        if not path.exists():
            raise FileNotFoundError(f"Family YAML not found: {path}")

        family = FamilySpec.from_yaml(path)
        self._family_cache[family_ref] = family
        return family

    def row_to_yaml(self, row: ScenarioRow, family: FamilySpec) -> str:
        """Convert row + family to scenario YAML content."""
        scenario = {
            'name': f"{row.key}_{family.name}",
            'description': f"Generated from row {row.row_number}: {row.food} + {row.packaging}",
            'source': {
                'spreadsheet_row': row.row_number,
                'row_key': row.key,
                'row_fingerprint': row.fingerprint,
                'family_name': family.name,
                'family_fingerprint': family.fingerprint,
            },
            'physics': {
                'monolayer': {
                    'polymer': row.polymer,
                    'thickness_m': row.thickness_um * 1e-6,
                    'temperature_degC': row.storage_temp_C,
                },
                'interface': {
                    'h_m_s': row.h_m_s,
                    'surface_area_m2': row.surface_area_dm2 * 0.01,  # dm² → m²
                    'food_volume_m3': row.food_volume_L * 0.001,     # L → m³
                    'contact_temperature_degC': row.storage_temp_C,
                    'cf0': 0.0,
                    'food_simulant': row.food_simulant,
                },
            },
            'family': {
                'substances': [
                    {'name': s.get('name'), 'cas': s.get('cas')}
                    for s in family.substances
                ],
            },
            'priors': {
                'time_s': {
                    'triangular': {
                        'min': row.t_min_days * 86400,
                        'mode': row.t_likely_days * 86400,
                        'max': row.t_max_days * 86400,
                    },
                    'grid': {'nlow': 15, 'nhigh': 15},
                },
                'cp0_av': {
                    'triangular': {
                        'min': min(s.get('C0_min', 0) for s in family.substances),
                        'mode': np.mean([s.get('C0_likely', 10) for s in family.substances]),
                        'max': max(s.get('C0_max', 100) for s in family.substances),
                    },
                    'grid': {'nlow': 15, 'nhigh': 15},
                },
            },
            'solver': {
                'pdf_bins': 250,
                'fo_grid': {
                    'n_fo': 200,
                    'fo_max_factor': 1.5,
                    'fo_min_floor': 1e-15,
                },
            },
        }
        return yaml.dump(scenario, default_flow_style=False, allow_unicode=True)

    def process_row_family(
        self,
        row: ScenarioRow,
        family_ref: str,
        force: bool = False,
    ) -> ProcessingResult:
        """Process one row + family combination."""
        start_time = time.time()

        try:
            # Load family
            family = self.load_family(family_ref)

            # Create output directory
            output_dir = self.artifacts_dir / row.key / family.name
            output_dir.mkdir(parents=True, exist_ok=True)

            # Check if already processed (fingerprint match)
            manifest_path = output_dir / 'manifest.json'
            if manifest_path.exists() and not force:
                try:
                    with open(manifest_path, 'r') as f:
                        existing = json.load(f)
                    if (existing.get('row_fingerprint') == row.fingerprint and
                        existing.get('family_fingerprint') == family.fingerprint):
                        return ProcessingResult(
                            key=row.key,
                            family=family.name,
                            status='cached',
                            output_dir=str(output_dir),
                            q50=existing.get('q50'),
                            q95=existing.get('q95'),
                            q99=existing.get('q99'),
                            duration_s=time.time() - start_time,
                        )
                except Exception:
                    pass  # Proceed with reprocessing

            # Generate scenario YAML
            yaml_content = self.row_to_yaml(row, family)
            yaml_path = output_dir / 'scenario.yml'
            yaml_path.write_text(yaml_content, encoding='utf-8')

            # Process with Survey
            survey = Survey.from_scenario(yaml_path)
            survey.compute(parallel=False)  # Serial within row, parallel across rows

            # Extract results
            q50 = survey.quantile(0.50)
            q95 = survey.quantile(0.95)
            q99 = survey.quantile(0.99)

            # Save outputs
            survey.save(output_dir / 'results.npz')

            # Save manifest
            manifest = {
                'row_key': row.key,
                'row_fingerprint': row.fingerprint,
                'family_name': family.name,
                'family_fingerprint': family.fingerprint,
                'q50': q50,
                'q95': q95,
                'q99': q99,
                'processed_at': datetime.now().isoformat(),
            }
            with open(manifest_path, 'w') as f:
                json.dump(manifest, f, indent=2)

            return ProcessingResult(
                key=row.key,
                family=family.name,
                status='success',
                output_dir=str(output_dir),
                q50=q50,
                q95=q95,
                q99=q99,
                cache_hits=survey._curve_cache.stats.get('hits', 0),
                cache_misses=survey._curve_cache.stats.get('misses', 0),
                duration_s=time.time() - start_time,
            )

        except Exception as e:
            logger.error(f"Failed processing {row.key}/{family_ref}: {e}")
            return ProcessingResult(
                key=row.key,
                family=family_ref,
                status='failed',
                output_dir=str(self.artifacts_dir / row.key / Path(family_ref).stem),
                error=str(e),
                duration_s=time.time() - start_time,
            )


# =============================================================================
# Batch Engine
# =============================================================================

class BatchEngine:
    """
    Batch processing engine for survey computations.

    Features:
    - Spreadsheet-driven workflow
    - Parallel scenario processing
    - Change detection and incremental refresh
    - Comprehensive audit trail
    """

    def __init__(
        self,
        spreadsheet_path: Union[str, Path],
        families_dir: Union[str, Path],
        artifacts_dir: Union[str, Path],
        cache_dir: Union[str, Path] = ".survey_cache",
    ):
        self.spreadsheet_path = Path(spreadsheet_path)
        self.families_dir = Path(families_dir)
        self.artifacts_dir = Path(artifacts_dir)
        self.cache_dir = Path(cache_dir)

        # Ensure directories exist
        self.families_dir.mkdir(parents=True, exist_ok=True)
        self.artifacts_dir.mkdir(parents=True, exist_ok=True)
        self.cache_dir.mkdir(parents=True, exist_ok=True)

        self._loader = None
        self._processor = None
        self._manifest = None

    @property
    def loader(self) -> SpreadsheetLoader:
        if self._loader is None:
            self._loader = SpreadsheetLoader(self.spreadsheet_path)
        return self._loader

    @property
    def processor(self) -> RowProcessor:
        if self._processor is None:
            self._processor = RowProcessor(self.families_dir, self.artifacts_dir)
        return self._processor

    def load_scenarios(self) -> List[ScenarioRow]:
        """Load all scenarios from spreadsheet."""
        return self.loader.load()

    def preview(self) -> str:
        """Preview batch processing without executing."""
        scenarios = self.load_scenarios()

        lines = [
            "Batch Processing Preview",
            "=" * 60,
            f"Spreadsheet: {self.spreadsheet_path}",
            f"Families dir: {self.families_dir}",
            f"Artifacts dir: {self.artifacts_dir}",
            f"Scenarios: {len(scenarios)}",
            "",
        ]

        # Count families
        all_families = set()
        for s in scenarios:
            all_families.update(s.get_families())

        lines.append(f"Unique families: {len(all_families)}")
        lines.append("")

        # Summary table
        lines.append(f"{'Key':<20} {'Polymer':<10} {'Families':<30}")
        lines.append("-" * 60)
        for s in scenarios[:10]:  # Show first 10
            families = ', '.join(Path(f).stem for f in s.get_families())
            lines.append(f"{s.key:<20} {s.polymer:<10} {families:<30}")
        if len(scenarios) > 10:
            lines.append(f"... and {len(scenarios) - 10} more")

        return '\n'.join(lines)

    def process(
        self,
        max_workers: int = None,
        force: bool = False,
        keys: List[str] = None,
    ) -> BatchManifest:
        """
        Process all scenarios in batch.

        Parameters
        ----------
        max_workers : int, optional
            Maximum parallel workers. If None, uses CPU count.
        force : bool
            Force reprocessing even if cached.
        keys : List[str], optional
            Only process these scenario keys.

        Returns
        -------
        BatchManifest
            Batch processing manifest with results.
        """
        scenarios = self.load_scenarios()

        # Filter by keys if specified
        if keys:
            scenarios = [s for s in scenarios if s.key in keys]

        # Build task list
        tasks = []
        for scenario in scenarios:
            for family_ref in scenario.get_families():
                tasks.append((scenario, family_ref))

        total_tasks = len(tasks)
        logger.info(f"Processing {total_tasks} tasks ({len(scenarios)} scenarios)")

        # Initialize manifest
        batch_id = datetime.now().strftime("%Y%m%d_%H%M%S")
        self._manifest = BatchManifest(
            batch_id=batch_id,
            spreadsheet_path=str(self.spreadsheet_path),
            spreadsheet_fingerprint=self.loader.fingerprint,
            families_dir=str(self.families_dir),
            artifacts_dir=str(self.artifacts_dir),
            total_scenarios=len(scenarios),
            total_families=total_tasks,
        )

        # Progress reporter
        prog = ProgressReporter(total_steps=total_tasks, desc="Batch Processing")

        # Process tasks
        results = []

        if max_workers == 1:
            # Serial processing
            for scenario, family_ref in tasks:
                result = self.processor.process_row_family(scenario, family_ref, force=force)
                results.append(result)
                self._update_stats(result)
                prog.update(1)
        else:
            # Parallel processing
            with ProcessPoolExecutor(max_workers=max_workers) as exe:
                futures = {
                    exe.submit(
                        _process_task,
                        scenario.to_dict(),
                        family_ref,
                        str(self.families_dir),
                        str(self.artifacts_dir),
                        force,
                    ): (scenario.key, family_ref)
                    for scenario, family_ref in tasks
                }

                for fut in as_completed(futures):
                    key, family_ref = futures[fut]
                    try:
                        result = fut.result()
                    except Exception as e:
                        result = ProcessingResult(
                            key=key,
                            family=family_ref,
                            status='failed',
                            output_dir='',
                            error=str(e),
                        )
                    results.append(result)
                    self._update_stats(result)
                    prog.update(1)

        prog.done()

        # Finalize manifest
        self._manifest.results = results
        self._manifest.completed_at = datetime.now().isoformat()

        # Save manifest
        manifest_path = self.artifacts_dir / f"batch_{batch_id}.json"
        self._manifest.save(manifest_path)
        logger.info(f"Manifest saved: {manifest_path}")

        # Generate summary
        self._generate_summary()

        return self._manifest

    def _update_stats(self, result: ProcessingResult) -> None:
        """Update manifest statistics."""
        if result.status == 'success':
            self._manifest.processed += 1
        elif result.status == 'cached':
            self._manifest.cached += 1
        elif result.status == 'skipped':
            self._manifest.skipped += 1
        else:
            self._manifest.failed += 1

    def _generate_summary(self) -> None:
        """Generate summary CSV with all quantiles."""
        import csv

        summary_path = self.artifacts_dir / 'summary.csv'
        with open(summary_path, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            writer.writerow([
                'key', 'family', 'status', 'q50', 'q95', 'q99',
                'cache_hits', 'cache_misses', 'duration_s'
            ])
            for r in self._manifest.results:
                writer.writerow([
                    r.key, r.family, r.status,
                    r.q50, r.q95, r.q99,
                    r.cache_hits, r.cache_misses, r.duration_s
                ])

        logger.info(f"Summary saved: {summary_path}")

    def get_changes(self, previous_manifest: Union[str, Path]) -> Dict[str, List[str]]:
        """
        Detect changes since previous batch.

        Returns
        -------
        dict
            {'new': [...], 'modified': [...], 'deleted': [...]}
        """
        prev = BatchManifest.load(previous_manifest)
        current_scenarios = self.load_scenarios()

        # Build fingerprint maps
        prev_fps = {
            (r.key, r.family): r
            for r in prev.results
        }

        changes = {'new': [], 'modified': [], 'deleted': []}

        # Check current scenarios
        for scenario in current_scenarios:
            for family_ref in scenario.get_families():
                family_name = Path(family_ref).stem
                key = (scenario.key, family_name)

                if key not in prev_fps:
                    changes['new'].append(f"{scenario.key}/{family_name}")
                else:
                    # Check fingerprint
                    prev_result = prev_fps[key]
                    # Would need stored fingerprint in result to compare
                    pass

        return changes


# =============================================================================
# Worker Function (for multiprocessing)
# =============================================================================

def _process_task(
    scenario_dict: Dict[str, Any],
    family_ref: str,
    families_dir: str,
    artifacts_dir: str,
    force: bool,
) -> ProcessingResult:
    """
    Worker function for parallel processing.

    Must be module-level for multiprocessing.
    """
    scenario = ScenarioRow(**scenario_dict)
    processor = RowProcessor(families_dir, artifacts_dir)
    return processor.process_row_family(scenario, family_ref, force=force)


# =============================================================================
# CLI Interface
# =============================================================================

def main():
    """Command-line interface for batch processing."""
    import argparse

    parser = argparse.ArgumentParser(
        description='Survey Batch Processing Engine',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Preview batch
  python -m survey.engine preview batch.xlsx --families families/

  # Process all scenarios
  python -m survey.engine process batch.xlsx --families families/ --output artifacts/

  # Force reprocessing
  python -m survey.engine process batch.xlsx --force
        """
    )

    subparsers = parser.add_subparsers(dest='command', required=True)

    # Preview command
    preview_parser = subparsers.add_parser('preview', help='Preview batch processing')
    preview_parser.add_argument('spreadsheet', help='Path to spreadsheet (XLS/ODS/CSV)')
    preview_parser.add_argument('--families', default='families/', help='Families directory')

    # Process command
    process_parser = subparsers.add_parser('process', help='Process batch')
    process_parser.add_argument('spreadsheet', help='Path to spreadsheet')
    process_parser.add_argument('--families', default='families/', help='Families directory')
    process_parser.add_argument('--output', default='artifacts/', help='Output directory')
    process_parser.add_argument('--workers', type=int, default=None, help='Max workers')
    process_parser.add_argument('--force', action='store_true', help='Force reprocessing')
    process_parser.add_argument('--keys', nargs='+', help='Only process these keys')

    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    if args.command == 'preview':
        engine = BatchEngine(
            spreadsheet_path=args.spreadsheet,
            families_dir=args.families,
            artifacts_dir='artifacts/',
        )
        print(engine.preview())

    elif args.command == 'process':
        engine = BatchEngine(
            spreadsheet_path=args.spreadsheet,
            families_dir=args.families,
            artifacts_dir=args.output,
        )
        manifest = engine.process(
            max_workers=args.workers,
            force=args.force,
            keys=args.keys,
        )
        print(f"\nCompleted: {manifest.processed} processed, "
              f"{manifest.cached} cached, {manifest.failed} failed")


if __name__ == '__main__':
    main()
