# CHANGELOG - Survey Module

All notable changes to the Survey module are documented in this file.

## [0.3.0] - 2026-01-10

### Fixed - CAS-based k/k0 Inference

#### survey/models.py
- **CAS-based migrant lookup**: `_get_migrant()` now prioritizes CAS numbers over substance names for PubChem queries
  - CAS lookup is more reliable than name lookup (substance names in databases often don't match PubChem)
  - Falls back to name lookup if CAS fails
  - Caches results by CAS or name to avoid redundant queries

- **`infer_k()` and `infer_k0()`**: Now accept optional `cas` parameter
  - Passes CAS to `_compute_k_from_layer()` and `_compute_k0_from_food()`
  - Enables reliable Flory-Huggins computation for substances with non-standard names

- **`infer_all()`**: Now passes both `substance.name` and `substance.cas` to inference methods

- **Polymer name normalization**: Added `.strip()` to handle trailing spaces in polymer codes (e.g., "PET " → "PET")

### Impact

| Metric | Before (name lookup) | After (CAS lookup) |
|--------|---------------------|-------------------|
| k=1.0 fallback | 50.3% | 5.3% |
| k0=1.0 fallback | 48.2% | 1.7% |
| k=k0=1.0 (both) | 48.2% | 1.7% |

### Remaining Limitations

Substances still getting k=k0=1.0 (~1.7%):
- **Missing logP in PubChem**: Organometallic compounds, salts, complex mixtures
- **Unsupported polymers**: Paper/paperboard (no FH monomer defined)
- **Example**: Dibutyltin dilaurate (CAS 77-58-7) - PubChem lacks XLogP data

See **Known Limitations** section in README.md for details.

---

## [0.2.0] - 2026-01-09

### Added - Risk Assessment Percentiles

#### survey/aggregation.py
- **Risk-relevant percentiles** (`compute_risk_percentiles`): Computes standard percentiles for regulatory risk assessment
  - P50: Median (central tendency)
  - P75: Third quartile
  - P90: High exposure
  - P95: Reasonable worst case (commonly used regulatory threshold)
  - P97.5: Used in EFSA/regulatory frameworks
  - P99: Extreme exposure
- **RISK_PERCENTILES** constant: `[50, 75, 90, 95, 97.5, 99]`
- **`quantile_from_cdf()`**: Extracts quantiles from integrated PDF histogram (for visualization consistency)

### Changed
- **`quantile_from_samples()`**: Enhanced documentation clarifying this is the preferred method for risk assessment (exact empirical CDF, no binning error)
- **`survey/__init__.py`**: Exports `compute_risk_percentiles`, `RISK_PERCENTILES`, `quantile_from_cdf`

### Technical Notes

#### Two Quantile Computation Methods

| Method | Use Case | Accuracy |
|--------|----------|----------|
| `quantile_from_samples()` | **Risk assessment** | Exact empirical CDF from weighted samples |
| `quantile_from_cdf()` | Visualization | Binned approximation from PDF histogram |

For **regulatory reporting**, always use `compute_risk_percentiles()` or `quantile_from_samples()` which compute exact values without histogram discretization error.

#### Dual-Resolution Output Strategy

Production workflow supports:
1. **High-resolution PDF/CDF** (250-500 bins): Stored in NPZ for detailed plotting
2. **Low-resolution percentiles** (50 points): For Excel aggregation (P2, P4, ..., P100)
3. **Risk percentiles** (6 points): P50, P75, P90, P95, P97.5, P99 for regulatory reporting

### API Usage

```python
from survey.aggregation import compute_risk_percentiles, RISK_PERCENTILES

# Compute risk-relevant percentiles
risk_pct = compute_risk_percentiles(cf_samples, weights)
print(f"P95 (regulatory threshold): {risk_pct['p95']:.4f} mg/kg")
print(f"P97.5 (EFSA framework): {risk_pct['p97_5']:.4f} mg/kg")

# All risk percentiles
for level in RISK_PERCENTILES:
    key = f"p{str(level).replace('.', '_')}"
    print(f"P{level}: {risk_pct[key]:.4f}")
```

---

## [0.1.0] - 2026-01-09

### Initial Release

#### Core Features
- **Survey class**: Production-grade migration estimation with caching
- **Deterministic uncertainty propagation**: Finite-difference quadrature on triangular priors (no Monte Carlo)
- **Parallel master curve computation**: Distributed across workers with persistent cache
- **Content-addressed caching**: SHA-256 fingerprinting prevents redundant computation
- **Multilayer support**: Functional barriers with reference layer selection

#### Aggregation Module (`survey/aggregation.py`)
- `combine_tensors()`: Tensor product for two independent components
- `combine_multiple_tensors()`: Sequential tensor product for N components
- `aggregate_components()`: Per-substance aggregation across components
- `aggregate_family_weighted()`: Occurrence-weighted family aggregation
- `compute_pdf_from_samples()`: Weighted histogram to PDF/CDF
- `compute_statistics()`: Mean, std, q50, q95, q99, max
- `compute_percentiles()`: Configurable resolution CDF discretization
- `aggregate_packaging_components()`: Multi-component food packaging aggregation
- `aggregate_food_exposure()`: Family-level exposure aggregation

#### Web Applications
- **Family Editor** (port 8000): Spreadsheet-style substance family management
- **Survey Simulator** (port 8001): Interactive batch simulation with PDF/CDF visualization

#### Supported Polymers
| Polymer | Code | Notes |
|---------|------|-------|
| Polypropylene | PP | Common packaging |
| Low-density polyethylene | LDPE | Films |
| High-density polyethylene | HDPE | Bottles |
| Polyethylene terephthalate | PET, gPET, wPET | Bottles |
| Polystyrene | PS | Foam packaging |
| Polyamide | PA | Barrier layers |
| Ethylene vinyl acetate | EVA | Sealing layers |
| Polyvinyl chloride | PVC | Various |
| Paper/Paperboard | Paper | Cellulose-based |

#### Supported Food Simulants
| Simulant | EU Reference | Description |
|----------|--------------|-------------|
| `water` | Simulant A | Aqueous foods |
| `water3aceticacid` | Simulant B | Acidic foods |
| `ethanol` | Simulant C | Alcoholic beverages |
| `oliveoil` | Simulant D | Fatty foods |
| `ethanol50` | Simulant D1 | Dairy products |
| `tenax` | Simulant E | Dry foods |

---

## License

MIT License. See LICENSE file in repository root.

## Author

**Olivier Vitrac, PhD, HDR**
- Email: olivier.vitrac@gmail.com
- Affiliation: INRAE / Generative Simulation
