# Survey Examples

This directory contains example files and output from the SFPPy Survey module.

## Directory Structure

```
examples/
├── README.md                 # This file
├── batch_scenarios/          # YAML scenario templates
├── batch_results/            # Output from YAML batch runs
├── batch_results_json/       # Output from JSON batch runs
├── families/                 # Family definition files
└── output/                   # General simulation output
```

## Subdirectories

### `batch_scenarios/`
Pre-built YAML scenario templates for common food packaging applications:
- `dairy_PET_pot.yml` — Dairy products in PET containers
- `water_PET_bottle.yml` — Water in PET bottles
- `water_HDPE_cap.yml` — Water contact with HDPE caps
- `fatty_PP_tray.yml` — Fatty foods in PP trays
- `hot_fill_PET.yml` — Hot-fill applications in PET

### `batch_results/`
Output from running `python survey/run_batch.py examples/batch_scenarios/`:
- `SUMMARY.md` — Markdown report with all results
- `SUMMARY.json` — Machine-readable summary
- Individual scenario folders with PDF/CDF plots and data

### `families/`
Family definition spreadsheets for use with Family Editor:
- `plasticizers.xlsx` — Common plasticizer families
- `antioxidants.xlsx` — Antioxidant substance families

## Usage

### Run All Example Scenarios
```bash
python survey/run_batch.py survey/examples/batch_scenarios/ -o survey/examples/batch_results/
```

### Import to Survey Simulator
1. Open Survey Simulator: http://127.0.0.1:8001
2. Click "Import JSON" and select `example_pf_jobs.json`
3. Click "Run Simulation"

---
*Part of SFPPy - Scientific Framework for Food Packaging & Migration Modeling*
