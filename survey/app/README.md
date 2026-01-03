# Survey Web Applications

FastAPI-based web applications for the SFPPy Survey module.

## Applications

### Family Editor (`main.py`) - Port 8000
Spreadsheet-style editor for substance families and packaging configurations.

**Features:**
- Create/edit/delete substance families
- PubChem integration for substance validation
- Import/Export XLSX spreadsheets
- Link to Survey Simulator

### Survey Simulator (`simulator.py`) - Port 8001
Interactive batch migration simulator with uncertainty propagation.

**Features:**
- PF Items (Packaged Food Items) management
- Pre-defined substance families
- JSON/XLSX import/export
- Parallel simulation execution
- Interactive PDF/CDF charts with zoom
- Fullscreen mode and bulk export

## Directory Structure

```
app/
├── README.md              # This file
├── main.py                # Family Editor application
├── simulator.py           # Survey Simulator application
├── default_config.yml     # Default simulation configuration
├── templates/             # Jinja2 HTML templates
│   ├── index.html         # Family Editor UI
│   └── simulator.html     # Survey Simulator UI
└── static/                # Static assets (if any)
```

## Running Applications

### Using Launcher (Recommended)
```bash
# Both applications
python survey/launcher.py

# Family Editor only
python survey/launcher.py --app editor

# Survey Simulator only
python survey/launcher.py --app simulator
```

### Direct Uvicorn
```bash
# Family Editor
python -m uvicorn survey.app.main:app --reload --port 8000

# Survey Simulator
python -m uvicorn survey.app.simulator:app --reload --port 8001
```

## API Endpoints

See `../API_REFERENCE.md` for complete API documentation.

### Key Endpoints

**Family Editor (8000):**
- `GET /` — Editor UI
- `GET /api/spreadsheet` — Get current spreadsheet
- `POST /api/spreadsheet/import` — Import XLSX
- `GET /api/spreadsheet/export` — Export XLSX

**Survey Simulator (8001):**
- `GET /simulator` — Simulator UI
- `POST /api/jobs/import-json` — Import PF Jobs
- `POST /api/jobs/run-batch` — Run simulation
- `POST /api/results/export` — Export results

## Configuration

Edit `default_config.yml` to customize:
- Parallel processing (workers)
- Monte Carlo samples
- Output formats
- PDF bins

---
*Part of SFPPy - Scientific Framework for Food Packaging & Migration Modeling*
