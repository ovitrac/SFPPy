# SFPPy Survey Web Applications - API Reference

## Overview

This document describes the REST API endpoints for the SFPPy Survey web applications:

1. **Family Editor** - Port 8000
2. **Survey Simulator** - Port 8001

---

## Family Editor (Port 8000)

Base URL: `http://127.0.0.1:8000`

### Health & System

#### GET /health
Check application health and status.

**Response:**
```json
{
    "status": "ok",
    "family_substances": 0,
    "working_folder": "/home/user/survey/examples/families"
}
```

#### GET /api/system-info
Get system configuration.

**Response:**
```json
{
    "working_folder": "/home/user/survey/examples/families",
    "yaml_files": ["plasticizers.yml", "antioxidants.yml"]
}
```

### Family Management

#### GET /api/families
List all family YAML files in working folder.

**Response:**
```json
{
    "success": true,
    "families": [
        {
            "filename": "plasticizers.yml",
            "name": "plasticizers",
            "description": "Common plasticizers for PVC",
            "polymer": "PVC",
            "n_substances": 3
        }
    ]
}
```

#### GET /api/family/{filename}
Get a specific family definition.

**Parameters:**
- `filename` (path): YAML filename (e.g., "plasticizers.yml")

**Response:**
```json
{
    "success": true,
    "family": {
        "name": "plasticizers",
        "description": "Common plasticizers for PVC",
        "polymer": "PVC",
        "thickness_um": 200,
        "temperature_C": 25,
        "contact_days": 10,
        "substances": {
            "77-90-7": {
                "name": "ATBC",
                "C0_ppm": 150,
                "SML_ppm": 60
            }
        }
    }
}
```

#### POST /api/family
Create or update a family.

**Request Body:**
```json
{
    "filename": "new_family.yml",
    "family": {
        "name": "new_family",
        "description": "Description here",
        "polymer": "LDPE",
        "thickness_um": 100,
        "temperature_C": 25,
        "contact_days": 10,
        "substances": {}
    }
}
```

**Response:**
```json
{
    "success": true,
    "message": "Family saved to new_family.yml"
}
```

#### DELETE /api/family/{filename}
Delete a family file.

**Response:**
```json
{
    "success": true,
    "message": "Family new_family.yml deleted"
}
```

### Substance Management

#### POST /api/family/{filename}/substance
Add a substance to a family.

**Request Body:**
```json
{
    "identifier": "77-90-7",
    "name": "ATBC",
    "C0_ppm": 150,
    "SML_ppm": 60
}
```

**Response:**
```json
{
    "success": true,
    "message": "Substance 77-90-7 added to plasticizers.yml"
}
```

#### DELETE /api/family/{filename}/substance/{identifier}
Remove a substance from a family.

**Response:**
```json
{
    "success": true,
    "message": "Substance 77-90-7 removed from plasticizers.yml"
}
```

### PubChem Integration

#### GET /api/pubchem/search/{query}
Search PubChem for substance information.

**Parameters:**
- `query` (path): CAS number, name, or CID

**Response:**
```json
{
    "success": true,
    "results": [
        {
            "cid": 6505,
            "name": "Acetyl tributyl citrate",
            "cas": "77-90-7",
            "molecular_weight": 402.5,
            "molecular_formula": "C20H34O8"
        }
    ]
}
```

---

## Survey Simulator (Port 8001)

Base URL: `http://127.0.0.1:8001`

### Health & System

#### GET /health
Check application health and system resources.

**Response:**
```json
{
    "status": "ok",
    "cpu_count": 20,
    "default_workers": 10,
    "active_jobs": 0,
    "spreadsheet_support": {
        "xlsx": true,
        "ods": false
    }
}
```

#### GET /api/system-info
Get detailed system information.

**Response:**
```json
{
    "cpu_count": 20,
    "default_workers": 10,
    "spreadsheet_support": {"xlsx": true, "ods": false},
    "upload_dir": "/tmp/sfppy_simulator/uploads",
    "output_dir": "/tmp/sfppy_simulator/output"
}
```

### Spreadsheet Operations

#### GET /api/download-template
Download XLSX template with example data.

**Response:** File download (`sfppy_survey_template.xlsx`)

**Content-Type:** `application/vnd.openxmlformats-officedocument.spreadsheetml.sheet`

#### POST /api/upload-spreadsheet
Upload and parse a spreadsheet file.

**Request:** `multipart/form-data`
- `file`: XLSX or ODS file

**Response:**
```json
{
    "success": true,
    "file_id": "a1b2c3d4",
    "filename": "my_families.xlsx",
    "file_path": "/tmp/sfppy_simulator/uploads/a1b2c3d4_my_families.xlsx",
    "families": 3,
    "substances": 7,
    "configs": [...]
}
```

### Editor Endpoints

#### POST /api/editor/new
Create a new empty spreadsheet in memory.

**Response:**
```json
{
    "success": true,
    "message": "New spreadsheet created",
    "families": 0,
    "substances": 0
}
```

#### POST /api/editor/load-demo
Load demo spreadsheet with example families.

**Response:**
```json
{
    "success": true,
    "message": "Demo spreadsheet loaded with 3 families and 7 substances",
    "families": 3,
    "substances": 7
}
```

**Demo Data Includes:**
| Family | Polymer | Substances |
|--------|---------|------------|
| plasticizers_pvc | PVC | ATBC (77-90-7), DEHA (103-23-1), DINCH (166412-78-8) |
| antioxidants_ldpe | LDPE | Irganox 1076 (2082-79-3), BHT (128-37-0) |
| uv_stabilizers_pp | PP | Tinuvin 326 (3896-11-5), Chimassorb 81 (1843-05-6) |

#### GET /api/editor/current
Get current spreadsheet data.

**Response:**
```json
{
    "success": true,
    "families": [
        {
            "name": "plasticizers_pvc",
            "description": "Common plasticizers for PVC",
            "polymer": "PVC",
            "thickness_um": 200.0,
            "temperature_C": 25.0,
            "contact_days": 10.0,
            "food_volume_ml": 1000.0,
            "food_density": 1.0
        }
    ],
    "substances": [
        {
            "family_name": "plasticizers_pvc",
            "identifier": "77-90-7",
            "C0_min": 50.0,
            "C0_likely": 150.0,
            "C0_max": 300.0,
            "unit": "mg/kg"
        }
    ],
    "configs": [...]
}
```

#### POST /api/editor/add-family
Add a family to the current spreadsheet.

**Request Body:**
```json
{
    "name": "new_family",
    "description": "Optional description",
    "polymer": "LDPE",
    "thickness_um": 100.0,
    "temperature_C": 25.0,
    "contact_days": 10.0,
    "food_volume_ml": 1000.0,
    "food_density": 1.0
}
```

**Response:**
```json
{
    "success": true,
    "message": "Family 'new_family' added",
    "family_count": 4
}
```

#### PUT /api/editor/update-family/{name}
Update an existing family.

**Parameters:**
- `name` (path): Current family name

**Request Body:** Same as `add-family`

**Response:**
```json
{
    "success": true,
    "message": "Family 'new_family' updated"
}
```

#### DELETE /api/editor/delete-family/{name}
Delete a family and all its substances.

**Response:**
```json
{
    "success": true,
    "message": "Family 'new_family' and its substances deleted"
}
```

#### POST /api/editor/add-substance
Add a substance to a family.

**Request Body:**
```json
{
    "family_name": "plasticizers_pvc",
    "identifier": "84-74-2",
    "C0_min": 50.0,
    "C0_likely": 100.0,
    "C0_max": 200.0,
    "unit": "mg/kg"
}
```

**Response:**
```json
{
    "success": true,
    "message": "Substance '84-74-2' added to 'plasticizers_pvc'",
    "substance_count": 4
}
```

#### PUT /api/editor/update-substance
Update an existing substance.

**Request:** `application/x-www-form-urlencoded`
- `family_name`: Family name
- `old_identifier`: Current identifier
- `identifier`: New identifier
- `C0_min`, `C0_likely`, `C0_max`: Concentration values
- `unit`: Unit (default: mg/kg)

**Response:**
```json
{
    "success": true,
    "message": "Substance updated"
}
```

#### DELETE /api/editor/delete-substance
Delete a substance.

**Query Parameters:**
- `family_name`: Family name
- `identifier`: Substance identifier

**Response:**
```json
{
    "success": true,
    "message": "Substance deleted"
}
```

#### POST /api/editor/save-xlsx
Export current data to XLSX file.

**Request:** `application/x-www-form-urlencoded`
- `filename`: Output filename (default: "survey_data.xlsx")

**Response:** File download

#### POST /api/editor/run-batch
Run batch simulation on current spreadsheet data.

**Request:** `application/x-www-form-urlencoded`
- `n_workers`: Number of parallel workers
- `n_samples`: Monte Carlo samples (default: 1000)

**Response:**
```json
{
    "success": true,
    "job_id": "a1b2c3d4",
    "message": "Batch job started with 10 workers"
}
```

#### POST /api/editor/load-from-upload
Load uploaded spreadsheet into editor.

**Request:** `application/x-www-form-urlencoded`
- `file_path`: Path to uploaded file

**Response:**
```json
{
    "success": true,
    "message": "Spreadsheet loaded into editor",
    "families": 3,
    "substances": 7
}
```

### Batch Job Management

#### POST /api/start-batch
Start a batch job from uploaded file.

**Request:** `application/x-www-form-urlencoded`
- `file_path`: Path to spreadsheet file
- `n_workers`: Number of workers
- `n_samples`: Number of samples (default: 1000)

**Response:**
```json
{
    "success": true,
    "job_id": "a1b2c3d4",
    "message": "Batch job started with 10 workers"
}
```

#### GET /api/job-status/{job_id}
Get status of a batch job.

**Response (Running):**
```json
{
    "success": true,
    "job_id": "a1b2c3d4",
    "status": "running",
    "progress": 45.5,
    "total_tasks": 3,
    "completed_tasks": 1,
    "failed_tasks": 0,
    "elapsed_s": 12.5,
    "results": [],
    "n_workers": 10,
    "n_samples": 1000
}
```

**Response (Completed):**
```json
{
    "success": true,
    "job_id": "a1b2c3d4",
    "status": "completed",
    "progress": 100.0,
    "total_tasks": 3,
    "completed_tasks": 3,
    "failed_tasks": 0,
    "elapsed_s": 35.2,
    "results": [
        {
            "task_id": "task_0000",
            "family_name": "plasticizers_pvc",
            "success": true,
            "duration_s": 11.5,
            "n_samples": 1000,
            "quantiles": {
                "q50": 0.0091,
                "q75": 0.0145,
                "q90": 0.0178,
                "q95": 0.0190,
                "q99": 0.0250,
                "mean": 0.0098,
                "std": 0.0052,
                "min": 0.0001,
                "max": 0.0307
            },
            "centers": [0.0, 0.001, 0.002, ...],
            "pdf": [0.0, 0.12, 0.25, ...],
            "cdf": [0.0, 0.05, 0.15, ...]
        }
    ]
}
```

#### GET /api/jobs
List all batch jobs.

**Response:**
```json
{
    "success": true,
    "jobs": [
        {
            "job_id": "a1b2c3d4",
            "status": "completed",
            "spreadsheet_name": "Editor Data",
            "created_at": "2025-12-18T17:43:19.543840",
            "n_workers": 10,
            "n_samples": 1000
        }
    ]
}
```

#### DELETE /api/job/{job_id}
Delete a job and its output files.

**Response:**
```json
{
    "success": true,
    "message": "Job a1b2c3d4 deleted"
}
```

### Results

#### GET /api/results/{job_id}/summary
Get summary statistics for completed job.

**Response:**
```json
{
    "success": true,
    "job_id": "a1b2c3d4",
    "total_families": 3,
    "successful": 3,
    "failed": 0,
    "total_duration_s": 35.2,
    "families": [
        {
            "name": "plasticizers_pvc",
            "n_samples": 1000,
            "duration_s": 11.5,
            "q50": 0.0091,
            "q95": 0.0190,
            "q99": 0.0250,
            "mean": 0.0098,
            "max": 0.0307
        }
    ]
}
```

#### GET /api/results/{job_id}/distributions
Get distribution data for visualization.

**Response:**
```json
{
    "success": true,
    "job_id": "a1b2c3d4",
    "distributions": [
        {
            "family_name": "plasticizers_pvc",
            "centers": [0.0, 0.001, 0.002, ...],
            "pdf": [0.0, 0.12, 0.25, ...],
            "cdf": [0.0, 0.05, 0.15, ...],
            "quantiles": {"q50": 0.0091, "q95": 0.0190, ...}
        }
    ]
}
```

### Configuration

#### GET /api/config
Get current configuration.

**Response:**
```json
{
    "success": true,
    "config": {
        "parallel": {"workers": null, "max_workers": 20},
        "sampling": {
            "n_samples": 1000,
            "presets": {
                "quick": 100,
                "fast": 500,
                "standard": 1000,
                "precise": 5000,
                "high_precision": 10000
            }
        },
        "physics": {
            "default_polymer": "LDPE",
            "default_thickness_um": 100,
            "default_temperature_C": 25,
            "default_contact_days": 10
        },
        "solver": {"pdf_bins": 250},
        "output": {"plot_dpi": 150, "quantiles": [0.50, 0.75, 0.90, 0.95, 0.99]},
        "polymers": ["LDPE", "HDPE", "PP", "PET", "PS", "PVC", "PA", "EVOH", "generic"],
        "food_simulants": ["ethanol10", "ethanol20", "ethanol50", "ethanol95", "oliveoil", "water", "acetic3"]
    }
}
```

#### GET /api/config/yaml
Get configuration as YAML string.

**Response:**
```json
{
    "success": true,
    "yaml": "parallel:\n  workers: null\n  max_workers: 20\n..."
}
```

#### POST /api/config
Update configuration from YAML.

**Request:** `application/x-www-form-urlencoded`
- `config_yaml`: YAML configuration string

**Response:**
```json
{
    "success": true,
    "message": "Configuration updated",
    "config": {...}
}
```

#### POST /api/config/reset
Reset configuration to defaults.

**Response:**
```json
{
    "success": true,
    "message": "Configuration reset to defaults",
    "config": {...}
}
```

#### GET /api/config/defaults
Get default configuration.

**Response:**
```json
{
    "success": true,
    "config": {...},
    "yaml": "..."
}
```

---

## Error Responses

All endpoints return errors in this format:

```json
{
    "success": false,
    "error": "Error message describing the issue"
}
```

**HTTP Status Codes:**
- `200` - Success
- `400` - Bad Request (invalid input)
- `404` - Not Found (resource doesn't exist)
- `500` - Internal Server Error

---

## Data Models

### Family

```typescript
interface Family {
    name: string;
    description?: string;
    polymer: string;          // LDPE, HDPE, PP, PET, PS, PVC, PA, EVOH, generic
    thickness_um: number;     // Layer thickness in µm
    temperature_C: number;    // Contact temperature in °C
    contact_days: number;     // Contact duration in days
    food_volume_ml?: number;  // Food volume in mL (default: 1000)
    food_density?: number;    // Food density in kg/L (default: 1.0)
}
```

### Substance

```typescript
interface Substance {
    family_name: string;      // Parent family name
    identifier: string;       // CAS number, PubChem CID, or name
    C0_min: number;          // Minimum initial concentration (mg/kg)
    C0_likely: number;       // Most likely concentration (mg/kg)
    C0_max: number;          // Maximum concentration (mg/kg)
    unit?: string;           // Unit (default: "mg/kg")
}
```

### SimulationResult

```typescript
interface SimulationResult {
    task_id: string;
    family_name: string;
    success: boolean;
    error?: string;
    duration_s: number;
    n_samples: number;
    quantiles: {
        q50: number;
        q75: number;
        q90: number;
        q95: number;
        q99: number;
        mean: number;
        std: number;
        min: number;
        max: number;
    };
    centers: number[];        // PDF bin centers
    pdf: number[];           // PDF values
    cdf: number[];           // CDF values
}
```

---

## Usage Examples

### cURL Examples

```bash
# Check health
curl http://127.0.0.1:8001/health

# Load demo data
curl -X POST http://127.0.0.1:8001/api/editor/load-demo

# Get current spreadsheet
curl http://127.0.0.1:8001/api/editor/current

# Add a family
curl -X POST http://127.0.0.1:8001/api/editor/add-family \
  -H "Content-Type: application/json" \
  -d '{"name":"test","polymer":"LDPE","thickness_um":100,"temperature_C":25,"contact_days":10}'

# Add a substance
curl -X POST http://127.0.0.1:8001/api/editor/add-substance \
  -H "Content-Type: application/json" \
  -d '{"family_name":"test","identifier":"77-90-7","C0_min":50,"C0_likely":100,"C0_max":200}'

# Run simulation
curl -X POST http://127.0.0.1:8001/api/editor/run-batch \
  -d "n_workers=4&n_samples=500"

# Check job status
curl http://127.0.0.1:8001/api/job-status/{job_id}

# Get results
curl http://127.0.0.1:8001/api/results/{job_id}/summary
```

### Python Examples

```python
import requests

BASE_URL = "http://127.0.0.1:8001"

# Load demo
r = requests.post(f"{BASE_URL}/api/editor/load-demo")
print(r.json())

# Run simulation
r = requests.post(f"{BASE_URL}/api/editor/run-batch",
                  data={"n_workers": 4, "n_samples": 1000})
job_id = r.json()["job_id"]

# Poll for completion
import time
while True:
    r = requests.get(f"{BASE_URL}/api/job-status/{job_id}")
    status = r.json()
    print(f"Progress: {status['progress']:.1f}%")
    if status["status"] in ["completed", "failed"]:
        break
    time.sleep(1)

# Get results
r = requests.get(f"{BASE_URL}/api/results/{job_id}/summary")
for fam in r.json()["families"]:
    print(f"{fam['name']}: Q95={fam['q95']:.4f} mg/kg")
```

---

## Authors

**Olivier Vitrac, PhD, HDR**
- Email: olivier.vitrac@gmail.com
- Affiliation: INRAE / Generative Simulation
