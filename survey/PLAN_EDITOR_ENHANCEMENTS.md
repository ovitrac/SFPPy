# Plan of Work: Survey Simulator Editor Enhancements

**Date:** 2025-12-19
**Author:** Claude Code (for Olivier Vitrac)
**Status:** Draft for Approval (v2 - Collapsible UI)

---

## 1. Objectives

Enhance the **SFPPy Survey Simulator** editor with:

1. **Collapsible rows** for each packaged food (keep full batch view)
2. **YAML syntax alignment** (matching `basic_monolayer.yml` structure)
3. **XLS import/export** for all jobs
4. **Substance family validation** (ensure families exist)
5. **Triangular priors** for contact time and concentration

---

## 2. YAML Structure Reference

From `survey/examples/basic_monolayer.yml`:

```yaml
name: basic_monolayer_example

physics:
  monolayer:
    polymer: PP
    thickness_m: 50e-6           # 50 µm
    temperature_degC: 40.0

  interface:
    h_m_s: 1e-7                  # Mass transfer coefficient
    surface_area_m2: 0.06        # 6 dm² (typical tray)
    food_volume_m3: 0.001        # 1 L
    contact_temperature_degC: 40.0
    cf0: 0.0
    food_simulant: oliveoil      # EU Simulant D

family:
  masses_g_mol: [100.0, 150.0, 200.0, 250.0]
  # OR:
  # substances:
  #   - name: limonene
  #   - name: benzaldehyde

priors:
  time_s:
    triangular:
      min: 0.0
      mode: 604800.0             # 7 days
      max: 2592000.0             # 30 days
    grid:
      nlow: 15
      nhigh: 15

  cp0_av:
    triangular:
      min: 0.0
      mode: 10.0
      max: 50.0
    grid:
      nlow: 15
      nhigh: 15
```

---

## 3. Collapsible UI Design

### 3.1 Master View (All Packaged Foods)

A compact table showing all jobs at a glance:

```
┌──────────────────────────────────────────────────────────────────────────────────────────┐
│ 📦 BATCH JOBS                                                    [+ Add Food] [📥 Import XLS] │
├──────────────────────────────────────────────────────────────────────────────────────────┤
│ ▶ │ Name              │ Polymer │ Simulant   │ Time (mode) │ Substances │ Status  │ 🗑️  │
├───┼───────────────────┼─────────┼────────────┼─────────────┼────────────┼─────────┼─────┤
│ ▼ │ yogurt_cup        │ PP      │ ethanol50  │ 30 days     │ 3          │ ✓ Valid │ 🗑️  │
│   └─────────────────────────────────────────────────────────────────────────────────────┤
│   │  [Expanded details panel - see section 3.2]                                         │
│   └─────────────────────────────────────────────────────────────────────────────────────┤
│ ▶ │ olive_oil_bottle  │ PET     │ oliveoil   │ 90 days     │ 5          │ ✓ Valid │ 🗑️  │
│ ▶ │ water_bottle      │ gPET    │ water      │ 180 days    │ 2          │ ⚠ Check │ 🗑️  │
│ ▶ │ meat_tray         │ PP      │ oliveoil   │ 7 days      │ 4          │ ✓ Valid │ 🗑️  │
└──────────────────────────────────────────────────────────────────────────────────────────┘
                                                              [📤 Export XLS] [▶ Run Batch]
```

**Columns:**
| Column | Description |
|--------|-------------|
| ▶/▼ | Expand/collapse toggle |
| Name | Job/food identifier |
| Polymer | Primary layer polymer |
| Simulant | Food simulant |
| Time (mode) | Most likely contact time |
| Substances | Number of substances in family |
| Status | Validation status |
| 🗑️ | Delete action |

### 3.2 Expanded Detail Panel (Per Food)

When a row is expanded, show editable details grouped by YAML sections:

```
┌─────────────────────────────────────────────────────────────────────────────────────────┐
│ ▼ yogurt_cup                                                                    [Save] │
├─────────────────────────────────────────────────────────────────────────────────────────┤
│                                                                                          │
│ ┌─ physics.monolayer ──────────────────────────────────────────────────────────────────┐ │
│ │ Polymer: [PP      ▼]   Thickness: [50   ] [µm▼]   Temperature: [25   ] °C           │ │
│ └──────────────────────────────────────────────────────────────────────────────────────┘ │
│                                                                                          │
│ ┌─ physics.interface ──────────────────────────────────────────────────────────────────┐ │
│ │ Simulant: [ethanol50 ▼]   h: [1e-7  ] m/s   cf0: [0     ]                           │ │
│ │ Surface Area: [600  ] [cm²▼]   Food Volume: [125  ] [mL ▼]                          │ │
│ └──────────────────────────────────────────────────────────────────────────────────────┘ │
│                                                                                          │
│ ┌─ priors.time_s (triangular) ─────────────────────────────────────────────────────────┐ │
│ │ Min: [0     ] [days▼]   Mode: [30    ] [days▼]   Max: [90    ] [days▼]              │ │
│ └──────────────────────────────────────────────────────────────────────────────────────┘ │
│                                                                                          │
│ ┌─ priors.cp0_av (triangular) ─────────────────────────────────────────────────────────┐ │
│ │ Min: [0     ] mg/kg   Mode: [10    ] mg/kg   Max: [50    ] mg/kg                    │ │
│ └──────────────────────────────────────────────────────────────────────────────────────┘ │
│                                                                                          │
│ ┌─ family.substances ──────────────────────────────────────────────────────────────────┐ │
│ │ [✓] ATBC (77-90-7)           C0: min=0, mode=20, max=100 mg/kg                       │ │
│ │ [✓] DEHA (103-23-1)          C0: min=0, mode=15, max=80 mg/kg                        │ │
│ │ [✓] Irganox 1076 (2082-79-3) C0: min=0, mode=5, max=30 mg/kg                         │ │
│ │ [+ Add Substance]                                                                     │ │
│ └──────────────────────────────────────────────────────────────────────────────────────┘ │
│                                                                                          │
└─────────────────────────────────────────────────────────────────────────────────────────┘
```

### 3.3 CSS Classes for Collapsible Rows

```css
/* Collapsible food rows */
.food-row {
    cursor: pointer;
    transition: background-color 0.2s;
}
.food-row:hover {
    background-color: rgba(34, 197, 94, 0.1);
}
.food-row.expanded {
    background-color: rgba(34, 197, 94, 0.05);
}

/* Detail panel */
.food-details {
    display: none;
    padding: 1rem;
    border-left: 3px solid #22c55e;
    background: linear-gradient(to right, rgba(34, 197, 94, 0.02), transparent);
}
.food-details.visible {
    display: block;
}

/* Section groups matching YAML */
.yaml-section {
    background: #1a1a1a;
    border: 1px solid #333;
    border-radius: 6px;
    padding: 0.75rem;
    margin-bottom: 0.5rem;
}
.yaml-section-header {
    font-family: 'Fira Code', monospace;
    color: #22c55e;
    font-size: 0.8rem;
    margin-bottom: 0.5rem;
}
```

---

## 4. Updated Data Models

### 4.1 Backend Models (Pydantic)

```python
from pydantic import BaseModel, validator
from typing import List, Optional
from enum import Enum

# Validated enums
class PolymerType(str, Enum):
    LDPE = "LDPE"
    HDPE = "HDPE"
    LLDPE = "LLDPE"
    PP = "PP"
    oPP = "oPP"
    PET = "PET"
    gPET = "gPET"
    rPET = "rPET"
    wPET = "wPET"
    PS = "PS"
    HIPS = "HIPS"
    pPVC = "pPVC"
    PVC = "PVC"
    PA6 = "PA6"
    PA66 = "PA6,6"

class SimulantType(str, Enum):
    oliveoil = "oliveoil"
    ethanol50 = "ethanol50"
    ethanol95 = "ethanol95"
    water = "water"
    water3aceticacid = "water3aceticacid"
    isooctane = "isooctane"
    tenax = "tenax"

class TimeUnit(str, Enum):
    s = "s"
    min = "min"
    h = "h"
    days = "days"
    months = "months"
    years = "years"


class TriangularPrior(BaseModel):
    """Triangular distribution specification."""
    min_val: float = 0.0
    mode: float
    max_val: float
    unit: str = ""  # For time: days, for conc: mg/kg

    @validator('max_val')
    def max_gte_mode(cls, v, values):
        if 'mode' in values and v < values['mode']:
            raise ValueError('max must be >= mode')
        return v

    @validator('mode')
    def mode_gte_min(cls, v, values):
        if 'min_val' in values and v < values['min_val']:
            raise ValueError('mode must be >= min')
        return v


class SubstanceSpec(BaseModel):
    """Substance with concentration prior."""
    identifier: str           # CAS or name
    c0_min: float = 0.0
    c0_mode: float = 10.0
    c0_max: float = 50.0
    unit: str = "mg/kg"
    validated: bool = False   # True if found in PubChem


class PackagedFoodJob(BaseModel):
    """Complete job specification aligned with YAML structure."""

    # Identity
    name: str
    description: str = ""

    # physics.monolayer
    polymer: PolymerType = PolymerType.LDPE
    thickness_value: float = 100.0
    thickness_unit: str = "µm"
    layer_temperature_C: float = 40.0

    # physics.interface
    simulant: SimulantType = SimulantType.oliveoil
    h_m_s: float = 1e-7
    surface_area_value: float = 600.0
    surface_area_unit: str = "cm²"
    food_volume_value: float = 1000.0
    food_volume_unit: str = "mL"
    contact_temperature_C: float = 40.0
    cf0: float = 0.0

    # priors.time_s
    time_prior: TriangularPrior

    # priors.cp0_av
    conc_prior: TriangularPrior

    # family.substances
    substances: List[SubstanceSpec] = []

    # Validation state
    is_valid: bool = True
    validation_errors: List[str] = []


class BatchJobsInput(BaseModel):
    """Collection of jobs for batch processing."""
    jobs: List[PackagedFoodJob]
```

### 4.2 Updated Spreadsheet Format

**Sheet: Families** (enhanced)

| Column | Type | YAML Path | Description |
|--------|------|-----------|-------------|
| name | str | name | Job identifier |
| description | str | - | Optional description |
| polymer | enum | physics.monolayer.polymer | From PolymerType |
| thickness_um | float | physics.monolayer.thickness_m | In µm |
| layer_temp_C | float | physics.monolayer.temperature_degC | Layer temp |
| simulant | enum | physics.interface.food_simulant | From SimulantType |
| h_m_s | float | physics.interface.h_m_s | Mass transfer coef |
| surface_area_cm2 | float | physics.interface.surface_area_m2 | In cm² |
| food_volume_ml | float | physics.interface.food_volume_m3 | In mL |
| contact_temp_C | float | physics.interface.contact_temperature_degC | Contact temp |
| time_min_days | float | priors.time_s.triangular.min | In days |
| time_mode_days | float | priors.time_s.triangular.mode | In days |
| time_max_days | float | priors.time_s.triangular.max | In days |

**Sheet: Substances** (unchanged + validation)

| Column | Type | Description |
|--------|------|-------------|
| family_name | str | Links to Families.name |
| identifier | str | CAS or substance name |
| C0_min | float | Minimum concentration |
| C0_likely | float | Most likely concentration |
| C0_max | float | Maximum concentration |
| unit | str | Concentration unit |
| **validated** | bool | NEW: True if found in PubChem |

---

## 5. Validation Requirements

### 5.1 Polymer Validation

```python
VALID_POLYMERS = set(Dpiringer.piringer_data.keys())

def validate_polymer(polymer: str) -> tuple[bool, str]:
    """Check polymer exists in Dpiringer model."""
    if polymer in VALID_POLYMERS:
        return True, ""
    return False, f"Unknown polymer: {polymer}. Valid: {', '.join(sorted(VALID_POLYMERS))}"
```

### 5.2 Simulant Validation

```python
VALID_SIMULANTS = {
    "oliveoil", "oil", "ethanol", "ethanol95", "ethanol50",
    "water", "water3aceticacid", "isooctane", "tenax",
    "acetonitrile", "methanol"
}

def validate_simulant(simulant: str) -> tuple[bool, str]:
    """Check simulant exists in food.py."""
    if simulant.lower() in VALID_SIMULANTS:
        return True, ""
    return False, f"Unknown simulant: {simulant}. Valid: {', '.join(sorted(VALID_SIMULANTS))}"
```

### 5.3 Substance Family Validation

```python
async def validate_substance(identifier: str) -> tuple[bool, dict]:
    """
    Validate substance exists and retrieve properties.

    Returns:
        (valid, info_dict) where info_dict contains:
        - name, cas, M (molar mass), formula
    """
    try:
        from patankar.loadpubchem import migrant
        m = migrant(identifier)
        return True, {
            "name": m.name,
            "cas": m.CAS,
            "M": m.M,
            "formula": m.formula,
        }
    except Exception as e:
        return False, {"error": str(e)}


def validate_family_substances(family_name: str, substances: List[SubstanceSpec]) -> List[str]:
    """
    Validate all substances in a family.

    Returns list of error messages (empty if all valid).
    """
    errors = []
    if not substances:
        errors.append(f"Family '{family_name}' has no substances")

    for sub in substances:
        valid, info = validate_substance(sub.identifier)
        if not valid:
            errors.append(f"Substance '{sub.identifier}' not found: {info.get('error', 'Unknown')}")
        else:
            sub.validated = True

    return errors
```

---

## 6. XLS Import/Export Specification

### 6.1 Import Flow

```
[📥 Import XLS]
      │
      ▼
┌─────────────────────────────────────┐
│  1. Parse Families sheet            │
│  2. Parse Substances sheet          │
│  3. Validate polymers               │
│  4. Validate simulants              │
│  5. Validate substances (PubChem)   │
│  6. Show validation report          │
│  7. Add valid jobs to batch         │
└─────────────────────────────────────┘
```

### 6.2 Export Flow

```
[📤 Export XLS]
      │
      ▼
┌─────────────────────────────────────┐
│  1. Collect all PackagedFoodJob     │
│  2. Convert to FamilyRow format     │
│  3. Convert to SubstanceRow format  │
│  4. Write to XLSX via write_xlsx()  │
│  5. Download file                   │
└─────────────────────────────────────┘
```

### 6.3 API Endpoints

```python
@app.post("/api/editor/import-xlsx")
async def import_xlsx(file: UploadFile = File(...)) -> JSONResponse:
    """
    Import jobs from XLSX spreadsheet.

    Returns:
        {
            "success": true,
            "jobs_imported": 5,
            "jobs_skipped": 1,
            "validation_errors": [
                {"family": "unknown_polymer", "error": "Unknown polymer: XXX"}
            ]
        }
    """
    pass


@app.get("/api/editor/export-xlsx")
async def export_xlsx() -> FileResponse:
    """
    Export all current jobs to XLSX file.

    Returns downloadable XLSX file.
    """
    pass


@app.post("/api/editor/validate-substance")
async def validate_substance_endpoint(identifier: str) -> JSONResponse:
    """
    Validate a single substance identifier.

    Returns:
        {
            "valid": true,
            "name": "Limonene",
            "cas": "138-86-3",
            "M": 136.23,
            "formula": "C10H16"
        }
    """
    pass


@app.post("/api/editor/validate-job")
async def validate_job(job: PackagedFoodJob) -> JSONResponse:
    """
    Validate a complete job specification.

    Returns:
        {
            "valid": true,
            "errors": [],
            "warnings": ["Substance X not found in cache, will query PubChem"]
        }
    """
    pass
```

---

## 7. JavaScript for Collapsible Rows

```javascript
// Toggle food row expansion
function toggleFoodRow(jobId) {
    const row = document.querySelector(`[data-job-id="${jobId}"]`);
    const details = document.querySelector(`[data-job-details="${jobId}"]`);
    const arrow = row.querySelector('.expand-arrow');

    row.classList.toggle('expanded');
    details.classList.toggle('visible');
    arrow.textContent = row.classList.contains('expanded') ? '▼' : '▶';
}

// Render all jobs in master table
function renderJobsTable(jobs) {
    const tbody = document.getElementById('jobs-table-body');
    tbody.innerHTML = '';

    jobs.forEach(job => {
        // Summary row
        const row = document.createElement('tr');
        row.className = 'food-row';
        row.dataset.jobId = job.name;
        row.onclick = () => toggleFoodRow(job.name);
        row.innerHTML = `
            <td class="expand-arrow">▶</td>
            <td>${job.name}</td>
            <td>${job.polymer}</td>
            <td>${job.simulant}</td>
            <td>${job.time_prior.mode} ${job.time_prior.unit}</td>
            <td>${job.substances.length}</td>
            <td>${job.is_valid ? '✓ Valid' : '⚠ Check'}</td>
            <td><button onclick="deleteJob('${job.name}'); event.stopPropagation();">🗑️</button></td>
        `;
        tbody.appendChild(row);

        // Details row (hidden by default)
        const detailsRow = document.createElement('tr');
        detailsRow.innerHTML = `
            <td colspan="8">
                <div class="food-details" data-job-details="${job.name}">
                    ${renderJobDetails(job)}
                </div>
            </td>
        `;
        tbody.appendChild(detailsRow);
    });
}

// Render expanded job details matching YAML sections
function renderJobDetails(job) {
    return `
        <div class="yaml-section">
            <div class="yaml-section-header">physics.monolayer</div>
            <div class="section-content">
                <label>Polymer:
                    <select onchange="updateJob('${job.name}', 'polymer', this.value)">
                        ${POLYMER_OPTIONS.map(p =>
                            `<option value="${p}" ${job.polymer === p ? 'selected' : ''}>${p}</option>`
                        ).join('')}
                    </select>
                </label>
                <label>Thickness:
                    <input type="number" value="${job.thickness_value}"
                           onchange="updateJob('${job.name}', 'thickness_value', this.value)">
                    <select onchange="updateJob('${job.name}', 'thickness_unit', this.value)">
                        <option value="µm" ${job.thickness_unit === 'µm' ? 'selected' : ''}>µm</option>
                        <option value="mm" ${job.thickness_unit === 'mm' ? 'selected' : ''}>mm</option>
                    </select>
                </label>
                <label>Temperature:
                    <input type="number" value="${job.layer_temperature_C}"
                           onchange="updateJob('${job.name}', 'layer_temperature_C', this.value)"> °C
                </label>
            </div>
        </div>

        <div class="yaml-section">
            <div class="yaml-section-header">physics.interface</div>
            <div class="section-content">
                <label>Simulant:
                    <select onchange="updateJob('${job.name}', 'simulant', this.value)">
                        ${SIMULANT_OPTIONS.map(s =>
                            `<option value="${s[0]}" ${job.simulant === s[0] ? 'selected' : ''}>${s[1]}</option>`
                        ).join('')}
                    </select>
                </label>
                <label>Surface Area:
                    <input type="number" value="${job.surface_area_value}"
                           onchange="updateJob('${job.name}', 'surface_area_value', this.value)">
                    <select onchange="updateJob('${job.name}', 'surface_area_unit', this.value)">
                        <option value="cm²">cm²</option>
                        <option value="dm²">dm²</option>
                        <option value="m²">m²</option>
                    </select>
                </label>
                <label>Food Volume:
                    <input type="number" value="${job.food_volume_value}"
                           onchange="updateJob('${job.name}', 'food_volume_value', this.value)">
                    <select onchange="updateJob('${job.name}', 'food_volume_unit', this.value)">
                        <option value="mL">mL</option>
                        <option value="L">L</option>
                    </select>
                </label>
            </div>
        </div>

        <div class="yaml-section">
            <div class="yaml-section-header">priors.time_s (triangular)</div>
            <div class="section-content">
                <label>Min: <input type="number" value="${job.time_prior.min_val}"></label>
                <label>Mode: <input type="number" value="${job.time_prior.mode}"></label>
                <label>Max: <input type="number" value="${job.time_prior.max_val}"></label>
                <select>
                    <option value="days" selected>days</option>
                    <option value="months">months</option>
                    <option value="years">years</option>
                </select>
            </div>
        </div>

        <div class="yaml-section">
            <div class="yaml-section-header">priors.cp0_av (triangular)</div>
            <div class="section-content">
                <label>Min: <input type="number" value="${job.conc_prior.min_val}"></label>
                <label>Mode: <input type="number" value="${job.conc_prior.mode}"></label>
                <label>Max: <input type="number" value="${job.conc_prior.max_val}"></label>
                <span>mg/kg</span>
            </div>
        </div>

        <div class="yaml-section">
            <div class="yaml-section-header">family.substances</div>
            <div class="section-content">
                ${job.substances.map(sub => `
                    <div class="substance-row ${sub.validated ? 'validated' : 'unvalidated'}">
                        <span class="status">${sub.validated ? '✓' : '?'}</span>
                        <span class="identifier">${sub.identifier}</span>
                        <span class="c0">C0: ${sub.c0_min}/${sub.c0_mode}/${sub.c0_max} ${sub.unit}</span>
                        <button onclick="removeSubstance('${job.name}', '${sub.identifier}')">🗑️</button>
                    </div>
                `).join('')}
                <button onclick="addSubstance('${job.name}')">+ Add Substance</button>
            </div>
        </div>
    `;
}
```

---

## 8. Implementation Tasks (Updated)

### Phase 1: Backend Updates

| Task | Priority | Description |
|------|----------|-------------|
| 1.1 | HIGH | Create `PackagedFoodJob` model with YAML alignment |
| 1.2 | HIGH | Add validation functions for polymer/simulant/substance |
| 1.3 | HIGH | Update `FamilyRow` with new fields |
| 1.4 | MEDIUM | Add triangular prior support to batch.py |

### Phase 2: API Endpoints

| Task | Priority | Description |
|------|----------|-------------|
| 2.1 | HIGH | `/api/editor/import-xlsx` with validation |
| 2.2 | HIGH | `/api/editor/export-xlsx` |
| 2.3 | HIGH | `/api/editor/validate-substance` |
| 2.4 | MEDIUM | `/api/editor/validate-job` |

### Phase 3: Frontend - Collapsible UI

| Task | Priority | Description |
|------|----------|-------------|
| 3.1 | HIGH | Master table with expand/collapse |
| 3.2 | HIGH | Detail panel with YAML section styling |
| 3.3 | HIGH | Polymer/simulant dropdowns |
| 3.4 | MEDIUM | Unit conversion on display |

### Phase 4: Validation & Testing

| Task | Priority | Description |
|------|----------|-------------|
| 4.1 | HIGH | Real-time substance validation (PubChem) |
| 4.2 | HIGH | Visual feedback for validation errors |
| 4.3 | MEDIUM | Test with existing YAML examples |

---

## 9. Approval Checklist

Please confirm:

- [ ] **Collapsible design**: Does the master table + expandable details meet your needs?
- [ ] **YAML alignment**: Are the section names (`physics.monolayer`, etc.) appropriate?
- [ ] **Validation**: Should invalid jobs be blocked or just flagged?
- [ ] **XLS format**: Is the updated spreadsheet format (with triangular priors) acceptable?
- [ ] **Substance validation**: Should we cache PubChem results locally?

---

*Plan v2 - Collapsible UI with YAML alignment*
