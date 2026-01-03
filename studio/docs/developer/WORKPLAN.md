# WORKPLAN: Session Save/Load and Script Conversion

## Overview

This workplan defines the design for:
1. **Session file format** - JSON schema to encode complete simulation configurations
2. **Save/Load functionality** - Persist and restore UI state
3. **Script-to-session conversion** - Convert existing Python examples to loadable sessions
4. **Simulation restart** - Continue from previous state with arbitrary CF values

---

## 1. Session File Format Specification

### 1.1 File Extension and Structure

- **Extension**: `.sfppy.json` (e.g., `example1.sfppy.json`)
- **Format**: JSON with strict schema validation
- **Version**: Include schema version for forward compatibility

### 1.2 JSON Schema (v1.0)

```json
{
  "$schema": "https://sfppy.io/session-schema/v1.0.json",
  "version": "1.0",
  "metadata": {
    "name": "string",
    "description": "string",
    "created": "ISO8601 datetime",
    "modified": "ISO8601 datetime",
    "author": "string (optional)",
    "source_script": "string (optional, e.g., example1.py)"
  },
  "geometry": { ... },
  "substances": [ ... ],
  "layers": [ ... ],
  "food": { ... },
  "contact_steps": [ ... ],
  "simulation_config": { ... },
  "results": { ... },
  "restart_state": { ... }
}
```

### 1.3 Detailed Schema Components

#### 1.3.1 Geometry

```json
"geometry": {
  "shape": "cylinder | bottle | box_container | pouch | sphere",
  "dimensions": {
    // Shape-specific, all with explicit units
    "radius": { "value": 50, "unit": "mm" },
    "height": { "value": 100, "unit": "mm" },
    // For bottle:
    "body_radius": { "value": 40, "unit": "mm" },
    "body_height": { "value": 200, "unit": "mm" },
    "neck_radius": { "value": 18, "unit": "mm" },
    "neck_height": { "value": 50, "unit": "mm" },
    // For box_container:
    "length": { "value": 190, "unit": "mm" },
    "width": { "value": 100, "unit": "mm" }
  },
  "computed": {
    "volume_m3": 0.000785,
    "surface_m2": 0.0471
  }
}
```

#### 1.3.2 Substances

```json
"substances": [
  {
    "id": "uuid-or-cid",
    "source": "pubchem | manual | internal",
    "lookup_name": "irganox 1076",  // For PubChem lookup
    "cid": 18725,                    // PubChem CID (if available)
    "properties": {
      "name": "Octadecyl 3-(3,5-di-tert-butyl-4-hydroxyphenyl)propionate",
      "cas": ["2082-79-3"],
      "mw": 530.9,
      "logP": 13.8,
      "formula": "C35H62O3"
    },
    "regulatory": {
      "EU": { "authorized": true, "SML": 6.0 },
      "US": { "fcs_listed": true },
      "CN": { "gb_listed": true }
    },
    "SML": 6.0,
    "color": "#3B82F6"  // For plotting
  }
]
```

#### 1.3.3 Layers (Assembly)

```json
"layers": [
  {
    "index": 1,                        // 1 = contact layer (food side)
    "polymer": "LDPE | PP | gPET | wPET | ...",
    "thickness": { "value": 100, "unit": "um" },
    "name": "Contact layer (optional)",
    "substances": [
      {
        "substance_id": "uuid-or-cid",
        "C0": { "value": 1000, "unit": "mg/kg" },
        "D_auto": true,                // Auto-compute diffusivity
        "D_override": null,            // Manual override (m2/s)
        "k_auto": true,                // Auto-compute partition
        "k_override": null             // Manual override
      }
    ]
  },
  {
    "index": 2,
    "polymer": "PP",
    "thickness": { "value": 300, "unit": "um" },
    "substances": [
      {
        "substance_id": "toluene-cid",
        "C0": { "value": 200, "unit": "mg/kg" }
      }
    ]
  }
]
```

#### 1.3.4 Food / Contact Medium

```json
"food": {
  "category": "realfood | simulant | setoff",
  "texture": "liquid | semisolid | solid",
  "affinity": "fat | aqueous | acidic | alcoholic",
  "simulant": "ethanol | ethanol95 | olive_oil | water | isooctane | null",
  "name": "sandwich (optional)",
  "k0_values": {
    "substance_id_1": { "value": 1.0, "auto": true },
    "substance_id_2": { "value": 0.5, "auto": false }
  }
}
```

#### 1.3.5 Contact Steps (Multi-step scenarios)

```json
"contact_steps": [
  {
    "index": 1,
    "name": "Cold storage",
    "type": "storage | hotfilling | setoff | ambient",
    "temperature": { "value": 7, "unit": "degC" },
    "duration": { "value": 10, "unit": "days" },
    "with_food_contact": true,
    "food_override": null  // Optional: override food properties for this step
  },
  {
    "index": 2,
    "name": "Warming before consumption",
    "type": "storage",
    "temperature": { "value": 25, "unit": "degC" },
    "duration": { "value": 4, "unit": "hours" },
    "with_food_contact": true
  }
]
```

#### 1.3.6 Simulation Configuration

```json
"simulation_config": {
  "solver": "senspatankar",
  "time_factor": 2.0,           // Simulate for 2*tcontact
  "n_time_points": 1000,
  "n_space_points": 200,
  "output_units": {
    "time": "days",
    "length": "um",
    "concentration": "mg/kg"
  }
}
```

#### 1.3.7 Results (Optional - for completed simulations)

```json
"results": {
  "computed_at": "ISO8601 datetime",
  "elapsed_seconds": 2.5,
  "substances": [
    {
      "substance_id": "...",
      "CF_at_tcontact": 5.23,
      "CF_equilibrium": 8.45,
      "SML": 6.0,
      "compliant": true,
      "margin_percent": 12.8
    }
  ],
  "time_series": {
    "t_days": [0, 0.1, ...],
    "CF_by_substance": {
      "substance_id_1": [0, 0.5, ...],
      "substance_id_2": [0, 0.3, ...]
    }
  },
  "concentration_profiles": {
    "x_um": [0, 1, 2, ...],
    "times_days": [0.1, 1, 3, ...],
    "Cx_by_substance": {
      "substance_id_1": [[...], [...], ...]
    }
  }
}
```

#### 1.3.8 Restart State (For continuing simulations)

```json
"restart_state": {
  "enabled": false,
  "from_step": 1,                    // Which contact step to restart from
  "initial_CF": {                    // Arbitrary CF values for restart
    "substance_id_1": 2.5,           // mg/kg in food at restart
    "substance_id_2": 1.0
  },
  "initial_Cx": {                    // Optional: layer concentration profiles
    "substance_id_1": {
      "layer_1": [1000, 999, 998, ...],  // Concentration at each x point
      "layer_2": [0, 0.1, 0.2, ...]
    }
  },
  "source_session": "previous_session.sfppy.json"
}
```

---

## 2. Example Mappings

### 2.1 example1.py (Monolayer, Fully Numeric)

**Key features to encode:**
- Cylinder geometry (19 cm height, 30 mm radius)
- Two substances: Irganox 1076, Irgafos 168
- Single LDPE layer (100 um, C0=5000 mg/kg)
- Food: semisolid + fat
- Single contact step: 10 days at 7C

```json
{
  "version": "1.0",
  "metadata": {
    "name": "LDPE Sandwich Migration",
    "source_script": "example1.py"
  },
  "geometry": {
    "shape": "cylinder",
    "dimensions": {
      "height": { "value": 19, "unit": "cm" },
      "radius": { "value": 30, "unit": "mm" }
    }
  },
  "substances": [
    { "id": "irganox-1076", "lookup_name": "irganox 1076", "SML": 6.0, "color": "#EF4444" },
    { "id": "irgafos-168", "lookup_name": "irgafos 168", "SML": 6.0, "color": "#3B82F6" }
  ],
  "layers": [
    {
      "index": 1,
      "polymer": "LDPE",
      "thickness": { "value": 100, "unit": "um" },
      "substances": [
        { "substance_id": "irganox-1076", "C0": { "value": 5000, "unit": "mg/kg" } },
        { "substance_id": "irgafos-168", "C0": { "value": 5000, "unit": "mg/kg" } }
      ]
    }
  ],
  "food": {
    "category": "realfood",
    "texture": "semisolid",
    "affinity": "fat",
    "simulant": "ethanol",
    "name": "sandwich"
  },
  "contact_steps": [
    {
      "index": 1,
      "temperature": { "value": 7, "unit": "degC" },
      "duration": { "value": 10, "unit": "days" },
      "with_food_contact": true
    }
  ]
}
```

### 2.2 example1_extensions.py (Named + Two-step)

**Additional features:**
- Two contact steps: 10d@7C + 4h@25C

```json
"contact_steps": [
  { "index": 1, "name": "Cold storage", "temperature": {"value": 7, "unit": "degC"}, "duration": {"value": 10, "unit": "days"} },
  { "index": 2, "name": "Warming", "temperature": {"value": 25, "unit": "degC"}, "duration": {"value": 4, "unit": "hours"} }
]
```

### 2.3 example2.py (Multilayer + Functional Barrier)

**Key features:**
- Bottle geometry (body + neck)
- PP layer (300 um) + PET functional barrier (30 um)
- Substance: toluene at 10 mg/kg
- Long contact: 450 days at 20C

```json
{
  "geometry": {
    "shape": "bottle",
    "dimensions": {
      "body_radius": { "value": 40, "unit": "mm" },
      "body_height": { "value": 200, "unit": "mm" },
      "neck_radius": { "value": 18, "unit": "mm" },
      "neck_height": { "value": 50, "unit": "mm" }
    }
  },
  "layers": [
    { "index": 1, "polymer": "wPET", "thickness": {"value": 30, "unit": "um"}, "substances": [{"substance_id": "toluene", "C0": {"value": 0, "unit": "mg/kg"}}] },
    { "index": 2, "polymer": "PP", "thickness": {"value": 300, "unit": "um"}, "substances": [{"substance_id": "toluene", "C0": {"value": 10, "unit": "mg/kg"}}] }
  ],
  "contact_steps": [
    { "index": 1, "temperature": {"value": 20, "unit": "degC"}, "duration": {"value": 450, "unit": "days"} }
  ]
}
```

### 2.4 example3.py (ABA Trilayer + Multi-step)

**Key features:**
- Box container geometry
- ABA structure: wPET(20um) + PP(500um) + gPET(20um)
- Limonene at 200 mg/kg in PP only
- Three contact steps: setoff + hotfilling + storage

```json
{
  "geometry": {
    "shape": "box_container",
    "dimensions": {
      "length": { "value": 19, "unit": "cm" },
      "width": { "value": 10, "unit": "cm" },
      "height": { "value": 8, "unit": "cm" }
    }
  },
  "layers": [
    { "index": 1, "polymer": "wPET", "thickness": {"value": 20, "unit": "um"}, "substances": [{"substance_id": "limonene", "C0": {"value": 0, "unit": "mg/kg"}}] },
    { "index": 2, "polymer": "PP", "thickness": {"value": 500, "unit": "um"}, "substances": [{"substance_id": "limonene", "C0": {"value": 200, "unit": "mg/kg"}}] },
    { "index": 3, "polymer": "gPET", "thickness": {"value": 20, "unit": "um"}, "substances": [{"substance_id": "limonene", "C0": {"value": 0, "unit": "mg/kg"}}] }
  ],
  "contact_steps": [
    { "index": 1, "name": "Setoff storage", "type": "setoff", "temperature": {"value": 20, "unit": "degC"}, "duration": {"value": 4, "unit": "months"}, "with_food_contact": false },
    { "index": 2, "name": "Hot filling", "type": "hotfilling", "temperature": {"value": 90, "unit": "degC"}, "duration": {"value": 30, "unit": "minutes"}, "with_food_contact": true },
    { "index": 3, "name": "Long storage", "type": "storage", "temperature": {"value": 30, "unit": "degC"}, "duration": {"value": 6, "unit": "months"}, "with_food_contact": true }
  ]
}
```

### 2.5 example3_shortvariant.py

Same as example3.py but demonstrates:
- Layer thickness variants: `[10e-6, 0.5e-3, 10e-6]`
- Substance swap: limonene -> toluene

**Encoding variants:**
```json
"variants": [
  {
    "name": "Reference",
    "layer_overrides": null,
    "substance_override": null
  },
  {
    "name": "Toluene variant",
    "substance_override": "toluene"
  },
  {
    "name": "Thin FB",
    "layer_overrides": {
      "1": { "thickness": {"value": 10, "unit": "um"} },
      "3": { "thickness": {"value": 10, "unit": "um"} }
    }
  }
]
```

---

## 3. Implementation Phases

### Phase 1: Core Session Format (Priority: HIGH)

**Tasks:**
1. Define JSON schema validator (jsonschema library)
2. Create `SessionModel` dataclass with Pydantic
3. Implement serialization/deserialization
4. Add unit conversion utilities

**Files:**
- `studio/app/models/session.py` - Session data models
- `studio/app/utils/session_io.py` - Save/Load functions
- `studio/app/schemas/session_v1.json` - JSON schema

### Phase 2: UI Integration (Priority: HIGH)

**Tasks:**
1. Add Save Session button (Config tab)
2. Add Load Session button (Config tab)
3. Implement state hydration from session file
4. Handle validation errors gracefully

**API Endpoints:**
- `POST /api/sessions/save` - Save current state
- `POST /api/sessions/load` - Load from file
- `GET /api/sessions/list` - List saved sessions
- `DELETE /api/sessions/{id}` - Delete session

### Phase 3: Script Conversion Tool (Priority: MEDIUM)

**Tasks:**
1. Create CLI tool: `python -m studio.convert example1.py`
2. Parse Python AST to extract configuration
3. Map SFPPy constructs to session JSON
4. Handle dynamic/computed values

**Approach:**
- Use `ast` module to parse Python scripts
- Extract `Packaging3D`, `migrant`, `layer.*` calls
- Map to session JSON structure
- Flag unsupported constructs for manual review

### Phase 4: Restart Functionality (Priority: MEDIUM)

**Tasks:**
1. Extend Patankar solver to accept initial CF
2. Add UI for "Continue from previous" option
3. Implement layer concentration profile restoration
4. Chain simulations from arbitrary state

**Technical Notes:**
- Patankar solver needs `C0_food` parameter (initial CF in food)
- May need to reconstruct layer Cx profile from previous simulation
- Contact step index determines which step to restart from

### Phase 5: Built-in Examples (Priority: LOW)

**Tasks:**
1. Convert example1-3 to `.sfppy.json` files
2. Bundle with Studio in `studio/examples/`
3. Add "Load Example" dropdown in UI
4. Validate against schema on startup

---

## 4. File Standards

### 4.1 Session Files

| Aspect | Standard |
|--------|----------|
| Extension | `.sfppy.json` |
| Encoding | UTF-8 |
| Indentation | 2 spaces |
| Numbers | Max 6 significant digits for floats |
| Units | Always explicit `{"value": x, "unit": "y"}` |
| IDs | UUID v4 or `{compound}-{cid}` format |

### 4.2 Unit Conventions

| Quantity | Allowed Units |
|----------|---------------|
| Length | nm, um, mm, cm, m |
| Time | s, min, hours, days, weeks, months, years |
| Temperature | degC, K |
| Concentration | mg/kg, ppm, g/kg, % |
| Diffusivity | m2/s (always SI) |
| Partition coefficient | dimensionless |

### 4.3 Validation Rules

1. **Required fields**: geometry, layers (>=1), contact_steps (>=1)
2. **Substances**: Must have either `lookup_name` or `properties.name`
3. **Layer references**: All `substance_id` must exist in `substances[]`
4. **Units**: All physical quantities must have explicit units
5. **Indices**: Must be consecutive integers starting at 1

---

## 5. Restart from Previous State

### 5.1 Use Cases

1. **Continue simulation**: Run 10 more days after initial 10 days
2. **What-if analysis**: Start from CF=2 mg/kg, vary temperature
3. **Multi-stage process**: Different conditions for each stage
4. **Validation**: Compare against experimental data at t>0

### 5.2 Technical Implementation

```python
# In Patankar solver, add initial_CF parameter
def senspatankar(layer, food, initial_CF=0, initial_Cx=None, ...):
    """
    initial_CF: float - Initial concentration in food (mg/kg)
    initial_Cx: array - Initial concentration profile in layers
    """
    # Adjust mass balance for non-zero initial CF
    # food.C0 = initial_CF
    # layer.Cx = initial_Cx if provided
```

### 5.3 Session Encoding

```json
"restart_state": {
  "enabled": true,
  "from_step": 2,
  "initial_CF": {
    "irganox-1076": 2.5,
    "irgafos-168": 1.2
  },
  "source_results": {
    "session_id": "abc123",
    "timestamp": "2025-01-15T10:30:00Z"
  }
}
```

---

## 6. Design Decisions (APPROVED)

| Question | Decision |
|----------|----------|
| **Variant storage** | Separate sessions (independent files) |
| **Cx profiles storage** | JSON in session, CSV export for verification |
| **Kinetics export** | CSV format for easy verification |
| **Batch processing** | Supported - queue multiple sessions |

### 6.1 CSV Export Format

**Kinetics (CF vs time):**
```csv
time_days,substance_1_CF_mg_kg,substance_2_CF_mg_kg,...
0.0,0.0,0.0
0.1,0.5,0.3
1.0,2.1,1.5
...
```

**Concentration profiles (Cx vs x at selected times):**
```csv
x_um,t=0.1_days,t=1.0_days,t=5.0_days,...
0.0,1000.0,950.2,850.1
1.0,999.8,948.5,845.3
...
```

### 6.2 Batch Processing

```json
{
  "batch": {
    "name": "FB thickness study",
    "sessions": [
      "fb_10um.sfppy.json",
      "fb_20um.sfppy.json",
      "fb_30um.sfppy.json"
    ],
    "parallel": true,
    "output_folder": "results/fb_study/"
  }
}
```

---

## 7. Testing Strategy

### 7.1 Phase 1 Tests: Core Session Format

**Test file:** `studio/tests/test_session_format.py`

```python
# Test cases:
def test_session_schema_validation():
    """Valid session passes schema validation."""

def test_session_schema_rejects_invalid():
    """Invalid sessions rejected with clear errors."""

def test_unit_conversion():
    """Units converted correctly (um->m, days->s, etc.)."""

def test_session_save_load_roundtrip():
    """Save then load produces identical session."""

def test_example1_encoding():
    """example1.py configuration encodes correctly."""

def test_example2_encoding():
    """example2.py (multilayer) encodes correctly."""

def test_example3_encoding():
    """example3.py (multi-step) encodes correctly."""
```

### 7.2 Phase 2 Tests: UI Integration

**Test file:** `studio/tests/test_session_api.py`

```python
def test_save_session_endpoint():
    """POST /api/sessions/save returns valid session file."""

def test_load_session_endpoint():
    """POST /api/sessions/load hydrates UI state correctly."""

def test_list_sessions_endpoint():
    """GET /api/sessions/list returns available sessions."""

def test_invalid_session_rejected():
    """Loading invalid session returns 400 with error details."""
```

### 7.3 Phase 3 Tests: Script Conversion

**Test file:** `studio/tests/test_script_conversion.py`

```python
def test_convert_example1():
    """example1.py converts to valid session JSON."""

def test_convert_example2():
    """example2.py converts to valid session JSON."""

def test_convert_example3():
    """example3.py converts to valid session JSON."""

def test_converted_session_produces_same_results():
    """Simulation from converted session matches original script output."""
```

### 7.4 Phase 4 Tests: Restart Functionality

**Test file:** `studio/tests/test_restart.py`

```python
def test_restart_with_initial_CF():
    """Restart from CF=2 mg/kg gives correct continuation."""

def test_restart_mass_balance():
    """Total mass conserved when restarting from previous state."""

def test_restart_from_step_2():
    """Multi-step restart from step 2 skips step 1 correctly."""
```

### 7.5 Phase 5 Tests: Batch Processing

**Test file:** `studio/tests/test_batch.py`

```python
def test_batch_queue_multiple_sessions():
    """Batch queues and runs multiple sessions."""

def test_batch_parallel_execution():
    """Parallel batch runs faster than sequential."""

def test_batch_csv_export():
    """Batch exports combined CSV results."""
```

### 7.6 Reference Test Data

Create reference fixtures from actual SFPPy runs:
- `studio/tests/fixtures/example1_reference.json` - Expected results
- `studio/tests/fixtures/example2_reference.json`
- `studio/tests/fixtures/example3_reference.json`

**Generate with:**
```bash
cd /home/olivi/natacha/python
python example1.py  # Save CF values at key times
python example2.py
python example3.py
```

---

## 8. Timeline (Suggested Order)

| Phase | Description | Dependencies |
|-------|-------------|--------------|
| 1 | Core session format | None |
| 2 | UI save/load | Phase 1 |
| 3 | Script conversion | Phase 1 |
| 4 | Restart functionality | Phase 1, Phase 2 |
| 5 | Built-in examples | Phase 3 |

---

## Approval Required

Before implementation:
1. Review session schema completeness
2. Confirm unit conventions
3. Validate example mappings
4. Agree on restart mechanism design
