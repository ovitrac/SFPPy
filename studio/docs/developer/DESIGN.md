# SFPPy Studio — Design Document

**Version:** 0.1 (Draft)
**Author:** Olivier Vitrac, PhD, HDR
**Organization:** INRAE + Generative Simulation
**Date:** 2025-12-24

---

## 0. Version & Traceability

### 0.1 Version Management

```python
# studio/version.py
__version__ = "0.1.0"
__version_info__ = (0, 1, 0)

# Import SFPPy version from utils
from utils import SFPPy_version  # e.g., "1.50"
```

### 0.2 Traceability Information (captured per job)

| Field | Source | Example |
|-------|--------|---------|
| `sfppy_version` | `utils.SFPPy_version` | "1.50" |
| `studio_version` | `studio.version.__version__` | "0.1.0" |
| `user_name` | User profile or system | "olivier.vitrac" |
| `user_email` | User profile | "olivier.vitrac@gmail.com" |
| `machine_name` | `socket.gethostname()` | "workstation-01" |
| `machine_os` | `platform.system()` | "Linux" |
| `python_version` | `sys.version` | "3.11.5" |
| `created_at` | `datetime.utcnow()` | "2025-12-24T15:30:00Z" |
| `modified_at` | Auto-updated | "2025-12-24T16:45:00Z" |

### 0.3 Job Metadata Structure

```python
@dataclass
class JobMetadata:
    # Version info
    sfppy_version: str
    studio_version: str

    # User info
    user_name: str
    user_email: Optional[str]
    user_organization: Optional[str]

    # Machine info
    machine_name: str
    machine_os: str
    python_version: str

    # Timestamps
    created_at: datetime
    modified_at: datetime

    # Job identity
    job_id: str
    job_name: str
    parent_job_id: Optional[str]  # If cloned/rerun from another job
```

---

## 1. Overview

**SFPPy Studio** is a comprehensive web application for food contact migration analysis. It provides a complete workflow from material/food definition through simulation to regulatory compliance reporting.

### 1.1 Design Principles

- **Multi-layer first**: All scenarios support 1-10 layers from the start
- **Abstract pipeline**: Predefined components (polymers, foods, conditions) assembled via operators
- **Notebook-compatible**: Workflow mirrors `gui.ipynb` / `demo.ipynb` patterns
- **Full traceability**: All jobs stored with edit/relaunch capability
- **Regulatory focus**: Built-in SML comparison and compliance assessment

### 1.2 Technology Stack

| Component | Technology |
|-----------|------------|
| Backend | FastAPI |
| Frontend | Jinja2 + Tailwind CSS + Chart.js |
| Solver | `senspatankar` from `patankar.migration` |
| Database | SQLite for job storage |
| Export | PDF (matplotlib), CSV, Excel |

---

## 2. Application Structure

```
studio/
├── app/
│   ├── __init__.py
│   ├── main.py              # FastAPI app entry point
│   ├── routes/
│   │   ├── assembly.py      # Material assembly endpoints
│   │   ├── substances.py    # Substance search/management
│   │   ├── simulation.py    # Simulation execution
│   │   ├── jobs.py          # Job history/management
│   │   └── reports.py       # Report generation
│   ├── models/
│   │   ├── assembly.py      # Layer, Material models
│   │   ├── food.py          # Food, Simulant models
│   │   ├── substance.py     # Substance models
│   │   ├── scenario.py      # Multi-step scenario
│   │   └── job.py           # Job persistence
│   ├── services/
│   │   ├── solver.py        # senspatankar wrapper
│   │   ├── estimators.py    # D, k estimation
│   │   ├── pubchem.py       # PubChem integration
│   │   └── compliance.py    # Regulatory checks
│   ├── templates/
│   │   ├── base.html        # Base template with tabs
│   │   ├── assembly.html    # Tab 1: Assembly
│   │   ├── food.html        # Tab 2: Food & Conditions
│   │   ├── substances.html  # Tab 3: Substances
│   │   ├── simulation.html  # Tab 4: Simulation
│   │   ├── results.html     # Tab 5: Results
│   │   ├── jobs.html        # Tab 6: Job History
│   │   └── config.html      # Tab 7: Configuration
│   └── static/
│       ├── css/
│       ├── js/
│       └── img/
├── data/
│   ├── polymers.yaml        # Predefined polymer library
│   ├── foods.yaml           # Predefined food types
│   ├── simulants.yaml       # Food simulants
│   ├── conditions.yaml      # Contact conditions presets
│   └── regulatory.yaml      # SML, FCM data
├── jobs/                    # Job storage (JSON/SQLite)
├── exports/                 # Generated reports
├── launcher.py              # Application launcher
├── DESIGN.md                # This document
└── README.md                # Quick start guide
```

---

## 3. Tab Structure

### Tab Navigation

```
┌────────────────────────────────────────────────────────────────────────────┐
│  🧪 SFPPy Studio — Food Contact Migration Analysis                         │
├────────────────────────────────────────────────────────────────────────────┤
│  📦 Assembly │ 🍽️ Food │ ⚗️ Substances │ 🔬 Simulate │ 📊 Results │      │
│              │         │               │             │            │ 📋 Jobs │ ⚙️ │
└────────────────────────────────────────────────────────────────────────────┘
```

---

## 4. Tab 1: 📦 Assembly (Material Layers)

### 4.1 Purpose
Define the multilayer material structure (1-10 layers).

### 4.2 Layout

```
┌─────────────────────────────────────────────────────────────────────────────┐
│ 📦 Material Assembly                                                        │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│  Number of layers: [◀] 3 [▶]     Total thickness: 380 µm                   │
│                                                                             │
│  ┌─────────────────────────────────────────────────────────────────────┐   │
│  │ FOOD ║ Layer 1    ║ Layer 2      ║ Layer 3      │                   │   │
│  │ ←──→ ║ PET 30µm   ║ rPP 300µm    ║ LDPE 50µm    │                   │   │
│  │      ║ [FB] ✓     ║              ║              │                   │   │
│  └─────────────────────────────────────────────────────────────────────┘   │
│                                                                             │
│  ┌─────────────────────────────────────────────────────────────────────┐   │
│  │ Layer 1 (Food Contact)                                    [Remove]  │   │
│  │ ┌─────────────┐ ┌─────────────┐ ┌─────────────┐ ┌────────────────┐ │   │
│  │ │ Polymer     │ │ Thickness   │ │ Temperature │ │ Options        │ │   │
│  │ │ [PET     ▼] │ │ [30] [µm ▼] │ │ [25] °C     │ │ ☐ Glassy (FB)  │ │   │
│  │ └─────────────┘ └─────────────┘ └─────────────┘ │ ☐ Plasticized  │ │   │
│  │                                                  └────────────────┘ │   │
│  │ Estimated D: 1.2e-14 m²/s (at 25°C)    k: 1.0                      │   │
│  └─────────────────────────────────────────────────────────────────────┘   │
│                                                                             │
│  ┌─────────────────────────────────────────────────────────────────────┐   │
│  │ Layer 2                                                   [Remove]  │   │
│  │ ... (similar structure)                                             │   │
│  └─────────────────────────────────────────────────────────────────────┘   │
│                                                                             │
│  [+ Add Layer]                                                              │
│                                                                             │
├─────────────────────────────────────────────────────────────────────────────┤
│ 💡 Help: Layer 1 is always in contact with food. Use [◀][▶] to adjust the  │
│    number of layers. Glassy polymers (PET, PA) act as functional barriers. │
└─────────────────────────────────────────────────────────────────────────────┘
```

### 4.3 Parameter Override Panel

**Critical for Example 1:** Users must be able to override computed D, k, k0 values.

```
┌─────────────────────────────────────────────────────────────────────────────┐
│ Layer 2 (rPP) — Parameter Overrides                                         │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│  ┌─────────────────────────────────────────────────────────────────────┐   │
│  │ Diffusivity (D)                                              [?]    │   │
│  │ ○ Use computed (Piringer model): 1.2e-14 m²/s                       │   │
│  │ ● Override: [1.5e-14    ] m²/s                                      │   │
│  │                                                                      │   │
│  │ 💡 Computed D is estimated from polymer type, temperature, and MW.  │   │
│  │    Override when you have experimental or literature values.        │   │
│  └─────────────────────────────────────────────────────────────────────┘   │
│                                                                             │
│  ┌─────────────────────────────────────────────────────────────────────┐   │
│  │ Partition coefficient (k)                                    [?]    │   │
│  │ ○ Use computed (FHP model): 1.0                                     │   │
│  │ ● Override: [0.8        ]                                           │   │
│  │                                                                      │   │
│  │ 💡 k = C_polymer / C_food at equilibrium. Override for specific     │   │
│  │    polymer/food combinations with known experimental values.        │   │
│  └─────────────────────────────────────────────────────────────────────┘   │
│                                                                             │
│  ┌─────────────────────────────────────────────────────────────────────┐   │
│  │ Reference partition (k0)                                     [?]    │   │
│  │ ○ Use default: 1.0                                                  │   │
│  │ ● Override: [1.2        ]                                           │   │
│  │                                                                      │   │
│  │ 💡 k0 is the reference partition coefficient for the simulant.      │   │
│  │    Used to scale k when switching from simulant to real food.       │   │
│  └─────────────────────────────────────────────────────────────────────┘   │
│                                                                             │
│  [Reset to Computed Values]                                                 │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

### 4.4 Data Model

```python
@dataclass
class ParameterOverride:
    use_computed: bool = True
    override_value: Optional[float] = None
    computed_value: Optional[float] = None  # For display
    model_used: Optional[str] = None        # "Piringer", "FHP", etc.

@dataclass
class Layer:
    index: int                    # 1 = food contact
    polymer: str                  # "PET", "LDPE", "PP", etc.
    thickness: float              # in meters (SI)
    thickness_unit: str           # display unit
    temperature: float            # °C
    is_glassy: bool              # functional barrier flag
    is_plasticized: bool

    # Parameter overrides (critical for flexibility)
    D_override: ParameterOverride = field(default_factory=ParameterOverride)
    k_override: ParameterOverride = field(default_factory=ParameterOverride)
    k0_override: ParameterOverride = field(default_factory=ParameterOverride)

    # Calculated (from substance assignment)
    substances: List[SubstanceAssignment] = field(default_factory=list)

    def get_D(self) -> float:
        """Return override value if set, else computed."""
        if not self.D_override.use_computed and self.D_override.override_value:
            return self.D_override.override_value
        return self.D_override.computed_value

    def get_k(self) -> float:
        """Return override value if set, else computed."""
        if not self.k_override.use_computed and self.k_override.override_value:
            return self.k_override.override_value
        return self.k_override.computed_value

@dataclass
class Assembly:
    name: str
    layers: List[Layer]          # index 1..n
    total_thickness_m: float     # computed
```

### 4.4 Predefined Polymers (`data/polymers.yaml`)

```yaml
polymers:
  - code: LDPE
    name: Low Density Polyethylene
    Tg: -120  # °C
    density: 920  # kg/m³
    icon: "♳"

  - code: HDPE
    name: High Density Polyethylene
    Tg: -120
    density: 960
    icon: "♴"

  - code: PP
    name: Polypropylene
    Tg: -10
    density: 905
    icon: "♵"

  - code: PET
    name: Polyethylene Terephthalate
    Tg: 75
    density: 1380
    icon: "♳"
    variants:
      - code: gPET
        name: Glassy PET
        description: Below Tg, functional barrier
      - code: wPET
        name: Wet/Plasticized PET
        description: Above Tg or plasticized

  - code: PS
    name: Polystyrene
    Tg: 100
    density: 1050
    icon: "♶"

  # ... more polymers
```

---

## 5. Tab 2: 🍽️ Food & Conditions

### 5.1 Purpose
Define food type, simulant, geometry, and contact conditions.

### 5.2 Layout

```
┌─────────────────────────────────────────────────────────────────────────────┐
│ 🍽️ Food & Contact Conditions                                               │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│  ┌─────────────────────────────────┐  ┌─────────────────────────────────┐  │
│  │ Food Type                       │  │ Packaging Geometry              │  │
│  │                                 │  │                                 │  │
│  │ Category: [Fatty foods    ▼]   │  │ Shape: [Cylinder         ▼]    │  │
│  │ Texture:  [Semisolid      ▼]   │  │                                 │  │
│  │ Process:  [Ambient        ▼]   │  │ Radius: [40] mm                 │  │
│  │                                 │  │ Height: [200] mm                │  │
│  │ 🥪 Sandwich (fatty, semisolid) │  │                                 │  │
│  │                                 │  │ Volume: 1005.3 cm³              │  │
│  │ Simulant: [Ethanol 50%    ▼]   │  │ Surface: 602.9 cm²              │  │
│  │ Override: ☐ Use real food      │  │ V/S ratio: 1.67 cm              │  │
│  └─────────────────────────────────┘  └─────────────────────────────────┘  │
│                                                                             │
│  ┌─────────────────────────────────────────────────────────────────────┐   │
│  │ Contact Conditions                                                   │   │
│  │                                                                      │   │
│  │ ○ Single step    ● Multi-step chain                                 │   │
│  │                                                                      │   │
│  │ Step 1: [Storage    ▼]  Temp: [25]°C  Duration: [4] [months▼]       │   │
│  │         ☐ Without food (set-off)                                    │   │
│  │                                                                      │   │
│  │ Step 2: [Hot-fill   ▼]  Temp: [85]°C  Duration: [30] [min▼]         │   │
│  │         ☑ With food                                                 │   │
│  │                                                                      │   │
│  │ Step 3: [Storage    ▼]  Temp: [4]°C   Duration: [6] [months▼]       │   │
│  │         ☑ With food                                                 │   │
│  │                                                                      │   │
│  │ [+ Add Step]  [Remove Last]                                         │   │
│  └─────────────────────────────────────────────────────────────────────┘   │
│                                                                             │
├─────────────────────────────────────────────────────────────────────────────┤
│ 💡 Help: Multi-step chains simulate real-world scenarios (set-off during    │
│    storage → hot-fill → long-term storage). Results will be merged.        │
└─────────────────────────────────────────────────────────────────────────────┘
```

### 5.3 Data Model

```python
@dataclass
class Food:
    category: str          # "fatty", "aqueous", "acidic", "alcoholic"
    texture: str           # "liquid", "semisolid", "solid"
    process: str           # "ambient", "hotfilled", "frozen"
    name: str              # Human readable
    simulant: str          # "ethanol", "ethanol50", "olive_oil", etc.
    use_real_food: bool

@dataclass
class Geometry:
    shape: str             # "cylinder", "bottle", "rectangle"
    dimensions: dict       # shape-specific
    volume_m3: float
    surface_m2: float

@dataclass
class ContactStep:
    index: int
    condition_type: str    # "storage", "hotfill", "transport"
    temperature_C: float
    duration: float
    duration_unit: str
    with_food: bool        # False = set-off (no food contact)

@dataclass
class Scenario:
    food: Food
    geometry: Geometry
    steps: List[ContactStep]
```

---

## 6. Tab 3: ⚗️ Substances

### 6.1 Purpose
Search, visualize, and assign substances to layers with initial concentrations.

### 6.2 Integration with loadpubchem.py

The backend uses `patankar/loadpubchem.py` which provides:
- **Local cache** in `patankar/cache.PubChem/` (JSON files per CID)
- **PubChem API** connection for new searches
- **Thumbnail images** in `cache.PubChem/thumbs/`
- **Regulatory data** from embedded databases

```python
# Backend service using loadpubchem
from patankar.loadpubchem import migrant, migrantToxtree

def search_substance(query: str) -> List[dict]:
    """Search PubChem with local cache."""
    m = migrant(query)  # Uses cache first, then PubChem API
    return {
        "name": m.name,
        "cid": m.cid,
        "cas": m.CAS,
        "mw": m.M,
        "logP": m.logP,
        "smiles": m.SMILES,
        "image_path": m.get_image_path(),  # Local cached thumbnail
        "regulatory": {
            "eu_sml": m.SML,
            "eu_fcm": m.FCM,
            "fda_fcn": m.FCN,
            "ttc": m.TTC
        }
    }
```

### 6.3 Layout

```
┌─────────────────────────────────────────────────────────────────────────────┐
│ ⚗️ Substances                                                               │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│  ┌───────────────────────────────────────────────────────────────────────┐ │
│  │ 🔍 Search Substance                                             [?]   │ │
│  │ ─────────────────────────────────────────────────────────────────────  │ │
│  │                                                                        │ │
│  │ Search by: ○ Name  ○ CAS number  ○ PubChem CID                        │ │
│  │                                                                        │ │
│  │ Query: [Irganox 1076                         ] [🔍 Search]            │ │
│  │                                                                        │ │
│  │ 💡 Search uses local cache first (fast), then PubChem API if needed.  │ │
│  │    Common additives are pre-cached for instant results.               │ │
│  └───────────────────────────────────────────────────────────────────────┘ │
│                                                                             │
│  ┌───────────────────────────────────────────────────────────────────────┐ │
│  │ Search Results                                              3 found   │ │
│  │ ─────────────────────────────────────────────────────────────────────  │ │
│  │                                                                        │ │
│  │ ┌─────────────────────────────────────────────────────────────────┐   │ │
│  │ │ ○ Irganox 1076                                                   │   │ │
│  │ │   ┌────────┐  CAS: 2082-79-3    CID: 91597                      │   │ │
│  │ │   │ [IMG]  │  Octadecyl 3,5-di-tert-butyl-4-                    │   │ │
│  │ │   │        │  hydroxyhydrocinnamate                              │   │ │
│  │ │   └────────┘  MW: 530.87   logP: 12.4                           │   │ │
│  │ │                                                     [+ Add]      │   │ │
│  │ └─────────────────────────────────────────────────────────────────┘   │ │
│  │                                                                        │ │
│  │ ┌─────────────────────────────────────────────────────────────────┐   │ │
│  │ │ ○ Irganox 1010                                                   │   │ │
│  │ │   ┌────────┐  CAS: 6683-19-8    CID: 71306                      │   │ │
│  │ │   │ [IMG]  │  Pentaerythritol tetrakis(3,5-di-tert-             │   │ │
│  │ │   │        │  butyl-4-hydroxyhydrocinnamate)                     │   │ │
│  │ │   └────────┘  MW: 1177.63  logP: 23.0                           │   │ │
│  │ │                                                     [+ Add]      │   │ │
│  │ └─────────────────────────────────────────────────────────────────┘   │ │
│  └───────────────────────────────────────────────────────────────────────┘ │
│                                                                             │
│  ┌───────────────────────────────────────────────────────────────────────┐ │
│  │ Assigned Substances (2)                                               │ │
│  │ ─────────────────────────────────────────────────────────────────────  │ │
│  │                                                                        │ │
│  │ ┌─────────────────────────────────────────────────────────────────┐   │ │
│  │ │ 🧪 Irganox 1076                                       [Remove]   │   │ │
│  │ │ ─────────────────────────────────────────────────────────────── │   │ │
│  │ │ ┌──────────┐                                                     │   │ │
│  │ │ │          │  CAS: 2082-79-3        MW: 530.87 g/mol            │   │ │
│  │ │ │  [2D     │  logP: 12.4            SMILES: CCCCCCCC...         │   │ │
│  │ │ │  struct] │                                                     │   │ │
│  │ │ │          │  ┌─────────────────────────────────────────────┐   │   │ │
│  │ │ └──────────┘  │ Regulatory Status                     [?]   │   │   │ │
│  │ │               │ ─────────────────────────────────────────── │   │   │ │
│  │ │               │ EU 10/2011: SML = 6 mg/kg  FCM #397   ✅    │   │   │ │
│  │ │               │ US FDA:     FCN 000XXX                 ✅    │   │   │ │
│  │ │               │ China GB:   Authorized                 ✅    │   │   │ │
│  │ │               │ ToxTree:    Class III  TTC=0.09 µg/d        │   │   │ │
│  │ │               └─────────────────────────────────────────────┘   │   │ │
│  │ │                                                                  │   │ │
│  │ │  Layer Assignment                                         [?]   │   │ │
│  │ │  ───────────────────────────────────────────────────────────── │   │ │
│  │ │  ☐ Layer 1 (PET 30µm)       C0: [________] mg/kg               │   │ │
│  │ │  ☑ Layer 2 (rPP 300µm)      C0: [1000    ] mg/kg               │   │ │
│  │ │  ☐ Layer 3 (LDPE 50µm)      C0: [________] mg/kg               │   │ │
│  │ │                                                                  │   │ │
│  │ │  💡 C0 = initial concentration in the layer before migration.   │   │ │
│  │ │     Only layers containing the substance should be checked.     │   │ │
│  │ └─────────────────────────────────────────────────────────────────┘   │ │
│  │                                                                        │ │
│  │ [+ Search More Substances]                                             │ │
│  └───────────────────────────────────────────────────────────────────────┘ │
│                                                                             │
├─────────────────────────────────────────────────────────────────────────────┤
│ 💡 Help: Substances are searched via PubChem with local caching for speed. │
│    SML = Specific Migration Limit (EU 10/2011). FCM = Food Contact         │
│    Material number. TTC = Threshold of Toxicological Concern.              │
└─────────────────────────────────────────────────────────────────────────────┘
```

### 6.4 Substance Visualization Modal

When clicking on substance image or [View Details]:

```
┌─────────────────────────────────────────────────────────────────────────────┐
│ 🧪 Irganox 1076 — Detailed View                                    [Close] │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│  ┌─────────────────────────┐  ┌─────────────────────────────────────────┐  │
│  │                         │  │ Identifiers                             │  │
│  │     [2D Structure]      │  │ Name: Octadecyl 3,5-di-tert-butyl-4-   │  │
│  │        300x300          │  │       hydroxyhydrocinnamate             │  │
│  │                         │  │ CAS: 2082-79-3                          │  │
│  │                         │  │ PubChem CID: 91597                      │  │
│  │                         │  │ InChIKey: SSDSCDGVMJFTEQ-UHFFFAOYSA-N  │  │
│  └─────────────────────────┘  └─────────────────────────────────────────┘  │
│                                                                             │
│  ┌─────────────────────────────────────────────────────────────────────┐   │
│  │ Physical Properties                                                  │   │
│  │ ─────────────────────────────────────────────────────────────────── │   │
│  │ Molecular Weight: 530.87 g/mol    Molecular Formula: C35H62O3       │   │
│  │ logP: 12.4                        Melting Point: 50-55°C            │   │
│  │ Density: ~1.0 g/cm³               Water Solubility: Insoluble       │   │
│  └─────────────────────────────────────────────────────────────────────┘   │
│                                                                             │
│  ┌─────────────────────────────────────────────────────────────────────┐   │
│  │ Regulatory Status                                                    │   │
│  │ ─────────────────────────────────────────────────────────────────── │   │
│  │                                                                      │   │
│  │ 🇪🇺 EU 10/2011                                                       │   │
│  │    FCM Number: 397                                                   │   │
│  │    SML: 6 mg/kg food                                                │   │
│  │    Status: ✅ Authorized                                             │   │
│  │                                                                      │   │
│  │ 🇺🇸 US FDA                                                            │   │
│  │    FCN: Listed                                                       │   │
│  │    Status: ✅ Authorized                                             │   │
│  │                                                                      │   │
│  │ 🇨🇳 China GB 9685-2016                                                │   │
│  │    Status: ✅ Authorized                                             │   │
│  │                                                                      │   │
│  │ ☠️ ToxTree Assessment                                                 │   │
│  │    Cramer Class: III (high)                                         │   │
│  │    TTC: 0.09 µg/day (1.5 µg/kg bw/day)                              │   │
│  └─────────────────────────────────────────────────────────────────────┘   │
│                                                                             │
│  ┌─────────────────────────────────────────────────────────────────────┐   │
│  │ SMILES                                                               │   │
│  │ CCCCCCCCCCCCCCCCCCOC(=O)CCc1cc(C(C)(C)C)c(O)c(C(C)(C)C)c1           │   │
│  └─────────────────────────────────────────────────────────────────────┘   │
│                                                                             │
│  [Open in PubChem ↗]  [Copy Data]  [Add to Scenario]                       │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

### 6.3 Data Model

```python
@dataclass
class Substance:
    name: str
    cas: str
    cid: int                    # PubChem CID
    mw: float                   # Molecular weight
    logP: float
    smiles: str
    image_url: str
    # Regulatory
    eu_sml: Optional[float]     # mg/kg
    eu_fcm: Optional[int]
    fda_fcn: Optional[str]
    ttc: Optional[float]        # µg/day

@dataclass
class SubstanceAssignment:
    substance: Substance
    layer_assignments: Dict[int, float]  # layer_index → C0 (mg/kg)
```

---

## 7. Tab 4: 🔬 Simulate

### 7.1 Purpose
Configure solver, run simulation, link to previous results.

### 7.2 Layout

```
┌─────────────────────────────────────────────────────────────────────────────┐
│ 🔬 Simulation                                                               │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│  ┌───────────────────────────────────────────────────────────────────────┐ │
│  │ Scenario Summary                                                       │ │
│  │ ─────────────────────────────────────────────────────────────────────  │ │
│  │ Assembly: 3 layers (PET 30µm + rPP 300µm + LDPE 50µm) = 380µm         │ │
│  │ Food: Sandwich (fatty, semisolid) / Ethanol 50%                        │ │
│  │ Conditions: 3 steps (4mo storage → 30min hot-fill → 6mo storage)      │ │
│  │ Substances: Irganox 1076 (layer 2), Irgafos 168 (layer 2)             │ │
│  │                                                                        │ │
│  │ Pipeline: m % food << geometry >> ABA >> food1 >> food2 >> food3      │ │
│  └───────────────────────────────────────────────────────────────────────┘ │
│                                                                             │
│  ┌───────────────────────────────────────────────────────────────────────┐ │
│  │ Previous Results (Optional)                                            │ │
│  │                                                                        │ │
│  │ ☐ Link to previous simulation result                                  │ │
│  │   Job: [Select previous job...              ▼]                        │ │
│  │   This allows continuation from a previous state.                     │ │
│  └───────────────────────────────────────────────────────────────────────┘ │
│                                                                             │
│  ┌───────────────────────────────────────────────────────────────────────┐ │
│  │ Solver Settings                                                        │ │
│  │                                                                        │ │
│  │ Spatial points: [1800]     Time steps: [1000]                         │ │
│  │ D estimator: [Piringer ▼]  k estimator: [FHP ▼]                       │ │
│  │ Workers: [6]               Timeout: [300] s                            │ │
│  │                                                                        │ │
│  │ [Use defaults]                                                         │ │
│  └───────────────────────────────────────────────────────────────────────┘ │
│                                                                             │
│  ┌───────────────────────────────────────────────────────────────────────┐ │
│  │ Job Name: [sandwich_3step_irganox_2024-12-24            ]             │ │
│  │                                                                        │ │
│  │                    [🚀 Run Simulation]                                 │ │
│  └───────────────────────────────────────────────────────────────────────┘ │
│                                                                             │
│  ┌───────────────────────────────────────────────────────────────────────┐ │
│  │ Progress                                                     [Cancel] │ │
│  │ ████████████████████░░░░░░░░░░  68%                                   │ │
│  │ Step 2/3: Hot-fill (85°C, 30min)  ETA: 12s                            │ │
│  └───────────────────────────────────────────────────────────────────────┘ │
│                                                                             │
├─────────────────────────────────────────────────────────────────────────────┤
│ 💡 Help: Review the scenario summary. Link to previous results for chained │
│    simulations. Adjust solver settings in ⚙️ Configuration tab.            │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

## 8. Tab 5: 📊 Results

### 8.1 Purpose
Visualize results, compare with SML, export reports.

### 8.2 Layout

```
┌─────────────────────────────────────────────────────────────────────────────┐
│ 📊 Results                                                                  │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│  Job: sandwich_3step_irganox_2024-12-24          Status: ✅ Completed      │
│                                                                             │
│  ┌───────────────────────────────────────────────────────────────────────┐ │
│  │ CF(t) — Concentration in Food vs Time                    [🔍 Zoom]    │ │
│  │ ┌─────────────────────────────────────────────────────────────────┐   │ │
│  │ │                                                                  │   │ │
│  │ │         📈 (Chart.js interactive chart)                         │   │ │
│  │ │                                                                  │   │ │
│  │ │    - Irganox 1076 (green)                                       │   │ │
│  │ │    - Irgafos 168 (blue)                                         │   │ │
│  │ │    - SML reference (dashed red)                                 │   │ │
│  │ │                                                                  │   │ │
│  │ └─────────────────────────────────────────────────────────────────┘   │ │
│  │ Scale: [Linear ▼]  Show: ☑ Merged steps  ☐ Individual steps          │ │
│  └───────────────────────────────────────────────────────────────────────┘ │
│                                                                             │
│  ┌─────────────────────────────────┐ ┌─────────────────────────────────┐  │
│  │ Compliance Summary              │ │ Key Values                      │  │
│  │                                 │ │                                 │  │
│  │ Irganox 1076:                   │ │ Irganox 1076:                   │  │
│  │   CF(final) = 0.023 mg/kg      │ │   CF(10d) = 0.018 mg/kg        │  │
│  │   SML = 6 mg/kg                │ │   CF(∞) = 0.089 mg/kg          │  │
│  │   Status: ✅ COMPLIANT          │ │   t(95%) = 45 days             │  │
│  │                                 │ │                                 │  │
│  │ Irgafos 168:                    │ │ Irgafos 168:                    │  │
│  │   CF(final) = 0.041 mg/kg      │ │   CF(10d) = 0.032 mg/kg        │  │
│  │   SML = 5 mg/kg                │ │   CF(∞) = 0.156 mg/kg          │  │
│  │   Status: ✅ COMPLIANT          │ │   t(95%) = 62 days             │  │
│  └─────────────────────────────────┘ └─────────────────────────────────┘  │
│                                                                             │
│  ┌───────────────────────────────────────────────────────────────────────┐ │
│  │ Export                                                                 │ │
│  │                                                                        │ │
│  │ [📄 PDF Report] [📊 Excel] [📋 CSV] [🖼️ PNG Charts] [🔄 Reuse Result] │ │
│  └───────────────────────────────────────────────────────────────────────┘ │
│                                                                             │
├─────────────────────────────────────────────────────────────────────────────┤
│ 💡 Help: CF(∞) = equilibrium concentration. t(95%) = time to reach 95% of  │
│    equilibrium. Click [🔄 Reuse Result] to use this as input for next job. │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

## 9. Tab 6: 📋 Job History

### 9.1 Purpose
Browse, edit, relaunch, compare previous jobs.

### 9.2 Layout

```
┌─────────────────────────────────────────────────────────────────────────────┐
│ 📋 Job History                                                              │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│  Search: [________________] [🔍]    Filter: [All ▼]    Sort: [Date ▼]      │
│                                                                             │
│  ┌───────────────────────────────────────────────────────────────────────┐ │
│  │ ☐ │ Job Name                      │ Date       │ Status │ Actions    │ │
│  ├───┼───────────────────────────────┼────────────┼────────┼────────────┤ │
│  │ ☐ │ sandwich_3step_irganox        │ 2024-12-24 │ ✅     │ [👁️][✏️][🔄]│ │
│  │ ☐ │ bottle_rPP_toluene            │ 2024-12-23 │ ✅     │ [👁️][✏️][🔄]│ │
│  │ ☐ │ tray_LDPE_antioxidants        │ 2024-12-22 │ ✅     │ [👁️][✏️][🔄]│ │
│  │ ☐ │ film_multilayer_test          │ 2024-12-21 │ ⚠️     │ [👁️][✏️][🔄]│ │
│  └───────────────────────────────────────────────────────────────────────┘ │
│                                                                             │
│  Selected: 0    [Compare Selected] [Export Selected] [Delete Selected]     │
│                                                                             │
│  Legend: 👁️ View  ✏️ Edit & Rerun  🔄 Use as Template                       │
│                                                                             │
├─────────────────────────────────────────────────────────────────────────────┤
│ 💡 Help: All simulations are saved for traceability. Use [✏️] to modify    │
│    parameters and rerun. Use [🔄] to create a new job from this template.  │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

## 10. Tab 7: ⚙️ Configuration

### 10.1 Purpose
Global settings for solver, estimators, user profile, and defaults.

### 10.2 Layout

```
┌─────────────────────────────────────────────────────────────────────────────┐
│ ⚙️ Configuration                                                            │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│  ┌───────────────────────────────────────────────────────────────────────┐ │
│  │ 👤 User Profile (for traceability & reports)                    [?]   │ │
│  │ ─────────────────────────────────────────────────────────────────────  │ │
│  │                                                                        │ │
│  │ Name:         [Olivier Vitrac                    ]                    │ │
│  │ Email:        [olivier.vitrac@gmail.com          ]                    │ │
│  │ Organization: [INRAE + Generative Simulation     ]                    │ │
│  │ Role:         [Researcher                     ▼]                      │ │
│  │                                                                        │ │
│  │ 💡 This information appears in exported reports and job metadata      │ │
│  │    for full traceability.                                             │ │
│  └───────────────────────────────────────────────────────────────────────┘ │
│                                                                             │
│  ┌───────────────────────────────────────────────────────────────────────┐ │
│  │ 🔧 Solver Settings                                              [?]   │ │
│  │ ─────────────────────────────────────────────────────────────────────  │ │
│  │                                                                        │ │
│  │ Spatial discretization:  [1800   ] points                             │ │
│  │ Time steps:              [1000   ] steps                              │ │
│  │ Convergence tolerance:   [1e-6   ]                                    │ │
│  │ Max step (adaptive):     [1e-3   ] (0 = auto)                         │ │
│  │                                                                        │ │
│  │ 💡 Higher values = more precision but slower. Default is balanced.    │ │
│  └───────────────────────────────────────────────────────────────────────┘ │
│                                                                             │
│  ┌───────────────────────────────────────────────────────────────────────┐ │
│  │ 📐 Estimator Models                                             [?]   │ │
│  │ ─────────────────────────────────────────────────────────────────────  │ │
│  │                                                                        │ │
│  │ Diffusivity (D):    [Piringer (EU standard)              ▼]          │ │
│  │   Options: Piringer, WLF, Experimental database                       │ │
│  │                                                                        │ │
│  │ Partition (k):      [FHP (Flory-Huggins-Piringer)        ▼]          │ │
│  │   Options: FHP, UNIFAC, Experimental database                         │ │
│  │                                                                        │ │
│  │ 💡 Piringer is the EU regulatory standard. Override per-layer if     │ │
│  │    you have specific experimental values.                             │ │
│  └───────────────────────────────────────────────────────────────────────┘ │
│                                                                             │
│  ┌───────────────────────────────────────────────────────────────────────┐ │
│  │ 🌍 Default Settings                                             [?]   │ │
│  │ ─────────────────────────────────────────────────────────────────────  │ │
│  │                                                                        │ │
│  │ Regulatory framework:  [EU 10/2011                        ▼]          │ │
│  │ Temperature unit:      [°C                                ▼]          │ │
│  │ Thickness unit:        [µm                                ▼]          │ │
│  │ Time unit:             [days                              ▼]          │ │
│  │ Concentration unit:    [mg/kg                             ▼]          │ │
│  │                                                                        │ │
│  │ 💡 These are defaults for new scenarios. Can be overridden per job.   │ │
│  └───────────────────────────────────────────────────────────────────────┘ │
│                                                                             │
│  ┌───────────────────────────────────────────────────────────────────────┐ │
│  │ 💾 Export Settings                                              [?]   │ │
│  │ ─────────────────────────────────────────────────────────────────────  │ │
│  │                                                                        │ │
│  │ Default export folder: [~/sfppy_exports                 ] [Browse]    │ │
│  │ PDF report template:   [Standard (with logo)             ▼]          │ │
│  │ Include in reports:    ☑ Scenario details  ☑ Charts  ☑ Data tables  │ │
│  │                        ☑ Regulatory comparison  ☐ Raw solver output  │ │
│  │                                                                        │ │
│  │ 💡 Reports include full traceability: user, date, versions, machine.  │ │
│  └───────────────────────────────────────────────────────────────────────┘ │
│                                                                             │
│  ┌───────────────────────────────────────────────────────────────────────┐ │
│  │ ⚡ Performance                                                  [?]   │ │
│  │ ─────────────────────────────────────────────────────────────────────  │ │
│  │                                                                        │ │
│  │ Max parallel workers:  [6      ] (available: 12 CPUs)                 │ │
│  │ Job timeout:           [300    ] seconds                              │ │
│  │ Cache PubChem results: ☑ Enabled                                      │ │
│  │                                                                        │ │
│  │ 💡 More workers = faster batch processing but higher memory usage.    │ │
│  └───────────────────────────────────────────────────────────────────────┘ │
│                                                                             │
│  ┌───────────────────────────────────────────────────────────────────────┐ │
│  │ ℹ️ System Information                                                  │ │
│  │ ─────────────────────────────────────────────────────────────────────  │ │
│  │                                                                        │ │
│  │ SFPPy Version:    1.50                                                │ │
│  │ Studio Version:   0.1.0                                               │ │
│  │ Python:           3.11.5                                              │ │
│  │ Machine:          workstation-01 (Linux)                              │ │
│  │ Cache location:   patankar/cache.PubChem/                             │ │
│  │ Jobs stored:      42 jobs (1.2 GB)                                    │ │
│  │                                                                        │ │
│  │ [Clear Cache]  [Export Config]  [Import Config]                       │ │
│  └───────────────────────────────────────────────────────────────────────┘ │
│                                                                             │
│  [Save Configuration]  [Reset to Defaults]                                 │
│                                                                             │
├─────────────────────────────────────────────────────────────────────────────┤
│ 💡 Help: Configuration is saved locally and persists across sessions.      │
│    User profile information appears in all exported reports.               │
└─────────────────────────────────────────────────────────────────────────────┘
```

### 10.3 Sections Summary

| Section | Settings | Help |
|---------|----------|------|
| **👤 User Profile** | Name, email, organization, role | For report headers and traceability |
| **🔧 Solver** | Spatial points, time steps, tolerance, max_step | Higher = more precise, slower |
| **📐 Estimators** | D model (Piringer/WLF), k model (FHP/UNIFAC) | Piringer = EU standard |
| **🌍 Defaults** | Regulatory framework, units (T, L, t, C) | Applied to new scenarios |
| **💾 Export** | Folder, template, report contents | Full traceability in reports |
| **⚡ Performance** | Workers, timeout, caching | Balance speed vs resources |
| **ℹ️ System** | Versions, machine, storage stats | Read-only info |

---

## 11. Help System

### 11.1 Design Principle

**Every action and parameter has clear, contextual help.** Help is provided at three levels:

1. **Inline hints** (💡) — Brief tips visible in the UI
2. **Tooltip help** ([?]) — Detailed explanation on hover/click
3. **Panel help** — Dynamic help panel at the bottom of each tab

### 11.2 Help Icon Behavior

```
┌─────────────────────────────────────┐
│ Diffusivity (D)             [?]    │  ← Click or hover for help
├─────────────────────────────────────┤
│ ┌─────────────────────────────────┐ │
│ │ What is Diffusivity?            │ │  ← Tooltip popup
│ │                                 │ │
│ │ D (m²/s) is the diffusion      │ │
│ │ coefficient that describes how │ │
│ │ fast molecules move through    │ │
│ │ the polymer matrix.            │ │
│ │                                 │ │
│ │ Higher D = faster migration    │ │
│ │ Typical: 1e-16 to 1e-12 m²/s  │ │
│ │                                 │ │
│ │ [📖 Learn more...]             │ │
│ └─────────────────────────────────┘ │
└─────────────────────────────────────┘
```

### 11.3 Help Content by Tab

#### Tab 1: 📦 Assembly

| Element | Help Text |
|---------|-----------|
| Layer count [◀ n ▶] | Number of layers in the packaging. Layer 1 is always in contact with food. |
| Polymer dropdown | Select the polymer type. Properties (D, k) are estimated from polymer characteristics. |
| Thickness | Layer thickness. Common units: µm (films), mm (sheets), cm (thick containers). |
| Glassy checkbox | Check if polymer is below its glass transition temperature (Tg). Glassy polymers act as functional barriers. |
| D override | Override the computed diffusivity with an experimental or literature value. Required when default models don't apply. |
| k override | Override the partition coefficient. Use experimental values for specific polymer/food combinations. |

#### Tab 2: 🍽️ Food

| Element | Help Text |
|---------|-----------|
| Food category | Select the food type. This determines which simulant to use for testing. |
| Texture | Liquid, semisolid, or solid. Affects mass transfer at the surface. |
| Geometry shape | Select 3D shape. Volume and surface area are calculated automatically. |
| Contact time | Duration of food contact. Typical: 10 days (EU test), shelf life for real assessment. |
| Contact temperature | Storage/use temperature. Higher T = faster migration. |
| Multi-step chain | Enable for complex scenarios (set-off → filling → storage). Each step can have different conditions. |
| Without food | Check for set-off simulation (no food, just material equilibration). |

#### Tab 3: ⚗️ Substances

| Element | Help Text |
|---------|-----------|
| Search | Search by name, CAS number, or PubChem CID. Local cache is checked first. |
| C0 (initial concentration) | Amount of substance in the layer before migration (mg/kg polymer). |
| SML | Specific Migration Limit from EU 10/2011. Maximum allowed in food. |
| FCM | Food Contact Material number from EU positive list. |
| TTC | Threshold of Toxicological Concern. Safety limit for non-listed substances. |
| Layer assignment | Check which layers contain this substance. Each layer can have different C0. |

#### Tab 4: 🔬 Simulate

| Element | Help Text |
|---------|-----------|
| Pipeline preview | Shows the symbolic representation of your scenario (substance % food >> layers >> food). |
| Link previous | Continue from a previous simulation's final state. For multi-step chains. |
| Job name | Unique name for this simulation. Used in exports and job history. |
| Solver settings | Advanced: spatial points, time steps, tolerance. Higher = more accurate but slower. |

#### Tab 5: 📊 Results

| Element | Help Text |
|---------|-----------|
| CF(t) | Concentration in food vs time. Main output for compliance assessment. |
| Cx(x) | Concentration profile through material layers at contact time. Shows gradients. |
| Q95 | 95th percentile concentration (probabilistic mode). |
| SML comparison | Red line shows regulatory limit. Below = compliant. |
| Merged steps | For multi-step scenarios, show combined result across all steps. |
| Log scale | Use logarithmic scale for wide concentration ranges. |

### 11.4 Dynamic Help Panel

Each tab has a collapsible help panel at the bottom:

```
┌─────────────────────────────────────────────────────────────────────────────┐
│ 💡 Help                                                          [Expand ▼]│
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│ Current context: Layer 2 (rPP) — Parameter Override                        │
│                                                                             │
│ You are configuring manual overrides for diffusion parameters. This is     │
│ useful when:                                                                │
│                                                                             │
│ • You have experimental D or k values from literature                      │
│ • The Piringer model doesn't apply to your specific polymer/additive       │
│ • You want to perform sensitivity analysis with different values           │
│                                                                             │
│ Recommended values for Irganox 1076 in PP at 25°C:                         │
│   D = 1.2e-14 m²/s (Piringer estimate)                                     │
│   k = 1.0 (FHP estimate)                                                   │
│                                                                             │
│ Reference: Piringer et al., Food Additives & Contaminants, 2008            │
│                                                                             │
│ [📖 Open Documentation]  [📧 Report Issue]                                 │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

### 11.5 Help API

```python
# Backend provides help content
GET /api/help/context/{tab}/{element}

# Example response
{
    "element": "D_override",
    "title": "Diffusivity Override",
    "short": "Override computed D with experimental value",
    "long": "The diffusion coefficient D (m²/s) describes...",
    "tips": ["Use scientific notation: 1.2e-14", "Typical range: 1e-16 to 1e-12"],
    "references": ["Piringer et al. (2008)", "EU JRC Technical Report"],
    "related": ["k_override", "piringer_model"]
}
```

---

## 12. Footer

```html
<footer class="opacity-20 hover:opacity-100 transition-opacity duration-300">
  <div class="text-xs text-gray-500 text-center py-2">
    <span>INRAE</span> · <span>Generative Simulation</span> ·
    <a href="mailto:olivier.vitrac@gmail.com">Olivier Vitrac</a>, PhD, HDR
  </div>
</footer>
```

**Behavior:** Footer is nearly invisible (20% opacity) but becomes fully visible on mouse hover.

---

## 12. Development Roadmap

### Phase 1: Core Foundation (MVP)
- [ ] Project structure and launcher
- [ ] Tab 1: Assembly (multilayer builder)
- [ ] Tab 2: Food & Conditions (single step)
- [ ] Tab 3: Substances (search + assignment)
- [ ] Tab 4: Simulate (basic solver wrapper)
- [ ] Tab 5: Results (CF plot + export)
- [ ] Base template with navigation

### Phase 2: Full Workflow
- [ ] Multi-step scenarios (Tab 2)
- [ ] Job persistence (Tab 6)
- [ ] Link to previous results
- [ ] PDF report generation
- [ ] Compliance assessment

### Phase 3: Advanced Features
- [ ] Tab 7: Configuration
- [ ] Curve fitting (D, k estimation from data)
- [ ] Sensitivity analysis
- [ ] Batch comparison
- [ ] Functional barrier wizard

---

## 13. API-First Architecture

### 13.1 Principle

**All functionalities are exposed as REST endpoints** that can be tested and debugged independently of the GUI. The web interface is a thin client that calls these APIs.

```
┌─────────────────────────────────────────────────────────────────┐
│                        Web Interface                            │
│                    (Jinja2 + JavaScript)                        │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼ REST calls (JSON)
┌─────────────────────────────────────────────────────────────────┐
│                         FastAPI Backend                         │
│  ┌─────────────┐ ┌─────────────┐ ┌─────────────┐ ┌───────────┐ │
│  │ /api/       │ │ /api/       │ │ /api/       │ │ /api/     │ │
│  │ assembly    │ │ substances  │ │ simulation  │ │ jobs      │ │
│  └─────────────┘ └─────────────┘ └─────────────┘ └───────────┘ │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼ Python calls
┌─────────────────────────────────────────────────────────────────┐
│                      SFPPy Core Library                         │
│            (patankar.migration, patankar.layer, etc.)           │
└─────────────────────────────────────────────────────────────────┘
```

### 13.2 Testing Strategy

Each API endpoint can be tested via:
- **curl** / **httpie** from command line
- **Swagger UI** at `/docs` (FastAPI auto-generated)
- **pytest** with `httpx.AsyncClient`
- **Postman** / **Insomnia** collections

### 13.3 Complete API Reference

#### Assembly APIs

| Method | Endpoint | Request | Response | Description |
|--------|----------|---------|----------|-------------|
| GET | `/api/polymers` | - | `{polymers: [...]}` | List predefined polymers |
| GET | `/api/polymers/{code}` | - | `{polymer: {...}}` | Get polymer details |
| POST | `/api/assembly/validate` | `{layers: [...]}` | `{valid, errors, warnings}` | Validate assembly |
| POST | `/api/assembly/estimate` | `{layers: [...], substance: ...}` | `{D, k, estimated_CF}` | Estimate properties |

#### Food & Conditions APIs

| Method | Endpoint | Request | Response | Description |
|--------|----------|---------|----------|-------------|
| GET | `/api/foods` | - | `{foods: [...]}` | List food categories |
| GET | `/api/simulants` | - | `{simulants: [...]}` | List food simulants |
| GET | `/api/conditions` | - | `{conditions: [...]}` | List condition presets |
| POST | `/api/geometry/calculate` | `{shape, dimensions}` | `{volume_m3, surface_m2}` | Calculate V, S |

#### Substance APIs

| Method | Endpoint | Request | Response | Description |
|--------|----------|---------|----------|-------------|
| POST | `/api/substances/search` | `{query, limit}` | `{results: [...]}` | Search PubChem |
| GET | `/api/substances/{cid}` | - | `{substance: {...}}` | Get substance by CID |
| GET | `/api/substances/cas/{cas}` | - | `{substance: {...}}` | Get substance by CAS |
| GET | `/api/substances/{cid}/image` | - | PNG image | Get structure image |
| GET | `/api/substances/{cid}/regulatory` | - | `{eu, us, cn}` | Get regulatory data |

#### Simulation APIs

| Method | Endpoint | Request | Response | Description |
|--------|----------|---------|----------|-------------|
| POST | `/api/simulation/create` | `{assembly, food, conditions, substances}` | `{job_id}` | Create job |
| POST | `/api/simulation/{job_id}/run` | - | `{status: "started"}` | Start simulation |
| GET | `/api/simulation/{job_id}/status` | - | `{progress, eta, current_step}` | Poll progress |
| POST | `/api/simulation/{job_id}/cancel` | - | `{status: "cancelled"}` | Cancel job |
| GET | `/api/simulation/{job_id}/results` | - | `{CF, Cx, times, ...}` | Get results |

#### Job Management APIs

| Method | Endpoint | Request | Response | Description |
|--------|----------|---------|----------|-------------|
| GET | `/api/jobs` | `?limit&offset&status` | `{jobs: [...], total}` | List all jobs |
| GET | `/api/jobs/{id}` | - | `{job: {...}}` | Get job details |
| PUT | `/api/jobs/{id}` | `{name, ...}` | `{job: {...}}` | Update job metadata |
| DELETE | `/api/jobs/{id}` | - | `{success: true}` | Delete job |
| POST | `/api/jobs/{id}/clone` | - | `{new_job_id}` | Clone job as template |
| POST | `/api/jobs/{id}/rerun` | `{modified_params}` | `{new_job_id}` | Rerun with changes |

#### Export APIs

| Method | Endpoint | Request | Response | Description |
|--------|----------|---------|----------|-------------|
| GET | `/api/jobs/{id}/export/pdf` | - | PDF file | Export report as PDF |
| GET | `/api/jobs/{id}/export/xlsx` | - | Excel file | Export data as Excel |
| GET | `/api/jobs/{id}/export/csv` | - | CSV file | Export data as CSV |
| GET | `/api/jobs/{id}/export/json` | - | JSON file | Export full job JSON |
| GET | `/api/jobs/{id}/charts/cf` | `?format=png|svg` | Image | Export CF chart |
| GET | `/api/jobs/{id}/charts/cx` | `?format=png|svg` | Image | Export Cx chart |

#### Configuration APIs

| Method | Endpoint | Request | Response | Description |
|--------|----------|---------|----------|-------------|
| GET | `/api/config` | - | `{solver, estimators, ...}` | Get all config |
| PUT | `/api/config` | `{solver, ...}` | `{config: {...}}` | Update config |
| GET | `/api/config/defaults` | - | `{defaults: {...}}` | Get factory defaults |
| POST | `/api/config/reset` | - | `{config: {...}}` | Reset to defaults |

#### Curve Fitting APIs (Advanced)

| Method | Endpoint | Request | Response | Description |
|--------|----------|---------|----------|-------------|
| POST | `/api/fitting/upload` | `{data: [[t, CF], ...]}` | `{data_id}` | Upload experimental data |
| POST | `/api/fitting/fit` | `{job_id, data_id, params}` | `{D_fitted, k_fitted, r2}` | Fit D and k |
| GET | `/api/fitting/{fit_id}` | - | `{results: {...}}` | Get fitting results |

### 13.4 Example: Testing Simulation via curl

```bash
# 1. Create a simulation job
curl -X POST http://localhost:8002/api/simulation/create \
  -H "Content-Type: application/json" \
  -d '{
    "assembly": {
      "layers": [
        {"polymer": "PET", "thickness_um": 30, "is_glassy": true},
        {"polymer": "PP", "thickness_um": 300}
      ]
    },
    "food": {
      "category": "fatty",
      "simulant": "ethanol50",
      "temperature_C": 25
    },
    "conditions": {
      "steps": [
        {"duration": 10, "unit": "days", "with_food": true}
      ]
    },
    "substances": [
      {"cid": 91597, "name": "Irganox 1076", "layers": {"2": 1000}}
    ]
  }'

# Response: {"success": true, "job_id": "abc123"}

# 2. Run the simulation
curl -X POST http://localhost:8002/api/simulation/abc123/run

# 3. Poll progress
curl http://localhost:8002/api/simulation/abc123/status

# 4. Get results
curl http://localhost:8002/api/simulation/abc123/results

# 5. Export PDF
curl http://localhost:8002/api/jobs/abc123/export/pdf -o report.pdf
```

### 13.5 Swagger Documentation

FastAPI auto-generates interactive documentation at:
- **Swagger UI:** `http://localhost:8002/docs`
- **ReDoc:** `http://localhost:8002/redoc`
- **OpenAPI JSON:** `http://localhost:8002/openapi.json`

---

## 14. Example Coverage Matrix

**Critical Requirement:** The UI must be able to replicate every operation demonstrated in the command-line examples.

### 14.1 Example 1: Monolayer Migration (`example1.py`)

| Step | CLI Code | UI Equivalent | Tab | API Endpoint |
|------|----------|---------------|-----|--------------|
| Define geometry | `Packaging3D('Cylinder', h, r)` | Geometry selector + dimensions | 🍽️ Food | `POST /api/geometry/calculate` |
| Get substance | `migrant("irganox 1076")` | PubChem search | ⚗️ Substances | `POST /api/substances/search` |
| Create layer | `polymer.LDPE(l=(100,"um"), substance=m, C0=5000)` | Layer card with polymer dropdown + thickness + C0 | 📦 Assembly | `POST /api/assembly/validate` |
| Define food | `class sandwich(food.realfood, food.semisolid, food.fat)` | Food category + texture dropdowns | 🍽️ Food | `GET /api/foods` |
| Set conditions | `contacttime=(10,"days"), contacttemperature=(7,"degC")` | Duration + temperature inputs | 🍽️ Food | - |
| Run solver | `senspatankar(LDPE_layer, FOOD_layer)` | [🚀 Run Simulation] button | 🔬 Simulate | `POST /api/simulation/run` |
| Interpolate | `simulation.interpolate_CF(tnew)` | Data table with custom times | 📊 Results | `GET /api/simulation/{id}/results?times=...` |
| Plot CF(t) | `simulation.plotCF()` | Interactive Chart.js plot | 📊 Results | `GET /api/jobs/{id}/charts/cf` |
| Plot Cx(x) | `simulation.plotCx()` | Concentration profile chart | 📊 Results | `GET /api/jobs/{id}/charts/cx` |
| Store results | `CFSimulationContainer.add(simulation)` | Auto-saved to job history | 📋 Jobs | `GET /api/jobs` |
| Export PDF | `print_figure(hfig, folder)` | [📄 PDF Report] button | 📊 Results | `GET /api/jobs/{id}/export/pdf` |
| Export Excel | `.save_as_excel(filename)` | [📊 Excel] button | 📊 Results | `GET /api/jobs/{id}/export/xlsx` |

### 14.2 Example 2: Multilayer with Functional Barrier (`example2.py`)

| Step | CLI Code | UI Equivalent | Tab | API Endpoint |
|------|----------|---------------|-----|--------------|
| Multi-part geometry | Body + neck bottle | Bottle shape with body/neck dimensions | 🍽️ Food | `POST /api/geometry/calculate` |
| Multilayer assembly | `PET + rPP` (2 layers) | Layer selector [◀ 2 ▶], add PET + rPP | 📦 Assembly | `POST /api/assembly/validate` |
| Functional barrier flag | `gPET` (glassy PET) | ☑ Glassy (FB) checkbox on layer 1 | 📦 Assembly | - |
| Compare with/without FB | Run twice, overlay results | Run multiple jobs, select for comparison | 📋 Jobs | `POST /api/jobs/compare` |
| Sensitivity analysis | Loop over FB thickness (2-60 µm) | Batch parameter sweep (Advanced) | ⚙️ Config | `POST /api/sensitivity/run` |
| Color-coded results | `plt.cm.viridis(i/n)` | Auto-assigned colors in comparison chart | 📊 Results | - |

### 14.3 Example 3: Multi-Step Chain (`example3.py`)

| Step | CLI Code | UI Equivalent | Tab | API Endpoint |
|------|----------|---------------|-----|--------------|
| Trilayer ABA | `A + B + A` | Layer selector [◀ 3 ▶], configure each | 📦 Assembly | `POST /api/assembly/validate` |
| Substance injection | `m % medium1` | Assign substance to layer(s) with C0 | ⚗️ Substances | - |
| Multi-step chain | `>> medium1 >> medium2 >> medium3` | ● Multi-step chain, add 3 steps | 🍽️ Food | - |
| Step without food | Set-off period | ☐ Without food (set-off) checkbox | 🍽️ Food | - |
| Different temperatures | 25°C → 85°C → 4°C | Temperature input per step | 🍽️ Food | - |
| Result merging | `sol1 + sol2 + sol3` | ☑ Merged steps toggle in results | 📊 Results | `GET /api/simulation/{id}/results?merged=true` |
| Pipeline preview | `m % food >> ABA >> food1 >> ...` | Pipeline summary text in Simulate tab | 🔬 Simulate | - |

### 14.4 Example 4: Curve Fitting (`example4.py`)

| Step | CLI Code | UI Equivalent | Tab | API Endpoint |
|------|----------|---------------|-----|--------------|
| Create layerLink | `D = layerLink("D", indices=0, values=Dref)` | (Internal - automatic) | - | - |
| Pseudo-experimental | `R.pseudoexperiment(npoints=30)` | [Generate Demo Data] button | ⚙️ Config | `POST /api/fitting/demo-data` |
| Upload real data | Load CSV/Excel | [📤 Upload Data] button | ⚙️ Config | `POST /api/fitting/upload` |
| Fit D and k | Optimization loop | [🔧 Fit Parameters] button | ⚙️ Config | `POST /api/fitting/fit` |
| Show fitted vs measured | Overlay plot | Comparison chart with markers | 📊 Results | `GET /api/fitting/{id}/chart` |
| R² and residuals | Error metrics | Fitting statistics panel | 📊 Results | `GET /api/fitting/{id}` |

### 14.5 Example 5: RAG Demo (`example5_rag.py`)

| Step | CLI Code | UI Equivalent | Tab | API Endpoint |
|------|----------|---------------|-----|--------------|
| Simple simulation | Toluene in ethanol | Standard workflow | All | - |
| Log-scale plot | `ax.set_yscale('log')` | Scale: [Log ▼] dropdown | 📊 Results | - |
| PDF export | `print_figure()` | [📄 PDF Report] | 📊 Results | `GET /api/jobs/{id}/export/pdf` |

### 14.6 Notebook: gui.ipynb

| Widget | CLI Equivalent | UI Equivalent | Tab |
|--------|----------------|---------------|-----|
| Food tree widget | `create_food_tree_widget()` | Food category dropdowns | 🍽️ Food |
| Packaging widget | `create_packaging_widget()` | Geometry shape selector | 🍽️ Food |
| Multi-layer widget | `create_multi_layer_widget()` | Layer cards with [◀ n ▶] | 📦 Assembly |
| Substance widget | `create_substance_widget()` | PubChem search + assignment | ⚗️ Substances |
| Simulation widget | `create_simulation_widget()` | Simulate tab controls | 🔬 Simulate |
| Plot widget | `create_plotmigration_widget()` | Results charts | 📊 Results |

### 14.7 Coverage Checklist

| Feature | Example | Covered | UI Element |
|---------|---------|---------|------------|
| Single layer | Ex1 | ✅ | Layer selector = 1 |
| Multi-layer (2-10) | Ex2, Ex3 | ✅ | Layer selector [◀ n ▶] |
| Functional barrier | Ex2 | ✅ | ☑ Glassy checkbox |
| Plasticized polymer | Ex2 | ✅ | ☑ Plasticized checkbox |
| Single substance | Ex1 | ✅ | Add 1 substance |
| Multiple substances | Ex1 | ✅ | Add multiple substances |
| Substance per layer | Ex3 | ✅ | Layer assignment checkboxes |
| Single step | Ex1, Ex2 | ✅ | ○ Single step radio |
| Multi-step chain | Ex3 | ✅ | ● Multi-step + [+ Add Step] |
| Step without food | Ex3 | ✅ | ☐ Without food checkbox |
| Various geometries | Ex1, Ex2 | ✅ | Shape dropdown |
| Custom dimensions | Ex1, Ex2 | ✅ | Dimension inputs |
| Temperature variation | Ex3 | ✅ | Temp input per step |
| CF(t) plot | All | ✅ | Main results chart |
| Cx(x) plot | Ex1 | ✅ | Profile chart tab |
| Log scale | Ex5 | ✅ | Scale dropdown |
| SML comparison | Ex1 | ✅ | Compliance panel |
| Result overlay | Ex2 | ✅ | Compare selected jobs |
| Result merging | Ex3 | ✅ | ☑ Merged steps toggle |
| Link previous result | Ex3 | ✅ | Job selector in Simulate |
| PDF export | All | ✅ | [📄 PDF Report] |
| Excel export | Ex1 | ✅ | [📊 Excel] |
| CSV export | Ex1 | ✅ | [📋 CSV] |
| Job history | - | ✅ | Jobs tab |
| Edit & rerun | - | ✅ | [✏️] button |
| Clone as template | - | ✅ | [🔄] button |
| Curve fitting | Ex4 | ✅ | Advanced tab |
| Sensitivity analysis | Ex2 | ✅ | Batch sweep (Advanced) |
| Demo data generation | Ex4 | ✅ | [Generate Demo Data] |

### 14.8 Verification Protocol

For each release, run these verification tests:

```bash
# 1. Verify Example 1 can be reproduced via API
curl -X POST http://localhost:8002/api/test/verify-example1

# 2. Verify Example 2 can be reproduced via API
curl -X POST http://localhost:8002/api/test/verify-example2

# 3. Verify Example 3 can be reproduced via API
curl -X POST http://localhost:8002/api/test/verify-example3

# 4. Verify Example 4 can be reproduced via API
curl -X POST http://localhost:8002/api/test/verify-example4

# 5. Full integration test
curl -X POST http://localhost:8002/api/test/verify-all
```

Each verification endpoint will:
1. Create the same scenario as the example
2. Run the simulation
3. Compare results with reference values (within tolerance)
4. Return pass/fail with diff if failed

---

## 15. References

- SFPP3 Design: https://sfpp3-simulation.contactalimentaire.fr/
- SFPPy Notebooks: `notebooks/gui.ipynb`, `notebooks/demo.ipynb`
- Survey Simulator: `survey/` (FastAPI + Jinja2 pattern)
