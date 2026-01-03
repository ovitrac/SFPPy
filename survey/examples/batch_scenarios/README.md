# Batch Scenarios

Pre-built YAML scenario templates for common food packaging migration simulations.

## Available Scenarios

| File | Polymer | Simulant | Use Case |
|------|---------|----------|----------|
| `dairy_PET_pot.yml` | gPET | ethanol50 | Dairy products (yogurt, cream) |
| `water_PET_bottle.yml` | gPET | water | Bottled water |
| `water_HDPE_cap.yml` | HDPE | water | Bottle caps, closures |
| `fatty_PP_tray.yml` | PP | oliveoil | Fatty foods (cheese, meat) |
| `hot_fill_PET.yml` | PET | water3aceticacid | Hot-filled beverages |

## Scenario Format

Each YAML file follows the `Survey.from_scenario()` format:

```yaml
name: "scenario_name"
description: "Scenario description"

packaging:
  polymer: "gPET"           # Polymer type
  thickness: 200            # µm
  temperature: 25           # °C
  surface_area: 600         # cm²
  food_volume: 1000         # mL

contact:
  simulant: "water"         # Food simulant
  temperature: 25           # °C
  time:
    min: 0
    mode: 864000            # seconds (10 days)
    max: 2592000            # seconds (30 days)

substances:
  - identifier: "112-62-9"  # CAS or name
    c0:
      min: 0
      mode: 100             # mg/kg
      max: 500
```

## Running Scenarios

### Single Scenario
```bash
python survey/run_batch.py survey/examples/batch_scenarios/water_PET_bottle.yml
```

### All Scenarios
```bash
python survey/run_batch.py survey/examples/batch_scenarios/
```

### With Custom Output
```bash
python survey/run_batch.py survey/examples/batch_scenarios/ -o my_results/ -w 4
```

## JSON Export Format

The `example_pf_jobs.json` file demonstrates the JSON format used by Survey Simulator:

```json
[
  {
    "name": "water_PET_bottle",
    "polymer": "gPET",
    "thickness_value": 180,
    "simulant": "water",
    "substances": [
      {"identifier": "112-62-9", "c0_mode": 100, "c0_max": 500}
    ]
  }
]
```

---
*Part of SFPPy Survey Module*
