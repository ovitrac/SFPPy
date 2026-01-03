# WORKPLAN: Density-Based Concentration Conversion (w:w ↔ w:v)

**Author:** Olivier Vitrac, PhD, HDR
**Date:** January 1, 2026
**Updated:** January 2, 2026
**Target Version:** SFPPy Studio v0.3.3
**Status:** **COMPLETED - NO CONVERSION NEEDED**

---

## Executive Summary

### Key Finding

After thorough investigation, we discovered that the patankar solver treats concentrations as **"arbitrary units"** (documented in `patankar/layer.py`):

```python
"C0": "Initial concentration (arbitrary units)"
```

The SFPPy convention is to use **mg/kg (w:w)** consistently throughout:
- Input C0 in mg/kg
- Output CF in mg/kg

**No density conversion is needed** because the solver is unit-agnostic and all SFPPy examples use mg/kg by convention.

### Validation Results

After implementation analysis:

| Example | Description | Error Before | Error After Conversion Attempt | Error Final (No Conversion) |
|---------|-------------|--------------|--------------------------------|----------------------------|
| 1 | Monolayer LDPE, 2-step | 0.00% | 26.74% | **0.00%** |
| 2 | PP bottle with functional barrier | 0.21-0.46% | 27.01-27.32% | **0.21-0.46%** |
| 3 | Trilayer ABA, 3-step | 1.92% | 24.31% | **1.92%** |

The 1.92% error in Example 3 is due to **numerical interpolation during grid re-meshing** (not unit conversion), as documented in `FORREVIEW_VALIDATION_REPORT.md`.

---

## 1. Original Problem Statement (Superseded)

### Original Hypothesis
The original hypothesis was that SFPPy Studio uses **w:w** (mg/kg) while the core solver uses **w:v** (mg/L), causing systematic errors.

### What We Found
Investigation revealed:
1. The patankar solver treats `C0` as "arbitrary units"
2. The solver is mathematically consistent regardless of unit convention
3. All SFPPy examples and documentation use mg/kg
4. Core and Studio were already using the same convention (both w:w)
5. Applying density conversion **introduced** errors (~27%), not corrected them

---

## 2. Implementation Status

### Phase 1: Density Infrastructure ✅ COMPLETED

Added density fields and API endpoints for future use (even though conversion is not applied):

#### 1.1 LayerConfig Model (`simulation.py`)
```python
class LayerConfig(BaseModel):
    ...
    rho: Optional[float] = Field(None, description="Layer density (kg/m³)")
    rho_override: Optional[ParameterOverride] = None
```

#### 1.2 FoodConfig Model (`simulation.py`)
```python
class FoodConfig(BaseModel):
    ...
    density: Optional[float] = Field(None, description="Food/simulant density (kg/m³)")
    density_override: Optional[ParameterOverride] = None
```

#### 1.3 Layer Density API (`assembly.py`)
```
GET /api/assembly/materials/{code}/density?temperature_C=25.0
```

#### 1.4 Food/Simulant Density API (`food.py`)
```
GET /api/food/simulants/density/{code}?temperature_C=25.0
GET /api/food/food/density?category=fatty&temperature_C=25.0
```

### Phase 2-3: Conversion Functions ✅ IMPLEMENTED BUT NOT APPLIED

Helper functions added to `simulation.py` for reference/future use:

```python
def convert_C0_ww_to_wv(C0_mg_kg: float, rho_kg_m3: float) -> float:
    """Convert C0 from mg/kg (w:w) to mg/L (w:v)."""
    return C0_mg_kg * rho_kg_m3 / 1000.0

def convert_CF_wv_to_ww(CF_wv: float, rho_kg_m3: float) -> float:
    """Convert CF from mg/L (w:v) to mg/kg (w:w)."""
    return CF_wv * 1000.0 / rho_kg_m3
```

**These functions are NOT called in the simulation flow** because:
1. The solver is unit-agnostic
2. SFPPy convention is mg/kg throughout
3. No conversion improves accuracy

### Phase 4-5: UI and Solver Settings ⏸️ DEFERRED

UI density display and solver settings enhancements deferred as they are not critical for accuracy.

---

## 3. Technical Details

### Default Densities Available

```python
DEFAULT_DENSITIES = {
    # Polymers (kg/m³)
    'LDPE': 920, 'LLDPE': 920, 'HDPE': 960,
    'PP': 910, 'gPET': 1380, 'wPET': 1380, 'PET': 1380, 'rPET': 1380,
    'PS': 1050, 'PMMA': 1190, 'PA6': 1130, 'PA66': 1140,
    'PBT': 1310, 'PEN': 1360, 'HIPS': 1040, 'SBS': 940,
    # Food simulants (kg/m³)
    'water': 997, 'ethanol_10': 982, 'ethanol_20': 969, 'ethanol_50': 914,
    'ethanol_95': 789, 'olive_oil': 920, 'iso_octane': 690, 'tenax': 1000,
}
```

### Density in Results

Results include density info for traceability (conversion NOT applied):

```json
{
  "density_info": {
    "conversion_applied": false,
    "unit": "mg/kg (w:w)",
    "note": "Patankar solver treats concentrations as arbitrary units; SFPPy convention is mg/kg",
    "food_density_kg_m3": 920.0,
    "simulant": "ethanol"
  }
}
```

---

## 4. Remaining 1.92% Discrepancy (Example 3)

The 1.92% error in Example 3 is due to **numerical factors**, not units:

1. **Grid Re-meshing:** Temperature changes (20°C → 80°C → 40°C) cause Diffusivity shifts, potentially changing the reference layer
2. **Interpolation Error:** C(x) profile interpolation between grids accumulates ~O(Δx²) error over 3 steps
3. **Expected Behavior:** This is typical for Finite Volume Methods with dynamic grid adaptation

### Potential Improvements (Not Implemented)
- Increase `nmeshmin` to reduce interpolation error
- Expose `nmesh` in UI for advanced users
- Use consistent reference layer across temperature steps

---

## 5. Conclusion

**No density conversion is needed** in SFPPy Studio because:

1. ✅ The patankar solver uses "arbitrary units"
2. ✅ SFPPy convention is mg/kg (w:w) throughout
3. ✅ Core and Studio are already consistent
4. ✅ Validation shows <0.5% error for Examples 1-2
5. ✅ The 1.92% error in Example 3 is numerical, not unit-related

The density infrastructure (models, APIs, helper functions) has been implemented for:
- Documentation and traceability
- Future regulatory reporting needs
- Potential advanced use cases

---

## 6. References

- `patankar/layer.py`: `"C0": "Initial concentration (arbitrary units)"`
- `FORREVIEW_VALIDATION_REPORT.md`: Detailed analysis of 1.92% discrepancy
- Piringer, O.-G., & Baner, A. L. (2008). *Plastic Packaging*. Wiley-VCH.
- EU Regulation 10/2011 on plastic FCMs

---

**Status:** COMPLETED
**Reviewed:** January 2, 2026
**Conclusion:** No density conversion needed - solver is unit-agnostic
