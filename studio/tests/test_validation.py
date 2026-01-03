#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SFPPy Studio Validation Tests

Compare Studio API simulation results against SFPPy core examples.

Test cases:
1. example1_extensions.py - Monolayer LDPE with Irganox 1076 (2-step)
2. example2.py - Multilayer PP bottle with toluene (functional barrier)
3. example3_shortvariant.py - Trilayer ABA with limonene (3-step multi-temp)

All dimensions are stored in SI units (meters) internally.

@author: SFPPy Studio Validation Suite
@date: January 2026
"""
import sys
import json
import numpy as np
from pathlib import Path

# Setup paths
ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "studio" / "app"))

# Results storage
RESULTS = {
    "example1_extensions": {"core": {}, "studio": {}, "error": {}},
    "example2": {"core": {}, "studio": {}, "error": {}},
    "example3_shortvariant": {"core": {}, "studio": {}, "error": {}},
}

TOLERANCE_PERCENT = 5.0  # Acceptable error threshold

# ==============================================================================
# Test 1: example1_extensions.py
# ==============================================================================
def test_example1_extensions():
    """
    Example 1 Extensions: Monolayer LDPE with Irganox 1076
    - Step 1: 10 days at 7°C
    - Step 2: +4 hours at 25°C (warming)
    """
    print("\n" + "=" * 70)
    print("TEST 1: example1_extensions.py - Monolayer LDPE, 2-step")
    print("=" * 70)

    # --- CORE SIMULATION ---
    from patankar.layer import LDPE
    from patankar.loadpubchem import migrant
    from patankar.geometry import Packaging3D
    import patankar.food as food

    # Geometry: Cylinder 19cm height, 30mm radius
    internal_volume, surface_area = Packaging3D(
        "Cylinder", length=(19, "cm"), radius=(30, "mm")
    ).get_volume_and_area()

    m = migrant("irganox 1076")
    film = LDPE(migrant=m, l=(100, "um"))

    class Sandwich(food.realfood, food.semisolid, food.fat):
        name = "sandwich"

    F = Sandwich(
        volume=internal_volume,
        surfacearea=surface_area,
        contacttime=(10, "days"),
        contacttemperature=(7, "C"),
        substance=m,
        simulant="ethanol"
    )

    # Step 1: 10 days @ 7°C
    kin = F.migration(film)
    CF1_core = float(kin.CFtarget)

    # Step 2: +4h @ 25°C
    F2 = F.copy().update(contacttime=(4, "hours"), contacttemperature=(25, "C"))
    kin2 = kin >> F2
    CF2_core = float(kin2.CFtarget)

    print(f"  CORE: Step 1 (10d @ 7°C):  CF = {CF1_core:.4f} mg/kg")
    print(f"  CORE: Step 2 (+4h @ 25°C): CF = {CF2_core:.4f} mg/kg")

    RESULTS["example1_extensions"]["core"] = {
        "step1_CF": CF1_core,
        "step2_CF": CF2_core,
    }

    # --- STUDIO SIMULATION ---
    # Run step 1 alone first
    from routes.simulation import execute_real_simulation, SimulationInput, get_sfppy_modules, apply_solver_settings

    sfppy = get_sfppy_modules()
    apply_solver_settings(sfppy)

    # Step 1 only
    config_step1 = SimulationInput(
        name="Example 1 - Step 1",
        layers=[{
            "index": 1,
            "polymer": "LDPE",
            "thickness": 100,
            "thickness_unit": "um",
        }],
        substances=[{
            "id": "irganox1076",
            "name": "irganox 1076",
            "SML": 6.0,
            "layer_assignments": {1: {"C0": 1000}},  # Default C0
        }],
        steps=[
            {"index": 1, "temperature_C": 7, "duration": 10, "duration_unit": "days", "with_food": True},
        ],
        geometry={"shape": "cylinder", "dimensions": {"radius": 30, "height": 190}},  # mm
        food={"category": "fatty", "simulant": "ethanol"},
    )

    result_step1 = execute_real_simulation(config_step1, sfppy=sfppy)
    CF1_studio = result_step1.get("substances", [{}])[0].get("CF_at_tcontact", 0) if result_step1.get("substances") else 0

    # Full 2-step simulation
    config_full = SimulationInput(
        name="Example 1 Extensions",
        layers=[{
            "index": 1,
            "polymer": "LDPE",
            "thickness": 100,
            "thickness_unit": "um",
        }],
        substances=[{
            "id": "irganox1076",
            "name": "irganox 1076",
            "SML": 6.0,
            "layer_assignments": {1: {"C0": 1000}},  # Default C0
        }],
        steps=[
            {"index": 1, "temperature_C": 7, "duration": 10, "duration_unit": "days", "with_food": True},
            {"index": 2, "temperature_C": 25, "duration": 4, "duration_unit": "hours", "with_food": True},
        ],
        geometry={"shape": "cylinder", "dimensions": {"radius": 30, "height": 190}},  # mm
        food={"category": "fatty", "simulant": "ethanol"},
    )

    result_full = execute_real_simulation(config_full, sfppy=sfppy)
    CF2_studio = result_full.get("substances", [{}])[0].get("CF_at_tcontact", 0) if result_full.get("substances") else 0

    print(f"  STUDIO: Step 1 (10d @ 7°C):  CF = {CF1_studio:.4f} mg/kg")
    print(f"  STUDIO: After Step 2 (+4h @ 25°C): CF = {CF2_studio:.4f} mg/kg")

    RESULTS["example1_extensions"]["studio"] = {
        "step1_CF": CF1_studio,
        "step2_CF": CF2_studio,
    }

    # Calculate errors
    err1 = abs(CF1_core - CF1_studio) / CF1_core * 100 if CF1_core > 0 else 0
    err2 = abs(CF2_core - CF2_studio) / CF2_core * 100 if CF2_core > 0 else 0

    RESULTS["example1_extensions"]["error"] = {
        "step1_error_pct": err1,
        "step2_error_pct": err2,
    }

    print(f"  ERROR: Step 1 = {err1:.2f}%, Step 2 = {err2:.2f}%")
    return max(err1, err2) < TOLERANCE_PERCENT


# ==============================================================================
# Test 2: example2.py
# ==============================================================================
def test_example2():
    """
    Example 2: PP bottle with toluene, functional barrier study
    - 300µm PP with 30µm PET functional barrier
    - 450 days at 20°C
    """
    print("\n" + "=" * 70)
    print("TEST 2: example2.py - PP bottle with functional barrier")
    print("=" * 70)

    # --- CORE SIMULATION ---
    from patankar.loadpubchem import migrant
    from patankar.geometry import Packaging3D
    import patankar.food as food
    import patankar.layer as polymer
    from patankar.migration import senspatankar as solver

    # Bottle geometry
    bottle = Packaging3D(
        "bottle",
        body_radius=(40, "mm"),
        body_height=(0.2, "m"),
        neck_radius=(1.8, "cm"),
        neck_height=0.05  # Default SI unit: meters
    )
    internalvolume, contactsurface = bottle.get_volume_and_area()

    surrogate = migrant("toluene")
    contactTemperature = (20, "degC")
    contactTime = (450, "days")
    maxConcentration = 10

    # PP walls with toluene
    PPwalls_with_toluene = polymer.PP(
        l=(300, "um"),
        substance=surrogate,
        C0=maxConcentration,
        T=contactTemperature
    )

    # Liquid food
    class liquidFood(food.realfood, food.liquid, food.fat):
        name = "liquidFood"

    FOODlayer = liquidFood(
        volume=internalvolume,
        surfacearea=contactsurface,
        contacttime=contactTime,
        contacttemperature=contactTemperature
    )

    # Without FB
    ref_simulation = solver(PPwalls_with_toluene, FOODlayer, name="bottle-rPP")
    # Use CFtarget for consistency with Studio (concentration at contact time)
    CF_noFB_core = float(ref_simulation.CFtarget)

    # With 30µm PET FB
    PET_functionalBarrier = polymer.wPET(
        l=(30, "um"),
        substance=surrogate,
        C0=0,
        T=contactTemperature
    )
    FBwalls = PET_functionalBarrier + PPwalls_with_toluene
    fb_simulation = solver(FBwalls, FOODlayer, name="bottleFB-PET-rPP")
    # Use CFtarget for consistency with Studio (concentration at contact time)
    CF_withFB_core = float(fb_simulation.CFtarget)

    print(f"  CORE: Without FB: CF = {CF_noFB_core:.4f} mg/kg")
    print(f"  CORE: With 30µm PET FB: CF = {CF_withFB_core:.4f} mg/kg")

    RESULTS["example2"]["core"] = {
        "without_FB_CF": CF_noFB_core,
        "with_FB_CF": CF_withFB_core,
    }

    # --- STUDIO SIMULATION (without FB) ---
    from routes.simulation import execute_real_simulation, SimulationInput, get_sfppy_modules, apply_solver_settings

    sfppy = get_sfppy_modules()
    apply_solver_settings(sfppy)

    # Without FB
    config_noFB = SimulationInput(
        name="Example 2 - No FB",
        layers=[{
            "index": 1,
            "polymer": "PP",
            "thickness": 300,
            "thickness_unit": "um",
        }],
        substances=[{
            "id": "toluene",
            "name": "toluene",
            "SML": 60.0,
            "layer_assignments": {1: {"C0": 10}},
        }],
        steps=[
            {"index": 1, "temperature_C": 20, "duration": 450, "duration_unit": "days", "with_food": True},
        ],
        geometry={"shape": "bottle", "dimensions": {
            "body_radius": 40, "body_height": 200, "neck_radius": 18, "neck_height": 50
        }},  # mm
        food={"category": "fatty liquid", "simulant": "ethanol"},
    )

    result_noFB = execute_real_simulation(config_noFB, sfppy=sfppy)
    CF_noFB_studio = result_noFB.get("substances", [{}])[0].get("CF_at_tcontact", 0) if result_noFB.get("substances") else 0

    # With FB
    config_withFB = SimulationInput(
        name="Example 2 - With FB",
        layers=[
            {"index": 1, "polymer": "wPET", "thickness": 30, "thickness_unit": "um"},
            {"index": 2, "polymer": "PP", "thickness": 300, "thickness_unit": "um"},
        ],
        substances=[{
            "id": "toluene",
            "name": "toluene",
            "SML": 60.0,
            "layer_assignments": {1: {"C0": 0}, 2: {"C0": 10}},
        }],
        steps=[
            {"index": 1, "temperature_C": 20, "duration": 450, "duration_unit": "days", "with_food": True},
        ],
        geometry={"shape": "bottle", "dimensions": {
            "body_radius": 40, "body_height": 200, "neck_radius": 18, "neck_height": 50
        }},
        food={"category": "fatty liquid", "simulant": "ethanol"},
    )

    result_withFB = execute_real_simulation(config_withFB, sfppy=sfppy)
    CF_withFB_studio = result_withFB.get("substances", [{}])[0].get("CF_at_tcontact", 0) if result_withFB.get("substances") else 0

    print(f"  STUDIO: Without FB: CF = {CF_noFB_studio:.4f} mg/kg")
    print(f"  STUDIO: With 30µm PET FB: CF = {CF_withFB_studio:.4f} mg/kg")

    RESULTS["example2"]["studio"] = {
        "without_FB_CF": CF_noFB_studio,
        "with_FB_CF": CF_withFB_studio,
    }

    # Calculate errors
    err_noFB = abs(CF_noFB_core - CF_noFB_studio) / CF_noFB_core * 100 if CF_noFB_core > 0 else 0
    err_withFB = abs(CF_withFB_core - CF_withFB_studio) / CF_withFB_core * 100 if CF_withFB_core > 0 else 0

    RESULTS["example2"]["error"] = {
        "without_FB_error_pct": err_noFB,
        "with_FB_error_pct": err_withFB,
    }

    print(f"  ERROR: Without FB = {err_noFB:.2f}%, With FB = {err_withFB:.2f}%")
    return max(err_noFB, err_withFB) < TOLERANCE_PERCENT


# ==============================================================================
# Test 3: example3_shortvariant.py
# ==============================================================================
def test_example3_shortvariant():
    """
    Example 3 Short Variant: Trilayer ABA with limonene
    - wPET(30µm) + PP(500µm) + gPET(30µm)
    - Step 1: Setoff 4 months @ 20°C
    - Step 2: Hot-fill (default ~2min @ 85°C)
    - Step 3: Storage 6 months @ 30°C
    """
    print("\n" + "=" * 70)
    print("TEST 3: example3_shortvariant.py - Trilayer ABA, 3-step")
    print("=" * 70)

    # --- CORE SIMULATION ---
    from patankar.loadpubchem import migrant
    from patankar.layer import gPET, wPET, PP
    from patankar.geometry import Packaging3D
    from patankar.food import setoff, ambient, hotfilled, realfood, liquid, fat

    # Geometry
    box = Packaging3D("box_container", length=(19, "cm"), width=(10, "cm"), height=(8, "cm"))
    Vint, Acontact = box.get_volume_and_area()

    # Migrant and layers
    m = migrant("limonene")
    A1 = wPET(l=(30, "um"), migrant=m, C0=0)
    B = PP(l=(0.5, "mm"), migrant=m, CP0=200)
    A2 = gPET(l=(30, "um"), migrant=m, C0=0)
    ABA = A1 + B + A2

    # Contact conditions
    class Step1(setoff, ambient):
        contacttime = (4, "months")
        contacttemperature = (20, "C")

    class Step2(hotfilled, realfood, liquid, fat):
        pass

    class Step3(ambient, realfood, liquid, fat):
        contacttime = (6, "months")
        contacttemperature = (30, "C")

    F1, F2, F3 = Step1(), Step2(), Step3()
    box >> F1 >> F2 >> F3  # Propagate geometry

    # Chained simulation
    F1 >> ABA >> F1 >> F2 >> F3
    sol_ref = F1.lastsimulation + F2.lastsimulation + F3.lastsimulation

    # Use CFtarget for consistency with Studio (concentration at contact time)
    CF_core = float(sol_ref.CFtarget)
    print(f"  CORE: Final CF = {CF_core:.4f} mg/kg")

    # Get per-step CFtarget values for consistency
    CF1_core = float(F1.lastsimulation.CFtarget) if hasattr(F1.lastsimulation, 'CFtarget') and F1.lastsimulation.CFtarget is not None else 0
    CF2_core = float(F2.lastsimulation.CFtarget) if hasattr(F2.lastsimulation, 'CFtarget') and F2.lastsimulation.CFtarget is not None else 0
    CF3_core = float(F3.lastsimulation.CFtarget) if hasattr(F3.lastsimulation, 'CFtarget') and F3.lastsimulation.CFtarget is not None else 0

    print(f"    Step 1 (setoff): CF = {CF1_core:.4f}")
    print(f"    Step 2 (hot-fill): CF = {CF2_core:.4f}")
    print(f"    Step 3 (storage): CF = {CF3_core:.4f}")

    RESULTS["example3_shortvariant"]["core"] = {
        "final_CF": CF_core,
        "step1_CF": CF1_core,
        "step2_CF": CF2_core,
        "step3_CF": CF3_core,
    }

    # --- STUDIO SIMULATION ---
    from routes.simulation import execute_real_simulation, SimulationInput, get_sfppy_modules, apply_solver_settings

    sfppy = get_sfppy_modules()
    apply_solver_settings(sfppy)

    config = SimulationInput(
        name="Example 3 Short Variant",
        layers=[
            {"index": 1, "polymer": "wPET", "thickness": 30, "thickness_unit": "um"},
            {"index": 2, "polymer": "PP", "thickness": 500, "thickness_unit": "um"},
            {"index": 3, "polymer": "gPET", "thickness": 30, "thickness_unit": "um"},
        ],
        substances=[{
            "id": "limonene",
            "name": "limonene",
            "SML": 60.0,
            "layer_assignments": {1: {"C0": 0}, 2: {"C0": 200}, 3: {"C0": 0}},
        }],
        steps=[
            {"index": 1, "temperature_C": 20, "duration": 4, "duration_unit": "months", "with_food": False, "setoff_type": "stacked"},
            {"index": 2, "temperature_C": 80, "duration": 20, "duration_unit": "minutes", "with_food": True},  # hotfilled default: 20min @ 80°C
            {"index": 3, "temperature_C": 30, "duration": 6, "duration_unit": "months", "with_food": True},
        ],
        geometry={"shape": "box_container", "dimensions": {"length": 190, "width": 100, "height": 80}},  # mm
        food={"category": "fatty liquid", "simulant": "ethanol"},
    )

    result = execute_real_simulation(config, sfppy=sfppy)
    substances = result.get("substances", [])

    if substances:
        CF_studio = substances[0].get("CF_at_tcontact", 0)
        steps = substances[0].get("steps", [])
        CF1_studio = steps[0]["CF_final"] if len(steps) > 0 else 0
        CF2_studio = steps[1]["CF_final"] if len(steps) > 1 else 0
        CF3_studio = steps[2]["CF_final"] if len(steps) > 2 else 0
    else:
        CF_studio = CF1_studio = CF2_studio = CF3_studio = 0

    print(f"  STUDIO: Final CF = {CF_studio:.4f} mg/kg")
    print(f"    Step 1 (setoff): CF = {CF1_studio:.4f}")
    print(f"    Step 2 (hot-fill): CF = {CF2_studio:.4f}")
    print(f"    Step 3 (storage): CF = {CF3_studio:.4f}")

    RESULTS["example3_shortvariant"]["studio"] = {
        "final_CF": CF_studio,
        "step1_CF": CF1_studio,
        "step2_CF": CF2_studio,
        "step3_CF": CF3_studio,
    }

    # Calculate error
    err_final = abs(CF_core - CF_studio) / CF_core * 100 if CF_core > 0 else 0

    RESULTS["example3_shortvariant"]["error"] = {
        "final_error_pct": err_final,
    }

    print(f"  ERROR: Final CF = {err_final:.2f}%")
    return err_final < TOLERANCE_PERCENT


# ==============================================================================
# Main
# ==============================================================================
def main():
    print("\n" + "=" * 70)
    print("SFPPy STUDIO VALIDATION SUITE")
    print("Comparing Studio API results against SFPPy Core examples")
    print("=" * 70)

    results = {}

    try:
        results["example1_extensions"] = test_example1_extensions()
    except Exception as e:
        print(f"  EXCEPTION: {e}")
        import traceback
        traceback.print_exc()
        results["example1_extensions"] = False

    try:
        results["example2"] = test_example2()
    except Exception as e:
        print(f"  EXCEPTION: {e}")
        import traceback
        traceback.print_exc()
        results["example2"] = False

    try:
        results["example3_shortvariant"] = test_example3_shortvariant()
    except Exception as e:
        print(f"  EXCEPTION: {e}")
        import traceback
        traceback.print_exc()
        results["example3_shortvariant"] = False

    # Summary
    print("\n" + "=" * 70)
    print("VALIDATION SUMMARY")
    print("=" * 70)

    all_pass = True
    for test_name, passed in results.items():
        status = "✅ PASS" if passed else "❌ FAIL"
        print(f"  {test_name}: {status}")
        if not passed:
            all_pass = False

    print("\n" + "-" * 70)
    if all_pass:
        print("  OVERALL: ✅ ALL TESTS PASSED")
    else:
        print("  OVERALL: ❌ SOME TESTS FAILED")

    # Output detailed results for documentation
    print("\n" + "=" * 70)
    print("DETAILED RESULTS (for README.md)")
    print("=" * 70)

    for test_name, data in RESULTS.items():
        print(f"\n### {test_name}")
        print(f"Core results: {data['core']}")
        print(f"Studio results: {data['studio']}")
        print(f"Errors: {data['error']}")

    return all_pass


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
