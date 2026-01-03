#!/usr/bin/env python3
"""
SFPPy Survey Module - Dependency Checker
=========================================

Check that all required and optional dependencies are installed.

Usage:
    python survey/check_dependencies.py
    python survey/check_dependencies.py --install  # Show install commands
    python survey/check_dependencies.py --test     # Run basic functionality tests

Author: Olivier Vitrac
"""

import sys
import importlib
from pathlib import Path
from typing import Tuple, Optional


def check_import(module_name: str) -> Tuple[bool, Optional[str]]:
    """Check if a module can be imported, return (success, version)."""
    try:
        mod = importlib.import_module(module_name)
        version = getattr(mod, "__version__", None)
        return True, version
    except ImportError:
        return False, None


def print_status(name: str, installed: bool, version: Optional[str] = None,
                 required: bool = True, description: str = "") -> None:
    """Print formatted status line."""
    if installed:
        ver_str = f" ({version})" if version else ""
        print(f"  [{'OK':^8}] {name}{ver_str}")
    else:
        status = "MISSING" if required else "OPTIONAL"
        desc = f" - {description}" if description else ""
        print(f"  [{status:^8}] {name}{desc}")


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Check SFPPy Survey dependencies")
    parser.add_argument("--install", action="store_true",
                        help="Show pip install commands for missing packages")
    parser.add_argument("--test", action="store_true",
                        help="Run basic functionality tests")
    parser.add_argument("--quiet", "-q", action="store_true",
                        help="Only show errors")
    args = parser.parse_args()

    if not args.quiet:
        print("=" * 70)
        print("SFPPy Survey Module - Dependency Check")
        print("=" * 70)

    missing_required = []
    missing_optional = []
    all_packages = []

    # =========================================================================
    # Core Python Dependencies
    # =========================================================================
    if not args.quiet:
        print("\n[1/4] Core Scientific Libraries")
        print("-" * 40)

    core_deps = [
        ("numpy", "numpy", True, "Numerical computing"),
        ("scipy", "scipy", True, "Scientific computing"),
        ("yaml", "pyyaml", True, "YAML parsing"),
    ]

    for import_name, pip_name, required, desc in core_deps:
        installed, version = check_import(import_name)
        all_packages.append((pip_name, installed, required))
        if not args.quiet:
            print_status(import_name, installed, version, required, desc)
        if not installed and required:
            missing_required.append(pip_name)

    # =========================================================================
    # Web Framework Dependencies
    # =========================================================================
    if not args.quiet:
        print("\n[2/4] Web Framework Libraries")
        print("-" * 40)

    web_deps = [
        ("fastapi", "fastapi", True, "Web API framework"),
        ("uvicorn", "uvicorn", True, "ASGI server"),
        ("pydantic", "pydantic", True, "Data validation"),
        ("jinja2", "jinja2", True, "HTML templating"),
        ("starlette", "starlette", True, "HTTP toolkit"),
    ]

    for import_name, pip_name, required, desc in web_deps:
        installed, version = check_import(import_name)
        all_packages.append((pip_name, installed, required))
        if not args.quiet:
            print_status(import_name, installed, version, required, desc)
        if not installed and required:
            missing_required.append(pip_name)

    # =========================================================================
    # Spreadsheet Dependencies
    # =========================================================================
    if not args.quiet:
        print("\n[3/4] Spreadsheet Support (Optional)")
        print("-" * 40)

    spreadsheet_deps = [
        ("openpyxl", "openpyxl", False, "XLSX file support"),
        ("odf", "odfpy", False, "ODS file support"),
    ]

    for import_name, pip_name, required, desc in spreadsheet_deps:
        installed, version = check_import(import_name)
        all_packages.append((pip_name, installed, required))
        if not args.quiet:
            print_status(import_name, installed, version, required, desc)
        if not installed and not required:
            missing_optional.append((pip_name, desc))

    # =========================================================================
    # Integration Dependencies
    # =========================================================================
    if not args.quiet:
        print("\n[4/4] External Integrations (Optional)")
        print("-" * 40)

    integration_deps = [
        ("requests", "requests", False, "PubChem API access"),
        ("aiofiles", "aiofiles", False, "Async file operations"),
    ]

    for import_name, pip_name, required, desc in integration_deps:
        installed, version = check_import(import_name)
        all_packages.append((pip_name, installed, required))
        if not args.quiet:
            print_status(import_name, installed, version, required, desc)
        if not installed and not required:
            missing_optional.append((pip_name, desc))

    # =========================================================================
    # Check SFPPy modules
    # =========================================================================
    if not args.quiet:
        print("\n" + "=" * 70)
        print("SFPPy Module Check")
        print("=" * 70)

    sfppy_modules = [
        ("survey.survey", "Survey class"),
        ("survey.models", "Data models"),
        ("survey.io", "I/O utilities"),
        ("survey.spreadsheet", "Spreadsheet handler"),
        ("survey.batch", "Batch runner"),
        ("survey.app.main", "Family Editor app"),
        ("survey.app.simulator", "Survey Simulator app"),
        ("patankar", "Patankar solver"),
    ]

    sfppy_ok = True
    for module_name, desc in sfppy_modules:
        installed, _ = check_import(module_name)
        if not args.quiet:
            status = "OK" if installed else "ERROR"
            print(f"  [{status:^8}] {module_name} - {desc}")
        if not installed:
            sfppy_ok = False

    # =========================================================================
    # Summary
    # =========================================================================
    print("\n" + "=" * 70)
    print("Summary")
    print("=" * 70)

    if missing_required:
        print(f"\n  MISSING REQUIRED: {len(missing_required)} packages")
        for pkg in missing_required:
            print(f"    - {pkg}")
    else:
        print("\n  All required dependencies are installed.")

    if missing_optional:
        print(f"\n  MISSING OPTIONAL: {len(missing_optional)} packages")
        for pkg, desc in missing_optional:
            print(f"    - {pkg}: {desc}")

    if not sfppy_ok:
        print("\n  WARNING: Some SFPPy modules failed to import.")
        print("  Make sure SFPPy is installed: pip install -e .")

    # =========================================================================
    # Install commands
    # =========================================================================
    if args.install or missing_required:
        print("\n" + "=" * 70)
        print("Installation Commands")
        print("=" * 70)

        if missing_required:
            cmd = f"pip install {' '.join(missing_required)}"
            print(f"\n  Required packages:\n    {cmd}")

        if missing_optional:
            pkgs = [pkg for pkg, _ in missing_optional]
            cmd = f"pip install {' '.join(pkgs)}"
            print(f"\n  Optional packages:\n    {cmd}")

        print("\n  All packages at once:")
        all_missing = missing_required + [pkg for pkg, _ in missing_optional]
        if all_missing:
            print(f"    pip install {' '.join(all_missing)}")
        else:
            print("    (nothing to install)")

    # =========================================================================
    # Functionality tests
    # =========================================================================
    if args.test:
        print("\n" + "=" * 70)
        print("Functionality Tests")
        print("=" * 70)

        tests_passed = 0
        tests_failed = 0

        # Test 1: Survey import
        print("\n  [TEST 1] Survey class import...")
        try:
            from survey.survey import Survey
            print("    PASS: Survey class imported")
            tests_passed += 1
        except Exception as e:
            print(f"    FAIL: {e}")
            tests_failed += 1

        # Test 2: Models import
        print("\n  [TEST 2] Data models import...")
        try:
            from survey.models import (
                LayerSpec, PackagingSpec, SubstanceSpec,
                PriorSpec, SurveyConfig
            )
            print("    PASS: All models imported")
            tests_passed += 1
        except Exception as e:
            print(f"    FAIL: {e}")
            tests_failed += 1

        # Test 3: Spreadsheet module
        print("\n  [TEST 3] Spreadsheet module...")
        try:
            from survey.spreadsheet import (
                SpreadsheetData, FamilyRow, SubstanceRow,
                check_dependencies as check_spreadsheet_deps
            )
            deps = check_spreadsheet_deps()
            print(f"    PASS: XLSX={deps['xlsx']}, ODS={deps['ods']}")
            tests_passed += 1
        except Exception as e:
            print(f"    FAIL: {e}")
            tests_failed += 1

        # Test 4: Batch runner
        print("\n  [TEST 4] Batch runner...")
        try:
            from survey.batch import (
                BatchRunner, SimulationTask, SimulationResult,
                get_cpu_count, get_default_workers
            )
            cpus = get_cpu_count()
            workers = get_default_workers()
            print(f"    PASS: CPUs={cpus}, DefaultWorkers={workers}")
            tests_passed += 1
        except Exception as e:
            print(f"    FAIL: {e}")
            tests_failed += 1

        # Test 5: FastAPI apps
        print("\n  [TEST 5] FastAPI applications...")
        try:
            from survey.app.main import app as editor_app
            from survey.app.simulator import app as simulator_app
            print(f"    PASS: Editor={editor_app.title}, Simulator={simulator_app.title}")
            tests_passed += 1
        except Exception as e:
            print(f"    FAIL: {e}")
            tests_failed += 1

        # Test 6: Patankar solver
        print("\n  [TEST 6] Patankar solver...")
        try:
            from patankar import dfood
            print("    PASS: Patankar solver available")
            tests_passed += 1
        except Exception as e:
            print(f"    FAIL: {e}")
            tests_failed += 1

        print("\n  " + "-" * 40)
        print(f"  Results: {tests_passed} passed, {tests_failed} failed")

    # =========================================================================
    # Final status
    # =========================================================================
    print("\n" + "=" * 70)

    if missing_required:
        print("STATUS: INCOMPLETE - Install missing required packages")
        return 1
    elif not sfppy_ok:
        print("STATUS: WARNING - Some SFPPy modules not available")
        return 1
    else:
        print("STATUS: READY - All required dependencies installed")
        return 0


if __name__ == "__main__":
    sys.exit(main())
