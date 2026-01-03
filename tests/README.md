# SFPPy Test Suite

Unit and integration tests for the SFPPy framework.

## Quick Start

```bash
# Run all tests
python -m pytest tests/ -v

# Run specific test file
python -m pytest tests/test_units.py -v

# Run specific test class
python -m pytest tests/test_layer.py::TestMultilayer -v

# Run with coverage report
python -m pytest tests/ --cov=patankar --cov-report=term-missing
```

## Test Structure

```
tests/
├── __init__.py          # Suite description
├── conftest.py          # Shared fixtures
├── test_units.py        # Unit conversion tests
├── test_layer.py        # Layer property tests
└── README.md            # This file
```

## Test Tiers

### Tier 1: Unit Tests (Current)

Critical function validation with fast execution.

| File | Tests | Scope |
|------|-------|-------|
| `test_units.py` | 30 | `_toSI()`, `check_units()`, unit consistency |
| `test_layer.py` | 22 | Layer creation, properties, multilayer, mesh |

### Tier 2: Integration Tests (Planned)

- Migration simulations with reference solutions
- Multilayer diffusion scenarios
- Food-packaging contact workflows

### Tier 3: Regression Tests (Planned)

- Full example validation against stored results
- Numerical stability checks

## Fixtures (conftest.py)

Reusable test objects defined in `conftest.py`:

```python
# Layers
ldpe_layer      # LDPE 100 µm, C0=1000, T=25°C
hdpe_layer      # HDPE 50 µm, C0=500, T=40°C
gpet_layer      # gPET 12 µm, C0=0, T=25°C

# Food simulants
ethanol_food    # Ethanol, 10 days, 25°C
water_food      # Water, 10 days, 25°C

# Migrants
limonene_migrant
toluene_migrant

# Tolerances
rtol            # 1e-6 (relative)
atol            # 1e-12 (absolute)
```

Usage in tests:

```python
def test_example(ldpe_layer, limonene_migrant):
    """Test using fixtures."""
    assert ldpe_layer.D > 0
    assert limonene_migrant.M > 0
```

## Writing New Tests

### 1. Create test file

Name: `test_<module>.py`

```python
"""
Tests for <module> functionality.

@project: SFPPy - Safe Food Packaging in Python
@author: Olivier Vitrac
@license: MIT
"""

import pytest
import numpy as np
from numpy.testing import assert_allclose


class TestFeatureName:
    """Tests for feature X."""

    def test_basic_case(self):
        """Description of what is tested."""
        from patankar.module import function
        result = function(input)
        assert_allclose(result, expected, rtol=1e-6)
```

### 2. Use fixtures for common objects

Add to `conftest.py`:

```python
@pytest.fixture
def my_fixture():
    """Create object for testing."""
    from patankar.module import MyClass
    return MyClass(param=value)
```

### 3. Test categories

- **Creation tests**: Object instantiation
- **Property tests**: Computed values
- **Operator tests**: `+`, `>>`, `%`, `@`
- **Error tests**: Invalid inputs raise exceptions

## Running Tests

### Basic

```bash
python -m pytest tests/ -v
```

### With output capture disabled (see print statements)

```bash
python -m pytest tests/ -v -s
```

### Stop on first failure

```bash
python -m pytest tests/ -v -x
```

### Run only failed tests from last run

```bash
python -m pytest tests/ -v --lf
```

### Parallel execution (requires pytest-xdist)

```bash
python -m pytest tests/ -v -n auto
```

## Dependencies

Required:
- `pytest`

Optional:
- `pytest-cov` (coverage reports)
- `pytest-xdist` (parallel execution)

Install with:

```bash
mamba install pytest pytest-cov -y
```

## Expected Output

```
tests/test_layer.py::TestLayerCreation::test_ldpe_creation PASSED
tests/test_layer.py::TestLayerCreation::test_hdpe_creation PASSED
...
tests/test_units.py::TestToSI::test_micrometers_to_meters PASSED
...
============================== 52 passed in 0.67s ==============================
```

## Troubleshooting

### Import errors

Ensure you're running from the project root:

```bash
cd /path/to/python
python -m pytest tests/ -v
```

### Missing fixtures

Check that `conftest.py` is present in `tests/` directory.

### Slow tests

Use `-x` to stop on first failure during development:

```bash
python -m pytest tests/ -v -x --tb=short
```
