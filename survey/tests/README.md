# Survey Tests

Unit and integration tests for the SFPPy Survey module.

## Running Tests

```bash
# Run all tests
pytest survey/tests/

# Run with verbose output
pytest survey/tests/ -v

# Run specific test file
pytest survey/tests/test_survey.py

# Run with coverage
pytest survey/tests/ --cov=survey
```

## Test Files

| File | Description |
|------|-------------|
| `test_survey.py` | Core Survey class tests |
| `test_models.py` | Data model tests |
| `test_batch.py` | Batch runner tests |
| `test_api.py` | API endpoint tests |

## Test Fixtures

Test fixtures are located in `../examples/` and include:
- YAML scenario files
- XLSX spreadsheet samples
- JSON PF job exports

## Writing Tests

```python
import pytest
from survey import Survey
from survey.models import PriorSpec

def test_survey_from_scenario():
    survey = Survey.from_scenario("examples/batch_scenarios/water_PET_bottle.yml")
    assert survey is not None
    assert len(survey.substances) > 0

def test_prior_spec():
    prior = PriorSpec(mode=100, max_val=500)
    assert prior.min_val == 0  # Always 0 by design
    assert prior.mode == 100
```

---
*Part of SFPPy Survey Module*
