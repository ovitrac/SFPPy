# Non-Pure Compound Cache for SFPPy

This directory stores metadata for **non-pure substances** — mixtures, polymers, and undefined materials that do NOT have a defined molecular structure in PubChem.

## Purpose

Many substances used in food contact materials are not pure compounds:
- **Mixtures:** Oils, waxes, resins (CAS 8xxx-xx-x series)
- **Polymers:** Polyethylene, polypropylene, etc.
- **Undefined:** Natural products, petroleum derivatives

These substances cannot be looked up in PubChem because they lack a single molecular structure. This cache provides:

1. **Fast failure** — Avoids wasted PubChem API calls for known non-pure substances
2. **Metadata storage** — Records what we know about these substances
3. **Explicit handling** — Users can opt-in to include non-pure substances via `pure_only=False`

## Directory Contents

| File | Description |
|------|-------------|
| `README.md` | This file |
| `nonpure_index.json` | CAS → metadata index |
| `cas{CAS}.json` | Individual substance records (optional) |

## File Format

### nonpure_index.json
Maps CAS numbers to substance metadata:
```json
{
  "8001-79-4": {
    "name": "Castor oil",
    "type": "mixture",
    "category": "vegetable oil",
    "source": "manual",
    "note": "Complex mixture of triglycerides",
    "date_added": "2026-01-09"
  },
  "9002-88-4": {
    "name": "Polyethylene",
    "type": "polymer",
    "category": "polyolefin",
    "source": "manual",
    "note": "No defined MW - polymer",
    "date_added": "2026-01-09"
  }
}
```

### Substance Types

| Type | Description | Examples |
|------|-------------|----------|
| `mixture` | Multi-component blend | Oils, waxes, resins |
| `polymer` | Macromolecule | PE, PP, PET, PS |
| `undefined` | Unknown composition | Natural extracts |
| `uvcb` | Unknown Variable composition | Petroleum fractions |

### Categories

Common categories for food contact materials:
- `vegetable oil` — Plant-derived oils
- `mineral oil` — Petroleum-derived oils
- `wax` — Natural and synthetic waxes
- `resin` — Polymeric resins
- `polyolefin` — PE, PP, etc.
- `polyester` — PET, PBT, etc.
- `other` — Miscellaneous

## Usage

### In migrant class (default: pure_only=True)
```python
from patankar.loadpubchem import migrant

# Default behavior: only pure compounds
m = migrant("8001-79-4")  # Raises error - castor oil is not pure

# Opt-in to include non-pure substances
m = migrant("8001-79-4", pure_only=False)  # Returns metadata from cache.NonPure
```

### Checking if substance is non-pure
```python
from patankar.loadpubchem import is_nonpure_substance

is_nonpure_substance("8001-79-4")  # True (castor oil)
is_nonpure_substance("50-00-0")    # False (formaldehyde - pure)
```

## CAS Number Patterns

Non-pure substances often follow recognizable CAS patterns:

| Pattern | Type | Examples |
|---------|------|----------|
| `8xxx-xx-x` | Natural mixtures | 8001-79-4 (castor oil) |
| `9xxx-xx-x` | Polymers | 9002-88-4 (polyethylene) |
| `6xxxx-xx-x` | Complex mixtures | 64742-54-7 (petroleum distillate) |

## Cache Management

### Add a new non-pure substance
```python
from patankar.loadpubchem import NonPureIndex

db = NonPureIndex()
db.add("8001-79-4", {
    "name": "Castor oil",
    "type": "mixture",
    "category": "vegetable oil",
    "source": "manual"
})
db.save()
```

### Clear this cache
```bash
# Safe to clear - will be rebuilt from known substances
rm patankar/cache.NonPure/nonpure_index.json
```

### View statistics
```bash
python3 -c "import json; d=json.load(open('patankar/cache.NonPure/nonpure_index.json')); print(f'{len(d)} non-pure substances')"
```

## Behavior with pure_only Parameter

| Substance | pure_only=True (default) | pure_only=False |
|-----------|--------------------------|-----------------|
| Pure compound in PubChem | Returns data | Returns data |
| Pure compound NOT cached | Queries PubChem | Queries PubChem |
| Non-pure in cache.NonPure | **Raises error** | Returns metadata |
| Non-pure NOT in any cache | Queries PubChem (fails) | Queries PubChem (fails) |

## Why Separate Caches?

1. **Clarity** — Pure vs non-pure substances have fundamentally different data
2. **Safety** — Default `pure_only=True` prevents accidental use of mixtures for property prediction
3. **Reversibility** — Easy to clear non-pure cache without losing PubChem data
4. **Auditability** — Clear provenance for each type of substance

## Related

- **`cache.PubChem/`** — Cache for pure compounds with molecular structures
- **`cache.ToxTree/`** — Toxicological classification (pure compounds only)
- **`loadpubchem.py`** — Main module managing both caches

---
*Part of SFPPy — Scientific Framework for Food Packaging & Migration Modeling*
