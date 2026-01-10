# PubChem Cache for SFPPy

This directory is managed by `patankar.loadpubchem` and serves as a **local cache** for **pure chemical compounds** retrieved from PubChem. It enables offline operation and reduces API calls to PubChem.

## Directory Contents

| File/Folder | Description |
|-------------|-------------|
| `pubchem_index.json` | Synonym → CID index (126k+ entries) |
| `.pubchem_rate.lock` | Inter-process rate limiter lock file |
| `cid{N}.full.json` | Complete PubChem record for CID N |
| `cid{N}.simple.json` | Lightweight record (name, CAS, M, SMILES, etc.) |
| `structure/{CID}.sdf` | 3D molecular structure files (for ToxTree) |
| `thumbs/{CID}.png` | 2D structure images (large, variable size) |
| `thumbs_tiny/{CID}.png` | 2D structure images (small, 100x100 for UI) |
| `sync_thumbs_tiny.py` | Script to synchronize tiny thumbnails |

## File Formats

### pubchem_index.json
Maps all known synonyms (names, CAS numbers, IUPAC names) to their PubChem CIDs:
```json
{
  "formaldehyde": [712],
  "50-00-0": [712],
  "methanal": [712]
}
```

### cid{N}.simple.json
Lightweight record with essential properties:
```json
{
  "CID": 712,
  "name": "formaldehyde",
  "synonyms": ["methanal", "formalin"],
  "CAS": ["50-00-0"],
  "M": 30.026,
  "formula": "CH2O",
  "SMILES": "C=O",
  "InChI": "InChI=1S/CH2O/c1-2/h1H2",
  "InChIKey": "WSFSSNUMVMOOMR-UHFFFAOYSA-N",
  "logP": 0.35,
  "date": "2026-01-09"
}
```

These properties are accessible from the `migrant` class:

```python
from loadpubchem import migrant
m = migrant("Irganox 1010")
print(m.cid)
print(m.compound)
m.image        # shows the molecule in IPython
m.structure["atoms"]    # atom coordinates (DataFrame)
m.structure["bonds"]    # bond information (DataFrame)
```

## How It Works

1. When a compound is requested, `loadpubchem` checks `pubchem_index.json` for local matches
2. If found, data is loaded from `cid{N}.simple.json` or `cid{N}.full.json`
3. If not found, a query is made to **PubChem** with rate limiting and exponential backoff
4. Results are stored locally and the index is updated with new synonyms

## Rate Limiting (v1.50+)

The `.pubchem_rate.lock` file coordinates rate limiting across parallel processes:

- **Target rate:** 3 requests/second (PubChem NCBI policy)
- **Lock expiry:** 5 seconds (handles crashed processes)
- **Retry strategy:** Exponential backoff with max 5 retries
- **Implementation:** Atomic writes via temp file + rename

## What Belongs Here vs cache.NonPure

| Type | Example | This Cache? |
|------|---------|-------------|
| Pure compounds | Formaldehyde (50-00-0) | **Yes** |
| Monomers | Styrene (100-42-5) | **Yes** |
| Additives | BHT (128-37-0) | **Yes** |
| Mixtures | Castor oil (8001-79-4) | No → `cache.NonPure/` |
| Polymers | Polyethylene | No → `cache.NonPure/` |
| Undefined | Petroleum wax (8002-74-2) | No → `cache.NonPure/` |

**Rule:** If it has a single, defined molecular structure (SMILES, InChI) in PubChem, it belongs here.

## Cache Management

### Refresh entire index
```python
from patankar.loadpubchem import CompoundIndex
db = CompoundIndex()
db.refresh_index()  # Rebuilds from all .full.json files
```

### Clear cache (use with caution)
```bash
rm -rf patankar/cache.PubChem/cid*.json
rm -rf patankar/cache.PubChem/structure/
rm -rf patankar/cache.PubChem/thumbs/
rm patankar/cache.PubChem/pubchem_index.json
```

### View statistics
```bash
# Count cached compounds
ls patankar/cache.PubChem/*.simple.json | wc -l

# Check index size
python3 -c "import json; d=json.load(open('patankar/cache.PubChem/pubchem_index.json')); print(f'{len(d)} entries')"
```

## Thumbnail Synchronization

The `sync_thumbs_tiny.py` script downloads missing tiny thumbnails (100x100) from PubChem to enable **offline operation of SFPPy Studio**.

### Usage
```bash
cd patankar/cache.PubChem

# Dry run - see what would be downloaded
python sync_thumbs_tiny.py --dry-run

# Download all missing thumbnails (~6 req/sec)
python sync_thumbs_tiny.py

# Limit downloads (useful for testing)
python sync_thumbs_tiny.py --limit 100

# Custom delay between requests
python sync_thumbs_tiny.py --delay 0.2

# Quiet mode (no per-file output)
python sync_thumbs_tiny.py --quiet

# Reset progress and re-attempt failed downloads
python sync_thumbs_tiny.py --reset
```

### Options
| Option | Description |
|--------|-------------|
| `--dry-run, -n` | Show what would be downloaded without downloading |
| `--limit N, -l N` | Limit to N downloads (0 = no limit) |
| `--delay SECONDS, -d` | Delay between requests (default: 0.17s = ~6/sec) |
| `--reset` | Clear progress file and retry all missing |
| `--quiet, -q` | Suppress per-file output |

### Features
- **Resumable**: Progress saved in `.sync_progress.txt`
- **Reliable**: Exponential backoff on 503/429 errors, 3 retries
- **Validated**: Verifies PNG magic bytes before saving
- **Progress**: Reports status every 50 downloads

### PubChem API Used
```
https://pubchem.ncbi.nlm.nih.gov/image/imgsrv.fcgi?cid={CID}&t=s
```
Where `t=s` = small (100x100), `t=l` = large thumbnail.

## Notes

- **Do not manually modify files** — the cache is auto-managed
- If the index becomes corrupted, it will be rebuilt automatically
- To force refresh, delete the index file and restart

## Related

- **`cache.NonPure/`** — Cache for mixtures, polymers, and undefined substances
- **`cache.ToxTree/`** — Toxicological classification results (Cramer classes)
- **`loadpubchem.py`** — Main module managing this cache

---
*Part of SFPPy — Scientific Framework for Food Packaging & Migration Modeling*
