# ToxTree Cache for SFPPy 🍏⏩🍎

This directory is managed by `patankar.private.loadpubchem.migrantToxtree` and serves as a **local cache** for toxicological assessments performed using **Toxtree**. It helps optimize repeated evaluations by storing results locally, allowing for faster access usage (even on machines without ToxTree).

## 📁 Folder Contents

- `cidXXXX.default.csv` 📄: **Base Toxtree evaluation** (always generated for each compound).
- `cidXXXX.engine.csv` 📑: Engine-specific Toxtree outputs, where `engine` corresponds to:
  - `cramer`
  - `cramer2`
  - `cramer3`
  - `kroes`
  - `dnabinding`
  - `skin`
  - `eye`
  - `Ames`
- `cidXXXX.engine.json` 📂: Processed **JSON versions** of the Toxtree outputs for faster lookup.

## 🔹 How It Works

1. When a compound is assessed, `migrantToxtree` first checks if **JSON results exist**.
2. If found, data is retrieved directly from JSON (fastest lookup).
3. If missing, `Toxtree` is executed to generate a **CSV file**.
4. The CSV results are then **converted into JSON** and stored for future use.
5. If structure files are not available, they are downloaded from **PubChem** and stored in `structure/`.

## ⚠️ Notes

- **Do not manually modify or delete files** unless necessary, as `migrantToxtree` automatically manages the cache.
- To force an update, use the `refresh=True` flag when initializing `migrantToxtree`.
- To fully regenerate all files, use `no_cache=True` to overwrite existing data.

🚀 **Using local caching improves performance and allows for seamless offline toxicological assessments in SFPPy.**