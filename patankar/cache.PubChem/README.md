# PubChem Cache for SFPPy 🍏⏩🍎

This directory is managed by `patankar.private.loadpubchem` and serves as a **local cache** for chemical data retrieved from PubChem. It helps optimize repeated queries by storing chemical properties locally, allowing for faster access and offline usage.

## 📁 Folder Contents
- `pubchem_index.json` 📖: A local index mapping compound names, CAS numbers, and synonyms to their **CID** (PubChem Compound ID).

- `cidXXXX.full.json` 📂: Full records containing **all available PubChem properties** for a compound.

- `cidXXXX.simple.json` 📑: Lightweight records with only the most essential fields, including:

  - CID
  - Name & Synonyms
  - CAS number
  - Molecular Formula
  - Molecular Weight
  - SMILES / InChI / InChIKey
  - LogP value (if available)
  - Additional metadata or temporary files created during data retrieval.

- `structure/cidXXXX.sdf` 🧱: Cached **SDF structure file** (Toxtree needs them).

- `thumbs/cidXXXX.png` 🎨: **Cropped PNG image** of the molecule structure.

  

  > These properties are easily accessible from migrant that compute additional quantities from these inputs.

  ```python
  from loadpubchem import migrant
  m = migrant("Irganox 1010")
  print(m.cid)
  print(m.compound)
  # etc.
  m.image # shows the molecule in IPython
  m.structure["atoms"] # show atom coordinates (DataFrame)
  m.structure["bonds"] # show bonds (DataFrame)
  m.structure["metadata"].T # show metadata (DataFrame)
  ```

## 🔹 How It Works
1. When a compound is requested, `loadpubchem` first checks `pubchem_index.json` to see if it exists locally.
2. If found, it retrieves data from `cidXXXX.simple.json` or `cidXXXX.full.json`, depending on the requested format.
3. If not found, a query is made to **PubChem**, and the results are stored locally in both formats for future use.
4. The index is automatically updated with new synonyms and compound mappings to improve efficiency.

## ⚠️ Notes
- **Do not manually modify or delete files** unless necessary, as the `loadpubchem` module automatically manages the cache.
- If the index becomes outdated or corrupted, it will be rebuilt when `loadpubchem` detects inconsistencies.
- To force an update, you can clear this folder and rerun queries.

🔬 **Using local caching improves performance and allows for seamless offline chemical property retrieval in SFPPy.**