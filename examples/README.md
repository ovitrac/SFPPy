# SFPPy Examples  🍏⏩🍎

Progressive tutorials demonstrating food packaging migration simulation — from basic monolayer analysis to advanced parameter fitting and AI-assisted regulatory queries.

## Quick Start

```bash
cd examples/
python example1.py
```

All examples include a **path bootstrap** and can be run directly without installing SFPPy.

## Learning Path

| Example | Level | Topic | Key Concepts |
|---------|-------|-------|--------------|
| [example1.py](#example-1-monolayer-migration) | Beginner | Monolayer migration | Geometry, layers, food, solver |
| [example1_extensions.py](#example-1-extensions-tiered-assessment) | Beginner | Tiered M2/M3 assessment | Equilibrium, kinetics, chaining |
| [example2.py](#example-2-functional-barriers) | Intermediate | Functional barriers | Multilayer `+`, parametric study |
| [example3.py](#example-3-chained-simulations) | Intermediate | Multi-step scenarios | Operators `>>`, `@`, `%`, variants |
| [example3_shortvariant.py](#example-3-short-potential-release) | Intermediate | Potential release (PR) | Risk decomposition by step/layer |
| [example4.py](#example-4-parameter-fitting) | Advanced | Parameter estimation | `layerLink`, fitting, sensitivity |
| [example5.py](#example-5-rag-knowledge-base) | Advanced | Regulatory AI queries | LlamaIndex, Ollama, RAG |

---

## Example 1: Monolayer Migration

**File:** `example1.py`

**Scenario:** Migration of antioxidants (Irganox 1076, Irgafos 168) from a 100 µm LDPE film into a fatty sandwich over 10 days at 7°C.

### What You Learn

- Define packaging geometry with `Packaging3D`
- Load migrant properties from PubChem with `migrant()`
- Create polymer layers (`LDPE`, `PP`, `gPET`, etc.)
- Define food contact conditions via multiple inheritance
- Run migration simulations with `senspatankar`
- Compare multiple migrants with `CFSimulationContainer`
- Plot and export results

### Key Code Patterns

```python
# Geometry
sandwich_geom = Packaging3D('Cylinder', length=(19, "cm"), radius=(30, "mm"))
V, A = sandwich_geom.get_volume_and_area()

# Migrant from PubChem
m = migrant("irganox 1076")

# Polymer layer with initial concentration
film = polymer.LDPE(l=(100, "um"), substance=m, C0=5000, T=(7, "degC"))

# Food definition (multiple inheritance)
class sandwich(food.realfood, food.semisolid, food.fat):
    name = "sandwich"

F = sandwich(volume=V, surfacearea=A, contacttime=(10, "days"), ...)

# Simulation
result = solver(film, F, name="LDPE-sandwich")
result.plotCF()  # migration kinetics
result.plotCx()  # concentration profiles
```

### Compact Syntax (Operators)

```python
# substance % food << geometry >> layer >> food
m % F << sandwich_geom >> film >> F
F.lastsimulation.plotCF()
```

---

## Example 1 Extensions: Tiered Assessment

**File:** `example1_extensions.py`

**Scenario:** Same as Example 1, but demonstrates the tiered approach used in regulatory assessments.

### What You Learn

- **Tier M2 (Equilibrium):** Compute partition coefficient $K_{F/P}$
- **Tier M3 (Kinetics):** Full diffusion simulation
- Chain simulations with `>>` for multi-step contact
- Calculate severity relative to SML

### Key Concepts

```python
# M2: Equilibrium partition coefficient
K_F_over_P = film.k / F.k0  # <1 means stronger polymer affinity

# M3: Kinetics simulation
kin = F.migration(film)
severity = 100 * kin.CFtarget / m.SML  # % of SML

# Chaining: add warming step (4h @ 25°C)
F2 = F.copy().update(contacttime=(4, "hours"), contacttemperature=(25, "C"))
kin2 = kin >> F2  # continues from previous state
(kin + kin2).plotCF()  # combined kinetics
```

---

## Example 2: Functional Barriers

**File:** `example2.py`

**Scenario:** Migration of toluene from a 300 µm recycled PP bottle into fatty liquid food. Evaluate the protective effect of a PET functional barrier (2–60 µm).

### What You Learn

- Create multilayer structures with `+` operator
- Define bottle geometry (body + neck)
- Run parametric studies (varying barrier thickness)
- Compare scenarios with color-coded plots

### Key Code Patterns

```python
# Multilayer: PET barrier + PP core
PP_core = polymer.PP(l=(300, "um"), substance=toluene, C0=10, T=(20, "degC"))
PET_barrier = polymer.wPET(l=(30, "um"), substance=toluene, C0=0, T=(20, "degC"))

# Assembly (food contacts left side)
bilayer = PET_barrier + PP_core

# Parametric study
for thickness in range(2, 61, 4):
    bilayer.l[0] = _toSI((thickness, "um")).item()
    sim = solver(bilayer, food, name=f"FB-{thickness}um")
    results.add(sim, f"FB = {thickness} µm", color)

results.plotCF()  # all curves overlaid
```

---

## Example 3: Chained Simulations

**File:** `example3.py`

**Scenario:** Limonene migration from a trilayer ABA container (PET/PP/PET) through three processing steps:
1. Setoff storage (4 months @ 20°C)
2. Hot-filling (90°C)
3. Long-term storage (6 months @ 30°C)

### What You Learn

- Build trilayer structures: `ABA = A + B + A`
- Chain multiple contact steps with `>>`
- Compare variants (different migrants, thicknesses)
- Export results to Excel/CSV

### Operator Reference

| Operator | Meaning | Example |
|----------|---------|---------|
| `+` | Assemble layers / merge results | `ABA = A + B + A` |
| `>>` | Propagate & simulate | `food >> layer >> food` |
| `@` | Thermalize (no simulation) | `food @ layer` |
| `%` | Inject substance | `m % food >> layer >> food` |

### Key Code Patterns

```python
# Trilayer assembly
A = gPET(l=(20, "um"), migrant=m, C0=0)
B = PP(l=(0.5, "mm"), migrant=m, C0=200)
ABA = A + B + A

# Define contact steps
class Step1(stacked, ambient): contacttime = (4, "months")
class Step2(hotfilled, realfood, liquid, fat): pass
class Step3(ambient, realfood, liquid, fat): contacttime = (6, "months")

# Propagate geometry
container >> step1 >> step2 >> step3

# Chain simulation through all steps
step1 >> ABA >> step1 >> step2 >> step3

# Merge results
full_kinetics = step1.lastsimulation + step2.lastsimulation + step3.lastsimulation
full_kinetics.plotCF()

# Variant: different migrant
m2 = migrant("toluene")
step1 @ ABA.update(solute=m2) >> step1 >> step2 >> step3
```

---

## Example 3 Short: Potential Release

**File:** `example3_shortvariant.py`

**Scenario:** Compact version of Example 3 with **Potential Release (PR)** analysis — decomposing migration risk by step and by layer.

### What You Learn

- Compute intrinsic PR per processing step
- Compute layer-based PR in multilayer systems
- Quantify the contribution of each step to total migration

### Key Concepts

```python
# Potential Release analysis
R123 = F1.potentialRelease(ABA) >> F2 >> F3
R23 = F2.potentialRelease(ABA) >> F3  # step 1 omitted
R13 = F1.potentialRelease(ABA) >> F3  # step 2 omitted

# Intrinsic PR per step
PR[0] = 1 - (1 - PRN) / (1 - PR23)  # contribution of step 1
PR[1] = 1 - (1 - PRN) / (1 - PR13)  # contribution of step 2
PR[2] = 1 - (1 - PRN) / (1 - PR12)  # contribution of step 3
```

---

## Example 4: Parameter Fitting

**File:** `example4.py`

**Scenario:** Estimate diffusivity (D) and partition coefficient (k) from experimental kinetic data using pseudo-experiments and optimization.

### What You Learn

- Use `layerLink` for dynamic parameter modification
- Generate pseudo-experimental data with noise
- Perform sensitivity analysis
- Fit parameters to recover original values
- Detect overfitting

### Key Code Patterns

```python
# Create layer with linkable parameters
P = layer(l=(100, "um"), D=(1e-10, "cm**2/s"), k=0.1, C0=1000)

# Create layerLink objects
D_link = layerLink("D", indices=0, values=P.D)
k_link = layerLink("k", indices=0, values=P.k)
P.Dlink, P.klink = D_link, k_link

# Simulation
R = F.migration(P)

# Generate pseudo-experiment (30 points, 1% noise)
E = R.pseudoexperiment(npoints=30, std_relative=0.01)

# Squared error function
d2 = R - E  # callable: d2() returns current error

# Sensitivity: vary D and k
for i in range(10):
    D_link[0] /= 1.1
    k_link[0] *= 1.1
    R.rerun(name=f"iter {i}")
    print(f"Error change: {100 * d2() / d2_orig - 100:.2f}%")

# Fit: optimize D and k
R.fit(E)
print(f"Fitted D = {D_link.values}, Original = {D_orig}")
```

---

## Example 5: RAG Knowledge Base

**File:** `example5.py`

**Scenario:** Query a local Markdown knowledge base (EU Regulation 10/2011) using Retrieval-Augmented Generation (RAG) with a local LLM.

### What You Learn

- Index Markdown documents with LlamaIndex
- Use local embeddings (HuggingFace BGE)
- Query with local LLM (Mistral via Ollama)
- Traceable, auditable responses with source citations

### Prerequisites

```bash
# Install RAG dependencies
pip install llama-index llama-index-embeddings-huggingface \
            llama-index-llms-ollama sentence-transformers chromadb

# Install Ollama and pull Mistral
# https://ollama.ai
ollama pull mistral
```

### Usage

```bash
# Interactive mode
python example5.py

# With specific query
python example5.py "What is a functional barrier?"

# Rebuild index
python example5.py --rebuild
```

### Key Architecture

```
docs/KB/          # Markdown knowledge base
docs/KBIndex/     # Persistent vector index (ChromaDB)
```

```python
# Minimal RAG pipeline
documents = SimpleDirectoryReader("docs/KB").load_data()
embed_model = HuggingFaceEmbedding(model_name="BAAI/bge-small-en-v1.5")
index = VectorStoreIndex.from_documents(documents, embed_model=embed_model)
llm = Ollama(model="mistral")
response = index.as_query_engine(llm=llm).query("What is migration?")
```

---

## Output Files

Examples generate outputs in `examples/tmp/`:

| File Type | Content |
|-----------|---------|
| `*.png` | High-resolution plots |
| `*.pdf` | Vector graphics |
| `*.xlsx` | Excel spreadsheets |
| `*.csv` | Raw data tables |

---

## Troubleshooting

### Import Errors

Ensure you're running from the `examples/` directory or project root:

```bash
cd /path/to/SFPPy
python examples/example1.py
```

### Missing Dependencies

```bash
# Core SFPPy
mamba install numpy scipy matplotlib pandas

# Example 5 (RAG)
pip install llama-index llama-index-embeddings-huggingface chromadb
```

### Slow First Run

First execution downloads:
- PubChem data (cached locally)
- Embedding models (~100 MB for Example 5)

Subsequent runs use cached data.

---

## Further Reading

- **Full documentation:** `docs/`
- **API reference:** `patankar/` module docstrings
- **Notebooks:** `notebooks/` (Jupyter tutorials)
- **Survey framework:** `survey/` (batch simulations)
- **Studio GUI:** `studio/` (interactive web interface)
