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
| [example6_cosmetic.py](#example-6-cosmetic-exposure) | Advanced | Cosmetic exposure & verdicts | `cosmetic`, dose, tiers, gates |
| [example7_cosmetic_populations.py](#example-7-populations-and-comparable-units) | Advanced | Populations & unit conversions | body weight, retention, SML→TDI |
| [example8_cosmetic_specification.py](#example-8-from-exposure-to-a-material-specification) | Advanced | **Reverse problem** — MAE → CP0max | `PR`, formats, `setCF`/`setMAE` |

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

## Example 6: Cosmetic Exposure

**File:** `example6_cosmetic.py`

**Scenario:** Substances migrating from a **recycled PET bottle** into a **body lotion** — assessed all the way to a dose and a verdict.

### Why this example is different

Every other example stops at $C_F$, the concentration in the contacting phase. For food that is enough: EU practice folds the whole exposure chain into the limit itself (6 dm²/kg food, 1 kg/day, 60 kg body weight), so the SML *is* the exposure model.

Regulation (EC) 1223/2009 carries no such convention for cosmetic packaging. The chain has to be modelled explicitly, per product category, and it ends in a **dose**, not a concentration:

```
C_product --(amount applied, retention)--> dose --(area, absorption)-->
    systemic [mg/kg bw/day] -> DNEL / TTC
    local    [µg/cm²]       -> site-of-contact threshold
```

### What You Learn

- Assess the **unidentified population** with `product.declare_unidentified()` — one T0 bound, no molecule named
- Assess a named substance: `substance % product`, then `product.assess()`
- Read a `cosmeticassessment` — dose, **binding threshold with its tier and origin**, margin, and the provenance of every input
- Run the standard pipeline unchanged: `substance % product << geometry >> packaging >> product`
- Climb tiers **only on failure**, and watch `refine()` refuse when nothing failed (gate G8)
- Invert the question: `maxconcentration()` and `required_decontamination()` turn a verdict into a supplier specification
- See **gates** refuse rather than return a plausible wrong number

### The idea it is built around

> We do not need to simulate *accurately*. We need to simulate **robustly**.
> Refinement is required only when a problem is found — and then the exit is a
> **branch, not an escalation**: refine one rung, **or** measure.

The example makes this concrete. In step 4 the crude one-line **M0 total-transfer bound already passes**, so the Patankar solver is never needed; it is run only to show that M0 bounds M3. Running it anyway would have replaced a parameter-free bound with an estimate needing $D$, $k$ and contact conditions — three more things to defend, for a conclusion already reached.

### The unidentified population

The commonest substance in a recyclate is the one nobody named. `declare_unidentified()` says so explicitly, and tier T0 answers it — the genotoxic TTC applies to *any* molecule, so **one bound covers them all** without a structure, a CAS or a database lookup.

The declaration is deliberate: leaving `substance` as `None` remains a refusal (G2), so a forgotten injection never becomes a silent T0 answer. And `unidentified % product` raises, so the sentinel can never reach the solver — tier M3 is structurally unreachable for a molecule with no structure, which is correct.

The example shows the same T0 bound moving **27×** between a leave-on body lotion and a rinse-off washing gel for a child. Same toxicology, same threshold; retention and body weight did that. A single food-style limit cannot express it.

### Running

```bash
python example6_cosmetic.py
```

No arguments. First run populates the PubChem/ToxTree caches for toluene and bisphenol A.

### Product categories

`bodylotion` (leave-on, adult, EtOH 95), `shampoo` (rinse-off, adult, EtOH 50), `washinggel` (rinse-off, **child**, EtOH 50). Each inherits its prescribed simulant, so the simulant doctrine is carried by the type rather than passed as an argument.

### Gates

A gate is the engine declining to return a number it cannot defend. `help_gates()` lists all twelve. The example fires five, including:

- **G1** — an unrecognised unit is refused, because source data spans a 1000× mg/µg divide
- **G7** — systemic and local are **not commensurable**; absence of a local threshold is an explicit state, never a silent default
- **G9** — in SFPPy the `r-` prefix means **RUBBERY**, not recycled

### Attribution

The exposure classes, simulant doctrine, three product categories and Tier 1 defaults come from the **CosPaTox consortium**. Any use must cite the April 2024 deliverables (Guideline, Substance List, Scientific Dossier, WCC Calculator) at <https://cospatox.com/publication/>. No DOI or consortium-mandated citation format exists — cite by name, consortium, month/year and URL.

---

## Example 7: Populations and Comparable Units

**File:** `example7_cosmetic_populations.py`

**Scenario:** The same substance, at the same concentration, assessed across **four populations × two usage types** — and reduced to one comparable unit.

### The unit that makes the comparison possible

Everything ends in **mg of substance per kg of body weight per day**. That is the unit of a TTC and of a TDI, so exposure and threshold can simply be divided. Getting every input into it is most of the work, and it is where the mistakes live.

**Three conversions:**

1. **Exposure → dose**

   $$E = C_F \times \underbrace{(\text{amount} \times \text{frequency} \times \text{retention})}_{\text{kg product/day}} \times \frac{\alpha}{\text{bodyweight}}$$

   `retention` separates leave-on (≈1) from rinse-off (≈0.01). It is *not* a dilution of the migrant in the formulation — it is the fraction of applied product that stays on the skin.

2. **SML → TDI** — an SML limits a concentration *in food*, derived from a dose through a convention. Running it backwards recovers the dose:

   $$TDI = SML \times \frac{1\ \mathrm{kg\ food/day}}{60\ \mathrm{kg\ adult}}$$

   `patankar.cosmetic` applies exactly this, and every threshold built that way **carries the convention in its note**.

3. **TTC → nothing.** Already mg/kg bw/day. No convention, no assumption about who is exposed.

### The asymmetry worth understanding

The 60 kg in conversion 2 is **not** an assumption about who uses the cosmetic. It only *undoes* the food convention buried inside the SML. What comes out is a dose per kg of body weight, which is intrinsically population-neutral — so applying it to a 5 kg infant needs no adjustment.

The population enters on the **exposure** side, through body weight, applied amount and exposed area. Keeping the two sides separate is what makes the comparison honest.

### What the grid shows

At $C_F = 5$ mg/kg of bisphenol A, against a T2 threshold of $8.33\times10^{-4}$ mg/kg bw/day:

```
usage       population   mg/kg bw/day    margin   verdict
leave-on    adult           6.667e-04      1.25     PASS
leave-on    child           1.000e-03     0.833     EXCEED
leave-on    toddler         1.250e-03     0.667     EXCEED
leave-on    infant          1.500e-03     0.556     EXCEED
rinse-off   adult           6.667e-06     125       PASS
rinse-off   infant          1.500e-05      55.6     PASS
```

**Identical chemistry in every row.** Same substance, same concentration, same threshold. The verdict is not — and where it flips is decided by who uses the product and how, **never by the transfer model**. A leave-on product on a small body weight is the demanding case, and no refinement of the migration model addresses it.

A safety assessment that reports one number for "the product" has silently chosen a population. Here the choice is explicit and auditable.

### Design answers, per population

Step 4 inverts it into maximum tolerable concentrations, including a column that needs **no substance at all** — the genotoxic TTC applied to the unidentified population via `declare_unidentified()`. For a leave-on product on an infant that is the tightest specification in the table, and the one a recyclate must meet before any identification work begins.

### ⚠️ The exposure figures are placeholders

The populations, applied amounts and exposed areas are **illustrative**, set explicitly per instance precisely to show they are inputs to be sourced — not constants of nature. The module's own class defaults are likewise declared as placeholders. Replace them with the **SCCS Notes of Guidance** and the published **CosPaTox** values before any external use.

---

## Example 8: From Exposure to a Material Specification

**File:** `example8_cosmetic_specification.py`

**Scenario:** The **reverse problem** — given what a consumer may be exposed to, what is the maximum concentration the incoming material may contain?

Examples 6 and 7 answer *"is this safe?"*. This one answers what a formulator and a recycler actually negotiate over, and its output is a **specification a decontamination process can be held to**.

### The inversion, and why it is exact

Transport is **linear in the initial concentration** — double the load in the wall and every concentration in the product doubles. So the potential release

$$PR = \frac{C_F}{C_{P0}}$$

does not depend on $C_{P0}$. One run at any convenient guess characterises the system:

$$MAE \;\rightarrow\; C_{F,\max} \qquad C_{P0,\text{guess}} \;\rightarrow\; C_{F,\text{guess}} \qquad \boxed{C_{P0,\max} = C_{P0,\text{guess}}\frac{C_{F,\max}}{C_{F,\text{guess}}} = \frac{C_{F,\max}}{PR}}$$

The example **demonstrates** the guess cancelling rather than asserting it, sweeping it over seven orders of magnitude for one invariant answer.

> `PR` here is on a **concentration basis**. SFPPy's solver also carries the mass-basis family in `SensPatankarResult.PR` ($PRE = m_{Feq}/m_0$, $PR(C_F) = C_F V_F/m_0$, $PRT = PR/PRE$). Both are linear in $C_{P0}$, so the inversion is identical — but the two must not be read against each other.

### API

```python
PR, cp0, cf, prov = product.potentialrelease(wall, tier="M0")   # or "M3"
r = product.maxinitialconcentration(wall, MAE=None, tier="M0", verify=True)
r.CP0max, r.CFmax, r.PR, r.CP0guess, r.threshold, r.provenance
```

`tier="M0"` is the analytic mass balance; `tier="M3"` uses a stored Patankar solution from the standard pipeline. **M0 yields the tighter specification**, so the cheap route is the conservative one.

### Feeding a decision without computing it

Most decisions come from a detection limit, an analytical result or a limit somebody set — not from a simulation.

```python
product.setCF(0.02, unit="mg/kg", origin="censored")        # a detection limit
product.setMAE(1e-4, unit="mg/kg bw/day", origin="internal limit", tier="T2")
```

Suggested origins: `measured`, `censored`, `MAE-derived`, `supplier`, `regulatory`, `assumed`. An explicitly set **MAE outranks the whole derived threshold ladder** — an internal specification, a customer requirement or a dossier DNEL is as binding as anything derivable, and often tighter. The origin travels with the result.

### Cosmetic formats

`SHAPE_REGISTRY` covers the common primary packs: `tube`, `tottle`, `airless`, `pump`, `aerosol`, `bidon`, `stick`, `roll_on`, `vial`, `ampoule`, `sachet`, `stickpack`, `pouch`, `compact`, `palette`, `cup`, `tray`.

```
format        A/V (1/m)
bidon              83.3
tube / airless    100.0
stick             191.7
compact           266.7
sachet            565.0     <- 6.8x the bidon
```

Migration does not care about the silhouette, only about how much wall faces how much product. **A single-dose sachet is the demanding format in almost any portfolio**, and that falls out of the geometry with no chemistry involved. Equally, `tube` and `airless` dimensioned to the same A/V carry the *same* specification despite looking nothing alike.

> ⚠️ **Varnished metal** (`aerosol`, `bidon`). The silhouette is a cylinder; the regime is not. A 5–10 µm coating on impermeable metal holds the **entire** reservoir in the coating — a finite-dose problem tending to near-total depletion, governed by $K$ rather than $D$. **Model the coating as the layer, never the metal.**

### Robustness

`verify=True` re-runs the transfer at ten times the guess and checks $PR$ is unchanged. Linearity is the assumption the whole inversion rests on, so it is **measured, not trusted**: a concentration-dependent $D$ or $k$ fires gate **G13-linearity** instead of returning a number that looks reasonable and is not. $PR = 0$ reports *unbounded* explicitly rather than dividing by zero; a missing $C_{P0}$, an absent M3 simulation and an unknown tier each refuse with a named remedy.

---

## Further Reading

- **Full documentation:** `docs/`
- **API reference:** `patankar/` module docstrings
- **Notebooks:** `notebooks/` (Jupyter tutorials)
- **Survey framework:** `survey/` (batch simulations)
- **Studio GUI:** `studio/` (interactive web interface)
