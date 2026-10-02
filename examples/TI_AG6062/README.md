# Case studies of Techniques de l'Ingénieur AG 6 062

Scripts of the seven case studies (C1–C7) of:

> Vitrac O., Nguyen P.-M. *Conception sûre des emballages – Cas d'étude et outils numériques*.
> Techniques de l'Ingénieur, AG 6 062 (2027, à paraître).

**Authors**

- Olivier Vitrac — Innovation Lab, Adservio Group, Puteaux, France
- Phuong-Mai Nguyen — LNE (Laboratoire national de métrologie et d'essais), Trappes, France

The article is the third of the series *Conception sûre des emballages* (Techniques de l'Ingénieur):

| Part | Reference | Title |
|---|---|---|
| 1/3 | AG 6 058 (2026) | Principes, modélisation et fondements |
| 2/3 | AG 6 059 | Tiers avancés M2–M3 : partage, diffusion et barrières |
| 3/3 | AG 6 062 | Cas d'étude et outils numériques |

## Scripts

Each case isolates one design lever of safe packaging design and quantifies its effect on the severity $\hat{s}$
(compliance when $\hat{s} \leq 100\,\%$).

| Script | Case (section of AG 6 062) | Lever | Tier |
|---|---|---|---|
| [case_C1_reduce_C0.py](case_C1_reduce_C0.py) | C1 (4.1) — Diminuer $C_0$ | initial concentration, blending virgin and recycled PP | M0 |
| [case_C2_change_contaminant.py](case_C2_change_contaminant.py) | C2 (4.2) — Changer de contaminant | toxicological axis T (TTC, Cramer, specific SML) | T0–T2 |
| [case_C3_fatty_to_aqueous.py](case_C3_fatty_to_aqueous.py) | C3 (4.3) — Aliment gras → aqueux | partition coefficient $K$ | M2 |
| [case_C4_time_temperature.py](case_C4_time_temperature.py) | C4 (4.4) — Parcours temps/température | $D(T)$, Fourier number; reproduces table 16 | M3 |
| [case_C5_change_additive.py](case_C5_change_additive.py) | C5 (4.5) — Changer d'additif | diffusivity $D(M)$ (Piringer) | M3 |
| [case_C6_diffusion_barrier.py](case_C6_diffusion_barrier.py) | C6 (4.6) — Barrière à effet diffusion | wPET barrier thickness, lag time | M3 + FB |
| [case_C7_functional_barrier.py](case_C7_functional_barrier.py) | C7 (4.7) — Diffusivité et solubilité couplées | barrier permeance $D_j K_{i,j/s}/l_j$; split of $D$ and $K$ with `layerLink` | M2+M3 + FB |

The scripts reproduce the migration and severity calculations of the article. Estimates quoted from the literature
(for instance the free-volume and experimental diffusivities of case C5) are not recomputed.

## Run

From the SFPPy root:

```bash
python examples/TI_AG6062/case_C7_functional_barrier.py
```

Each script prints its results and a key message; none writes files. The values reported in the article were
computed with **SFPPy 1.9.3**.

## Cite

```text
Vitrac O., Nguyen P.-M. Conception sûre des emballages – Cas d'étude et outils numériques.
Techniques de l'Ingénieur, AG 6 062 (2027).
```

License: MIT, as SFPPy.
