## **🔬 🌊Predicting Partition Coefficients from $logP$ Using the Polarity Index $P'$ **  🍏⏩🍎

### **📌 Synopsis**  
Directly substituting partition coefficients or Henry-like constants ($k$) with **$logP$ values** is ❌ **theoretically unsound**, as no established framework supports such an approach.  

This document presents a **✅ generalized method** based on **Flory-Huggins theory** 📖, introducing a **robust polarity index ($P'$) 🔢**. By deriving **$P'$ from $logP$ and molecular mass ($M$) ⚖️**, this approach provides a **🔍 more consistent and predictive** framework for estimating partition coefficients **across different systems 🌎**.

***

[TOC]

***

### Objective

The goal is to establish a qualitative correlation between an approximate “polarity index” ($P’$) and the octanol-water partition coefficient ($logP$) for small organic molecules ($i$), using the Flory-Huggins (FH) theory of mixtures. The polarity index $P’$ is introduced as an alternative to the square root of molar cohesive energies $\sqrt{E}$ in the approximation of FH interaction parameters for a solute $i$ in a fluid or polymer phase $k$:

$$
         \chi_{i+k} \approx 0.34 \times (\sqrt{E_i} - \sqrt{E_k}) ^ 2
$$

where $E_i$ and $E_k$ are the cohesive energies of $i$ and $k$, respectively. The empirical coefficient 0.34 is derived from experimental data.

>> References
>> 
>> -	Guillaume Gillet, Olivier Vitrac, and Stéphane Desobry
>>	Prediction of Solute Partition Coefficients between Polyolefins and Alcohols Using a Generalized Flory−Huggins Approach
>>	*Industrial & Engineering Chemistry Research* **2009** 48 (11), 5285-5301
>>	https://doi.org/10.1021/ie801141h
>>
>> - Thomas Brouwer and Boelo Schuur
>>	Model Performances Evaluated for Infinite Dilution Activity Coefficients Prediction at 298.15 K
>>	*Industrial & Engineering Chemistry Research* **2019** 58 (20), 8903-8914
>>	https://doi.org/10.1021/acs.iecr.9b00727

***

### Derivation of the Polarity Index ($P’$) from logP
Since $logP$ values include a significant entropic contribution, particularly for large substances, this effect is removed using the same FH theory.
For solute $i$, $logP$ is related to the octanol-water partition coefficient $lnK_{i,o/w}$ at infinite dulution:

$$
logP = \frac{lnK_{i,o/w}}{ln10}
$$

where:

$$
lnK_{i,o/w} = \gamma_{i,w} - \gamma_{i,o}
$$

At infinite dilution in $k = o, w$, the FH theory gives the following expression for the activity coefficient of $i$ in phase $k$:

$$
\gamma_{i,k} \approx \chi_{i+k} +1 - r_{i,k}
$$

where $r_{i,k} \approx \frac{V_i}{V_k}$ is the molar volume ratio between $i$ and $k$. As a result:

$$
ln10 \times logP \approx  \Delta\chi_{i}^{scale:ow} + S_i^{correction:ow}
$$

where the correction term accounts for entropic effects:

$$
S_i^{correction:ow} = -\left( r_{i,w} - r_{i,o} \right)
$$

By generalizing the FH coefficient approximation with a scaling parameter $\alpha$ and substituting cohesive energies with their corresponding polarity indices $P’$, we obtain:

$$
\Delta\chi_{i}^{scale:ow} = \chi_{i,w} - \chi_{i,o}
   \approx \alpha \times \left[ (P'_i - P'_w) ^ 2 - (P'_i - P'_o) ^ 2 \right] \\
$$

which expands to:

$$
= \alpha \times \left[P'_i \cdot (P'_w - P'_o) + (P'_w + P'_o)(P'_w - P'_o) \\ \right] \\
   = \alpha\cdot(P'_w - P'_o)\times\left(P'_i + \beta \right) \\
   = A\cdot P'_i + B
$$

where $A$ and $B$ are constants.

> Since $P’$ is defined based on differences, we can set $B = 0$ without affecting predictions.

***

### Definition and Range of $P’$
-	$P’$ is designed to be a bounded measure of polarity, ranging from 0 (highly hydrophobic substances) to 10.2 (water, the most polar reference).
-	It follows that:

$$
A \cdot P' \approx ln10\times logP - S_i^{correction:ow}(V_i)  \approx U\left(logP,V_i\right)
$$

where choosing $A = 1$ leads to:

$$
\alpha = \frac{1}{(P'_w - P'_o)}
$$

> Empirical fitting suggests $\alpha \approx 0.14$.
> To ensure a smooth and bounded prediction of $P’$, we approximate it using a second-degree polynomial regression:

$$
\hat{P’}(M,V) = a \cdot  \left[U\left(logP,V_i=f(M)\right)\right]^2 + b \cdot  U\left(logP,V_i=f(M)\right) + c
$$

where $V_i$ is estimated using Miller’s empirical formula:

$$
V_i = f\left(M_i\right) \approx 0.997 \cdot M_i^{1.03}
$$

with $M_i$ in $g \cdot mol^{-1}$ and $V_i$ in $cm^3 \cdot mol^{-1}$.

***

### Fitting Methodology
- The parameters $a$, $b$, and $c$ are optimized based on reported $logP$ and $P’$ values from open databases.
- The upper bound of $P’$ is determined by water’s polarity, which depends on database values for $logP$ and parameter choices.

***

### Application to Activity Coefficients
Once $P’$ is derived, it can be used to estimate activity coefficients:

$$
\gamma_{i,k} \approx \chi_{i+k} +1 - \left(r_{i,k} - n(r_{i,k})\right) \\
\approx \alpha \times [\hat{P’}(logP_i,V_i) - \hat{P’}(logP_k,V_k)] ^ 2  + 1 - \left(r_{i,k} - n(r_{i,k})\right) 
$$

The correction term $n(r_{i,k})$ accounts for the presence of surrounding molecules within the effective volume of the solute. This effect is significant for large, flexible solutes and negligible for small, rigid molecules.

***

### Special Cases
- For polymers ($k=P$):
Since the solute is considered infinitely small compared to the polymer, $r_{i,k}  \rightarrow 0$ and $n(r_{i,k}) \rightarrow 0$.
- For food simulants ($k=F$): For a generalized approximation applicable to various substances, we propose:

$$
n(r) = \frac{r-5}{r}
$$

>> References
>
>>		Guillaume Gillet, Olivier Vitrac, and Stéphane Desobry
>>	Prediction of Partition Coefficients of Plastic Additives between Packaging Materials and Food Simulants
>>	*Industrial & Engineering Chemistry Research* **2010** 49 (16), 7263-7280
>>	https://doi.org/10.1021/ie9010595

***

### Summary
✔ The methodology is theoretically well-founded and aligns with FH theory.

✔ It provides a practical way to estimate activity coefficients without requiring direct cohesive energy measurements.

✔ It accounts for entropy corrections, making $P’$ a more transferable descriptor.

✔ The polynomial fit is simple and efficient, but its generalization across all chemical systems remains to be tested.

- $P’$ is introduced as a polarity index replacing cohesive energy terms in FH calculations.
- It is directly linked to $logP$ and fitted using second-degree polynomial regression.
	•	The method provides a universal approach to predicting activity coefficients in various media, including polymers and food simulants.
	•	A correction term $n(r)$ accounts for molecular interactions in complex environments.

***
### Extensions to $k_{i,k}^H$ (Henry-like coefficients)

Calculated activity coefficients, $\gamma_{i,k}$ are associated with volume fractions $\phi_{i,k}$ in the amorphous phase. As a result they are not practical to expression mass fluxes within an effective medium. In volume-averaged transport equations, Henry-like partition coefficients, $k_{i,k}^H$ are preferred. They express the linear relationship at infinite dilution between the partial pressure $p_{i}$ (or fugacity) in phase/component $k$ (pressure of a theoretical gas phase in equilibrium with the local volume) and the volume-averaged concentration of the species $i$ in $k$:
$$
p_i = k_{i,k}^HC_{i,k}\cdot\left(1-c_k\right)\cdot\left(1-\epsilon_k\right)=\gamma_{i,k}\phi_{i,k}P_{i,sat}^{(T)}
$$
where $c_k$ and $\epsilon_k$ are the crystallinity and porosity of the effective medium $k$, and $P_{i,sat}^{(T)}$ the vapor saturation pressure of $i$.

Finally by noticing that $C_{i,k}\cdot\left(1-c\right)\cdot\left(1-\epsilon\right)\times V_i^m=\phi_{i,k}$ with $V_i^m$ the molar volume of $i$, one gets: 
$$
k_{i,k} =\frac{P_{i,sat}^{(T)}V_i^m}{\left(1-c_k\right)\cdot\left(1-\epsilon_k\right)}\times\gamma_{i,k}
$$
The only assumption made here is that the sorption isotherm is linear and that the volatility of the substances is low to neglect the amount of substances in the gas phase.

> This model is fully implemented in `SFPPy` within the `kFH` model. For non-volatile substances and mass transfer through dense phases, it is convenient to scale all $k_{i,k}^H$ and consequently all fugacities with $P_{i,sat}^{(T)}=1$.
> The correct $P_{i,sat}^{(T)}$ value needs to be assigned in the presence of an air layer or of a porous medium since $k_{i,air}$ is by definition $\frac{1}{RT}$.


***

### Outlook
⚠ More validation is needed, especially for extreme cases (very high or low $logP$) and polymers.

⚠ Alternative functional forms for $P’$ may improve accuracy.

⚠ The empirical scaling factor ($\alpha$) and correction function ($n(r)$) should be tested for robustness.

***

<div style="border: 2px solid #4CAF50; border-radius: 8px; padding: 10px; background: linear-gradient(to right, #4CAF50, #FF4D4D); color: white; text-align: center; font-weight: bold;">
  <span style="font-size: 20px;">🍏⏩🍎 <strong>SFPPy for Food Contact Compliance and Risk Assessment</strong></span><br>
  Contact <a href="mailto:olivier.vitrac@gmail.com" style="color: #fff; text-decoration: underline;">Olivier Vitrac</a> for questions |
  <a href="https://github.com/ovitrac/SFPPy" style="color: #fff; text-decoration: underline;">Website</a> |
  <a href="https://ovitrac.github.io/SFPPy/" style="color: #fff; text-decoration: underline;">Documentation</a>
</div>