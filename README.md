# Effect of Al short-range order on the stiffness and strength of Ti-7Al and estimating the order with Bayesian inference

Ti-7Al is the α phase of titanium alloys such as Ti-6Al-4V. When modeling such system for molecular dynamics (MD) simulations, Al atoms are often placed at random. However, the ordering varies depending on factors like ageing, temperature or processes the alloy has been subjected to. In Ti-7Al, the Al atoms generally have Ti as their nearest neighbors rather than other Al atoms. This is called short-range order (SRO). SRO is measured using Warren Cowley parameter. Warren Cowley parameter value of zero indicates disorder, while the more negative it is the higher the degree of order. Ageing in moderate temperature increases the degree of order, so its important to know how mechanical properties in Ti-7Al varies with varying SRO.
<br> <br>
Same mechanical parts made from the same Ti-7Al but aged differently can have different strength. Models that assume random Al may overpredict the strength of aged material. In this project we tried to quantify much does Al ordering change the stiffness and strength of Ti-7Al. We explored the uncertainity in these properties and where they come from. 
<br> <br>
Another problem happens to be that SRO is very hard to measure directly. So, we attempted to make a model that will be able to infer how ordered the alloy is based on the mechanical property values of that alloy. This inverse model will be able an ordinary stiffness measurement (for example ultrasound) into an estimate of how ordered the alloy is and help us understanding its physical behavior. 


---

## Methodology

```
Build boxes with a chosen amount of order (α₁)
        ↓
MD simulations (LAMMPS): relaxation, elastic constants, tensile test
        ↓
Statistics: distributions, trends, tests on mean and spread, noise
        ↓
Forward UQ: single-crystal uncertainty to a polycrystal 
        ↓
Inverse UQ: estimate α₁ from a measured property
```

### Warren-Cowley parameter α₁

α₁ describes the first shell of neighbours around each atom.

- **α₁ = 0:** Disorder
- **α₁ < 0:** Order
Four α₁ levels were used for this study: **0, -0.04, -0.08 and -0.12**.

### Simulation setup

- **Material:** Ti-7Al (11.78 at% Al), HCP.
- **Potential:** 2NN MEAM (Kim et al., 2016).
- **Boxes:** 20 independent boxes of 4032 atoms at each level, 80 in total. Each box has its own random seeds. Ordering was introduced by swap Monte Carlo until α₁ was within 0.001 of the target.
- **Relaxation:** energy minimization, then 5 ps NVT and 50 ps anisotropic NPT at 300 K and 0 bar, with a 1 fs time step.
- **Elastic constants:** at 0 K, from ±0.5% strains in all six directions.
- **Tension:** along [2-1-10] at 300 K and 10¹⁰ s⁻¹ up to 15% strain.
- **Extra checks:**
  - three more strain rates (5×10¹⁰, 2×10⁹ and 4×10⁸ s⁻¹);
  - 8 large boxes of 32,256 atoms.

The 20 boxes at one level have the same α₁ but different atom positions. Their scatter is the uncertainty we study.

---

## Results

### 1. Distributions at each level

<img width="2153" height="874" alt="fig1_distributions" src="https://github.com/user-attachments/assets/5c5ee2b9-b8d6-4741-835a-8a9b1bc3aeb6" />



| Level | Yield stress (GPa) | G (GPa) |
|---|---|---|
| α₁ = 0 | 6.90 ± 0.07 | 44.57 ± 0.11 |
| α₁ = -0.04 | 6.64 ± 0.06 | 43.97 ± 0.12 |
| α₁ = -0.08 | 6.30 ± 0.10 | 43.13 ± 0.08 |
| α₁ = -0.12 | 5.95 ± 0.08 | 42.02 ± 0.09 |

Normality was checked with 48 Shapiro-Wilk tests: 12 properties at 4 levels. 47 of them passed. A normal (Gaussian) model is therefore used in the UQ steps.

### 2. Effect of ordering on the properties

<img width="2202" height="1415" alt="fig2_trends" src="https://github.com/user-attachments/assets/097e0865-4b85-4bb4-a43d-e029ece443b0" />


*All 80 boxes against α₁, with a straight-line fit.*

| Property | Change per -0.1 in α₁ | R² |
|---|---|---|
| Yield stress | -11.6% | 0.95 |
| C₆₆ | -7.4% | 0.98 |
| C₄₄ | -5.0% | 0.96 |
| G | -4.8% | 0.97 |
| E from tension | -3.7% | 0.34 |
| B | +1.0% | 1.00 |

Ordering mainly lowers the **resistance to shear**. Resistance to compression (B) is almost unchanged.
Between random and α₁ = -0.12, the means move by 12 to 27 box-to-box standard deviations (Welch t-test p < 10⁻³⁰; ANOVA across all four levels p ≈ 10⁻⁵⁰). So the shift is not scatter. A Brown-Forsythe test finds no change for E, yield, G and C₆₆ (p = 0.38 to 0.81). C₄₄ is the only exception: its scatter shrinks, and the inverse model accounts for this. So, the size of the scatter hardly changes with ordering.

### 3. Analysis of the scatter

<img width="2217" height="870" alt="fig3_noise_vs_arrangement" src="https://github.com/user-attachments/assets/d38aeea6-5110-4fd2-a88c-798d454696cf" />

*The same boxe’s properties were measured two ways. One is Young's modulus from the tension test and another from the elastic constants.* Young's modulus from tension has an R² of only 0.34, so its scatter was examined further.

| Level | SD, tension E (GPa) | SD, E from elastic constants (GPa) | Share of tension variance that is noise |
|---|---|---|---|
| α₁ = 0 | 2.21 | 0.40 | 97% |
| α₁ = -0.04 | 2.78 | 0.31 | 99% |
| α₁ = -0.08 | 2.39 | 0.33 | 98% |
| α₁ = -0.12 | 3.10 | 0.45 | 98% |

The two measurements do not correlate from box to box (r ≈ 0). So the tension-test scatter is dominated by **thermal noise** from MD, not an effect of where the atoms sit. So, for the UQ steps, stiffness was taken from the **elastic constants**, which are about 7 times less scattered.

### 4. From single crystal to polycrystal

<img width="3032" height="901" alt="fig4_polycrystal_uq" src="https://github.com/user-attachments/assets/d9c507a4-1f9a-47a9-9345-eadd9c0f34b5" />


*(a) Polycrystal Young's modulus with 95% intervals. (b) Share of variance from crystal and texture. (c) Monte Carlo SD compared with a simple formula.*

A real part is a polycrystal made of many grains with different orientations (texture). Its modulus was estimated as follows:

- **Crystal uncertainty:** a multivariate normal distribution fitted to the 20 sets of elastic constants at each level.
- **Texture:** 50 grain orientations with weights that sum to one. Two textures were assumed, not measured:
  - *random*, where all orientations are equally likely;
  - *basal*, where the c-axis is mostly perpendicular to the load.
- **Texture uncertainty:** 5% random noise on each weight, renormalized so the weights still sum to one.
- **Monte Carlo:** 100,000 samples per case.

| Level | Texture | Polycrystal E (GPa) | Share of variance from crystal |
|---|---|---|---|
| α₁ = 0 | random | 118.15 ± 0.34 | 96% |
| α₁ = 0 | basal | 116.50 ± 0.37 | 99% |
| α₁ = -0.12 | random | 112.01 ± 0.24 | 90% |
| α₁ = -0.12 | basal | 109.18 ± 0.36 | 100% |

Ordering lowers polycrystal E by **6 to 7 GPa** for both textures, about 20 standard deviations. Almost all of the spread (89 to 100%) comes from the single-crystal elastic constants. A 5% texture uncertainty adds very little.
A simple formula that ignores the sum-to-one rule for the texture weights overestimates the SD by 2 to 4 times.

### 5. Inverse model: how ordered is the alloy?

<img width="2154" height="1712" alt="fig8_inverse" src="https://github.com/user-attachments/assets/b3fdb457-f9b2-4ea5-b784-7eb617f2a1c2" />


*(a) Forward model for G with 90% band. (b) Estimated α₁ for three measured G values. (c) Width of the 90% interval against measurement error. (d) Leave-one-box-out check.*

The idea is to measure a property such as the shear modulus G, and let the model return the likely range of α₁.

| Part | Choice |
|---|---|
| Forward model | quadratic fit of each property against α₁ (80 boxes) |
| Scatter | box-to-box SD; constant, or changing with α₁ where Brown-Forsythe says so (only C₄₄) |
| Prior | α₁ equally likely anywhere between -0.134 and +0.01 |
| Likelihood | Gaussian: model scatter + fit uncertainty + measurement error |
| Result | posterior of α₁ on a fine grid, reported as a 90% interval |

**Example: measured G with 1% error**

| Measured G | Estimated α₁ (90% interval) | Meaning |
|---|---|---|
| 44.5 GPa | -0.044 to +0.010 | random or nearly random |
| 43.3 GPa | -0.101 to -0.035 | partly ordered |
| 42.1 GPa | -0.134 to -0.095 | strongly ordered |


| Property used | 0% error | 1% error | 2% error |
|---|---|---|---|
| G | 0.016 | 0.073 | 0.109 |
| C₆₆ | 0.016 | 0.046 | 0.087 |
| Yield stress | 0.034 | 0.043 | 0.063 |
| G + C₄₄ + C₆₆ | 0.014 | 0.035 | 0.068 |
| B | 0.004 | 0.127 | 0.129 |
| E from tension | 0.114 | 0.116 | 0.120 |

Shear properties work well when measured to about 1% (possible with ultrasound). Combining G, C₄₄ and C₆₆ halves the interval. B changes very little and tension E is too noisy. So, E from tensile test is not useful.

<p float="left">
<img width="2154" height="874" alt="fig5_size_check" src="https://github.com/user-attachments/assets/fece24aa-550b-4a52-ae71-d5fc6ed459d3" />
<img width="1001" height="863" alt="fig6_strain_rate" src="https://github.com/user-attachments/assets/9c92a9d7-35ff-4f15-bb98-d2126608e1ca" />
<img width="1076" height="881" alt="fig7_stress_strain" src="https://github.com/user-attachments/assets/59187187-b491-4741-9436-a51610453a9b" />
</p>

---
## Conclusions

- Going from random to strongly ordered Al lowers the **yield stress by about 14%** and the **shear stiffness by 6 to 9%**, while the bulk modulus barely changes.
- These changes are far larger than the scatter between simulations, so they are real.
- A Bayesian model can estimate the degree of order from a measured shear modulus, with an error bar, if the measurement is accurate to about 1%.
