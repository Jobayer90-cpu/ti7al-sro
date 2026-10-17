<div align="center">

# ti7al-sro

**Effect of Al short-range order on the stiffness and strength of Ti-7Al and estimating the order with Bayesian inference**

[![LAMMPS](https://img.shields.io/badge/LAMMPS-MEAM%20package-darkgreen)](https://www.lammps.org)
[![Python](https://img.shields.io/badge/Python-3.8%2B-blue)](https://www.python.org)
[![NumPy](https://img.shields.io/badge/NumPy-SciPy-orange)](https://scipy.org)
[![Matplotlib](https://img.shields.io/badge/Matplotlib-figures-yellow)](https://matplotlib.org)
[![OVITO](https://img.shields.io/badge/OVITO-visualization-purple)](https://www.ovito.org)
[![Status](https://img.shields.io/badge/status-ongoing-lightgrey)](#note)

</div>


## Table of Contents

- [Project Overview](#project-overview)
- [Clone the Repository](#clone-the-repository)
- [Workflow](#workflow)
- [Simulation Setup](#simulation-setup)
- [Running the Simulations](#running-the-simulations)
- [Running the Analysis](#running-the-analysis)
- [Results](#results)
- [Conclusions](#conclusions)
- [Project Organization](#project-organization)
- [Contact Information](#contact-information)
- [References](#references)
- [Citation](#citation)
- [Note](#note)

## Project Overview

Ti-7Al is a model alloy for the α phase of titanium alloys such as Ti-6Al-4V, which is used in jet engines for highly stressed rotating parts. In molecular dynamics (MD) simulations of such alloys, Al atoms are usually placed at random. In real alloys, however, the arrangement depends on ageing, temperature and processing history: Al atoms tend to avoid each other and have fewer Al neighbours than random placement would give. This is called **short-range order (SRO)**.

SRO is measured with the Warren–Cowley parameter α₁. A value of zero means a random arrangement, and more negative values mean more order. Ageing at moderate temperatures increases the order, so parts made from the same Ti-7Al but aged differently can have different strength, and models that assume random Al may overpredict the strength of aged material.

This project asks three questions:

1. How much does Al ordering change the stiffness and strength of Ti-7Al?
2. How large is the uncertainty in these properties, and where does it come from?
3. SRO is very hard to measure directly. Can the degree of order be estimated from an ordinary stiffness measurement, such as ultrasound?

## Clone the Repository

```bash
git clone https://github.com/Jobayer90-cpu/ti7al-sro.git
cd ti7al-sro
```

Requirements:

- [LAMMPS](https://www.lammps.org) with the MEAM package
- Python 3.8+ with NumPy, SciPy and Matplotlib

## Workflow

<p align="center">
  <img width="100%" alt="Workflow" src="workflow.png" />
</p>

## Simulation Setup

| Item | Setting |
|---|---|
| Material | Ti-7Al (11.78 at% Al), HCP |
| Potential | 2NN MEAM (Kim et al., 2016) |
| Ordering levels | α₁ = 0, -0.04, -0.08, -0.12 (limit at this composition: -0.134) |
| Boxes | 20 independent boxes of 4032 atoms per level, 80 in total, each with its own random seeds |
| Ordering method | swap Monte Carlo until α₁ is within 0.001 of the target |
| Relaxation | energy minimization, 5 ps NVT, 50 ps anisotropic NPT at 300 K and 0 bar, 1 fs time step |
| Elastic constants | 0 K, ±0.5% strain in all six directions |
| Tension | along [2-1-10] at 300 K and 10¹⁰ s⁻¹ up to 15% strain |
| Extra checks | three more strain rates (5×10¹⁰, 2×10⁹, 4×10⁸ s⁻¹) and 8 large boxes of 32,256 atoms |

The 20 boxes at one level have the same α₁ but different atom positions. Their scatter is the uncertainty studied here.

## Running the Simulations

1. Build the ordered structures and the LAMMPS run folders:

   ```bash
   python python/1_build/build_structures.py
   ```

2. Relax all boxes (Linux or WSL):

   ```bash
   cd sims && bash run_relax_all.sh
   ```

3. Check the relaxations and set up the mechanical tests:

   ```bash
   python python/3_check_relax/check_relax.py
   python python/4_prepare_mech/prepare_mech.py
   ```

4. Run the elastic constant and tension simulations:

   ```bash
   cd sims && bash run_mech_all.sh
   ```

A single box can also be run by hand:

```bash
cd sims/disordered/B1_relax/box01
mpirun -np 4 lmp -in in.relax
```

## Running the Analysis

These steps only need the csv files in `results/`, so they can be run without LAMMPS.

```bash
python python/6_extract/extract_results.py      # collect results from all runs
python python/7_analyze/analyze.py              # statistics and polycrystal UQ
python python/7b_inverse/infer_alpha.py         # Bayesian inverse model
python python/8_figures/make_figures_clean.py   # figures
```

## Results

### 1. Distributions at each level

<img width="100%" alt="fig1_distributions" src="https://github.com/user-attachments/assets/5c5ee2b9-b8d6-4741-835a-8a9b1bc3aeb6" />

_Figure 1_: Yield strength and tensile Young's modulus of all 80 boxes. Black: mean with 95% bootstrap interval.

| Level | Yield stress (GPa) | G (GPa) |
|---|---|---|
| α₁ = 0 | 6.90 ± 0.07 | 44.57 ± 0.11 |
| α₁ = -0.04 | 6.64 ± 0.06 | 43.97 ± 0.12 |
| α₁ = -0.08 | 6.30 ± 0.10 | 43.13 ± 0.08 |
| α₁ = -0.12 | 5.95 ± 0.08 | 42.02 ± 0.09 |

Normality was checked with 48 Shapiro–Wilk tests (12 properties at 4 levels), and 47 passed. A normal (Gaussian) model is therefore used in the UQ steps.

### 2. Effect of ordering on the properties

<img width="100%" alt="fig2_trends" src="https://github.com/user-attachments/assets/097e0865-4b85-4bb4-a43d-e029ece443b0" />

_Figure 2_: All 80 boxes against α₁, with a straight-line fit.

| Property | Change per -0.1 in α₁ | R² |
|---|---|---|
| Yield stress | -11.6% | 0.95 |
| C₆₆ | -7.4% | 0.98 |
| C₄₄ | -5.0% | 0.96 |
| G | -4.8% | 0.97 |
| E from tension | -3.7% | 0.34 |
| B | +1.0% | 1.00 |

- Ordering mainly lowers the **resistance to shear**. The resistance to compression (B) is almost unchanged.
- Between random and α₁ = -0.12, the means move by 12 to 27 box-to-box standard deviations (Welch t-test p < 10⁻³⁰; ANOVA across all four levels p ≈ 10⁻⁵⁰), so the shift is not scatter.
- The size of the scatter hardly changes with ordering. A Brown–Forsythe test finds no change for E, yield, G and C₆₆ (p = 0.38 to 0.81). C₄₄ is the only exception: its scatter shrinks, and the inverse model accounts for this.

> **In plain terms**
> - **Yield stress** is how hard you must pull before the metal deforms permanently. It drops by about 14% from random to strongly ordered.
> - **Shear stiffness** (G, C₄₄, C₆₆) is the resistance to being twisted or slid sideways. It drops by 6 to 9%.
> - **Bulk modulus** (B) is the resistance to being squeezed from all sides. It barely changes.
> - Boxes at the same ordering level still differ a little. That spread is the *scatter*. The jump between levels is 12 to 27 times larger than the scatter, so it comes from ordering, not from chance.

### 3. Where the scatter comes from

<img width="100%" alt="fig3_noise_vs_arrangement" src="https://github.com/user-attachments/assets/d38aeea6-5110-4fd2-a88c-798d454696cf" />

_Figure 3_: The same boxes measured two ways: Young's modulus from the tension test and from the elastic constants.

| Level | SD, tension E (GPa) | SD, E from elastic constants (GPa) | Share of tension variance that is noise |
|---|---|---|---|
| α₁ = 0 | 2.21 | 0.40 | 97% |
| α₁ = -0.04 | 2.78 | 0.31 | 99% |
| α₁ = -0.08 | 2.39 | 0.33 | 98% |
| α₁ = -0.12 | 3.10 | 0.45 | 98% |

The two measurements do not correlate from box to box (r ≈ 0), so the tension-test scatter is dominated by **thermal noise** from MD, not by where the atoms sit. For the UQ steps, stiffness was therefore taken from the **elastic constants**, which are about 7 times less scattered.

> **In plain terms: how the noise was measured**
>
> Young's modulus was obtained from the same 20 boxes in two ways:
> - **Tension test at 300 K:** the atoms vibrate, so the stress reading jitters. Box-to-box spread: 2.21 GPa.
> - **Elastic constants at 0 K:** no vibration, so almost no jitter. Spread: 0.40 GPa, which comes only from where the Al atoms sit.
>
> The tension test sees both effects. Removing the part caused by atom positions leaves the jitter (spreads add by their squares):
>
> noise² = 2.21² − 0.40², so noise ≈ 2.17 GPa, which is **97%** of the tension-test variance.
>
> The tension test mostly measures vibration noise, so the elastic constants were used instead.

### 4. From single crystal to polycrystal

<img width="100%" alt="fig4_polycrystal_uq" src="https://github.com/user-attachments/assets/d9c507a4-1f9a-47a9-9345-eadd9c0f34b5" />

_Figure 4_: (a) Polycrystal Young's modulus with 95% intervals. (b) Share of variance from crystal and texture. (c) Monte Carlo SD compared with a simple formula.

A real part is a polycrystal made of many grains with different orientations (texture). Its modulus was estimated as follows:

- **Crystal uncertainty:** a multivariate normal distribution fitted to the 20 sets of elastic constants at each level.
- **Texture:** 50 grain orientations with weights that sum to one. Two textures were assumed, not measured: *random*, where all orientations are equally likely, and *basal*, where the c-axis is mostly perpendicular to the load.
- **Texture uncertainty:** 5% random noise on each weight, renormalized so the weights still sum to one.
- **Monte Carlo:** 100,000 samples per case.

| Level | Texture | Polycrystal E (GPa) | Share of variance from crystal |
|---|---|---|---|
| α₁ = 0 | random | 118.15 ± 0.34 | 96% |
| α₁ = 0 | basal | 116.50 ± 0.37 | 99% |
| α₁ = -0.12 | random | 112.01 ± 0.24 | 90% |
| α₁ = -0.12 | basal | 109.18 ± 0.36 | 100% |

- Ordering lowers polycrystal E by **6 to 7 GPa** for both textures, about 20 standard deviations.
- Almost all of the spread (89 to 100%) comes from the single-crystal elastic constants. A 5% texture uncertainty adds very little.
- A simple formula that ignores the sum-to-one rule for the texture weights overestimates the SD by 2 to 4 times.

> **In plain terms: how a polycrystal value is obtained from a single crystal**
>
> No polycrystal was simulated. A titanium crystal is stiffer in some directions than in others: about 157 GPa when pulled along its c-axis and about 120 GPa when pulled in the basal plane. A real part contains millions of grains pointing in different directions, so:
> 1. the crystal's stiffness is computed for 50 grain directions;
> 2. these values are averaged, each weighted by how common that direction is (the texture);
> 3. to get the uncertainty, this calculation is repeated many times (Monte Carlo). Each time the elastic constants get a small random change matching the box-to-box scatter, and the texture weights get 5% random noise.
>
> The average of all repeats is the expected modulus, and their spread is its uncertainty.

<details>
<summary><b>Minimal example code (pure Python, runs in any online compiler)</b></summary>

```python
import math
import random

# ---------------------------------------------------------------
# 1. Single-crystal elastic constants from MD (GPa), random Al
#    mean and box-to-box SD of the 20 boxes
# ---------------------------------------------------------------
mean = {"C11": 177.74, "C12": 93.57, "C13": 79.39, "C33": 203.88, "C44": 40.51}
sd   = {"C11": 0.16,   "C12": 0.20,  "C13": 0.15,  "C33": 0.12,   "C44": 0.16}


def E_direction(C, theta):
    """Young's modulus of a hexagonal crystal pulled at angle theta
    from its c-axis (0 = along c, 90 deg = in the basal plane)."""
    C11, C12, C13, C33, C44 = C["C11"], C["C12"], C["C13"], C["C33"], C["C44"]
    D = C33 * (C11 + C12) - 2 * C13 ** 2
    S11 = (C11 * C33 - C13 ** 2) / ((C11 - C12) * D)
    S13 = -C13 / D
    S33 = (C11 + C12) / D
    S44 = 1 / C44
    s, c = math.sin(theta), math.cos(theta)
    return 1 / (S11 * s**4 + S33 * c**4 + (2 * S13 + S44) * s**2 * c**2)


# ---------------------------------------------------------------
# 2. Texture: 50 grain orientations and how common each one is
#    random texture -> every orientation equally likely
# ---------------------------------------------------------------
N = 50
angles = [math.acos(1 - (i + 0.5) / N) for i in range(N)]   # spread evenly
weights = [1 / N] * N                                       # sum = 1


def polycrystal_E(C, w):
    """weighted average of E over all grain orientations"""
    return sum(wi * E_direction(C, a) for wi, a in zip(w, angles))


print("Single crystal, along c-axis    :", round(E_direction(mean, 0), 1), "GPa")
print("Single crystal, in basal plane  :", round(E_direction(mean, math.pi / 2), 1), "GPa")
print("Polycrystal (one calculation)   :", round(polycrystal_E(mean, weights), 2), "GPa")

# ---------------------------------------------------------------
# 3. Monte Carlo: repeat with uncertain inputs
# ---------------------------------------------------------------
random.seed(1)
results = []
for run in range(20000):
    C = {k: random.gauss(mean[k], sd[k]) for k in mean}            # crystal scatter
    w = [wi * (1 + random.gauss(0, 0.05)) for wi in weights]       # 5 % texture noise
    total = sum(w)
    w = [wi / total for wi in w]                                   # weights sum to 1
    results.append(polycrystal_E(C, w))

m = sum(results) / len(results)
s = math.sqrt(sum((r - m) ** 2 for r in results) / (len(results) - 1))
print("Polycrystal, Monte Carlo        :", round(m, 2), "+/-", round(s, 2), "GPa")
```

Output:

```
Single crystal, along c-axis    : 157.4 GPa
Single crystal, in basal plane  : 120.1 GPa
Polycrystal (one calculation)   : 118.16 GPa
Polycrystal, Monte Carlo        : 118.16 +/- 0.23 GPa
```

This simplified version varies each elastic constant independently. The full analysis in `python/7_analyze/analyze.py` also uses the correlations between them, which gives a slightly larger spread (± 0.34 GPa).

</details>

### 5. Inverse model: how ordered is the alloy?

<img width="100%" alt="fig8_inverse" src="https://github.com/user-attachments/assets/b3fdb457-f9b2-4ea5-b784-7eb617f2a1c2" />

_Figure 5_: (a) Fitted curve predicting G from α₁ (the forward model); the grey band shows where 90% of boxes are expected to fall. (b) Estimated α₁ for three measured G values. (c) Width of the 90% interval against measurement error. (d) Leave-one-box-out check.

The idea is to measure a property such as the shear modulus G and let the model return the likely range of α₁.

| Part | Choice |
|---|---|
| Forward model | quadratic fit of each property against α₁ (80 boxes) |
| Scatter | box-to-box SD; constant, or changing with α₁ where Brown–Forsythe says so (only C₄₄) |
| Prior | α₁ equally likely anywhere between -0.134 and +0.01 |
| Likelihood | Gaussian: model scatter + fit uncertainty + measurement error |
| Result | posterior of α₁ on a fine grid, reported as a 90% interval |

**Example: measured G with 1% error**

| Measured G | Estimated α₁ (90% interval) | Meaning |
|---|---|---|
| 44.5 GPa | -0.044 to +0.010 | random or nearly random |
| 43.3 GPa | -0.101 to -0.035 | partly ordered |
| 42.1 GPa | -0.134 to -0.095 | strongly ordered |

**Which measurement is useful?** Width of the 90% interval for α₁ (smaller is better; with no information it would be 0.13).

| Property used | 0% error | 1% error | 2% error |
|---|---|---|---|
| G | 0.016 | 0.073 | 0.109 |
| C₆₆ | 0.016 | 0.046 | 0.087 |
| Yield stress | 0.034 | 0.043 | 0.063 |
| G + C₄₄ + C₆₆ | 0.014 | 0.035 | 0.068 |
| B | 0.004 | 0.127 | 0.129 |
| E from tension | 0.114 | 0.116 | 0.120 |

- Shear properties work well when measured to about 1%, which is possible with ultrasound. Combining G, C₄₄ and C₆₆ halves the interval.
- B changes too little and tension E is too noisy to be useful.
- In leave-one-box-out tests, the 90% interval contained the true α₁ in 85 to 90% of cases.

> **In plain terms: why the measurement must be accurate to about 1%**
>
> "Accurate to 1%" means the measured value is off by at most about 1%. For G ≈ 44 GPa that is about ±0.45 GPa. The whole effect of ordering on G is only about 6% (44.6 → 42.0 GPa), so the measurement error has to be clearly smaller than that window.
>
> | Error in G | ± in G | 90% interval for α₁ | What you learn |
> |---|---|---|---|
> | 0% | 0 | 0.016 | almost the exact α₁ |
> | 0.5% | ±0.22 GPa | 0.038 | very useful |
> | **1%** | **±0.45 GPa** | **0.073** | **random vs. partly vs. strongly ordered** |
> | 2% | ±0.9 GPa | 0.109 | only a rough hint |
> | 5% | ±2.2 GPa | ≈0.13 | nothing; the error is as big as the whole effect |
>
> 1% is not a fixed rule. It is roughly the largest error that still gives a useful answer. Ultrasonic methods (for example resonant ultrasound spectroscopy) can reach it; a normal tensile test cannot.
>
> **"Prior only"** in Figure 5(c) is what is known before any measurement: α₁ could be anywhere between -0.134 (fully ordered) and +0.01 (random), so the 90% range is about 0.13 wide. A measurement whose interval is as wide as this line teaches nothing. The further a curve lies below it, the more useful that measurement is.

### 6. Checks

<img width="100%" alt="fig5_size_check" src="https://github.com/user-attachments/assets/fece24aa-550b-4a52-ae71-d5fc6ed459d3" />

_Figure 6_: Box size. Large boxes (32,256 atoms, open diamonds) give the same means as the 4032-atom boxes (filled circles, mean ± SD), so only the scatter depends on box size.

<p align="center">
  <img width="45%" alt="fig6_strain_rate" src="https://github.com/user-attachments/assets/9c92a9d7-35ff-4f15-bb98-d2126608e1ca" />
  &nbsp;&nbsp;
  <img width="48%" alt="fig7_stress_strain" src="https://github.com/user-attachments/assets/59187187-b491-4741-9436-a51610453a9b" />
</p>

_Figure 7_: (left) Strain rate: the ordered alloy is weaker at every rate tested. (right) Stress–strain curves, average of 20 boxes per level.

## Conclusions

- Going from random to strongly ordered Al lowers the yield stress by about 14% and the shear stiffness by 6 to 9%, while the bulk modulus barely changes.
- These changes are far larger than the scatter between simulations, so they come from ordering and not from noise.
- Young's modulus from MD tension tests is dominated by thermal noise. Elastic constants are the better input for uncertainty analysis.
- At the polycrystal scale, almost all of the uncertainty comes from the single-crystal properties rather than the texture.
- A Bayesian model can estimate the degree of order from a measured shear modulus, with an error bar, if the measurement is accurate to about 1%.

## Project Organization

```
├── README.md                    <- README for users.
├── potentials                   <- 2NN MEAM potential files (library.meam, AlTi.meam).
├── templates                    <- LAMMPS input templates (relax, elastic constants, tension).
├── python
│   ├── common                   <- Shared settings and SRO tools (Warren–Cowley, swap Monte Carlo).
│   ├── 1_build                  <- Build ordered structures and run folders.
│   ├── 3_check_relax            <- Check the relaxations.
│   ├── 4_prepare_mech           <- Set up elastic constant and tension runs.
│   ├── 6_extract                <- Collect results from all runs.
│   ├── 7_analyze                <- Statistics and polycrystal uncertainty propagation.
│   ├── 7b_inverse               <- Bayesian inverse model.
│   └── 8_figures                <- Figures.
├── sims                         <- Batch scripts for LAMMPS runs.
└── results
    ├── *.csv                    <- Results tables.
    └── figures_clean            <- Figures shown in this README.
```

## Contact Information

- Md Jobayer Aziz: [md-jobayer.aziz1@louisiana.edu](mailto:md-jobayer.aziz1@louisiana.edu)

## References

- Interatomic potential: Y.-K. Kim, H.-K. Kim, W.-S. Jung, B.-J. Lee, *Atomistic modeling of the Ti–Al binary system*, Comput. Mater. Sci. 119 (2016) 1–8. [https://doi.org/10.1016/j.commatsci.2016.03.039](https://doi.org/10.1016/j.commatsci.2016.03.039)
- Random-alloy MD and polycrystal UQ for Ti-7Al: M.M. Billah, P. Acar, Acta Mater. 270 (2024) 119879. [https://doi.org/10.1016/j.actamat.2024.119879](https://doi.org/10.1016/j.actamat.2024.119879)
- SRO and planar slip in Ti–Al: R. Zhang et al., Sci. Adv. 5 (2019) eaax2799. [https://doi.org/10.1126/sciadv.aax2799](https://doi.org/10.1126/sciadv.aax2799)
- [LAMMPS](https://www.lammps.org): A.P. Thompson et al., Comput. Phys. Commun. 271 (2022) 108171. [https://doi.org/10.1016/j.cpc.2021.108171](https://doi.org/10.1016/j.cpc.2021.108171)
- [LAMMPS documentation](https://docs.lammps.org/Manual.html)

## Citation

If you use this code or data, please cite:

```bibtex
@misc{aziz2026ti7alsro,
  author       = {Aziz, Md Jobayer},
  title        = {Effect of Al short-range order on the stiffness and strength of Ti-7Al and estimating the order with Bayesian inference},
  year         = {2026},
  howpublished = {\url{https://github.com/Jobayer90-cpu/ti7al-sro}}
}
```

## Note

- Claude (Anthropic) was used as an assistant for code development, data analysis, figure preparation and editing of the text. All simulations, results and conclusions were checked by the author.
- This project is ongoing.
