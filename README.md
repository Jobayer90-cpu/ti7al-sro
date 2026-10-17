# ti7al-sro

**Effect of Al short-range order on the stiffness and strength of Ti-7Al and estimating the order with Bayesian inference**

[![LAMMPS](https://img.shields.io/badge/LAMMPS-MEAM-darkgreen)](https://www.lammps.org)
[![Python](https://img.shields.io/badge/Python-3.10-blue)](https://www.python.org)
[![Status](https://img.shields.io/badge/status-ongoing-lightgrey)](#note)

## Why this project

Ti-7Al is titanium with 7 wt% aluminium. It is often used as a simple version of Ti-6Al-4V, the alloy in jet engine parts.

In MD simulations, the Al atoms are usually placed at random. In a real alloy they are not. When the alloy is aged (heat treated around 500 to 600 °C), Al atoms start to avoid each other. This is called short-range order (SRO). The crystal stays the same, only the arrangement of Ti and Al atoms changes.

I measure SRO with the Warren–Cowley parameter α₁. α₁ = 0 means random. More negative means more ordered.

I wanted to know:

1. How much does SRO change stiffness and strength?
2. How much do the results scatter, and why?
3. Can we estimate SRO from a stiffness measurement?

## Main findings

- From random (α₁ = 0) to ordered (α₁ = -0.12), yield strength drops by about 14% and shear stiffness by 6 to 9%. Bulk modulus barely changes.
- These changes are much bigger than the scatter between simulations, so they are real.
- Young's modulus from the tensile test is mostly thermal noise. Getting it from the elastic constants is about 7 times more precise.
- For a polycrystal, Young's modulus drops by 6 to 7 GPa with ordering.
- If the shear modulus is measured within about 1%, a Bayesian model can estimate the degree of order with an error bar.

## What I did

<p align="center">
  <img width="100%" alt="Workflow" src="workflow.png" />
</p>

- **Structures:** 4 levels of order (α₁ = 0, -0.04, -0.08, -0.12). 20 boxes of 4032 atoms at each level. Same order, different atom arrangement. Order was made by swapping Ti and Al atoms.
- **Simulations (LAMMPS, 2NN MEAM potential):**
  - relaxation at 300 K and 0 bar (5 ps NVT + 50 ps NPT)
  - elastic constants at 0 K (±0.5% strain)
  - tensile test at 300 K, 10¹⁰ s⁻¹, up to 15% strain
- **Checks:** three more strain rates, and 8 bigger boxes (32,256 atoms).
- **Polycrystal:** averaged the crystal stiffness over many grain directions and used Monte Carlo to get the uncertainty.
- **Inverse model:** Bayesian model that goes from a measured stiffness back to α₁.

## Results

### 1. Effect of ordering

<img width="100%" alt="fig2_trends" src="https://github.com/user-attachments/assets/097e0865-4b85-4bb4-a43d-e029ece443b0" />

All 80 boxes against α₁. Ordering increases to the right.

| Property | Change per -0.1 in α₁ |
|---|---|
| Yield strength | -11.6% |
| C₆₆ (shear) | -7.4% |
| C₄₄ (shear) | -5.0% |
| Shear modulus G | -4.8% |
| Bulk modulus B | +1.0% |

Ordering makes the alloy weaker and easier to shear, but not easier to compress.

The simulated yield strength (about 6 to 7 GPa) is much higher than in experiments because the crystals are perfect and the strain rate is very high. The change between levels is what matters.

### 2. Scatter and noise

<img width="100%" alt="fig3_noise_vs_arrangement" src="https://github.com/user-attachments/assets/d38aeea6-5110-4fd2-a88c-798d454696cf" />

Young's modulus of the same boxes found in two ways.

I got Young's modulus in two ways:

- **Tensile test at 300 K:** atoms vibrate, so the stress signal is noisy. SD = 2.21 GPa.
- **Elastic constants at 0 K:** no vibration. SD = 0.40 GPa, coming only from where the Al atoms sit.

Spreads add as squares, so the noise is √(2.21² − 0.40²) ≈ 2.17 GPa. That is 97 to 99% of the tensile-test scatter. So I used the elastic constants for the rest of the work.

### 3. From single crystal to polycrystal

<img width="100%" alt="fig4_polycrystal_uq" src="https://github.com/user-attachments/assets/d9c507a4-1f9a-47a9-9345-eadd9c0f34b5" />

(a) Polycrystal Young's modulus. (b) Where the uncertainty comes from. (c) Monte Carlo vs a simple formula.

I did not simulate a polycrystal. A Ti crystal is stiffer along its c-axis (157 GPa) than across it (120 GPa). A real part has grains in all directions, so I:

1. calculated the crystal stiffness in 50 directions,
2. averaged them using texture weights (how common each direction is),
3. repeated this 100,000 times with small random changes in the elastic constants and 5% noise in the weights.

| | Random texture | Basal texture |
|---|---|---|
| α₁ = 0 | 118.15 ± 0.34 GPa | 116.50 ± 0.37 GPa |
| α₁ = -0.12 | 112.01 ± 0.24 GPa | 109.18 ± 0.36 GPa |

Almost all of the uncertainty (89 to 100%) comes from the crystal, not the texture. The textures were assumed, not measured.

### 4. Estimating order from a measurement

<img width="100%" alt="fig8_inverse" src="https://github.com/user-attachments/assets/b3fdb457-f9b2-4ea5-b784-7eb617f2a1c2" />

(a) G against α₁. (b) Estimated α₁ for three measured G values. (c) Effect of measurement error. (d) Test on boxes the model had not seen.

The model fits G against α₁, then works backwards: for a measured G, it gives the likely range of α₁.

| Measured G (±1%) | Likely α₁ (90% range) | Meaning |
|---|---|---|
| 44.5 GPa | -0.044 to +0.010 | random |
| 43.3 GPa | -0.101 to -0.035 | partly ordered |
| 42.1 GPa | -0.134 to -0.095 | strongly ordered |

G only changes by about 6% over the whole range, so the measurement must be accurate. At 1% error (±0.45 GPa) the model still works. At 5% it tells you nothing. Ultrasonic testing can reach about 1%, a normal tensile test cannot.

To test it, I hid each box and estimated its α₁ from the other 79. The true value was inside the 90% range in 85 to 90% of cases.

### 5. Checks

<img width="100%" alt="fig5_size_check" src="https://github.com/user-attachments/assets/fece24aa-550b-4a52-ae71-d5fc6ed459d3" />

Bigger boxes give the same averages. Only the scatter depends on box size.

<p align="center">
  <img width="45%" alt="fig6_strain_rate" src="https://github.com/user-attachments/assets/9c92a9d7-35ff-4f15-bb98-d2126608e1ca" />
  &nbsp;&nbsp;
  <img width="48%" alt="fig7_stress_strain" src="https://github.com/user-attachments/assets/59187187-b491-4741-9436-a51610453a9b" />
</p>

Left: the ordered alloy is weaker at every strain rate. Right: average stress–strain curve at each level.

## Conclusion

- SRO makes Ti-7Al weaker and easier to shear. Yield strength drops by about 14% from random to ordered.
- Simulations that assume random Al may overpredict the strength of aged material.
- For uncertainty work, take stiffness from the elastic constants, not from the MD tensile test.
- A Bayesian inverse model can estimate the degree of order from a measured shear modulus, with an error bar. It works when the measurement is accurate to about 1%, and its error bars held up in tests on unseen boxes.

## Limitations

- One interatomic potential.
- Elastic constants at 0 K only.
- Order was made by swapping atoms, not by simulating real ageing.
- Textures were assumed, and grain boundaries are ignored.
- The inverse model is tested on simulated data only.

## How to run

Needs Python 3.8+ with NumPy, SciPy and Matplotlib. LAMMPS with the MEAM package is only needed to rerun the simulations.

```bash
git clone https://github.com/Jobayer90-cpu/ti7al-sro.git
cd ti7al-sro
pip install -r requirements.txt
```

**Redo the analysis** (about one minute, uses the saved results in `results/`):

```bash
python python/7_analyze/analyze.py
python python/7b_inverse/infer_alpha.py
python python/8_figures/make_figures_clean.py
```

**Redo the simulations** (needs a cluster):

```bash
python python/1_build/build_structures.py
cd sims && bash run_relax_all.sh && cd ..
python python/3_check_relax/check_relax.py
python python/4_prepare_mech/prepare_mech.py
python python/4_prepare_mech/prepare_mech.py --batch B4_rate B5b_tension80 B5c_cij80
cd sims && bash run_mech_all.sh && cd ..
python python/6_extract/extract_results.py
```

## Folders

```
potentials/   MEAM potential files
templates/    LAMMPS input files
python/       scripts, run in the order of the folder numbers
sims/         batch scripts for LAMMPS
results/      results.csv, tension.csv, cij.csv and figures
```

## References

1. Kim et al., Comput. Mater. Sci. 119 (2016) 1–8. (potential)
2. Billah and Acar, Acta Mater. 270 (2024) 119879.
3. Zhang et al., Sci. Adv. 5 (2019) eaax2799.

## Contact

Md Jobayer Aziz, md-jobayer.aziz1@louisiana.edu

## Note

This project is ongoing. Claude (Anthropic) was used as an assistant for coding, analysis, figures and editing. All results were checked by me.

