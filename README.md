# Simulation Code: A Background-Free Search for Physics Beyond the Standard Model Using Atom Interferometry
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18896187.svg)](https://doi.org/10.5281/zenodo.18896187)

This repository contains the advanced MATLAB simulation suite used to generate the data, projections, and figures for the paper:

**A Background-Free Search for Physics Beyond the Standard Model Using Atom Interferometry**
*Ivan Cagnani (Linnaeus University)*

## Overview

Extracting a μrad-scale Beyond Standard Model (BSM) phase shift from an atom interferometer requires overcoming gigaradian-scale gravitational backgrounds, massive vibration noise, and the fundamental kinematic degeneracy between scalar BSM couplings ($\propto \phi_g \phi_e$) and the quadratic Stark shift.

This code rigorously validates the experimental protocol proposed in the paper. To avoid the superficiality of simplified algebraic models and the prohibitive $\mathcal{O}(N^3)$ computational bottleneck of massive matrix exponentiation, this simulator is built on a high-speed, two-stage architecture: an **ab initio physics engine** followed by a **statistical noise engine**.

## Methodology

The script executes two primary computational tasks:

### Part 1: The Fast SSFM Physics Engine
Because the $k$-odd parity of the BSM and Stark operators is strictly dependent on the physical spatial separation of the atomic wavepackets, it must be evaluated quantum mechanically. 
* The engine natively solves the 1D time-dependent Schrödinger equation using the **Split-Step Fourier Method (SSFM)**.
* To prevent catastrophic momentum aliasing from the macroscopic gravitational acceleration ($g = 9.81$ m/s²), the simulation utilizes a massive spatial grid of **$N = 1,048,576$ ($2^{20}$) points**. 
* By alternating between position and momentum space via Fast Fourier Transforms (FFTs), it explicitly resolves the physical spatial separation ($\Delta z \propto k_{\text{eff}}$) and establishes the true theoretical baseline differential phases.

### Part 2: Statistical Monte Carlo Noise Engine
Using the true quantum phases established by the SSFM, the Monte Carlo engine simulates the signal extraction under severe experimental noise conditions over 25,000 complete 4-point measurement cycles (representing a 17-day integration time).
* **Mid-Fringe Spin-Squeezing:** It injects massive $2\pi$ common-mode vibration noise onto $N=10^6$ optimally spin-squeezed atoms. It utilizes active phase sweeping and a $\pi/2$ optical offset to geometrically open the covariance trace into a perfect circle, mathematically defeating *fringe-peak rectification bias* (ellipse-collapse bias).
* **Lock-In Quadrature Validation:** It injects a severe macroscopic low-frequency differential drift (~357 μrad). The simulation proves that the Double-Difference lock-in successfully algebraically annihilates the gigaradian gravity background and rejects the wandering drift to perfectly isolate the target BSM/Stark parity signal ($>10\sigma$ confidence) at the quantum limit.

## Output Visualizations

![Atom Interferometer Dashboard](BSM_Search_Plots_Ivan_Cagnani_2026.png)

Running the script automatically generates a publication-ready **3-Panel Dashboard**:
1. **Ab Initio Spatial Separation:** Visualizes the probability density $|\psi(x)|^2$ at $t=T$, proving the fundamental $k$-odd parity of the spatial potentials.
2. **Noise Rejection (Mid-Fringe):** Plots the active phase sweep of the simultaneous dual-isotope sequence, demonstrating unbiased covariance extraction in high-noise regimes.
3. **Lock-In Quadrature Isolation:** A histogram of the 25,000 lock-in outputs, proving the successful isolation of the target nanoradian signal from the annihilated gravitational background.

## Code Files

**`BSM_Search_Simulator_Ivan_Cagnani_2026.m`**
* **This is the primary, public-facing script for the project.**
* It is a plain-text MATLAB file, fully commented, and self-contained (phase extraction functions are embedded).
* Highly optimized, it completes a massive quantum simulation and 25,000 statistical cycles in under a few minutes on standard hardware.

## How to Run

1. Open the `BSM_Search_Simulator_Ivan_Cagnani_2026.m` file in MATLAB.
2. Click the **Run** button.
3. The SSFM physics engine will compute the true phases (expect 30 to 120 seconds depending on the processor).
4. The Monte Carlo engine will execute the 25,000 cycles (expect 3 to 30 seconds depending on the processor).
5. The final sensitivity calculations will output to the command window, and the 3-panel dashboard will render.

### Dependencies

This code relies on native MATLAB FFT optimizations and requires a full desktop version of MATLAB (R2016b or newer). **It is not compatible with GNU Octave.**

It requires the following toolbox:
* **Statistics and Machine Learning Toolbox™** (specifically for the `fitdist` and `pdf` functions used to render the Gaussian fit on the final histogram).
