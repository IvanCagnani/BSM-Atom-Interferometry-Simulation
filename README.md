# Simulation Code: A Background-Free Search for Physics Beyond the Standard Model Using Atom Interferometry
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18896187.svg)](https://doi.org/10.5281/zenodo.18896187)

This repository contains the advanced MATLAB simulation suite used to generate the data, statistical proofs, and figures for the paper:

**A Background-Free Search for Physics Beyond the Standard Model Using Atom Interferometry**
*Ivan Cagnani (Linnaeus University)*

## Overview

Extracting a μrad-scale Beyond Standard Model (BSM) phase shift from an atom interferometer requires overcoming gigaradian-scale gravitational backgrounds, massive vibration noise, and the fundamental kinematic degeneracy between scalar BSM couplings and the spatial gradient of the quadratic Stark shift.

This code rigorously validates the experimental protocol proposed in the paper. To avoid the superficiality of simplified algebraic models and the prohibitive computational bottleneck of massive matrix exponentiation, this simulator is built on a high-speed, two-stage architecture: an **ab initio physics engine** followed by a **statistical noise engine**.

## Methodology

The script executes two primary computational tasks:

### Part 1: The Fast SSFM Physics Engine
Because the $k$-odd parity of the BSM and Stark operators is strictly dependent on the physical spatial separation of the atomic wavepackets, it must be evaluated quantum mechanically. 
* The engine natively solves the 1D time-dependent Schrödinger equation using the **Split-Step Fourier Method (SSFM)**.
* To prevent catastrophic momentum aliasing from the macroscopic gravitational acceleration ($g = 9.81$ m/s²), the simulation utilizes a massive spatial grid of **$N = 1,048,576$ ($2^{20}$) points**. 
* By alternating between position and momentum space via Fast Fourier Transforms (FFTs), it explicitly resolves the physical spatial separation and establishes the true theoretical baseline differential phases.

### Part 2: Statistical Monte Carlo Noise Engine
Using the true quantum phases established by the SSFM, the Monte Carlo engine simulates the signal extraction under severe experimental noise conditions over 25,000 complete 4-point measurement cycles (representing a 17-day integration time).
* **Mid-Fringe Spin-Squeezing:** It injects massive $2\pi$ common-mode vibration noise onto $N=10^6$ optimally spin-squeezed atoms. It utilizes active phase sweeping and a $\pi/2$ optical offset to geometrically open the covariance trace into a perfect circle, mathematically defeating *fringe-peak rectification bias*.
* **Lock-In Quadrature Validation:** It injects a severe macroscopic low-frequency differential drift (~357 μrad). The simulation proves that the Double-Difference lock-in successfully algebraically annihilates the gigaradian gravity background and rejects the wandering drift to isolate the target parity signal at the fundamental quantum limit.

## Output Visualizations

Running the script automatically generates **6 publication-ready figures** and an exact console printout for LaTeX manuscript integration:

1. **Fig 1: Spatial Separation:** Visualizes the probability density at $t=T$, proving the fundamental $k$-odd parity of the spatial potentials.
2. **Fig 2: Mid-Fringe Covariance:** Plots the active phase sweep of the simultaneous dual-isotope sequence, demonstrating unbiased covariance extraction in high-noise regimes.
3. **Fig 3: Lock-In Extraction:** A histogram of the 25,000 lock-in outputs, proving the successful isolation of the target nanoradian signal from the annihilated gravitational background.
4. **Fig 4: Sigma Distribution:** A histogram of the 120 independent Monte Carlo true-seed extractions, confirming a mean statistical significance of 12.6$\sigma$.
5. **Fig 5: Integration Convergence:** Plots the simulated phase error against the ideal $1/\sqrt{N}$ Standard Quantum Limit (SQL), proving metrological stability over 25,000 cycles.
6. **Fig 6: Q-Q Plot:** A Quantile-Quantile plot visually validating the strict normal distribution of the extracted phase.

## Code Files

This repository contains two versions of the simulation suite:

**1. `BSM_Search_Simulator_Ivan_Cagnani_2026.m` (Standard Version)**
* **The primary, fast-running script.**
* Executes a single 25,000-cycle Monte Carlo run using MATLAB's default pseudorandom number generator.
* Ideal for quick verification and generating the core physical visualizations. Highly optimized, it completes the massive quantum simulation and statistical cycles in under two minutes on standard hardware.

**2. `BSM_Search_Simulator_Ivan_Cagnani_2026_with_Statistical_Validation.m` (True-Seed Validation Version)**
* **The rigorous statistical proof script.**
* Wraps the Monte Carlo engine in a 120-run master loop, explicitly initializing MATLAB's Mersenne Twister with 120 cryptographically secure, true random integers derived from atmospheric noise (via random.org).
* Mathematically proves the isolated signal (12.6$\sigma$ mean confidence) is structurally sound and completely free from algorithmic PRNG artifacting using native Anderson-Darling and Lilliefors normality tests.
* *Note:* Because it executes 3,000,000 total statistical cycles (120 runs × 25,000 cycles), expect a longer execution time (approx. 15 to 30 minutes depending on processor speed).

## How to Run

1. Open either `.m` file in MATLAB.
2. Click the **Run** button.
3. The SSFM physics engine will compute the true phases (expect 30 to 120 seconds).
4. **If running the Standard version:** The Monte Carlo engine will execute the 25,000 cycles (~5 to 30 seconds).
5. **If running the Validation version:** The script will iterate through the 120 true-seed runs and output the Anderson-Darling and Lilliefors normality results directly to the command window.
6. The final sensitivity calculations will print to the console (formatted for LaTeX), and the 6 figures will render.

### Dependencies

This code relies on native MATLAB FFT optimizations and requires a full desktop version of MATLAB (R2020a or newer). **It is not compatible with GNU Octave**.

It requires the following toolbox:
* **Statistics and Machine Learning Toolbox™**
  * Used for the `fitdist` and `pdf` functions to render the Gaussian fits.
  * Used for the `adtest` (Anderson-Darling) and `lillietest` (Lilliefors) functions in the Validation version.
