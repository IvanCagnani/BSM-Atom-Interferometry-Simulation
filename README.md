# **Monte Carlo Simulation for "A Background-Free Search for Physics Beyond the Standard Model Using Atom Interferometry"**

This repository contains the MATLAB simulation code used to generate the figures and sensitivity projections in the paper:

A Background-Free Search for Physics Beyond the Standard Model Using Atom Interferometry  
Ivan Cagnani (Linnaeus University)  
\[Link to your paper (e.g., on arXiv) will go here\]

## **Code Files**

This simulation is provided in two formats. The .m file is the primary, recommended version.

1. **BSM\_Monte\_Carlo\_Simulator.m**  
   * **This is the primary, public-facing script for the project.**  
   * It is a plain-text MATLAB file, fully commented, and can be viewed directly in your browser with syntax highlighting.  
   * This file will run on most versions of MATLAB (R2016b or newer).  
2. **BSM\_Monte\_Carlo\_Simulator.mlx**  
   * This is the original MATLAB Live Script *development* file. It is not formatted for presentation.  
   * For the cleanest, most readable version of the code, please use the .m file.

## **Methodology**

This script validates the paper's 4-point differential quadrature technique. It performs a Monte Carlo simulation to test the experiment's ability to isolate a target BSM signal from both quantum noise (QPN) and large, injected systematic errors.

The code performs two main tasks:

1. **Single-Shot Noise Analysis:** It runs 1,000,000 trials of a single 4-point measurement to determine the 1-sigma noise floor (sigma\_1). This generates the histograms in Figure 1\.  
2. **Sensitivity Projection:** It scales this sigma\_1 noise floor by 1/sqrt(N) to project the total integration time required for a 5-sigma discovery. This generates the plot in Figure 2\.

## **How to Run**

1. Open the BSM\_Monte\_Carlo\_Simulator.m file in MATLAB.  
2. Click the "Run" button.  
3. The script will execute the Monte Carlo simulation (which may take a few seconds) and then output the final sensitivity calculations to the command window and generate the two figures from the paper.

### **Dependencies**

This code is **not** compatible with GNU Octave. It requires a full version of MATLAB and the following two toolboxes:

* **Parallel Computing Toolbox™** (for the parfor loop)  
* **Statistics and Machine Learning Toolbox™** (for the fitdist function used in plotting)
