# 3D Heliocentric Low-Thrust Interplanetary Transfer Optimization (ICLOCS2 + IPOPT)

[![MATLAB](https://img.shields.io/badge/MATLAB-R2021b%2B-blue.svg)](https://www.mathworks.com/products/matlab.html)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](#author--license)
[![Optimization: IPOPT](https://img.shields.io/badge/Solver-IPOPT%20%2F%20ICLOCS2-orange.svg)](https://imperialcollegelondon.github.io/ICLOCS/)

## Overview

This repository contains an advanced MATLAB-based framework for solving **3D Heliocentric Low-Thrust Interplanetary Transfers** with simultaneous optimization of the **departure epoch**, **arrival epoch**, and **time-of-flight (TOF)**. 

This work represents an extension and refactoring of the direct collocation framework originally developed for launcher ascent trajectory optimization (using ICLOCS2 and IPOPT). The formulation is adapted to handle continuous low-thrust propulsion, dual-boundary ephemeris matching via normalized polynomials, and singularity-free local Cartesian steering controls in interplanetary space.

## Key Features

- **Direct Collocation Architecture:** Formulates the interplanetary trajectory optimization as a Non-Linear Programming (NLP) problem using trapezoidal transcription, solved efficiently via IPOPT.
- **Simultaneous Launch & Arrival Window Optimization:** Optimizes both the time of flight ($t_f$) and the arrival epoch shift ($p_1 = \Delta t_f$) within prescribed bounds, automatically determining the optimal departure epoch ($\Delta t_0$) and total revolutions ($n_{\text{rev}}$).
- **Smooth $C^\infty$ Ephemeris Parser:** Automatically parses JPL Horizons vector tables (`.txt`/`.csv`) for Earth and target bodies, evaluating continuous boundary conditions at both departure and arrival through centered and scaled 4th-degree polynomial approximations (Horner's scheme) in $[-1, 1]$.
- **Singularity-Free Cartesian Control Formulation:** Employs a 3D heliocentric spherical state representation $[r, \theta, \phi, v, \gamma, \chi, m]$ paired with a local Cartesian unit thrust direction vector $[u_x, u_y, u_z]$ (Up, North, East) subject to $u_x^2 + u_y^2 + u_z^2 - 1 = 0$ and engine throttling $\tau \in [0, 1]$, ensuring a smooth and well-conditioned Hessian matrix.
- **Multi-Grid Homotopy Continuation:** Implements a robust continuation strategy starting from an energy-optimal solution on a fast 150-node grid and transitioning to a fuel-optimal bang-bang profile on a dense 300-node grid with warm-started multipliers.
- **Rigorous Verification Pipeline:** Includes open-loop verification via MATLAB's `ode45` integrator with strict tolerances ($\text{RelTol} = \text{AbsTol} = 10^{-10}$), ensuring full consistency between collocation nodes and continuous dynamics.
- **Advanced Post-Processing:** Automated reconstruction of pitch and yaw pointing angles ($\alpha, \beta - \chi$), extraction of Cartesian velocity/position costates (Primer Vector theory), switching functions, Hamiltonian invariance verification, and 3D HECI trajectory plots.

## Repository Structure

- `main.m` — Main orchestration script configuring heliocentric canonical units, spacecraft propulsion parameters, running the two-step energy-to-fuel homotopy solver, verifying via `ode45` simulation, and generating kinematic/costate/3D plots.
- `low_thrust.m` — Optimal Control Problem (OCP) definition file containing state/input bounds, 12-equation boundary conditions (`b_unscaled`), unit-norm path constraints, and scaled/unscaled cost functions for ICLOCS2.
- `low_thrust_dynamics.m` — 3D heliocentric inertial equations of motion in spherical coordinates with local Cartesian thrust projection.
- `get_boundary_conditions.m` — Ephemeris parser, automatic multi-revolution estimator, and smooth 4th-degree polynomial generator for Earth departure and target arrival windows.
- `settings_low_thrust.m` — Configuration file for ICLOCS2 (trapezoidal transcription) and tuned IPOPT solver settings.
- `ephem_earth.txt` & `ephem_mars.txt` — JPL Horizons heliocentric state vector tables for Earth and Mars (target).

## Prerequisites & Dependencies

To run this code, you need:
1. **MATLAB** (R2021b or newer recommended).
2. **ICLOCS2 Toolbox** (Imperial College London Optimal Control Software v2). Ensure the `src` and `usr` directories are properly added to your MATLAB path.
3. **IPOPT** NLP solver (interfaced via MEX files within ICLOCS2).
4. **ADiGator** *(optional)* for automatic differentiation (numeric finite-difference derivatives are configured by default in `settings_low_thrust.m`).
5. **JPL Horizons Ephemeris Files** (`ephem_earth.txt` and `ephem_mars.txt` placed in the root directory or inside an `ephemeris/` subfolder).

## Usage

1. Clone or download this repository and open it in MATLAB.
2. *(Optional)* Replace `ephem_earth.txt` and `ephem_mars.txt` or adjust `tf_epoch`, `arr_window_days`, and `TOF_window_days` in `main.m` to explore different missions or launch windows.
3. Update the `addpath` lines at the top of `main.m` to point to your local installation of `ICLOCS2`.
4. Run the main script:
   ```matlab
   run('main.m');
   ```

---

## Academic Context & Acknowledgments

This code is developed solely for academic research, advanced astrodynamics studies, and educational purposes at **Alma Mater Studiorum - Università di Bologna**.

This framework interfaces with third-party open-source optimization tools, which remain the property of their respective authors under their own licensing terms:
- **ICLOCS2:** Y. Nie, O. Faqir, and E. C. Kerrigan (2018), *"ICLOCS2: Try this Optimal Control Problem Solver Before you Try the Rest"*, Imperial College London (MIT License).
- **IPOPT:** A. Wächter and L. T. Biegler (2006), *"On the implementation of an interior-point filter line-search algorithm for large-scale nonlinear programming"*, Mathematical Programming (EPL License).
- **Ephemeris Data:** Provided by the NASA/JPL Solar System Dynamics group via the Horizons System.

---

## Author & License

- **Author:** Gabriele Borghini (2026)
- **Copyright:** © 2026 Gabriele Borghini. All Rights Reserved.
- **License (MIT):** Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files, to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, and sublicense copies of the Software, subject to the inclusion of the original copyright and permission notice in all copies or substantial portions of the Software.
- **Disclaimer:** This software is provided "as is", without warranty of any kind, express or implied, including but not limited to the warranties of accuracy, numerical stability, or fitness for a particular purpose.
