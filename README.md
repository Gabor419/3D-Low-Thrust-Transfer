# 3D Heliocentric Low-Thrust Interplanetary Transfer Optimization (ICLOCS2 + IPOPT)

[![MATLAB](https://img.shields.io/badge/MATLAB-R2021b%2B-blue.svg)](https://www.mathworks.com/products/matlab.html)
[![License: Proprietary](https://img.shields.io/badge/License-Proprietary-red.svg)](#author--license)
[![Optimization: IPOPT](https://img.shields.io/badge/Solver-IPOPT%20%2F%20ICLOCS2-orange.svg)](https://imperialcollegelondon.github.io/ICLOCS/)

## Overview

This repository contains an advanced MATLAB-based framework for solving **3D Heliocentric Low-Thrust Interplanetary Transfers** with a prescribed or free time-of-flight (TOF). 

This work represents a natural extension and refactoring of the direct collocation framework originally developed for launcher ascent trajectory optimization (using ICLOCS2 and IPOPT). The formulation is adapted to handle continuous low-thrust propulsion, ephemeris-driven boundary matching, and variable steering controls in interplanetary space.

## Key Features

- **Direct Collocation Architecture:** Formulates the interplanetary trajectory optimization as a Non-Linear Programming (NLP) problem, solved efficiently via IPOPT.
- **Free-TOF & Ephemeris Parser:** Automatically parses JPL Horizons vector tables (`.txt`/`.csv`) for Earth and target bodies, evaluating continuous boundary conditions through smooth $C^\infty$ local polynomial approximations (Horner's method).
- **Singularity-Free Control Formulation:** Employs an inertial spherical frame with steering angles ($\alpha, \beta_{\text{rel}}$) and engine throttling ($\tau$), ensuring a smooth and well-conditioned Hessian matrix for rapid convergence.
- **Homotopy Continuation Method:** Implements a multi-step continuation strategy transitioning from an energy-optimal smooth profile to a fuel-optimal bang-bang-like trajectory.
- **Rigorous Verification Pipeline:** Includes open-loop verification via MATLAB's `ode45` integrator with strict tolerances ($\leq 10^{-8}$), ensuring precise consistency between the collocation nodes and continuous integration.
- **Advanced Post-Processing:** Automated extraction of Cartesian costates (Primer Vector theory), switching functions, Hamiltonian evaluation, and 3D heliocentric trajectory visualizations with thrust vectors.

## Repository Structure

- `main.m` — Main orchestration script configuring heliocentric canonical units, running the energy-to-fuel homotopy solver, verifying via ODE simulation, and generating kinematic/costate plots.
- `low_thrust.m` — Optimal Control Problem (OCP) definition file containing boundary conditions, state/input bounds, cost functions, and constraints for ICLOCS2.
- `low_thrust_dynamics.m` — 3D heliocentric inertial equations of motion in spherical coordinates, optimized for ADiGator compatibility.
- `get_boundary_conditions.m` — Ephemeris parser and smooth polynomial generator for target rendezvous at a variable final time ($t_f$).
- `settings_low_thrust.m` — Configuration file for ICLOCS2 (trapezoidal transcription) and bespoke IPOPT solver settings.
- `ephem_earth.txt` & `ephem_target.txt` — Pre-extracted JPL Horizons heliocentric state vectors (2027-01-01 to 2030-12-31).

## Prerequisites & Dependencies

To run this code, you need:
1. **MATLAB** (latest version recommended).
2. **ICLOCS2 Toolbox** (Imperial College London Optimal Control Software v2). Ensure the `src` and `usr` directories are properly added to your MATLAB path.
3. **IPOPT** NLP solver (interfaced via MEX files within ICLOCS2).
4. **ADiGator** (optional, for automatic differentiation and analytical derivative generation).
5. **JPL Horizons Ephemeris Files** (`ephem_earth.txt` and `ephem_target.txt` from 2027-01-01 to 2030-12-31 are already included in the repository).

## Usage

1. Clone or download this repository and add it to your MATLAB workspace.
2. *(Optional)* Replace `ephem_earth.txt` and `ephem_target.txt` in the working directory if targeting a different celestial body or launch window.
3. Update the path configurations in `main.m` to point to your local installation of `ICLOCS2`.
4. Run the `main.m` script:
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
- **Usage Rights:** The source code, mathematical formulations, and custom scripts contained in this repository are the intellectual property of the author. They may not be copied, modified, redistributed, or used for commercial purposes without the express written consent of the author.
- **Disclaimer:** This software is provided "as is", without warranty of any kind, express or implied, including but not limited to the warranties of accuracy, numerical stability, or fitness for a particular purpose.