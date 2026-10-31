# Thermal Stability of Wigner Crystals in 2D Semiconductors

## Project Overview

Developed a Python-based stochastic simulation to investigate the thermal stability of two-dimensional Wigner crystals, where electrons form ordered lattice structures due to Coulomb repulsion. The project studies how temperature-driven particle fluctuations affect crystal stability and potentially contribute to lattice breakdown.

## Theory

The model describes interacting electrons confined by a harmonic potential, subject to Coulomb repulsion and stochastic thermal perturbations. Particle dynamics follow a Langevin-inspired equation of motion:

$$m\frac{d^2\mathbf{r}_i}{dt^2}=-\mathbf{r}_i+\sum_{j\neq i}\frac{\mathbf{r}_i-\mathbf{r}_j}{|\mathbf{r}_i-\mathbf{r}_j|^3}+\boldsymbol{\eta}_i(t)$$

where the stochastic term represents Gaussian thermal fluctuations controlled by temperature \(T\).

Crystal stability is assessed through particle displacement trajectories and positional variance, with larger fluctuations indicating reduced structural stability.

## How It Was Built

Developed a Python stochastic simulation of 26 interacting electrons over 10,000 time steps using finite-difference integration of coupled equations of motion. Modeled nonlinear Coulomb interactions and harmonic confinement using NumPy, incorporating Gaussian-distributed thermal perturbations into particle dynamics. Implemented modular indexing to calculate pairwise electrostatic forces while avoiding self-interactions. Applied time-series analysis and variance estimation to individual electron trajectories, examining positional fluctuations and the average variance of selected inner-shell particles to assess lattice stability.

## Results

Simulated temperature-driven fluctuations in a 26-electron Wigner crystal, tracking particle positions, velocities, and accelerations across 10,000 time steps. Analyzed positional variance and particle trajectories to quantify thermal perturbations and investigate lattice stability under stochastic forcing. Generated time-series and two-dimensional trajectory visualizations to characterize electron displacement from initial lattice configurations.

## Technologies

Python, NumPy, Matplotlib, Stochastic Differential Equations, Langevin Dynamics, Gaussian Noise, Finite-Difference Methods, Monte Carlo Simulation, Time-Series Analysis, Variance Analysis, Computational Statistical Mechanics.
