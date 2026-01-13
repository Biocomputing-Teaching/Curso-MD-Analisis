---
layout: default
title: Cookbook
permalink: /cookbook/
---

This chapter complements the classic episodes with short recipes that explain the latest OpenMM Cookbook tutorials. Each section summarizes the essential content and links to the original notebook so students can go deeper from an English starting point.

## 1. Load files and report results
The **Loading Input Files and Reporting Results** tutorial explains how to use `OpenMM` to read inputs generated with AmberTools (prmtop+inpcrd) and, similarly, with CHARMM, GROMACS, or Tinker files. It shows how to easily log energy, pressure, and coordinates using reporters built in the application layer, and how to define custom reporters that trigger every integration step. [Open in the Cookbook](https://openmm.github.io/openmm-cookbook/latest/notebooks/tutorials/loading_and_reporting.html)

## 2. Build systems from scratch
In **Building Systems from Scratch**, OpenMM Setup is set aside to work directly with the library layer. A binary Lennard-Jones system is created, mixing rules are defined (Lorentz-Berthelot), and the base forces (`Nonbonded`, `HarmonicBond`, etc.) are built manually. This approach is ideal for experimenting with simple models or inspecting `System` objects before using them in biomolecular simulations. [Open in the Cookbook](https://openmm.github.io/openmm-cookbook/latest/notebooks/tutorials/building_systems.html)

## 3. Simulation parameters
The **Selecting Values for Simulation Parameters** tutorial dives into the trade-offs between speed and accuracy. It details how to choose the integration step size based on the fastest relevant mode, how to work with constraints (`HBonds`, `AllBonds`, `HAngles`), how to apply hydrogen mass repartitioning, and which thermostat/barostat settings suit different ensembles. It is a roadmap for tuning numerical accuracy even when working with previously generated scripts. [Open in the Cookbook](https://openmm.github.io/openmm-cookbook/latest/notebooks/tutorials/simulation_parameters.html)

## 4. Alchemical free energies
**Alchemical Free Energy Calculations** shows how to use `CustomForce` and global parameters to represent intermediate states between two systems, and how to combine it with `OpenMMTools` and `PyMBAR` to estimate free energy differences. It includes recommendations for installing dependencies (`openmmtools`, `pymbar`) and links to the reference `YANK` code. This is ideal for extending the protein-ligand focused episodes. [Open in the Cookbook](https://openmm.github.io/openmm-cookbook/latest/notebooks/tutorials/Alchemical_free_energy_calculations.html)

## 5. Coarse-grained polymers
The **Implementing a Coarse-Grained Polymer Force Field** chapter teaches how to build polymer topologies with Lennard-Jones beads, assign mass and realistic physical parameters, and define the force both via API and XML files. It combines redox modeling techniques and manual parameterization to show how protein simulations can be adapted to coarse-grained frameworks. [Open in the Cookbook](https://openmm.github.io/openmm-cookbook/latest/notebooks/tutorials/coarse_grained_polymer.html)

## 6. Replica exchange with solute tempering (REST)
**Running a REST simulation** details the creation of a REST system: copy a vanilla `System`, define REST atoms, add `CustomBondForce`, `CustomAngleForce`, and `CustomTorsionForce` with global factors, and scale the `NonbondedForce` for modified interactions. It explains running replicas and coordinating with `OpenMMTools`. This material fits perfectly after the MSM/slow processes episode. [Open in the Cookbook](https://openmm.github.io/openmm-cookbook/latest/notebooks/tutorials/Running_a_REST_simulation.html)

## 7. Umbrella sampling
Finally, **Umbrella Sampling** summarizes the theoretical foundations (CVs, bias potentials, and free energy profiles), and explains how to apply Gaussian biases to explore high-energy regions. It includes the relation \(F(x)=-k_B T \log(p(x))\) and how to reconstruct the profile by combining histograms. [Open in the Cookbook](https://openmm.github.io/openmm-cookbook/latest/notebooks/tutorials/umbrella_sampling.html)
