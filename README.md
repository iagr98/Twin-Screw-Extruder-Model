# MATLAB Model — Reactive Extrusion in a Corotating Twin‑Screw Extruder

This repository contains a MATLAB implementation of a **reactive extrusion** model for the ring‑opening polymerization of lactide (PLA) in a **corotating, tightly‑meshing twin‑screw extruder**. The model couples a reactor‑network representation of the extruder with reaction kinetics, pressure/throughput calculations, and property estimation (e.g., viscosity links, number‑average molar mass). It is designed for fast simulation and exploratory studies, not for turnkey plant operation.

> The model structure and equations summarized here are based on the accompanying overview document included with this project. Key ideas: overall coupling of apparatus and reaction model, pressure build‑up via restrictive elements, and a simplified kinetic scheme with catalyst activation + propagation. 

---

## Contents

- `src/` — MATLAB functions and classes implementing the model
- `results/` — results of all ran simulations
- `overview.pdf` — document explaining how the model is built
- `README.md` — this file
- `LICENSE` — MIT


---

## What the model does

### 1) Apparatus model (extruder as a reactor network)
The corotating twin‑screw extruder is represented as a sequence of **ideal reactors** (CSTR‑cascade). Each element (screw (RSE/LSE), kneading block (KB), Die) contributes:
- **Drag flow** and **pressure‑driven flow**, with restrictions in **non‑ or counter‑conveying** elements or due to the **die**. Only **fully filled** elements can build pressure.
- Global pressure field solved from a **linear system** `A·p = B`, assembled from element‑to‑element pressure/drag couplings; solved in MATLAB (e.g., `linsolve`). 

### 2) Reaction model (simplified kinetics)
A reduced kinetic scheme is used:
- **Catalyst activation / deactivation** between species `C, A` and initiator/chain species `ROH, I, R, D` with forward/backward rates `ka1, ka2`.
- **Propagation** between active chains `R` and monomer `M` with `kp`, reversible deactivation `kd`. 

For a (local) batch‑type description within each reactor, the ODE system is:  
`dM/dt = −kp·M·(I+R) + kd·R`, …, up to `dD/dt = −ka1·C·D + ka2·A·R` 

### 3) Viscosity model
This model uses a viscosity model which depends on Monomer concentration, shear rate within the extruder and temperature in the substance.

### 4) Simulation run nature
- First a simulation (**flow**) solving mass balance and energy (Temperature) balance is ran while assuming constant monomer concentration for viscosity calculation. Level and temperature profile along the extruder is determined.
- Simulation is coupled with a Reaction model simulation for determining species concentration along extruder.
- Viscosity is updated with respect to Monomer concentration, temperature and shear rate.
- A second **flow** simulaton is ran for determination of final state of extruder.

---

## Getting started

### Requirements
- MATLAB R2021b or newer (tested); base toolboxes (no Simulink required for the core model).

### Installation and Running of the Model
1. Clone the repository  
   ```bash
   git clone https://github.com/iagr98/Twin-Screw-Extruder-Model.git   
   ```
2. Add `src/` to the MATLAB path

3. Run `RTD_simulation.m`

4. Simulation's Workspace saved in folder `results`
---


## Citation

If you use this code or the model structure in academic work, please cite the accompanying overview and this repository.

```
@misc{reactive_extrusion_matlab,
  title  = {MATLAB Model — Reactive Extrusion in a Corotating Twin‑Screw Extruder},
  author = {Niclas Conen, Filip Latz},
  year   = {2025},
  url    = {https://github.com/iagr98/Twin-Screw-Extruder-Model.git}
}
```

---

## License

MIT
