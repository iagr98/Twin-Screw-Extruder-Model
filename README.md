# MATLAB Model — Reactive Extrusion in a Corotating Twin‑Screw Extruder

This repository contains a MATLAB implementation of a **reactive extrusion** model for the ring‑opening polymerization of lactide (PLA) in a **corotating, tightly‑meshing twin‑screw extruder**. The model couples a reactor‑network representation of the extruder with reaction kinetics, pressure/throughput calculations, and property estimation (e.g., viscosity links, number‑average molar mass). It is designed for fast simulation and exploratory studies, not for turnkey plant operation.

> The model structure and equations summarized here are based on the accompanying overview document included with this project. Key ideas: overall coupling of apparatus and reaction model, pressure build‑up via restrictive elements, and a simplified kinetic scheme with catalyst activation + propagation. fileciteturn0file2

---

## Contents

- `src/` — MATLAB functions and classes implementing the model
- `results/` — results of all ran simulations
- `overview.pdf` — document explaining how the model is built
- `README.md` — this file
- `LICENSE` — MIT

> You may rename or restructure folders; the layout above is a suggested starting point.

---

## What the model does

### 1) Apparatus model (extruder as a reactor network)
The corotating twin‑screw extruder is represented as a sequence of **ideal reactors** (CSTR‑like cells). Each element (screw, kneading block, left/right conveying) contributes:
- **Drag flow** and **pressure‑driven flow**, with restrictions in **non‑ or counter‑conveying** elements or due to the **die**. Only **fully filled** elements can build pressure. fileciteturn0file6
- Global pressure field solved from a **linear system** `A·p = B`, assembled from element‑to‑element pressure/drag couplings; solved in MATLAB (e.g., `linsolve`). fileciteturn0file7 fileciteturn0file3

### 2) Reaction model (simplified kinetics)
A reduced kinetic scheme is used:
- **Catalyst activation / deactivation** between species `C, A` and initiator/chain species `ROH, I, R, D` with forward/backward rates `ka1, ka2`.
- **Propagation** between active chains `R` and monomer `M` with `kp`, reversible deactivation `kd`. fileciteturn0file0

For a (local) batch‑type description within each reactor, the ODE system is:  
`dM/dt = −kp·M·(I+R) + kd·R`, …, up to `dD/dt = −ka1·C·D + ka2·A·R` (Eqs. 10–16). fileciteturn0file0

### 3) Mass balances in each reactor
Each cell uses a component mass balance with **flow** + **reaction** terms. The general CSTR form (Eq. 17) is implemented; **isochoric** reaction implies `dV/dt` follows the **overall mass balance** (Eq. 18). fileciteturn0file1

### 4) Property estimation (Mn)
Number‑average molar mass `Mn` is updated from reaction state. In the simplest case (no catalyst‑initiated chains) and with quasi‑instantaneous activation,  
`dMn/dt = − 2·M_rep · (dM/dt) / ROH0` (Eqs. 19–23). A two‑branch scheme avoids numerical issues when catalyst also starts chains (Eq. 24). fileciteturn0file5 fileciteturn0file4

---

## Assumptions & scope (read before using)

- **Extruder type:** corotating, tightly‑meshing twin‑screw only. Element taxonomy covers conveying screws and kneading blocks; left‑handed (upstream) or >90° kneading disks are treated as **restrictive** elements. fileciteturn0file6  
- **Filling & pressure:** only **fully filled** cells accumulate pressure; partially filled cells are at ambient pressure. Pressure‑driven flow between neighboring cells follows Eq. (3). The global system uses `A·p=B` (Eq. 4). fileciteturn0file6  
- **Reactions:** simplified scheme (activation + propagation) without chain‑length resolution. Suitable for **fast studies**; for detailed MWD predictions, extend with method‑of‑moments. fileciteturn0file0  
- **Thermal model:** energy balance hooks are included conceptually; if not provided, temperature is treated as given per zone (isothermal cells). fileciteturn0file2  
- **Rheology:** viscosity can be coupled to conversion/Mn and shear; provide a closure or correlation as needed (placeholder provided). fileciteturn0file2

---

## Getting started

### Requirements
- MATLAB R2021b or newer (tested); base toolboxes (no Simulink required for the core model).
- Optional: Optimization Toolbox (for parameter estimation).

### Installation
1. Clone the repository  
   ```bash
   git clone https://github.com/<your-org>/<your-repo>.git
   cd <your-repo>
   ```
2. Add `src/` to the MATLAB path:
   ```matlab
   addpath(genpath(fullfile(pwd,'src')));
   ```

### Quick start (minimal example)
```matlab
% 1) Load or define extruder layout and operating point
layout   = examples.layout_kneading_demo();   % cells/elements, left/right screws, KB angles
op       = examples.operating_point_demo();   % N (rpm), feed rate, die resistance, T(z)

% 2) Set kinetic parameters & initial composition
kin      = params.kinetics_default();         % ka1, ka2, kp, kd, etc.
state0   = params.initial_state_demo();       % M0, ROH0, catalyst C0, etc.

% 3) Run the simulation
results  = rextruder.simulate(layout, op, kin, state0);

% 4) Inspect outputs
plots.overview(results);                      % p(z), fill(z), M conversion, Mn(z), throughput
```

**Outputs** typically include per‑cell: pressure, fill fraction, drag/pressure flows, monomer conversion, `Mn`, and (optionally) temperature.

---

## How it works (equations implemented)

- **Filling logic & flows:** Drag flows `ṁ_D` from element geometry + screw speed; pressure flows `ṁ_P = (k_P·ρ/(dΩ·η))·Δp`. **Partially filled:** only `ṁ_D`; **fully filled:** mass balance uses both and pressure unknowns are solved globally. (Eqs. 1–4). fileciteturn0file6  
- **Pressure system:** For internal cells `i` the balance yields a tri‑diagonal contribution to `A`; boundaries pinned to ambient (`p_U`). Solved with `linsolve`. (Eqs. 5–9). fileciteturn0file8  
- **Reactor balances:** Component CSTR balance with reaction source from the kinetic ODEs (Eq. 17), with **isochoric** link to volume change (Eq. 18). fileciteturn0file1  
- **Mn update:** Eq. (19–24) including safeguarded branch for catalyst‑initiated chains. fileciteturn0file5 fileciteturn0file4

---

## Repository structure (suggested)

```
src/
  + rextruder/                 % main API (simulate, assemble, solvePressure)
  + apparatus/                 % element types, geometry, drag/pressure flow
  + kinetics/                  % ODE system for Eqs. (10–16)
  + balances/                  % CSTR balance (Eq. 17), volume link (Eq. 18)
  + props/                     % viscosity, Mn update (Eqs. 19–24)
  + utils/                     % helpers, linsolve wrapper, checks
examples/
  layout_kneading_demo.m
  operating_point_demo.m
  run_minimal_demo.m
params/
  kinetics_default.m
  initial_state_demo.m
plots/
  overview.m
data/
  geometry/
  temperature_profiles/
docs/
  overview.pdf
```

---

## Inputs & parameters

- **Extruder layout:** element list (type, length, pitch, staggering angle), die resistance; which cells are **restrictive** (L‑screw, ≥90° kneading). fileciteturn0file6  
- **Operating point:** screw speed `N`, feed rate, zone temperatures `T(z)`.
- **Kinetics:** `ka1, ka2, kp, kd`, initial species (`M0, ROH0, C0, A0, R0, D0, I0`). fileciteturn0file0  
- **Physical props:** density `ρ`, viscosity model (η vs. shear, T, conversion/Mn). fileciteturn0file2

---

## Extending the model

- **Energy balance:** implement heat transfer per zone, viscous dissipation, and temperature‑dependent kinetics and viscosity. fileciteturn0file2  
- **MWD / moments:** upgrade kinetics to a **method‑of‑moments** formulation if you need full distributions. (See references in the overview document.) fileciteturn0file4  
- **Non‑ideal mixing:** replace CSTR assumption with PFR/CSTR networks or RTD‑based coupling.

---

## References (from the overview)

- Jacobsen et al., 2000 — single‑step reactive extrusion of PLLA. fileciteturn0file2  
- Zapata‑González & Saldívar‑Guerra, 2023 — method of moments overview. fileciteturn0file4  
- Additional references listed at the end of the document. fileciteturn0file4

---

## Citation

If you use this code or the model structure in academic work, please cite the accompanying overview and this repository.

```
@misc{reactive_extrusion_matlab,
  title  = {MATLAB Model — Reactive Extrusion in a Corotating Twin‑Screw Extruder},
  author = {<Your Name>},
  year   = {2025},
  url    = {https://github.com/<your-org>/<your-repo>}
}
```

---

## License

Add a license file before publishing (MIT/BSD/Apache recommended).
