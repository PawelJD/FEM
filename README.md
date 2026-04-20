# 🌡️ 2D Heat Transfer Simulation using Finite Element Method

A C++ implementation of a **2D transient heat transfer simulation** using the Finite Element Method (FEM). The program calculates how temperature distributes and evolves over time across a 2D mesh, using 4-node quadrilateral elements (Q4).

---

## Overview

This simulator solves the **transient heat conduction problem** — given a 2D mesh with defined material properties and boundary conditions, it computes the temperature at every node at each time step. The governing equation is discretized using FEM and integrated over time using the implicit **Backward Euler** method, which ensures numerical stability.

---

## Theory

The simulation is grounded in the following concepts:

- **Heat conduction equation** (Fourier's law) discretized over a 2D domain
- **Q4 shape functions** - bilinear interpolation over 4-node quadrilateral elements
- **Numerical integration** via Gauss quadrature (2×2 or higher)
- **FEM assembly** - local element matrices are assembled into global system matrices
- **Time integration** - implicit Backward Euler scheme: `[H + C/Δt]{T}ⁿ⁺¹ = {P} + [C/Δt]{T}ⁿ`

---

## ✨ Features

- 2D transient (time-dependent) heat conduction
- 4-node quadrilateral elements (Q4)
- Gaussian numerical integration
- Computation of all key FEM matrices:
  - **H** - conductivity matrix
  - **C** - heat capacity matrix
  - **Hbc** - boundary convection matrix
  - **P** - load vector
- Implicit time stepping (Backward Euler)
- Global matrix assembly from local element contributions
- System of equations solved via Gaussian elimination
- Outputs minimum and maximum temperature at each time step

---

## 📁 Project Structure

```
FEM/
├── Final/
│   ├── Main.cpp               # Entry point, simulation loop
│   ├── GaussIntegration.cpp   # Gauss quadrature points and weights
│   └── ...                    # Element matrix computations, assembly
├── Mesh4x4.txt                # Sample mesh — 16 nodes, 9 elements (uniform)
├── Mesh4x4_Mixgrid.txt        # Sample mesh — 16 nodes, 9 elements (non-uniform)
├── Mesh31x31.txt              # Sample mesh — 961 nodes, 900 elements (large)
└── README.md
```

---

## Input Format

The simulation reads from a plain text input file. Below is an example of a valid input:

```
SimulationTime      500
SimulationStepTime  50
Conductivity        25
Alfa                300
Tot                 1200
InitialTemp         100
Density             7800
SpecificHeat        700
Nodes               16
Elements            9
*Node
  1, 0.0, 0.0
  2, 0.1, 0.0
  ...
*Element
  1, 1, 2, 5, 4
  ...
```

| Parameter            | Description                           |
|----------------------|---------------------------------------|
| `SimulationTime`     | Total simulation duration [s]         |
| `SimulationStepTime` | Time step size Δt [s]                 |
| `Conductivity`       | Thermal conductivity k [W/(m·K)]      |
| `Alfa`               | Convection coefficient α [W/(m²·K)]   |
| `Tot`                | Ambient temperature T∞ [K]            |
| `InitialTemp`        | Initial temperature of all nodes [K]  |
| `Density`            | Material density ρ [kg/m³]            |
| `SpecificHeat`       | Specific heat capacity c [J/(kg·K)]   |

---

## 📐 Sample Meshes

Three ready-to-use mesh files are included for testing:

| File                  | Nodes | Elements | Description                               |
|-----------------------|-------|----------|-------------------------------------------|
| `Mesh4x4.txt`         | 16    | 9        | Uniform square grid                       |
| `Mesh4x4_Mixgrid.txt` | 16    | 9        | Non-uniform (irregular) node distribution |
| `Mesh31x31.txt`       | 961   | 900      | Large uniform grid for performance testing|

---

## ⚙️ Algorithm

```
1. Load mesh file (nodes, elements, simulation parameters)
2. Initialize global matrices H, C, P (zeroed)
3. For each element:
   a. Compute Jacobian at each Gauss point
   b. Calculate shape function derivatives (dN/dx, dN/dy)
   c. Compute local matrices: H, C, Hbc, P
4. Assemble local matrices into global system
5. Apply boundary conditions
6. Time loop:
   a. Build system: [H + C/Δt]{T}ⁿ⁺¹ = {P} + [C/Δt]{T}ⁿ
   b. Solve with Gaussian elimination
   c. Output min/max temperature for current time step
```

---
## 📊 Output

At each time step the program prints the minimum and maximum nodal temperature:

```
Time:  50s  |  Tmin =  100.00 K  |  Tmax =  312.45 K
Time: 100s  |  Tmin =  100.00 K  |  Tmax =  489.12 K
...
