# Chaos-Quantized PRNG: Lorenz System Implementation in Julia

[![Language](https://img.shields.io/badge/language-Julia-9558b2.svg)](https://julialang.org/)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

This repository contains an implementation of a Chaos-based Pseudo-Random Number Generator (PRNG). It uses the **3D Lorenz System** as an entropy source and applies a sophisticated **Quantization Pipeline** to transform deterministic chaotic flows into uncorrelated bitstreams.

The theoretical framework for this implementation is detailed in the accompanying paper: *[Quantization Techniques in Chaos-Based PRNGs](./docs/Quantization_Techniques_in_Pseudo_Random_Numbers_Generators.pdf)*.

---

## 🔬 Scientific Background

Chaotic systems are inherently continuous and deterministic. To be used in digital cryptography, they must undergo **Quantization**. This project implements a robust three-stage pipeline to address the core challenges of chaotic PRNGs:

1.  **The Correlation Problem:** Breaking the continuity of the Lorenz flow.
2.  **Dynamic Degradation:** Preventing the "collapse" into limit cycles caused by finite 64-bit float precision.
3.  **The Bias Trap:** Ensuring a perfectly uniform 50/50 distribution of bits.



---

## 🛠 Quantization Pipeline

The generator processes the Lorenz state $(x, y, z)$ through the following mathematical chain:

* **Magnification:** States are multiplied by $10^{14}$ to shift high-entropy microscopic fluctuations into the integer domain.
* **LSB Extraction:** The 8 least significant bits are extracted via bitwise masking (`& 0xFF`).
* **XOR Mixing:** The bits from all three dimensions are combined ($Q = x_{lsb} \oplus y_{lsb} \oplus z_{lsb}$) to eliminate structural symmetry and maximize diffusion.


---

## ⚡ Why Julia?

This project was implemented in **Julia** to achieve the "best of both worlds":
* **Performance:** Execution speed matching C/C++ for the heavy numerical integration (RK4).
* **Ergonomics:** Syntax as clean and readable as Python.
* **Bitwise Efficiency:** Highly efficient handling of large-integer scaling and bitwise XOR operations.

---

## 📂 Repository Structure

```text
.
├── docs/
│   ├── Quantization_Techniques_in_Pseudo_Random_Numbers_Generators.pdf # The paper
├── output/
│   └── plot1_waveform.png
|   |-- plot2_phase_raster.png
|   |-- plot3_acf.png
|   |-- plot4_returnmap.png
├── scripts/
│   ├── main.jl    # Core simulation code
│   ├── required-packages.jl       # Loads the required packages
└── README.md