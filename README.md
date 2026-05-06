# Riemann‑ID

**Validation of the Real‑Imaginary Duality Principle: from the Riemann Hypothesis to prime density generation**

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20054082.svg)](https://doi.org/10.5281/zenodo.20054082)

This repository contains the complete source code and generated figures for the paper *“An Information‑Dynamics Proof of the Non‑Trivial Zeros of the Riemann Zeta Function: From Dynamic Real‑Imaginary Coupling to Prime Density Generation”* (Kai Huang, 2026).

## What is inside

| File | Description |
|------|-------------|
| `prime_density.py` | Numerical integration of the generalized Ginzburg–Landau (GL) equation driven by the coupling matrix $K(x)=\sqrt{2x\ln x}$. Evolves the complex information field from random noise and outputs **fig1_prime_density.pdf**. |
| `prime_density_binned_validation.py` | Binned $\pi(x)$ comparison. Divides [2, 100] into 10 sub‑intervals and compares the model‑predicted prime count with the exact count (using `sympy.primepi`). Outputs `output_prime_density_binned_validation.txt`. |
| `GUE_spacings.py` | Generates a 1500×1500 Gaussian Unitary Ensemble (GUE) random matrix, computes normalized nearest‑neighbour spacings, and compares them with the first 100 Riemann zeros (Odlyzko data). Reports the Kolmogorov–Smirnov statistic. Outputs **fig3_gue_spacings.pdf**. |
| `coupling_matrix.py` | Plots the dynamic coupling matrix $K(x)$ together with the prime density $1/\ln x$. Outputs **fig2_coupling_matrix.pdf**. |
| `zero_mapping.py` | Illustrates the schematic mapping from eigenvalues $\lambda_n^2$ of the squared operator $\Delta_K$ to critical‑line zeros via $s(1-s)=\lambda_n^2$. Outputs **fig4_zero_mapping.pdf**. |
| `fig1_prime_density.pdf` | Self‑organised steady‑state density vs. target $1/\ln x$ and relative error (log scale). |
| `fig2_coupling_matrix.pdf` | Coupling matrix $K(x)$ and prime density $1/\ln x$. |
| `fig3_gue_spacings.pdf` | Nearest‑neighbour spacing distribution: GUE model vs. first 100 Riemann zeros. |
| `fig4_zero_mapping.pdf` | Schematic of the eigenvalue‑to‑zero mapping on the critical line. |

## Reproducing the results

### Requirements
- Python ≥ 3.8
- NumPy, SciPy, Matplotlib, SymPy

Install dependencies:
```bash
pip install numpy scipy matplotlib sympy
```

### Generate the figures
```bash
python prime_density.py
python prime_density_binned_validation.py
python GUE_spacings.py
python coupling_matrix.py
python zero_mapping.py
```

All PDF figures will be written to the current working directory.

## Citation
If you use this code or data in your research, please cite:

**Kai Huang.** *An Information‑Dynamics Proof of the Non‑Trivial Zeros of the Riemann Zeta Function: From Dynamic Real‑Imaginary Coupling to Prime Density Generation*. Zenodo, 2026.  
DOI: [10.5281/zenodo.20054082](https://doi.org/10.5281/zenodo.20054082)

```bibtex
@misc{huang2026riemann,
  author       = {Kai Huang},
  title        = {An Information‑Dynamics Proof of the Non‑Trivial Zeros of the 
                  Riemann Zeta Function: From Dynamic Real‑Imaginary Coupling 
                  to Prime Density Generation},
  year         = 2026,
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.20054082},
}
```

## Contact
Kai Huang – hkaiopen@foxmail.com  
Project repository: [https://github.com/hkaiopen/Riemann-ID](https://github.com/hkaiopen/Riemann-ID)
