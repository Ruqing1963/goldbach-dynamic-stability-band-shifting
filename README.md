# Dynamic Stability of the Goldbach Locus

**Conductor Orbit Propagation and the Band Shifting Law in GSp(4)**

[![DOI](https://img.shields.io/badge/DOI-pending-blue)]()
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## Overview

This repository accompanies the paper:

> R. Chen, *Dynamic Stability of the Goldbach Locus: Conductor Orbit Propagation and the Band Shifting Law in GSp(4)*, February 2026.

We study how the conductor distribution of the Goldbach–Frey Jacobian family evolves as the even integer $N$ varies from $2^{10}$ to $2^{14}$. Despite violent discontinuities in the static conduit $\operatorname{rad}_{\mathrm{odd}}(N/2)$, the Goldbach locus exhibits **dynamic stability**: its $\rho$-band shifts rigidly according to a near-exact linear law while preserving its internal structure.

### Key Results

| Result | Value |
|--------|-------|
| **Band Shifting Law** | $\langle\rho\rangle_{\mathrm{GB}} = 1.994\,\xi + 3.424$, $R^2 = 0.9974$ |
| **Slope convergence** | $\alpha \to 2$ (algebraically exact from discriminant) |
| **Jump–conduit correlation** | $\operatorname{Corr}(\Delta\xi, \Delta\rho_{\min}) = 0.981$ |
| **Bandwidth constancy** | $\mathrm{BW}_{\mathrm{GB}} = 1.09 \pm 0.20$ (translation-invariant) |
| **Non-emptiness** | 513/513 even $N \in [1024, 2048]$ have Goldbach pairs |

### Figures

| | |
|---|---|
| ![BSL](figures/fig1_BSL.png) | ![Propagation](figures/fig2_propagation.png) |
| **Figure 1.** Band Shifting Law | **Figure 2.** Conductor orbit propagation |

![Envelope](figures/fig3_envelope.png)
**Figure 3.** Goldbach band vs. composite envelope

## Repository Structure

```
├── README.md
├── LICENSE
├── .gitignore
├── paper/
│   ├── Dynamic_Stability_Goldbach_Locus_Band_Shifting.tex
│   └── Dynamic_Stability_Goldbach_Locus_Band_Shifting.pdf
├── figures/
│   ├── fig1_BSL.pdf / .png
│   ├── fig2_propagation.pdf / .png
│   └── fig3_envelope.pdf / .png
├── scripts/
│   ├── dynamic_stability_scan.py      # Main scan: BSL, correlations, figures
│   └── sensitivity_analysis.py        # R² across N-ranges, bandwidth vs density
└── data/
    └── propagation_data.csv           # 513-row dataset (N, rad_M, ρ statistics)
```

## Quick Start

```bash
# Reproduce all results and figures
python3 scripts/dynamic_stability_scan.py

# Run sensitivity analysis (extended to N=16384)
python3 scripts/sensitivity_analysis.py
```

**Dependencies:** Python ≥ 3.10, NumPy ≥ 1.24, Matplotlib ≥ 3.7. No other packages required.

## Series Context

This is **Paper #10** in the Titan Project conductor rigidity series:

| # | Paper | Repository | Zenodo |
|---|-------|-----------|--------|
| 7 | The Goldbach Mirror | [goldbach-mirror-conductor-rigidity](https://github.com/Ruqing1963/goldbach-mirror-conductor-rigidity) | [10.5281/zenodo.18684892](https://zenodo.org/records/18684892) |
| 8 | The Goldbach Mirror II | [goldbach-mirror-II-geometric-foundations](https://github.com/Ruqing1963/goldbach-mirror-II-geometric-foundations) | [10.5281/zenodo.18719056](https://zenodo.org/records/18719056) |
| 9 | The Algebraic Vacuum | [goldbach-algebraic-vacuum-zero-ramification](https://github.com/Ruqing1963/goldbach-algebraic-vacuum-zero-ramification) | [10.5281/zenodo.18720040](https://zenodo.org/records/18720040) |
| **10** | **Dynamic Stability** | **this repository** | pending |

## Citation

```bibtex
@misc{chen2026dynamicstability,
  author       = {Ruqing Chen},
  title        = {Dynamic Stability of the Goldbach Locus: Conductor Orbit
                  Propagation and the Band Shifting Law in {GSp}(4)},
  year         = {2026},
  publisher    = {Zenodo},
  howpublished = {\url{https://github.com/Ruqing1963/goldbach-dynamic-stability-band-shifting}}
}
```

## License

[MIT](LICENSE)
