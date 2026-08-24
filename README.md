# Influence of local connectivity on the Anderson transition beyond bandwidth broadening

This repository contains numerical codes and processed data associated with the study:

**Influence of local connectivity on the Anderson transition beyond bandwidth broadening**

The work investigates Anderson localization in modified three-dimensional cubic lattices with coordination numbers (7) and (8). The localization transition is characterized using the consecutive level-spacing ratio and the inverse participation ratio (IPR).

For the IPR finite-size scaling analysis, the **typical IPR**

[
\mathrm{IPR}_{\mathrm{typ}}
===========================

\exp\left[
\left\langle
\ln(\mathrm{IPR})
\right\rangle
\right]
]

is used.

Python is used to construct the lattice Hamiltonians and perform the numerical diagonalization and localization calculations. MATLAB scripts are provided for plotting the processed data and reproducing the corresponding finite-size scaling analyses.

---

## Models

The single-particle Anderson Hamiltonian is

[
H=-t\sum_{\langle i,j\rangle}
\left(c_i^\dagger c_j+c_j^\dagger c_i\right)
+\sum_i \epsilon_i c_i^\dagger c_i ,
]

where:

* (t=1) is the hopping amplitude;
* (\langle i,j\rangle) denotes connected site pairs (hopping bonds);
* (\epsilon_i) is the random on-site energy;
* (\epsilon_i) is sampled independently from the uniform distribution
  ([-W/2,W/2]);
* (W) is the disorder strength;
* periodic boundary conditions are used.

Two modified cubic lattice realizations are considered.

### 7NN lattice

The 7NN lattice is obtained by supplementing the ordinary simple-cubic hopping network with an additional alternating pattern of diagonal hopping bonds.

The resulting lattice has a uniform coordination number

[
\mathrm{CN}=7
]

under periodic boundary conditions.

### 8NN lattice

The 8NN lattice is constructed by introducing two additional diagonal hopping connections according to a checkerboard-type parity pattern.

The resulting lattice has a uniform coordination number

[
\mathrm{CN}=8
]

under periodic boundary conditions.

The lattice-generation scripts verify that the Hamiltonian is symmetric and that every site has the required coordination number.

---

## Localization quantities

### Inverse participation ratio

For a normalized eigenstate (\psi_n), the inverse participation ratio is

[
\mathrm{IPR}_n
==============

\sum_i |\psi_n(i)|^4.
]

For the band-center finite-size scaling analysis, the logarithmic IPR is averaged over the selected eigenstates and disorder realizations. The typical IPR is then calculated as

[
\mathrm{IPR}_{\mathrm{typ}}
===========================

\exp\left[
\left\langle
\ln(\mathrm{IPR})
\right\rangle
\right].
]

The corresponding scaling form is

[
\mathrm{IPR}_{\mathrm{typ}}(W,L)
================================

L^{-D_2}
F\left[(W-W_c)L^{1/\nu}\right],
]

where (W_c) is the critical disorder strength, (\nu) is the correlation-length exponent, and (D_2) is the correlation (multifractal) dimension.

The processed data therefore contain both

[
\left\langle\ln(\mathrm{IPR})\right\rangle
]

and

[
\mathrm{IPR}_{\mathrm{typ}}
===========================

\exp\left[
\left\langle\ln(\mathrm{IPR})\right\rangle
\right].
]

### Consecutive level-spacing ratio

For ordered eigenvalues (E_n), the consecutive level spacings are

[
\delta_n=E_{n+1}-E_n.
]

The consecutive level-spacing ratio is

[
r_n=
\frac{\min(\delta_n,\delta_{n+1})}
{\max(\delta_n,\delta_{n+1})}.
]

The disorder-averaged quantity

[
\langle r\rangle
]

is used to distinguish extended and localized spectral statistics.

The finite-size scaling variable is

[
(W-W_c)L^{1/\nu}.
]

---

## Repository structure

```text
anderson-localization-modified-lattices/
│
├── README.md
├── .gitignore
│
├── 7NN/
│   ├── python/
│   │   ├── lattice_7nn.py
│   │   ├── run_ipr_7nn.py
│   │   └── run_level_spacing_7nn.py
│   │
│   ├── matlab/
│   │   ├── plot_ipr_7nn.m
│   │   ├── collapse_ipr_7nn.m
│   │   ├── plot_level_spacing_7nn.m
│   │   └── collapse_level_spacing_7nn.m
│   │
│   └── data/
│       ├── ipr_7nn_L6.csv
│       ├── level_spacing_7nn_L6.csv
│       └── paper/
│           ├── ipr_7nn_paper.log
│           └── level_spacing_7nn_paper.log
│
├── 8NN/
│   ├── python/
│   │   ├── lattice_8nn.py
│   │   ├── run_ipr_8nn.py
│   │   └── run_level_spacing_8nn.py
│   │
│   ├── matlab/
│   │   ├── plot_ipr_8nn.m
│   │   ├── collapse_ipr_8nn.m
│   │   ├── plot_level_spacing_8nn.m
│   │   └── collapse_level_spacing_8nn.m
│   │
│   └── data/
│       ├── ipr_8nn_L6.csv
│       ├── level_spacing_8nn_L6.csv
│       └── paper/
│           ├── ipr_8nn_paper.log
│           └── level_spacing_8nn_paper.log
│
└── shared/
    └── python/
        └── localization_utils.py
```

---

## Python requirements

The Python calculations require:

* Python 3
* NumPy
* SciPy

Install the required packages with

```bash
python -m pip install numpy scipy
```

On Windows, the Python launcher may alternatively be used:

```bash
py -m pip install numpy scipy
```

---

## Running the Python calculations

Run the commands from the root directory of the repository.

### 7NN typical IPR

```bash
py 7NN/python/run_ipr_7nn.py
```

### 7NN level-spacing ratio

```bash
py 7NN/python/run_level_spacing_7nn.py
```

### 8NN typical IPR

```bash
py 8NN/python/run_ipr_8nn.py
```

### 8NN level-spacing ratio

```bash
py 8NN/python/run_level_spacing_8nn.py
```

For macOS or Linux, replace `py` with `python3`.

The Python scripts calculate eigenvalues and eigenvectors close to the band center, (E=0), using sparse diagonalization.

---

## Important note about calculation parameters

The default parameters included in the Python runner scripts are small test settings intended primarily to verify that the programs execute correctly.

These demonstration calculations are not expected to reproduce the statistically converged numerical results reported in the manuscript.

The production calculations used in the study employed:

* larger lattice sizes;
* denser disorder grids;
* more eigenstates near the band center;
* many independent disorder realizations;
* substantially greater computing time and memory.

The processed numerical results used for the 7NN and 8NN band-center scaling analyses are provided in the corresponding `data/paper/` directories.

---

## MATLAB plotting and scaling scripts

The MATLAB scripts automatically locate the required data files relative to their own locations. The MATLAB working directory therefore does not need to be manually changed to the data directory.

### 7NN analysis

Run:

```text
7NN/matlab/plot_ipr_7nn.m
7NN/matlab/collapse_ipr_7nn.m
7NN/matlab/plot_level_spacing_7nn.m
7NN/matlab/collapse_level_spacing_7nn.m
```

### 8NN analysis

Run:

```text
8NN/matlab/plot_ipr_8nn.m
8NN/matlab/collapse_ipr_8nn.m
8NN/matlab/plot_level_spacing_8nn.m
8NN/matlab/collapse_level_spacing_8nn.m
```

The IPR plotting scripts reproduce the averaged logarithmic IPR and typical-IPR curves contained in the repository.

The level-spacing scripts reproduce the corresponding disorder-dependent level-spacing-ratio curves.

The collapse scripts reproduce the finite-size scaling analyses using the critical parameters specified in the corresponding MATLAB files.

---

## Data organization

Two types of numerical data are included.

### Demonstration data

Files such as

```text
ipr_7nn_L6.csv
level_spacing_7nn_L6.csv
ipr_8nn_L6.csv
level_spacing_8nn_L6.csv
```

were generated using small test calculations.

They are included to demonstrate the output format and verify that the Python workflow functions correctly. They should not be interpreted as the final production data reported in the manuscript.

### Paper data

The files inside

```text
7NN/data/paper/
8NN/data/paper/
```

contain the processed numerical data used for the 7NN and 8NN band-center level-spacing and typical-IPR plots and finite-size scaling analyses reported in the manuscript.

For the IPR calculations, these files contain the averaged logarithmic IPR,

[
\langle\ln(\mathrm{IPR})\rangle,
]

and the corresponding typical IPR,

[
\mathrm{IPR}_{\mathrm{typ}}
===========================

\exp[\langle\ln(\mathrm{IPR})\rangle].
]

Individual Hamiltonian matrices, eigenvectors, and results from every individual disorder realization are not included because of their substantially larger storage requirements.

---

## Reproducing the band-center scaling analyses

To reproduce one of the supplied plots or scaling analyses:

1. Clone or download the complete repository.
2. Keep the directory structure unchanged.
3. Open the corresponding MATLAB script.
4. Run the script.
5. The script will automatically read the required processed data from the appropriate `data/paper/` directory.

For example, the 7NN typical-IPR scaling collapse is generated using

```text
7NN/matlab/collapse_ipr_7nn.m
```

with the processed data stored in

```text
7NN/data/paper/ipr_7nn_paper.log
```

---

## Scope of the archived data

The repository is intended primarily to document and reproduce the numerical workflow associated with the **7NN and 8NN band-center finite-size scaling analyses**.

It contains:

* lattice-generation routines for the modified 7NN and 8NN lattices;
* typical-IPR calculations near the band center;
* consecutive level-spacing-ratio calculations near the band center;
* processed production data for the corresponding finite-size scaling analyses;
* MATLAB scripts for plotting and scaling these quantities.

Large raw datasets, individual disorder realizations, complete eigenvector datasets, and Hamiltonian matrices from the production calculations are not included because of their storage requirements.

Additional numerical data associated with the study may be obtained from the authors upon reasonable request.

---

## Numerical considerations

Sparse diagonalization near the band center can occasionally experience convergence difficulties, particularly for clean or nearly degenerate systems.

The Python scripts use shift-invert sparse diagonalization with a small nonzero spectral shift and suitable iteration parameters to improve numerical convergence.

Large production calculations should preferably be performed using high-performance computing resources.

---

## Citation

This repository accompanies the article:

> **Influence of local connectivity on the Anderson transition beyond bandwidth broadening**

**Mohammed Zahid Malik and Raja Ghosh**

*Physical Review B* (2026).

The complete volume, article number, and final DOI will be added after publication.

When using the codes or processed numerical data, please cite the associated article and this repository.

---

## Authors

Repository maintained by **Mohammed Zahid Malik**.

The associated research article is authored by:

* Mohammed Zahid Malik
* Raja Ghosh

---

## License

No reuse license has currently been assigned.

Until a license is added, the code and data remain protected by default copyright and should not be redistributed or reused without permission.
