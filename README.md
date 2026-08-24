# Influence of local connectivity on the Anderson transition beyond bandwidth broadening

This repository contains the numerical codes and processed data used in the study:

**Influence of local connectivity on the Anderson transition beyond bandwidth broadening**

The project investigates Anderson localization in modified three-dimensional cubic lattices with coordination numbers 7 and 8. The numerical analysis uses the typical inverse participation ratio and the consecutive level-spacing ratio to characterize the transition between extended and localized states.

Python is used to construct the Hamiltonians and calculate the localization quantities. MATLAB is used to reproduce the raw-data plots and finite-size scaling-collapse figures.

---

## Models

The single-particle Anderson Hamiltonian is

\[
H=-t\sum_{\langle i,j\rangle}
\left(c_i^\dagger c_j+c_j^\dagger c_i\right)
+\sum_i \epsilon_i c_i^\dagger c_i ,
\]

where:

- \(t=1\) is the hopping amplitude;
- \(\epsilon_i\) is the on-site disorder;
- \(\epsilon_i\) is sampled from a uniform box distribution
  \([-W/2,W/2]\);
- \(W\) is the disorder strength;
- periodic boundary conditions are used.

Two modified cubic lattices are considered.

### 7NN lattice

Each lattice site has seven nearest-neighbour connections. The ordinary cubic-lattice hopping is supplemented with an additional diagonal hopping pattern.

### 8NN lattice

Each lattice site has eight nearest-neighbour connections. Two additional diagonal connections are introduced through a checkerboard-type pattern.

The lattice-generation scripts verify that the Hamiltonian is symmetric and that every lattice site has the required coordination number.

---

## Localization quantities

### Typical inverse participation ratio

For a normalized eigenstate \(\psi_n\), the inverse participation ratio is

\[
\mathrm{IPR}_n=\sum_i |\psi_n(i)|^4.
\]

The typical IPR is calculated as

\[
\mathrm{IPR}_{\mathrm{typ}}
=
\exp\left[
\left\langle
\ln(\mathrm{IPR})
\right\rangle
\right].
\]

The finite-size scaling form used for the IPR collapse is

\[
\mathrm{IPR}_{\mathrm{typ}}(W,L)
=
L^{-\tau}
F\left[(W-W_c)L^{1/\nu}\right].
\]

### Consecutive level-spacing ratio

For ordered eigenvalues \(E_n\), the level spacings are

\[
\delta_n=E_{n+1}-E_n.
\]

The consecutive spacing ratio is

\[
r_n=
\frac{\min(\delta_n,\delta_{n+1})}
{\max(\delta_n,\delta_{n+1})}.
\]

The disorder-averaged value \(\langle r\rangle\) distinguishes extended and localized spectral statistics.

The finite-size scaling variable is

\[
(W-W_c)L^{1/\nu}.
\]

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

- Python 3
- NumPy
- SciPy

Install the required packages with

```bash
python -m pip install numpy scipy
```

On Windows, the Python launcher may be used:

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

The Python scripts calculate eigenvalues and eigenvectors close to the band center, \(E=0\), using sparse diagonalization.

---

## Important note about calculation parameters

The default parameters currently included in the Python runner scripts are small test settings intended to verify that the programs run correctly.

For example, the demonstration files use a small lattice and a small number of disorder realizations. These test calculations are not expected to reproduce the statistically converged results reported in the manuscript.

Large-scale production calculations require:

- larger lattice sizes;
- a denser disorder grid;
- more eigenvalues near the band center;
- many independent disorder realizations;
- substantially greater computing time and memory.

The processed numerical data used for the manuscript figures are provided in the corresponding `data/paper/` directories.

---

## MATLAB plotting scripts

The MATLAB scripts automatically locate the required data files relative to their own locations. Therefore, the MATLAB working directory does not need to be manually set to the data folder.

### 7NN figures

Run:

```text
7NN/matlab/plot_ipr_7nn.m
7NN/matlab/collapse_ipr_7nn.m
7NN/matlab/plot_level_spacing_7nn.m
7NN/matlab/collapse_level_spacing_7nn.m
```

### 8NN figures

Run:

```text
8NN/matlab/plot_ipr_8nn.m
8NN/matlab/collapse_ipr_8nn.m
8NN/matlab/plot_level_spacing_8nn.m
8NN/matlab/collapse_level_spacing_8nn.m
```

The plotting scripts reproduce the raw disorder-dependent curves. The collapse scripts apply the critical parameters specified near the beginning of each MATLAB file.

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

were generated using small test calculations. They are included to demonstrate the output format and verify that the Python workflow functions correctly.

### Paper data

The files inside

```text
7NN/data/paper/
8NN/data/paper/
```

contain the processed disorder-averaged data used directly for the manuscript plots and finite-size scaling analysis.

Individual Hamiltonian matrices, eigenvectors, and data from every disorder realization are not included because of their substantially larger storage requirements.

---

## Reproducing the figures

To reproduce a figure:

1. Clone or download the complete repository.
2. Keep the directory structure unchanged.
3. Open the corresponding MATLAB script.
4. Run the script.
5. The script will read the appropriate file from the `data/paper/` directory.

For example, the 7NN IPR-collapse figure is generated by

```text
7NN/matlab/collapse_ipr_7nn.m
```

using

```text
7NN/data/paper/ipr_7nn_paper.log
```

---

## Numerical considerations

Sparse diagonalization near the band center can occasionally experience convergence difficulties, especially for clean or nearly degenerate systems.

The Python scripts use shift-invert diagonalization with a small nonzero spectral shift and increased iteration parameters to improve numerical convergence.

Large production calculations should preferably be performed on a high-performance computing system.

---

## Citation

This repository accompanies the manuscript:

> > **Influence of local connectivity on the Anderson transition beyond bandwidth broadening**

The complete journal citation and DOI will be added after publication.

When using the codes or processed data, please cite the associated manuscript and this repository.

---

## Authors

Repository maintained by **Zahid Malik**.

Additional author and institutional information can be added after publication.

---

## License

No reuse license has currently been assigned.

Until a license is added, the code and data remain protected by default copyright and should not be redistributed or reused without permission.
