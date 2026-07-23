"""
Common localization utilities used by the 7NN and 8NN calculations.
"""

import numpy as np
from scipy.sparse import csr_matrix


def apply_box_disorder(
    clean_hamiltonian: csr_matrix,
    disorder_strength: float,
    seed: int | None = None,
) -> csr_matrix:
    """
    Add box-distributed on-site disorder to a clean Hamiltonian.

    The site energies are sampled uniformly from

        [-W/2, W/2],

    where W is the disorder strength.

    Parameters
    ----------
    clean_hamiltonian
        Clean sparse Hamiltonian.
    disorder_strength
        Disorder strength W.
    seed
        Random seed for reproducibility.

    Returns
    -------
    scipy.sparse.csr_matrix
        Hamiltonian containing random diagonal disorder.
    """

    if disorder_strength < 0:
        raise ValueError("Disorder strength must be non-negative.")

    rng = np.random.default_rng(seed)
    number_of_sites = clean_hamiltonian.shape[0]

    disorder = rng.uniform(
        low=-disorder_strength / 2.0,
        high=disorder_strength / 2.0,
        size=number_of_sites,
    )

    disordered_hamiltonian = clean_hamiltonian.copy().tolil()
    disordered_hamiltonian.setdiag(disorder)

    return disordered_hamiltonian.tocsr()


def compute_ipr(eigenvectors: np.ndarray) -> np.ndarray:
    """
    Calculate the inverse participation ratio of each eigenstate.

    Eigenvectors are assumed to be stored column-wise, as returned
    by scipy.sparse.linalg.eigsh.

    Parameters
    ----------
    eigenvectors
        Array whose columns contain normalized eigenvectors.

    Returns
    -------
    numpy.ndarray
        IPR value for each eigenstate.
    """

    if eigenvectors.ndim != 2:
        raise ValueError("eigenvectors must be a two-dimensional array.")

    return np.sum(np.abs(eigenvectors) ** 4, axis=0)


def compute_mean_log_ipr(eigenvectors: np.ndarray) -> float:
    """
    Calculate the mean logarithmic IPR over the selected eigenstates.

    The typical IPR is obtained later as

        IPR_typ = exp(mean_log_ipr).
    """

    ipr_values = compute_ipr(eigenvectors)

    if np.any(ipr_values <= 0):
        raise ValueError("All IPR values must be positive.")

    return float(np.mean(np.log(ipr_values)))