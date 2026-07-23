"""
Calculate the mean consecutive level-spacing ratio for the 7NN lattice.

The script:
1. Builds the clean 7NN Hamiltonian.
2. Adds box-distributed on-site disorder.
3. Computes eigenvalues closest to the band center E = 0.
4. Calculates the consecutive level-spacing ratio.
5. Averages over disorder realizations.
6. Saves the results in CSV format.

The parameters below are small test values. The production parameters
used in the research calculation will be added after the test succeeds.
"""

from __future__ import annotations

import multiprocessing as mp
import sys
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import eigsh

from lattice_7nn import build_7nn_hamiltonian


# ------------------------------------------------------
# Import shared localization utilities
# ------------------------------------------------------
CURRENT_DIRECTORY = Path(__file__).resolve().parent
REPOSITORY_ROOT = CURRENT_DIRECTORY.parents[1]
SHARED_PYTHON_DIRECTORY = REPOSITORY_ROOT / "shared" / "python"

sys.path.insert(0, str(SHARED_PYTHON_DIRECTORY))

from localization_utils import apply_box_disorder  # noqa: E402


# ======================================================
# TEST PARAMETERS
# ======================================================
SYSTEM_SIZES = [6]

HOPPING = 1.0

NUMBER_OF_REALIZATIONS = 2

NUMBER_OF_DISORDER_VALUES = 3
DISORDER_VALUES = np.linspace(
    1.0,
    25.0,
    NUMBER_OF_DISORDER_VALUES,
)

NUMBER_OF_EIGENVALUES = 20

NUMBER_OF_CORES = 1

BASE_RANDOM_SEED = 24680


# ======================================================
# LEVEL-SPACING RATIO
# ======================================================
def calculate_spacing_ratios(
    eigenvalues: np.ndarray,
) -> np.ndarray:
    """
    Calculate consecutive level-spacing ratios.

    For ordered eigenvalues E_i, the spacings are

        s_i = E_{i+1} - E_i,

    and the consecutive spacing ratio is

        r_i = min(s_i, s_{i+1}) / max(s_i, s_{i+1}).
    """

    sorted_eigenvalues = np.sort(np.real(eigenvalues))
    spacings = np.diff(sorted_eigenvalues)

    smaller_spacings = np.minimum(
        spacings[:-1],
        spacings[1:],
    )

    larger_spacings = np.maximum(
        spacings[:-1],
        spacings[1:],
    )

    # Avoid division by zero in the unlikely case of exactly
    # degenerate numerical eigenvalues.
    valid = larger_spacings > 0.0

    if not np.any(valid):
        raise RuntimeError(
            "No valid level-spacing ratios could be calculated."
        )

    return smaller_spacings[valid] / larger_spacings[valid]


# ======================================================
# SINGLE DISORDER REALIZATION
# ======================================================
def calculate_single_realization(
    task: tuple,
) -> float:
    """
    Calculate the mean spacing ratio for one disorder realization.
    """

    (
        clean_hamiltonian,
        disorder_strength,
        number_of_eigenvalues,
        seed,
    ) = task

    disordered_hamiltonian = apply_box_disorder(
        clean_hamiltonian=clean_hamiltonian,
        disorder_strength=disorder_strength,
        seed=seed,
    )

    matrix_dimension = disordered_hamiltonian.shape[0]

    eigenvalues = eigsh(
        disordered_hamiltonian,
        k=number_of_eigenvalues,
        sigma=0.0,
        which="LM",
        return_eigenvectors=False,
        ncv=min(
            matrix_dimension - 1,
            max(2 * number_of_eigenvalues + 1, 60),
        ),
        maxiter=20000,
        tol=1.0e-8,
    )

    spacing_ratios = calculate_spacing_ratios(eigenvalues)

    return float(np.mean(spacing_ratios))


# ======================================================
# MAIN CALCULATION
# ======================================================
def main() -> None:
    """Run the 7NN level-spacing calculation."""

    output_directory = REPOSITORY_ROOT / "7NN" / "data"
    output_directory.mkdir(parents=True, exist_ok=True)

    seed_sequence = np.random.SeedSequence(BASE_RANDOM_SEED)

    print(
        "\n========== Starting 7NN level-spacing "
        "calculation =========="
    )

    for lattice_size in SYSTEM_SIZES:
        print(f"\n========== L = {lattice_size} ==========")

        clean_hamiltonian = build_7nn_hamiltonian(
            nx=lattice_size,
            ny=lattice_size,
            nz=lattice_size,
            hopping=HOPPING,
        )

        number_of_sites = clean_hamiltonian.shape[0]

        if NUMBER_OF_EIGENVALUES >= number_of_sites:
            raise ValueError(
                "NUMBER_OF_EIGENVALUES must be smaller than "
                f"the Hamiltonian dimension {number_of_sites}."
            )

        output_rows = []

        for disorder_strength in DISORDER_VALUES:
            child_sequences = seed_sequence.spawn(
                NUMBER_OF_REALIZATIONS
            )

            realization_seeds = [
                int(sequence.generate_state(1)[0])
                for sequence in child_sequences
            ]

            tasks = [
                (
                    clean_hamiltonian,
                    float(disorder_strength),
                    NUMBER_OF_EIGENVALUES,
                    seed,
                )
                for seed in realization_seeds
            ]

            if NUMBER_OF_CORES == 1:
                realization_mean_ratios = [
                    calculate_single_realization(task)
                    for task in tasks
                ]
            else:
                with mp.Pool(processes=NUMBER_OF_CORES) as pool:
                    realization_mean_ratios = pool.map(
                        calculate_single_realization,
                        tasks,
                    )

            mean_ratio = float(
                np.mean(realization_mean_ratios)
            )

            if NUMBER_OF_REALIZATIONS > 1:
                standard_error = float(
                    np.std(
                        realization_mean_ratios,
                        ddof=1,
                    )
                    / np.sqrt(NUMBER_OF_REALIZATIONS)
                )
            else:
                standard_error = 0.0

            output_rows.append(
                [
                    disorder_strength,
                    mean_ratio,
                    standard_error,
                ]
            )

            print(
                f"L = {lattice_size}, "
                f"W = {disorder_strength:.2f}, "
                f"<r> = {mean_ratio:.6f}"
            )

        output_array = np.asarray(output_rows, dtype=float)

        output_file = (
            output_directory
            / f"level_spacing_7nn_L{lattice_size}.csv"
        )

        np.savetxt(
            output_file,
            output_array,
            delimiter=",",
            header="W,mean_spacing_ratio,standard_error",
            comments="",
        )

        print(f"Saved: {output_file}")

    print("\n========== Calculation completed ==========")


if __name__ == "__main__":
    mp.freeze_support()
    main()