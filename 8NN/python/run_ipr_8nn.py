"""
Calculate the typical inverse participation ratio for the 8NN lattice.

The script:
1. Builds the clean 8NN Hamiltonian.
2. Adds independent box-distributed on-site disorder.
3. Computes eigenstates closest to E = 0.
4. Calculates the mean logarithmic IPR.
5. Saves data in CSV format for later MATLAB plotting and scaling analysis.

The parameters below are small test values. They will be replaced with
the actual production parameters after the script is tested.
"""

from __future__ import annotations

import multiprocessing as mp
import sys
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import eigsh

from lattice_8nn import build_8nn_hamiltonian


# ------------------------------------------------------
# Import the shared localization utilities
# ------------------------------------------------------
CURRENT_DIRECTORY = Path(__file__).resolve().parent
REPOSITORY_ROOT = CURRENT_DIRECTORY.parents[1]
SHARED_PYTHON_DIRECTORY = REPOSITORY_ROOT / "shared" / "python"

sys.path.insert(0, str(SHARED_PYTHON_DIRECTORY))

from localization_utils import (  # noqa: E402
    apply_box_disorder,
    compute_mean_log_ipr,
)


# ======================================================
# TEST PARAMETERS
# ======================================================
SYSTEM_SIZES = [6]

HOPPING = 1.0

NUMBER_OF_REALIZATIONS = 2

NUMBER_OF_DISORDER_VALUES = 3
DISORDER_VALUES = np.linspace(0.0, 4.0, NUMBER_OF_DISORDER_VALUES)

NUMBER_OF_EIGENSTATES = 20

NUMBER_OF_CORES = 1

BASE_RANDOM_SEED = 12345


# ======================================================
# SINGLE DISORDER REALIZATION
# ======================================================
def calculate_single_realization(
    task: tuple,
) -> float:
    """
    Calculate the mean logarithmic IPR for one disorder realization.
    """

    clean_hamiltonian, disorder_strength, number_of_eigenstates, seed = task

    disordered_hamiltonian = apply_box_disorder(
        clean_hamiltonian=clean_hamiltonian,
        disorder_strength=disorder_strength,
        seed=seed,
    )

    _, eigenvectors = eigsh(
        disordered_hamiltonian,
        k=number_of_eigenstates,
        sigma=1.0e-8,
        which="LM",
        ncv=min(
            disordered_hamiltonian.shape[0] - 1,
            max(2 * number_of_eigenstates + 1, 60),
        ),
        maxiter=20000,
        tol=1.0e-8,
)

    return compute_mean_log_ipr(eigenvectors)


# ======================================================
# MAIN CALCULATION
# ======================================================
def main() -> None:
    """Run the 8NN IPR calculation and save the results."""

    output_directory = REPOSITORY_ROOT / "8NN" / "data"
    output_directory.mkdir(parents=True, exist_ok=True)

    seed_sequence = np.random.SeedSequence(BASE_RANDOM_SEED)

    print("\n========== Starting 8NN IPR calculation ==========")

    for lattice_size in SYSTEM_SIZES:
        print(f"\n========== L = {lattice_size} ==========")

        clean_hamiltonian = build_8nn_hamiltonian(
            nx=lattice_size,
            ny=lattice_size,
            nz=lattice_size,
            hopping=HOPPING,
        )

        number_of_sites = clean_hamiltonian.shape[0]

        if NUMBER_OF_EIGENSTATES >= number_of_sites:
            raise ValueError(
                "NUMBER_OF_EIGENSTATES must be smaller than "
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
                    NUMBER_OF_EIGENSTATES,
                    seed,
                )
                for seed in realization_seeds
            ]

            if NUMBER_OF_CORES == 1:
                mean_log_ipr_values = [
                    calculate_single_realization(task)
                    for task in tasks
                ]
            else:
                with mp.Pool(processes=NUMBER_OF_CORES) as pool:
                    mean_log_ipr_values = pool.map(
                        calculate_single_realization,
                        tasks,
                    )

            mean_log_ipr = float(np.mean(mean_log_ipr_values))
            typical_ipr = float(np.exp(mean_log_ipr))

            if NUMBER_OF_REALIZATIONS > 1:
                standard_error = float(
                    np.std(mean_log_ipr_values, ddof=1)
                    / np.sqrt(NUMBER_OF_REALIZATIONS)
                )
            else:
                standard_error = 0.0

            output_rows.append(
                [
                    disorder_strength,
                    mean_log_ipr,
                    typical_ipr,
                    standard_error,
                ]
            )

            print(
                f"W = {disorder_strength:.2f}, "
                f"mean ln(IPR) = {mean_log_ipr:.6f}, "
                f"IPR_typ = {typical_ipr:.6e}"
            )

        output_array = np.asarray(output_rows, dtype=float)

        output_file = (
            output_directory
            / f"ipr_8nn_L{lattice_size}.csv"
        )

        np.savetxt(
            output_file,
            output_array,
            delimiter=",",
            header=(
                "W,mean_log_ipr,typical_ipr,"
                "standard_error_mean_log_ipr"
            ),
            comments="",
        )

        print(f"Saved: {output_file}")

    print("\n========== Calculation completed ==========")


if __name__ == "__main__":
    mp.freeze_support()
    main()