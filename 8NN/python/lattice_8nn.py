"""
Construction of the clean 8NN checkerboard lattice Hamiltonian.

The lattice has periodic boundary conditions. Both diagonals are added
on alternating xy plaquettes according to the parity condition
(i + j + k) % 2 == 0.
"""

from scipy.sparse import csr_matrix, lil_matrix


def build_8nn_hamiltonian(
    nx: int,
    ny: int,
    nz: int,
    hopping: float = 1.0,
) -> csr_matrix:
    """
    Build the clean 8NN tight-binding Hamiltonian.

    Parameters
    ----------
    nx, ny, nz
        Number of lattice sites along the x, y, and z directions.
        Even dimensions are required for the periodic checkerboard pattern.
    hopping
        Hopping amplitude for all bonds.

    Returns
    -------
    scipy.sparse.csr_matrix
        Sparse Hamiltonian of dimension
        (nx * ny * nz) x (nx * ny * nz).
    """

    if nx < 2 or ny < 2 or nz < 2:
        raise ValueError("Each lattice dimension must be at least 2.")

    if nx % 2 != 0 or ny % 2 != 0 or nz % 2 != 0:
        raise ValueError(
            "nx, ny, and nz must be even for the periodic "
            "8NN checkerboard pattern."
        )

    site_index = {}
    index = 0

    for k in range(nz):
        for i in range(nx):
            for j in range(ny):
                site_index[(i, j, k)] = index
                index += 1

    number_of_sites = nx * ny * nz
    adjacency = {
        site: set() for site in range(number_of_sites)
    }

    # Add the six simple-cubic nearest-neighbor bonds.
    # Only the positive directions are generated here, and every
    # connection is then inserted in both directions.
    for k in range(nz):
        for i in range(nx):
            for j in range(ny):
                site = site_index[(i, j, k)]

                positive_neighbors = [
                    ((i + 1) % nx, j, k),
                    (i, (j + 1) % ny, k),
                    (i, j, (k + 1) % nz),
                ]

                for neighbor_position in positive_neighbors:
                    neighbor = site_index[neighbor_position]

                    adjacency[site].add(neighbor)
                    adjacency[neighbor].add(site)

    # Add both diagonals on alternating xy plaquettes.
    for k in range(nz):
        for i in range(nx):
            for j in range(ny):
                if (i + j + k) % 2 == 0:
                    bottom_left = site_index[(i, j, k)]

                    bottom_right = site_index[
                        ((i + 1) % nx, j, k)
                    ]

                    top_left = site_index[
                        (i, (j + 1) % ny, k)
                    ]

                    top_right = site_index[
                        (
                            (i + 1) % nx,
                            (j + 1) % ny,
                            k,
                        )
                    ]

                    # First diagonal: bottom-left to top-right.
                    adjacency[bottom_left].add(top_right)
                    adjacency[top_right].add(bottom_left)

                    # Second diagonal: bottom-right to top-left.
                    adjacency[bottom_right].add(top_left)
                    adjacency[top_left].add(bottom_right)

    # Verify uniform coordination number CN = 8.
    coordination_numbers = [
        len(neighbors) for neighbors in adjacency.values()
    ]

    if any(cn != 8 for cn in coordination_numbers):
        raise RuntimeError(
            "The generated lattice does not have uniform CN = 8."
        )

    hamiltonian = lil_matrix(
        (number_of_sites, number_of_sites),
        dtype=float,
    )

    for site, neighbors in adjacency.items():
        for neighbor in neighbors:
            hamiltonian[site, neighbor] = -hopping

    return hamiltonian.tocsr()