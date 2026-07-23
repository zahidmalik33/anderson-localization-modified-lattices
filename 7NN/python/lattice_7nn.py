"""
Construction of the clean 7NN modified cubic lattice Hamiltonian.

The lattice has periodic boundary conditions and one additional
diagonal hopping connection per site.
"""

from scipy.sparse import csr_matrix, lil_matrix


def build_7nn_hamiltonian(
    nx: int,
    ny: int,
    nz: int,
    hopping: float = 1.0,
) -> csr_matrix:
    """
    Build the clean 7NN tight-binding Hamiltonian.

    Parameters
    ----------
    nx, ny, nz
        Number of lattice sites along the x, y, and z directions.
    hopping
        Nearest-neighbor and diagonal hopping amplitude.

    Returns
    -------
    scipy.sparse.csr_matrix
        Sparse Hamiltonian of dimension
        (nx * ny * nz) x (nx * ny * nz).
    """

    if nx < 2 or ny < 2 or nz < 2:
        raise ValueError("Each lattice dimension must be at least 2.")

    # The alternating diagonal pattern requires an even y dimension.
    if ny % 2 != 0:
        raise ValueError(
            "ny must be even for the periodic 7NN diagonal pattern."
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

    # Six nearest-neighbor connections with periodic boundaries.
    nearest_neighbor_shifts = [
        (-1, 0, 0),
        (1, 0, 0),
        (0, -1, 0),
        (0, 1, 0),
        (0, 0, -1),
        (0, 0, 1),
    ]

    for (i, j, k), site in site_index.items():
        for dx, dy, dz in nearest_neighbor_shifts:
            neighbor_position = (
                (i + dx) % nx,
                (j + dy) % ny,
                (k + dz) % nz,
            )

            neighbor = site_index[neighbor_position]
            adjacency[site].add(neighbor)

    # Add alternating diagonal connections in each xy plane.
    for k in range(nz):
        for i in range(nx):
            for j in range(0, ny, 2):
                site_a = site_index[((i + 1) % nx, j, k)]
                site_b = site_index[(i, (j + 1) % ny, k)]

                adjacency[site_a].add(site_b)
                adjacency[site_b].add(site_a)

    # Verify that every site has exactly seven neighbors.
    coordination_numbers = [
        len(neighbors) for neighbors in adjacency.values()
    ]

    if any(cn != 7 for cn in coordination_numbers):
        raise RuntimeError(
            "The generated lattice does not have uniform CN = 7."
        )

    hamiltonian = lil_matrix(
        (number_of_sites, number_of_sites),
        dtype=float,
    )

    for site, neighbors in adjacency.items():
        for neighbor in neighbors:
            hamiltonian[site, neighbor] = -hopping

    return hamiltonian.tocsr()