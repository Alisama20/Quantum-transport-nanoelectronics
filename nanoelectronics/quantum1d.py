"""1D time-independent Schrödinger equation via finite differences.

The Hamiltonian is discretized on a uniform grid of n interior points in [0, L]:

    H_ij = (2t₀ + V_i) δ_ij  −  t₀ (δ_{i,j−1} + δ_{i,j+1})

where t₀ = ℏ²/(2m·Δx²) is the hopping energy.  Dirichlet BCs ψ(0) = ψ(L) = 0.

Energy units: eV.  Length units: nm.
"""

import numpy as np
from scipy.linalg import eigh

HBAR = 197.327    # eV·nm  (ℏc in natural units)
M_E  = 0.511e6    # eV/c²  (electron rest mass)


def build_hamiltonian(V: np.ndarray, L: float = 10.0) -> np.ndarray:
    """Tridiagonal finite-difference Hamiltonian matrix.

    Parameters
    ----------
    V : potential on the n interior grid points [eV]
    L : physical domain length [nm]
    """
    n  = len(V)
    dx = L / (n + 1)                       # grid spacing [nm]
    t0 = HBAR**2 / (2.0 * M_E * dx**2)    # hopping energy [eV]
    diag = 2.0 * t0 + V
    off  = np.full(n - 1, -t0)
    return np.diag(diag) + np.diag(off, -1) + np.diag(off, 1)


def solve(V: np.ndarray, L: float = 10.0, n_states: int = 8):
    """Eigenvalues and eigenvectors of the 1D Hamiltonian.

    Parameters
    ----------
    V        : potential [eV], length n  (interior points)
    L        : domain length [nm]
    n_states : number of lowest eigenstates to return

    Returns
    -------
    E   : eigenvalues [eV], shape (n_states,)
    psi : physically normalised eigenvectors [nm^-½], shape (n, n_states)
          satisfying ∫|ψ_k|² dx = 1
    x   : grid positions [nm], shape (n,)
    """
    n  = len(V)
    dx = L / (n + 1)
    H  = build_hamiltonian(V, L)
    ns = min(n_states, n)
    E, vecs = eigh(H, subset_by_index=[0, ns - 1])
    psi = vecs / np.sqrt(dx)   # normalise to continuous limit
    x   = np.linspace(dx, L - dx, n)
    return E, psi, x


def particle_in_box_exact(n: int, L: float = 10.0) -> float:
    """Exact energy level E_n = n²π²ℏ²/(2mL²) [eV] (n = 1, 2, 3, …)."""
    return (n * np.pi * HBAR / L)**2 / (2.0 * M_E)


def make_barrier(n: int, V0: float, x_lo: float = 3.0,
                 x_hi: float = 7.0, L: float = 10.0) -> np.ndarray:
    """Rectangular potential barrier of height V0 in [x_lo, x_hi]."""
    x  = np.linspace(L / (n + 1), L - L / (n + 1), n)
    V  = np.zeros(n)
    V[(x >= x_lo) & (x <= x_hi)] = V0
    return V


def make_double_well(n: int, V0: float = 2.0, L: float = 10.0) -> np.ndarray:
    """Symmetric double well: two barriers flanking a central well."""
    x = np.linspace(L / (n + 1), L - L / (n + 1), n)
    V = np.zeros(n)
    V[(x >= 2.5) & (x <= 3.5)] = V0
    V[(x >= 6.5) & (x <= 7.5)] = V0
    return V
