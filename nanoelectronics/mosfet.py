"""Ballistic MOSFET model with 2D electron gas in the channel.

Current formula (Lundstrom ballistic FET model):

    I_DS = I₀ · √m · ∫_{E_C+U}^∞ √(E − E_C − U) · [f_S(E) − f_D(E)] dE

where
    I₀  = q · W / (π² · ℏ²)
    m   : electron effective mass
    U   : self-consistent gate-induced band shift [eV]

Self-consistent band shift (linear screening):
    U = −V_GS·(C_G/C_Σ) − V_DS·(C_D/C_Σ) + q·(N − N₀)/C_Σ

with N = ∫ g₀ · [f_S(E) + f_D(E)] dE,  g₀ = m·L·W/(2π·ℏ²)  [2D DOS × area].
"""

import numpy as np

K_B  = 8.617333262e-5    # eV/K
Q_E  = 1.602176634e-19   # C
HBAR = 1.054571817e-34   # J·s
M_E  = 9.109e-31         # kg (free electron mass)


def _fermi(E_eV: np.ndarray, mu_eV: float, T: float) -> np.ndarray:
    return 1.0 / (np.exp(np.clip((E_eV - mu_eV) / (K_B * T), -500, 500)) + 1.0)


def iv_ballistic(E_F: float = -5.0,
                 E_C: float = -4.7,
                 L: float = 40e-9,
                 W: float = 120e-9,
                 C_G: float = 0.1e-15,
                 C_S: float = 0.0,
                 C_D: float = 0.0,
                 T: float = 5.0,
                 V_GS: float = 0.5,
                 V_min: float = 0.0,
                 V_max: float = 0.5,
                 n_V: int = 51,
                 n_E: int = 100_000):
    """Ballistic I_DS vs V_DS at fixed V_GS.

    Parameters
    ----------
    E_F, E_C   : Fermi energy and conduction band edge [eV]
    L, W       : channel length and width [m]
    C_G, C_S, C_D : gate, source, drain capacitances [F]
    T          : temperature [K]
    V_GS       : gate-source voltage [V]

    Returns V_arr [V], I_arr [μA].
    """
    C_ES = C_S + C_D + C_G
    I_0  = Q_E * W / (np.pi**2 * HBAR**2)       # prefactor [A / (kg^0.5 J^1.5)]
    g_0  = 0.5 * M_E * L * W / (np.pi * HBAR**2) # 2D DOS [J^-1]

    # Equilibrium carrier density N₀ (V_DS = 0, U = 0)
    E_max = E_C + 5.0                # upper integration limit [eV]
    E_eq  = np.linspace(E_C, E_max, n_E)
    dE_eq = E_eq[1] - E_eq[0]
    N_0   = g_0 * dE_eq * Q_E * np.sum(_fermi(E_eq, E_F, T))

    V_arr = np.linspace(V_min, V_max, n_V)
    I_arr = np.zeros(n_V)

    for idx, V_DS in enumerate(V_arr):
        mu_s = E_F
        mu_d = E_F - V_DS

        # Self-consistent band shift U
        U = 0.0
        for _ in range(2000):
            U_old = U
            E_sc  = np.linspace(E_C + U, E_max, n_E)
            dE    = E_sc[1] - E_sc[0]
            N_S   = g_0 * dE * Q_E * np.sum(_fermi(E_sc, mu_s, T))
            N_D   = g_0 * dE * Q_E * np.sum(_fermi(E_sc, mu_d, T))
            N     = N_S + N_D
            U_calc = (-V_GS * (C_G / C_ES)
                      - V_DS * (C_D / C_ES)
                      + Q_E * (N - N_0) / C_ES)
            U  = U_old + 0.01 * (U_calc - U_old)
            if abs(U - U_old) < 1e-5:
                break

        # Current integral: I = I₀ √m ∫ √(E−E_C−U) [f_S−f_D] dE_J
        E_sc  = np.linspace(E_C + U, E_max, n_E)
        dE    = E_sc[1] - E_sc[0]
        E_above_J = (E_sc - (E_C + U)) * Q_E          # energy above band edge [J]
        f_diff    = _fermi(E_sc, mu_s, T) - _fermi(E_sc, mu_d, T)
        suma      = np.sum(np.sqrt(E_above_J) * f_diff)
        I_arr[idx] = I_0 * np.sqrt(M_E) * dE * suma * Q_E

    return V_arr, I_arr * 1e6    # V, μA


def iv_semiclassical(V_GS: float = 0.5,
                     V_T: float = 0.3,
                     mu_n: float = 300e-4,
                     L: float = 40e-9,
                     W: float = 120e-9,
                     C_G: float = 0.1e-15,
                     V_min: float = 0.0,
                     V_max: float = 0.5,
                     n_V: int = 51):
    """Long-channel MOSFET (drift-diffusion / gradual channel approximation).

    Returns V_arr [V], I_arr [μA].
    """
    C_OX  = C_G / (W * L)     # oxide capacitance density [F/m²]
    beta  = mu_n * C_OX * (W / L)  # transconductance parameter [A/V²]
    V_arr = np.linspace(V_min, V_max, n_V)
    I_arr = np.zeros(n_V)
    if V_GS > V_T:
        V_dsat = V_GS - V_T
        for i, V_DS in enumerate(V_arr):
            if V_DS < V_dsat:
                I_arr[i] = beta * ((V_GS - V_T) * V_DS - 0.5 * V_DS**2)
            else:
                I_arr[i] = 0.5 * beta * (V_GS - V_T)**2
    return V_arr, I_arr * 1e6  # V, μA
