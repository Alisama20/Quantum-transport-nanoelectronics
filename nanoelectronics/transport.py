"""Landauer-Büttiker quantum transport through a single-level molecular junction.

Model (after Datta, Quantum Transport, Ch. 2):
  - Single molecular orbital (LUMO) at energy ε₀ = E_H + U
  - Lorentzian transmission broadened by coupling to source Γ_S and drain Γ_D
  - Non-equilibrium occupation weighted by transit times τ_S = ℏ/Γ_S, τ_D = ℏ/Γ_D
  - Self-consistent Hartree shift U = (q²/C_Σ)·(N − N₀)  [Coulomb blockade]

Current:
    I = q ∫ D(E−U, ε₀, Γ) · (f_S − f_D) / (τ_S + τ_D) dE

Occupation:
    N = ∫ D(E−U, ε₀, Γ) · [τ_D·f_S(E) + τ_S·f_D(E)] / (τ_S + τ_D) dE
"""

import numpy as np

K_B     = 8.617333262e-5    # eV/K
Q_E     = 1.602176634e-19   # C
HBAR_EV = 6.582119569e-16   # eV·s


def fermi(E: np.ndarray, mu: float, T: float) -> np.ndarray:
    return 1.0 / (np.exp(np.clip((E - mu) / (K_B * T), -500, 500)) + 1.0)


def lorentzian(E: np.ndarray, E_H: float, gamma: float) -> np.ndarray:
    """Lorentzian (Breit-Wigner) DOS normalised to 1 eV⁻¹ per level."""
    return (gamma / np.pi) / ((E - E_H)**2 + (gamma / 2.0)**2)


def iv_curve(E_F: float = -5.0,
             E_H: float = -4.7,
             gamma_S: float = 0.1,
             gamma_D: float = 0.1,
             qqC: float = 0.0,
             T: float = 298.0,
             V_min: float = -4.0,
             V_max: float = 4.0,
             n_V: int = 80,
             eta: float = 0.5,
             n_E: int = 10_000):
    """Compute I-V characteristic and differential conductance.

    Parameters
    ----------
    E_F           : equilibrium Fermi level [eV]
    E_H           : LUMO/HOMO bare energy [eV]
    gamma_S/D     : source/drain level broadening [eV]  (Γ = Γ_S + Γ_D)
    qqC           : Coulomb charging energy q²/C_Σ [eV/electron]
                    0 → non-interacting; larger values → Coulomb blockade
    T             : temperature [K]
    eta           : capacitive lever arm (fraction of bias at source junction)
    n_E           : integration grid points

    Returns
    -------
    V_arr : bias voltages [V]
    I_arr : currents [A]
    G_arr : differential conductance dI/dV [A/V]
    """
    gamma  = gamma_S + gamma_D
    tau_S  = HBAR_EV / gamma_S
    tau_D  = HBAR_EV / gamma_D
    tau    = tau_S + tau_D

    span   = 10.0 * gamma
    E_arr  = np.linspace(E_H - span, E_H + span, n_E)
    dE     = E_arr[1] - E_arr[0]

    # Equilibrium occupation N₀ (V = 0)
    D0  = lorentzian(E_arr, E_H, gamma)
    N_0 = np.trapz(D0 * fermi(E_arr, E_F, T), E_arr)

    V_arr = np.linspace(V_min, V_max, n_V)
    I_arr = np.zeros(n_V)

    for idx, V in enumerate(V_arr):
        mu_S = E_F + eta * V
        mu_D = E_F - (1.0 - eta) * V
        f_S  = fermi(E_arr, mu_S, T)
        f_D  = fermi(E_arr, mu_D, T)

        # Self-consistent Hartree shift U
        U = 0.0
        for _ in range(500):
            U_old = U
            D  = lorentzian(E_arr, E_H + U, gamma)
            N  = np.trapz(D * (tau_D * f_S + tau_S * f_D) / tau, E_arr)
            U  = U_old + 0.01 * (qqC * (N - N_0) - U_old)
            if abs(U - U_old) < 1e-5:
                break

        D = lorentzian(E_arr, E_H + U, gamma)
        I_arr[idx] = Q_E / tau * np.trapz(D * (f_S - f_D), E_arr)

    G_arr = np.gradient(I_arr, V_arr)
    return V_arr, I_arr, G_arr
