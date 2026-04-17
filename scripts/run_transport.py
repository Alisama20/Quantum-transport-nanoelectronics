"""Single-level molecular junction: I-V and conductance.

Two scenarios:
  1. Non-interacting (qqC = 0) vs Coulomb-blockaded junction — shows how
     charging shifts and suppresses the current.
  2. Temperature dependence — sharp conductance peak broadens with T.

Outputs:
    transport_iv.png           — I-V with and without charging
    transport_conductance.png  — dI/dV vs V at several temperatures
"""

import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import numpy as np
import matplotlib.pyplot as plt

from nanoelectronics import iv_curve

FIGDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "figures"))
os.makedirs(FIGDIR, exist_ok=True)

E_F = -5.0     # eV  (equilibrium Fermi level)
E_H = -4.7     # eV  (LUMO level)


def plot_iv_charging() -> None:
    """I-V curves: non-interacting vs Coulomb-blockaded junction."""
    cases = [
        {"qqC": 0.0,  "T": 298, "label": "Non-interacting ($q^2/C_\\Sigma = 0$)",  "color": "royalblue"},
        {"qqC": 0.5,  "T": 298, "label": "$q^2/C_\\Sigma = 0.5$ eV",              "color": "tomato"},
        {"qqC": 1.5,  "T": 298, "label": "$q^2/C_\\Sigma = 1.5$ eV",              "color": "seagreen"},
    ]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))

    for c in cases:
        print(f"  qqC={c['qqC']:.1f}, T={c['T']} K ...")
        V, I, G = iv_curve(E_F=E_F, E_H=E_H, gamma_S=0.1, gamma_D=0.1,
                           qqC=c["qqC"], T=c["T"], V_min=-4, V_max=4, n_V=80)
        ax1.plot(V, I * 1e6, color=c["color"], lw=1.6, label=c["label"])
        ax2.plot(V, G * 1e6, color=c["color"], lw=1.6, label=c["label"])

    for ax, ylabel, title in [
        (ax1, "$I$  [$\\mu$A]", "Current-voltage characteristic"),
        (ax2, "$dI/dV$  [$\\mu$A/V]", "Differential conductance"),
    ]:
        ax.axhline(0, color="k", lw=0.5, alpha=0.5)
        ax.axvline(0, color="k", lw=0.5, alpha=0.5)
        ax.set_xlabel("$V$ [V]", fontsize=11)
        ax.set_ylabel(ylabel, fontsize=11)
        ax.set_title(title, fontsize=11)
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.25)

    fig.suptitle(
        f"Single-level molecular junction  ($E_F = {E_F}$ eV,  "
        f"$\\epsilon_0 = {E_H}$ eV,  $\\Gamma = 0.2$ eV,  $T = 298$ K)",
        fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "transport_iv.png"), dpi=150)
    plt.close(fig)
    print("  saved transport_iv.png")


def plot_conductance_temperature() -> None:
    """Differential conductance for several temperatures (non-interacting)."""
    temps  = [77, 150, 298, 500]
    colors = ["navy", "royalblue", "tomato", "darkorange"]

    fig, ax = plt.subplots(figsize=(8, 4.5))

    for T, color in zip(temps, colors):
        print(f"  T={T} K ...")
        V, I, G = iv_curve(E_F=E_F, E_H=E_H, gamma_S=0.1, gamma_D=0.1,
                           qqC=0.0, T=T, V_min=-1.5, V_max=1.5, n_V=120)
        ax.plot(V, G * 1e6, color=color, lw=1.6, label=f"$T = {T}$ K")

    ax.axhline(0, color="k", lw=0.5, alpha=0.5)
    ax.axvline(0, color="k", lw=0.5, alpha=0.5)
    ax.set_xlabel("$V$ [V]", fontsize=11)
    ax.set_ylabel("$dI/dV$  [$\\mu$A/V]", fontsize=11)
    ax.set_title(
        f"Conductance peak vs temperature — non-interacting junction\n"
        f"($\\epsilon_0 - E_F = {E_H - E_F:.1f}$ eV,  $\\Gamma = 0.2$ eV)",
        fontsize=11)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "transport_conductance.png"), dpi=150)
    plt.close(fig)
    print("  saved transport_conductance.png")


def main():
    print("=== Molecular junction transport ===")
    plot_iv_charging()
    plot_conductance_temperature()
    print("figures written to", FIGDIR)


if __name__ == "__main__":
    main()
