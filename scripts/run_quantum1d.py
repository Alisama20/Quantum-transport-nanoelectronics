"""Solve the 1D Schrödinger equation and generate figures.

Two scenarios:
  1. Particle in a box (V=0) — compare numerical eigenvalues with analytical,
     plot the first four wave functions.
  2. Double-well potential — show tunnelling-split doublets.

Outputs:
    quantum1d_particle_in_box.png   — wave functions + energy levels
    quantum1d_double_well.png       — tunnelling doublets in a symmetric well
"""

import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from nanoelectronics import solve, particle_in_box_exact, make_double_well, make_barrier

FIGDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "figures"))
os.makedirs(FIGDIR, exist_ok=True)

N_GRID = 200
L      = 10.0   # nm


def plot_particle_in_box() -> None:
    n_show = 4
    V = np.zeros(N_GRID)
    E, psi, x = solve(V, L=L, n_states=n_show)

    # Analytical comparison
    E_ana = np.array([particle_in_box_exact(n + 1, L) for n in range(n_show)])
    print("  Particle in a box — first 4 levels:")
    for k in range(n_show):
        print(f"    n={k+1}: E_num={E[k]:.4f} eV  E_ana={E_ana[k]:.4f} eV  "
              f"err={abs(E[k]-E_ana[k])/E_ana[k]*100:.4f}%")

    colors = ["royalblue", "tomato", "seagreen", "darkorange"]
    scale  = 0.06   # scale wave functions for overlay

    fig, ax = plt.subplots(figsize=(8, 5.5))

    # Potential walls (infinite at 0 and L)
    ax.axvline(0,  color="k", lw=2)
    ax.axvline(L, color="k", lw=2)
    ax.set_xlim(-0.5, L + 0.5)

    for k in range(n_show):
        ax.axhline(E[k], color=colors[k], ls="--", lw=0.8, alpha=0.6)
        # Plot |ψ|² centred on energy level
        prob = psi[:, k]**2
        prob /= prob.max()
        ax.plot(x, E[k] + scale * prob, color=colors[k], lw=1.6,
                label=f"$n={k+1}$: $E_{k+1}$ = {E[k]:.3f} eV")
        # Fill
        ax.fill_between(x, E[k], E[k] + scale * prob,
                        color=colors[k], alpha=0.15)

    ax.set_xlabel("$x$ [nm]", fontsize=11)
    ax.set_ylabel("Energy  [eV]", fontsize=11)
    ax.set_title(f"Particle in a box — first {n_show} eigenstates  ($L = {L}$ nm, $N = {N_GRID}$)",
                 fontsize=11)
    ax.legend(fontsize=9, loc="upper left")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "quantum1d_particle_in_box.png"), dpi=150)
    plt.close(fig)
    print("  saved quantum1d_particle_in_box.png")


def plot_double_well() -> None:
    V0_list = [1.5, 3.0, 6.0]
    n_show  = 4

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5), sharey=False)

    for ax, V0 in zip(axes, V0_list):
        V = make_double_well(N_GRID, V0=V0, L=L)
        E, psi, x = solve(V, L=L, n_states=n_show)

        colors = ["royalblue", "tomato", "seagreen", "darkorange"]
        scale  = 0.12 * V0

        # Potential profile
        ax.plot(x, V, "k-", lw=1.2, label="$V(x)$")

        for k in range(n_show):
            ax.axhline(E[k], color=colors[k], ls="--", lw=0.8, alpha=0.6)
            prob = psi[:, k]**2
            prob /= prob.max()
            ax.plot(x, E[k] + scale * prob, color=colors[k], lw=1.4,
                    label=f"$E_{k+1}={E[k]:.2f}$ eV")
            ax.fill_between(x, E[k], E[k] + scale * prob,
                            color=colors[k], alpha=0.15)

        # Show tunnelling splitting of first doublet
        if n_show >= 2:
            delta = E[1] - E[0]
            ax.text(0.97, 0.05,
                    f"$\\Delta E_{{1,2}} = {delta*1e3:.1f}$ meV",
                    ha="right", va="bottom", transform=ax.transAxes,
                    fontsize=8, color="navy")

        ax.set_xlabel("$x$ [nm]", fontsize=10)
        ax.set_ylabel("Energy  [eV]", fontsize=10)
        ax.set_title(f"$V_0 = {V0}$ eV", fontsize=10)
        ax.legend(fontsize=7, loc="upper center")
        ax.grid(True, alpha=0.2)
        ax.set_xlim(0, L)

    fig.suptitle("Double-well potential — tunnelling doublets", fontsize=12)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "quantum1d_double_well.png"), dpi=150)
    plt.close(fig)
    print("  saved quantum1d_double_well.png")


def main():
    print("=== 1D Schrodinger equation ===")
    plot_particle_in_box()
    plot_double_well()
    print("figures written to", FIGDIR)


if __name__ == "__main__":
    main()
