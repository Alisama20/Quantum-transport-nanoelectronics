"""Ballistic MOSFET: I-V family and comparison with semiclassical model.

Outputs:
    mosfet_iv_family.png    — I_DS vs V_DS for several V_GS values (ballistic)
    mosfet_comparison.png   — ballistic vs long-channel MOSFET at V_GS = 0.5 V
"""

import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import numpy as np
import matplotlib.pyplot as plt

from nanoelectronics import iv_ballistic, iv_semiclassical

FIGDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "figures"))
os.makedirs(FIGDIR, exist_ok=True)

# Device parameters (match the C++ source)
PARAMS = dict(
    E_F=-5.0, E_C=-4.7,
    L=40e-9, W=120e-9,
    C_G=0.1e-15, C_S=0.0, C_D=0.0,
    T=5.0,
    V_min=0.0, V_max=0.5, n_V=51,
)


def plot_iv_family() -> None:
    """Family of I_DS vs V_DS curves for different V_GS."""
    vgs_list = [0.3, 0.4, 0.5, 0.6, 0.7]
    colors   = plt.cm.plasma(np.linspace(0.1, 0.85, len(vgs_list)))

    fig, ax = plt.subplots(figsize=(8, 5))

    for V_GS, color in zip(vgs_list, colors):
        print(f"  V_GS = {V_GS:.1f} V ...")
        V, I = iv_ballistic(**PARAMS, V_GS=V_GS)
        ax.plot(V, I, color=color, lw=1.8,
                label=f"$V_{{GS}} = {V_GS:.1f}$ V")

    ax.set_xlabel("$V_{DS}$  [V]", fontsize=11)
    ax.set_ylabel("$I_{DS}$  [$\\mu$A]", fontsize=11)
    ax.set_title(
        f"Ballistic MOSFET — $I_{{DS}}$–$V_{{DS}}$ family  ($T = {PARAMS['T']:.0f}$ K, "
        f"$L = {round(PARAMS['L']*1e9)}$ nm, $W = {round(PARAMS['W']*1e9)}$ nm)",
        fontsize=11,
    )
    ax.legend(fontsize=9, title="$V_{GS}$")
    ax.set_xlim(0, 0.5)
    ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "mosfet_iv_family.png"), dpi=150)
    plt.close(fig)
    print("  saved mosfet_iv_family.png")


def plot_comparison() -> None:
    """Ballistic vs semiclassical I_DS vs V_DS at V_GS = 0.5 V."""
    V_GS = 0.5
    V_T  = PARAMS["E_C"] - PARAMS["E_F"]   # = -4.7 - (-5.0) = +0.3 V

    # Room-temperature ballistic for fair comparison with the classical model
    params_rt = {**PARAMS, "T": 298}
    print(f"  Ballistic model (T = 298 K) ...")
    V_b, I_b = iv_ballistic(**params_rt, V_GS=V_GS)

    print(f"  Ballistic model (T = {PARAMS['T']:.0f} K) ...")
    V_b5, I_b5 = iv_ballistic(**PARAMS, V_GS=V_GS)

    print("  Semiclassical model ...")
    V_c, I_c = iv_semiclassical(
        V_GS=V_GS, V_T=V_T, mu_n=300e-4,
        L=PARAMS["L"], W=PARAMS["W"], C_G=PARAMS["C_G"],
        V_min=PARAMS["V_min"], V_max=PARAMS["V_max"], n_V=PARAMS["n_V"],
    )

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(V_b,  I_b,  "bo-", ms=4, lw=1.4, label="Ballistic — $T = 298$ K")
    ax.plot(V_b5, I_b5, "cs-", ms=4, lw=1.4, label=f"Ballistic — $T = {PARAMS['T']:.0f}$ K")
    ax.plot(V_c,  I_c,  "r-",  lw=2.0,       label="Semiclassical (drift-diffusion, $T = 298$ K)")

    ax.set_xlabel("$V_{DS}$  [V]", fontsize=11)
    ax.set_ylabel("$I_{DS}$  [$\\mu$A]", fontsize=11)
    ax.set_title(
        f"Ballistic vs semiclassical MOSFET  ($V_{{GS}} = {V_GS}$ V,  $V_T = {V_T:.1f}$ V)\n"
        f"$L = {int(PARAMS['L']*1e9)}$ nm,  $W = {int(PARAMS['W']*1e9)}$ nm",
        fontsize=11,
    )
    ax.legend(fontsize=10)
    ax.set_xlim(0, 0.5)
    ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "mosfet_comparison.png"), dpi=150)
    plt.close(fig)
    print("  saved mosfet_comparison.png")


def main():
    print("=== Ballistic MOSFET ===")
    plot_iv_family()
    plot_comparison()
    print("figures written to", FIGDIR)


if __name__ == "__main__":
    main()
