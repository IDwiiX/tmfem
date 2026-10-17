#!/usr/bin/env python3
"""Phase 3 figures for the README.

Inputs (all written by `ctest --test-dir build -L phase3`):
- results/phase3/tensile_fit_curves.csv, tensile_fit_data.csv (P3-a)
- results/phase3/material_point_tabulated.csv (P3-b)
- results/phase3/newton_history_*.csv (P3-d)
- results/phase3/kirsch_plastic_plane_{stress,strain}_l3_gauss.csv (P3-f)

Run:  python3 scripts/plot_phase3.py
Writes docs/figures/phase3_{tensile_fit,newton_history,plastic_zone}.png.
"""
import os
import sys

# Same environment workaround as plot_phase1.py: the numpy in ~/.local breaks the system matplotlib.
if not sys.flags.no_user_site:
    os.execv(sys.executable, [sys.executable, "-s"] + sys.argv)

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results", "phase3")
FIG = os.path.join(ROOT, "docs", "figures")

BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
    "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 2.0,
    "legend.frameon": False, "font.size": 10, "axes.titlesize": 11, "axes.titleweight": "bold",
})

E_MPA = 207423.0   # material_tensile.json
A_M = 0.010        # hole radius, cases/kirsch/kirsch.geo


def load(name):
    with open(os.path.join(RES, name)) as f:
        lines = [ln for ln in f if not ln.startswith("#")]
    names = lines[0].strip().split(",")
    data = np.array([[float(v) for v in ln.split(",")] for ln in lines[1:]])
    return {n: data[:, i] for i, n in enumerate(names)}


def tensile_fit():
    c, d, mp = load("tensile_fit_curves.csv"), load("tensile_fit_data.csv"), load("material_point_tabulated.csv")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))
    a = c["alpha"]
    ax1.fill_between(a, c["sigma_lsq_MPa"] - c["band95_halfwidth_MPa"], c["sigma_lsq_MPa"] + c["band95_halfwidth_MPa"],
                     color=BLUE, alpha=0.18, linewidth=0, label="95% confidence band (LSQ)")
    ax1.plot(a, c["sigma_lsq_MPa"], color=BLUE, label="Ludwik, least squares (3 free)")
    ax1.plot(a, c["sigma_source_MPa"], color=ORANGE, linestyle="--", label="Ludwik, source (310.86, 559.27, 0.4077)")
    ax1.plot(a, c["sigma_loglog_MPa"], color=INK2, linestyle=":", label="Ludwik, log-log (sigma_y fixed)")
    ax1.plot(d["alpha"], d["sigma_MPa"], "o", color=INK, markersize=6, label="true-curve data")
    ax1.set(xlabel="equivalent plastic strain alpha = eps_true - sigma/E", ylabel="true stress (MPa)",
            title="P3-a  Ludwik fits to the tensile curve")
    ax1.legend(loc="lower right", fontsize=8.5)

    # Uniaxial curves in total strain: eps = alpha + sigma / E.
    ax2.plot(a + c["sigma_lsq_MPa"] / E_MPA, c["sigma_lsq_MPa"], color=BLUE, label="Ludwik LSQ")
    ax2.plot(a + c["sigma_source_MPa"] / E_MPA, c["sigma_source_MPa"], color=ORANGE, linestyle="--",
             label="Ludwik source")
    ax2.plot(mp["eps_yy"], mp["sigma_yy_MPa"], "s", color=AQUA, markersize=9, markerfacecolor="none",
             markeredgewidth=2, label="material point, tabulated (P3-b)")
    eps_data = d["alpha"] + d["sigma_MPa"] / E_MPA
    ax2.plot(eps_data, d["sigma_MPa"], "o", color=INK, markersize=4, label="true-curve data")
    ax2.annotate("457.24 MPa: eps_yy = 0.03158\neps_xx = eps_zz = -0.01535", xy=(0.031581, 457.24),
                 xytext=(0.06, 400), fontsize=8.5, color=INK2, arrowprops={"arrowstyle": "-", "color": INK2})
    ax2.set(xlabel="true strain eps_yy", ylabel="true stress (MPa)", title="Uniaxial material point",
            xlim=(0, 0.185), ylim=(280, 600))
    ax2.legend(loc="lower right", fontsize=8.5)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "phase3_tensile_fit.png"), dpi=150)


def newton_history():
    series = [("bar_ludwik", "bar, plane strain (asserted)", BLUE, "o", "-"),
              ("bar_plane_stress_ludwik", "bar, plane stress (asserted)", AQUA, "s", "-"),
              ("kirsch_ludwik", "hole plate, plane strain (diagnostic)", ORANGE, "^", "--")]
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    for name, label, color, marker, ls in series:
        h = load(f"newton_history_{name}.csv")
        # Iterate 0 is the start of the step; for the displacement-controlled bars it is the old equilibrium
        # (round-off), so every series starts at the first correction.
        k = h["iteration"] >= 1
        ax.semilogy(h["iteration"][k], h["residual_over_Fref_final"][k], marker=marker, color=color, linestyle=ls,
                    markersize=8, label=label)
    ax.axhline(1e-11, color=INK2, linewidth=1, linestyle=":")
    ax.text(4.6, 2.0e-11, "rate floor 1e-11 F_ref", color=INK2, fontsize=8.5)
    ax.set(xlabel="Newton correction in the plastic step", ylabel="||R_free|| / F_ref",
           title="P3-d  Newton residuals, Ludwik hardening")
    ax.legend(loc="upper right", fontsize=8.5)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "phase3_newton_history.png"), dpi=150)


def plastic_zone():
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), sharey=True)
    vmax = max(load(f"kirsch_plastic_plane_{p}_l3_gauss.csv")["alpha"].max() for p in ("stress", "strain"))
    for ax, plane in zip(axes, ("stress", "strain")):
        g = load(f"kirsch_plastic_plane_{plane}_l3_gauss.csv")
        x, y, al = g["x_m"] / A_M, g["y_m"] / A_M, g["alpha"]
        win = (x < 1.45) & (y < 0.55)
        el = win & (al <= 0)
        pl = win & (al > 0)
        ax.scatter(x[el], y[el], s=0.6, color=GRID, linewidths=0)
        sc = ax.scatter(x[pl], y[pl], c=al[pl], s=4, cmap="Blues", vmin=0, vmax=vmax, linewidths=0)
        t = np.linspace(0, np.pi / 2, 200)
        ax.plot(np.cos(t), np.sin(t), color=INK, linewidth=1.5)
        r_p = np.hypot(x[pl], y[pl]).max()
        ax.set(xlim=(0.6, 1.45), ylim=(0, 0.55), aspect="equal", xlabel="x / a", title=f"plane {plane}")
        ax.text(1.02, 0.47, f"zone reaches r = {r_p:.3f} a\nmax alpha = {al.max():.2e}", fontsize=8.5, color=INK2)
    axes[0].set_ylabel("y / a")
    fig.suptitle("P3-f  plastic zone at the hole edge, S = 140 MPa along y (gray: elastic points)",
                 fontweight="bold", fontsize=11)
    cb = fig.colorbar(sc, ax=axes, shrink=0.85, pad=0.02)
    cb.set_label("equivalent plastic strain alpha")
    fig.savefig(os.path.join(FIG, "phase3_plastic_zone.png"), dpi=150, bbox_inches="tight")


if __name__ == "__main__":
    os.makedirs(FIG, exist_ok=True)
    tensile_fit()
    newton_history()
    plastic_zone()
    print("wrote docs/figures/phase3_{tensile_fit,newton_history,plastic_zone}.png")
