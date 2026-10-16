#!/usr/bin/env python3
"""Phase 1 figures for the README, from the CSVs the phase1 tests write to results/phase1/.

Run after `ctest --test-dir build -L phase1`:  python3 scripts/plot_phase1.py
Writes docs/figures/phase1_{erf_profile,thin_film_decay,steel_plate_inner_face,convergence}.png.
"""
import os
import sys

# The numpy in ~/.local (2.x) breaks the system matplotlib (built for numpy 1.x): restart without
# user site-packages so the system numpy and matplotlib are used together.
if not sys.flags.no_user_site:
    os.execv(sys.executable, [sys.executable, "-s"] + sys.argv)

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker  # noqa: E402
import numpy as np  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results", "phase1")
FIG = os.path.join(ROOT, "docs", "figures")

# Reference categorical palette, slots 1-3 in fixed order; text and grid in neutral inks.
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
    "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 2.0,
    "legend.frameon": False, "font.size": 10, "axes.titlesize": 11, "axes.titleweight": "bold",
})


def load(name):
    path = os.path.join(RES, name)
    if not os.path.exists(path):
        sys.exit(f"missing {path}: run `ctest --test-dir build -L phase1` first")
    with open(path) as f:
        lines = [line for line in f if not line.startswith("#")]  # names=True would read a comment as header
    return np.genfromtxt(lines, delimiter=",", names=True)


def save(fig, name):
    os.makedirs(FIG, exist_ok=True)
    out = os.path.join(FIG, name)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


def erf_profile():
    d = load("carburizing_profile.csv")
    x_mm = d["x_m"] * 1e3
    keep = x_mm <= 2.0
    fe_marks = keep & (np.arange(len(x_mm)) % 10 == 0)
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    for col, label, color in [("t3600s", "t = 1 h", BLUE), ("t25370s", "t = 7.05 h", ORANGE)]:
        ax.plot(x_mm[keep], d["C_exact_" + col][keep], color=color, label=f"erf solution, {label}")
        ax.plot(x_mm[fe_marks], d["C_fe_" + col][fe_marks], "o", ms=4.5, mfc="none", mew=1.2, color=color,
                label=f"FE (lumped BE), {label}")
    ax.plot([0.5], [0.8], "s", color=INK, ms=6)
    ax.annotate("0.8 wt% at 0.5 mm (P1-b)", (0.5, 0.8), xytext=(0.75, 0.95), color=INK2,
                arrowprops=dict(arrowstyle="-", color=INK2, lw=0.8))
    ax.set(xlabel="depth x [mm]", ylabel="carbon C [wt%]",
           title="Carburizing: FE vs erf solution (D = 1.6e-11 m²/s)")
    ax.legend(loc="upper right", fontsize=8.5)
    save(fig, "phase1_erf_profile.png")


def thin_film_decay():
    d = load("thin_film_decay.csv")
    t = d["t_s"]
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(9.6, 3.8))
    ax.loglog(t, d["C0_exact_per_m3"], color=BLUE, label="analytic M / sqrt(pi D t)")
    step = max(1, len(t) // 40)
    ax.loglog(t[::step], d["C0_fe_per_m3"][::step], "o", ms=4.5, mfc="none", mew=1.2, color=ORANGE,
              label="FE (BE, dt = 0.01 t)")
    ax.set(xlabel="time t [s]", ylabel="surface concentration C(0, t) [atoms/m³]",
           title="Thin film: peak decays as t^-1/2")
    ax.legend(fontsize=8.5)
    drift = np.abs(d["dose_rel_drift"])
    ax2.loglog(t[1:], np.maximum(drift[1:], 1e-16), color=BLUE)  # 1e-16 floor: exact zeros on a log axis
    ax2.axhline(1e-10, color=INK2, lw=1, ls="--")
    ax2.text(t[1], 1e-10, " P1-h tolerance 1e-10", va="bottom", color=INK2, fontsize=8.5)
    ax2.set(xlabel="time t [s]", ylabel="|dose drift| (relative)", title="Dose conservation")
    fig.tight_layout()
    save(fig, "phase1_thin_film_decay.png")


def steel_plate():
    d = load("steel_plate_synthetic.csv")
    t_h = d["t_s"] / 3600.0
    c = lambda T: T - 273.15  # noqa: E731
    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(7.2, 5.6), sharex=True, gridspec_kw={"height_ratios": [2, 1]})
    ax.plot(t_h, c(d["T_sa_K"]), color=INK2, lw=1.2, ls="--", label="sol-air T_sa (synthetic)")
    ax.plot(t_h, c(d["T_lumped_K"]), color=BLUE, label="lumped ODE (RK4, 1 s)")
    ax.plot(t_h[::6], c(d["T_inner_fe_K"])[::6], "o", ms=4, mfc="none", mew=1.2, color=ORANGE,
            label="1D FE inner face (BE, 300 s)")
    ax.set(ylabel="temperature [°C]", title="Steel plate 5 mm: inner face, FE vs lumped model")
    ax.legend(loc="lower right", fontsize=8.5)
    ax2.plot(t_h, d["T_inner_fe_K"] - d["T_lumped_K"], color=ORANGE, label="FE inner − lumped (P1-j, tol 0.05 K)")
    ax2.plot(t_h, d["T_outer_fe_K"] - d["T_inner_fe_K"], color=AQUA, label="FE outer − inner (P1-k, 0.0255 K)")
    ax2.axhline(0, color=INK2, lw=0.8)
    ax2.set(xlabel="time [h]", ylabel="difference [K]")
    ax2.set_xticks(np.arange(0, 73, 12))
    ax2.legend(loc="lower right", fontsize=8.5)
    fig.tight_layout()
    save(fig, "phase1_steel_plate_inner_face.png")


def slope_guide(ax, x, y0, order, label):
    x = np.array([x[0], x[-1]])
    ax.loglog(x, y0 * (x / x[0]) ** order, color=INK2, lw=1, ls=":")
    ax.text(x[-1], y0 * (x[-1] / x[0]) ** order, f"  slope {label}", color=INK2, fontsize=8.5, va="center")


def convergence():
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(9.6, 3.9))
    for name, label, color, marker, order in [("convergence_time_theta1.csv", "backward Euler", BLUE, "o", 1),
                                               ("convergence_time_theta0.5.csv", "Crank-Nicolson", ORANGE, "s", 2)]:
        d = load(name)
        ax.loglog(d["dt_s"], d["max_nodal_error"], marker=marker, ms=6, color=color, label=label)
        slope_guide(ax, d["dt_s"], d["max_nodal_error"][0] * 1.6, order, str(order))
    ax.set(xlabel="time step dt", ylabel="max nodal time error", title="Time order (P1-m)")
    ax.legend(fontsize=8.5, loc="lower right")
    for name, label, color, marker in [("convergence_space_consistent.csv", "consistent mass", BLUE, "o"),
                                       ("convergence_space_lumped.csv", "lumped mass", ORANGE, "s")]:
        d = load(name)
        ax2.loglog(d["h_m"], d["l2_error"], marker=marker, ms=6, color=color, label=label)
    d = load("convergence_space_consistent.csv")
    slope_guide(ax2, d["h_m"], d["l2_error"][0] * 1.6, 2, "2")
    ax2.set(xlabel="element size h", ylabel="L2 error", title="Space order (P1-m), CN dt = 1e-4")
    for a in (ax, ax2):
        a.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax2.legend(fontsize=8.5, loc="lower right")
    fig.tight_layout()
    save(fig, "phase1_convergence.png")


def main():
    erf_profile()
    thin_film_decay()
    steel_plate()
    convergence()


if __name__ == "__main__":
    main()
