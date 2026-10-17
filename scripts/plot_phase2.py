#!/usr/bin/env python3
"""Phase 2 figures for the README.

Inputs:
- results/kirsch/kirsch_t6_nodal.csv, written by `./build/src/tmfem cases/kirsch/kirsch_t6.json`
- results/phase2/convergence_t{3,6}.csv, written by `ctest --test-dir build -L P2-f`

Run:  python3 scripts/plot_phase2.py
Writes docs/figures/phase2_{kirsch_von_mises,kirsch_profiles,convergence}.png.
The nodal stresses are element averages at the nodes (plots only); the validated point values come
from the tests.
"""
import os
import sys

# Same environment workaround as plot_phase1.py: the numpy in ~/.local breaks the system matplotlib.
if not sys.flags.no_user_site:
    os.execv(sys.executable, [sys.executable, "-s"] + sys.argv)

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.tri as mtri  # noqa: E402
import numpy as np  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")
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

S_MPA, A_M = 140.0, 0.010  # cases/kirsch/kirsch_t6.json


def load(path):
    with open(path) as f:
        lines = [ln for ln in f if not ln.startswith("#")]
    names = lines[0].strip().split(",")
    data = np.array([[float(v) for v in ln.split(",")] for ln in lines[1:]])
    return {n: data[:, i] for i, n in enumerate(names)}


def kirsch_polar(r, t):
    """Kirsch stresses in MPa, theta from the load axis y (docs/learning/phase2.md section 9)."""
    q2 = (A_M / r) ** 2
    q4 = q2 * q2
    srr = S_MPA / 2 * (1 - q2) + S_MPA / 2 * (1 - 4 * q2 + 3 * q4) * np.cos(2 * t)
    stt = S_MPA / 2 * (1 + q2) - S_MPA / 2 * (1 + 3 * q4) * np.cos(2 * t)
    return srr, stt


def von_mises_map(d):
    x, y = d["x_m"] / A_M, d["y_m"] / A_M
    tri = mtri.Triangulation(x, y)
    cx, cy = x[tri.triangles].mean(axis=1), y[tri.triangles].mean(axis=1)
    tri.set_mask(np.hypot(cx, cy) < 1.0)  # Delaunay fills the hole; remove it
    fig, ax = plt.subplots(figsize=(5.6, 4.6))
    cs = ax.tricontourf(tri, d["von_mises_MPa"], levels=24, cmap="viridis")
    fig.colorbar(cs, ax=ax, label="von Mises stress [MPa]")
    ax.set_aspect("equal")
    ax.grid(False)
    ax.set_xlabel("x / a")
    ax.set_ylabel("y / a  (load axis)")
    ax.set_title("Kirsch plate, S = 140 MPa along y, finest T6 mesh (level 4)")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "phase2_kirsch_von_mises.png"), dpi=150)
    plt.close(fig)


def profiles(d):
    x, y = d["x_m"], d["y_m"]
    r = np.hypot(x, y)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 3.8))
    on_axis = np.abs(y) < 1e-12  # y = 0: theta = 90 deg from the load axis, sigma_tt = sigma_yy
    order = np.argsort(x[on_axis])
    xs = x[on_axis][order]
    rr = np.linspace(A_M, 5 * A_M, 400)
    a1.plot(rr / A_M, kirsch_polar(rr, np.pi / 2)[1], color=INK2, lw=1.2, label="Kirsch (exact)")
    a1.plot(xs / A_M, d["syy_MPa"][on_axis][order], "o", ms=3, color=BLUE, label="tmfem T6 (nodal average)")
    a1.set_xlabel("x / a  on y = 0")
    a1.set_ylabel(r"$\sigma_{yy}$ [MPa]")
    a1.set_title("Stress decay from the hole edge")
    a1.legend()

    on_hole = np.abs(r - A_M) < 1e-9 * A_M
    t = np.arctan2(x[on_hole], y[on_hole])
    sxx, syy, sxy = d["sxx_MPa"][on_hole], d["syy_MPa"][on_hole], d["sxy_MPa"][on_hole]
    et = np.stack([np.cos(t), -np.sin(t)])
    stt = et[0] ** 2 * sxx + et[1] ** 2 * syy + 2 * et[0] * et[1] * sxy
    o = np.argsort(t)
    tt = np.linspace(0, np.pi / 2, 300)
    a2.plot(np.degrees(tt), kirsch_polar(A_M, tt)[1], color=INK2, lw=1.2, label="Kirsch (exact)")
    a2.plot(np.degrees(t[o]), stt[o], "o", ms=3, color=ORANGE, label="tmfem T6 (nodal average)")
    a2.set_xlabel(r"$\theta$ from the load axis [deg]")
    a2.set_ylabel(r"$\sigma_{\theta\theta}$ at r = a [MPa]")
    a2.set_title("Hoop stress around the hole (Kt = 3 at 90 deg)")
    a2.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "phase2_kirsch_profiles.png"), dpi=150)
    plt.close(fig)


def convergence():
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 3.8))
    for t, color in ((3, BLUE), (6, ORANGE)):
        c = load(os.path.join(RES, "phase2", f"convergence_t{t}.csv"))
        for plane, style, name in ((0, "-o", "plane stress"), (1, "--s", "plane strain")):
            m = c["plane_strain"] == plane
            h = c["h_m"][m] * 1e3
            fit = c["in_fit"][m] == 1  # level 0 is pre-asymptotic and is not in the asserted fit
            for ax, key in ((a1, "L2_rel"), (a2, "energy_rel")):
                slope = np.polyfit(np.log(h[fit]), np.log(c[key][m][fit]), 1)[0]
                ax.loglog(h, c[key][m], style, color=color, ms=4, lw=1.5, label=f"T{t} {name} (slope {slope:.2f})")
    for ax, title, orders in ((a1, "Displacement, relative L2 error", (2, 3)),
                              (a2, "Relative energy-norm error", (1, 2))):
        ax.set_xlabel("h = sqrt(area / cells) [mm]  (slopes fitted on levels 1-4)")
        ax.set_title(title + f" (theory: T3 {orders[0]}, T6 {orders[1]})")
        ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "phase2_convergence.png"), dpi=150)
    plt.close(fig)


def main():
    os.makedirs(FIG, exist_ok=True)
    nodal = load(os.path.join(RES, "kirsch", "kirsch_t6_nodal.csv"))
    von_mises_map(nodal)
    profiles(nodal)
    convergence()
    print("wrote docs/figures/phase2_{kirsch_von_mises,kirsch_profiles,convergence}.png")


if __name__ == "__main__":
    main()
