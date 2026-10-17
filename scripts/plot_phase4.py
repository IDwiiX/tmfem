#!/usr/bin/env python3
"""Phase 4 figures for the README.

Inputs:
- written by `ctest --test-dir build -L phase4`:
  results/phase4/disk_profile.csv (P4-c/d), gaussian_{variance,profile}_t{3,6}.csv (P4-e),
  disk_convergence_t{3,6}.csv (P4-g), gaussian_convergence_t{3,6}.csv (P4-h)
- written by `./build/src/tmfem cases/beam_block/{copper,graphite}.json` (demo D4):
  results/phase4/beam_block/<material>/sweep.csv and s1/end_of_pulse_nodal.csv

Run:  python3 scripts/plot_phase4.py
Writes docs/figures/phase4_{thin_disk,gaussian,convergence,beam_block_maps,beam_block_sweep}.png.
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
RES = os.path.join(ROOT, "results", "phase4")
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

MATERIALS = (("copper", BLUE, "copper (OFE, annealed)"), ("graphite", ORANGE, "graphite (R4550)"))


def load(*parts):
    with open(os.path.join(RES, *parts)) as f:
        lines = [ln for ln in f if not ln.startswith("#")]
    names = lines[0].strip().split(",")
    data = np.array([[float(v) for v in ln.split(",")] for ln in lines[1:]])
    return {n: data[:, i] for i, n in enumerate(names)}


def slope(h, e):
    return np.polyfit(np.log(h[1:]), np.log(e[1:]), 1)[0]  # levels 1-4, as asserted


def thin_disk():
    d = load("disk_profile.csv")
    rr = np.linspace(0, 1, 200)
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    ax.plot(rr, -(1 - rr**2) / 4, color=INK2, lw=1.2, label="Timoshenko and Goodier (exact)")
    ax.plot(rr, -(1 - 3 * rr**2) / 4, color=INK2, lw=1.2)
    ax.plot(d["r_over_b"], d["sr_fe"], "o", ms=4, color=BLUE, label=r"tmfem T6, $\sigma_r$")
    ax.plot(d["r_over_b"], d["st_fe"], "s", ms=4, color=ORANGE, label=r"tmfem T6, $\sigma_\theta$")
    ax.axhline(0, color=INK2, lw=0.6)
    ax.set_xlabel("r / b")
    ax.set_ylabel(r"stress / (E $\alpha$ $\Delta T$)")
    ax.set_title(r"P4-c/d: thin disk, T = $\Delta T$(1 - r$^2$/b$^2$), plane stress")
    ax.legend(loc="upper left")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "phase4_thin_disk.png"), dpi=150)


def gaussian():
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.0))
    for t, color, marker, ms, fill in ((3, BLUE, "o", 4, BLUE), (6, ORANGE, "s", 9, "none")):  # T6 hollow, over T3
        p = load(f"gaussian_profile_t{t}.csv")
        o = np.argsort(p["x_m"])
        sel = o[:: max(1, len(o) // 40)]
        a2.plot(p["x_m"][sel] * 1e3, p["T_fe_K"][sel], marker, ms=ms, color=color, mfc=fill, mew=1.5, linestyle="none",
                label=f"tmfem T{t}")
        v = load(f"gaussian_variance_t{t}.csv")
        a1.plot(v["t_s"] * 1e3, v["variance_m2"] * 1e6, marker, ms=ms + 2, color=color, mfc=fill, mew=1.5,
                linestyle="none", label=f"tmfem T{t}")
    p = load("gaussian_profile_t6.csv")
    o = np.argsort(p["x_m"])
    a2.plot(p["x_m"][o] * 1e3, p["T_exact_K"][o], color=INK2, lw=1.2, label="exact")
    x = np.linspace(0, 12, 200)
    a2.plot(x, 100 * np.exp(-x**2 / 2), color=INK2, lw=1.0, linestyle=":", label="initial (t = 0)")
    tt = np.linspace(0, v["t_s"].max(), 50)
    a1.plot(tt * 1e3, (1e-6 + 2 * (400 / 3.44e6) * tt) * 1e6, color=INK2, lw=1.2, label=r"$s_0^2 + 2at$")
    a1.set_xlabel("t [ms]")
    a1.set_ylabel(r"variance per coordinate [mm$^2$]")
    a1.set_title("P4-e: variance of the spreading Gaussian")
    a1.legend()
    a2.set_xlabel("x [mm] (y = 0)")
    a2.set_ylabel("T [K]")
    a2.set_title("P4-e: profile at t_end (s = 2 s0)")
    a2.set_xlim(0, 8)
    a2.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "phase4_gaussian.png"), dpi=150)


def convergence():
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.2))
    for t, color in ((3, BLUE), (6, ORANGE)):
        c = load(f"disk_convergence_t{t}.csv")
        h = c["h_m"] * 1e3
        a1.loglog(h, c["L2_rel"], "-o", color=color, ms=4, lw=1.5, label=f"T{t} L2 (slope {slope(h, c['L2_rel']):.2f})")
        a1.loglog(h, c["energy_rel"], "--s", color=color, ms=4, lw=1.5,
                  label=f"T{t} energy (slope {slope(h, c['energy_rel']):.2f})")
        g = load(f"gaussian_convergence_t{t}.csv")
        h = g["h_m"] * 1e3
        a2.loglog(h, g["L2_rel"], "-o", color=color, ms=4, lw=1.5, label=f"T{t} L2 (slope {slope(h, g['L2_rel']):.2f})")
    a1.set_title("P4-g: thin disk, thermal stress")
    a2.set_title("P4-h: Gaussian heat spreading (CN, dt ~ h^((p+1)/2))")
    for ax in (a1, a2):
        ax.set_xlabel("h [mm]")
        ax.set_ylabel("relative error")
        ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "phase4_convergence.png"), dpi=150)


def beam_maps():
    fig, axes = plt.subplots(2, 2, figsize=(10, 8.2))
    for row, (mat, _, label) in enumerate(MATERIALS):
        d = load("beam_block", mat, "s1", "end_of_pulse_nodal.csv")
        x, y = d["x_m"] * 1e3, d["y_m"] * 1e3
        tri = mtri.Triangulation(x, y)
        keep = (x < 4.5) & (y < 4.5)
        tri.set_mask(~np.all(keep[tri.triangles], axis=1))
        for col, (key, cmap, unit) in enumerate((("T_K", "Oranges", "T [K]"), ("von_mises_MPa", "Blues", "von Mises [MPa]"))):
            ax = axes[row, col]
            vmin = 293.15 if key == "T_K" else 0.0  # T0; the field has no undershoot (min - T0 > -1e-9 K)
            cs = ax.tricontourf(tri, d[key], levels=np.linspace(vmin, d[key][keep].max(), 21), cmap=cmap)
            fig.colorbar(cs, ax=ax, label=unit)
            ax.set_aspect("equal")
            ax.set_xlim(0, 4)
            ax.set_ylim(0, 4)
            ax.grid(False)
            ax.set_xlabel("x [mm]")
            ax.set_ylabel("y [mm]")
            ax.set_title(f"{label}, s = 1 mm: {unit.split(' [')[0]}")
    fig.suptitle("Demo D4 (not a validation): end of a 10 us, 1 kJ/m pulse, plane strain, quarter model", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "phase4_beam_block_maps.png"), dpi=150)


def beam_sweep():
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.2))
    for mat, color, label in MATERIALS:
        s = load("beam_block", mat, "sweep.csv")
        a1.loglog(s["s_mm"], s["vm_centre_MPa"], "o", ms=7, color=color, label=f"{label}, tmfem")
        a1.loglog(s["s_mm"], s["vm_centre_estimate_MPa"], "-", lw=1.2, color=color, label=f"{label}, infinite-plane estimate")
        name = "von Mises / yield" if mat == "copper" else "max principal / flexural (optimistic)"
        a2.loglog(s["s_mm"], s["strength_ratio"], "-o", ms=6, color=color, label=f"{label}: {name}")
        if mat == "graphite":
            a2.loglog(s["s_mm"], s["compressive_ratio"], "--s", ms=6, color=color, label=f"{label}: |min principal| / compressive")
    a2.axhline(1.0, color=INK2, lw=1.0)
    a2.annotate("strength limit", (0.52, 1.08), color=INK2, fontsize=9)
    a1.set_xlabel("beam rms size s [mm]")
    a1.set_ylabel("von Mises at the centre [MPa]")
    a1.set_title("Centre stress at the end of the pulse")
    a2.set_xlabel("beam rms size s [mm]")
    a2.set_ylabel("peak stress / strength")
    a2.set_title("Strength check (elastic; above 1 copper would yield)")
    a1.legend(fontsize=8)
    a2.legend(fontsize=8)
    fig.suptitle("Demo D4 (not a validation): equal deposit E_l = 1 kJ/m in both materials", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "phase4_beam_block_sweep.png"), dpi=150)


if __name__ == "__main__":
    os.makedirs(FIG, exist_ok=True)
    thin_disk()
    gaussian()
    convergence()
    beam_maps()
    beam_sweep()
    print("wrote docs/figures/phase4_*.png")
