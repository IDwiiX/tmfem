#!/usr/bin/env python3
"""Recompute every reference value in docs/REFERENCE_VALUES.md from first principles.

Run:  python3 references/verify_reference_values.py
Writes references/data/kirsch_sigma140MPa.csv and prints a check report.
Needs only Python 3 and numpy.
"""
from math import erf, sqrt, pi, exp, log, cos, sin, radians
import csv
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
R = 8.314


def erfinv(y, lo=0.0, hi=6.0):
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if erf(mid) < y:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def kirsch(S, r_over_a, theta_deg):
    """Infinite plate, hole radius a, remote stress S along y; theta from the load axis (y)."""
    a2 = 1.0 / r_over_a**2
    a4 = a2 * a2
    t = radians(theta_deg)
    srr = S / 2 * (1 - a2) + S / 2 * (1 - 4 * a2 + 3 * a4) * cos(2 * t)
    stt = S / 2 * (1 + a2) - S / 2 * (1 + 3 * a4) * cos(2 * t)
    srt = -S / 2 * (1 + 2 * a2 - 3 * a4) * sin(2 * t)
    return srr, stt, srt


def read_csv(name):
    with open(os.path.join(DATA, name)) as f:
        rows = [r for r in csv.reader(f) if r and not r[0].startswith("#")]
    return np.array([[float(v) for v in r] for r in rows[1:]])


def main():
    print("== Phase 1: diffusion and heat ==")
    J = -3e-11 * (0.8 - 1.2) / (10e-3 - 5e-3)
    print(f"ex1 flux J = {J:.3e} kg/m2s   (doc 2.4e-12: off by 1000)")
    z = erfinv(1 - (0.8 - 0.25) / (1.2 - 0.25))
    t = (0.5e-3 / (2 * z)) ** 2 / 1.6e-11
    print(f"ex2 time = {t:.0f} s = {t/3600:.3f} h   (doc 7.025 h)")
    print(f"ex3 t(500C) = {5.3e-13 * 10 / 4.8e-14:.1f} h")
    D1 = 1.8e-5 * exp(-152000 / (R * 1000))
    D2 = 1.8e-5 * exp(-152000 / (R * 800))
    print(f"ex4 D1 = {D1:.3e}, D2 = {D2:.3e}, x2 = {15.5 * sqrt(D2 / D1):.3f} mm")
    M = 2 * 3e20 * sqrt(4e-13 * 360 / pi)
    Cs2 = M / sqrt(pi * 3e-12 * 360)
    T1, T2, T3 = 1373.15, 1473.15, 1573.15
    Q = R * log(3e-12 / 4e-13) / (1 / T1 - 1 / T2)
    D0 = 4e-13 * exp(Q / (R * T1))
    D3 = D0 * exp(-Q / (R * T3))
    print(f"ex5 M = {M:.3e} at/cm2, Cs(1200C) = {Cs2:.3e} at/cm3 (doc 2.8e19)")
    print(f"ex5 Q = {Q/1000:.1f} kJ/mol, D0 = {D0:.3f} cm2/s, D(1300C) = {D3:.3e} cm2/s")

    k, cp, rho, L, he, hi = 60.0, 434.0, 7850.0, 0.005, 22.0, 9.0
    print(f"steel plate: Bi = {he*L/k:.5f}, tau = {rho*cp*L/(he+hi):.1f} s, "
          f"L^2/alpha = {L*L*rho*cp/k:.2f} s, T_in ~ {he/(he+hi):.3f} T_sa + {hi/(he+hi)*26:.2f} C")

    print("\n== Phase 2: Kirsch, S = 140 MPa ==")
    thetas = [0, 25, 50, 90, 125, 150, 200, 240, 310, 360]
    out = os.path.join(DATA, "kirsch_sigma140MPa.csv")
    with open(out, "w", newline="") as f:
        w = csv.writer(f)
        f.write("# Recomputed Kirsch solution, remote 140 MPa along y, theta from the load axis (y).\n")
        f.write("# Source tables at r/a = 1, 2, 4 match; the source table at r/a = 1.5 is wrong.\n")
        w.writerow(["r_over_a", "theta_deg", "sigma_rr_MPa", "sigma_tt_MPa", "sigma_rt_MPa"])
        for r in [1, 1.5, 2, 4]:
            for th in thetas:
                w.writerow([r, th] + [f"{v:.4f}" for v in kirsch(140.0, r, th)])
    print(f"wrote {out}")
    print("peak sigma_tt at r = a, theta = 90:", f"{kirsch(140.0, 1, 90)[1]:.2f} MPa (Kt = 3)")
    print("r = 1.5a, theta = 0:", ["%.2f" % v for v in kirsch(140.0, 1.5, 0)], "(doc 2100, -3080, 0)")

    print("\n== Phase 3: tensile ==")
    d = read_csv("tensile_true.csv")
    eps, sig = d[:, 0], d[:, 1]
    E = 207423.0
    sy = 310.5611
    ep = eps - sig / E
    mask = (sig > sy) & (ep > 1e-6)
    m, c = np.polyfit(np.log(ep[mask]), np.log(sig[mask] - sy), 1)
    print(f"log-log fit with sigma_y = {sy}: m = {m:.4f}, K = {exp(c):.1f} MPa (doc 0.4077, 559.27)")
    sA, nu = 457.2405, 0.3
    e_el = sA / E
    e_pl = 0.031581 - e_el
    print(f"point A: eps_el = {e_el:.5f}, eps_pl = {e_pl:.5f}, eps_xx = {-nu*e_el - 0.5*e_pl:.5f} (doc -0.01535)")
    for label, (sy_, K_, m_) in {"our fit": (sy, exp(c), m), "doc fit": (310.862, 559.27, 0.4077)}.items():
        s_pred = sy_ + K_ * e_pl**m_
        ep_pred = ((sA - sy_) / K_) ** (1 / m_)
        print(f"Ludwik {label}: sigma at eps_p_A = {s_pred:.1f} MPa ({(s_pred/sA-1)*100:+.1f}%), "
              f"eps_p at 457.24 MPa = {ep_pred:.5f} ({(ep_pred/e_pl-1)*100:+.1f}%)")

    print("\n== Phase 4: thin disk, T = dT (1 - r^2/b^2), plane stress ==")
    b, n = 1.0, 20001
    rr = np.linspace(0, b, n)
    T = 1.0 - rr**2 / b**2
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (T[1:] * rr[1:] + T[:-1] * rr[:-1]) * np.diff(rr))])
    tot = cum[-1] / b**2
    for frac in [0.25, 0.5, 1.0]:
        i = int(frac * (n - 1))
        r = rr[i]
        s_r = tot - cum[i] / r**2
        s_t = tot + cum[i] / r**2 - T[i]
        print(f"r/b = {frac}: sigma_r/(E a dT) = {s_r:.5f} vs {-(1 - r*r)/4:.5f}, "
              f"sigma_t/(E a dT) = {s_t:.5f} vs {-(1 - 3*r*r)/4:.5f}")


if __name__ == "__main__":
    main()
