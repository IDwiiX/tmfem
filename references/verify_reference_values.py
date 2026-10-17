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

    print("\n== Phase 4: held bar (P4-a) and free body (P4-b), uniform dT, from 3D Hooke's law ==")
    # sigma = C (eps - alpha dT I) with eps = diag(exx, eyy, ezz); unknown strains solved from the zero stresses.
    for E, nu, al, dT in [(200e9, 0.3, 1.2e-5, 100.0), (117e9, 0.34, 16.9e-6, 46.2), (11.5e9, 0.2, 4.2e-6, 123.0)]:
        C = E / ((1 + nu) * (1 - 2 * nu)) * np.array([[1 - nu, nu, nu], [nu, 1 - nu, nu], [nu, nu, 1 - nu]])
        th = C @ (al * dT * np.ones(3))  # thermal stress of fully prevented expansion

        def solve(fixed, free):
            """Strains in `fixed` are 0; stresses in `free` are 0. Returns (strain, stress)."""
            e = np.zeros(3)
            e[free] = np.linalg.solve(C[np.ix_(free, free)], th[free])
            return e, C @ e - th

        ead = E * al * dT
        cases = {
            "P4-a plane stress (exx = 0; syy = szz = 0)": (solve([0], [1, 2]), {"sxx/(-E a dT)": (0, -ead)}),
            "P4-a plane strain (exx = ezz = 0; syy = 0)": (
                solve([0, 2], [1]),
                {"sxx/(-E a dT/(1-nu))": (0, -ead / (1 - nu)), "szz/(-E a dT/(1-nu))": (2, -ead / (1 - nu))}),
            "P4-b plane stress (all free)": (solve([], [0, 1, 2]), {}),
            "P4-b plane strain (ezz = 0; sxx = syy = 0)": (solve([2], [0, 1]), {"szz/(-E a dT)": (2, -ead)}),
        }
        print(f"E = {E/1e9:.1f} GPa, nu = {nu}, alpha dT = {al*dT:.3e}")
        for name, ((e, s), ratios) in cases.items():
            parts = [f"{k} = {s[i]/ref:.12f}" for k, (i, ref) in ratios.items()]
            zero = [i for i in range(3) if i not in [j for j, _ in ratios.values()]]
            parts.append(f"max |zero stress|/(E a dT) = {max(abs(s[i]) for i in zero)/ead:.1e}")
            if name.startswith("P4-a"):
                eyy_ref = (1 + nu) * al * dT if "stress" in name else (1 + nu) / (1 - nu) * al * dT
                parts.append(f"eyy/target = {e[1]/eyy_ref:.12f}")
            else:
                ex_ref = al * dT if "stress" in name else (1 + nu) * al * dT
                parts.append(f"exx/alpha_eff dT = {e[0]/ex_ref:.12f}")
            print(f"  {name}: " + ", ".join(parts))


    print("\n== Phase 4: P4-f (i) steady slab, P4-f (iii) amplification, P4-e/h Gaussian ==")
    L, q, h, k, T_inf = 0.012, 1e6, 100.0, 1.0, 300.0
    T = lambda x: T_inf + q * L / (2 * h) + q * x * (L - x) / (2 * k)
    dx = 1e-6
    for x in (0.0, 0.003, 0.006, 0.011):
        xc = min(max(x, 2 * dx), L - 2 * dx)
        print(f"slab: x = {x:.3f}: k T'' + q = {k * (T(xc + dx) - 2 * T(xc) + T(xc - dx)) / dx**2 + q:+.2e} W/m3")
    print(f"slab: Robin at x = 0: k T'(0) - h (T(0) - T_inf) = {k * (T(dx) - T(-dx)) / (2 * dx) - h * (T(0) - T_inf):+.2e}; "
          f"at x = L: -k T'(L) - h (T(L) - T_inf) = {-k * (T(L + dx) - T(L - dx)) / (2 * dx) - h * (T(L) - T_inf):+.2e}")
    lam_dt = 0.5
    for theta in (1.0, 0.5):
        G = (1 - (1 - theta) * lam_dt) / (1 + theta * lam_dt)
        y = 1.0
        for _ in range(20):  # scalar theta step for y' = -lam y
            y = (y - (1 - theta) * lam_dt * y) / (1 + theta * lam_dt)
        print(f"G^n: theta = {theta}: recurrence {y:.15e}, G^20 = {G**20:.15e}")
    a, s0, T0 = 400 / 3.44e6, 1e-3, 100.0
    Tg = lambda x, y, t: T0 * s0**2 / (s0**2 + 2 * a * t) * np.exp(-(x * x + y * y) / (2 * (s0**2 + 2 * a * t)))
    t, d, dt = 4e-3, 2e-5, 1e-7
    for (x, y) in ((0.0, 0.0), (1e-3, 5e-4), (2.5e-3, 1e-3)):
        lap = (Tg(x + d, y, t) + Tg(x - d, y, t) + Tg(x, y + d, t) + Tg(x, y - d, t) - 4 * Tg(x, y, t)) / d**2
        res = (Tg(x, y, t + dt) - Tg(x, y, t - dt)) / (2 * dt) - a * lap
        print(f"Gaussian: (x, y) = ({x:.1e}, {y:.1e}): (dT/dt - a lap T)/(a |lap T| + 1e-30) = {res / (a * abs(lap)):+.1e}")
    g = np.linspace(-0.02, 0.02, 2001)
    X, Y = np.meshgrid(g, g)
    for t in (0.0, 4.3e-3, 1.29e-2):
        F = Tg(X, Y, t)
        var = (X**2 * F).sum() / F.sum()
        E = F.sum() * (g[1] - g[0])**2
        print(f"Gaussian: t = {t:.2e} s: variance {var:.6e} vs s0^2 + 2at {s0**2 + 2 * a * t:.6e}; "
              f"integral {E:.9e} vs 2 pi s0^2 T0 {2 * pi * s0**2 * T0:.9e}")


if __name__ == "__main__":
    main()
