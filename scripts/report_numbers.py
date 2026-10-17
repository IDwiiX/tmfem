#!/usr/bin/env python3
"""Collect every computed number quoted in report/main.tex (and reused in the README prose), with its source.

Inputs (all written by a full run, see README "Build and test"):
  results/validation/*.csv         one row per asserted or reported value (tests/support/validation_record.hpp)
  results/ctest_full.log           ctest -V output; values the tests print but do not record
  results/phase*/..., results/phase4/beam_block/*/{sweep,history}.csv   solver and demo outputs
  cases/beam_block/*.json, references/data/janaf_graphite_cp.csv        inputs for derived demo numbers
Outputs:
  report/generated/numbers.tex        one macro per number (\\N<Name>), the report quotes nothing else
  report/generated/results_table.tex  condensed results table, one row per reference ID
  docs/PROVENANCE.md                  every number, its value, source file and producing test or command

Run after ctest -V and the demo cases:  python3 scripts/report_numbers.py
Any missing input or unmatched log pattern is an error, so a renamed printf cannot silently drop a number.
"""
import csv
import glob
import json
import math
import os
import re
import subprocess
import sys
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
import build_results_table as brt  # noqa: E402  (shares BENCHMARK, id_key, load order)

RESULTS = os.path.join(ROOT, "results")
LOG = os.path.join(RESULTS, "ctest_full.log")
GEN = os.path.join(ROOT, "report", "generated")
PROVENANCE = os.path.join(ROOT, "docs", "PROVENANCE.md")
DEMO = os.path.join(RESULTS, "phase4", "beam_block")


class Numbers:
    """Ordered registry: macro name -> (LaTeX value, plain value, source, producer)."""

    def __init__(self):
        self.items = []
        self.names = set()

    def add(self, name, value, source, producer, sig=4, fmt=None, unit=""):
        if name in self.names:
            raise SystemExit(f"duplicate macro {name}")
        if not re.fullmatch(r"[A-Za-z]+", name):
            raise SystemExit(f"macro name must be letters only: {name}")
        self.names.add(name)
        plain = fmt.format(value) if fmt else f"{value:.{sig}g}"
        self.items.append((name, plain, unit, source, producer))
        return value


# ----------------------------------------------------------------------------------------------------------
# inputs

def validation_rows():
    rows = []
    for path in glob.glob(os.path.join(RESULTS, "validation", "*.csv")):
        with open(path, newline="") as f:
            for r in csv.DictReader(f):
                r["file"] = os.path.relpath(path, ROOT)
                rows.append(r)
    if not rows:
        raise SystemExit("no results/validation/*.csv: run ctest first")
    return rows


def find_row(rows, rid, quantity_part):
    hits = [r for r in rows if r["id"] == rid and r["quantity"] == quantity_part]
    hits = hits or [r for r in rows if r["id"] == rid and quantity_part in r["quantity"]]
    if len(hits) != 1:
        raise SystemExit(f"{rid} '{quantity_part}': expected one row, found {len(hits)}")
    return hits[0]


def read_log():
    if not os.path.exists(LOG):
        raise SystemExit(f"{LOG} missing: run  ctest --test-dir build -V > results/ctest_full.log 2>&1")
    with open(LOG) as f:
        return [re.sub(r"^\d+: ", "", line.rstrip("\n")) for line in f]


def log_match(lines, pattern, after=None):
    """First match of pattern; with `after`, only in the lines following the first line matching `after`."""
    start = 0
    if after:
        start = next((i for i, l in enumerate(lines) if re.search(after, l)), None)
        if start is None:
            raise SystemExit(f"log: no line matches {after!r}")
    for line in lines[start:]:
        m = re.search(pattern, line)
        if m:
            return m
    raise SystemExit(f"log: no line matches {pattern!r} (after {after!r})")


def read_csv_commented(path):
    with open(path) as f:
        return list(csv.DictReader(l for l in f if not l.startswith("#")))


def git_head():
    try:
        return subprocess.check_output(["git", "-C", ROOT, "rev-parse", "--short", "HEAD"], text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


# ----------------------------------------------------------------------------------------------------------
# tolerance use: |error| / tolerance in the tolerance's own units (None for report-only rows)

def tolerance_use(r):
    tol, ref, val = r["tolerance"].strip(), float(r["reference"]), float(r["tmfem"])
    if tol in ("report only", "not asserted") or math.isnan(ref):
        return None
    diff = abs(val - ref)
    rel = diff / abs(ref) if ref != 0 else diff
    m = re.match(r">=\s*([\d.]+)", tol)
    if m:
        return float(m.group(1)) / val
    m = re.match(r"\+/-\s*([\d.eE+-]+)", tol)
    if m:
        return diff / float(m.group(1))
    m = re.match(r"([\d.eE+-]+)%", tol)
    if m:
        return rel * 100 / float(m.group(1))
    m = re.match(r"([\d.eE+-]+)\s*relative", tol)
    if m:
        return rel / float(m.group(1))
    m = re.match(r"([\d.eE+-]+)\s*(K|MPa|absolute)", tol)
    if m:
        return diff / float(m.group(1))
    m = re.fullmatch(r"([\d.eE+-]+)", tol)
    if m:
        return diff / float(m.group(1))
    raise SystemExit(f"cannot parse tolerance {tol!r} ({r['id']} {r['quantity']})")


# ----------------------------------------------------------------------------------------------------------
# number groups

def validation_numbers(N, rows):
    """Values recorded by the tests. Macro names are listed here explicitly so the report can rely on them."""
    spec = [
        ("PbTimeFE", "P1-b", "time finite element", 6), ("PbTimeRef", "P1-b", "time analytic", 6),
        ("PhDrift", "P1-h", "dose drift", 2), ("PhSlope", "P1-h", "peak decay slope", 6),
        ("PiBiot", "P1-i", "Biot", 3), ("PiTau", "P1-i", "time constant", 4),
        ("PjDiff", "P1-j", "inner face", 2),
        ("PkDiff", "P1-k", "through-thickness", 3),
        ("PmBE", "P1-m", "time order BE", 3), ("PmCN", "P1-m", "time order CN", 3),
        ("PmBEf", "P1-m", "forced time order BE consistent", 3), ("PmCNf", "P1-m", "forced time order CN consistent", 3),
        ("PmSpace", "P1-m", "space order consistent", 3), ("PmSpaceLumped", "P1-m", "space order lumped", 3),
        ("PaPatchTthree", "P2-a", "T3 max", 2), ("PaPatchTsix", "P2-a", "T6 max", 2),
        ("KtTsix", "P2-c", "Kt finest T6", 6), ("KtTthree", "P2-c", "Kt finest T3", 5),
        ("PdErr", "P2-d", "max abs stress", 2), ("PeErr", "P2-e", "max abs stress", 2),
        ("PfTthreeLtwo", "P2-f", "T3 plane stress L2", 3), ("PfTthreeEn", "P2-f", "T3 plane stress energy", 3),
        ("PfTsixLtwo", "P2-f", "T6 plane stress L2", 3), ("PfTsixEn", "P2-f", "T6 plane stress energy", 3),
        ("PfTsixLtwoPe", "P2-f", "T6 plane strain L2", 3),
        ("FitSy", "P3-a", "LSQ sigma_y0", 4), ("FitK", "P3-a", "LSQ K", 4), ("FitM", "P3-a", "LSQ m", 3),
        ("FitS", "P3-a", "residual standard deviation", 3),
        ("LogLogM", "P3-a", "log-log m (sigma_y fixed)", 3), ("LogLogK", "P3-a", "log-log K [MPa] (sigma_y fixed)", 3),
        ("PbEyy", "P3-b", "eps_yy", 5), ("PbExx", "P3-b", "eps_xx", 5),
        ("PcZ", "P3-c", "driver - data", 3),
        ("PdRatePe", "P3-d", "plane-strain bar", 3), ("PdRatePs", "P3-d", "plane-stress bar", 3),
        ("PeBar", "P3-e", "FE vs driver", 2),
        ("PsEpsFE", "P3-f", "plane stress, T6 level 4: hole-edge eps_yy", 4),
        ("PsAlphaFE", "P3-f", "plane stress, T6 level 4: max alpha", 3),
        ("PsSigFE", "P3-f", "plane stress, T6 level 4: sigma_yy at the point nearest", 4),
        ("PsZone", "P3-f", "plane stress, T6 level 4: plastic zone", 3),
        ("PsSyyMax", "P3-f", "plane stress, T6 level 4: max sigma_yy", 4),
        ("PeZone", "P3-f", "plane strain, T6 level 4: plastic zone", 3),
        ("PeSyyEdge", "P3-f", "plane strain, T6 level 4: sigma_yy at the point nearest", 4),
        ("DiskSrZero", "P4-c", "r/b = 0.0", 5), ("DiskSrHalf", "P4-c", "r/b = 0.5", 5),
        ("DiskSrOne", "P4-c", "r/b = 1.0", 2),
        ("DiskStZero", "P4-d", "r/b = 0.0", 5), ("DiskStHalf", "P4-d", "r/b = 0.5", 5),
        ("DiskStOne", "P4-d", "r/b = 1.0", 6),
        ("GaussVarTthree", "P4-e", "T3 max rel variance", 2), ("GaussVarTsix", "P4-e", "T6 max rel variance", 2),
        ("GaussEnTsix", "P4-e", "T6 energy drift", 2),
        ("PgTthreeLtwo", "P4-g", "T3 plane stress L2", 3), ("PgTthreeEn", "P4-g", "T3 plane stress energy", 3),
        ("PgTsixLtwo", "P4-g", "T6 plane stress L2", 3), ("PgTsixEn", "P4-g", "T6 plane stress energy", 3),
        ("PhTthree", "P4-h", "T3 L2", 3), ("PhTsix", "P4-h", "T6 L2", 3),
    ]
    orders = {"PmBE", "PmCN", "PmBEf", "PmCNf", "PmSpace", "PmSpaceLumped", "PfTthreeLtwo", "PfTthreeEn",
              "PfTsixLtwo", "PfTsixEn", "PfTsixLtwoPe", "PdRatePe", "PdRatePs", "PgTthreeLtwo", "PgTthreeEn",
              "PgTsixLtwo", "PgTsixEn", "PhTthree", "PhTsix"}
    for name, rid, part, sig in spec:
        r = find_row(rows, rid, part)
        # Measured orders keep two decimals, so 2.0005 does not print as an exact "2".
        N.add(name, float(r["tmfem"]), r["file"], f"test `{r['test']}`", sig=sig, fmt="{:.2f}" if name in orders else None)
    # Neuber references of the plane-stress plate (recorded as the reference of the report-only rows).
    for name, part in (("NeuberEps", "plane stress, T6 level 4: hole-edge eps_yy"),
                       ("NeuberAlpha", "plane stress, T6 level 4: max alpha"),
                       ("NeuberSig", "plane stress, T6 level 4: sigma_yy at the point nearest")):
        r = find_row(rows, "P3-f", part)
        N.add(name, float(r["reference"]), r["file"] + " (reference column)", f"test `{r['test']}`", sig=4)
    for name, part in (("PsEpsDiff", "hole-edge eps_yy"), ("PsAlphaDiff", "max alpha"),
                       ("PsSigDiff", "sigma_yy at the point nearest")):
        r = find_row(rows, "P3-f", "plane stress, T6 level 4: " + part)
        d = (float(r["tmfem"]) / float(r["reference"]) - 1) * 100
        N.add(name, d, r["file"], f"FE / Neuber - 1, test `{r['test']}`", fmt="{:+.1f}", unit="%")
    pf = [r for r in rows if r["id"] == "P2-f"]
    diff = max(abs(float(r["tmfem"]) - float(find_row(rows, "P2-f", r["quantity"].replace("plane strain", "plane stress"))["tmfem"]))
               for r in pf if "plane strain" in r["quantity"])
    N.add("PfPlaneDiff", diff, "results/validation/P2-f__*.csv", "max |plane strain - plane stress| observed order", sig=2)
    asserted = [r for r in rows if tolerance_use(r) is not None]
    N.add("Rows", len(rows), "results/validation/*.csv", "count of recorded rows", fmt="{:d}")
    N.add("Asserted", len(asserted), "results/validation/*.csv", "rows with a tolerance", fmt="{:d}")
    N.add("Ids", len({r["id"] for r in rows}), "results/validation/*.csv", "distinct reference IDs", fmt="{:d}")
    worst = max((r for r in asserted if not r["tolerance"].startswith(">=")), key=tolerance_use)
    N.add("WorstUse", 100 * tolerance_use(worst), worst["file"], f"largest |error|/tolerance among two-sided rows, {worst['id']}",
          fmt="{:.0f}", unit="%")


def log_numbers(N, lines):
    """Values the tests print to stdout but do not record."""
    m = log_match(lines, r"P3-a least squares on stress, n = \d+ plastic points, 3 parameters free \(dof (\d+)")
    dof = int(m.group(1))
    N.add("FitDof", dof, "results/ctest_full.log", "test `P3-a Ludwik least-squares fit with uncertainty`", fmt="{:d}")
    N.add("BandMax", math.sqrt(dof), "results/ctest_full.log", "sqrt(dof): bound on |residual|/s", sig=3)
    m = log_match(lines, r"P1-j max .*?= ([\d.]+) K \(dt = 150 s: ([\d.]+) K\)")
    N.add("PjDiffHalf", float(m.group(2)), "results/ctest_full.log", "test `P1-j steel plate 1D FE vs lumped ODE`", sig=2)
    m = log_match(lines, r"Tests? passed, (\d+) tests? failed out of (\d+)|(\d+)% tests passed, (\d+) tests failed out of (\d+)")
    total = int(m.group(5) or m.group(2))
    failed = int(m.group(4) or m.group(1))
    N.add("Tests", total, "results/ctest_full.log", "ctest summary line", fmt="{:d}")
    N.add("Failed", failed, "results/ctest_full.log", "ctest summary line", fmt="{:d}")
    skipped = sum(1 for l in lines if l.rstrip().endswith("(Skipped)"))
    N.add("Skipped", skipped, "results/ctest_full.log", "lines ending in (Skipped)", fmt="{:d}")
    # Plane-strain hole plate, level 4: sigma_zz at the Gauss point nearest the edge (biaxial notch root).
    m = log_match(lines, r"level 4 \(.*szz = ([-\d.]+)", after=r"P3-f elastoplastic hole plate, plane strain")
    N.add("PeSzz", float(m.group(1)), "results/ctest_full.log", "test `P3-f elastoplastic hole plate, plane strain, level 4 (report)`", sig=3)
    for t in (3, 6):
        m = log_match(lines, rf"P4-h T{t} level 2: .*?, (\d+)% of the error")
        N.add(f"PhTimePart{'Tthree' if t == 3 else 'Tsix'}", float(m.group(1)), "results/ctest_full.log",
              "test `P4-h Gaussian spreading convergence`", fmt="{:.0f}", unit="%")
    m = log_match(lines, r"P4-f \(ii\) T6: max \|F - A T_inf 1\|_1 / \(sum\|A_ij\| T_inf\) = ([\d.e+-]+)")
    N.add("PfLoadIdentity", float(m.group(1)), "results/ctest_full.log",
          "test `P4-f (ii) discrete energy balance with source, curved Robin rim and face sink`", sig=2)
    m = log_match(lines, r"zone r - a = ([\d.]+) a", after=r"P3-f elastoplastic hole plate, plane stress, S")
    N.add("PsZoneLevelOne", float(m.group(1)), "results/ctest_full.log", "test `P3-f elastoplastic hole plate, plane stress (report)`", sig=3)


def demo_numbers(N):
    mats = {}
    for name in ("copper", "graphite"):
        with open(os.path.join(ROOT, "cases", "beam_block", f"{name}.json")) as f:
            mats[name] = json.load(f)["material"]
    sweep = {name: {float(r["s_mm"]): r for r in read_csv_commented(os.path.join(DEMO, name, "sweep.csv"))}
             for name in mats}
    src = {name: os.path.relpath(os.path.join(DEMO, name, "sweep.csv"), ROOT) for name in mats}
    cmd = {name: f"`./build/src/tmfem cases/beam_block/{name}.json`" for name in mats}
    tag = {"copper": "Cu", "graphite": "Gr"}
    stag = {0.5: "Half", 1.0: "One", 2.0: "Two", 4.0: "Four"}
    for name, rows in sweep.items():
        for s, r in rows.items():
            t = tag[name] + stag[s]
            N.add(f"Demo{t}dT", float(r["dT_centre_K"]), src[name], cmd[name], fmt="{:.3g}")
            N.add(f"Demo{t}vm", float(r["vm_centre_MPa"]), src[name], cmd[name], fmt="{:.2g}" if float(r["vm_centre_MPa"]) < 10 else "{:.3g}")
            N.add(f"Demo{t}Tmax", float(r["T_max_end_of_pulse_K"]), src[name], cmd[name], fmt="{:.0f}")
            N.add(f"Demo{t}ratio", float(r["strength_ratio"]), src[name], cmd[name], sig=3)
            N.add(f"Demo{t}gps", float(r["gps_shift_percent"]), src[name], cmd[name], sig=2, unit="%")
            if name == "graphite":
                N.add(f"Demo{t}comp", float(r["compressive_ratio"]), src[name], cmd[name], sig=2)
    # FE centre values against the closed-form infinite-plane estimates the demo writes next to them.
    gap_dT = max(abs(float(r["dT_centre_K"]) / float(r["dT0_estimate_K"]) - 1) for rows in sweep.values() for r in rows.values())
    gap_vm = max(abs(float(r["vm_centre_MPa"]) / float(r["vm_centre_estimate_MPa"]) - 1) for rows in sweep.values() for r in rows.values())
    N.add("DemoGapdT", 100 * gap_dT, "results/phase4/beam_block/*/sweep.csv", "max |FE / estimate - 1| over both materials and all s", sig=2, unit="%")
    N.add("DemoGapvm", 100 * gap_vm, "results/phase4/beam_block/*/sweep.csv", "same, von Mises", sig=2, unit="%")
    # Material comparison from the case files: rho c ratio and the thermal-shock merit E alpha / ((1 - nu) rho c).
    rc = {n: m["rho_kg_m3"] * m["c_J_kgK"] for n, m in mats.items()}
    fom = {n: m["E_GPa"] * 1e9 * m["alpha_per_K"] / ((1 - m["nu"]) * rc[n]) for n, m in mats.items()}
    cases = "cases/beam_block/{copper,graphite}.json"
    N.add("HeatRatio", rc["copper"] / rc["graphite"], cases, "rho c (Cu) / rho c (graphite)", sig=3)
    N.add("FomCu", fom["copper"], cases, "E alpha / ((1 - nu) rho c) [Pa m^3/J]", sig=2)
    N.add("FomGr", fom["graphite"], cases, "E alpha / ((1 - nu) rho c) [Pa m^3/J]", sig=2)
    N.add("StressRatio", float(sweep["copper"][1.0]["vm_centre_MPa"]) / float(sweep["graphite"][1.0]["vm_centre_MPa"]),
          "results/phase4/beam_block/*/sweep.csv", "von Mises at the centre, Cu / graphite, s = 1 mm", sig=2)
    N.add("CuYield", mats["copper"]["strength"]["yield_MPa"], cases, "case input", fmt="{:g}")
    # Stress relaxation after the pulse, copper s = 1 mm (history.csv: first row = end of pulse).
    hist_path = os.path.join(DEMO, "copper", "s1", "history.csv")
    hist = read_csv_commented(hist_path)
    vm0 = float(hist[0]["vm_centre_MPa"])
    for t_s, name in ((1e-2, "RelaxTenms"), (1e-1, "RelaxHundredms")):
        row = min(hist, key=lambda r: abs(math.log(float(r["t_s"]) / t_s)))
        if abs(float(row["t_s"]) / t_s - 1) > 1e-6:
            raise SystemExit(f"history.csv has no output at t = {t_s} s")
        N.add(name, 100 * float(row["vm_centre_MPa"]) / vm0, os.path.relpath(hist_path, ROOT),
              f"von Mises at the centre at t = {t_s:g} s / end of pulse", sig=2, unit="%")
    janaf_numbers(N, mats["graphite"], sweep["graphite"][0.5])


def janaf_numbers(N, mat, row):
    """Graphite at s = 0.5 mm: temperature reached with the JANAF c(T) for the same deposit, against constant c."""
    path = os.path.join(ROOT, "references", "data", "janaf_graphite_cp.csv")
    tab = [(float(r["T_K"]), float(r["Cp_J_molK"])) for r in read_csv_commented(path)]
    M_kg_mol = 12.011e-3
    with open(os.path.join(ROOT, "cases", "beam_block", "graphite.json")) as f:
        case = json.load(f)
    T0 = case["T0_K"]
    s_m = 0.5e-3
    q_J_m3 = case["beam"]["E_per_length_J_m"] / (2 * math.pi * s_m ** 2)  # centre of the Gaussian deposit
    target = q_J_m3 / mat["rho_kg_m3"] * M_kg_mol  # J/mol to absorb

    def cp(T):  # piecewise-linear in the table
        for (t1, c1), (t2, c2) in zip(tab, tab[1:]):
            if t1 <= T <= t2:
                return c1 + (c2 - c1) * (T - t1) / (t2 - t1)
        raise SystemExit(f"T = {T} K outside the JANAF table")

    T, H, dT = T0, 0.0, 0.01
    while H < target:  # midpoint rule on a 0.01 K grid; error far below the quoted digits
        H += cp(T + dT / 2) * dT
        T += dT
    T -= (H - target) / cp(T)
    dT_const = q_J_m3 / (mat["rho_kg_m3"] * mat["c_J_kgK"])
    src = "references/data/janaf_graphite_cp.csv, cases/beam_block/graphite.json"
    N.add("JanafT", T, src, "integration of JANAF c(T) for the centre deposit, s = 0.5 mm", fmt="{:.0f}")
    N.add("JanafdT", T - T0, src, "same", fmt="{:.0f}")
    N.add("ConstdT", dT_const, src + ", results/.../graphite/sweep.csv", "q / (rho c), constant c", fmt="{:.0f}")
    N.add("JanafOver", 100 * (dT_const / (T - T0) - 1), src, "constant-c overestimate of dT", fmt="{:.0f}", unit="%")
    N.add("CpRatio", cp(800.0) / cp(298.15), src, "Cp(800 K) / Cp(298.15 K)", sig=2)
    if abs(dT_const / float(row["dT0_estimate_K"]) - 1) > 1e-9:
        raise SystemExit("constant-c dT disagrees with the demo's own estimate")


# ----------------------------------------------------------------------------------------------------------
# outputs

def tex_escape(s):
    rep = {"\\": r"\textbackslash{}", "&": r"\&", "%": r"\%", "$": r"\$", "#": r"\#", "_": r"\_",
           "{": r"\{", "}": r"\}", "~": r"\textasciitilde{}", "^": r"\textasciicircum{}", "|": r"\textbar{}",
           "<": r"\textless{}", ">": r"\textgreater{}"}
    return "".join(rep.get(c, c) for c in s)


def tex_num(plain):
    return rf"\num{{{plain}}}" if re.fullmatch(r"[+-]?[\d.]+(e[+-]?\d+)?", plain) else plain


def write_numbers(N, head):
    os.makedirs(GEN, exist_ok=True)
    with open(os.path.join(GEN, "numbers.tex"), "w") as f:
        f.write(f"% Generated by scripts/report_numbers.py from the run at commit {head}. Do not edit.\n")
        f.write("% Sources: docs/PROVENANCE.md\n")
        for name, plain, unit, _, _ in N.items:
            val = tex_num(plain) + (r"\,\%" if unit == "%" else "")
            f.write(f"\\newcommand{{\\N{name}}}{{{val}}}\n")


def write_results_table(rows, head):
    by_id = defaultdict(list)
    for r in rows:
        by_id[r["id"]].append(r)
    lines = [f"% Generated by scripts/report_numbers.py from the run at commit {head}. Do not edit.",
             r"\begin{tabular}{@{}l p{2.8cm} p{4.1cm} r r p{1.8cm} r@{}}", r"\toprule",
             r"ID & Benchmark & Quantity (largest use of tolerance) & Reference & tmfem & Tolerance & Use \\",
             r"\midrule"]
    for rid in sorted(by_id, key=lambda i: brt.id_key({"id": i, "quantity": ""})):
        group = by_id[rid]
        scored = [(tolerance_use(r), r) for r in group if tolerance_use(r) is not None]
        bench = tex_escape(brt.BENCHMARK.get(rid, ""))
        if not scored:
            lines.append(rf"{rid} & {bench} & \multicolumn{{5}}{{l}}{{report only ({len(group)} values), see \S\ref{{sec:plasticity}}}} \\")
            continue
        use, r = max(scored, key=lambda x: x[0])
        lines.append(rf"{rid} & {bench} & {tex_escape(r['quantity'])} & {tex_num(brt.fmt(r['reference']))} & "
                     rf"{tex_num(brt.fmt(r['tmfem']))} & {tex_escape(r['tolerance'])} & {100 * use:.0f}\,\% \\")
    for rid, (reason, _) in brt.SKIPPED.items():
        if rid not in by_id:
            lines.append(rf"{rid} & {tex_escape(brt.BENCHMARK[rid])} & \multicolumn{{5}}{{l}}{{{tex_escape(reason)}}} \\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    with open(os.path.join(GEN, "results_table.tex"), "w") as f:
        f.write("\n".join(lines) + "\n")


def write_provenance(N, head):
    out = ["# Provenance of quoted numbers", "",
           f"Generated by `python3 scripts/report_numbers.py` from the run at commit `{head}`. Do not edit by hand.",
           "", "Every computed number in `report/main.tex` is one of these macros (`\\N<Name>`) or a row of the",
           "generated results table. Inputs and reference values (material data, loads, the Abaqus strains, the",
           "source's Ludwik parameters and the source-corrections table) are quoted from the case files and",
           "`docs/REFERENCE_VALUES.md`, not listed here. The computed numbers in the README prose are written",
           "from this list into its marked blocks by the same script. Numbers from runs before a disclosed change are labelled *historical* where",
           "they appear (docs/VALIDATION.md) and are not listed here.", "",
           "Producers: `test` = Catch2 test (`ctest --test-dir build -V > results/ctest_full.log`),",
           "`./build/src/tmfem <case>` = case run, other text = the calculation done in this script.", "",
           "| Macro | Value | Source file | Producer |", "| --- | --- | --- | --- |"]
    for name, plain, unit, source, producer in N.items:
        out.append(f"| `\\N{name}` | {plain}{' %' if unit == '%' else ''} | `{source}` | {producer} |")
    with open(PROVENANCE, "w") as f:
        f.write("\n".join(out) + "\n")


README = os.path.join(ROOT, "README.md")


def readme_blocks(N):
    """Write the README's computed prose numbers from the same registry as the report."""
    v = {name: plain for name, plain, _, _, _ in N.items}
    pct = lambda name: f"{100 * float(v[name]):.2g}%"
    demo = (f"A Gaussian beam deposits 1 kJ per metre of path in 10 µs into a long 40 mm x 40 mm block of copper or\n"
            f"graphite, solved in plane strain (properties and sources: [cases/beam_block/README.md](cases/beam_block/README.md)).\n"
            f"At equal deposited energy, graphite has a {v['HeatRatio']} times larger temperature rise than copper\n"
            f"({v['DemoGrOnedT']} K vs {v['DemoCuOnedT']} K at s = 1 mm) yet carries {v['StressRatio']} times less von Mises stress\n"
            f"({v['DemoGrOnevm']} vs {v['DemoCuOnevm']} MPa at the centre), because the stress scales with Eα/((1−ν)ρc). Annealed\n"
            f"copper exceeds yield for s ≤ 1 mm (von Mises / yield = {v['DemoCuOneratio']} at 1 mm); graphite reaches\n"
            f"{pct('DemoGrHalfratio')} of its flexural strength at s = 0.5 mm, an optimistic tensile check. The demo illustrates\n"
            f"graphite's thermal-shock resistance only; its lower energy density for the same beam is outside it. These are\n"
            f"elastic, constant-property estimates; the Gaussian stands in for a FLUKA or Geant4 energy map. Full\n"
            f"interpretation in [docs/VALIDATION.md](docs/VALIDATION.md#demo-d4-beam-heated-block-not-a-validation) and in\n"
            f"section 4 of the report. Numbers generated by `scripts/report_numbers.py` (sources: docs/PROVENANCE.md).")
    limits = (f"- **Demo:** plane strain instead of generalised plane strain (a σ_zz shift of {v['DemoCuHalfgps']}% to\n"
              f"  {v['DemoGrFourgps']}% of the peak |σ_zz|, twice that fraction in the centre von Mises); elastic and\n"
              f"  constant-property (graphite at s = 0.5 mm overestimates its temperature rise by {v['JanafOver']}%); assumed\n"
              f"  graphite ν and c; one-way, quasi-static coupling.")
    with open(README) as f:
        text = f.read()
    for tag, body in (("demo", demo), ("demo-limits", limits)):
        begin, end = f"<!-- {tag}:begin -->", f"<!-- {tag}:end -->"
        if begin not in text or end not in text:
            raise SystemExit(f"README.md lacks the {begin} / {end} markers")
        head, rest = text.split(begin, 1)
        text = f"{head}{begin}\n{body}\n{end}{rest.split(end, 1)[1]}"
    with open(README, "w") as f:
        f.write(text)


def main():
    rows = validation_rows()
    lines = read_log()
    N = Numbers()
    validation_numbers(N, rows)
    log_numbers(N, lines)
    demo_numbers(N)
    head = git_head()
    write_numbers(N, head)
    write_results_table(rows, head)
    write_provenance(N, head)
    readme_blocks(N)
    print(f"wrote {len(N.items)} numbers, results table ({len({r['id'] for r in rows})} IDs) and {PROVENANCE}")


if __name__ == "__main__":
    main()
