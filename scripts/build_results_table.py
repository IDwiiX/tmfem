#!/usr/bin/env python3
"""Rebuild the README results table from results/validation/*.csv.

Each row is written by the Catch2 test that asserts it (tests/support/validation_record.hpp), so
every number in the table traces to a test. Run after ctest:  python3 scripts/build_results_table.py
The table replaces the text between the results markers in README.md.
"""
import csv
import glob
import math
import os
import re

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
VALIDATION = os.path.join(ROOT, "results", "validation")
README = os.path.join(ROOT, "README.md")
REFERENCE_VALUES = os.path.join(ROOT, "docs", "REFERENCE_VALUES.md")
BEGIN, END = "<!-- results:begin -->", "<!-- results:end -->"

# Reference rows that no test can produce yet, with the reason shown in the table. Any other reference ID
# without a recorded row is an error, so the table always covers docs/REFERENCE_VALUES.md.
SKIPPED = {
    "P1-l": ("skipped: references/data/ttab.csv and the ANSYS curve are not provided",
             "P1-l steel plate vs ANSYS inner face"),
}

BENCHMARK = {
    "P1-a": "Steady diffusion flux (ex. 1)",
    "P1-b": "Carburizing, erf solution (ex. 2)",
    "P1-c": "Equivalent time (ex. 3)",
    "P1-d": "Arrhenius D(T) (ex. 4)",
    "P1-e": "Boron predeposition dose (ex. 5)",
    "P1-f": "Thin-film surface concentration (ex. 5)",
    "P1-g": "Arrhenius fit (ex. 5)",
    "P1-h": "Thin film on a half-space",
    "P1-i": "Steel plate parameters",
    "P1-j": "Steel plate, FE vs lumped ODE",
    "P1-k": "Steel plate, through-thickness difference",
    "P1-l": "Steel plate, real ttab vs ANSYS",
    "P1-m": "1D heat convergence orders",
    "P2-a": "Patch test, distorted mesh",
    "P2-b": "Uniform bar in tension",
    "P2-c": "Kirsch plate, stress concentration",
    "P2-d": "Kirsch plate, r/a = 2 and 4",
    "P2-e": "Kirsch plate, r/a = 1.5",
    "P2-f": "Kirsch convergence orders",
    "P3-a": "Ludwik fit to the true tensile curve",
    "P3-b": "Material point, uniaxial, tabulated hardening",
    "P3-c": "Material point, Ludwik along the true curve",
    "P3-d": "Newton convergence rate",
    "P3-e": "Plane-strain bar FE vs material point",
    "P3-f": "Elastoplastic hole plate, S = 140 MPa (report only)",
    "P4-a": "Bar held in x, uniform dT",
    "P4-b": "Free body, uniform dT",
    "P4-c": "Thin disk, parabolic T: sigma_r",
    "P4-d": "Thin disk, parabolic T: sigma_theta",
    "P4-e": "2D Gaussian heat spreading",
    "P4-f": "2D heat with source and convection",
    "P4-g": "Thin disk convergence orders",
    "P4-h": "2D Gaussian heat convergence orders",
}


def id_key(row):
    m = re.match(r"P(\d+)-([a-z])", row["id"])
    return (int(m.group(1)), m.group(2), row["quantity"]) if m else (99, row["id"], row["quantity"])


def fmt(v):
    v = float(v)
    return "n/a" if math.isnan(v) else f"{v:.6g}"


def error_text(row):
    """Error against the reference. Report-only rows compare with a value that is not a target (an
    estimate or the source's parameters), so their difference is shown in parentheses; no reference: n/a."""
    ref, val = float(row["reference"]), float(row["tmfem"])
    if math.isnan(ref):
        return "n/a"
    if ref == 0.0 or "+/-" in row["tolerance"] or row["tolerance"].endswith(" K"):
        text = f"{val - ref:+.3g}"
    else:
        text = f"{(val - ref) / abs(ref) * 100:+.3g}%"
    return f"({text})" if row["tolerance"] == "report only" else text


def reference_ids():
    """Row IDs of the phase tables in docs/REFERENCE_VALUES.md (stretch goals excluded)."""
    with open(REFERENCE_VALUES) as f:
        return sorted(set(re.findall(r"^\| (P\d-[a-z]) \|", f.read(), re.M)))


def load_rows():
    rows = []
    for path in glob.glob(os.path.join(VALIDATION, "*.csv")):
        with open(path, newline="") as f:
            rows.extend(csv.DictReader(f))
    return rows


def main():
    rows = load_rows()
    recorded = {r["id"] for r in rows}
    missing = [i for i in reference_ids() if i not in recorded and i not in SKIPPED]
    if missing:
        raise SystemExit(f"reference rows with no recorded result: {', '.join(missing)} (run the full ctest)")
    for rid, (reason, test) in SKIPPED.items():
        if rid not in recorded:
            rows.append({"id": rid, "quantity": reason, "reference": "nan", "tmfem": "nan",
                         "tolerance": "report only", "test": test})
    rows.sort(key=id_key)
    lines = ["| ID | Benchmark | Quantity | Reference | tmfem | Error | Tolerance | Test |",
             "| --- | --- | --- | --- | --- | --- | --- | --- |"]
    for r in rows:
        lines.append(f"| {r['id']} | {BENCHMARK.get(r['id'], '')} | {r['quantity']} | {fmt(r['reference'])} | "
                     f"{fmt(r['tmfem'])} | {error_text(r)} | {r['tolerance']} | `{r['test']}` |")
    table = "\n".join(lines)
    with open(README) as f:
        text = f.read()
    if BEGIN not in text or END not in text:
        raise SystemExit(f"README.md lacks the {BEGIN} / {END} markers")
    head, rest = text.split(BEGIN, 1)
    _, tail = rest.split(END, 1)
    with open(README, "w") as f:
        f.write(f"{head}{BEGIN}\n{table}\n{END}{tail}")
    print(f"wrote {len(rows)} rows to {README}")


if __name__ == "__main__":
    main()
