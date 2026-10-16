#!/usr/bin/env python3
"""Rebuild the README results table from results/validation/*.csv.

Each row is written by the Catch2 test that asserts it (tests/support/validation_record.hpp), so
every number in the table traces to a test. Run after ctest:  python3 scripts/build_results_table.py
The table replaces the text between the results markers in README.md.
"""
import csv
import glob
import os
import re

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
VALIDATION = os.path.join(ROOT, "results", "validation")
README = os.path.join(ROOT, "README.md")
BEGIN, END = "<!-- results:begin -->", "<!-- results:end -->"

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
    "P1-m": "1D heat convergence orders",
}


def id_key(row):
    m = re.match(r"P(\d+)-([a-z])", row["id"])
    return (int(m.group(1)), m.group(2), row["quantity"]) if m else (99, row["id"], row["quantity"])


def fmt(v):
    return f"{float(v):.6g}"


def error_text(row):
    ref, val = float(row["reference"]), float(row["tmfem"])
    if ref == 0.0 or "+/-" in row["tolerance"] or row["tolerance"].endswith(" K"):
        return f"{val - ref:+.3g}"
    return f"{(val - ref) / abs(ref) * 100:+.3g}%"


def main():
    rows = []
    for path in glob.glob(os.path.join(VALIDATION, "*.csv")):
        with open(path, newline="") as f:
            rows.extend(csv.DictReader(f))
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
