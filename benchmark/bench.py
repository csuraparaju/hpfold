#!/usr/bin/env python3
"""Benchmark mc_hp_pt across easy, stress, and scaling workloads.

Usage:
    python bench.py [suite] [--binary PATH] [--csv FILE]

Suites:
    easy            small problems, quick sanity check (~30s)
    stress          published HP benchmark sequences (~minutes)
    scale-replicas  wall time vs replica count
    scale-steps     wall time vs steps_per_swap
    scale-swaps     wall time vs swap round count
    all             everything (long)
"""

import argparse
import csv
import re
import statistics
import subprocess
import sys
import time
from dataclasses import dataclass
from typing import Optional


# Canonical 2D HP benchmark sequences with published ground-state energies.
# References: Unger & Moult 1993, Lesh et al. 2003, etc.
# These are the standard test sequences used in HP-folding literature.
BENCH = {
    "hp20": ("HPHPPHHPHPPHPHHPPHPH",                                                  -9),
    "hp24": ("HHPPHPPHPPHPPHPPHPPHPPHH",                                              -9),
    "hp25": ("PPHPPHHPPPPHHPPPPHHPPPPHH",                                             -8),
    "hp36": ("PPPHHPPHHPPPPPHHHHHHHPPHHPPPPHHPPHPP",                                 -14),
    "hp48": ("PPHPPHHPPHHPPPPPHHHHHHHHHHPPPPPPHHPPHHPPHPPHHHHH",                     -23),
    "hp50": ("HHPHPHPHPHHHHPHPPPHPPPHPPPPHPPPHPPPHPHHHHPHPHPHPHH",                   -21),
    "hp60": ("PPHHHPHHHHHHHHPPPHHHHHHHHHHPHPPPHHHHHHHHHHHHPPPPHHHHHHPHHPHP",         -36),
}


@dataclass
class Case:
    name: str
    sequence: str
    replicas: int = 128
    t_min: float = 0.5
    t_max: float = 20.0
    steps_per_swap: int = 5000
    n_swaps: int = 200
    repeats: int = 1
    known_optimum: Optional[int] = None


@dataclass
class Run:
    wall_time: float
    best_energy: int
    swap_acceptance: float


ENERGY_RE = re.compile(r"Best chain \(energy (-?\d+)\)")
ACCEPT_RE = re.compile(r"swap acceptance:\s*([\d.]+)")


def run_one(binary: str, c: Case) -> Run:
    """Invoke the binary once and parse one run's results."""
    cmd = [
        binary,
        "-s", c.sequence,
        "-r", str(c.replicas),
        "-l", str(c.t_min),
        "-u", str(c.t_max),
        "-m", str(c.steps_per_swap),
        "-n", str(c.n_swaps),
    ]
    t0 = time.perf_counter()
    proc = subprocess.run(cmd, capture_output=True, text=True)
    elapsed = time.perf_counter() - t0

    if proc.returncode != 0:
        raise RuntimeError(
            f"{binary} exited {proc.returncode}: {proc.stderr.strip() or '(no stderr)'}"
        )

    em = ENERGY_RE.search(proc.stdout)
    am = ACCEPT_RE.search(proc.stdout)
    if not em or not am:
        raise RuntimeError(f"could not parse output:\n{proc.stdout}")

    return Run(elapsed, int(em.group(1)), float(am.group(1)))


def run_case(binary: str, c: Case, verbose: bool = True) -> list:
    """Run a single case (possibly multiple repeats) and return all Run results."""
    if verbose:
        print(
            f"[{c.name}] n={len(c.sequence)} R={c.replicas} "
            f"steps={c.steps_per_swap} swaps={c.n_swaps} "
            f"T=[{c.t_min:g},{c.t_max:g}] repeats={c.repeats}"
        )
    runs = []
    for i in range(c.repeats):
        r = run_one(binary, c)
        runs.append(r)
        if verbose:
            tag = ""
            if c.known_optimum is not None:
                gap = r.best_energy - c.known_optimum
                tag = f" (opt={c.known_optimum}, gap={gap:+d})"
            print(
                f"  run {i+1}/{c.repeats}: E={r.best_energy}{tag}  "
                f"accept={r.swap_acceptance:.3f}  t={r.wall_time:.2f}s"
            )
    return runs


def summarize(c: Case, runs: list) -> dict:
    e = [r.best_energy for r in runs]
    t = [r.wall_time for r in runs]
    a = [r.swap_acceptance for r in runs]
    return {
        "name": c.name,
        "n": len(c.sequence),
        "replicas": c.replicas,
        "steps_per_swap": c.steps_per_swap,
        "n_swaps": c.n_swaps,
        "repeats": len(runs),
        "e_best": min(e),
        "e_mean": statistics.mean(e),
        "t_mean": statistics.mean(t),
        "t_std": statistics.stdev(t) if len(t) > 1 else 0.0,
        "accept": statistics.mean(a),
        "opt": c.known_optimum,
        "gap": (min(e) - c.known_optimum) if c.known_optimum is not None else None,
    }


def print_summary(rows: list) -> None:
    if not rows:
        return
    print()
    print(
        f"{'name':<14} {'n':>3} {'R':>3} {'E_best':>7} {'E_mean':>7} "
        f"{'t(s)':>14} {'accept':>7} {'gap':>5}"
    )
    print("-" * 72)
    for r in rows:
        gap = f"{r['gap']:+d}" if r['gap'] is not None else "   - "
        t_field = f"{r['t_mean']:>5.2f}±{r['t_std']:.2f}"
        print(
            f"{r['name']:<14} {r['n']:>3} {r['replicas']:>3} "
            f"{r['e_best']:>7d} {r['e_mean']:>7.2f} "
            f"{t_field:>14} {r['accept']:>7.3f} {gap:>5}"
        )
    print()


# ---------------- Test suites ----------------

def easy_suite():
    """Quick sanity tests, ~30s total."""
    return [
        Case("smoke", "HPHPPHHPHPP", replicas=4, n_swaps=20, repeats=3),
        Case("hp20", BENCH["hp20"][0], replicas=8, n_swaps=50, repeats=3,
             known_optimum=BENCH["hp20"][1]),
        Case("hp24", BENCH["hp24"][0], replicas=8, n_swaps=100, repeats=3,
             known_optimum=BENCH["hp24"][1]),
    ]


def stress_suite():
    """Published HP benchmarks; checks both speed and solution quality."""
    return [
        Case("hp20_thorough", BENCH["hp20"][0], replicas=16, n_swaps=500, repeats=3,
             known_optimum=BENCH["hp20"][1]),
        Case("hp25", BENCH["hp25"][0], replicas=16, n_swaps=500, repeats=2,
             known_optimum=BENCH["hp25"][1]),
        Case("hp36", BENCH["hp36"][0], replicas=16, n_swaps=1000, repeats=2,
             known_optimum=BENCH["hp36"][1]),
        Case("hp48", BENCH["hp48"][0], replicas=32, n_swaps=2000, repeats=2,
             known_optimum=BENCH["hp48"][1]),
        Case("hp60", BENCH["hp60"][0], replicas=64, n_swaps=3000, repeats=1,
             known_optimum=BENCH["hp60"][1]),
    ]


def scale_suite(dim: str):
    """One-axis scaling sweep on a fixed sequence."""
    seq, opt = BENCH["hp24"]
    if dim == "replicas":
        return [Case(f"R={r:02d}", seq, replicas=r, n_swaps=200, repeats=2,
                     known_optimum=opt)
                for r in (2, 4, 8, 16, 32, 64)]
    if dim == "steps":
        return [Case(f"m={m}", seq, replicas=8, steps_per_swap=m, n_swaps=100,
                     repeats=2, known_optimum=opt)
                for m in (500, 1000, 5000, 10000, 25000, 50000)]
    if dim == "swaps":
        return [Case(f"n={n}", seq, replicas=8, n_swaps=n, repeats=2,
                     known_optimum=opt)
                for n in (50, 100, 250, 500, 1000, 2000)]
    raise ValueError(f"unknown scaling dim: {dim}")


SUITES = {
    "easy":           ("Easy tests",            easy_suite),
    "stress":         ("Stress tests",          stress_suite),
    "scale-replicas": ("Scaling: replica count", lambda: scale_suite("replicas")),
    "scale-steps":    ("Scaling: steps/swap",    lambda: scale_suite("steps")),
    "scale-swaps":    ("Scaling: swap rounds",   lambda: scale_suite("swaps")),
}


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("suite", choices=list(SUITES) + ["all"], nargs="?", default="easy")
    p.add_argument("--binary", default="./mc_hp_pt", help="path to mc_hp_pt binary")
    p.add_argument("--csv", help="write summary to CSV file")
    args = p.parse_args()

    if args.suite == "all":
        suites_to_run = list(SUITES.values())
    else:
        suites_to_run = [SUITES[args.suite]]

    all_rows = []
    for title, build in suites_to_run:
        print("=" * 72)
        print(title)
        print("=" * 72)
        rows = []
        for c in build():
            try:
                runs = run_case(args.binary, c)
                rows.append(summarize(c, runs))
            except Exception as e:
                print(f"  ERROR: {e}", file=sys.stderr)
        print_summary(rows)
        all_rows.extend(rows)

    if args.csv and all_rows:
        with open(args.csv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(all_rows[0].keys()))
            w.writeheader()
            w.writerows(all_rows)
        print(f"summary written to {args.csv}")


if __name__ == "__main__":
    main()