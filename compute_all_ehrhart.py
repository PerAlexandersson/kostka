#!/usr/bin/env python3
"""
Compute dimension, Ehrhart polynomial, h*-vector for all GT(λ, w)
with |λ| ≤ MAX_N and λ ≥ w in dominance order (non-skew only).

Uses the kostka binary's table subcommand for batch computation per λ.
Output: gt_ehrhart_data.csv

Flag --kkt-counterexample: halt immediately if any Ehrhart polynomial
has a negative coefficient, and print the counterexample.
"""

import argparse
import subprocess
import json
import csv
import re
import time
import sys
import os

MAX_N = 9
KOSTKA = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                      "target", "release", "kostka")
OUTPUT = "gt_ehrhart_data.csv"
MAX_STATES = 5_000_000   # prevent OOM; raise if you have plenty of RAM
TIMEOUT = 600             # seconds per λ


def partitions(n, max_part=None):
    """Generate all integer partitions of n (reverse lex order)."""
    if max_part is None:
        max_part = n
    if n == 0:
        yield ()
        return
    for first in range(min(n, max_part), 0, -1):
        for rest in partitions(n - first, first):
            yield (first,) + rest


def dominates(lam, mu):
    """True if λ ≥ μ in dominance order."""
    s1 = s2 = 0
    for i in range(max(len(lam), len(mu))):
        s1 += lam[i] if i < len(lam) else 0
        s2 += mu[i] if i < len(mu) else 0
        if s1 < s2:
            return False
    return True


def has_negative_coefficient(poly_str):
    """Check if the polynomial string has any negative coefficient.

    The polynomial comes from the kostka binary as e.g.:
      (1/6)n^3 + n^2 + (11/6)n + 1
      (1/6)n^3 - (1/2)n^2 + (1/3)n
    Negative coefficients appear as ' - ' between terms or as a leading '-'.
    """
    if not poly_str or poly_str in ("0", "1"):
        return False
    # Split on ' + ' and ' - ', keeping the sign
    # A term is negative if preceded by ' - ' or starts with '-'
    s = poly_str.strip()
    # Normalize: replace ' - ' with ' + -' so we can split on ' + '
    s = s.replace(" - ", " + -")
    terms = [t.strip() for t in s.split(" + ") if t.strip()]
    for t in terms:
        if t.startswith("-") or t.startswith("(-"):
            return True
    return False


def main():
    parser = argparse.ArgumentParser(description="Batch Ehrhart computation")
    parser.add_argument("--kkt-counterexample", action="store_true",
                        help="Halt on first Ehrhart polynomial with negative coefficient")
    parser.add_argument("--max-n", type=int, default=MAX_N,
                        help=f"Maximum partition size (default {MAX_N})")
    args = parser.parse_args()

    max_n = args.max_n
    check_negative = args.kkt_counterexample

    env = os.environ.copy()
    env["PATH"] = os.path.expanduser("~/.cargo/bin") + ":" + env.get("PATH", "")

    print("Building release binary...", flush=True)
    subprocess.run(["cargo", "build", "--release"], check=True, env=env,
                   capture_output=True)
    print("Build complete.\n", flush=True)
    if check_negative:
        print("*** --kkt-counterexample active: will halt on negative "
              "Ehrhart coefficient ***\n", flush=True)

    total_pairs = 0
    skipped = []

    with open(OUTPUT, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "n", "lambda", "weight", "kostka", "degree",
            "polynomial", "hstar", "lambda_time_s"
        ])

        for n in range(1, max_n + 1):
            plist = list(partitions(n))
            print(f"\n{'='*60}", flush=True)
            print(f"|λ| = {n}: {len(plist)} partitions", flush=True)
            print(f"{'='*60}", flush=True)
            n_start = time.time()

            for lam in plist:
                lam_str = ",".join(map(str, lam))

                try:
                    t0 = time.time()
                    result = subprocess.run(
                        [KOSTKA, "table",
                         "--lambda", lam_str,
                         "--all-weights", "--ehrhart", "--hstar",
                         "--weight-as-partition",
                         "--max-states", str(MAX_STATES),
                         "--format", "json"],
                        capture_output=True, text=True,
                        timeout=TIMEOUT, env=env
                    )
                    dt = time.time() - t0
                except subprocess.TimeoutExpired:
                    print(f"  λ=({lam_str}): TIMEOUT ({TIMEOUT}s)", flush=True)
                    skipped.append((n, lam_str, "timeout"))
                    continue

                if result.returncode != 0:
                    dt = time.time() - t0
                    err = result.stderr.strip().split("\n")[-1] if result.stderr else "?"
                    print(f"  λ=({lam_str}): SKIPPED ({err})", flush=True)
                    skipped.append((n, lam_str, err))
                    continue

                count = 0
                for line in result.stdout.strip().split("\n"):
                    if not line:
                        continue
                    try:
                        d = json.loads(line)
                    except json.JSONDecodeError:
                        continue

                    wt = tuple(d["weight"])
                    if not dominates(lam, wt):
                        continue

                    poly = d.get("polynomial", "") or ""
                    hs = d.get("hstar")

                    # Check for negative Ehrhart coefficients
                    if check_negative and poly and has_negative_coefficient(poly):
                        wt_str = ",".join(map(str, wt))
                        print(f"\n{'!'*60}", flush=True)
                        print(f"COUNTEREXAMPLE FOUND!", flush=True)
                        print(f"  λ  = ({lam_str})", flush=True)
                        print(f"  w  = ({wt_str})", flush=True)
                        print(f"  K  = {d.get('value', '?')}", flush=True)
                        print(f"  deg= {d.get('degree', '?')}", flush=True)
                        print(f"  P  = {poly}", flush=True)
                        print(f"  h* = {hs}", flush=True)
                        print(f"{'!'*60}", flush=True)
                        # Still write this row before exiting
                        writer.writerow([
                            n, f"({lam_str})",
                            f"({wt_str})",
                            d.get("value", ""),
                            d.get("degree", ""),
                            poly,
                            str(hs) if hs is not None else "",
                            f"{dt:.3f}"
                        ])
                        f.flush()
                        sys.exit(1)

                    writer.writerow([
                        n,
                        f"({lam_str})",
                        f"({','.join(map(str, wt))})",
                        d.get("value", ""),
                        d.get("degree", ""),
                        poly,
                        str(hs) if hs is not None else "",
                        f"{dt:.3f}"
                    ])
                    count += 1
                    total_pairs += 1

                print(f"  λ=({lam_str}): {count} weights, {dt:.3f}s", flush=True)
                f.flush()

            n_elapsed = time.time() - n_start
            print(f"  --- |λ|={n} done: {n_elapsed:.1f}s, "
                  f"{total_pairs} pairs so far ---", flush=True)

    print(f"\n{'='*60}")
    print(f"COMPLETE: {total_pairs} pairs written to {OUTPUT}")
    if check_negative:
        print("No negative Ehrhart coefficients found.")
    if skipped:
        print(f"\nSkipped {len(skipped)} lambdas:")
        for sn, sl, reason in skipped:
            print(f"  |λ|={sn} λ=({sl}): {reason}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
