#!/usr/bin/env python3
"""Benchmark: kostka lr (Rust DP) vs lrcalc (Buch's C library).

Scales c^(6,4,3,2,1)_{(3,2,1),(4,3,2,1)} = 8 by factors t = 1..5.
"""

import subprocess
import time
import json
import sys

sys.path.insert(0, "/tmp/lrcalc_bench/lib/python3.12/site-packages")
import lrcalc

KOSTKA_BIN = "./target/release/kostka"

# Base case: c^(6,4,3,2,1)_{(3,2,1),(4,3,2,1)} = 8
LAM0 = [6, 4, 3, 2, 1]
MU0  = [3, 2, 1]
NU0  = [4, 3, 2, 1]

FACTORS = [1, 2, 3, 4, 5]


def fmt_part(p):
    return ",".join(str(x) for x in p)


def scale(p, t):
    return [x * t for x in p]


def bench_lrcalc(lam, mu, nu, target_secs=0.5):
    """Time lrcalc.lrcoef with adaptive repeats."""
    val = lrcalc.lrcoef(lam, mu, nu)
    # Probe
    t0 = time.perf_counter()
    lrcalc.lrcoef(lam, mu, nu)
    probe = time.perf_counter() - t0
    repeats = max(1, min(1000000, int(target_secs / max(probe, 1e-9))))
    t0 = time.perf_counter()
    for _ in range(repeats):
        lrcalc.lrcoef(lam, mu, nu)
    elapsed = (time.perf_counter() - t0) / repeats
    return int(val), elapsed


def bench_kostka(lam, mu, nu):
    """Run kostka lr --format json and extract internal DP timing."""
    cmd = [
        KOSTKA_BIN, "lr",
        "--lambda", fmt_part(lam),
        "--mu", fmt_part(mu),
        "--nu", fmt_part(nu),
        "--format", "json",
    ]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
    except subprocess.TimeoutExpired:
        return None, None
    if r.returncode != 0:
        return None, None
    data = json.loads(r.stdout)
    return int(data["lr"]), data["dp_ms"] / 1000.0


def fmt_time(t):
    if t is None or t < 0:
        return "—"
    if t < 0.001:
        return f"{t*1e6:.1f}µs"
    if t < 1:
        return f"{t*1000:.2f}ms"
    return f"{t:.2f}s"


def main():
    subprocess.run(["cargo", "build", "--release"], capture_output=True)

    print(f"{'t':>2}  {'|λ|':>5}  {'c':>8}  {'lrcalc':>12}  {'kostka DP':>12}  {'ratio':>8}")
    print("-" * 60)

    for t in FACTORS:
        lam = scale(LAM0, t)
        mu = scale(MU0, t)
        nu = scale(NU0, t)

        val_lr, t_lr = bench_lrcalc(lam, mu, nu)
        val_k, t_k = bench_kostka(lam, mu, nu)

        if val_k is None:
            match = "TIMEOUT"
            ratio_s = "—"
        elif val_lr != val_k:
            match = "MISMATCH"
            ratio_s = "—"
        else:
            match = ""
            if t_lr > 0 and t_k is not None and t_k > 0:
                r = t_k / t_lr
                ratio_s = f"{r:.1f}x" if r >= 1 else f"1/{1/r:.1f}x"
            else:
                ratio_s = "~0"

        sz = sum(lam)
        print(
            f"{t:>2}  {sz:>5}  {val_lr:>8}  {fmt_time(t_lr):>12}  {fmt_time(t_k):>12}  {ratio_s:>8}  {match}"
        )


if __name__ == "__main__":
    main()
