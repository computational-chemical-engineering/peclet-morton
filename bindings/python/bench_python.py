"""Python-side benchmark: vectorised Morton arithmetic vs the decode/encode
round trip, both over NumPy arrays. Run after building the native library:

    cmake -S . -B build && cmake --build build --target mortonarith_c
    PYTHONPATH=bindings/python python3 bindings/python/bench_python.py
"""
import time

import numpy as np

import mortonarith as ma

N = 1 << 22
rng = np.random.default_rng(0)
x = rng.integers(0, 1 << 30, size=N).astype(np.uint32)
y = rng.integers(0, 1 << 30, size=N).astype(np.uint32)
codes = ma.encode(x, y, bits=32)


def timed(label, fn, reps=5):
    fn()  # warmup
    best = min((lambda: (t0 := time.perf_counter(), fn(), time.perf_counter())[::2])()
               for _ in range(reps))
    secs = best[1] - best[0]
    print(f"  {label:<44} {N / secs / 1e6:8.1f} Mops/s")


print(f"== Python / NumPy, N={N} ==")
timed("encode", lambda: ma.encode(x, y, bits=32))
timed("decode", lambda: ma.decode(codes, dims=2, bits=32))
timed("shift +1 in x (arithmetic, no decode)",
      lambda: ma.shift(codes, axis=0, delta=1, dims=2, bits=32))


def roundtrip_shift():
    xx, yy = ma.decode(codes, dims=2, bits=32)
    return ma.encode((xx + 1).astype(np.uint32), yy, bits=32)


timed("shift +1 in x (decode + encode)", roundtrip_shift)
