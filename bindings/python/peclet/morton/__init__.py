"""peclet.morton - fast Morton (Z-order) codes with arithmetic, for NumPy.

This is a thin ``ctypes`` wrapper over the C++ ``morton`` library. Every
operation is vectorised: it runs over whole NumPy arrays in compiled code, so
there is no per-element Python overhead.

Supported configurations (dimensions x bits-per-axis):

    2D: 32, 16        3D: 21, 16

Codes are always returned as ``uint64``.

Example
-------
>>> import numpy as np
>>> from peclet.morton import encode, decode, shift
>>> x = np.array([1, 2, 3], dtype=np.uint32)
>>> y = np.array([4, 5, 6], dtype=np.uint32)
>>> codes = encode(x, y, bits=32)
>>> xb = shift(codes, axis=0, delta=+1, dims=2, bits=32)  # move +1 in x, no decode
>>> decode(xb, dims=2, bits=32)
(array([2, 3, 4], dtype=uint32), array([4, 5, 6], dtype=uint32))
"""

from __future__ import annotations

import ctypes
import glob
import os

import numpy as np

__all__ = ["encode", "decode", "shift", "box_zorder", "box_count"]

_HERE = os.path.dirname(os.path.abspath(__file__))


def _load_library() -> ctypes.CDLL:
    candidates = glob.glob(os.path.join(_HERE, "libmortonarith_c*"))
    if not candidates:
        raise ImportError(
            "mortonarith native library not found. Build it first:\n"
            "    cmake -S . -B build && cmake --build build --target mortonarith_c"
        )
    return ctypes.CDLL(candidates[0])


_lib = _load_library()

# (dims, bits) -> (numpy coord dtype, C suffix)
_CONFIG = {
    (2, 32): (np.uint32, "u32"),
    (2, 16): (np.uint16, "u16"),
    (3, 21): (np.uint32, "u32"),
    (3, 16): (np.uint16, "u16"),
}

_u64 = ctypes.c_uint64
_sz = ctypes.c_size_t


def _ptr(arr):
    return arr.ctypes.data_as(ctypes.c_void_p)


def _cfg(dims: int, bits: int):
    key = (dims, bits)
    if key not in _CONFIG:
        raise ValueError(
            f"unsupported configuration dims={dims}, bits={bits}; "
            f"supported: {sorted(_CONFIG)}"
        )
    return _CONFIG[key]


def encode(*coords: np.ndarray, bits: int) -> np.ndarray:
    """Interleave coordinate arrays into Morton codes (uint64).

    Pass 2 or 3 equally sized integer arrays (one per axis).
    """
    dims = len(coords)
    dtype, suffix = _cfg(dims, bits)
    arrs = [np.ascontiguousarray(c, dtype=dtype) for c in coords]
    n = arrs[0].size
    for a in arrs:
        if a.size != n:
            raise ValueError("coordinate arrays must have the same length")
    out = np.empty(n, dtype=np.uint64)
    fn = getattr(_lib, f"mortonarith_encode{dims}_{suffix}")
    args = [_ptr(a) for a in arrs] + [_ptr(out), _sz(n)]
    fn(*args)
    return out


def decode(codes: np.ndarray, *, dims: int, bits: int) -> tuple:
    """Decode Morton codes back into a tuple of `dims` coordinate arrays."""
    dtype, suffix = _cfg(dims, bits)
    codes = np.ascontiguousarray(codes, dtype=np.uint64)
    n = codes.size
    outs = [np.empty(n, dtype=dtype) for _ in range(dims)]
    fn = getattr(_lib, f"mortonarith_decode{dims}_{suffix}")
    fn(_ptr(codes), *[_ptr(o) for o in outs], _sz(n))
    return tuple(outs)


def shift(codes: np.ndarray, *, axis: int, delta: int, dims: int, bits: int) -> np.ndarray:
    """Add `delta` to one axis of each code, in Morton space (no decode/encode).

    `delta` may be negative. Coordinates wrap modulo 2**bits.
    """
    _, suffix = _cfg(dims, bits)
    codes = np.ascontiguousarray(codes, dtype=np.uint64)
    out = np.empty_like(codes)
    fn = getattr(_lib, f"mortonarith_add{dims}_{suffix}")
    fn(_ptr(codes), _ptr(out), _sz(codes.size), ctypes.c_uint(axis), ctypes.c_int64(delta))
    return out


def box_count(lo, hi, *, bits: int) -> int:
    """Number of cells in the inclusive box [lo, hi]."""
    dims = len(lo)
    dtype, suffix = _cfg(dims, bits)
    lo = np.ascontiguousarray(lo, dtype=dtype)
    hi = np.ascontiguousarray(hi, dtype=dtype)
    fn = getattr(_lib, f"mortonarith_box_count{dims}_{suffix}")
    fn.restype = _u64
    return int(fn(_ptr(lo), _ptr(hi)))


def box_zorder(lo, hi, *, bits: int) -> np.ndarray:
    """All Morton codes in the inclusive box [lo, hi], sorted in Z-order."""
    dims = len(lo)
    dtype, suffix = _cfg(dims, bits)
    lo = np.ascontiguousarray(lo, dtype=dtype)
    hi = np.ascontiguousarray(hi, dtype=dtype)
    n = box_count(lo, hi, bits=bits)
    out = np.empty(n, dtype=np.uint64)
    if n:
        fn = getattr(_lib, f"mortonarith_box_zorder{dims}_{suffix}")
        fn(_ptr(lo), _ptr(hi), _ptr(out))
    return out
