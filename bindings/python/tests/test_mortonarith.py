import numpy as np
import pytest

import peclet.morton as ma


def ref_encode(coords, bits):
    dims = len(coords)
    n = len(coords[0])
    out = np.zeros(n, dtype=np.uint64)
    for b in range(bits):
        for d in range(dims):
            bit = ((coords[d].astype(np.uint64) >> np.uint64(b)) & np.uint64(1))
            out |= bit << np.uint64(b * dims + d)
    return out


@pytest.mark.parametrize("dims,bits", [(2, 32), (2, 16), (3, 21), (3, 16)])
def test_encode_decode_roundtrip(dims, bits):
    rng = np.random.default_rng(0)
    dtype = np.uint32 if bits > 16 else np.uint16
    hi = (1 << bits) - 1
    coords = [rng.integers(0, hi + 1, size=10000).astype(dtype) for _ in range(dims)]

    codes = ma.encode(*coords, bits=bits)
    assert codes.dtype == np.uint64
    np.testing.assert_array_equal(codes, ref_encode(coords, bits))

    back = ma.decode(codes, dims=dims, bits=bits)
    for d in range(dims):
        np.testing.assert_array_equal(back[d].astype(np.uint64), coords[d].astype(np.uint64))


def test_shift_matches_decode_modify_encode():
    rng = np.random.default_rng(1)
    x = rng.integers(0, 1 << 20, size=5000).astype(np.uint32)
    y = rng.integers(0, 1 << 20, size=5000).astype(np.uint32)
    codes = ma.encode(x, y, bits=32)

    shifted = ma.shift(codes, axis=0, delta=+5, dims=2, bits=32)
    xb, yb = ma.decode(shifted, dims=2, bits=32)
    np.testing.assert_array_equal(xb, (x + 5).astype(np.uint32))
    np.testing.assert_array_equal(yb, y)

    shifted = ma.shift(codes, axis=1, delta=-3, dims=2, bits=32)
    xb, yb = ma.decode(shifted, dims=2, bits=32)
    np.testing.assert_array_equal(xb, x)
    np.testing.assert_array_equal(yb, (y - 3).astype(np.uint32))


def test_shift_wraps():
    codes = ma.encode(np.array([255], np.uint16), np.array([0], np.uint16), bits=16)
    # 16-bit axis wraps mod 2^16, so 65535 + 1 -> 0
    big = ma.encode(np.array([65535], np.uint16), np.array([10], np.uint16), bits=16)
    wrapped = ma.shift(big, axis=0, delta=+1, dims=2, bits=16)
    x, y = ma.decode(wrapped, dims=2, bits=16)
    assert x[0] == 0 and y[0] == 10
    assert codes is not None


def test_box_zorder_matches_brute_force():
    lo = (3, 5)
    hi = (12, 9)
    bits = 16
    codes = ma.box_zorder(lo, hi, bits=bits)
    assert ma.box_count(lo, hi, bits=bits) == codes.size

    expect = []
    for yy in range(lo[1], hi[1] + 1):
        for xx in range(lo[0], hi[0] + 1):
            expect.append(int(ma.encode(np.array([xx], np.uint16),
                                        np.array([yy], np.uint16), bits=bits)[0]))
    expect.sort()
    assert np.all(np.diff(codes) > 0)  # strictly increasing (Z-order)
    np.testing.assert_array_equal(codes, np.array(expect, dtype=np.uint64))


def test_box_zorder_3d():
    lo = (1, 2, 0)
    hi = (5, 4, 3)
    codes = ma.box_zorder(lo, hi, bits=21)
    assert codes.size == (5 - 1 + 1) * (4 - 2 + 1) * (3 - 0 + 1)
    assert np.all(np.diff(codes) > 0)


def test_unsupported_config_raises():
    with pytest.raises(ValueError):
        ma.encode(np.array([1], np.uint32), np.array([2], np.uint32), bits=24)
