# This is a module to test the code.
# How to run: pytest -q
#
# Adjust this import to point to your Grafic implementation:
# from your_module import Grafic
from grafic import Grafic
import os
import numpy as np
import pytest


def make_float_data(shape, seed=1234):
    rng = np.random.default_rng(seed)
    return rng.random(shape, dtype=np.float32) * 100 - 50  # some negative/positive values


def make_int64_data(shape, seed=5678):
    rng = np.random.default_rng(seed)
    return rng.integers(low=-2**31, high=2**31 - 1, size=shape, dtype=np.int64)


def make_int32_data(shape, seed=9012):
    rng = np.random.default_rng(seed)
    return rng.integers(low=-2**31, high=2**31 - 1, size=shape, dtype=np.int32)


def headers_close(h1, h2, rtol=0, atol=1e-6):
    # h = [n1, n2, n3, dx, x1, x2, x3, f1, f2, f3, f4]
    if int(h1[0]) != int(h2[0]):
        return False
    if int(h1[1]) != int(h2[1]):
        return False
    if int(h1[2]) != int(h2[2]):
        return False
    a1 = np.array(h1[3:], dtype=float)
    a2 = np.array(h2[3:], dtype=float)
    return np.allclose(a1, a2, rtol=rtol, atol=atol)


@pytest.mark.parametrize("layout", ["sliced", "single"])
def test_roundtrip_float32_little_endian(tmp_path, layout):
    shape = (7, 5, 3)
    arr = make_float_data(shape)
    out = tmp_path / f"le_{layout}.dat"

    g = Grafic(endian="<")
    g.set_data(arr)
    g.make_header(box_size_cu=1.0, offsets=(0.1, 0.2, 0.3), extras=(1.0, 2.0, 3.0, 4.0))
    orig_header = list(g.header)
    g.write_float(str(out), layout=layout)

    # Check record marker endianness: first 4 bytes = 44 as little-endian
    with open(out, "rb") as f:
        first4 = f.read(4)
        assert int.from_bytes(first4, "little") == 44
        assert int.from_bytes(first4, "big") != 44

    g2 = Grafic(endian="<")
    g2.read(str(out))
    assert g2.data.shape == shape
    assert np.allclose(g2.data, arr, rtol=0, atol=0)
    assert headers_close(orig_header, g2.header)


@pytest.mark.parametrize("layout", ["sliced", "single"])
def test_roundtrip_int64_big_endian(tmp_path, layout):
    shape = (6, 4, 2)
    arr = make_int64_data(shape)
    out = tmp_path / f"be_{layout}.dat"

    g = Grafic(endian=">")
    g.set_data(arr)
    g.make_header(box_size_cu=2.0, offsets=(10.0, -5.0, 3.5), extras=(0.0, 0.0, 0.0, 0.0))
    orig_header = list(g.header)
    g.write(str(out), dtype_out=np.int64, layout=layout)

    # Check record marker endianness: first 4 bytes = 44 as big-endian
    with open(out, "rb") as f:
        first4 = f.read(4)
        assert int.from_bytes(first4, "big") == 44
        assert int.from_bytes(first4, "little") != 44

    g2 = Grafic(endian=">")
    g2.read(str(out))
    assert g2.data.dtype == np.int64
    assert np.array_equal(g2.data, arr)
    assert headers_close(orig_header, g2.header)


def test_read_header_only(tmp_path):
    shape = (5, 5, 5)
    arr = make_float_data(shape, seed=42)
    out = tmp_path / "header_only.dat"

    g = Grafic(endian="<")
    g.set_data(arr)
    g.make_header(box_size_cu=5.0, offsets=(1.0, 2.0, 3.0), extras=(4.0, 5.0, 6.0, 7.0))
    orig_header = list(g.header)
    g.write_float(str(out), layout="sliced")

    g2 = Grafic(endian="<")
    hdr = g2.read_header_only(str(out))
    assert headers_close(orig_header, hdr)
    # data should still be None
    assert g2.data is None


def test_int32_requires_explicit_dtype(tmp_path):
    shape = (4, 3, 2)
    arr = make_int32_data(shape)
    out = tmp_path / "int32_sliced.dat"

    g = Grafic(endian="<")
    g.set_data(arr)
    g.make_header(box_size_cu=1.0)
    g.write(str(out), dtype_out=np.int32, layout="sliced")

    # Default read (dtype inference) will assume float32 for 4-byte elements
    g_default = Grafic(endian="<")
    g_default.read(str(out))
    assert g_default.data.dtype == np.float32  # default assumption

    # Explicit dtype read should match exactly
    g_int = Grafic(endian="<")
    g_int.read(str(out), dtype=np.int32)
    assert g_int.data.dtype == np.int32
    assert np.array_equal(g_int.data, arr)


def test_wrong_endianness_read_either_raises_or_autodetects(tmp_path):
    shape = (8, 6, 4)
    arr = make_float_data(shape)
    out = tmp_path / "le_file.dat"

    # Write little-endian file
    g = Grafic(endian="<")
    g.set_data(arr)
    g.make_header(box_size_cu=3.0)
    g.write_float(str(out), layout="sliced")

    # Try reading as big-endian. Depending on implementation this may either raise
    # or succeed if the reader auto-detects endianness. Accept both behaviors.
    g_wrong = Grafic(endian=">")
    try:
        g_wrong.read(str(out))
    except Exception:
        return
    # If no exception, data should be correct
    assert np.allclose(g_wrong.data, arr)


def test_no_temp_files_left(tmp_path):
    shape = (3, 3, 3)
    arr = make_float_data(shape)
    out = tmp_path / "clean.dat"

    g = Grafic(endian="<")
    g.set_data(arr)
    g.make_header(box_size_cu=1.0)
    g.write_float(str(out), layout="sliced")

    # Ensure only the target file exists; no temp files with our prefix
    leftovers = [p for p in os.listdir(tmp_path) if str(p).startswith(".grafic_tmp_")]
    assert leftovers == []

    # File is readable and correct
    g2 = Grafic(endian="<")
    g2.read(str(out))
    assert np.allclose(g2.data, arr)


def test_layout_detection_both_modes(tmp_path):
    shape = (5, 4, 3)
    arr = make_float_data(shape)

    # Sliced
    out_sliced = tmp_path / "sliced.dat"
    gs = Grafic(endian="<")
    gs.set_data(arr)
    gs.make_header(box_size_cu=1.0)
    gs.write_float(str(out_sliced), layout="sliced")
    rs = Grafic(endian="<")
    rs.read(str(out_sliced))
    assert np.allclose(rs.data, arr)

    # Single
    out_single = tmp_path / "single.dat"
    g1 = Grafic(endian="<")
    g1.set_data(arr)
    g1.make_header(box_size_cu=1.0)
    g1.write_float(str(out_single), layout="single")
    r1 = Grafic(endian="<")
    r1.read(str(out_single))
    assert np.allclose(r1.data, arr)


def test_header_content_values(tmp_path):
    shape = (9, 7, 5)
    arr = make_float_data(shape)
    out = tmp_path / "hdr_content.dat"

    offsets = (12.5, -7.25, 3.0)
    extras = (0.125, 256.0, -1.5, 42.0)
    box_size = 18.0

    g = Grafic(endian=">")
    g.set_data(arr)
    g.make_header(box_size_cu=box_size, offsets=offsets, extras=extras)
    g.write_float(str(out), layout="sliced")

    r = Grafic(endian=">")
    r.read(str(out))
    n1, n2, n3, dx, x1, x2, x3, f1, f2, f3, f4 = r.header

    assert (n1, n2, n3) == shape
    # dx is computed as box_size / n1
    assert np.isclose(dx, box_size / shape[0], rtol=0, atol=1e-6)
    assert np.allclose([x1, x2, x3], offsets, rtol=0, atol=1e-6)
    assert np.allclose([f1, f2, f3, f4], extras, rtol=0, atol=1e-6)


def test_truncated_file_errors(tmp_path):
    shape = (4, 4, 2)
    arr = make_float_data(shape)
    out = tmp_path / "trunc.dat"

    g = Grafic(endian="<")
    g.set_data(arr)
    g.make_header(box_size_cu=1.0)
    g.write_float(str(out), layout="sliced")

    # Truncate the file to cut in the middle of a record
    with open(out, "rb") as f:
        data = f.read()
    with open(out, "wb") as f:
        f.write(data[:-3])  # remove last 3 bytes

    r = Grafic(endian="<")
    with pytest.raises(Exception):
        r.read(str(out))
