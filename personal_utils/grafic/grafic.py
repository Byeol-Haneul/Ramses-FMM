import numpy as np
import struct
import os
import tempfile

class Grafic:
  """
  Read/write GRAFIC-style Fortran sequential unformatted files.
  
  Endianness:
    - Use endian="<" (little, default) or ">" (big). Header, record markers,
      and data are written with the same endianness.
  Layout:
    - "sliced": header + n3 records, each (n1 x n2) in Fortran order.
    - "single": header + 1 record, full (n1 x n2 x n3) in Fortran order.
  
  After read():
    - self.header = [n1, n2, n3, dx, xoff1, xoff2, xoff3, f1, f2, f3, f4]
    - self.data   = ndarray (n1, n2, n3), native endianness.
  """
  
  def __init__(self, endian: str = "<"):
      if endian not in ("<", ">"):
          raise ValueError("endian must be '<' (little) or '>' (big)")
      self.endian = endian
      self._rec_fmt = endian + "i"        # 32-bit record marker
      self._head_fmt = endian + "3i8f"    # header layout (3 ints, 8 floats)
  
      self.data = None
      self.header = None
  
  # ---------- Header utilities ----------
  
  def make_header(self, box_size_cu, offsets=(0.0, 0.0, 0.0), extras=(0.0, 0.0, 0.0, 0.0)):
      """
      Build header from current data and physical box size in code units.
  
      Assumes a cubic grid and isotropic spacing dx = box_size_cu / n1.
      """
      if self.data is None:
          raise ValueError("No data array set for Grafic object")
      if self.data.ndim != 3:
          raise ValueError("data must be 3D")
  
      n1, n2, n3 = map(int, self.data.shape)
      if not (n1 == n2 == n3):
          # Still allow, but dx is based on n1 to match legacy behavior
          pass
  
      dx = float(box_size_cu) / n1
      x1, x2, x3 = map(float, offsets)
      f1, f2, f3, f4 = map(float, extras)
      # f1 = float(box_size_cu)
      self.header = [np.int32(n1), np.int32(n2), np.int32(n3),
                     np.float32(dx), np.float32(x1), np.float32(x2), np.float32(x3),
                     np.float32(f1), np.float32(f2), np.float32(f3), np.float32(f4)]
  
  def _pack_header(self):
      if self.header is None:
          raise ValueError("Header is not set")
      # Unpack to Python types in the expected order for struct
      n1, n2, n3, dx, x1, x2, x3, f1, f2, f3, f4 = self.header
      return struct.pack(self._head_fmt, int(n1), int(n2), int(n3),
                         float(dx), float(x1), float(x2), float(x3),
                         float(f1), float(f2), float(f3), float(f4))
  
  # ---------- File I/O primitives ----------
  
  def _write_fortran_record(self, fh, payload_bytes: bytes):
      n = len(payload_bytes)
      fh.write(struct.pack(self._rec_fmt, n))
      fh.write(payload_bytes)
      fh.write(struct.pack(self._rec_fmt, n))
  
  def _read_fortran_record(self, fh):
      hdr = fh.read(4)
      if len(hdr) == 0:
          raise EOFError("Unexpected EOF while reading record marker")
      if len(hdr) != 4:
          raise IOError("Truncated record marker")
      n = struct.unpack(self._rec_fmt, hdr)[0]
      payload = fh.read(n)
      if len(payload) != n:
          raise IOError("Truncated record payload")
      tail = fh.read(4)
      if len(tail) != 4:
          raise IOError("Truncated trailing record marker")
      n2 = struct.unpack(self._rec_fmt, tail)[0]
      if n2 != n:
          raise IOError("Fortran record marker mismatch")
      return payload, n
  
  # ---------- Public API ----------
  
  def set_data(self, arr):
      arr = np.asarray(arr)
      if arr.ndim != 3:
          raise ValueError("data must be a 3D array")
      self.data = arr
  
  def write(self, output_name, dtype_out=np.float32, layout="sliced"):
      """
      Write data with dtype_out and layout ("sliced" or "single").
      """
      if self.data is None:
          raise ValueError("No data to write")
      if self.header is None:
          raise ValueError("Header not set; call make_header() first")
      if layout not in ("sliced", "single"):
          raise ValueError("layout must be 'sliced' or 'single'")
  
      # Prepare destination temp file in the same directory for atomic replace
      out_dir = os.path.dirname(os.path.abspath(output_name)) or "."
      fd, tmp_path = tempfile.mkstemp(prefix=".grafic_tmp_", dir=out_dir)
      os.close(fd)  # we'll reopen with buffered IO
  
      dtype_file = np.dtype(dtype_out).newbyteorder(self.endian)
  
      try:
          with open(tmp_path, "wb") as f:
              # Header
              self._write_fortran_record(f, self._pack_header())
  
              # Data
              if layout == "single":
                  arrF = np.asfortranarray(self.data)
                  if arrF.dtype != dtype_file:
                      arrF = arrF.astype(dtype_file, copy=False)
                  self._write_fortran_record(f, arrF.tobytes(order="F"))
              else:
                  # RAMSES/GRAFIC often expect z-slices (k index) as separate records
                  n1, n2, n3 = self.data.shape
                  for k in range(n3):
                      sli = np.asfortranarray(self.data[:, :, k])
                      if sli.dtype != dtype_file:
                          sli = sli.astype(dtype_file, copy=False)
                      self._write_fortran_record(f, sli.tobytes(order="F"))
  
          # Atomic replace
          os.replace(tmp_path, output_name)
      except Exception:
          # Clean up temp file on error
          try:
              os.remove(tmp_path)
          except OSError:
              pass
          raise
  
  def write_float(self, output_name, layout="sliced"):
      self.write(output_name, dtype_out=np.float32, layout=layout)
  
  def write_int64(self, output_name, layout="sliced"):
      self.write(output_name, dtype_out=np.int64, layout=layout)
  
  def read(self, filename, dtype=None):
    """
    Read a GRAFIC file (supports both 'sliced' and 'single' layouts).

    dtype:
      - None: infer from first data record size.
               4 bytes/elem -> float32 (ambiguous int32/float32; default to float32)
               8 bytes/elem -> int64
      - Or explicitly set, e.g., np.float32, np.int32, np.int64.

    Uses the instance's endianness (self.endian: '<' or '>') and header fmt (self._head_fmt).
    Requires self._read_fortran_record to be defined: returns (bytes, length).
    """
    import struct
    import numpy as np

    with open(filename, "rb") as f:
        # Header
        head_bytes, _ = self._read_fortran_record(f)
        try:
            n1, n2, n3, dx, x1, x2, x3, f1, f2, f3, f4 = struct.unpack(self._head_fmt, head_bytes)
        except struct.error as e:
            raise ValueError(f"Header unpack failed: {e}")
        self.header = [n1, n2, n3, dx, x1, x2, x3, f1, f2, f3, f4]

        # Peek first data record
        first_bytes, first_len = self._read_fortran_record(f)

        # Infer dtype and layout if needed
        layout = None

        def match_len(elem_bytes):
            if first_len == n1 * n2 * elem_bytes:
                return "sliced"
            if first_len == n1 * n2 * n3 * elem_bytes:
                return "single"
            return None

        if dtype is None:
            # Try 8-byte elements first (int64), then 4-byte (float32 default).
            layout = match_len(8)
            if layout is not None:
                base_dtype = np.int64
                itemsize = 8
            else:
                layout = match_len(4)
                if layout is not None:
                    base_dtype = np.float32  # default for 4-byte ambiguity
                    itemsize = 4
                else:
                    raise ValueError(
                        f"Cannot infer dtype/layout: first record len {first_len} "
                        f"not matching grid ({n1},{n2},{n3})."
                    )
        else:
            dt = np.dtype(dtype)
            base_dtype = dt.type
            itemsize = dt.itemsize
            layout = match_len(itemsize)
            if layout is None:
                raise ValueError(
                    f"Record length {first_len} inconsistent with dtype {dt} and grid ({n1},{n2},{n3})."
                )

        dtype_file = np.dtype(base_dtype).newbyteorder(self.endian)
        dtype_native = np.dtype(base_dtype).newbyteorder('=')

        if layout == "single":
            arrF = np.frombuffer(first_bytes, dtype=dtype_file).reshape((n1, n2, n3), order="F")
            arr = arrF.astype(dtype_native, copy=False)
        else:
            # Sliced: first record is one (n1, n2) slab in Fortran order
            slab0 = np.frombuffer(first_bytes, dtype=dtype_file).reshape((n1, n2), order="F")
            # Allocate native-endian output (n1, n2, n3)
            arr = np.empty((n1, n2, n3), dtype=dtype_native)
            arr[:, :, 0] = slab0.astype(dtype_native, copy=False)
            for k in range(1, n3):
                rec_bytes, rec_len = self._read_fortran_record(f)
                if rec_len != first_len:
                    raise ValueError(f"Inconsistent slice record length at k={k}: {rec_len} != {first_len}")
                slab = np.frombuffer(rec_bytes, dtype=dtype_file).reshape((n1, n2), order="F")
                arr[:, :, k] = slab.astype(dtype_native, copy=False)

    self.data = arr
    return arr
    
  def read_header_only(self, filename):
      """
      Read only the header of a GRAFIC file.
      """
      import struct
  
      with open(filename, "rb") as f:
          head_bytes, _ = self._read_fortran_record(f)
          try:
              n1, n2, n3, dx, x1, x2, x3, f1, f2, f3, f4 = struct.unpack(self._head_fmt, head_bytes)
          except struct.error as e:
              raise ValueError(f"Header unpack failed: {e}")
          self.header = [n1, n2, n3, dx, x1, x2, x3, f1, f2, f3, f4]
      return self.header
