"""
Convert an AADR 'tgeno' (transpose_packed) .geno file to the older
non-transposed PACKEDANCESTRYMAP ('GENO') format that the pipeline reads.

Why: AADR v66+ ships the 1240K genotypes in the newer transposed layout
(one record per individual). The analysis pipeline reads SNP rows, which are
scattered across the whole file in the transposed layout (unworkable over R2
range requests). This tool rewrites the matrix SNP-major.

Format facts (verified against AADR v62/v66 + AdmixTools mcio.c):
  tgeno  : 48-byte header "TGENO n m h h"; then n_indiv rows of ceil(n_snp/4)
           bytes; 4 genotypes/byte, MSB-first (first element in bits 6-7).
  packed : header max(48, ceil(n_indiv/4)) bytes "GENO n m h h"; then n_snp
           rows of ceil(n_indiv/4) bytes; same MSB-first packing.
  values : 0=hom-ref 1=het 2=hom-alt 3=missing.

Usage:
  python scripts/convert_tgeno_to_packed.py IN.tgeno.geno OUT.geno
  python scripts/convert_tgeno_to_packed.py --self-test
"""
from __future__ import annotations

import argparse
import mmap
import sys
import time
from pathlib import Path

import numpy as np

TGENO_HEADER = 48
BLOCK_SNPS = 8192          # multiple of 4; ~200 MB working set at n_indiv~23k


def _ceil_div4(n: int) -> int:
    return (n + 3) // 4


def unpack_msb(cols: np.ndarray) -> np.ndarray:
    """(rows, nbytes) uint8 -> (rows, nbytes*4) int8 genotypes, MSB-first."""
    g = np.empty((cols.shape[0], cols.shape[1] * 4), dtype=np.int8)
    g[:, 0::4] = (cols >> 6) & 3
    g[:, 1::4] = (cols >> 4) & 3
    g[:, 2::4] = (cols >> 2) & 3
    g[:, 3::4] = cols & 3
    return g


def pack_msb(gt: np.ndarray, rlen: int) -> np.ndarray:
    """(rows, n_elem) int8 genotypes -> (rows, rlen) uint8, MSB-first packing
    of elements along axis 1 (padding to rlen*4 with 0)."""
    rows, n = gt.shape
    padded = np.zeros((rows, rlen * 4), dtype=np.uint8)
    padded[:, :n] = gt.astype(np.uint8)
    q = padded.reshape(rows, rlen, 4)
    out = (q[:, :, 0] << 6) | (q[:, :, 1] << 4) | (q[:, :, 2] << 2) | q[:, :, 3]
    return out.astype(np.uint8)


def read_tgeno_header(mm: mmap.mmap):
    raw = mm[:64].split(b"\x00")[0].decode("ascii", errors="replace")
    parts = raw.split()
    if not parts or parts[0] != "TGENO":
        raise ValueError(f"Not a TGENO file (header={raw!r})")
    n_indiv, n_snp = int(parts[1]), int(parts[2])
    h1 = parts[3] if len(parts) > 3 else "0"
    h2 = parts[4] if len(parts) > 4 else "0"
    return n_indiv, n_snp, h1, h2


def convert(in_path: Path, out_path: Path) -> None:
    fh = open(in_path, "rb")
    mm = mmap.mmap(fh.fileno(), 0, access=mmap.ACCESS_READ)
    n_indiv, n_snp, h1, h2 = read_tgeno_header(mm)
    trlen = _ceil_div4(n_snp)                 # bytes per individual row (input)
    rlen = _ceil_div4(n_indiv)                # bytes per SNP row (output)
    out_header = max(48, rlen)

    expected_in = TGENO_HEADER + n_indiv * trlen
    actual_in = mm.size()
    print(f"input : {in_path.name}  TGENO {n_indiv} indiv x {n_snp} snp  "
          f"trlen={trlen}  size={actual_in:,} (expected {expected_in:,})")
    if actual_in < expected_in:
        raise ValueError("Input file is truncated.")

    # Individual-major view of the body: body[i, b] = byte b of individual i.
    body = np.frombuffer(mm, dtype=np.uint8, count=n_indiv * trlen,
                         offset=TGENO_HEADER).reshape(n_indiv, trlen)

    print(f"output: {out_path.name}  GENO  rlen={rlen}  header={out_header}  "
          f"est size={out_header + n_snp * rlen:,}")

    hdr = f"GENO {n_indiv} {n_snp} {h1} {h2}".encode("ascii")
    hdr = hdr + b"\x00" * (out_header - len(hdr))

    t0 = time.time()
    with open(out_path, "wb") as out:
        out.write(hdr)
        for s0 in range(0, n_snp, BLOCK_SNPS):
            s1 = min(s0 + BLOCK_SNPS, n_snp)
            b_lo = s0 // 4
            b_hi = _ceil_div4(s1)
            cols = body[:, b_lo:b_hi]                     # (n_indiv, nbytes)
            g = unpack_msb(cols)[:, : (s1 - s0)]          # (n_indiv, block)
            packed = pack_msb(g.T, rlen)                  # (block, rlen)
            out.write(packed.tobytes())
            del cols, g, packed
            if (s0 // BLOCK_SNPS) % 20 == 0:
                pct = 100 * s1 / n_snp
                print(f"  {s1:>9,}/{n_snp:,} SNPs ({pct:5.1f}%)  "
                      f"{time.time()-t0:6.1f}s", flush=True)
    del body                       # release the mmap view before closing
    mm.close(); fh.close()
    size = out_path.stat().st_size
    print(f"done: wrote {size:,} bytes in {time.time()-t0:.1f}s "
          f"(expected {out_header + n_snp * rlen:,})")


def self_test() -> None:
    """Round-trip a random matrix: build a tgeno file, convert, read with the
    pipeline's GenoFile, and assert genotypes match exactly."""
    sys.path.insert(0, str(Path(__file__).parent))
    from utils.parsers import GenoFile

    rng = np.random.default_rng(0)
    n_indiv, n_snp = 37, 9001               # deliberately non-multiples of 4
    truth = rng.integers(0, 4, size=(n_indiv, n_snp), dtype=np.int8)  # incl 3=missing
    trlen = _ceil_div4(n_snp)

    # Encode tgeno: per individual row, pack SNPs MSB-first.
    body = np.zeros((n_indiv, trlen), dtype=np.uint8)
    for i in range(n_indiv):
        body[i] = pack_msb(truth[i][None, :], trlen)[0]
    import tempfile, os
    d = Path(tempfile.mkdtemp())
    tg = d / "t.tgeno"; gg = d / "t.geno"
    hdr = f"TGENO {n_indiv} {n_snp} 0 0".encode(); hdr += b"\x00" * (TGENO_HEADER - len(hdr))
    with open(tg, "wb") as f:
        f.write(hdr); f.write(body.tobytes())

    convert(tg, gg)
    g = GenoFile.open(gg)
    assert g.n_indiv == n_indiv and g.n_snps == n_snp, (g.n_indiv, g.n_snps)
    ok = True
    for s in range(n_snp):
        row = g.read_snp_row(s)
        if not np.array_equal(row, truth[:, s]):
            ok = False
            print("MISMATCH at snp", s); break
    g.close()
    for p in (tg, gg): os.remove(p)
    os.rmdir(d)
    print("SELF-TEST:", "PASS ✅" if ok else "FAIL ❌")
    if not ok:
        sys.exit(1)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("input", nargs="?", help="input .tgeno/.geno (TGENO magic)")
    ap.add_argument("output", nargs="?", help="output .geno (GENO magic)")
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args()
    if a.self_test:
        self_test()
    elif a.input and a.output:
        convert(Path(a.input), Path(a.output))
    else:
        ap.error("provide input and output, or --self-test")
