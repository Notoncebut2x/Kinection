"""
R2-backed EIGENSTRAT PACKGENO reader.

Reads SNP rows from a GENO file stored in Cloudflare R2 using HTTP range
requests — one request per chunk — instead of local memory-mapping.

GENO binary format (EIGENSOFT PACKGENO):
  Header:  64-byte ASCII text — "GENO   {n_indiv} {n_snps} {hash1} {hash2}\\0..."
  Body:    n_snps rows × ceil(n_indiv/4) bytes per row
  Packing: 2 bits per genotype, LSB first within each byte:
             bits 1-0  → individual j*4 + 0
             bits 3-2  → individual j*4 + 1
             bits 5-4  → individual j*4 + 2
             bits 7-6  → individual j*4 + 3
  Values:  0=hom_ref  1=het  2=hom_alt  3=missing

Range request strategy:
  For a sorted array of geno_indices spanning [min, max], we fetch
  bytes [header + min*bps, header + (max+1)*bps) in one request and
  extract the requested rows from the buffer. At 32% SNP overlap density
  this reads at most ~3× the data we need per chunk — acceptable given
  R2 has no egress fees.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np

from . import r2_client

log = logging.getLogger(__name__)

# Bytes fetched to parse the header string; the real data offset is the
# record length, padded to >= 48 bytes (see packed_header_size).
GENO_HEADER_READ = 64


def packed_header_size(bytes_per_record: int) -> int:
    return max(48, bytes_per_record)


@dataclass
class R2GenoFile:
    """
    Reads EIGENSTRAT PACKGENO files from Cloudflare R2 via byte-range requests.

    Drop-in replacement for the local GenoFile class. Implements the same
    interface: open(), read_snp_row(), read_chunk(), close().

    Usage:
        geno = R2GenoFile.open(r2_client.GENO_KEY)
        rows = geno.read_chunk(chunk_geno_indices)   # (n, n_indiv) int8
        geno.close()
    """
    r2_key: str
    n_indiv: int
    n_snps: int
    bytes_per_snp: int
    header_size: int = 48

    @classmethod
    def open(cls, r2_key: str) -> 'R2GenoFile':
        """Fetch the header from R2 and return a configured instance."""
        header = r2_client.get_object_bytes(r2_key, byte_range=(0, GENO_HEADER_READ - 1))
        header_str = header.split(b'\x00')[0].decode('ascii', errors='replace')
        parts = header_str.split()
        if not parts or parts[0] != 'GENO':
            raise ValueError(f"Not a PACKGENO file (header={header_str!r}): {r2_key}")
        n_indiv       = int(parts[1])
        n_snps        = int(parts[2])
        bytes_per_snp = (n_indiv + 3) // 4
        header_size   = packed_header_size(bytes_per_snp)
        log.info(
            'R2GenoFile: %d individuals × %d SNPs, %d bytes/SNP, %d-byte header  [key=%s]',
            n_indiv, n_snps, bytes_per_snp, header_size, r2_key,
        )
        return cls(r2_key=r2_key, n_indiv=n_indiv, n_snps=n_snps,
                   bytes_per_snp=bytes_per_snp, header_size=header_size)

    def read_snp_row(self, snp_index: int) -> np.ndarray:
        """Read a single SNP row. Prefer read_chunk() for bulk reads."""
        return self.read_chunk(np.array([snp_index], dtype=np.int32))[0]

    def read_chunk(self, geno_indices: np.ndarray | list) -> np.ndarray:
        """
        Fetch a set of SNP rows from R2 in one HTTP range request.

        geno_indices must be sorted ascending. Fetches the byte span
        [min_idx, max_idx] inclusive, then extracts the requested rows.

        Returns: (n, n_indiv) int8 — 0=hom_ref 1=het 2=hom_alt 3=missing
        """
        geno_indices = np.asarray(geno_indices, dtype=np.int32)
        if len(geno_indices) == 0:
            return np.empty((0, self.n_indiv), dtype=np.int8)

        min_idx = int(geno_indices[0])
        max_idx = int(geno_indices[-1])
        bps     = self.bytes_per_snp

        byte_start = self.header_size + min_idx * bps
        byte_end   = self.header_size + (max_idx + 1) * bps - 1

        raw = r2_client.get_object_bytes(self.r2_key, byte_range=(byte_start, byte_end))
        buf = np.frombuffer(raw, dtype=np.uint8)

        result = np.empty((len(geno_indices), self.n_indiv), dtype=np.int8)
        for out_i, g_idx in enumerate(geno_indices):
            row_start = (int(g_idx) - min_idx) * bps
            row_bytes = buf[row_start : row_start + bps]
            # MSB-first: first individual of each group is in the high bits.
            gt = np.empty(len(row_bytes) * 4, dtype=np.int8)
            gt[0::4] = (row_bytes >> 6)  & 0x03
            gt[1::4] = (row_bytes >> 4)  & 0x03
            gt[2::4] = (row_bytes >> 2)  & 0x03
            gt[3::4] =  row_bytes        & 0x03
            result[out_i] = gt[: self.n_indiv]

        return result

    def close(self) -> None:
        pass  # no persistent connection to close

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.close()
