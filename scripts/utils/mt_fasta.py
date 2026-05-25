"""
Parser for AADR's compressed mitochondrial DNA repository.

The mt repo format (Mallick / Patterson "cTools"-style compression) is:

    >ID1\\tcompressed_seq1
    >ID2\\tcompressed_seq2
    ...

Compressed-sequence character meanings:
    Q       same as the rCRS reference at this position
    A C G T actual base call (differs from rCRS)
    n       no data (missing)
    -       padding (when the sample's mt is shorter than rCRS)

To reconstruct a sample's mt sequence, for each position i:
    - if compressed[i] == 'Q' → use rcrs[i]
    - else → use compressed[i]
Then pad with '-' to the full rCRS length (16,569).

Our parser returns sequences as plain Python strings, length 16,569,
where 'n' and '-' both mean "missing data at this position" and any of
A/C/G/T means a definite base call. Comparing two such strings at a
given position is a straightforward equality check (skip if either is
missing).

The rCRS reference itself is bundled in `scripts/data/mt_rcrs.txt`
(extracted from AADR's mtdna_uncompress_v66.py helper).
"""
from __future__ import annotations

import gzip
import logging
from functools import lru_cache
from io import TextIOWrapper
from pathlib import Path
from typing import IO, Iterator

log = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parent.parent
MT_RCRS_PATH = ROOT / "data" / "mt_rcrs.txt"
MT_GENOME_LEN = 16_569


@lru_cache(maxsize=1)
def rcrs() -> str:
    """The revised Cambridge Reference Sequence (rCRS) — 16,569 bp."""
    seq = MT_RCRS_PATH.read_text().strip()
    if len(seq) != MT_GENOME_LEN:
        raise ValueError(
            f"rCRS file is {len(seq)} chars, expected {MT_GENOME_LEN}"
        )
    return seq


def _open_maybe_gz(path: Path) -> IO[str]:
    if path.suffix == ".gz":
        return TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8")
    return open(path, "r")


def _decompress(compressed: str, ref: str) -> str:
    """Substitute Q→ref, pad with '-' to MT_GENOME_LEN. Lowercase 'n' kept as-is."""
    # Pad first so zip yields exactly MT_GENOME_LEN pairs
    if len(compressed) < MT_GENOME_LEN:
        compressed = compressed + ("-" * (MT_GENOME_LEN - len(compressed)))
    elif len(compressed) > MT_GENOME_LEN:
        compressed = compressed[:MT_GENOME_LEN]
    return "".join(r if c == "Q" else c for c, r in zip(compressed, ref))


def iter_mt_records(path: Path) -> Iterator[tuple[str, str]]:
    """
    Yield (genetic_id, mt_sequence) tuples from the AADR mt repo file.
    Reads a `.gz` file transparently. Skips malformed lines with a warn.
    """
    ref = rcrs()
    with _open_maybe_gz(path) as fh:
        for line_no, raw in enumerate(fh, start=1):
            line = raw.rstrip("\n")
            if not line:
                continue
            if not line.startswith(">"):
                log.warning("%s line %d: expected '>' prefix, skipping", path, line_no)
                continue
            try:
                header, compressed = line[1:].split("\t", 1)
            except ValueError:
                log.warning("%s line %d: missing TAB separator, skipping", path, line_no)
                continue
            yield header.strip(), _decompress(compressed, ref)


def load_mt_repo(path: Path) -> dict[str, str]:
    """
    Load the entire mt repo into memory: {genetic_id: 16569-char sequence}.
    For 7k+ individuals at 16,569 chars each this is ~115 MB of Python
    strings — acceptable for a one-shot pipeline run.
    """
    out = dict(iter_mt_records(path))
    log.info("Loaded %d mt sequences from %s", len(out), path.name)
    return out


def is_called(base: str) -> bool:
    """A position is callable if it's a single real base (not missing/padding)."""
    return len(base) == 1 and base in "ACGT"


def count_pairwise_diffs(
    seq_a: str,
    seq_b: str,
    positions: list[int] | None = None,
) -> tuple[int, int]:
    """
    Count (k, L) where k = sites differing, L = sites called in both.
    `positions` is an optional 0-indexed list to restrict the comparison
    to (e.g. only the modern AncestryDNA-covered mt positions). When None,
    iterates over the full 16,569.
    """
    k = l = 0
    if positions is None:
        positions = range(MT_GENOME_LEN)
    for i in positions:
        a, b = seq_a[i], seq_b[i]
        if is_called(a) and is_called(b):
            l += 1
            if a != b:
                k += 1
    return k, l
