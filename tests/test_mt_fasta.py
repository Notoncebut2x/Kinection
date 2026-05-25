"""
Tests for scripts.utils.mt_fasta — AADR mt FASTA parser.

The compression scheme is non-standard (Q = "same as rCRS"), so wrong-by-one
errors in decompression would silently produce nonsense TMRCAs. These tests
pin down the format contract.
"""
from __future__ import annotations

import gzip
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

from utils.mt_fasta import (  # noqa: E402
    MT_GENOME_LEN,
    _decompress,
    count_pairwise_diffs,
    is_called,
    iter_mt_records,
    load_mt_repo,
    rcrs,
)


# ---------------------------------------------------------------------------
# rCRS reference
# ---------------------------------------------------------------------------

def test_rcrs_length():
    assert len(rcrs()) == MT_GENOME_LEN == 16_569


def test_rcrs_starts_with_known_sequence():
    # rCRS position 1-50 starts: GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCAT
    assert rcrs()[:50] == "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCAT"


def test_rcrs_only_acgtn():
    bad = set(rcrs()) - set("ACGTN")
    assert not bad, f"rCRS contains unexpected chars: {bad}"


# ---------------------------------------------------------------------------
# _decompress
# ---------------------------------------------------------------------------

def test_decompress_q_substitutes_rcrs():
    ref = "ACGTACGT"
    compressed = "QQQQQQQQ"
    assert _decompress(compressed, ref) == ref


def test_decompress_keeps_explicit_bases():
    ref = "ACGTACGT"
    compressed = "TGCAQQQQ"
    # Positions 0-3 are explicit, 4-7 are Q → use ref
    assert _decompress(compressed, ref) == "TGCAACGT"


def test_decompress_pads_short_sequences():
    ref = "ACGTACGT"
    compressed = "TG"
    # Padded with '-' to length of ref
    assert _decompress(compressed, ref) == "TG------"


def test_decompress_truncates_long_sequences():
    ref = "ACGT"
    compressed = "TGCATGCA"
    assert _decompress(compressed, ref) == "TGCA"


def test_decompress_preserves_missing_n():
    ref = "ACGTACGT"
    compressed = "nnnnQQQQ"
    assert _decompress(compressed, ref) == "nnnnACGT"


# ---------------------------------------------------------------------------
# iter_mt_records and load_mt_repo
# ---------------------------------------------------------------------------

def _write_gz(path: Path, content: str) -> None:
    with gzip.open(path, "wb") as fh:
        fh.write(content.encode("utf-8"))


def test_iter_mt_records_parses_id_and_seq(tmp_path):
    path = tmp_path / "mini.fa.gz"
    _write_gz(path, ">SAMPLE1\tQQQQ\n>SAMPLE2\tAAAA\n")
    records = list(iter_mt_records(path))
    assert len(records) == 2
    assert records[0][0] == "SAMPLE1"
    assert records[1][0] == "SAMPLE2"
    # Both are length MT_GENOME_LEN (rCRS-padded for SAMPLE1, base-padded for SAMPLE2)
    assert len(records[0][1]) == MT_GENOME_LEN
    assert len(records[1][1]) == MT_GENOME_LEN


def test_iter_mt_records_skips_malformed_lines(tmp_path):
    path = tmp_path / "broken.fa.gz"
    _write_gz(path, ">SAMPLE1\tQQQQ\nNOT A HEADER\n>SAMPLE2 missing tab\n>SAMPLE3\tAAAA\n")
    records = list(iter_mt_records(path))
    ids = [r[0] for r in records]
    assert ids == ["SAMPLE1", "SAMPLE3"]


def test_load_mt_repo_returns_dict(tmp_path):
    path = tmp_path / "two.fa.gz"
    _write_gz(path, ">A\tQQQQ\n>B\tTTTT\n")
    repo = load_mt_repo(path)
    assert set(repo.keys()) == {"A", "B"}
    assert all(len(s) == MT_GENOME_LEN for s in repo.values())


# ---------------------------------------------------------------------------
# Difference counting
# ---------------------------------------------------------------------------

def test_count_pairwise_diffs_zero_when_identical():
    seq = "A" * MT_GENOME_LEN
    k, l = count_pairwise_diffs(seq, seq)
    assert k == 0
    assert l == MT_GENOME_LEN


def test_count_pairwise_diffs_skips_missing():
    a = "ACGTnn--"
    b = "ACGT----"
    k, l = count_pairwise_diffs(a, b, positions=range(8))
    # Only positions 0-3 are called in both; all match.
    assert k == 0
    assert l == 4


def test_count_pairwise_diffs_counts_mismatches():
    a = "ACGTACGT"
    b = "TCGAACGT"
    k, l = count_pairwise_diffs(a, b, positions=range(8))
    assert k == 2  # pos 0 A vs T, pos 3 T vs A
    assert l == 8


def test_count_pairwise_diffs_position_subset():
    a = "ACGTACGT"
    b = "TGCATGCA"  # all 8 positions differ
    k, l = count_pairwise_diffs(a, b, positions=[0, 2, 4])
    assert k == 3
    assert l == 3


# ---------------------------------------------------------------------------
# is_called
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("base,expected", [
    ("A", True), ("C", True), ("G", True), ("T", True),
    ("N", False), ("n", False), ("-", False), ("Q", False), ("", False),
])
def test_is_called(base, expected):
    assert is_called(base) is expected
