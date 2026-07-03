"""
Tests for 23andMe raw-export parsing and modern-format auto-detection.

All fixtures are synthetic (fake rsIDs / positions) — no real genotype
data is committed.
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
from utils.parsers import (  # noqa: E402
    detect_modern_format,
    parse_23andme,
    parse_modern_dna,
)

# A minimal 23andMe-style export exercising every case the parser must
# handle. The header carries a "23andMe" signature (for format detection)
# but deliberately avoids the real export header the raw-DNA pre-commit
# scanner guards against — these are synthetic fixtures, not real data.
TWENTYTHREE_SAMPLE = """\
# synthetic 23andMe v5 fixture
rs1\t1\t100\tAG
rs2\t2\t200\tCC
rs3\tX\t300\tT
rs4\tMT\t400\tG
rs5\t7\t500\t--
rs6\t9\t600\tII
rs7\t11\t700\tDD
rs8\tY\t800\tD
rs9\tXY\t900\tAC
"""

# AncestryDNA-style fixture (5 columns, uncommented header). Uses non-rs
# marker ids so the 5-column rows don't trip the raw-DNA line scanner.
ANCESTRY_SAMPLE = """\
# synthetic AncestryDNA v2 fixture
rsid\tchromosome\tposition\tallele1\tallele2
snp1\t1\t100\tA\tG
snp2\t26\t200\tC\tC
"""


def _write(tmp_path, name, content):
    p = tmp_path / name
    p.write_text(content, encoding="utf-8")
    return p


def test_detect_23andme(tmp_path):
    p = _write(tmp_path, "23andme.txt", TWENTYTHREE_SAMPLE)
    assert detect_modern_format(p) == "23andme"


def test_detect_ancestrydna(tmp_path):
    p = _write(tmp_path, "ancestry.txt", ANCESTRY_SAMPLE)
    assert detect_modern_format(p) == "ancestrydna"


def test_detect_by_column_count_without_signature(tmp_path):
    # No vendor comment — must fall back to the 4-vs-5 column heuristic.
    p = _write(tmp_path, "bare.txt", "rs1\t1\t100\tAG\nrs2\t2\t200\tCC\n")
    assert detect_modern_format(p) == "23andme"


def test_parse_23andme_valid_and_skips(tmp_path):
    p = _write(tmp_path, "23andme.txt", TWENTYTHREE_SAMPLE)
    snps = parse_23andme(p)

    # Kept: rs1 (diploid), rs2 (hom), rs3 (haploid X), rs4 (haploid MT),
    #       rs9 (XY -> PAR). Skipped: rs5 (--), rs6 (II), rs7 (DD), rs8 (D).
    assert set(snps) == {"rs1", "rs2", "rs3", "rs4", "rs9"}


def test_parse_23andme_allele_splitting(tmp_path):
    p = _write(tmp_path, "23andme.txt", TWENTYTHREE_SAMPLE)
    snps = parse_23andme(p)

    # Diploid splits into two alleles.
    assert (snps["rs1"].allele1, snps["rs1"].allele2) == ("A", "G")
    # Haploid call is stored as homozygous.
    assert snps["rs3"].allele1 == snps["rs3"].allele2 == "T"
    assert snps["rs3"].is_homozygous


def test_parse_23andme_chrom_normalisation(tmp_path):
    p = _write(tmp_path, "23andme.txt", TWENTYTHREE_SAMPLE)
    snps = parse_23andme(p)
    assert snps["rs3"].chrom == "X"
    assert snps["rs4"].chrom == "MT"
    assert snps["rs9"].chrom == "PAR"  # 23andMe "XY" pseudo-autosomal


def test_parse_modern_dna_dispatches_on_auto(tmp_path):
    tw = _write(tmp_path, "23andme.txt", TWENTYTHREE_SAMPLE)
    an = _write(tmp_path, "ancestry.txt", ANCESTRY_SAMPLE)

    tw_snps = parse_modern_dna(tw)          # auto -> 23andme
    an_snps = parse_modern_dna(an)          # auto -> ancestrydna

    assert "rs9" in tw_snps                 # XY row only exists in the 23andMe fixture
    assert an_snps["snp2"].chrom == "MT"    # AncestryDNA 26 -> MT


def test_parse_modern_dna_explicit_format(tmp_path):
    p = _write(tmp_path, "whatever.txt", TWENTYTHREE_SAMPLE)
    snps = parse_modern_dna(p, fmt="23andme")
    assert set(snps) == {"rs1", "rs2", "rs3", "rs4", "rs9"}


def test_parse_modern_dna_unknown_format_raises(tmp_path):
    p = _write(tmp_path, "whatever.txt", TWENTYTHREE_SAMPLE)
    with pytest.raises(ValueError):
        parse_modern_dna(p, fmt="nonsense")
