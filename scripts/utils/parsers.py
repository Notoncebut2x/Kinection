"""
Parsers for all data formats used in the DNA lineage analysis pipeline.

Handles:
  - AncestryDNA raw export (V1.0 and V2.0 arrays)
  - EIGENSTRAT PACKGENO binary (.geno)
  - EIGENSTRAT .ind individual list
  - Allen aDNA Resource .anno annotation file
"""

from __future__ import annotations

import csv
import struct
import logging
import mmap
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator

import numpy as np

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# CHROMOSOME MAP: AncestryDNA numeric codes → standard labels
# ---------------------------------------------------------------------------
ANCESTRY_CHROM_MAP = {
    str(i): str(i) for i in range(1, 23)
}
ANCESTRY_CHROM_MAP.update({
    "23": "X",
    "24": "Y",
    "25": "PAR",   # pseudo-autosomal region
    "26": "MT",
    "chromosome": None,  # header row artifact — skip
})

COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def complement(allele: str) -> str:
    return allele.translate(COMPLEMENT)


# ---------------------------------------------------------------------------
# Data containers
# ---------------------------------------------------------------------------

@dataclass
class SNP:
    rsid: str
    chrom: str          # standardised (1-22, X, Y, MT, PAR)
    position: int
    allele1: str        # forward-strand allele 1
    allele2: str        # forward-strand allele 2

    @property
    def is_homozygous(self) -> bool:
        return self.allele1 == self.allele2

    @property
    def genotype_str(self) -> str:
        return f"{self.allele1}{self.allele2}"


@dataclass
class Individual:
    genetic_id: str
    sex: str            # M / F / U
    population: str     # population label from .ind file
    index: int          # 0-based column index in .geno file


@dataclass
class AnnoRecord:
    genetic_id: str
    group_id: str
    locality: str
    political_entity: str
    lat: float | None
    lon: float | None
    date_bp: float | None       # years before 1950
    date_std: float | None
    date_range: str
    molecular_sex: str
    ydna_terminal: str          # terminal Y haplogroup (Yfull notation)
    ydna_isogg: str             # ISOGG v15.73 Y haplogroup
    ydna_manual: str            # manually curated if different
    mtdna_haplogroup: str
    snps_1240k: int | None      # SNPs hit on 1240k autosomal set
    assessment: str

    @property
    def date_ce(self) -> float | None:
        """Return approximate year in CE (negative = BCE)."""
        if self.date_bp is None:
            return None
        return 1950 - self.date_bp

    @property
    def best_y_haplogroup(self) -> str:
        """Return the best available Y haplogroup call."""
        for field in (self.ydna_manual, self.ydna_terminal, self.ydna_isogg):
            if field and field not in ("..", "", "n/a", "n/a  (sex unknown)"):
                return field
        return ""

    @property
    def valid_mtdna(self) -> str:
        """Return mtDNA haplogroup if available."""
        if self.mtdna_haplogroup and self.mtdna_haplogroup not in ("..", "", "n/a"):
            return self.mtdna_haplogroup
        return ""


# ---------------------------------------------------------------------------
# AncestryDNA parser
# ---------------------------------------------------------------------------

def parse_ancestry_dna(path: Path) -> dict[str, SNP]:
    """
    Parse an AncestryDNA raw export file (V1.0 or V2.0).

    Returns a dict keyed by rsID → SNP.
    Chromosome codes are normalised (24→Y, 26→MT, etc.).
    SNPs with unknown chromosome labels or indels ('I'/'D') are skipped.
    """
    path = Path(path)
    snps: dict[str, SNP] = {}
    skipped = 0

    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 5:
                continue
            rsid, raw_chrom, raw_pos, a1, a2 = (
                parts[0], parts[1], parts[2], parts[3], parts[4]
            )
            if rsid == "rsid":
                continue  # header row

            chrom = ANCESTRY_CHROM_MAP.get(raw_chrom)
            if chrom is None:
                skipped += 1
                continue

            # Skip indels — not useful for SNP matching
            if "I" in (a1, a2) or "D" in (a1, a2):
                skipped += 1
                continue
            # Skip 0-allele missing calls
            if a1 == "0" or a2 == "0":
                skipped += 1
                continue

            try:
                position = int(raw_pos)
            except ValueError:
                skipped += 1
                continue

            snps[rsid] = SNP(
                rsid=rsid,
                chrom=chrom,
                position=position,
                allele1=a1.upper(),
                allele2=a2.upper(),
            )

    logger.info(
        "Parsed %d SNPs from %s (skipped %d)", len(snps), path.name, skipped
    )
    return snps


# ---------------------------------------------------------------------------
# EIGENSTRAT .ind parser
# ---------------------------------------------------------------------------

def parse_ind_file(path: Path) -> list[Individual]:
    """
    Parse an EIGENSTRAT .ind file.

    Format: whitespace-delimited, columns are: ID  sex  population
    Returns a list ordered by their 0-based index (= column in .geno).
    """
    individuals: list[Individual] = []
    path = Path(path)

    with open(path, encoding="utf-8") as fh:
        for idx, line in enumerate(fh):
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 3:
                logger.warning("Malformed .ind line %d: %r", idx, line)
                continue
            individuals.append(Individual(
                genetic_id=parts[0],
                sex=parts[1],
                population=" ".join(parts[2:]),
                index=idx,
            ))

    logger.info("Parsed %d individuals from %s", len(individuals), path.name)
    return individuals


# ---------------------------------------------------------------------------
# Allen aDNA Resource .anno parser
# ---------------------------------------------------------------------------

# Column names are very long — map to short keys
_ANNO_COLUMN_MAP = {
    "Genetic ID (suffixes: \".DG\" is a high coverage shotgun genome with diploid genotype calls, \".AG\" is shotgun data with each position in the genome represented by a randomly chosen sequence, \".HO\" is Affymetrix Human Origins genotype data)": "genetic_id",
    "Group ID": "group_id",
    "Locality": "locality",
    "Political Entity": "political_entity",
    "Lat.": "lat",
    "Long.": "lon",
    "Date mean in BP in years before 1950 CE [OxCal mu for a direct radiocarbon date, and average of range for a contextual date]": "date_bp",
    "Date standard deviation in BP [OxCal sigma for a direct radiocarbon date, and standard deviation of the uniform distribution between the two bounds for a contextual date]": "date_std",
    "Full Date One of two formats. (Format 1) 95.4% CI calibrated radiocarbon age (Conventional Radiocarbon Age BP, Lab number) e.g. 2624-2350 calBCE (3990±40 BP, Ua-35016). (Format 2) Archaeological context range, e.g. 2500-1700 BCE": "date_range",
    "Molecular Sex": "molecular_sex",
    "Y haplogroup in terminal mutation notation automatically called based on Yfull with the software described in Lazaridis et al. Science 2022": "ydna_terminal",
    "Y haplogroup  in ISOGG v15.73 notation automatically called based on Yfull with the software described in Lazaridis et al. Science 2022": "ydna_isogg",
    "Y haplogroup manually called if different from automatic": "ydna_manual",
    "mtDNA haplogroup if >2x or published": "mtdna_haplogroup",
    "SNPs hit on autosomal targets (Computed using easystats on 1240k snpset)": "snps_1240k",
    "ASSESSMENT": "assessment",
}


def _safe_float(val: str) -> float | None:
    try:
        return float(val)
    except (ValueError, TypeError):
        return None


def _safe_int(val: str) -> int | None:
    try:
        return int(val)
    except (ValueError, TypeError):
        return None


def parse_anno_file(path: Path) -> dict[str, AnnoRecord]:
    """
    Parse the Allen aDNA Resource .anno annotation file.

    Returns a dict keyed by genetic_id → AnnoRecord.
    """
    path = Path(path)
    records: dict[str, AnnoRecord] = {}

    with open(path, encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader(fh, delimiter="\t")

        # Build a column → short-key mapping for whatever columns exist
        col_map: dict[str, str] = {}
        for full_name in (reader.fieldnames or []):
            stripped = full_name.strip()
            if stripped in _ANNO_COLUMN_MAP:
                col_map[full_name] = _ANNO_COLUMN_MAP[stripped]

        for row in reader:
            def get(key: str, default: str = "") -> str:
                for full, short in col_map.items():
                    if short == key:
                        return (row.get(full) or "").strip()
                return default

            genetic_id = get("genetic_id")
            if not genetic_id:
                continue

            rec = AnnoRecord(
                genetic_id=genetic_id,
                group_id=get("group_id"),
                locality=get("locality"),
                political_entity=get("political_entity"),
                lat=_safe_float(get("lat")),
                lon=_safe_float(get("lon")),
                date_bp=_safe_float(get("date_bp")),
                date_std=_safe_float(get("date_std")),
                date_range=get("date_range"),
                molecular_sex=get("molecular_sex"),
                ydna_terminal=get("ydna_terminal"),
                ydna_isogg=get("ydna_isogg"),
                ydna_manual=get("ydna_manual"),
                mtdna_haplogroup=get("mtdna_haplogroup"),
                snps_1240k=_safe_int(get("snps_1240k")),
                assessment=get("assessment"),
            )
            records[genetic_id] = rec

    logger.info("Parsed %d anno records from %s", len(records), path.name)
    return records


# ---------------------------------------------------------------------------
# EIGENSTRAT PACKGENO binary .geno parser
# ---------------------------------------------------------------------------

GENO_HEADER_SIZE = 64  # bytes — fixed in EIGENSOFT source


@dataclass
class GenoFile:
    """
    Memory-mapped interface to an EIGENSTRAT PACKGENO binary .geno file.

    Format:
      - 64-byte text header:  "GENO   {n_indiv} {n_snps} {hash1} {hash2}\\0..."
      - Then n_snps rows, each ceil(n_indiv / 4) bytes
      - Each byte packs 4 genotypes as 2-bit values (LSB first):
          0 = homozygous ref (AA)
          1 = heterozygous (Aa)
          2 = homozygous alt (aa)
          3 = missing
    """
    path: Path
    n_indiv: int
    n_snps: int
    bytes_per_snp: int
    _mm: mmap.mmap = field(default=None, repr=False)
    _fh: object = field(default=None, repr=False)

    @classmethod
    def open(cls, path: Path) -> "GenoFile":
        path = Path(path)
        with open(path, "rb") as fh:
            raw_header = fh.read(GENO_HEADER_SIZE)

        header_str = raw_header.split(b"\x00")[0].decode("ascii", errors="replace")
        parts = header_str.split()
        if parts[0] != "GENO":
            raise ValueError(f"Not a PACKGENO file (magic={parts[0]!r}): {path}")

        n_indiv = int(parts[1])
        n_snps = int(parts[2])
        bytes_per_snp = (n_indiv + 3) // 4

        obj = cls(
            path=path,
            n_indiv=n_indiv,
            n_snps=n_snps,
            bytes_per_snp=bytes_per_snp,
        )

        # Open memory-mapped file for efficient random access
        fh = open(path, "rb")
        obj._fh = fh
        obj._mm = mmap.mmap(fh.fileno(), 0, access=mmap.ACCESS_READ)

        expected_size = GENO_HEADER_SIZE + n_snps * bytes_per_snp
        actual_size = obj._mm.size()
        if actual_size < expected_size:
            missing = expected_size - actual_size
            missing_snps = missing // bytes_per_snp
            logger.warning(
                "GENO file is %d bytes short (~%d SNPs truncated). "
                "The download may be incomplete.",
                missing, missing_snps,
            )

        logger.info(
            "Opened GENO file: %d individuals, %d SNPs, %d bytes/SNP",
            n_indiv, n_snps, bytes_per_snp,
        )
        return obj

    def close(self) -> None:
        if self._mm:
            self._mm.close()
        if self._fh:
            self._fh.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    def read_snp_row(self, snp_index: int) -> np.ndarray:
        """
        Return genotypes for all individuals at a given SNP (0-based index).

        Output: int8 array of length n_indiv
            0 = hom ref, 1 = het, 2 = hom alt, 3 = missing
        """
        # Cast to int64 before multiply — snp_index * bytes_per_snp can exceed int32 max
        offset = GENO_HEADER_SIZE + int(snp_index) * self.bytes_per_snp
        raw = self._mm[offset: offset + self.bytes_per_snp]
        raw_np = np.frombuffer(raw, dtype=np.uint8)

        # Unpack 4 genotypes per byte (2 bits each, LSB first)
        gt = np.empty(len(raw_np) * 4, dtype=np.int8)
        gt[0::4] = raw_np & 0x03
        gt[1::4] = (raw_np >> 2) & 0x03
        gt[2::4] = (raw_np >> 4) & 0x03
        gt[3::4] = (raw_np >> 6) & 0x03

        return gt[: self.n_indiv]

    def read_individual_column(self, indiv_index: int) -> np.ndarray:
        """
        Return genotypes for one individual across all SNPs.

        This reads the entire file sequentially — use read_snp_row for
        random access to individual SNPs.

        Output: int8 array of length n_snps
        """
        byte_col = indiv_index // 4
        bit_shift = (indiv_index % 4) * 2
        mask = 0x03

        genotypes = np.empty(self.n_snps, dtype=np.int8)
        for snp_i in range(self.n_snps):
            offset = GENO_HEADER_SIZE + int(snp_i) * self.bytes_per_snp + byte_col
            if offset >= self._mm.size():
                genotypes[snp_i] = 3  # missing — truncated file
                continue
            byte_val = self._mm[offset]
            genotypes[snp_i] = (byte_val >> bit_shift) & mask

        return genotypes

    def iter_snp_rows(self, snp_indices: list[int] | None = None) -> Iterator[tuple[int, np.ndarray]]:
        """
        Iterate over SNP rows, yielding (snp_index, genotype_array) tuples.

        If snp_indices is None, iterates all SNPs in order.
        """
        indices = snp_indices if snp_indices is not None else range(self.n_snps)
        for i in indices:
            yield i, self.read_snp_row(i)
