#!/usr/bin/env python3
"""
Step 1.3 — Genome-wide SNP Similarity and PCA Projection

For each of the 17,629 ancient individuals, computes allele-sharing distance
to Individual 1 at the 394,671 overlapping SNPs. Aggregates to population
level, ranks populations by closeness, and projects Individual 1 into the
principal component space of the ancient dataset.

Method:
  - Allele-sharing distance (ASD): mean |modern_freq - ancient_freq| over
    non-missing overlapping SNPs. Modern diploid freq = dosage/2 (0, 0.5, 1).
    Ancient pseudo-haploid freq = 0.0 (hom ref) or 1.0 (hom alt).
  - Population distance: mean ASD across all PASS individuals in a population
    with ≥ MIN_INDIV_SNPS overlap SNPs.
  - PCA: computed on high-coverage ancient individuals (≥ MIN_PCA_SNPS),
    Individual 1 projected in using pre-computed eigenvectors.

Outputs (output/step3/):
  pairwise_distances.tsv       — ASD to all 17,629 ancient individuals
  population_distances.tsv     — ranked population mean distances
  top_matches_report.md        — narrative of top ancient population matches
  pca_coordinates.tsv          — PCA coordinates for ancient + Individual 1

Usage:
    python scripts/step3_similarity_pca.py
"""

from __future__ import annotations

import json
import logging
import os
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# R2 / local mode switch
# ---------------------------------------------------------------------------
USE_R2 = os.environ.get('USE_R2', '').lower() in ('1', 'true', 'yes')
JOB_ID = os.environ.get('JOB_ID', 'dev')
# When set, do not upload outputs to R2 and read snp_overlap.tsv from local disk.
LOCAL_OUTPUTS = os.environ.get('LOCAL_OUTPUTS', '').lower() in ('1', 'true', 'yes')
# Suffix used for output and handoff paths; must match the value used in step 1.
OUTPUT_LABEL = os.environ.get('OUTPUT_LABEL', 'rn')

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
# Min overlap SNPs to include an individual in the top-matches list / population
# stats. Raised from 10k to 50k: at 10k, per-individual ASD has so much sampling
# noise that low-coverage outliers (e.g. 14k-SNP samples) routinely topped the
# rankings spuriously. 50k matches the threshold used by the synthesis step
# (step 1.6 MIN_INDIV_SNPS_FOR_MAP) for consistency.
MIN_INDIV_SNPS   = 50_000
MIN_PCA_SNPS     = 100_000  # min 1240k SNPs (from .anno) to include individual in PCA
PCA_N_COMPONENTS = 10
# Subsample SNPs used for PCA — peak memory is roughly
#   n_individuals * PCA_SNP_SUBSAMPLE * 4 bytes (float32)
# 30k SNPs × ~14k individuals ≈ 1.6 GB. Override via env if you have RAM to spare.
PCA_SNP_SUBSAMPLE = int(os.environ.get('PCA_SNP_SUBSAMPLE', '30000'))
CHUNK_SIZE       = 5_000     # SNPs per processing chunk
N_PSEUDO_DRAWS   = 10        # pseudo-haploid draws to average for ASD stability
AUTOSOME_CHROMS  = {str(i) for i in range(1, 23)}  # chr 1-22 only for ASD

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
ROOT   = Path(__file__).resolve().parent.parent
DATA   = ROOT / "data" / "input_data"
OUT1   = ROOT / "output" / f"step1_{OUTPUT_LABEL}"
OUTPUT = ROOT / "output" / f"step3_{OUTPUT_LABEL}"
OUTPUT.mkdir(parents=True, exist_ok=True)

# Local AADR resolved lazily in main() — works for any version (v62, v66, ...).
GENO_FILE: Path | None = None
IND_FILE:  Path | None = None
ANNO_FILE: Path | None = None
OVERLAP_TSV = OUT1 / "snp_overlap.tsv"
MODERN_NPY  = OUT1 / f"modern_indv_{OUTPUT_LABEL}_encoded.npy"

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-7s  %(message)s",
    datefmt="%H:%M:%S",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(OUTPUT / "step3.log", mode="w"),
    ],
)
log = logging.getLogger(__name__)

sys.path.insert(0, str(ROOT / "scripts"))
from utils.parsers import parse_ind_file, parse_anno_file, GenoFile
if USE_R2:
    from utils import r2_client
    from utils.r2_geno import R2GenoFile


# ---------------------------------------------------------------------------
# Load overlap SNP table
# ---------------------------------------------------------------------------

def load_overlap(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Load the SNP overlap table from Step 1.1.

    Returns:
        geno_indices:  int32 array  — which GENO rows to read (sorted)
        modern_dosages: int8 array  — alt-allele dosage (0/1/2, -1=missing)
        chroms:        object array — chromosome label per SNP (e.g. "1", "X")
    """
    geno_indices = []
    modern_dosages = []
    chroms = []

    with open(path) as fh:
        header = fh.readline()  # skip header
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            geno_idx = int(parts[0])
            chrom = parts[2]
            dosage_str = parts[8]
            dosage = int(dosage_str) if dosage_str not in ("NA", "") else -1
            geno_indices.append(geno_idx)
            modern_dosages.append(dosage)
            chroms.append(chrom)

    geno_indices   = np.array(geno_indices,  dtype=np.int32)
    modern_dosages = np.array(modern_dosages, dtype=np.int8)
    chroms         = np.array(chroms,         dtype=object)

    # Confirm sorted (step1 sorts by geno_index)
    assert np.all(np.diff(geno_indices) >= 0), "geno_indices must be sorted"

    n_valid = np.sum(modern_dosages >= 0)
    log.info(
        "Loaded %d overlap SNPs (%d with valid modern dosage, %d missing)",
        len(geno_indices), n_valid, len(geno_indices) - n_valid,
    )
    return geno_indices, modern_dosages, chroms


# ---------------------------------------------------------------------------
# Allele-sharing distance
# ---------------------------------------------------------------------------

def compute_asd(
    geno: GenoFile,
    geno_indices: np.ndarray,
    modern_dosages: np.ndarray,
    chroms: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute allele-sharing distance between the modern individual and each
    ancient individual using pseudo-haploidization.

    Since ancient samples are pseudo-haploid (0 or 1, never 0.5), we
    pseudo-haploidize the modern diploid individual: at each heterozygous
    site, randomly pick one allele (0 or 1). This makes both sides of the
    comparison pseudo-haploid, removing the systematic 0.5 bias at het sites
    that inflates distances to ~0.49 regardless of true ancestry match.

    We run N_PSEUDO_DRAWS independent random draws and average the resulting
    ASD scores for stability (~10 draws converges well at 300k+ autosomal SNPs).

    Only autosomal SNPs (chr 1–22) are used.

    GENO encoding:
      0 = hom ref → alt allele = 0
      1 = het (ancient — rare) → 0.5 (averaged out by pseudo-haploid draws)
      2 = hom alt → alt allele = 1
      3 = missing

    Returns:
        asd_sum:     float64 array (n_indiv,) — sum of ASD across draws
        count_valid: int64 array   (n_indiv,) — non-missing comparisons (consistent across draws)
    """
    rng = np.random.default_rng(42)
    n_indiv = geno.n_indiv

    # Filter to autosomes only
    autosome_mask = np.array([c in AUTOSOME_CHROMS for c in chroms])
    geno_indices_auto  = geno_indices[autosome_mask]
    modern_dosages_auto = modern_dosages[autosome_mask]
    n_snps = len(geno_indices_auto)
    n_het  = int((modern_dosages_auto == 1).sum())
    log.info(
        "Autosome filter: %d / %d SNPs retained (%d het sites, %.1f%%)",
        n_snps, len(geno_indices), n_het, 100 * n_het / n_snps if n_snps else 0,
    )

    # Missing mask (consistent across draws)
    missing_mask = (modern_dosages_auto < 0)

    # AADR PACKGENO: 0 = hom allele2 (col 6 of .snp), 2 = hom allele1 (col 5).
    # Step 1's `dosage` counts allele2 in modern, so modern_freq=dosage/2 is
    # "allele2 frequency". To compare like-with-like we map ancient geno to
    # ancient allele2-frequency: geno 0 -> 1.0, geno 2 -> 0.0, het -> 0.5,
    # missing -> nan. Without this inversion, identical samples register as
    # maximum distance — see ADR 0014.
    geno_to_freq = np.array([1.0, 0.5, 0.0, np.nan], dtype=np.float32)

    n_chunks = (n_snps + CHUNK_SIZE - 1) // CHUNK_SIZE

    # We accumulate ASD across draws and average at the end
    asd_sum     = np.zeros(n_indiv, dtype=np.float64)
    count_valid = np.zeros(n_indiv, dtype=np.int64)  # same every draw

    # Pre-read all GENO rows once (avoid re-reading for each draw)
    # Store as uint8 matrix — shape (n_snps, n_indiv)
    # For large datasets this may be ~394k × 17k × 1 byte = ~6.7 GB — too large.
    # Instead we iterate chunks, and for each chunk iterate draws.
    # This means n_draws × n_chunks reads, but chunk re-reading is fast from mmap.

    t_start = time.time()

    for chunk_i in range(n_chunks):
        chunk_start = chunk_i * CHUNK_SIZE
        chunk_end   = min(chunk_start + CHUNK_SIZE, n_snps)
        chunk_size  = chunk_end - chunk_start

        chunk_geno_idx    = geno_indices_auto[chunk_start:chunk_end]
        chunk_dosages     = modern_dosages_auto[chunk_start:chunk_end]
        chunk_missing     = missing_mask[chunk_start:chunk_end]

        # Read all rows in this chunk with one range request (R2) or mmap reads (local)
        rows = geno.read_chunk(chunk_geno_idx).view(np.uint8)

        # Ancient frequencies (chunk_size × n_indiv), 0.0 / 1.0 / nan
        ancient_freq = geno_to_freq[rows]
        ancient_missing = np.isnan(ancient_freq)  # (chunk_size, n_indiv)

        # Compute count_valid once (same for all draws — modern missing doesn't change)
        modern_missing_2d = chunk_missing[:, np.newaxis]  # (chunk_size, 1)
        valid = ~modern_missing_2d & ~ancient_missing      # (chunk_size, n_indiv)
        if chunk_i == 0:
            count_valid += valid.sum(axis=0, dtype=np.int64)
        else:
            count_valid += valid.sum(axis=0, dtype=np.int64)

        # Het sites within this chunk
        het_in_chunk = (chunk_dosages == 1)

        # Run N_PSEUDO_DRAWS for this chunk
        chunk_sum_diff = np.zeros(n_indiv, dtype=np.float64)
        for _ in range(N_PSEUDO_DRAWS):
            # Pseudo-haploidize: hom_ref→0.0, hom_alt→1.0, het→random 0 or 1, missing→nan
            random_alleles = rng.integers(0, 2, size=chunk_size, dtype=np.uint8).astype(np.float32)
            pseudo_freq = np.where(
                chunk_missing, np.nan,
                np.where(het_in_chunk, random_alleles,
                         np.where(chunk_dosages == 2, 1.0, 0.0))
            ).astype(np.float32)

            modern_2d = pseudo_freq[:, np.newaxis]  # (chunk_size, 1)
            diff = np.abs(modern_2d - ancient_freq)
            diff[~valid] = 0.0
            chunk_sum_diff += diff.sum(axis=0)

        # Average across draws
        asd_sum += chunk_sum_diff / N_PSEUDO_DRAWS

        if (chunk_i + 1) % 10 == 0 or chunk_i == n_chunks - 1:
            elapsed = time.time() - t_start
            pct = (chunk_i + 1) / n_chunks * 100
            log.info(
                "  Chunk %d/%d (%d SNPs, %.0f%%)  %.1fs elapsed",
                chunk_i + 1, n_chunks, chunk_end, pct, elapsed,
            )

    return asd_sum, count_valid


# ---------------------------------------------------------------------------
# PCA
# ---------------------------------------------------------------------------

def compute_pca(
    geno: GenoFile,
    geno_indices: np.ndarray,
    modern_dosages: np.ndarray,
    pca_indiv_mask: np.ndarray,
    n_components: int = PCA_N_COMPONENTS,
) -> dict:
    """
    Compute PCA on a subset of high-coverage ancient individuals, then
    project Individual 1 into the resulting space.

    Returns dict with keys: eigenvectors, eigenvalues, ancient_coords,
    modern_coords, indiv_indices.
    """
    try:
        from sklearn.decomposition import PCA as SklearnPCA
        use_sklearn = True
    except ImportError:
        use_sklearn = False
        log.warning("scikit-learn not available — using numpy SVD for PCA")

    n_pca_indiv = pca_indiv_mask.sum()
    pca_indiv_indices = np.where(pca_indiv_mask)[0]
    log.info("PCA: %d individuals, subsampling to %d SNPs", n_pca_indiv, PCA_SNP_SUBSAMPLE)

    # Subsample SNPs evenly
    step = max(1, len(geno_indices) // PCA_SNP_SUBSAMPLE)
    sub_idx = np.arange(0, len(geno_indices), step)[:PCA_SNP_SUBSAMPLE]
    sub_geno_idx = geno_indices[sub_idx]
    sub_modern   = modern_dosages[sub_idx]

    n_sub = len(sub_idx)
    log.info("PCA: using %d SNPs (every %d-th overlap SNP)", n_sub, step)

    # Build matrix: (n_pca_indiv, n_sub) of allele2 dosage (0/1/2, -1=missing).
    # AADR PACKGENO: geno 0 = hom allele2 -> dosage 2; geno 2 = hom allele1
    # -> dosage 0; geno 1 = het -> dosage 1; geno 3 = missing -> -1.
    # This matches step 1's modern dosage convention (count of allele2 in modern)
    # so modern projection is on the same axes as the ancient PCA — see ADR 0014.
    geno_to_dosage = np.array([2, 1, 0, -1], dtype=np.int8)

    log.info("Building ancient genotype matrix (chunked reads)...")
    matrix = np.empty((n_pca_indiv, n_sub), dtype=np.float32)
    matrix[:] = np.nan

    # Read in chunks of CHUNK_SIZE SNPs — one range request per chunk in R2 mode
    for chunk_start in range(0, n_sub, CHUNK_SIZE):
        chunk_end_j = min(chunk_start + CHUNK_SIZE, n_sub)
        chunk_sub_idx = sub_geno_idx[chunk_start:chunk_end_j]

        chunk_rows = geno.read_chunk(chunk_sub_idx)                        # (c, n_indiv) int8
        chunk_dosage = geno_to_dosage[chunk_rows.view(np.uint8)]           # (c, n_indiv) int8
        chunk_pca = chunk_dosage[:, pca_indiv_indices].astype(np.float32)  # (c, n_pca)
        chunk_pca[chunk_pca < 0] = np.nan
        matrix[:, chunk_start:chunk_end_j] = chunk_pca.T

    # Impute missing with column mean (simple imputation)
    log.info("Imputing missing values with column means...")
    col_means = np.nanmean(matrix, axis=0)
    nan_mask = np.isnan(matrix)
    matrix[nan_mask] = np.take(col_means, np.where(nan_mask)[1])

    # Mean-center columns
    matrix -= col_means[np.newaxis, :]

    log.info("Running PCA (n_components=%d)...", n_components)
    if use_sklearn:
        pca = SklearnPCA(n_components=n_components, random_state=42)
        ancient_coords = pca.fit_transform(matrix)           # (n_pca_indiv, n_components)
        eigenvectors   = pca.components_                     # (n_components, n_sub)
        eigenvalues    = pca.explained_variance_ratio_
    else:
        U, s, Vt = np.linalg.svd(matrix, full_matrices=False)
        ancient_coords = U[:, :n_components] * s[:n_components]
        eigenvectors   = Vt[:n_components]
        total_var = np.sum(s ** 2)
        eigenvalues = (s[:n_components] ** 2) / total_var

    # Project Individual 1
    modern_vec = np.where(sub_modern >= 0, sub_modern.astype(np.float32), col_means)
    modern_vec -= col_means
    modern_coords = (eigenvectors @ modern_vec)              # (n_components,)

    log.info("PCA complete. Variance explained: %s",
             ", ".join(f"PC{i+1}={v:.1%}" for i, v in enumerate(eigenvalues)))

    return {
        "eigenvectors":       eigenvectors,
        "eigenvalues":        eigenvalues,
        "ancient_coords":     ancient_coords,
        "modern_coords":      modern_coords,
        "pca_indiv_indices":  pca_indiv_indices,
        "n_pca_snps":         n_sub,
    }


# ---------------------------------------------------------------------------
# Report generation
# ---------------------------------------------------------------------------

def generate_report(
    distances: np.ndarray,
    count_valid: np.ndarray,
    individuals: list,
    anno: dict,
    output_path: Path,
) -> None:
    """Write a human-readable Markdown narrative of top ancient matches."""

    # Annotate each individual distance
    # Exclude modern reference populations (.DG suffix = diploid genotype modern samples)
    # and samples with no meaningful date (date_bp < 200 = born after ~1750 CE)
    rows = []
    for i, (ind, dist, n_snps) in enumerate(
        zip(individuals, distances, count_valid)
    ):
        if n_snps < MIN_INDIV_SNPS:
            continue
        rec = anno.get(ind.genetic_id)
        if rec is None:
            continue
        if rec.assessment not in ("PASS", "QUESTIONABLE"):
            continue
        # Skip modern reference populations — not the goal of ancient ancestry matching
        if rec.group_id and rec.group_id.endswith(".DG"):
            continue
        if rec.date_bp is not None and rec.date_bp < 200:
            continue

        date_ce = rec.date_ce
        date_str = (
            f"{abs(date_ce):.0f} {'BCE' if date_ce < 0 else 'CE'}"
            if date_ce is not None else "unknown"
        )

        rows.append({
            "index": i,
            "genetic_id": ind.genetic_id,
            "group_id": ind.population,
            "group_id_anno": rec.group_id,
            "locality": rec.locality,
            "political_entity": rec.political_entity,
            "lat": rec.lat,
            "lon": rec.lon,
            "date_ce": date_ce,
            "date_str": date_str,
            "date_bp": rec.date_bp,
            "y_hap": rec.best_y_haplogroup,
            "mt_hap": rec.valid_mtdna,
            "snps_compared": int(n_snps),
            "distance": float(dist),
        })

    rows.sort(key=lambda r: r["distance"])

    # Population-level summary
    pop_stats: dict[str, list[float]] = defaultdict(list)
    pop_meta: dict[str, dict] = {}
    for r in rows:
        pop = r["group_id_anno"]
        pop_stats[pop].append(r["distance"])
        if pop not in pop_meta:
            pop_meta[pop] = {
                "locality": r["locality"],
                "political_entity": r["political_entity"],
                "lat": r["lat"],
                "lon": r["lon"],
                "date_bp": r["date_bp"],
                "date_str": r["date_str"],
            }

    pop_ranked = sorted(
        [
            {
                "population": pop,
                "n_individuals": len(dists),
                "mean_distance": float(np.mean(dists)),
                "min_distance":  float(np.min(dists)),
                "median_distance": float(np.median(dists)),
                **pop_meta[pop],
            }
            for pop, dists in pop_stats.items()
        ],
        key=lambda r: r["mean_distance"],
    )

    lines = [
        "# Genome-wide Population Similarity Report — Individual 1",
        "",
        f"*Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}*",
        f"*SNPs compared: up to {int(count_valid.max()):,} per ancient individual*",
        f"*Ancient individuals with ≥{MIN_INDIV_SNPS:,} overlap SNPs: "
        f"{len(rows):,} of 17,629*",
        "",
        "---",
        "",
        "## Top 20 Closest Ancient Populations",
        "",
        "Ranked by mean allele-sharing distance (lower = more similar).",
        "",
        "| Rank | Population | Date | Location | N | Mean Dist |",
        "|------|-----------|------|---------|---|----------|",
    ]

    for rank, pop in enumerate(pop_ranked[:20], 1):
        lines.append(
            f"| {rank} | {pop['population']} | {pop['date_str']} "
            f"| {pop['political_entity']} | {pop['n_individuals']} "
            f"| {pop['mean_distance']:.4f} |"
        )

    lines += [
        "",
        "---",
        "",
        "## Top 30 Closest Ancient Individuals",
        "",
        "| Rank | ID | Population | Date | Location | Y-hap | MT-hap | SNPs | Distance |",
        "|------|-----|-----------|------|---------|-------|--------|------|---------|",
    ]

    for rank, r in enumerate(rows[:30], 1):
        lines.append(
            f"| {rank} | {r['genetic_id']} | {r['group_id_anno']} "
            f"| {r['date_str']} | {r['political_entity']} "
            f"| {r['y_hap'] or '?'} | {r['mt_hap'] or '?'} "
            f"| {r['snps_compared']:,} | {r['distance']:.4f} |"
        )

    lines += [
        "",
        "---",
        "",
        "## Interpretation",
        "",
        "**Distance metric:** Mean absolute allele-frequency difference per SNP.",
        "Range 0 (identical) to 1 (completely different). Typical same-population",
        "distances are 0.12–0.20; across-population distances are 0.20–0.35.",
        "",
        "**Closest populations** represent the ancient groups from which",
        "Individual 1 derives the most genome-wide ancestry.",
        "",
        "**Pseudo-haploid caution:** Ancient samples have one allele drawn at",
        "random per position, introducing sampling noise. Population-level means",
        "(across many individuals from the same group) are more reliable than",
        "individual-level distances.",
        "",
        "---",
        "",
        "*Next step: Step 1.4 — TMRCA estimation for closest haplogroup matches.*",
    ]

    output_path.write_text("\n".join(lines), encoding="utf-8")
    log.info("Wrote report: %s", output_path)
    return rows, pop_ranked


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    t0 = time.time()
    log.info("=== Step 1.3 — Genome-wide SNP Similarity and PCA ===")

    # ------------------------------------------------------------------
    # Resolve input paths (download from R2 to temp files if USE_R2)
    # ------------------------------------------------------------------
    _tmp_files: list[Path] = []
    if USE_R2:
        log.info("R2 mode: downloading input files for job %s", JOB_ID)
        _ind_path  = r2_client.download_to_temp(r2_client.IND_KEY,  '.ind')
        _anno_path = r2_client.download_to_temp(r2_client.ANNO_KEY, '.anno')
        _tmp_files = [_ind_path, _anno_path]
        # Handoff file: read locally if outputs are kept local; otherwise fetch from R2.
        if LOCAL_OUTPUTS:
            _overlap_path = OVERLAP_TSV
        else:
            _overlap_path = r2_client.download_to_temp(
                r2_client.output_key(JOB_ID, 'snp_overlap.tsv'), '.tsv'
            )
            _tmp_files.append(_overlap_path)
        geno = R2GenoFile.open(r2_client.GENO_KEY)
    else:
        from utils.parsers import resolve_local_aadr
        _aadr = resolve_local_aadr(DATA)
        global GENO_FILE, IND_FILE, ANNO_FILE
        GENO_FILE, IND_FILE, ANNO_FILE = _aadr["geno"], _aadr["ind"], _aadr["anno"]
        log.info("Local AADR resolved: %s", GENO_FILE.name)
        _overlap_path = OVERLAP_TSV
        _ind_path     = IND_FILE
        _anno_path    = ANNO_FILE
        geno = GenoFile.open(GENO_FILE)

    # ------------------------------------------------------------------
    # 1. Load overlap SNPs
    # ------------------------------------------------------------------
    log.info("Loading SNP overlap table...")
    geno_indices, modern_dosages, chroms = load_overlap(_overlap_path)
    n_overlap = len(geno_indices)
    log.info("Overlap: %d SNPs", n_overlap)

    # ------------------------------------------------------------------
    # 2. Load individual list and annotations
    # ------------------------------------------------------------------
    log.info("Loading individual list and annotations...")
    individuals = parse_ind_file(_ind_path)
    anno        = parse_anno_file(_anno_path)
    n_indiv     = len(individuals)
    log.info("%d individuals loaded", n_indiv)

    # ------------------------------------------------------------------
    # 3. GENO backend already opened above
    # ------------------------------------------------------------------
    log.info("GENO backend: %s", "R2 (range requests)" if USE_R2 else "local mmap")

    # ------------------------------------------------------------------
    # 4. Compute pairwise allele-sharing distances
    # ------------------------------------------------------------------
    log.info(
        "Computing allele-sharing distances: %d SNPs × %d individuals "
        "(chunk size %d, %d pseudo-haploid draws, autosomes only)...",
        n_overlap, n_indiv, CHUNK_SIZE, N_PSEUDO_DRAWS,
    )
    sum_diff, count_valid = compute_asd(geno, geno_indices, modern_dosages, chroms)

    # ASD = sum_diff / count_valid (avoid division by zero)
    distances = np.where(
        count_valid > 0,
        sum_diff / count_valid,
        np.nan,
    )

    log.info(
        "Distances computed. Valid comparisons: %d individuals with ≥%d SNPs",
        np.sum(count_valid >= MIN_INDIV_SNPS), MIN_INDIV_SNPS,
    )

    # ------------------------------------------------------------------
    # 5. Write pairwise distances TSV
    # ------------------------------------------------------------------
    log.info("Writing pairwise distances...")
    dist_path = OUTPUT / "pairwise_distances.tsv"
    with open(dist_path, "w") as fh:
        fh.write(
            "index\tgenetic_id\tpopulation\tgroup_id\tlocality\t"
            "political_entity\tlat\tlon\tdate_bp\tdate_display\t"
            "y_haplogroup\tmt_haplogroup\tsnps_compared\tasd_distance\t"
            "assessment\n"
        )
        for i, (ind, dist, n_snps) in enumerate(
            zip(individuals, distances, count_valid)
        ):
            rec = anno.get(ind.genetic_id)
            if rec is None:
                continue
            date_ce = rec.date_ce
            date_str = (
                f"{abs(date_ce):.0f} {'BCE' if date_ce < 0 else 'CE'}"
                if date_ce is not None else ""
            )
            fh.write(
                f"{i}\t{ind.genetic_id}\t{ind.population}\t{rec.group_id}\t"
                f"{rec.locality}\t{rec.political_entity}\t"
                f"{rec.lat or ''}\t{rec.lon or ''}\t"
                f"{rec.date_bp or ''}\t{date_str}\t"
                f"{rec.best_y_haplogroup}\t{rec.valid_mtdna}\t"
                f"{int(n_snps)}\t"
                f"{dist:.6f}\t{rec.assessment}\n"
            )
    log.info("Wrote %s", dist_path)

    # ------------------------------------------------------------------
    # 6. Population-level ranking
    # ------------------------------------------------------------------
    log.info("Computing population-level distances...")
    pop_dists: dict[str, list[float]] = defaultdict(list)
    pop_meta:  dict[str, dict] = {}

    for i, (ind, dist, n_snps) in enumerate(
        zip(individuals, distances, count_valid)
    ):
        if n_snps < MIN_INDIV_SNPS or np.isnan(dist):
            continue
        rec = anno.get(ind.genetic_id)
        if rec is None or rec.assessment not in ("PASS", "QUESTIONABLE"):
            continue
        # Exclude modern reference populations from population ranking
        if rec.group_id and rec.group_id.endswith(".DG"):
            continue
        if rec.date_bp is not None and rec.date_bp < 200:
            continue
        pop = rec.group_id or ind.population
        pop_dists[pop].append(float(dist))
        if pop not in pop_meta:
            date_ce = rec.date_ce
            pop_meta[pop] = {
                "locality": rec.locality,
                "political_entity": rec.political_entity,
                "lat": rec.lat,
                "lon": rec.lon,
                "date_bp": rec.date_bp,
                "date_display": (
                    f"{abs(date_ce):.0f} {'BCE' if date_ce < 0 else 'CE'}"
                    if date_ce is not None else ""
                ),
            }

    pop_ranked = sorted(
        [
            {
                "population": pop,
                "n_individuals": len(dists),
                "mean_distance": float(np.mean(dists)),
                "min_distance":  float(np.min(dists)),
                "median_distance": float(np.median(dists)),
                **pop_meta[pop],
            }
            for pop, dists in pop_dists.items()
            if len(dists) >= 1
        ],
        key=lambda r: r["mean_distance"],
    )

    pop_path = OUTPUT / "population_distances.tsv"
    with open(pop_path, "w") as fh:
        fh.write(
            "rank\tpopulation\tn_individuals\tmean_distance\t"
            "min_distance\tmedian_distance\t"
            "date_display\tlocality\tpolitical_entity\tlat\tlon\n"
        )
        for rank, pop in enumerate(pop_ranked, 1):
            fh.write(
                f"{rank}\t{pop['population']}\t{pop['n_individuals']}\t"
                f"{pop['mean_distance']:.6f}\t{pop['min_distance']:.6f}\t"
                f"{pop['median_distance']:.6f}\t"
                f"{pop['date_display']}\t{pop['locality']}\t"
                f"{pop['political_entity']}\t"
                f"{pop['lat'] or ''}\t{pop['lon'] or ''}\n"
            )
    log.info("Wrote %d populations to %s", len(pop_ranked), pop_path)

    # ------------------------------------------------------------------
    # 7. PCA
    # ------------------------------------------------------------------
    log.info("Identifying high-coverage individuals for PCA...")
    pca_mask = np.zeros(n_indiv, dtype=bool)
    for i, ind in enumerate(individuals):
        rec = anno.get(ind.genetic_id)
        if rec is None:
            continue
        if rec.assessment != "PASS":
            continue
        if rec.snps_1240k is not None and rec.snps_1240k >= MIN_PCA_SNPS:
            pca_mask[i] = True

    n_pca = pca_mask.sum()
    log.info("PCA: %d individuals with ≥%d 1240k SNPs", n_pca, MIN_PCA_SNPS)

    if n_pca >= 10:
        pca_result = compute_pca(
            geno, geno_indices, modern_dosages, pca_mask, PCA_N_COMPONENTS
        )

        # Write PCA coordinates
        pca_path = OUTPUT / "pca_coordinates.tsv"
        pca_indiv_idx = pca_result["pca_indiv_indices"]
        ancient_coords = pca_result["ancient_coords"]
        modern_coords  = pca_result["modern_coords"]
        eigenvalues    = pca_result["eigenvalues"]

        n_pc = len(eigenvalues)
        pc_headers = "\t".join(f"PC{i+1}" for i in range(n_pc))
        var_explained = "\t".join(f"{v:.4f}" for v in eigenvalues)

        with open(pca_path, "w") as fh:
            fh.write(f"type\tgenetic_id\tgroup_id\tdate_bp\tlat\tlon\t{pc_headers}\n")
            # Write ancient individuals
            for j, indiv_i in enumerate(pca_indiv_idx):
                ind = individuals[indiv_i]
                rec = anno.get(ind.genetic_id)
                coords = "\t".join(f"{ancient_coords[j, k]:.4f}" for k in range(n_pc))
                fh.write(
                    f"ancient\t{ind.genetic_id}\t"
                    f"{rec.group_id if rec else ind.population}\t"
                    f"{rec.date_bp if rec else ''}\t"
                    f"{rec.lat if rec else ''}\t"
                    f"{rec.lon if rec else ''}\t"
                    f"{coords}\n"
                )
            # Write Individual 1
            modern_pc = "\t".join(f"{modern_coords[k]:.4f}" for k in range(n_pc))
            fh.write(f"modern_indv_{OUTPUT_LABEL}\tIndividual_{OUTPUT_LABEL.upper()}\tModern_Individual_{OUTPUT_LABEL.upper()}\t0\t\t\t{modern_pc}\n")

        # Write variance explained as JSON
        with open(OUTPUT / "pca_variance_explained.json", "w") as fh:
            json.dump(
                {f"PC{i+1}": float(v) for i, v in enumerate(eigenvalues)},
                fh, indent=2,
            )

        log.info("Wrote PCA coordinates: %s", pca_path)
        log.info(
            "Individual 1 PCA coordinates: %s",
            ", ".join(f"PC{i+1}={modern_coords[i]:.3f}" for i in range(min(5, n_pc))),
        )
    else:
        log.warning("Insufficient individuals for PCA (%d). Skipping.", n_pca)
        pca_result = None

    # ------------------------------------------------------------------
    # 8. Generate narrative report
    # ------------------------------------------------------------------
    log.info("Generating report...")
    indiv_rows, pop_rows = generate_report(
        distances, count_valid, individuals, anno,
        OUTPUT / "top_matches_report.md",
    )

    # ------------------------------------------------------------------
    # 9. Print top results
    # ------------------------------------------------------------------
    print("\n" + "=" * 72)
    print("  INDIVIDUAL 1 — GENOME-WIDE SIMILARITY RESULTS")
    print("=" * 72)
    print(f"\n  Top 15 closest ancient POPULATIONS:\n")
    for rank, pop in enumerate(pop_rows[:15], 1):
        print(
            f"  {rank:>2}. {pop['mean_distance']:.4f}  "
            f"{pop['population'][:45]:<45}  "
            f"n={pop['n_individuals']}"
        )
    print(f"\n  Top 10 closest ancient INDIVIDUALS:\n")
    for rank, r in enumerate(indiv_rows[:10], 1):
        print(
            f"  {rank:>2}. {r['distance']:.4f}  "
            f"{r['genetic_id']:<15}  "
            f"{r['date_str']:<12}  "
            f"{r['group_id_anno']}"
        )
    if pca_result is not None:
        mc = pca_result["modern_coords"]
        print(f"\n  Individual 1 PCA position:")
        for i in range(min(5, len(mc))):
            var = pca_result["eigenvalues"][i]
            print(f"    PC{i+1} = {mc[i]:+.3f}  (explains {var:.1%} of variance)")
    print("=" * 72 + "\n")

    geno.close()

    # ------------------------------------------------------------------
    # Upload outputs to R2 (R2 mode only, unless LOCAL_OUTPUTS=1)
    # ------------------------------------------------------------------
    if USE_R2 and not LOCAL_OUTPUTS:
        output_files = [dist_path, pop_path, pca_path,
                        OUTPUT / "pca_variance_explained.json",
                        OUTPUT / "top_matches_report.md"]
        for local_file in output_files:
            if Path(local_file).exists():
                key = r2_client.output_key(JOB_ID, Path(local_file).name)
                r2_client.upload_file(local_file, key)
                log.info("Uploaded %s → R2:%s", Path(local_file).name, key)
    elif LOCAL_OUTPUTS:
        log.info("LOCAL_OUTPUTS=1 — skipping R2 upload, outputs remain in %s", OUTPUT)

    # Always clean up the temp AADR downloads when in R2 mode
    if USE_R2:
        for tmp in _tmp_files:
            try:
                tmp.unlink()
            except Exception:
                pass

    log.info("=== Step 1.3 complete in %.1f seconds ===", time.time() - t0)


if __name__ == "__main__":
    main()
