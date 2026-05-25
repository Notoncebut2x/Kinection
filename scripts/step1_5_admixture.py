"""
Step 1.5 — Admixture Decomposition

Decomposes the modern individual's genome into proportions of ancient
source populations using constrained non-negative least squares (NNLS).

The classic 4-way / 5-way ancient European model:

  target_allele_freq  ≈  α_WHG · f_WHG  +  α_EEF · f_EEF
                       + α_Steppe · f_Steppe
                       + α_Levant · f_Levant
                       + α_Iran · f_Iran

subject to:  Σ αᵢ = 1   and   αᵢ ≥ 0

Sources are configured in SOURCES below. Block-bootstrap by chromosome
gives 95 % CIs for each proportion.

This is NOT a full qpAdm replacement — qpAdm uses outgroup f₄-statistics
and is robust to confounding. NNLS is interpretable and fast, but treats
sources as if they were homogeneous and ignores LD between SNPs.
Useful as a headline estimate; cite qpAdm if publishing.

Inputs:
  output/step1_rn/modern_indv_rn_encoded.npy   (modern dosages)
  output/step1_rn/snp_overlap.tsv              (SNP positions + geno_index)
  AADR reference data (local or R2)

Outputs:
  output/step1_5_rn/admixture_decomposition.json
  output/step1_5_rn/admixture_report.md
  output/step1_5_rn/source_coverage.tsv        (per-source diagnostics)
"""
from __future__ import annotations

import json
import logging
import os
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

USE_R2 = os.environ.get('USE_R2', '').lower() in ('1', 'true', 'yes')
JOB_ID = os.environ.get('JOB_ID', 'dev')
LOCAL_OUTPUTS = os.environ.get('LOCAL_OUTPUTS', '').lower() in ('1', 'true', 'yes')
# Suffix used for output and handoff paths; must match the value used in step 1.
OUTPUT_LABEL = os.environ.get('OUTPUT_LABEL', 'rn')

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data" / "input_data"
OUT1 = ROOT / "output" / f"step1_{OUTPUT_LABEL}"
OUTPUT = ROOT / "output" / f"step1_5_{OUTPUT_LABEL}"
OUTPUT.mkdir(parents=True, exist_ok=True)

# Local AADR resolved lazily in main() — works for any version (v62, v66, ...).
GENO_FILE: Path | None = None
IND_FILE:  Path | None = None
ANNO_FILE: Path | None = None

# Bootstrap config
N_BOOTSTRAP = int(os.environ.get('ADMIX_BOOTSTRAP', '200'))
MIN_SNPS_PER_SOURCE = 50_000   # SNPs must be present in this many source individuals

# Source population manifest — patterns are substring-matched against group_id
# from the AADR .anno file. Edit to add/remove sources.
SOURCES = {
    "WHG": {
        "description": "Western European Hunter-Gatherer (Mesolithic Europe)",
        "patterns": [
            "Germany_Mesolithic",
            "France_Mesolithic",
            "Spain_Mesolithic",
            "Netherlands_Doggerland_Mesolithic",
            "Belgium_Mesolithic",
            "Sweden_Mesolithic",
        ],
    },
    "EHG": {
        "description": "Eastern European Hunter-Gatherer",
        "patterns": [
            "Russia_YuzhniyOleniyOstrov_Mesolithic",
            "Russia_Minino_Mesolithic",
            "Russia_Mesolithic_Veretye",
            "Karelia_HG",
        ],
    },
    "EEF": {
        "description": "Anatolian / Early European Farmer",
        "patterns": [
            "Turkey_Marmara_Barcin_N",
            "Turkey_Central_Catalhoyuk_N_lc",
            "Turkey_Marmara_Mentese_N",
            "Turkey_Southeast_Cayonu_PPN",
        ],
    },
    "Steppe": {
        "description": "Steppe Pastoralist (Yamnaya / Afanasievo, Bronze Age)",
        "patterns": [
            "Russia_Afanasievo",
            "Russia_Samara_EBA_Yamnaya",
            "Russia_Caucasus_EBA_Yamnaya",
            "Ukraine_EBA_Yamnaya",
            "Russia_Remontnoye_EBA_Yamnaya",
            "Russia_Kalmykia_EBA_Yamnaya",
        ],
    },
    "Levant_N": {
        "description": "Levantine Neolithic Farmer",
        "patterns": [
            "Jordan_PPNB",
            "Israel_Natufian",
            "Israel_PPNB",
            "Cyprus_PPNB",
        ],
    },
    "Iran_N": {
        "description": "Iranian Neolithic / Caucasus Hunter-Gatherer",
        "patterns": [
            "Iran_GanjDareh_N",
            "Georgia_Kotias",
            "Georgia_Satsurblia",
        ],
    },
}

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-7s  %(message)s",
    datefmt="%H:%M:%S",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(OUTPUT / "step1_5.log", mode="w"),
    ],
)
log = logging.getLogger(__name__)

sys.path.insert(0, str(ROOT / "scripts"))
from utils.parsers import parse_ind_file, parse_anno_file, GenoFile
if USE_R2:
    from utils import r2_client
    from utils.r2_geno import R2GenoFile


# ---------------------------------------------------------------------------
# Source-population resolution
# ---------------------------------------------------------------------------

@dataclass
class Source:
    name: str
    description: str
    individuals: list  # list of utils.parsers.Individual
    geno_indices: np.ndarray  # column indices into the GENO matrix


def find_source_individuals(
    inds: list, anno: dict
) -> dict[str, Source]:
    """For each configured source, find AADR individuals whose group_id matches."""
    by_id = {ind.genetic_id: ind for ind in inds}
    sources: dict[str, Source] = {}

    for src_name, cfg in SOURCES.items():
        matched = []
        for genetic_id, rec in anno.items():
            if rec.group_id and any(p in rec.group_id for p in cfg["patterns"]):
                ind = by_id.get(genetic_id)
                if ind is not None:
                    matched.append(ind)
        if not matched:
            log.warning("Source %s: no individuals matched any pattern", src_name)
            continue
        indices = np.array([ind.index for ind in matched], dtype=np.int32)
        sources[src_name] = Source(
            name=src_name,
            description=cfg["description"],
            individuals=matched,
            geno_indices=indices,
        )
        log.info("Source %-8s : %d individuals", src_name, len(matched))
    return sources


# ---------------------------------------------------------------------------
# Allele frequency computation
# ---------------------------------------------------------------------------

def compute_source_freqs(
    geno_backend,
    sources: dict[str, Source],
    overlap_geno_indices: np.ndarray,
    chunk_size: int = 5_000,
) -> np.ndarray:
    """
    Read SNP rows for the overlap, compute per-source allele frequencies.
    Returns array of shape (n_sources, n_snps) with values in [0, 1] or NaN.
    """
    n_snps = len(overlap_geno_indices)
    src_names = list(sources.keys())
    n_src = len(src_names)
    freqs = np.full((n_src, n_snps), np.nan, dtype=np.float32)

    log.info("Computing source allele frequencies over %d SNPs ...", n_snps)
    t0 = time.time()

    for chunk_start in range(0, n_snps, chunk_size):
        chunk_end = min(chunk_start + chunk_size, n_snps)
        chunk_geno_idx = overlap_geno_indices[chunk_start:chunk_end]

        # (chunk_size, n_indiv_total) int8
        rows = geno_backend.read_chunk(chunk_geno_idx)
        # 0/1/2 = data, 3 = missing in EIGENSTRAT PACKGENO
        rows_u8 = rows.view(np.uint8)

        for si, src_name in enumerate(src_names):
            cols = sources[src_name].geno_indices
            sub = rows_u8[:, cols].astype(np.float32)
            sub[sub == 3] = np.nan  # missing
            # EIGENSTRAT PACKGENO encodes the count of the REFERENCE allele,
            # but step 1's snp_overlap.dosage counts the ALT allele. Invert
            # the source so both target and sources count the same allele:
            #   source_alt_freq = 1 - mean(geno)/2
            with np.errstate(invalid="ignore"):
                freqs[si, chunk_start:chunk_end] = 1.0 - (np.nanmean(sub, axis=1) / 2.0)

        if (chunk_start // chunk_size) % 10 == 0:
            log.info("  chunk %d/%d  (%.0fs elapsed)",
                     chunk_start // chunk_size + 1,
                     (n_snps + chunk_size - 1) // chunk_size,
                     time.time() - t0)

    log.info("Source allele frequencies computed in %.1fs", time.time() - t0)
    return freqs


# ---------------------------------------------------------------------------
# NNLS with sum-to-1 constraint
# ---------------------------------------------------------------------------

def fit_admixture(A: np.ndarray, b: np.ndarray) -> np.ndarray:
    """
    Solve: minimize ||A·α - b||²  subject to  α ≥ 0,  Σα = 1.

    A: (n_snps, n_sources) source allele frequencies
    b: (n_snps,)           target allele frequencies (modern individual)
    Returns α of length n_sources.
    """
    from scipy.optimize import minimize
    k = A.shape[1]

    def loss(alpha):
        return np.sum((A @ alpha - b) ** 2)

    def loss_grad(alpha):
        return 2.0 * A.T @ (A @ alpha - b)

    result = minimize(
        loss,
        x0=np.full(k, 1.0 / k),
        jac=loss_grad,
        method="SLSQP",
        bounds=[(0.0, 1.0)] * k,
        constraints={"type": "eq", "fun": lambda x: x.sum() - 1.0},
        options={"maxiter": 500, "ftol": 1e-8},
    )
    if not result.success:
        log.warning("Optimiser did not converge: %s", result.message)
    return np.clip(result.x, 0.0, 1.0)


def bootstrap_admixture(
    A: np.ndarray,
    b: np.ndarray,
    chroms: np.ndarray,
    n_iter: int,
    seed: int = 42,
) -> np.ndarray:
    """
    Block-bootstrap by chromosome. Returns (n_iter, n_sources) of fits.
    """
    rng = np.random.default_rng(seed)
    unique_chroms = np.unique(chroms)
    n_chroms = len(unique_chroms)
    k = A.shape[1]

    out = np.zeros((n_iter, k), dtype=np.float32)
    for i in range(n_iter):
        picks = rng.choice(unique_chroms, size=n_chroms, replace=True)
        # Build the resampled index by concatenating SNP indices for picked chroms
        idx_list = [np.where(chroms == c)[0] for c in picks]
        idx = np.concatenate(idx_list)
        out[i] = fit_admixture(A[idx], b[idx])
        if (i + 1) % 50 == 0:
            log.info("  bootstrap %d/%d", i + 1, n_iter)
    return out


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def write_report(
    proportions: dict[str, float],
    cis: dict[str, tuple[float, float]],
    sources: dict[str, Source],
    n_snps: int,
    residual: float,
    out_path: Path,
) -> None:
    lines = []
    lines.append("# Admixture Decomposition Report")
    lines.append("")
    lines.append(f"*SNPs used: {n_snps:,}  |  Residual: {residual:.4f}*")
    lines.append("")
    lines.append("## Estimated ancient ancestry components")
    lines.append("")
    lines.append("| Source | Proportion | 95% CI | Description |")
    lines.append("|---|---:|---:|---|")
    for name in sorted(proportions, key=lambda n: -proportions[n]):
        p = proportions[name]
        lo, hi = cis[name]
        src = sources[name]
        lines.append(
            f"| **{name}** | {p:.1%} | [{lo:.1%}, {hi:.1%}] | "
            f"{src.description} ({len(src.individuals)} individuals) |"
        )
    lines.append("")

    lines.append("## How to read this")
    lines.append("")
    lines.append(
        "Each percentage is the fraction of your genome best modelled as descending from "
        "that ancient source population. The five sources together capture the canonical "
        "Western Eurasian ancestry decomposition described in Lazaridis et al. and "
        "follow-up papers."
    )
    lines.append("")
    lines.append(
        "**Important caveats:**"
    )
    lines.append("")
    lines.append(
        "- NNLS is a simplified alternative to the qpAdm framework used in published "
        "ancient-DNA papers. Headline numbers will be in the same ballpark; exact "
        "decimals will differ."
    )
    lines.append(
        "- The Levant_N component captures *deep* Levantine farmer ancestry shared by "
        "many Mediterranean and Near Eastern populations — not Ashkenazi-specific "
        "ancestry per se. Modern Ashkenazi Jewish ancestry is roughly half-Levantine, "
        "half-European; the European half here is decomposed into WHG/EEF/Steppe."
    )
    lines.append(
        "- 95 % CIs are from chromosome-block bootstrap and reflect statistical noise, "
        "not biological uncertainty in the model itself."
    )
    lines.append("")
    out_path.write_text("\n".join(lines))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def load_overlap(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read snp_overlap.tsv → geno_indices, target dosages, chroms (autosomes only)."""
    geno_idx = []
    dosages = []
    chroms = []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        col = {name: i for i, name in enumerate(header)}
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            chrom = parts[col["chrom"]]
            if chrom not in [str(i) for i in range(1, 23)]:
                continue  # autosomes only
            d = parts[col["dosage"]]
            if d in ("", "-1", "NA"):
                continue
            geno_idx.append(int(parts[col["geno_index"]]))
            dosages.append(int(d))
            chroms.append(int(chrom))
    return (
        np.array(geno_idx, dtype=np.int32),
        np.array(dosages, dtype=np.int8),
        np.array(chroms, dtype=np.int8),
    )


def main() -> None:
    t0 = time.time()
    log.info("=== Step 1.5 — Admixture Decomposition ===")

    # Inputs
    _tmp_files: list[Path] = []
    if USE_R2:
        log.info("R2 mode: downloading reference files for job %s", JOB_ID)
        _ind_path  = r2_client.download_to_temp(r2_client.IND_KEY,  '.ind')
        _anno_path = r2_client.download_to_temp(r2_client.ANNO_KEY, '.anno')
        _tmp_files = [_ind_path, _anno_path]
        if LOCAL_OUTPUTS:
            _overlap_path = OUT1 / "snp_overlap.tsv"
        else:
            _overlap_path = r2_client.download_to_temp(
                r2_client.output_key(JOB_ID, "snp_overlap.tsv"), ".tsv"
            )
            _tmp_files.append(_overlap_path)
        geno = R2GenoFile.open(r2_client.GENO_KEY)
    else:
        from utils.parsers import resolve_local_aadr
        _aadr = resolve_local_aadr(DATA)
        global GENO_FILE, IND_FILE, ANNO_FILE
        GENO_FILE, IND_FILE, ANNO_FILE = _aadr["geno"], _aadr["ind"], _aadr["anno"]
        log.info("Local AADR resolved: %s", GENO_FILE.name)
        _ind_path = IND_FILE
        _anno_path = ANNO_FILE
        _overlap_path = OUT1 / "snp_overlap.tsv"
        geno = GenoFile.open(GENO_FILE)

    # Load overlap (target genotypes)
    log.info("Loading SNP overlap and target genotypes ...")
    geno_indices, target_dosages, chroms = load_overlap(_overlap_path)
    n_snps = len(geno_indices)
    log.info("Overlap: %d autosomal SNPs with valid target genotypes", n_snps)
    target_freqs = target_dosages.astype(np.float32) / 2.0  # [0, 1]

    # Resolve source individuals
    log.info("Resolving source individuals from AADR ...")
    inds = parse_ind_file(_ind_path)
    anno = parse_anno_file(_anno_path)
    sources = find_source_individuals(inds, anno)
    if len(sources) < 3:
        sys.exit("Too few source populations available — aborting.")

    # Compute source allele frequencies
    src_freqs = compute_source_freqs(geno, sources, geno_indices)

    # Filter SNPs with complete data across all sources
    valid_mask = ~np.isnan(src_freqs).any(axis=0)
    n_valid = int(valid_mask.sum())
    log.info("SNPs with complete source coverage: %d / %d", n_valid, n_snps)

    if n_valid < MIN_SNPS_PER_SOURCE:
        sys.exit(f"Too few SNPs ({n_valid}) with complete source coverage.")

    A = src_freqs[:, valid_mask].T.astype(np.float32)  # (n_snps, n_src)
    b = target_freqs[valid_mask]
    chroms_valid = chroms[valid_mask]

    # Fit
    log.info("Fitting NNLS admixture model ...")
    alpha = fit_admixture(A, b)
    residual = float(np.sqrt(np.mean((A @ alpha - b) ** 2)))

    log.info("Bootstrap (n=%d, block by chromosome) ...", N_BOOTSTRAP)
    boots = bootstrap_admixture(A, b, chroms_valid, n_iter=N_BOOTSTRAP)
    cis_lo = np.percentile(boots, 2.5, axis=0)
    cis_hi = np.percentile(boots, 97.5, axis=0)

    # Package outputs
    src_names = list(sources.keys())
    proportions = {n: float(a) for n, a in zip(src_names, alpha)}
    cis = {n: (float(lo), float(hi)) for n, lo, hi in zip(src_names, cis_lo, cis_hi)}

    log.info("=== Admixture proportions ===")
    for name in sorted(proportions, key=lambda n: -proportions[n]):
        p = proportions[name]
        lo, hi = cis[name]
        log.info("  %-10s  %5.1f%%  (95%% CI: %4.1f%% – %4.1f%%)",
                 name, p * 100, lo * 100, hi * 100)
    log.info("Model residual: %.4f", residual)

    # Write outputs
    json_path = OUTPUT / "admixture_decomposition.json"
    json_path.write_text(json.dumps({
        "proportions": proportions,
        "ci95": cis,
        "residual": residual,
        "n_snps_used": int(n_valid),
        "n_bootstrap": N_BOOTSTRAP,
        "sources": {name: {
            "description": src.description,
            "n_individuals": len(src.individuals),
        } for name, src in sources.items()},
    }, indent=2))
    log.info("Wrote %s", json_path)

    report_path = OUTPUT / "admixture_report.md"
    write_report(proportions, cis, sources, n_valid, residual, report_path)
    log.info("Wrote %s", report_path)

    # Diagnostics: per-source coverage
    cov_path = OUTPUT / "source_coverage.tsv"
    with open(cov_path, "w") as fh:
        fh.write("source\tdescription\tn_individuals\tindividual_ids\n")
        for name, src in sources.items():
            ids = ",".join(i.genetic_id for i in src.individuals)
            fh.write(f"{name}\t{src.description}\t{len(src.individuals)}\t{ids}\n")
    log.info("Wrote %s", cov_path)

    geno.close()

    # Upload to R2 unless LOCAL_OUTPUTS
    if USE_R2 and not LOCAL_OUTPUTS:
        for local_file in [json_path, report_path, cov_path]:
            key = r2_client.output_key(JOB_ID, local_file.name)
            r2_client.upload_file(local_file, key)
            log.info("Uploaded %s → R2:%s", local_file.name, key)
    elif LOCAL_OUTPUTS:
        log.info("LOCAL_OUTPUTS=1 — skipping R2 upload, outputs remain in %s", OUTPUT)

    if USE_R2:
        for tmp in _tmp_files:
            try:
                tmp.unlink()
            except Exception:
                pass

    log.info("=== Step 1.5 complete in %.1f seconds ===", time.time() - t0)


if __name__ == "__main__":
    main()
