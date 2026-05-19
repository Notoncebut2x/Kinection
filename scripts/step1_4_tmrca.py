#!/usr/bin/env python3
"""
Step 1.4 — TMRCA Estimation for Y-DNA haplogroup matches.

For the top Y-DNA haplogroup matches identified in step 2, estimate the time
to most recent common ancestor (TMRCA) between the modern individual and each
ancient sample, by counting pairwise Y-SNP differences and converting via the
Y-chromosome substitution rate.

Method:
  - Y is haploid in males. AADR PACKGENO encodes: 0 = hom for allele2
    (column 6 of .snp); 2 = hom for allele1 (column 5); 1 = het (noise on Y,
    treated as missing); 3 = missing. Step 1's "dosage" counts col-6 alleles
    in the modern sample, so modern dosage 0 corresponds to ancient geno 2
    (both carry allele1) — the encodings are inverted between the two. We
    compare in a normalised "allele1 or allele2" space.
  - Primary metric: k/L = fraction of called Y SNPs where the two differ.
    Comparable across matches; robust to rate-calibration uncertainty.
  - Secondary metric: approximate TMRCA in years. We use a calibrated per-
    panel-SNP rate mu_panel = 7e-6 / panel-SNP / year, anchored to the
    known R1b coalescence age (~20 ky) — the workplan's per-bp rate of
    0.74e-9 cannot be used directly because the 1240k Y panel is ascertained
    for polymorphic sites and has ~10,000x higher effective mutation rate
    per panel-SNP than a random Y bp. The calibration is an order-of-
    magnitude estimate; absolute TMRCAs are reliable only to within a
    factor of ~2.
  - TMRCA = k / (2 * mu_panel * L) years.  Poisson exact 95% CI on k.
  - Archaeological-date floor: TMRCA cannot be more recent than the ancient
    sample's date; flagged when violated.

mtDNA TMRCA is intentionally skipped: the 1240k panel contains zero mtDNA
SNPs, so we have no overlap between the modern array and AADR mt data here.
Future work: ingest AADR's separate mt-capture dataset.

Inputs (per OUTPUT_LABEL):
  output/step1_<label>/snp_overlap.tsv               (Y geno_indices, modern Y dosages)
  output/step2_<label>/ydna_haplogroup.json          (confidence gate)
  output/step2_<label>/ancient_haplogroup_matches.tsv (top Y matches)
  v62 .geno + .ind (local file or R2 via R2GenoFile)

Outputs (output/step1_4_<label>/):
  ydna_tmrca.tsv          — per-match TMRCA point estimate + 95% CI
  mtdna_tmrca.tsv         — stub (mt not available from 1240k)
  tmrca_timeline.json     — structured data for timeline visualisation
  tmrca_report.md         — human-readable narrative

Usage:
  python scripts/step1_4_tmrca.py
"""

from __future__ import annotations

import json
import logging
import os
import sys
import time
from pathlib import Path

import numpy as np
from scipy.stats import chi2

# ---------------------------------------------------------------------------
# Env config (matches step3 conventions)
# ---------------------------------------------------------------------------
USE_R2        = os.environ.get('USE_R2', '').lower() in ('1', 'true', 'yes')
JOB_ID        = os.environ.get('JOB_ID', 'dev')
LOCAL_OUTPUTS = os.environ.get('LOCAL_OUTPUTS', '').lower() in ('1', 'true', 'yes')
OUTPUT_LABEL  = os.environ.get('OUTPUT_LABEL', 'rn')

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
N_TOP_MATCHES   = 10          # how many top Y-haplogroup matches to compute TMRCA for
MIN_Y_OVERLAP   = 50          # minimum sites called in both to attempt TMRCA
# Calibrated per-panel-SNP rate, anchored to known R1b coalescence (~20 ky).
# The 1240k Y panel is ascertained for polymorphic sites, so this rate is
# ~10,000x larger than the per-bp Y substitution rate of 0.74e-9. Order-of-
# magnitude estimate only — absolute TMRCAs accurate within a factor of ~2.
Y_PANEL_RATE    = 7.0e-6      # per panel-SNP per year
CURRENT_YEAR    = 2026        # for converting BP -> calendar
# Sites with ancient het (geno==1) on Y are pseudo-het artefacts, treated as missing.

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data" / "input_data"
OUT1 = ROOT / "output" / f"step1_{OUTPUT_LABEL}"
OUT2 = ROOT / "output" / f"step2_{OUTPUT_LABEL}"
OUTPUT = ROOT / "output" / f"step1_4_{OUTPUT_LABEL}"
OUTPUT.mkdir(parents=True, exist_ok=True)

GENO_FILE   = DATA / "v62.0_1240k_public.geno"
IND_FILE    = DATA / "v62.0_1240k_public.ind"
ANNO_FILE   = DATA / "v62.0_1240k_public.anno"
OVERLAP_TSV = OUT1 / "snp_overlap.tsv"
YDNA_JSON   = OUT2 / "ydna_haplogroup.json"
MATCHES_TSV = OUT2 / "ancient_haplogroup_matches.tsv"

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-7s  %(message)s",
    datefmt="%H:%M:%S",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(OUTPUT / "step1_4.log", mode="w"),
    ],
)
log = logging.getLogger(__name__)

sys.path.insert(0, str(ROOT / "scripts"))
from utils.parsers import parse_ind_file, GenoFile
if USE_R2:
    from utils import r2_client
    from utils.r2_geno import R2GenoFile


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_modern_y(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """
    Read Y-chromosome rows from snp_overlap.tsv.

    Returns:
        y_geno_indices: int32 array, sorted, of .geno row indices for Y SNPs
        modern_y_norm: int8 array — normalised modern allele in AADR's
            allele1/allele2 space. 0 = allele1 (col 5), 1 = allele2 (col 6),
            -1 = missing/het.

    Encoding mapping (see module docstring):
        step1 dosage = count of allele2 in modern; dosage 0 = hom allele1,
        dosage 2 = hom allele2. To compare with AADR PACKGENO (where
        geno 2 = hom allele1, geno 0 = hom allele2) we map both into a
        shared {0=allele1, 1=allele2} space.
    """
    indices, alleles = [], []
    with open(path) as fh:
        fh.readline()  # header
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if parts[2] != "Y":
                continue
            geno_idx = int(parts[0])
            dosage_str = parts[8]
            if dosage_str in ("NA", "", "1"):
                allele = -1
            else:
                d = int(dosage_str)
                # dosage 0 -> hom allele1 -> 0; dosage 2 -> hom allele2 -> 1
                allele = 0 if d == 0 else 1 if d == 2 else -1
            indices.append(geno_idx)
            alleles.append(allele)
    idx = np.array(indices, dtype=np.int32)
    al  = np.array(alleles, dtype=np.int8)
    assert np.all(np.diff(idx) >= 0), "Y geno_indices must be sorted"
    return idx, al


def poisson_ci(k: int, alpha: float = 0.05) -> tuple[float, float]:
    """Exact Poisson 95% CI on a count k (Garwood/chi-square method)."""
    if k == 0:
        lo = 0.0
    else:
        lo = 0.5 * chi2.ppf(alpha / 2.0,     2 * k)
    hi = 0.5 * chi2.ppf(1 - alpha / 2.0, 2 * (k + 1))
    return float(lo), float(hi)


def tmrca_from_k_l(k: int, l: int) -> tuple[float, float, float]:
    """
    Returns (point, lo_95, hi_95) TMRCA in years given k diffs over L panel SNPs.
    Returns (nan, nan, nan) if L < MIN_Y_OVERLAP.
    """
    if l < MIN_Y_OVERLAP:
        return float("nan"), float("nan"), float("nan")
    denom = 2.0 * Y_PANEL_RATE * l
    k_lo, k_hi = poisson_ci(k)
    return k / denom, k_lo / denom, k_hi / denom


def write_stub(reason: str) -> None:
    """When Y haplogroup confidence is too low to attempt TMRCA."""
    log.warning("Skipping Y-DNA TMRCA: %s", reason)
    (OUTPUT / "ydna_tmrca.tsv").write_text(
        "genetic_id\tpopulation\tdate_bp\tn_y_sites\tn_diff\t"
        "tmrca_yr\ttmrca_lo_95\ttmrca_hi_95\tnote\n"
        f"#\t#\t#\t0\t0\tNA\tNA\tNA\t{reason}\n"
    )
    (OUTPUT / "tmrca_timeline.json").write_text(json.dumps({
        "label": OUTPUT_LABEL, "skipped": True, "reason": reason, "matches": []
    }, indent=2))
    write_mt_stub()
    (OUTPUT / "tmrca_report.md").write_text(
        f"# TMRCA Estimation — {OUTPUT_LABEL}\n\n"
        f"_Step 1.4 was skipped: **{reason}**_\n\n"
        f"Y-DNA TMRCA requires a confident haplogroup assignment from step 2 "
        f"to identify meaningful ancient comparators. mtDNA TMRCA is not yet "
        f"available because the 1240k panel contains no mtDNA SNPs.\n"
    )


def write_mt_stub() -> None:
    (OUTPUT / "mtdna_tmrca.tsv").write_text(
        "# mtDNA TMRCA not computed: the 1240k panel contains 0 mtDNA SNPs.\n"
        "# Future work: ingest AADR's separate mt-capture dataset.\n"
    )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    log.info("=== Step 1.4 — Y-DNA TMRCA Estimation ===")
    log.info("Label: %s", OUTPUT_LABEL)

    # ── 1. Confidence gate ──────────────────────────────────────────────
    if not YDNA_JSON.exists():
        write_stub("step 2 Y-DNA output not found")
        return
    ydna = json.loads(YDNA_JSON.read_text())
    hg, conf = ydna.get("haplogroup", ""), ydna.get("confidence", "low")
    log.info("Step 2 Y-haplogroup: %s (confidence: %s)", hg, conf)
    if conf == "low" or hg in ("", "Unknown"):
        write_stub(f"Y-haplogroup confidence too low (haplogroup={hg!r}, confidence={conf!r})")
        return

    # ── 2. Load top Y-matches from step 2 ───────────────────────────────
    if not MATCHES_TSV.exists():
        write_stub("step 2 haplogroup_matches.tsv not found")
        return
    matches: list[dict] = []
    with open(MATCHES_TSV) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            row = dict(zip(header, line.rstrip("\n").split("\t")))
            if "Y" not in row.get("match_type", ""):
                continue
            try:
                row["_y_score"] = int(row["y_proximity_score"])
                row["_combined"] = int(row["combined_score"])
                row["_date_bp"]  = float(row["date_bp"]) if row["date_bp"] else None
            except (KeyError, ValueError):
                continue
            matches.append(row)
    matches.sort(key=lambda r: (-r["_y_score"], -r["_combined"]))
    top = matches[:N_TOP_MATCHES]
    log.info("Loaded %d Y-matches from step 2, taking top %d", len(matches), len(top))
    if not top:
        write_stub("no Y-haplogroup matches in step 2 output")
        return

    # ── 3. Load modern Y SNPs ───────────────────────────────────────────
    y_indices, modern_y = load_modern_y(OVERLAP_TSV)
    log.info("Modern Y SNPs in overlap: %d (called: %d)",
             len(y_indices), int((modern_y >= 0).sum()))
    if len(y_indices) < MIN_Y_OVERLAP:
        write_stub(f"only {len(y_indices)} Y SNPs in overlap (need ≥{MIN_Y_OVERLAP})")
        return

    # ── 4. Open .ind + .geno ────────────────────────────────────────────
    _tmp_files: list[Path] = []
    if USE_R2:
        log.info("R2 mode: downloading AADR .ind ...")
        _ind_path = r2_client.download_to_temp(r2_client.IND_KEY, '.ind')
        _tmp_files.append(_ind_path)
        geno = R2GenoFile.open(r2_client.GENO_KEY)
    else:
        _ind_path = IND_FILE
        geno = GenoFile.open(GENO_FILE)
    individuals = parse_ind_file(_ind_path)
    id_to_col = {ind.genetic_id: i for i, ind in enumerate(individuals)}
    log.info("Loaded %d ancient individuals from .ind", len(individuals))

    # ── 5. Read Y rows for all ancients in one chunk ────────────────────
    log.info("Reading %d Y SNP rows for all %d individuals...",
             len(y_indices), geno.n_indiv)
    t0 = time.time()
    y_rows = geno.read_chunk(y_indices).view(np.uint8)  # (n_y, n_indiv)
    log.info("Read complete in %.1fs", time.time() - t0)
    geno.close()
    for p in _tmp_files:
        try:
            p.unlink()
        except OSError:
            pass

    # ── 6. Compute TMRCA per match ──────────────────────────────────────
    # Both modern_y and ancient_norm are in {0, 1} = {allele1, allele2} space.
    # AADR PACKGENO: geno 2 = hom allele1 -> 0; geno 0 = hom allele2 -> 1;
    # geno 1 (het) and geno 3 (missing) -> treated as missing on Y.
    modern_mask = (modern_y >= 0)
    results: list[dict] = []
    for m in top:
        gid = m["genetic_id"]
        col = id_to_col.get(gid)
        if col is None:
            log.warning("  %s not found in .ind, skipping", gid)
            continue
        ancient_col = y_rows[:, col]
        ancient_called = (ancient_col == 0) | (ancient_col == 2)
        # geno 2 -> allele1 (0); geno 0 -> allele2 (1)
        ancient_norm = np.where(ancient_col == 2, 0, 1).astype(np.int8)

        both = modern_mask & ancient_called
        L = int(both.sum())
        k = int(((modern_y != ancient_norm) & both).sum())
        diff_rate = (k / L) if L > 0 else float("nan")
        T, T_lo, T_hi = tmrca_from_k_l(k, L)

        # Archaeological-date floor: TMRCA cannot be more recent than the sample.
        date_bp = m.get("_date_bp")
        floor_violation = (
            date_bp is not None and not np.isnan(T) and T < date_bp
        )

        rec = {
            "genetic_id": gid,
            "population": m.get("group_id", ""),
            "ancient_y_haplogroup": m.get("ancient_y_haplogroup", ""),
            "locality": m.get("locality", ""),
            "lat": float(m["lat"]) if m.get("lat") else None,
            "lon": float(m["lon"]) if m.get("lon") else None,
            "date_bp": date_bp,
            "date_display": m.get("date_display", ""),
            "n_y_sites": L,
            "n_diff": k,
            "diff_rate": diff_rate if not np.isnan(diff_rate) else None,
            "tmrca_yr": T if not np.isnan(T) else None,
            "tmrca_lo_95": T_lo if not np.isnan(T_lo) else None,
            "tmrca_hi_95": T_hi if not np.isnan(T_hi) else None,
            "below_sample_age": floor_violation,
        }
        results.append(rec)
        if np.isnan(T):
            log.info("  %-22s  L=%4d  k=%3d  rate=NA   TMRCA=NA  (need ≥%d sites)",
                     gid, L, k, MIN_Y_OVERLAP)
        else:
            log.info("  %-22s  L=%4d  k=%3d  rate=%.4f  TMRCA≈%6.0f y  CI=[%.0f, %.0f]%s",
                     gid, L, k, diff_rate, T, T_lo, T_hi,
                     "  ⚠ below sample age" if floor_violation else "")

    # ── 7. Write outputs ────────────────────────────────────────────────
    log.info("Writing outputs...")

    # TSV
    tsv = OUTPUT / "ydna_tmrca.tsv"
    with open(tsv, "w") as fh:
        fh.write("genetic_id\tpopulation\tancient_y_haplogroup\tdate_bp\t"
                 "n_y_sites\tn_diff\tdiff_rate\ttmrca_yr\ttmrca_lo_95\ttmrca_hi_95\tnote\n")
        for r in results:
            note = "below_sample_age" if r["below_sample_age"] else ""
            date_str = str(r['date_bp']) if r['date_bp'] is not None else "NA"
            common = (
                f"{r['genetic_id']}\t{r['population']}\t{r['ancient_y_haplogroup']}\t"
                f"{date_str}\t{r['n_y_sites']}\t{r['n_diff']}\t"
            )
            if r["tmrca_yr"] is not None:
                fh.write(
                    common +
                    f"{r['diff_rate']:.4f}\t{r['tmrca_yr']:.0f}\t"
                    f"{r['tmrca_lo_95']:.0f}\t{r['tmrca_hi_95']:.0f}\t{note}\n"
                )
            else:
                fh.write(common + "NA\tNA\tNA\tNA\tinsufficient_sites\n")
    log.info("Wrote %s", tsv)

    # JSON for timeline viz
    timeline = {
        "label":      OUTPUT_LABEL,
        "skipped":    False,
        "method":     f"pairwise Y-SNP differences on 1240k panel, "
                      f"calibrated mu_panel={Y_PANEL_RATE:.1e}/panel-SNP/year, Poisson 95% CI",
        "modern_y_haplogroup":      hg,
        "modern_y_confidence":      conf,
        "n_modern_y_sites_called":  int((modern_y >= 0).sum()),
        "matches": results,
    }
    (OUTPUT / "tmrca_timeline.json").write_text(json.dumps(timeline, indent=2, default=str))
    log.info("Wrote %s", OUTPUT / "tmrca_timeline.json")

    write_mt_stub()

    # Markdown report
    report_lines = [
        f"# TMRCA Estimation — {OUTPUT_LABEL}",
        "",
        f"**Modern Y-haplogroup:** {hg} (confidence: {conf})",
        f"**Method:** Pairwise Y-SNP difference rate (`k/L`) on the 1240k panel "
        f"vs each ancient match. Approximate TMRCA in years uses a calibrated "
        f"per-panel-SNP rate of {Y_PANEL_RATE:.0e}/year, anchored to known "
        f"R1b coalescence (~20 ky).",
        f"**Modern Y SNPs called:** {int((modern_y >= 0).sum())} of "
        f"{len(modern_y)} Y positions in overlap.",
        "",
        "## Top Y-haplogroup matches",
        "",
        "| Ancient | Population | Y-hg | Date | Y sites | Diffs | k/L | TMRCA ≈ (yr) | 95% CI |",
        "|---|---|---|---|---:|---:|---:|---:|---|",
    ]
    for r in results:
        date_str = r["date_display"] or (f"{int(r['date_bp'])} BP" if r["date_bp"] else "—")
        if r["tmrca_yr"] is None:
            rate_str = "—"
            tmrca_str = "—"
            ci_str = "—"
        else:
            rate_str = f"{r['diff_rate']:.3f}"
            tmrca_str = f"{r['tmrca_yr']:,.0f}"
            ci_str = f"{r['tmrca_lo_95']:,.0f}–{r['tmrca_hi_95']:,.0f}"
            if r["below_sample_age"]:
                tmrca_str += " ⚠"
        report_lines.append(
            f"| {r['genetic_id']} | {r['population']} | {r['ancient_y_haplogroup']} | "
            f"{date_str} | {r['n_y_sites']} | {r['n_diff']} | {rate_str} | "
            f"{tmrca_str} | {ci_str} |"
        )
    report_lines += [
        "",
        "## How to read these numbers",
        "",
        "- **`k/L` (difference rate)** is the more robust metric — it directly",
        "  compares how often two Y-haplotypes differ across the called panel SNPs.",
        "  Lower = more closely related. This is comparable across matches and",
        "  doesn't depend on rate calibration.",
        "- **TMRCA in years** uses a per-panel-SNP mutation rate calibrated to known",
        "  R1b coalescence (~20 ky). It is an *order-of-magnitude* estimate;",
        "  absolute values are accurate to within roughly a factor of 2.",
        "- The 95% CI reflects Poisson noise only — not rate-calibration uncertainty.",
        "  Add another ~30–50% uncertainty if quoting absolute ages.",
        "",
        "## Caveats",
        "",
        "- AncestryDNA arrays cover only ~900 Y-SNPs out of the 32,670 in the 1240k",
        "  panel, so confidence intervals are wide. Whole-Y sequencing (Big-Y, YFull)",
        "  gives much tighter TMRCAs.",
        "- The 1240k panel is *ascertained* for polymorphic sites, so the standard",
        "  per-bp Y mutation rate (0.74 × 10⁻⁹/bp/yr) cannot be applied directly —",
        "  hence the calibrated per-panel-SNP rate used here.",
        "- Matches flagged ⚠ have a point estimate younger than the ancient sample's",
        "  archaeological date, which is impossible. Almost always this means k≈0",
        "  and the upper CI is the meaningful number.",
        "- mtDNA TMRCA is not yet available: the 1240k SNP panel contains zero mtDNA",
        "  sites. Will require ingesting AADR's separate mt-capture dataset.",
    ]
    (OUTPUT / "tmrca_report.md").write_text("\n".join(report_lines) + "\n")
    log.info("Wrote %s", OUTPUT / "tmrca_report.md")

    log.info("=== Step 1.4 complete ===")


if __name__ == "__main__":
    main()
