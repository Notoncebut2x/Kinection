#!/usr/bin/env python3
"""
Step 1.2 — Haplogroup Assignment

Assigns Y-DNA (paternal) and mtDNA (maternal) haplogroups to Individual 1
using their AncestryDNA SNP data, then searches the ancient dataset's .anno
file to find ancient individuals who share or are ancestral to those haplogroups.

Method (mirrors Deep Maniot study):
  - Y-DNA: check derived state at ISOGG-defined marker positions; walk
    haplogroup tree from root to tip to find deepest supported assignment.
  - mtDNA: check key PhyloTree B17 positions; score haplogroup candidates
    by number of required derived markers observed.
  - Ancient matching: search .anno haplogroup columns for prefix matches,
    compute phylogenetic proximity score, and rank by closeness + date.

Outputs (written to output/step2/):
  ydna_haplogroup.json         — Y-DNA haplogroup assignment + evidence
  mtdna_haplogroup.json        — mtDNA haplogroup assignment + evidence
  ancient_haplogroup_matches.tsv — ancient individuals sharing haplogroup
  haplogroup_report.md         — human-readable narrative

Usage:
    python scripts/step2_haplogroup.py
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
# Paths (used in local mode; ignored when USE_R2=1)
# ---------------------------------------------------------------------------
ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data" / "input_data"
OUTPUT = ROOT / "output" / f"step2_{OUTPUT_LABEL}"
OUTPUT.mkdir(parents=True, exist_ok=True)

MODERN_INDV1  = DATA / "AncestryDNA_rn.txt"
ANNO_FILE     = DATA / "v62.0_1240k_public.anno"
YDNA_MARKERS  = ROOT / "scripts" / "data" / "ydna_markers.json"
MTDNA_MARKERS = ROOT / "scripts" / "data" / "mtdna_markers.json"
GENO_FILE     = DATA / "v62.0_1240k_public.geno"
IND_FILE      = DATA / "v62.0_1240k_public.ind"
OVERLAP_TSV   = ROOT / "output" / f"step1_{OUTPUT_LABEL}" / "snp_overlap.tsv"

# Min SNPs required for a Y or MT distance to be reported
MIN_Y_SNPS  = 10
MIN_MT_SNPS = 10

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-7s  %(message)s",
    datefmt="%H:%M:%S",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(OUTPUT / "step2.log", mode="w"),
    ],
)
log = logging.getLogger(__name__)

sys.path.insert(0, str(ROOT / "scripts"))
from utils.parsers import parse_ancestry_dna, parse_anno_file, parse_ind_file, complement, GenoFile
if USE_R2:
    from utils import r2_client
    from utils.r2_geno import R2GenoFile


# ---------------------------------------------------------------------------
# Utility: haplogroup prefix matching
# ---------------------------------------------------------------------------

def haplogroup_prefix_match(query: str, target: str) -> int:
    """
    Return the length of the common haplogroup prefix between query and target.

    E.g.:
      query="R1b1a1b", target="R1b1a1b1a1a" → 7 (full query matched)
      query="R1b",     target="R1b1a1b"     → 3
      query="J2a",     target="J2b"         → 2 (J2)

    Higher score = closer match.
    """
    if not query or not target:
        return 0
    # Strip notation artefacts and compare character-by-character
    q = query.strip().lstrip("-")
    t = target.strip().lstrip("-")
    # Find longest common prefix
    length = 0
    for a, b in zip(q, t):
        if a == b:
            length += 1
        else:
            break
    return length


def strip_haplogroup_prefix(hap: str) -> str:
    """
    Remove everything up to and including the first '-' in haplogroup notation.
    E.g. 'R-M269' → 'R', 'I-M423' → 'I', 'R1b' → 'R1b'
    """
    if "-" in hap:
        return hap.split("-")[0].strip()
    return hap.strip()


# ---------------------------------------------------------------------------
# Y-DNA haplogroup caller
# ---------------------------------------------------------------------------

class YDNACaller:
    """
    Calls Y-DNA haplogroup from AncestryDNA SNP data using ISOGG markers.

    Strategy:
      1. Load all Y chromosome SNPs (AncestryDNA chr 24).
      2. For each marker in the ISOGG marker database, check if the individual
         carries the derived allele (or its complement, accounting for strand).
      3. Walk the haplogroup tree from root to tip, collecting all derived
         markers. Assign the deepest haplogroup with at least one confirmed
         derived marker.

    Limitations:
      - AncestryDNA covers a limited set of Y markers — assignment may be
        at a broad haplogroup level (e.g., R1b rather than R-L21).
      - Pseudo-haploid data: Y calls should be homozygous; het calls are flagged.
      - Some marker rsIDs may not be on the array — position matching is used
        as a fallback.
    """

    def __init__(self, marker_db: dict):
        self.markers = marker_db["markers"]
        self.tree = marker_db.get("tree", {})

        # Build lookup: position → list of marker definitions
        self._pos_lookup: dict[int, list[dict]] = defaultdict(list)
        self._rsid_lookup: dict[str, dict] = {}
        for m in self.markers:
            self._pos_lookup[m["position"]].append(m)
            if m.get("rsid"):
                self._rsid_lookup[m["rsid"]] = m

    def call(self, y_snps: dict[str, tuple]) -> dict:
        """
        y_snps: dict mapping rsid → (position, allele1, allele2)
        Returns a result dict with haplogroup, confidence, and evidence.
        """
        # Build lookup structures for modern Y SNPs
        modern_by_rsid: dict[str, tuple] = {}    # rsid → (pos, a1, a2)
        modern_by_pos: dict[int, tuple] = {}      # pos → (rsid, a1, a2)
        for rsid, (pos, a1, a2) in y_snps.items():
            modern_by_rsid[rsid] = (pos, a1, a2)
            modern_by_pos[pos] = (rsid, a1, a2)

        log.info("Y-DNA: checking %d markers against %d modern Y SNPs",
                 len(self.markers), len(y_snps))

        derived_haplogroups: dict[str, list[dict]] = defaultdict(list)
        ancestral_haplogroups: dict[str, list[dict]] = defaultdict(list)
        not_found: list[dict] = []

        for marker in self.markers:
            rsid = marker.get("rsid", "")
            pos = marker["position"]
            ancestral = marker["ancestral"].upper()
            derived = marker["derived"].upper()
            hap = marker["haplogroup"]

            # Try rsID match first, then position fallback
            match = None
            if rsid and rsid in modern_by_rsid:
                pos_m, a1, a2 = modern_by_rsid[rsid]
                match = (a1, a2, "rsid")
            elif pos in modern_by_pos:
                rsid_m, a1, a2 = modern_by_pos[pos]
                match = (a1, a2, "position")

            if match is None:
                not_found.append(marker)
                continue

            a1, a2, match_type = match
            # Y chromosome: should be homozygous (haploid)
            obs_alleles = {a1.upper(), a2.upper()}

            # Y-DNA: check direct allele only (no complement flip).
            # Y markers are defined on the GRCh37 forward strand and AncestryDNA
            # reports Y SNPs on the same strand. Including complements causes
            # false-positive derived calls (e.g., C = complement of derived G
            # would incorrectly call G-derived haplogroups).
            is_derived = bool(obs_alleles & {derived})
            is_ancestral = bool(obs_alleles & {ancestral})

            evidence = {
                "marker": marker["name"],
                "haplogroup": hap,
                "rsid": rsid,
                "position": pos,
                "observed": f"{a1}/{a2}",
                "expected_derived": derived,
                "match_type": match_type,
                "notes": marker.get("notes", ""),
            }

            if is_derived and not is_ancestral:
                derived_haplogroups[hap].append(evidence)
                log.debug("  DERIVED  %-12s  %-8s  obs=%s/%s  derived=%s",
                          marker["name"], hap, a1, a2, derived)
            elif is_ancestral:
                ancestral_haplogroups[hap].append(evidence)
                log.debug("  ANCESTRAL %-12s %-8s  obs=%s/%s  derived=%s",
                          marker["name"], hap, a1, a2, derived)

        log.info(
            "Y-DNA markers found: %d derived calls across %d haplogroups; "
            "%d markers not on array",
            sum(len(v) for v in derived_haplogroups.values()),
            len(derived_haplogroups),
            len(not_found),
        )

        # Determine deepest haplogroup with derived evidence
        if not derived_haplogroups:
            return {
                "haplogroup": "Unknown",
                "confidence": "low",
                "derived_markers": [],
                "ancestral_markers": [],
                "markers_checked": len(self.markers) - len(not_found),
                "markers_not_on_array": len(not_found),
                "notes": "No derived markers detected. The individual may have "
                         "a haplogroup not covered by this marker set, or the "
                         "AncestryDNA array may not cover sufficient Y markers.",
            }

        # Score haplogroups: prefer deepest (longest haplogroup name = more specific)
        # Also consider number of confirming markers
        scored = []
        for hap, evidence_list in derived_haplogroups.items():
            score = len(hap) * 10 + len(evidence_list)
            scored.append((score, hap, evidence_list))
        scored.sort(reverse=True)

        best_hap = scored[0][1]
        best_evidence = scored[0][2]

        # Confidence based on marker coverage
        n_derived = sum(len(v) for v in derived_haplogroups.values())
        if n_derived >= 3:
            confidence = "high"
        elif n_derived == 2:
            confidence = "medium"
        else:
            confidence = "low"

        all_derived = [e for ev_list in derived_haplogroups.values() for e in ev_list]
        all_ancestral = [e for ev_list in ancestral_haplogroups.values() for e in ev_list]

        log.info("Y-DNA haplogroup: %s  (confidence=%s, %d derived markers)",
                 best_hap, confidence, n_derived)

        return {
            "haplogroup": best_hap,
            "confidence": confidence,
            "all_derived_haplogroups": {k: v for k, v in derived_haplogroups.items()},
            "derived_markers": all_derived,
            "ancestral_markers": all_ancestral[:20],  # truncate for readability
            "markers_checked": len(self.markers) - len(not_found),
            "markers_not_on_array": len(not_found),
        }


# ---------------------------------------------------------------------------
# mtDNA haplogroup caller
# ---------------------------------------------------------------------------

class MtDNACaller:
    """
    Calls mtDNA haplogroup from AncestryDNA mtDNA SNP data (chr 26).

    Uses PhyloTree B17 key positions. Scores candidate haplogroups by
    the proportion of required derived markers observed.
    """

    def __init__(self, marker_db: dict):
        self.markers = marker_db["markers"]
        self.haplogroup_logic = marker_db.get("haplogroup_logic", {})

        # Build position → marker lookup
        self._pos_lookup: dict[int, dict] = {}
        for m in self.markers:
            self._pos_lookup[m["position"]] = m

    def call(self, mt_snps: dict[str, tuple]) -> dict:
        """
        mt_snps: dict mapping rsid → (position, allele1, allele2)
        """
        # Build position lookup for modern MT SNPs
        modern_by_pos: dict[int, tuple] = {}
        for rsid, (pos, a1, a2) in mt_snps.items():
            modern_by_pos[pos] = (rsid, a1, a2)

        log.info("mtDNA: checking %d positions against %d modern MT SNPs",
                 len(self.markers), len(mt_snps))

        # Determine derived state at each known position
        derived_positions: set[int] = set()
        ancestral_positions: set[int] = set()
        observed_evidence: list[dict] = []

        for marker in self.markers:
            pos = marker["position"]
            ancestral = marker["ancestral"].upper()
            derived = marker["derived"].upper()

            if pos not in modern_by_pos:
                continue

            rsid, a1, a2 = modern_by_pos[pos]
            obs = a1.upper()  # mtDNA is haploid — a1 == a2

            derived_set = {derived, complement(derived)}
            ancestral_set = {ancestral, complement(ancestral)}

            is_derived = obs in derived_set
            is_ancestral = obs in ancestral_set

            evidence = {
                "position": pos,
                "rsid": rsid,
                "observed": obs,
                "ancestral": ancestral,
                "derived": derived,
                "state": "derived" if is_derived else ("ancestral" if is_ancestral else "unknown"),
                "haplogroup_relevance": marker["haplogroup"],
                "notes": marker.get("notes", ""),
            }
            observed_evidence.append(evidence)

            if is_derived:
                derived_positions.add(pos)
            elif is_ancestral:
                ancestral_positions.add(pos)

        log.info(
            "mtDNA: %d positions observed; %d derived, %d ancestral",
            len(observed_evidence), len(derived_positions), len(ancestral_positions),
        )

        # Score each haplogroup in the logic table
        haplogroup_scores: list[tuple[float, str, list[int], list[int]]] = []

        for hap, rule in self.haplogroup_logic.items():
            if hap.startswith("_"):
                continue
            required = rule["required"]
            n_required = len(required)
            n_found = sum(1 for pos in required if pos in derived_positions)
            n_ancestral = sum(1 for pos in required if pos in ancestral_positions)
            score = n_found / n_required if n_required > 0 else 0

            haplogroup_scores.append((
                score,
                hap,
                [pos for pos in required if pos in derived_positions],
                [pos for pos in required if pos in ancestral_positions],
            ))

        haplogroup_scores.sort(reverse=True)

        # Best assignment: highest score with at least one required marker found
        best_hap = None
        best_score = 0.0
        best_found = []
        best_missing = []

        for score, hap, found_pos, ancestral_pos in haplogroup_scores:
            if score > 0:
                best_hap = hap
                best_score = score
                required = self.haplogroup_logic[hap]["required"]
                best_found = found_pos
                best_missing = [p for p in required if p not in derived_positions]
                break

        if best_hap is None:
            best_hap = "L3 or pre-L3"
            best_score = 0.0
            log.warning("mtDNA: No haplogroup-specific markers detected.")
        else:
            log.info(
                "mtDNA haplogroup: %s  (score=%.2f, %d/%d required markers found)",
                best_hap, best_score,
                len(best_found),
                len(self.haplogroup_logic[best_hap]["required"]),
            )

        if best_score >= 1.0:
            confidence = "high"
        elif best_score >= 0.5:
            confidence = "medium"
        else:
            confidence = "low"

        # Top 5 candidates for reporting
        candidates = [
            {"haplogroup": h, "score": round(s, 3),
             "derived_markers_found": len(f),
             "required": self.haplogroup_logic[h]["required"]}
            for s, h, f, _ in haplogroup_scores
            if s > 0
        ][:5]

        return {
            "haplogroup": best_hap,
            "confidence": confidence,
            "score": round(best_score, 3),
            "required_markers_found": best_found,
            "required_markers_missing": best_missing,
            "all_candidates": candidates,
            "observed_positions": [e for e in observed_evidence if e["state"] != "unknown"],
            "positions_checked": len(observed_evidence),
        }


# ---------------------------------------------------------------------------
# Ancient haplogroup matcher
# ---------------------------------------------------------------------------

def match_ancient_haplogroups(
    y_haplogroup: str,
    mt_haplogroup: str,
    anno_records: dict,
    top_n: int = 50,
) -> list[dict]:
    """
    Search the ancient .anno file for individuals whose Y or mtDNA haplogroup
    matches (prefix) the modern individual's haplogroup assignment.

    Returns a list of match dicts sorted by:
      1. Match type (Y+MT > Y-only > MT-only)
      2. Haplogroup proximity score (longer prefix match = closer)
      3. Date (oldest to most recent)
    """
    matches = []

    # Normalise query haplogroups — strip notation like 'R-M269' → 'R'
    y_query_isogg = strip_haplogroup_prefix(y_haplogroup) if y_haplogroup else ""
    mt_query = mt_haplogroup.strip() if mt_haplogroup else ""

    for genetic_id, rec in anno_records.items():
        if rec.assessment == "IGNORE":
            continue

        # Get the best Y haplogroup for this ancient individual
        ancient_y_raw = rec.best_y_haplogroup
        ancient_y = strip_haplogroup_prefix(ancient_y_raw)
        ancient_mt = rec.valid_mtdna

        y_score = 0
        mt_score = 0
        y_match = False
        mt_match = False

        if y_query_isogg and ancient_y:
            y_score = haplogroup_prefix_match(y_query_isogg, ancient_y)
            if y_score >= 1:
                y_match = True

        if mt_query and ancient_mt:
            mt_score = haplogroup_prefix_match(mt_query, ancient_mt)
            if mt_score >= 1:
                mt_match = True

        if not y_match and not mt_match:
            continue

        match_type = (
            "Y+MT" if (y_match and mt_match)
            else ("Y" if y_match else "MT")
        )
        combined_score = y_score * 10 + mt_score  # Y weighted more

        date_ce = rec.date_ce
        date_display = (
            f"{abs(date_ce):.0f} {'BCE' if date_ce < 0 else 'CE'}"
            if date_ce is not None else "unknown"
        )

        matches.append({
            "genetic_id": genetic_id,
            "group_id": rec.group_id,
            "locality": rec.locality,
            "political_entity": rec.political_entity,
            "lat": rec.lat,
            "lon": rec.lon,
            "date_bp": rec.date_bp,
            "date_display": date_display,
            "molecular_sex": rec.molecular_sex,
            "ancient_y_haplogroup": ancient_y_raw,
            "ancient_mt_haplogroup": ancient_mt,
            "match_type": match_type,
            "y_proximity_score": y_score,
            "mt_proximity_score": mt_score,
            "combined_score": combined_score,
            "assessment": rec.assessment,
            "snps_1240k": rec.snps_1240k,
        })

    # Sort: match type priority, then combined score, then date
    match_priority = {"Y+MT": 0, "Y": 1, "MT": 2}
    matches.sort(key=lambda m: (
        match_priority.get(m["match_type"], 9),
        -m["combined_score"],
        m["date_bp"] if m["date_bp"] is not None else 0,
    ))

    log.info(
        "Ancient haplogroup matches: %d total "
        "(Y+MT=%d, Y-only=%d, MT-only=%d)",
        len(matches),
        sum(1 for m in matches if m["match_type"] == "Y+MT"),
        sum(1 for m in matches if m["match_type"] == "Y"),
        sum(1 for m in matches if m["match_type"] == "MT"),
    )

    return matches[:top_n]


# ---------------------------------------------------------------------------
# Report generator
# ---------------------------------------------------------------------------

def generate_report(
    y_result: dict,
    mt_result: dict,
    matches: list[dict],
    output_path: Path,
) -> None:
    """Write a human-readable Markdown report."""

    y_hap = y_result.get("haplogroup", "Unknown")
    mt_hap = mt_result.get("haplogroup", "Unknown")
    y_conf = y_result.get("confidence", "unknown")
    mt_conf = mt_result.get("confidence", "unknown")

    y_matches = [m for m in matches if m["match_type"] in ("Y", "Y+MT")]
    mt_matches = [m for m in matches if m["match_type"] in ("MT", "Y+MT")]
    both_matches = [m for m in matches if m["match_type"] == "Y+MT"]

    lines = [
        "# Haplogroup Assignment Report — Individual 1",
        "",
        f"*Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}*",
        "",
        "---",
        "",
        "## Summary",
        "",
        f"| | Haplogroup | Confidence | Markers Found |",
        f"|---|---|---|---|",
        f"| **Paternal (Y-DNA)** | {y_hap} | {y_conf} | {len(y_result.get('derived_markers', []))} derived |",
        f"| **Maternal (mtDNA)** | {mt_hap} | {mt_conf} | score={mt_result.get('score', 0):.2f} |",
        "",
        "---",
        "",
        "## Paternal Lineage (Y-DNA)",
        "",
        f"**Assigned haplogroup: {y_hap}**  (confidence: {y_conf})",
        "",
    ]

    if y_result.get("derived_markers"):
        lines += [
            "### Derived markers observed:",
            "",
            "| Marker | Haplogroup | Observed | Match Type |",
            "|--------|-----------|---------|------------|",
        ]
        for ev in y_result["derived_markers"][:15]:
            lines.append(
                f"| {ev['marker']} | {ev['haplogroup']} "
                f"| {ev['observed']} | {ev['match_type']} |"
            )
        lines.append("")

    # Ancient Y matches
    if y_matches:
        lines += [
            f"### Closest ancient Y-DNA matches (top {min(10, len(y_matches))}):",
            "",
            "| Ancient ID | Group | Date | Y Haplogroup | Proximity |",
            "|-----------|-------|------|-------------|-----------|",
        ]
        for m in y_matches[:10]:
            lines.append(
                f"| {m['genetic_id']} | {m['group_id']} "
                f"| {m['date_display']} | {m['ancient_y_haplogroup']} "
                f"| {m['y_proximity_score']} |"
            )
        lines += ["", "---", ""]
    else:
        lines += [
            "No ancient Y-DNA matches found in dataset with current haplogroup assignment.",
            "",
            "---",
            "",
        ]

    lines += [
        "## Maternal Lineage (mtDNA)",
        "",
        f"**Assigned haplogroup: {mt_hap}**  (confidence: {mt_conf})",
        "",
    ]

    if mt_result.get("required_markers_found"):
        found_positions = mt_result["required_markers_found"]
        lines += [
            f"Required defining positions found: {', '.join(str(p) for p in found_positions)}",
            "",
        ]

    if mt_result.get("all_candidates"):
        lines += [
            "### Top haplogroup candidates by marker score:",
            "",
            "| Haplogroup | Score | Markers Found |",
            "|-----------|-------|--------------|",
        ]
        for c in mt_result["all_candidates"]:
            lines.append(
                f"| {c['haplogroup']} | {c['score']:.3f} | {c['derived_markers_found']} |"
            )
        lines.append("")

    if mt_matches:
        lines += [
            f"### Closest ancient mtDNA matches (top {min(10, len(mt_matches))}):",
            "",
            "| Ancient ID | Group | Date | mtDNA Haplogroup | Proximity |",
            "|-----------|-------|------|-----------------|-----------|",
        ]
        for m in mt_matches[:10]:
            lines.append(
                f"| {m['genetic_id']} | {m['group_id']} "
                f"| {m['date_display']} | {m['ancient_mt_haplogroup']} "
                f"| {m['mt_proximity_score']} |"
            )
        lines += ["", "---", ""]

    if both_matches:
        lines += [
            "## Individuals Matching Both Y-DNA and mtDNA",
            "",
            "These ancient individuals share both the paternal and maternal "
            "haplogroup lineage with Individual 1 — the strongest possible match.",
            "",
            "| Ancient ID | Group | Date | Location | Y-hap | MT-hap |",
            "|-----------|-------|------|---------|-------|--------|",
        ]
        for m in both_matches[:10]:
            lines.append(
                f"| {m['genetic_id']} | {m['group_id']} "
                f"| {m['date_display']} | {m['political_entity']} "
                f"| {m['ancient_y_haplogroup']} | {m['ancient_mt_haplogroup']} |"
            )
        lines.append("")

    lines += [
        "---",
        "",
        "## Interpretation Notes",
        "",
        "- **Haplogroup proximity score**: number of haplogroup characters shared "
          "as a prefix. Score of 3 = first 3 characters match (e.g., R1b matches R1b1a).",
        "- **Y-DNA confidence** depends on the number of ISOGG markers present on "
          "the AncestryDNA V2.0 array. Many deep-branch markers require Big-Y or "
          "dedicated Y-sequencing for definitive sub-haplogroup assignment.",
        "- **mtDNA confidence** is based on PhyloTree B17 key positions. The "
          "AncestryDNA array covers ~195 mtDNA positions — sufficient for major "
          "haplogroup assignment but not full HVR1/HVR2 resolution.",
        "- Ancient individuals are pseudo-haploid (one allele drawn randomly), "
          "which reduces per-individual reliability but is robust in aggregate.",
        "",
        "---",
        "",
        "*Next step: Step 1.3 — Genome-wide SNP similarity and PCA projection.*",
        "(Requires valid v62.0_1240k_public.snp file — see step1 output for status.)",
    ]

    output_path.write_text("\n".join(lines), encoding="utf-8")
    log.info("Wrote report: %s", output_path)


# ---------------------------------------------------------------------------
# Chromosome-specific allele-sharing distance (Y and MT)
# ---------------------------------------------------------------------------

def load_chrom_snps(overlap_tsv: Path, chrom: str) -> tuple[np.ndarray, np.ndarray]:
    """
    Read the step1 SNP overlap table and return GENO indices + modern dosages
    for a single chromosome (e.g. "Y" or "MT").

    Y and MT are haploid in the modern individual — dosage should be 0 or 2
    (hom_ref or hom_alt). Heterozygous calls are retained but treated as 0.5,
    which is unusual for haploid sites and may reflect array noise.

    Returns:
        geno_indices:   int32 array of GENO row indices
        modern_dosages: int8 array (0/1/2, -1=missing)
    """
    geno_indices = []
    modern_dosages = []

    with open(overlap_tsv) as fh:
        fh.readline()  # skip header
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if parts[2] != chrom:
                continue
            geno_indices.append(int(parts[0]))
            d = parts[8]
            modern_dosages.append(int(d) if d not in ("NA", "") else -1)

    return (
        np.array(geno_indices,   dtype=np.int32),
        np.array(modern_dosages, dtype=np.int8),
    )


def compute_chrom_asd(
    geno: GenoFile,
    geno_indices: np.ndarray,
    modern_dosages: np.ndarray,
    indiv_mask: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute allele-sharing distance at a set of SNPs (Y or MT) between the
    modern individual and every ancient individual.

    Both sides are treated as pseudo-haploid:
      modern dosage 0  → allele 0.0 (hom ref)
      modern dosage 2  → allele 1.0 (hom alt)
      modern dosage 1  → allele 0.5 (het — rare on haploid chr, treated as 0.5)
      modern dosage -1 → missing (nan)

    Ancient GENO values: 0→0.0, 1→0.5, 2→1.0, 3→nan

    indiv_mask: optional boolean array (n_indiv,) — if given, only those
    individuals are included; others get count_valid=0.

    Returns:
        sum_diff:    float64 (n_indiv,) — accumulated |diff|
        count_valid: int64   (n_indiv,) — non-missing comparisons
    """
    n_indiv = geno.n_indiv
    n_snps  = len(geno_indices)

    sum_diff    = np.zeros(n_indiv, dtype=np.float64)
    count_valid = np.zeros(n_indiv, dtype=np.int64)

    geno_to_freq = np.array([0.0, 0.5, 1.0, np.nan], dtype=np.float32)
    modern_freq  = np.where(
        modern_dosages >= 0, modern_dosages / 2.0, np.nan
    ).astype(np.float32)

    # Read all SNP rows in one chunk (Y/MT counts are small — a few hundred SNPs)
    all_rows = geno.read_chunk(geno_indices)           # (n_snps, n_indiv) int8

    for i, row in enumerate(all_rows):
        ancient_freq = geno_to_freq[row.view(np.uint8) if row.dtype == np.int8 else row]

        mf = modern_freq[i]
        if np.isnan(mf):
            continue

        valid = ~np.isnan(ancient_freq)
        if indiv_mask is not None:
            valid &= indiv_mask

        diff = np.abs(mf - ancient_freq)
        diff[~valid] = 0.0

        sum_diff    += diff
        count_valid += valid.astype(np.int64)

    return sum_diff, count_valid


def rank_chrom_distances(
    sum_diff: np.ndarray,
    count_valid: np.ndarray,
    individuals: list,
    anno: dict,
    min_snps: int,
    chrom: str,
) -> list[dict]:
    """
    Compute ASD = sum_diff / count_valid, filter by min_snps and PASS/QUESTIONABLE
    assessment, and return rows sorted by distance ascending.
    """
    rows = []
    for i, ind in enumerate(individuals):
        n = int(count_valid[i])
        if n < min_snps:
            continue
        rec = anno.get(ind.genetic_id)
        if rec is None:
            continue
        if rec.assessment not in ("PASS", "QUESTIONABLE"):
            continue
        dist = float(sum_diff[i]) / n

        date_ce = rec.date_ce
        date_str = (
            f"{abs(date_ce):.0f} {'BCE' if date_ce < 0 else 'CE'}"
            if date_ce is not None else "unknown"
        )
        rows.append({
            "genetic_id":      ind.genetic_id,
            "group_id":        rec.group_id,
            "locality":        rec.locality,
            "political_entity": rec.political_entity,
            "date_bp":         rec.date_bp,
            "date_str":        date_str,
            "molecular_sex":   rec.molecular_sex,
            "y_haplogroup":    rec.best_y_haplogroup,
            "mt_haplogroup":   rec.valid_mtdna,
            "snps_compared":   n,
            "asd_distance":    dist,
        })

    rows.sort(key=lambda r: r["asd_distance"])
    return rows


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    t0 = time.time()
    log.info("=== Step 1.2 — Haplogroup Assignment ===")

    # ------------------------------------------------------------------
    # Resolve input paths (download from R2 to temp files if USE_R2)
    # ------------------------------------------------------------------
    _tmp_files: list[Path] = []
    # Modern individual DNA file always read from local disk — never stored in R2.
    # MODERN_DNA env var overrides the default path so different individuals can be analysed.
    _modern_path = Path(os.environ['MODERN_DNA']) if os.environ.get('MODERN_DNA') else MODERN_INDV1

    if USE_R2:
        log.info("R2 mode: downloading AADR reference files for job %s", JOB_ID)
        _anno_path = r2_client.download_to_temp(r2_client.ANNO_KEY, '.anno')
        _ind_path  = r2_client.download_to_temp(r2_client.IND_KEY,  '.ind')
        _tmp_files = [_anno_path, _ind_path]
        # Handoff file: read locally if outputs are kept local; otherwise fetch from R2.
        if LOCAL_OUTPUTS:
            _overlap_path = OVERLAP_TSV
        else:
            _overlap_path = r2_client.download_to_temp(
                r2_client.output_key(JOB_ID, 'snp_overlap.tsv'), '.tsv'
            )
            _tmp_files.append(_overlap_path)
        geno_backend = R2GenoFile.open(r2_client.GENO_KEY)
    else:
        _anno_path    = ANNO_FILE
        _ind_path     = IND_FILE
        _overlap_path = OVERLAP_TSV
        geno_backend  = None  # opened later below

    # ------------------------------------------------------------------
    # 1. Load marker databases
    # ------------------------------------------------------------------
    log.info("Loading haplogroup marker databases...")
    with open(YDNA_MARKERS) as fh:
        ydna_db = json.load(fh)
    with open(MTDNA_MARKERS) as fh:
        mtdna_db = json.load(fh)

    log.info("Y-DNA markers loaded: %d", len(ydna_db["markers"]))
    log.info("mtDNA markers loaded: %d", len(mtdna_db["markers"]))

    # ------------------------------------------------------------------
    # 2. Parse modern individual
    # ------------------------------------------------------------------
    log.info("Parsing Individual 1 AncestryDNA file...")
    modern_snps = parse_ancestry_dna(_modern_path)

    # Separate Y and MT SNPs
    y_snps: dict[str, tuple] = {}
    mt_snps: dict[str, tuple] = {}

    for rsid, snp in modern_snps.items():
        if snp.chrom == "Y":
            y_snps[rsid] = (snp.position, snp.allele1, snp.allele2)
        elif snp.chrom == "MT":
            mt_snps[rsid] = (snp.position, snp.allele1, snp.allele2)

    log.info("Individual 1: %d Y SNPs, %d MT SNPs", len(y_snps), len(mt_snps))

    # Also build position-keyed MT lookup (markers are position-based)
    mt_by_pos: dict[int, tuple] = {}
    for rsid, (pos, a1, a2) in mt_snps.items():
        mt_by_pos[pos] = (rsid, a1, a2)

    mt_snps_by_pos: dict[str, tuple] = {
        rsid: (pos, a1, a2)
        for rsid, (pos, a1, a2) in mt_snps.items()
    }

    # ------------------------------------------------------------------
    # 3. Call Y-DNA haplogroup
    # ------------------------------------------------------------------
    log.info("Calling Y-DNA haplogroup...")
    y_caller = YDNACaller(ydna_db)
    y_result = y_caller.call(y_snps)
    log.info("Y-DNA result: %s  (confidence=%s)",
             y_result["haplogroup"], y_result["confidence"])

    y_path = OUTPUT / "ydna_haplogroup.json"
    with open(y_path, "w") as fh:
        json.dump(y_result, fh, indent=2)
    log.info("Wrote %s", y_path)

    # ------------------------------------------------------------------
    # 4. Call mtDNA haplogroup
    # ------------------------------------------------------------------
    log.info("Calling mtDNA haplogroup...")
    mt_caller = MtDNACaller(mtdna_db)
    mt_result = mt_caller.call(mt_snps_by_pos)
    log.info("mtDNA result: %s  (confidence=%s)",
             mt_result["haplogroup"], mt_result["confidence"])

    mt_path = OUTPUT / "mtdna_haplogroup.json"
    with open(mt_path, "w") as fh:
        json.dump(mt_result, fh, indent=2)
    log.info("Wrote %s", mt_path)

    # ------------------------------------------------------------------
    # 5. Parse ancient annotations and find haplogroup matches
    # ------------------------------------------------------------------
    log.info("Parsing ancient individual annotations...")
    anno = parse_anno_file(_anno_path)

    log.info(
        "Searching %d ancient individuals for haplogroup matches...",
        len(anno),
    )
    matches = match_ancient_haplogroups(
        y_haplogroup=y_result["haplogroup"],
        mt_haplogroup=mt_result["haplogroup"],
        anno_records=anno,
        top_n=100,
    )

    # ------------------------------------------------------------------
    # 6. Write ancient match table
    # ------------------------------------------------------------------
    tsv_path = OUTPUT / "ancient_haplogroup_matches.tsv"
    if matches:
        header = (
            "genetic_id\tgroup_id\tlocality\tpolitical_entity\t"
            "lat\tlon\tdate_bp\tdate_display\tmolecular_sex\t"
            "ancient_y_haplogroup\tancient_mt_haplogroup\t"
            "match_type\ty_proximity_score\tmt_proximity_score\t"
            "combined_score\tassessment\tsnps_1240k\n"
        )
        with open(tsv_path, "w") as fh:
            fh.write(header)
            for m in matches:
                fh.write(
                    f"{m['genetic_id']}\t{m['group_id']}\t{m['locality']}\t"
                    f"{m['political_entity']}\t"
                    f"{m['lat'] or ''}\t{m['lon'] or ''}\t"
                    f"{m['date_bp'] or ''}\t{m['date_display']}\t"
                    f"{m['molecular_sex']}\t"
                    f"{m['ancient_y_haplogroup']}\t{m['ancient_mt_haplogroup']}\t"
                    f"{m['match_type']}\t{m['y_proximity_score']}\t"
                    f"{m['mt_proximity_score']}\t{m['combined_score']}\t"
                    f"{m['assessment']}\t{m['snps_1240k'] or ''}\n"
                )
        log.info("Wrote %d matches to %s", len(matches), tsv_path)
    else:
        log.warning("No ancient haplogroup matches found.")
        tsv_path.write_text("No matches found.\n")

    # ------------------------------------------------------------------
    # 7. Generate Markdown report
    # ------------------------------------------------------------------
    report_path = OUTPUT / "haplogroup_report.md"
    generate_report(y_result, mt_result, matches, report_path)

    # ------------------------------------------------------------------
    # 8. Y-chromosome allele-sharing distance (paternal lineage)
    # ------------------------------------------------------------------
    log.info("Computing Y-chromosome ASD (paternal lineage)...")
    y_geno_idx, y_modern_dos = load_chrom_snps(_overlap_path, "Y")
    log.info("Y chromosome: %d SNPs in overlap", len(y_geno_idx))

    individuals = parse_ind_file(_ind_path)
    geno = geno_backend if USE_R2 else GenoFile.open(GENO_FILE)

    # Only compare to ancient males for Y
    male_mask = np.array(
        [anno.get(ind.genetic_id, None) is not None
         and anno[ind.genetic_id].molecular_sex == "M"
         for ind in individuals],
        dtype=bool,
    )
    log.info("Ancient males available for Y comparison: %d", int(male_mask.sum()))

    y_sum, y_count = compute_chrom_asd(geno, y_geno_idx, y_modern_dos, indiv_mask=male_mask)
    y_rows = rank_chrom_distances(y_sum, y_count, individuals, anno, MIN_Y_SNPS, "Y")

    y_dist_path = OUTPUT / "ydna_distances.tsv"
    with open(y_dist_path, "w") as fh:
        fh.write(
            "rank\tgenetic_id\tgroup_id\tlocality\tpolitical_entity\t"
            "date_bp\tdate_display\ty_haplogroup\tsnps_compared\tasd_distance\tassessment\n"
        )
        for rank, r in enumerate(y_rows, 1):
            rec = anno.get(r["genetic_id"])
            fh.write(
                f"{rank}\t{r['genetic_id']}\t{r['group_id']}\t{r['locality']}\t"
                f"{r['political_entity']}\t{r['date_bp'] or ''}\t{r['date_str']}\t"
                f"{r['y_haplogroup']}\t{r['snps_compared']}\t"
                f"{r['asd_distance']:.6f}\t{rec.assessment if rec else ''}\n"
            )
    log.info("Wrote %d Y-ASD rows to %s", len(y_rows), y_dist_path)

    # ------------------------------------------------------------------
    # 9. Mitochondrial allele-sharing distance (maternal lineage)
    # ------------------------------------------------------------------
    log.info("Computing mtDNA ASD (maternal lineage)...")
    mt_geno_idx, mt_modern_dos = load_chrom_snps(_overlap_path, "MT")
    log.info("MT chromosome: %d SNPs in overlap", len(mt_geno_idx))

    mt_sum, mt_count = compute_chrom_asd(geno, mt_geno_idx, mt_modern_dos)
    mt_rows = rank_chrom_distances(mt_sum, mt_count, individuals, anno, MIN_MT_SNPS, "MT")

    mt_dist_path = OUTPUT / "mtdna_distances.tsv"
    with open(mt_dist_path, "w") as fh:
        fh.write(
            "rank\tgenetic_id\tgroup_id\tlocality\tpolitical_entity\t"
            "date_bp\tdate_display\tmt_haplogroup\tsnps_compared\tasd_distance\tassessment\n"
        )
        for rank, r in enumerate(mt_rows, 1):
            rec = anno.get(r["genetic_id"])
            fh.write(
                f"{rank}\t{r['genetic_id']}\t{r['group_id']}\t{r['locality']}\t"
                f"{r['political_entity']}\t{r['date_bp'] or ''}\t{r['date_str']}\t"
                f"{r['mt_haplogroup']}\t{r['snps_compared']}\t"
                f"{r['asd_distance']:.6f}\t{rec.assessment if rec else ''}\n"
            )
    log.info("Wrote %d MT-ASD rows to %s", len(mt_rows), mt_dist_path)

    geno.close()

    # ------------------------------------------------------------------
    # Upload outputs to R2 (R2 mode only, unless LOCAL_OUTPUTS=1)
    # ------------------------------------------------------------------
    if USE_R2 and not LOCAL_OUTPUTS:
        for local_file in [y_path, mt_path, tsv_path, report_path, y_dist_path, mt_dist_path]:
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

    # ------------------------------------------------------------------
    # 10. Print top results to console
    # ------------------------------------------------------------------
    print("\n" + "=" * 70)
    print(f"  INDIVIDUAL 1 — HAPLOGROUP RESULTS")
    print("=" * 70)
    print(f"  Paternal (Y-DNA):  {y_result['haplogroup']}  [{y_result['confidence']} confidence]")
    print(f"  Maternal (mtDNA):  {mt_result['haplogroup']}  [{mt_result['confidence']} confidence]")
    print()
    if matches:
        print(f"  Top ancient matches ({len(matches)} found):")
        for m in matches[:10]:
            print(
                f"    [{m['match_type']:4s}] {m['ancient_y_haplogroup'] or '?':12s} / "
                f"{m['ancient_mt_haplogroup'] or '?':10s}  "
                f"{m['date_display']:15s}  {m['group_id']}"
            )
    print()
    print(f"  Top 10 paternal (Y-chr) matches — {len(y_rows)} ancient males compared:\n")
    for rank, r in enumerate(y_rows[:10], 1):
        print(
            f"    {rank:>2}. {r['asd_distance']:.4f}  "
            f"{r['genetic_id']:<18}  {r['date_str']:<12}  "
            f"{r['y_haplogroup']:<12}  {r['group_id']}"
        )
    print()
    print(f"  Top 10 maternal (mtDNA) matches — {len(mt_rows)} individuals compared:\n")
    for rank, r in enumerate(mt_rows[:10], 1):
        print(
            f"    {rank:>2}. {r['asd_distance']:.4f}  "
            f"{r['genetic_id']:<18}  {r['date_str']:<12}  "
            f"{r['mt_haplogroup']:<12}  {r['group_id']}"
        )
    print("=" * 70 + "\n")

    log.info("=== Step 1.2 complete in %.1f seconds ===", time.time() - t0)


if __name__ == "__main__":
    main()
