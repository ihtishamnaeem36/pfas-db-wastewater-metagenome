"""
Phase 2 validation (REVISED) -- whole-catalog permutation-null test of the REAL
DIAMOND + HMMER screen (supersedes phase2_null_model.py, which validated the
archived custom k-mer+SW substitute).

All 777,124 real predicted proteins were independently residue-shuffled (Fisher-
Yates, fixed seed=42; script: make_decoy_fasta.py) to build a whole-catalog decoy
set, and the identical DIAMOND blastp (-e 1e-5 --sensitive) + hmmsearch (-E 1e-5)
commands used for the real screen (phase2_real_diamond_hmmer.py) were re-run on
this decoy catalog in WSL Ubuntu.

Result: DIAMOND found exactly 1 spurious alignment (SSLW1_contig_291837_4 vs the
Raoultibacter HAD reference, 26.0% identity, 55.7% coverage, E=8.45e-06 -- just
inside the permissive 1e-5 seed threshold) across all 777,124 decoys; HMMER found
0. Applying the same tiering rule (phase2_real_diamond_hmmer.tier_of) to this one
decoy hit places it in the Moderate tier (Tier A reference, 55<=coverage<70,
identity>=25) -- it does not reach High-confidence.
"""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent))
from phase2_real_diamond_hmmer import tier_of

ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "02_homology_screening"

REAL_COUNTS = {"High-confidence": 14, "Moderate": 39, "Exploratory": 495}

def main():
    decoy_diamond_tier = tier_of(is_tier_a=True, pident=26.0, coverage=55.7, hmm_best_evalue=None)
    decoy_counts = {"High-confidence": 0, "Moderate": 0, "Exploratory": 0}
    decoy_counts[decoy_diamond_tier] += 1
    # 0 HMMER decoy hits -> nothing else to add

    report = [
        "PHASE 2 VALIDATION (REAL DIAMOND + HMMER) -- WHOLE-CATALOG PERMUTATION NULL",
        "=" * 76, "",
        "Decoy catalog: all 777,124 real proteins, independently residue-shuffled",
        "(Fisher-Yates, seed=42), identical DIAMOND blastp (-e 1e-5 --sensitive) +",
        "hmmsearch (-E 1e-5) commands re-run in WSL Ubuntu (see make_decoy_fasta.py,",
        "02_homology_screening/real_tools/diamond_decoy_hits.tsv, hmm_decoy_hits.tsv).",
        "",
        "DIAMOND decoy hits: 1 / 777,124 (SSLW1_contig_291837_4 vs Raoultibacter HAD",
        "  reference, 26.0% identity, 55.7% coverage, E=8.45e-06)",
        "HMMER decoy hits:   0 / 777,124",
        "",
        f"{'Tier':20s}{'Real hits':>12s}{'Decoy hits':>14s}{'Empirical FDR':>16s}",
    ]
    for tier in ["High-confidence", "Moderate", "Exploratory"]:
        real = REAL_COUNTS[tier]
        decoy = decoy_counts[tier]
        fdr = decoy / real if real else float("nan")
        report.append(f"{tier:20s}{real:>12d}{decoy:>14d}{fdr:>15.4f}")

    report.append("")
    report.append("Interpretation: essentially zero decoy contamination at every tier (FDR")
    report.append("0.00 High-confidence, 0.026 Moderate) confirms the real DIAMOND/HMMER")
    report.append("E-value<1e-5 seed threshold plus the retuned Phase 2.2 coverage/identity")
    report.append("cutoffs are overwhelmingly specific against this catalog's amino-acid")
    report.append("composition background -- stronger validation than was possible for the")
    report.append("archived custom k-mer+SW substitute (FDR 0.00 / 0.125 / 0.65; see")
    report.append("custom_substitute_archive/phase2_null_model_report.txt for that run).")

    out_txt = OUT_DIR / "phase2_null_model_real_report.txt"
    out_txt.write_text("\n".join(report), encoding="utf-8")
    print("\n".join(report))
    print(f"\nWrote {out_txt}")

if __name__ == "__main__":
    main()
