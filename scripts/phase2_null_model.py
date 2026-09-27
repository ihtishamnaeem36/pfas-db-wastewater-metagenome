"""
Phase 2 validation -- permutation-null false-positive-rate estimate

Reviewers of a homology-screening discovery paper will (rightly) ask: at the chosen
k-mer seed threshold (k=5, >=5 shared pentapeptides) and the Phase 2.2 tier
thresholds (coverage/identity), what fraction of hits could arise by chance alone,
given the real amino-acid composition of this catalog?

Method: every one of the 777,124 real predicted proteins is residue-shuffled in
place (fixed seed=42, per-sequence Fisher-Yates shuffle -- preserves exact length
and amino-acid composition per sequence, destroys any real sequence homology/motif
structure) and the *entire* Phase 2 pipeline (k-mer prefilter -> Smith-Waterman
scoring -> tiering) is rerun on this decoy catalog, identically to the real run.
This is a whole-catalog permutation null, not a subsample -- the strongest version
of this control.

Empirical FDR at each tier = (decoy hits reaching that tier) / (real hits reaching
that tier). Reported in 02_homology_screening/phase2_null_model_report.txt and
folded into the Methods/Results text.
"""
import re
import time
import random
import csv
from pathlib import Path
from collections import defaultdict
from Bio import Align
from Bio.Align import substitution_matrices

ROOT = Path(__file__).resolve().parents[1]
REF_FASTA = ROOT / "01_reference_database" / "search_reference.fasta"
CATALOG_FASTA = ROOT.parent / "Phase0_corrected" / "NR.protein.corrected.fa"
OUT_DIR = ROOT / "02_homology_screening"
K = 5
MIN_SHARED_KMERS = 5
SEED = 42

def read_fasta(path):
    d = {}
    name, seq = None, []
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    d[name] = "".join(seq)
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line.strip())
    if name is not None:
        d[name] = "".join(seq)
    return d

def kmers(seq, k=K):
    return {seq[i:i+k] for i in range(len(seq) - k + 1)}

aligner = Align.PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
aligner.mode = "local"
aligner.open_gap_score = -11
aligner.extend_gap_score = -1

def align_stats(a, b):
    aln = aligner.align(a, b)[0]
    score = aln.score
    t, q = aln.aligned
    aln_len = sum(e - s for s, e in t)
    ident = 0
    for (ts, te), (qs, qe) in zip(t, q):
        for i in range(te - ts):
            if a[ts + i] == b[qs + i]:
                ident += 1
    pident = 100.0 * ident / aln_len if aln_len else 0.0
    return score, pident, aln_len

def tier_of(ref_id, pident, coverage):
    is_tier_a = ref_id.startswith("TierA")
    if is_tier_a and coverage >= 70 and pident >= 30:
        return "High-confidence"
    if (is_tier_a and 40 <= coverage < 70) or (not is_tier_a and coverage >= 70 and pident >= 30):
        return "Moderate"
    return "Exploratory"

def main():
    rng = random.Random(SEED)
    refs = read_fasta(REF_FASTA)
    ref_kmers = {name: kmers(seq) for name, seq in refs.items()}
    ref_lens = {name: len(seq) for name, seq in refs.items()}

    t0 = time.time()
    n_scanned = 0
    decoy_hits = []  # (decoy_id, ref_id, target_seq)
    print("Scanning residue-shuffled decoy catalog (whole-catalog permutation null)...")
    with open(CATALOG_FASTA, encoding="utf-8", errors="replace") as f:
        name, seq_lines = None, []
        def process(nm, sq):
            nonlocal n_scanned
            n_scanned += 1
            if len(sq) < 30:
                return
            chars = list(sq)
            rng.shuffle(chars)
            shuffled = "".join(chars)
            tk = kmers(shuffled)
            for rname, rk in ref_kmers.items():
                if len(tk & rk) >= MIN_SHARED_KMERS:
                    decoy_hits.append((nm, rname, shuffled))
            if n_scanned % 200000 == 0:
                print(f"  scanned {n_scanned:,} / 777,124  ({time.time()-t0:.0f}s, "
                      f"{len(decoy_hits)} decoy k-mer hits so far)")
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    process(name, "".join(seq_lines))
                name = line[1:].split()[0]
                seq_lines = []
            else:
                seq_lines.append(line.strip())
        if name is not None:
            process(name, "".join(seq_lines))

    print(f"\nDecoy k-mer prefilter: {len(decoy_hits)} (gene, ref) pairs from "
          f"{n_scanned:,} shuffled decoys in {time.time()-t0:.0f}s "
          f"(real catalog: 74 pairs from 777,124 real proteins)")

    tier_counts = defaultdict(int)
    best_per_gene = {}
    for gid, rid, sq in decoy_hits:
        score, pident, aln_len = align_stats(sq, refs[rid])
        cov = 100.0 * aln_len / ref_lens[rid] if ref_lens[rid] else 0
        tier = tier_of(rid, pident, cov)
        if gid not in best_per_gene or score > best_per_gene[gid][0]:
            best_per_gene[gid] = (score, pident, cov, tier, rid)

    for gid, (score, pident, cov, tier, rid) in best_per_gene.items():
        tier_counts[tier] += 1

    real_counts = {"High-confidence": 6, "Moderate": 8, "Exploratory": 43}

    report = ["PHASE 2 VALIDATION -- WHOLE-CATALOG PERMUTATION NULL MODEL", "=" * 65,
              f"\nDecoy catalog: all {n_scanned:,} real proteins, per-sequence residue-",
              f"shuffled (seed={SEED}), same k-mer (k={K}, seed>={MIN_SHARED_KMERS}) + ",
              "Smith-Waterman + tiering pipeline as the real Phase 2 run.\n",
              f"{'Tier':20s}{'Real hits':>12s}{'Decoy hits':>14s}{'Empirical FDR':>16s}"]
    for tier in ["High-confidence", "Moderate", "Exploratory"]:
        real = real_counts[tier]
        decoy = tier_counts.get(tier, 0)
        fdr = decoy / real if real else float("nan")
        report.append(f"{tier:20s}{real:>12d}{decoy:>14d}{fdr:>15.3f}")

    report.append(f"\nTotal decoy candidate genes reaching any tier: {len(best_per_gene)} "
                  f"(vs. 57 real)")
    report.append("\nInterpretation: an empirical FDR near 0 at the High-confidence/Moderate "
                  "tiers indicates the coverage+identity thresholds (Phase 2.2) are doing "
                  "real discriminative work beyond the permissive k-mer seed stage, not "
                  "just passing through composition-driven noise.")

    out_txt = OUT_DIR / "phase2_null_model_report.txt"
    out_txt.write_text("\n".join(report), encoding="utf-8")
    print("\n".join(report))
    print(f"\nWrote {out_txt}")

if __name__ == "__main__":
    main()
