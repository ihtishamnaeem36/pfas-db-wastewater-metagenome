"""
Phase 2.1a -- fast k-mer prefilter (DIAMOND/HMMER substitute, step 1 of 2)

No DIAMOND/HMMER/BLAST binary is installable in this sandboxed Windows environment
(no admin/conda/WSL toolchain available). This substitutes a standard seed-and-extend
strategy implemented in pure Python:
  1. (this script) k-mer (k=5) shared-seed prefilter across the full 777,124-protein
     catalog against the 8 curated reference sequences -- fast, coarse recall pass,
     equivalent in spirit to DIAMOND's seed index / HMMER's MSV filter.
  2. (phase2_align_score.py) Biopython Smith-Waterman (BLOSUM62) local alignment of
     every candidate that survives step 1 against its best-matching reference --
     gives real percent identity / coverage / alignment score, equivalent to the
     DIAMOND blastp output columns the plan calls for.
This two-stage design is what makes exhaustive pairwise alignment of 8 queries against
777k targets computationally tractable in pure Python (full SW would be ~1e11 cell
updates; the k-mer stage cuts the SW stage down to a few thousand candidate pairs).

Threshold: sequences sharing >=6 five-mers with a reference are kept (permissive,
recall-oriented, matching the plan's "goal at this stage is recall" instruction).
"""
import re
import time
from pathlib import Path
from collections import defaultdict

ROOT = Path(__file__).resolve().parents[1]
REF_FASTA = ROOT / "01_reference_database" / "search_reference.fasta"
CATALOG_FASTA = ROOT.parent / "Phase0_corrected" / "NR.protein.corrected.fa"
OUT = ROOT / "02_homology_screening" / "kmer_prefilter_hits.tsv"

K = 5
MIN_SHARED_KMERS = 5

def read_fasta(path):
    name, seq = None, []
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line.strip())
    if name is not None:
        yield name, "".join(seq)

def kmers(seq, k=K):
    return {seq[i:i+k] for i in range(len(seq) - k + 1)}

def main():
    t0 = time.time()
    refs = list(read_fasta(REF_FASTA))
    ref_kmers = {name: kmers(seq) for name, seq in refs}
    ref_lens = {name: len(seq) for name, seq in refs}
    print(f"Loaded {len(refs)} reference sequences: {[n for n,_ in refs]}")

    n_scanned = 0
    n_hits = 0
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", encoding="utf-8") as out:
        out.write("gene_id\tref_id\tshared_kmers\ttarget_len\tref_len\n")
        for gid, seq in read_fasta(CATALOG_FASTA):
            n_scanned += 1
            if len(seq) < 30:
                continue
            tk = kmers(seq)
            for rname, rk in ref_kmers.items():
                shared = len(tk & rk)
                if shared >= MIN_SHARED_KMERS:
                    out.write(f"{gid}\t{rname}\t{shared}\t{len(seq)}\t{ref_lens[rname]}\n")
                    n_hits += 1
            if n_scanned % 100000 == 0:
                elapsed = time.time() - t0
                print(f"  scanned {n_scanned:,} / 777,124  ({elapsed:.0f}s elapsed, "
                      f"{n_scanned/elapsed:.0f} seq/s, {n_hits} hits so far)")

    elapsed = time.time() - t0
    print(f"\nDone. Scanned {n_scanned:,} proteins in {elapsed:.0f}s.")
    print(f"{n_hits} candidate (gene, reference) pairs with >= {MIN_SHARED_KMERS} shared "
          f"{K}-mers written to {OUT}")

if __name__ == "__main__":
    main()
