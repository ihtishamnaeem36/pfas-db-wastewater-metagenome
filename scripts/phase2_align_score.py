"""
Phase 2.1b -- local alignment scoring of k-mer-prefiltered candidates
(step 2 of 2; substitutes DIAMOND blastp output columns)

For every (gene, reference) pair that survived the k-mer prefilter, compute a real
Smith-Waterman local alignment (BLOSUM62, existence=-11, extend=-1 -- BLASTP defaults)
via Biopython's PairwiseAligner, and report percent identity, coverage of the reference
length, and raw alignment score, i.e. the quantities the plan's DIAMOND step needs.

Also scores every candidate against the two curated family alignments (family_FAcD,
family_HAD) as a simple profile-consistency check -- substituting hmmscan/hmmsearch
(no HMMER binary available in this environment; see README note in this script's
docstring tail for exactly what is and is not equivalent to real HMM scoring).

HMMER-substitute caveat: a true profile HMM scores position-specific conservation
across the whole family alignment with match/insert/delete states and trained
background frequencies. What is computed here is the max pairwise BLOSUM62 local
alignment score against each individual family member (i.e. best-in-family, not a
pooled profile score). This is a reasonable proxy for "does this candidate look like a
member of this family" but is NOT as sensitive to distant/divergent homologs as a real
HMM -- consistent with the caveat already flagged in the Phase 1 README for these thin
(n=3 per family) alignments. Documented as a methodological limitation, not concealed.
"""
import csv
from pathlib import Path
from Bio import Align
from Bio.Align import substitution_matrices

ROOT = Path(__file__).resolve().parents[1]
REF_FASTA = ROOT / "01_reference_database" / "search_reference.fasta"
FAM_FACD = ROOT / "01_reference_database" / "PFAS_Phase1_CuratedReferenceSet" / "family_FAcD.fasta"
FAM_HAD = ROOT / "01_reference_database" / "PFAS_Phase1_CuratedReferenceSet" / "family_HAD.fasta"
CATALOG_FASTA = ROOT.parent / "Phase0_corrected" / "NR.protein.corrected.fa"
HITS_IN = ROOT / "02_homology_screening" / "kmer_prefilter_hits.tsv"
OUT = ROOT / "02_homology_screening" / "phase2_candidate_hits.tsv"

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

aligner = Align.PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
aligner.mode = "local"
aligner.open_gap_score = -11
aligner.extend_gap_score = -1

def align_stats(a, b):
    aln = aligner.align(a, b)[0]
    score = aln.score
    t, q = aln.aligned  # aligned blocks, coordinates into a (target) and b (query/ref)
    aln_len = sum(e - s for s, e in t)
    ident = 0
    for (ts, te), (qs, qe) in zip(t, q):
        for i in range(te - ts):
            if a[ts + i] == b[qs + i]:
                ident += 1
    pident = 100.0 * ident / aln_len if aln_len else 0.0
    return score, pident, aln_len

def main():
    refs = read_fasta(REF_FASTA)
    fam_facd = read_fasta(FAM_FACD)
    fam_had = read_fasta(FAM_HAD)

    wanted_genes = set()
    pairs = []
    with open(HITS_IN, encoding="utf-8") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            wanted_genes.add(row["gene_id"])
            pairs.append((row["gene_id"], row["ref_id"]))

    print(f"Loading sequences for {len(wanted_genes)} candidate genes from catalog...")
    cand_seqs = {}
    with open(CATALOG_FASTA, encoding="utf-8", errors="replace") as f:
        name, seq = None, []
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name in wanted_genes:
                    cand_seqs[name] = "".join(seq)
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line.strip())
        if name in wanted_genes:
            cand_seqs[name] = "".join(seq)
    print(f"Recovered {len(cand_seqs)} / {len(wanted_genes)} candidate sequences.")

    rows = []
    for gid, rid in pairs:
        gseq = cand_seqs.get(gid)
        rseq = refs.get(rid)
        if not gseq or not rseq:
            continue
        score, pident, aln_len = align_stats(gseq, rseq)
        rlen = len(rseq)
        coverage = 100.0 * aln_len / rlen if rlen else 0.0

        fam_best_facd = max((align_stats(gseq, s)[0] for s in fam_facd.values()), default=0)
        fam_best_had = max((align_stats(gseq, s)[0] for s in fam_had.values()), default=0)

        rows.append(dict(
            gene_id=gid, ref_id=rid, ref_len=rlen, target_len=len(gseq),
            sw_score=round(score, 1), pident=round(pident, 1),
            aln_len=aln_len, coverage_pct=round(coverage, 1),
            hmm_proxy_best_FAcD=round(fam_best_facd, 1),
            hmm_proxy_best_HAD=round(fam_best_had, 1),
        ))
        print(f"  {gid:30s} vs {rid:55s} score={score:6.1f} pident={pident:5.1f}% cov={coverage:5.1f}%")

    rows.sort(key=lambda r: -r["sw_score"])
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"\nWrote {len(rows)} scored candidate pairs to {OUT}")

if __name__ == "__main__":
    main()
