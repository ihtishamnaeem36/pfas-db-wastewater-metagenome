"""
Phase 6a -- Sequence-based active-site residue conservation check (interim layer,
pending real Phase 6 structural modeling; see 06_structural_modeling/README.md)

WHY THIS EXISTS: Phase 6 proper (AlphaFold/ColabFold + TM-align, per the plan) requires
GPU compute that is not available in this environment. A live attempt to use the free
hosted ESMFold single-sequence structure API (api.esmatlas.com, no local GPU needed)
was made and failed (HTTP 504 / "Endpoint request timed out" on 3 consecutive tries,
including a 150 s timeout) -- logged in 06_structural_modeling/esmfold_attempt_log.txt.
This is reported as an infrastructure fact, not worked around by fabricating structures.

What this script does instead, honestly labeled as sequence-level (1-D), not
structural (3-D), evidence: DehH1/FA1 (UniProt Q1JU72) is the only Tier A reference
with an explicitly numbered catalytic triad in the curated citation table (Asp104-
His271-Asp128, numbering verified directly against the Q1JU72 sequence in this run:
position 104=D, 271=H, 128=D -- confirmed correct before use). For every candidate
whose best/co-scored reference is DehH1 or another FAcD-family Tier A hit, the
existing Phase 2 pairwise alignment is used to map these three reference positions
onto the candidate sequence and report whether a chemically compatible residue
(identical, or conservative substitution within {D,E} for the acidic pair, or a
non-conservative mismatch) occupies the aligned position.

This is NOT a substitute for real active-site pocket geometry (a 1-D alignment
position is not proof of 3-D proximity, especially below ~50% identity) -- it is
explicitly a lower tier of evidence than the Phase 6 structural check, and is
reported as such throughout.
"""
import csv
from pathlib import Path
from Bio import Align
from Bio.Align import substitution_matrices

ROOT = Path(__file__).resolve().parents[1]
REF_FASTA = ROOT / "01_reference_database" / "search_reference.fasta"
CATALOG_FASTA = ROOT.parent / "Phase0_corrected" / "NR.protein.corrected.fa"
CANDIDATES = ROOT / "02_homology_screening" / "phase2_ranked_candidates.tsv"
OUT = ROOT / "06_structural_modeling" / "phase6a_active_site_seqcheck.tsv"

DEHH1_ID = "TierA|FAcD|DehH1_FA1|Burkholderia_sp_FA1|Q1JU72"
CATALYTIC_POSITIONS = {104: "D", 271: "H", 128: "D"}  # verified against Q1JU72 above
ACIDIC = {"D", "E"}

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

def map_positions(cand_seq, ref_seq, ref_positions):
    """Align candidate (target) to ref (query), map 1-based ref positions -> candidate residue."""
    aln = aligner.align(cand_seq, ref_seq)[0]
    t_blocks, q_blocks = aln.aligned
    out = {}
    for pos in ref_positions:
        ref_idx0 = pos - 1
        found = None
        for (ts, te), (qs, qe) in zip(t_blocks, q_blocks):
            if qs <= ref_idx0 < qe:
                offset = ref_idx0 - qs
                cand_idx0 = ts + offset
                found = cand_seq[cand_idx0]
                break
        out[pos] = found  # None if position fell in a gap/unaligned region
    return out

def classify(expected, observed):
    if observed is None:
        return "not aligned (gap/outside alignment)"
    if observed == expected:
        return "IDENTICAL"
    if expected in ACIDIC and observed in ACIDIC:
        return "conservative (acidic-for-acidic)"
    if expected == "H" and observed in ("H", "K", "R"):
        return "conservative (basic-for-His)"
    return "NON-CONSERVATIVE MISMATCH"

def main():
    refs = read_fasta(REF_FASTA)
    ref_seq = refs[DEHH1_ID]
    print("Verifying reference catalytic-residue numbering against Q1JU72 sequence:")
    for pos, exp in CATALYTIC_POSITIONS.items():
        actual = ref_seq[pos - 1]
        status = "OK" if actual == exp else "MISMATCH -- DO NOT TRUST"
        print(f"  position {pos}: expected {exp}, sequence has {actual}  [{status}]")
        assert actual == exp, "Reference numbering verification failed"

    with open(CANDIDATES, encoding="utf-8") as f:
        # "FAcD" catches our original curated FAcD references; "acetate_F" catches the
        # Toolkit's EC 3.8.1.3 fluoroacetate dehalogenase entries (Q6NAM1, 3UMB/Q8XZN3)
        # -- same fold/mechanism family as FAcD, just not labeled "FAcD" in their header.
        # Other new high-spec families (reductive dehalogenase BCR, chlorobenzoyl-CoA
        # dehalogenase, triazine chlorohydrolase ATZA) are mechanistically distinct
        # folds -- mapping the DehH1 alpha/beta-hydrolase triad onto them would not be
        # meaningful, so they are intentionally excluded here, not overlooked. Must
        # exclude "HMM:family_FAcD.aln" pseudo-refs (no real alignment to map from).
        cands = [r for r in csv.DictReader(f, delimiter="\t")
                 if r["confidence_tier"] in ("High-confidence", "Moderate")
                 and ("FAcD" in r["best_ref"] or "acetate_F" in r["best_ref"])
                 and not r["best_ref"].startswith("HMM:")]

    wanted = {c["gene_id"] for c in cands}
    seqs = {}
    name, seq = None, []
    with open(CATALOG_FASTA, encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name in wanted:
                    seqs[name] = "".join(seq)
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line.strip())
        if name in wanted:
            seqs[name] = "".join(seq)

    rows = []
    for c in cands:
        gid = c["gene_id"]
        cseq = seqs.get(gid)
        if not cseq:
            continue
        mapped = map_positions(cseq, ref_seq, CATALYTIC_POSITIONS)
        classes = {pos: classify(CATALYTIC_POSITIONS[pos], mapped[pos]) for pos in CATALYTIC_POSITIONS}
        n_identical = sum(1 for v in classes.values() if v == "IDENTICAL")
        n_conservative_or_better = sum(1 for v in classes.values()
                                        if v in ("IDENTICAL",) or v.startswith("conservative"))
        rows.append(dict(
            gene_id=gid, best_ref=c["best_ref"], confidence_tier=c["confidence_tier"],
            pident_to_ref=c["pident"],
            Asp104_candidate_residue=mapped[104], Asp104_class=classes[104],
            His271_candidate_residue=mapped[271], His271_class=classes[271],
            Asp128_candidate_residue=mapped[128], Asp128_class=classes[128],
            n_identical_of_3=n_identical,
            n_conservative_or_better_of_3=n_conservative_or_better,
        ))

    OUT.parent.mkdir(parents=True, exist_ok=True)
    if rows:
        with open(OUT, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(rows)
        print(f"\nWrote {len(rows)} sequence-based active-site checks to {OUT}\n")
        for r in rows:
            print(f"  {r['gene_id']:28s} vs DehH1 triad: "
                  f"D104->{r['Asp104_candidate_residue']}({r['Asp104_class']}), "
                  f"H271->{r['His271_candidate_residue']}({r['His271_class']}), "
                  f"D128->{r['Asp128_candidate_residue']}({r['Asp128_class']})  "
                  f"[{r['n_conservative_or_better_of_3']}/3 conservative-or-better]")
    else:
        print("No FAcD-family High-confidence/Moderate candidates to check.")

if __name__ == "__main__":
    main()
