"""
Phase 6 prep -- extract High-confidence candidate sequences to a FASTA ready for
AlphaFold/ColabFold submission (see 06_structural_modeling/README.md for why the
folding step itself is not run in this environment), and write an empty, correctly
headered structural-evidence-table template for whoever runs Phase 6 to fill in.
"""
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CANDIDATES = ROOT / "02_homology_screening" / "phase2_ranked_candidates.tsv"
CATALOG_FASTA = ROOT.parent / "Phase0_corrected" / "NR.protein.corrected.fa"
OUT_FASTA = ROOT / "06_structural_modeling" / "candidates_for_alphafold.fasta"
OUT_TEMPLATE = ROOT / "06_structural_modeling" / "structural_evidence_table.tsv"

def main():
    with open(CANDIDATES, encoding="utf-8") as f:
        rows = [r for r in csv.DictReader(f, delimiter="\t")
                if r["confidence_tier"] == "High-confidence"]
    wanted = {r["gene_id"]: r for r in rows}

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

    OUT_FASTA.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_FASTA, "w", encoding="utf-8") as f:
        for gid, s in seqs.items():
            ref = wanted[gid]["best_ref"]
            f.write(f">{gid} best_ref={ref} pident={wanted[gid]['pident']} cov={wanted[gid]['coverage_pct']}\n")
            for i in range(0, len(s), 60):
                f.write(s[i:i+60] + "\n")

    with open(OUT_TEMPLATE, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["gene_id", "reference_structure_compared", "mean_pLDDT",
                    "active_site_pLDDT", "TM_score", "RMSD_angstrom",
                    "active_site_residue_conservation", "notes"])
        for gid in wanted:
            w.writerow([gid, "", "", "", "", "", "", ""])

    print(f"Wrote {len(seqs)} sequences to {OUT_FASTA}")
    print(f"Wrote empty template to {OUT_TEMPLATE}")

if __name__ == "__main__":
    main()
