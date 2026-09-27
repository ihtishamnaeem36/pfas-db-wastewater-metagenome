"""
Phase 5 -- Genomic context / neighborhood analysis

For every High-confidence/Moderate candidate, pull all genes sharing its contig
(sample_contig_N_*) from the corrected NR-protein FASTA, and annotate each neighbor
with its NR best-hit description (from NR.protein.taxonomy.txt, when available) to
check for co-located halogenated-compound-metabolism/transport/stress-response genes.
Contig length (bp, from total residues of all its genes as a proxy) and gene count are
reported as an explicit completeness/confidence criterion per the plan.
"""
import csv
import re
from pathlib import Path
from collections import defaultdict

ROOT = Path(__file__).resolve().parents[1]
P5 = ROOT.parents[2] / "project 5" / "Corrected Data"
CATALOG_FASTA = ROOT.parent / "Phase0_corrected" / "NR.protein.corrected.fa"
TAXO = P5 / "NR.protein.taxonomy.txt"
CANDIDATES = ROOT / "02_homology_screening" / "phase2_ranked_candidates.tsv"
OUT = ROOT / "05_genomic_neighborhood" / "neighborhood_context.tsv"

def contig_key(gene_id):
    m = re.match(r"^([A-Z]+\d)_contig_(\d+)_(\d+)$", gene_id)
    if not m:
        return None, None
    sample, contig_num, gene_num = m.groups()
    return f"{sample}_contig_{contig_num}", int(gene_num)

def main():
    with open(CANDIDATES, encoding="utf-8") as f:
        cands = [r for r in csv.DictReader(f, delimiter="\t")
                 if r["confidence_tier"] in ("High-confidence", "Moderate")]

    target_contigs = {}
    for c in cands:
        ck, gnum = contig_key(c["gene_id"])
        target_contigs[ck] = c["gene_id"]

    print(f"Scanning catalog for {len(target_contigs)} target contigs...")
    neighbors = defaultdict(list)  # contig -> [(gene_id, gene_num, seq_len)]
    name, seqlen = None, 0
    with open(CATALOG_FASTA, encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    ck, gnum = contig_key(name)
                    if ck in target_contigs:
                        neighbors[ck].append((name, gnum, seqlen))
                name = line[1:].split()[0]
                seqlen = 0
            else:
                seqlen += len(line.strip())
        if name is not None:
            ck, gnum = contig_key(name)
            if ck in target_contigs:
                neighbors[ck].append((name, gnum, seqlen))

    print("Loading NR best-hit descriptions for all neighbor genes...")
    wanted_genes = {gid for genes in neighbors.values() for gid, _, _ in genes}
    descs = {}
    with open(TAXO, encoding="utf-8") as f:
        next(f)
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if parts[0] in wanted_genes:
                descs[parts[0]] = parts[14] if len(parts) > 14 else "NA"
                if len(descs) == len(wanted_genes):
                    break

    HALOGEN_KEYWORDS = re.compile(
        r"dehalogen|halo|fluor|chlor|hydrolase|oxidoreductase|transport|efflux|"
        r"stress|resist|abc transporter|permease|reductase", re.I)

    rows = []
    for c in cands:
        ck, gnum = contig_key(c["gene_id"])
        genes = sorted(neighbors.get(ck, []), key=lambda x: x[1])
        total_len_aa = sum(g[2] for g in genes)
        n_genes = len(genes)
        neighbor_hits = []
        for gid, gn, glen in genes:
            if gid == c["gene_id"]:
                continue
            d = descs.get(gid, "NA")
            flag = "RELEVANT" if HALOGEN_KEYWORDS.search(d) else ""
            neighbor_hits.append(f"{gid}(gene#{gn},{glen}aa,{d}{' <'+flag+'>' if flag else ''})")
        rows.append(dict(
            gene_id=c["gene_id"], contig=ck, confidence_tier=c["confidence_tier"],
            n_genes_on_contig=n_genes,
            contig_length_proxy_aa=total_len_aa,
            n_relevant_neighbors=sum(1 for g in genes if g[0] != c["gene_id"]
                                      and HALOGEN_KEYWORDS.search(descs.get(g[0], ""))),
            neighbors="; ".join(neighbor_hits) if neighbor_hits else "(no other genes on contig)",
        ))

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    print(f"\nWrote neighborhood context for {len(rows)} candidates to {OUT}\n")
    for r in rows:
        print(f"  {r['gene_id']:28s} contig has {r['n_genes_on_contig']} genes, "
              f"{r['contig_length_proxy_aa']} aa total, "
              f"{r['n_relevant_neighbors']} halogen/transport/stress-relevant neighbor(s)")

if __name__ == "__main__":
    main()
