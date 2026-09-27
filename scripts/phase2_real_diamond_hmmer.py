"""
Phase 2 (REVISED) -- real DIAMOND + real HMMER screening

SUPERSEDES the custom k-mer + Smith-Waterman substitute (archived in
02_homology_screening/custom_substitute_archive/). DIAMOND and HMMER were not
installable via the normal package manager in the original Windows sandbox (no
admin rights), but this machine turned out to have WSL2 (Ubuntu) available, and
both tools could be obtained there WITHOUT root access:
  - diamond: official static Linux binary from the GitHub release, chmod +x, run.
  - hmmer: bioconda's linux-64 hmmer-3.4 build (hdbdd923_2, the libgcc/libstdc++-only
    variant with no MPI/GSL dependency) extracted directly from the .tar.bz2 conda
    package (apt-get download used to fetch a standalone zstd binary for a different,
    MPI-linked build first, then abandoned once the simpler no-MPI build was found via
    the Anaconda package-files API) -- no conda/mamba installation, no root, no
    internet-facing package manager beyond plain curl/apt-get download of individual
    .deb files for offline extraction.

Commands actually run (WSL Ubuntu, 12 CPU threads):
  diamond makedb --in search_reference.fasta -d pfas_refs
  diamond blastp -q NR.protein.corrected.fa -d pfas_refs -e 1e-5 --sensitive -k 0 \
      -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send \
             evalue bitscore qlen slen -p 12
    -> 679 pairwise alignments, 394 unique query genes, 33 s wall time
  hmmsearch --cpu 12 -E 1e-5 --tblout hmm_hits.tsv \
      PFAS_TierA_families.hmm NR.protein.corrected.fa
    -> 434 hits (family_FAcD.aln + family_HAD.aln profiles), 1.6 s wall time

This script merges both real outputs, assigns confidence tiers, and writes the
result in the exact schema of the original phase2_ranked_candidates.tsv so every
downstream script (phase4/5/6a/7/8/9) picks it up unchanged.
"""
import csv
from pathlib import Path
from collections import defaultdict

ROOT = Path(__file__).resolve().parents[1]
REF_FASTA = ROOT / "01_reference_database" / "search_reference.fasta"
DIAMOND = ROOT / "02_homology_screening" / "real_tools" / "diamond_hits.tsv"
HMM = ROOT / "02_homology_screening" / "real_tools" / "hmm_hits.tsv"
TPM = ROOT / "data" / "gene_TPM_corrected.tsv"
OUT = ROOT / "02_homology_screening" / "phase2_ranked_candidates.tsv"
OUT_RAW = ROOT / "02_homology_screening" / "phase2_candidate_hits.tsv"

TIER_A_PREFIX = "TierA"

def read_fasta_lens(path):
    lens = {}
    name, n = None, 0
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    lens[name] = n
                name = line[1:].split()[0]
                n = 0
            else:
                n += len(line.strip())
    if name is not None:
        lens[name] = n
    return lens

def load_diamond():
    """qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"""
    hits = defaultdict(list)
    with open(DIAMOND, encoding="utf-8") as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            qid, sid = p[0], p[1]
            pident, alnlen = float(p[2]), int(p[3])
            evalue, bitscore = float(p[10]), float(p[11])
            slen = int(p[13])
            coverage = 100.0 * alnlen / slen if slen else 0.0
            hits[qid].append(dict(ref_id=sid, pident=pident, aln_len=alnlen,
                                   coverage_pct=coverage, evalue=evalue, bitscore=bitscore,
                                   ref_len=slen))
    return hits

def load_hmm():
    """hmmsearch --tblout: target name, query name (family), E-value(full), score(full)"""
    hits = defaultdict(list)
    with open(HMM, encoding="utf-8") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            p = line.split()
            target, family = p[0], p[2]
            evalue, score = float(p[4]), float(p[5])
            hits[target].append(dict(family=family, evalue=evalue, score=score))
    return hits

def tier_of(is_tier_a, pident, coverage, hmm_best_evalue):
    # Thresholds recalibrated for REAL DIAMOND/HMMER output (much higher recall than the
    # archived custom k-mer+SW substitute): the naive Phase-2.2-draft thresholds (coverage
    # >=40%, HMM E<1e-10) produced 198 "Moderate" candidates from real DIAMOND's far more
    # sensitive search -- too broad to serve as the plan's "central results table" and
    # inconsistent with the plan's own instruction to apply STRICTER thresholds at this
    # step. Retuned against the actual coverage/identity/E-value distributions observed
    # (see 02_homology_screening/tier_threshold_tuning_notes.txt for the percentile
    # analysis behind these specific cutoffs).
    if is_tier_a and coverage >= 70 and pident >= 30:
        return "High-confidence"
    if hmm_best_evalue is not None and hmm_best_evalue < 1e-30:
        return "High-confidence"
    if is_tier_a and 55 <= coverage < 70 and pident >= 25:
        return "Moderate"
    if not is_tier_a and coverage >= 70 and pident >= 30:
        return "Moderate"
    if hmm_best_evalue is not None and hmm_best_evalue < 1e-15:
        return "Moderate"
    return "Exploratory"

def get_tpm_rows(gene_ids):
    wanted = set(gene_ids)
    found = {}
    with open(TPM, encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
        samples = header[1:]
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if parts[0] in wanted:
                found[parts[0]] = dict(zip(samples, map(float, parts[1:])))
    return found, samples

def main():
    ref_lens = read_fasta_lens(REF_FASTA)
    diamond_hits = load_diamond()
    hmm_hits = load_hmm()

    all_genes = set(diamond_hits) | set(hmm_hits)
    print(f"DIAMOND: {len(diamond_hits)} genes with hits (e<1e-5, sensitive mode)")
    print(f"HMMER:   {len(hmm_hits)} genes with hits (e<1e-5, family profiles)")
    print(f"Union:   {len(all_genes)} candidate genes total")

    raw_rows = []
    best_per_gene = {}
    for gid in all_genes:
        d_hits = sorted(diamond_hits.get(gid, []), key=lambda h: h["evalue"])
        h_hits = sorted(hmm_hits.get(gid, []), key=lambda h: h["evalue"])
        best_d = d_hits[0] if d_hits else None
        best_h = h_hits[0] if h_hits else None

        for dh in d_hits:
            raw_rows.append(dict(gene_id=gid, ref_id=dh["ref_id"], ref_len=dh["ref_len"],
                                  pident=dh["pident"], aln_len=dh["aln_len"],
                                  coverage_pct=round(dh["coverage_pct"], 1),
                                  evalue=dh["evalue"], bitscore=dh["bitscore"], method="DIAMOND"))
        for hh in h_hits:
            raw_rows.append(dict(gene_id=gid, ref_id=hh["family"], ref_len="", pident="",
                                  aln_len="", coverage_pct="", evalue=hh["evalue"],
                                  bitscore=hh["score"], method="HMMER"))

        if best_d:
            is_tier_a = best_d["ref_id"].startswith(TIER_A_PREFIX)
            pident, coverage = best_d["pident"], best_d["coverage_pct"]
            ref_id = best_d["ref_id"]
            sw_score = best_d["bitscore"]
        else:
            is_tier_a = True  # HMM-only hit is against the Tier A family profiles
            pident, coverage = 0.0, 0.0
            ref_id = f"HMM:{best_h['family']}"
            sw_score = best_h["score"] if best_h else 0.0

        hmm_best_evalue = best_h["evalue"] if best_h else None
        tier = tier_of(is_tier_a, pident, coverage, hmm_best_evalue)

        best_per_gene[gid] = dict(
            gene_id=gid, best_ref=ref_id, ref_tier="A" if is_tier_a else "B",
            sw_score=round(sw_score, 1), pident=round(pident, 1), coverage_pct=round(coverage, 1),
            confidence_tier=tier,
            hmm_proxy_best_FAcD=round(next((h["score"] for h in h_hits if h["family"].startswith("family_FAcD")), 0.0), 1),
            hmm_proxy_best_HAD=round(next((h["score"] for h in h_hits if h["family"].startswith("family_HAD")), 0.0), 1),
            diamond_evalue=best_d["evalue"] if best_d else None,
            hmm_evalue=hmm_best_evalue,
        )

    tpm_data, samples = get_tpm_rows(best_per_gene.keys())

    out_rows = []
    for gid, row in best_per_gene.items():
        tpm = tpm_data.get(gid, {})
        total_tpm = sum(tpm.values())
        out_rows.append(dict(
            gene_id=gid, best_ref=row["best_ref"], ref_tier=row["ref_tier"],
            sw_score=row["sw_score"], pident=row["pident"], coverage_pct=row["coverage_pct"],
            confidence_tier=row["confidence_tier"],
            hmm_proxy_best_FAcD=row["hmm_proxy_best_FAcD"],
            hmm_proxy_best_HAD=row["hmm_proxy_best_HAD"],
            total_tpm_across_18_samples=round(total_tpm, 4),
            n_samples_detected=sum(1 for v in tpm.values() if v > 0),
            **{f"TPM_{s}": round(tpm.get(s, 0.0), 4) for s in samples},
        ))

    order = {"High-confidence": 0, "Moderate": 1, "Exploratory": 2}
    out_rows.sort(key=lambda r: (order[r["confidence_tier"]], -float(r["sw_score"])))

    with open(OUT, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)

    with open(OUT_RAW, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(raw_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(raw_rows)

    from collections import Counter
    c = Counter(r["confidence_tier"] for r in out_rows)
    print(f"\nWrote {len(out_rows)} ranked candidates to {OUT}")
    print(f"Tier breakdown: {dict(c)}")
    print(f"\n(For comparison, the archived custom k-mer+SW substitute found: "
          f"{{'High-confidence': 6, 'Moderate': 8, 'Exploratory': 43}} from 57 total)")
    for r in out_rows:
        if r["confidence_tier"] != "Exploratory":
            print(f"  [{r['confidence_tier']:16s}] {r['gene_id']:28s} -> {r['best_ref'][:50]:50s} "
                  f"pident={r['pident']:.1f}% cov={r['coverage_pct']:.1f}% "
                  f"TPM_sum={r['total_tpm_across_18_samples']:.3f} "
                  f"n_detected={r['n_samples_detected']}/18")

if __name__ == "__main__":
    main()
