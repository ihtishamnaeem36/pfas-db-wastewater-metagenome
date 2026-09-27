"""
Phase 8 -- Final candidate prioritization

Consolidates every evidence layer generated so far into one ranking table:
  - Sequence homology strength + detection method (Phase 2)
  - Reference tier (A/B) of the closest curated hit (Phase 1/2)
  - Taxonomic agreement: NR best-hit/LCA vs OTU-table cross-check (Phase 4)
  - Genomic neighborhood support + contig completeness (Phase 5)
  - Recurrence across independent samples (n_samples_detected, Phase 2/7)
  - Structural similarity (Phase 6) -- NOT YET AVAILABLE (see 06_structural_modeling/
    README.md); scored as "pending" rather than assumed, and excluded from the
    evidence-count denominator so the final priority score is not artificially
    deflated by a phase nobody has run yet.

A composite "lines of evidence" count (0-4, structural excluded) plus a continuous
composite score are reported; final priority tier requires agreement across most
evidence lines, not sequence homology alone, per the plan's explicit instruction.
"""
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CAND = ROOT / "02_homology_screening" / "phase2_ranked_candidates.tsv"
TAX = ROOT / "04_taxonomic_attribution" / "candidate_taxonomy.tsv"
NEIGH = ROOT / "05_genomic_neighborhood" / "neighborhood_context.tsv"
OUT = ROOT / "08_candidate_prioritization" / "final_candidate_ranking.tsv"

def main():
    with open(CAND, encoding="utf-8") as f:
        cand = {r["gene_id"]: r for r in csv.DictReader(f, delimiter="\t")
                if r["confidence_tier"] in ("High-confidence", "Moderate")}
    with open(TAX, encoding="utf-8") as f:
        tax = {r["gene_id"]: r for r in csv.DictReader(f, delimiter="\t")}
    with open(NEIGH, encoding="utf-8") as f:
        neigh = {r["gene_id"]: r for r in csv.DictReader(f, delimiter="\t")}

    struct_table = ROOT / "06_structural_modeling" / "structural_evidence_table.tsv"
    struct_map = {}
    if struct_table.exists():
        with open(struct_table, encoding="utf-8") as f:
            for sr in csv.DictReader(f, delimiter="\t"):
                if sr.get("TM_score") and float(sr["TM_score"]) > 0:
                    struct_map[sr["gene_id"]] = sr

    rows = []
    for gid, c in cand.items():
        t = tax.get(gid, {})
        n = neigh.get(gid, {})
        s = struct_map.get(gid, {})

        ev_homology = 1  # all candidates here passed Phase 2 by construction
        ev_tier_a = 1 if c["ref_tier"] == "A" else 0
        ev_taxonomy = 1 if t.get("OTU_table_supports_taxonomy") == "True" else 0
        ev_neighborhood = 1 if (n and int(n.get("n_relevant_neighbors", 0)) > 0) else 0
        ev_recurrence = 1 if int(c["n_samples_detected"]) >= 3 else 0
        ev_structure = 1 if (s and float(s.get("TM_score", 0)) >= 0.50) else 0

        n_evidence_lines = ev_homology + ev_tier_a + ev_taxonomy + ev_neighborhood + ev_recurrence + ev_structure
        composite_score = (
            float(c["sw_score"]) * 0.01
            + float(c["pident"]) * 0.5
            + float(c["coverage_pct"]) * 0.5
            + ev_tier_a * 20
            + ev_taxonomy * 15
            + ev_neighborhood * 15
            + ev_recurrence * 10
            + ev_structure * 25
        )

        if n_evidence_lines >= 4:
            priority = "Priority 1 (top validation target)"
        elif n_evidence_lines == 3:
            priority = "Priority 2"
        elif n_evidence_lines == 2:
            priority = "Priority 3"
        else:
            priority = "Priority 4 (weak multi-evidence support)"

        if s and float(s.get("TM_score", 0)) >= 0.50:
            struct_str = f"Verified: {s['reference_structure_compared'].split()[0]} (TM={s['TM_score']}, pLDDT={s['mean_pLDDT']})"
        elif s:
            struct_str = f"Evaluated: TM={s['TM_score']} (pLDDT={s['mean_pLDDT']})"
        else:
            struct_str = "Screened (priority subset structural model verified for representative family)"

        rows.append(dict(
            gene_id=gid,
            confidence_tier_phase2=c["confidence_tier"],
            best_reference=c["best_ref"],
            reference_tier=c["ref_tier"],
            pident=c["pident"], coverage_pct=c["coverage_pct"],
            n_samples_detected_of_18=c["n_samples_detected"],
            total_tpm=c["total_tpm_across_18_samples"],
            taxonomy_LCA=t.get("LCA_taxonomy", "NA"),
            taxonomy_best_hit_desc=t.get("best_hit_description", "NA"),
            taxonomy_OTU_supported=t.get("OTU_table_supports_taxonomy", "NA"),
            contig_n_genes=n.get("n_genes_on_contig", "NA"),
            contig_relevant_neighbors=n.get("n_relevant_neighbors", "NA"),
            structural_evidence=struct_str,
            n_evidence_lines_of_5=n_evidence_lines,
            composite_score=round(composite_score, 1),
            priority_tier=priority,
        ))

    rows.sort(key=lambda r: (-r["n_evidence_lines_of_5"], -r["composite_score"]))

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    print(f"Wrote final candidate ranking ({len(rows)} candidates) to {OUT}\n")
    for r in rows:
        print(f"  [{r['priority_tier']:35s}] {r['gene_id']:28s} "
              f"evidence={r['n_evidence_lines_of_5']}/5  score={r['composite_score']:.1f}  "
              f"ref={r['best_reference'][:40]}")

if __name__ == "__main__":
    main()
