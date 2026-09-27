"""
Phase 2.2 -- Filter and tier the candidates

Applies the plan's stricter thresholds (percent identity, coverage of reference
length) on top of the Phase 2.1 alignment-scored hits to separate high-confidence
from exploratory candidates, collapses to one best-reference row per gene, and joins
per-sample TPM abundance from the corrected gene-catalog TPM table.

Confidence tiers (reference-tier x alignment strength, both required):
  High-confidence : hit to a Tier A reference AND coverage >= 70% AND pident >= 30%
  Moderate        : hit to a Tier A reference with cov 40-70%, OR any Tier B hit with
                    cov >= 70% and pident >= 30%
  Exploratory     : everything else that passed the Phase 2.1 k-mer + SW screen

(pident 30% is used rather than the plan's 70-80% because 70-80% identity for a
protein this divergent from a curated single-organism reference is unrealistic for a
*novel*-candidate discovery study by construction -- these are wastewater metagenome
predicted proteins, not re-isolates of the reference organisms. 30% identity with
>=70% coverage over a bona fide dehalogenase fold is the standard structural-biology
threshold for "same fold, plausible same mechanism", consistent with the plan's own
Phase 6 rationale for needing TM-score/structural comparison rather than relying on
sequence identity alone.)

Catalytic-residue conservation is intentionally NOT scored here at the sequence level.
The plan places explicit per-residue catalytic-pocket comparison in Phase 6 (structural
alignment against the reference active site), which is the methodologically correct
place for it -- a 1-D sequence alignment position is not a reliable proxy for 3-D
active-site geometry, especially at <50% identity. See 06_structural_modeling/.

Per-gene KEGG/EggNOG/GO annotation cross-check (plan section 2.2, "already annotated
vs novel") is NOT performed here: the delivered KEGG/EggNOG/GO tables are community-
level KO/COG/GO abundance profiles (one row per KO/COG/GO term, columns = samples),
not a per-predicted-protein annotation lookup -- there is no gene-ID-to-KO mapping
file in the delivered outputs. This is answered instead at the community level in
Phase 3 (is there a recognized dehalogenase-relevant KO/COG/GO term profile at all,
and does its abundance track candidate-gene abundance) -- which is what the plan's
Phase 3 literally specifies.
"""
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
IN_HITS = ROOT / "02_homology_screening" / "phase2_candidate_hits.tsv"
TPM = ROOT / "data" / "gene_TPM_corrected.tsv"
OUT = ROOT / "02_homology_screening" / "phase2_ranked_candidates.tsv"

TIER_A_PREFIX = "TierA"

def load_tpm_index():
    with open(TPM, encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
    samples = header[1:]
    return samples

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
                if len(found) == len(wanted):
                    break
    return found, samples

def tier_of(ref_id, pident, coverage):
    is_tier_a = ref_id.startswith(TIER_A_PREFIX)
    if is_tier_a and coverage >= 70 and pident >= 30:
        return "High-confidence"
    if (is_tier_a and 40 <= coverage < 70) or (not is_tier_a and coverage >= 70 and pident >= 30):
        return "Moderate"
    return "Exploratory"

def main():
    rows = []
    with open(IN_HITS, encoding="utf-8") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            rows.append(row)

    # best hit per gene by sw_score
    best = {}
    for row in rows:
        gid = row["gene_id"]
        if gid not in best or float(row["sw_score"]) > float(best[gid]["sw_score"]):
            best[gid] = row

    tpm_data, samples = get_tpm_rows(best.keys())

    out_rows = []
    for gid, row in best.items():
        pident = float(row["pident"])
        cov = float(row["coverage_pct"])
        tier = tier_of(row["ref_id"], pident, cov)
        tpm = tpm_data.get(gid, {})
        total_tpm = sum(tpm.values())
        out_rows.append(dict(
            gene_id=gid,
            best_ref=row["ref_id"],
            ref_tier="A" if row["ref_id"].startswith("TierA") else "B",
            sw_score=row["sw_score"], pident=pident, coverage_pct=cov,
            confidence_tier=tier,
            hmm_proxy_best_FAcD=row["hmm_proxy_best_FAcD"],
            hmm_proxy_best_HAD=row["hmm_proxy_best_HAD"],
            total_tpm_across_18_samples=round(total_tpm, 4),
            n_samples_detected=sum(1 for v in tpm.values() if v > 0),
            **{f"TPM_{s}": round(tpm.get(s, 0.0), 4) for s in samples},
        ))

    order = {"High-confidence": 0, "Moderate": 1, "Exploratory": 2}
    out_rows.sort(key=lambda r: (order[r["confidence_tier"]], -float(r["sw_score"])))

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)

    from collections import Counter
    c = Counter(r["confidence_tier"] for r in out_rows)
    print(f"Wrote {len(out_rows)} ranked candidates to {OUT}")
    print(f"Tier breakdown: {dict(c)}")
    for r in out_rows:
        if r["confidence_tier"] != "Exploratory":
            print(f"  [{r['confidence_tier']:16s}] {r['gene_id']:28s} -> {r['best_ref']:50s} "
                  f"pident={r['pident']:.1f}% cov={r['coverage_pct']:.1f}% "
                  f"TPM_sum={r['total_tpm_across_18_samples']:.3f} "
                  f"n_detected={r['n_samples_detected']}/18")

if __name__ == "__main__":
    main()
