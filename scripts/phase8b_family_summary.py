"""
Phase 8b -- Family-level summary of the 386-candidate ranked list

With 386 non-exploratory candidates, a flat gene-by-gene table is not the most
informative way to present the result -- most candidates are members of a
small number of reference enzyme families, and recurrence of a FAMILY across
many independent genes/samples/cities is itself a stronger discovery signal
than any single gene (consistent with the manuscript's own Discussion logic
for single-gene recurrence). This aggregates by best-matching reference family.
"""
import csv
import re
from pathlib import Path
from collections import defaultdict

ROOT = Path(__file__).resolve().parents[1]
CAND = ROOT / "02_homology_screening" / "phase2_ranked_candidates.tsv"
OUT = ROOT / "08_candidate_prioritization" / "family_level_summary.tsv"

def family_of(ref_id):
    if ref_id.startswith("TierA|"):
        parts = ref_id.split("|")
        return f"OriginalCurated:{parts[1]}"  # e.g. OriginalCurated:FAcD, OriginalCurated:HAD, OriginalCurated:DEF2
    if ref_id.startswith("Toolkit|"):
        parts = ref_id.split("|")
        other_id = parts[3]
        # collapse individual accessions to family by stripping trailing numeric/chain suffixes
        fam = re.sub(r"_?\d*$", "", other_id.split("_")[0]) if "_" in other_id else other_id
        return f"Toolkit:{fam}"
    if ref_id.startswith("HMM:"):
        return f"HMMprofile:{ref_id[4:]}"
    return "Other"

def city_of(gene_id):
    m = re.match(r"^([MPS])", gene_id)
    return {"M": "Mardan", "P": "Peshawar", "S": "Swat"}.get(m.group(1)) if m else "?"

def main():
    with open(CAND, encoding="utf-8") as f:
        rows = [r for r in csv.DictReader(f, delimiter="\t")
                if r["confidence_tier"] in ("High-confidence", "Moderate")]

    fam_stats = defaultdict(lambda: dict(n_genes=0, n_high=0, n_mod=0, cities=set(),
                                          samples=set(), tpm_sum=0.0, best_pident=0.0,
                                          best_coverage=0.0))
    tpm_cols_cache = None
    for r in rows:
        fam = family_of(r["best_ref"])
        s = fam_stats[fam]
        s["n_genes"] += 1
        if r["confidence_tier"] == "High-confidence":
            s["n_high"] += 1
        else:
            s["n_mod"] += 1
        s["cities"].add(city_of(r["gene_id"]))
        s["tpm_sum"] += float(r["total_tpm_across_18_samples"])
        s["best_pident"] = max(s["best_pident"], float(r["pident"]))
        s["best_coverage"] = max(s["best_coverage"], float(r["coverage_pct"]))
        if tpm_cols_cache is None:
            tpm_cols_cache = [c for c in r.keys() if c.startswith("TPM_")]
        for c in tpm_cols_cache:
            if float(r[c]) > 0:
                s["samples"].add(c[4:])

    out_rows = []
    for fam, s in fam_stats.items():
        out_rows.append(dict(
            reference_family=fam, n_candidate_genes=s["n_genes"],
            n_high_confidence=s["n_high"], n_moderate=s["n_mod"],
            n_cities=len(s["cities"]), cities=",".join(sorted(s["cities"])),
            n_samples_with_detection=len(s["samples"]),
            total_tpm=round(s["tpm_sum"], 2),
            best_pident=round(s["best_pident"], 1), best_coverage=round(s["best_coverage"], 1),
        ))
    out_rows.sort(key=lambda r: -r["n_candidate_genes"])

    with open(OUT, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)

    print(f"Wrote {len(out_rows)} reference families to {OUT}\n")
    print(f"{'Family':35s}{'nGenes':>8s}{'nHC':>5s}{'nMod':>6s}{'nCities':>9s}{'nSamples':>10s}{'TotalTPM':>10s}")
    for r in out_rows:
        print(f"{r['reference_family']:35s}{r['n_candidate_genes']:>8d}{r['n_high_confidence']:>5d}"
              f"{r['n_moderate']:>6d}{r['n_cities']:>9d}{r['n_samples_with_detection']:>10d}"
              f"{r['total_tpm']:>10.1f}")

if __name__ == "__main__":
    main()
