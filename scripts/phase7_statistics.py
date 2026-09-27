"""
Phase 7 -- Comparative and statistical analysis

- Kruskal-Wallis (source type, city) + pairwise Wilcoxon on summed candidate-gene
  (High-confidence + Moderate tier) TPM abundance per sample, n=18, non-parametric
  given small n and compositional data.
- PCoA (classical MDS on Bray-Curtis dissimilarity, implemented directly with numpy
  since scikit-bio/scikit-learn are not installed in this environment) on the OTU
  table, colored by source/city, to see whether candidate-gene abundance tracks
  broader community compositional pattern.
- Explicit small-n caveat: n=2 per city x source cell -> exploratory/descriptive only.
"""
import csv
import re
import json
from pathlib import Path
from collections import defaultdict
import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
CANDIDATES = ROOT / "02_homology_screening" / "phase2_ranked_candidates.tsv"
OTU = ROOT.parent / "Phase0_corrected" / "All.Taxa.OTU.corrected.xls"
OUT = ROOT / "07_statistics"
OUT.mkdir(parents=True, exist_ok=True)

def parse_sample_meta(s):
    m = re.match(r"^([MPS])(HW|CW|SLW)(\d)$", s)
    city = {"M": "Mardan", "P": "Peshawar", "S": "Swat"}[m.group(1)]
    source = {"HW": "Hospital", "CW": "Community", "SLW": "Slaughterhouse"}[m.group(2)]
    return city, source, m.group(3)

def load_candidate_abundance():
    with open(CANDIDATES, encoding="utf-8") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    samples = [c for c in rows[0].keys() if c.startswith("TPM_")]
    sample_names = [c[4:] for c in samples]
    totals = defaultdict(float)
    for r in rows:
        if r["confidence_tier"] in ("High-confidence", "Moderate"):
            for c, s in zip(samples, sample_names):
                totals[s] += float(r[c])
    return sample_names, totals

def kruskal_wallis_and_pairwise(sample_names, totals, group_fn, label):
    groups = defaultdict(list)
    for s in sample_names:
        groups[group_fn(s)].append(totals[s])
    lines = [f"\n--- Kruskal-Wallis: candidate-gene abundance by {label} ---"]
    for g, vals in groups.items():
        lines.append(f"  {g:15s} n={len(vals)}  values={[round(v,3) for v in vals]}")
    stat, p = stats.kruskal(*groups.values())
    lines.append(f"  H = {stat:.3f}, p = {p:.4f}  (n per group: "
                 f"{[len(v) for v in groups.values()]}; SMALL N -- exploratory only)")
    lines.append(f"  Pairwise Wilcoxon rank-sum (uncorrected, exploratory):")
    keys = list(groups.keys())
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            try:
                s2, p2 = stats.ranksums(groups[keys[i]], groups[keys[j]])
                lines.append(f"    {keys[i]} vs {keys[j]}: stat={s2:.3f} p={p2:.4f}")
            except Exception as e:
                lines.append(f"    {keys[i]} vs {keys[j]}: n/a ({e})")
    return lines

def bray_curtis(mat):
    n = mat.shape[0]
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            num = np.abs(mat[i] - mat[j]).sum()
            den = (mat[i] + mat[j]).sum()
            d = num / den if den > 0 else 0.0
            D[i, j] = D[j, i] = d
    return D

def classical_pcoa(D):
    n = D.shape[0]
    D2 = D ** 2
    J = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * J @ D2 @ J
    eigvals, eigvecs = np.linalg.eigh(B)
    idx = np.argsort(eigvals)[::-1]
    eigvals, eigvecs = eigvals[idx], eigvecs[:, idx]
    pos = eigvals > 1e-8
    coords = eigvecs[:, pos] * np.sqrt(eigvals[pos])
    var_explained = eigvals[pos] / eigvals[pos].sum()
    return coords, var_explained

def main():
    report = ["PHASE 7 -- STATISTICS", "=" * 60]
    sample_names, totals = load_candidate_abundance()
    report.append(f"\nSummed High-confidence + Moderate candidate-gene TPM per sample:")
    for s in sample_names:
        c, src, tp = parse_sample_meta(s)
        report.append(f"  {s:8s} ({c:10s} {src:15s} T{tp})  {totals[s]:8.3f}")

    report += kruskal_wallis_and_pairwise(
        sample_names, totals, lambda s: parse_sample_meta(s)[1], "source type")
    report += kruskal_wallis_and_pairwise(
        sample_names, totals, lambda s: parse_sample_meta(s)[0], "city")

    report.append("\nNOTE: with n=2 per city x source cell, ALL comparisons above are "
                   "exploratory/descriptive, not confirmatory, per the plan.")

    # PCoA on OTU table
    with open(OTU, encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
    otu_samples = header[1:-1]
    otu_rows = []
    with open(OTU, encoding="utf-8") as f:
        next(f)
        for line in f:
            parts = line.rstrip("\n").split("\t")
            otu_rows.append([float(x) for x in parts[1:-1]])
    mat = np.array(otu_rows).T  # samples x OTUs
    D = bray_curtis(mat)
    coords, var_exp = classical_pcoa(D)
    np.save(OUT / "pcoa_coords.npy", coords)
    with open(OUT / "pcoa_results.json", "w", encoding="utf-8") as f:
        json.dump(dict(
            samples=otu_samples,
            coords=coords[:, :4].tolist(),
            var_explained_pct=[round(100 * v, 2) for v in var_exp[:6]],
        ), f, indent=2)
    report.append(f"\n--- PCoA (Bray-Curtis, OTU table) ---")
    report.append(f"Variance explained, first 6 axes (%): {[round(100*v,2) for v in var_exp[:6]]}")

    # correlate candidate-gene abundance with community structure axis 1 and with
    # dehalogenase-relevant KO/GO abundance (phase 3 output)
    phase3_path = ROOT / "03_functional_context" / "dehalogenase_relevant_term_abundance_by_sample.tsv"
    if phase3_path.exists():
        with open(phase3_path, encoding="utf-8") as f:
            p3 = {row["sample"]: row for row in csv.DictReader(f, delimiter="\t")}
        common = [s for s in sample_names if s in p3]
        x = np.array([totals[s] for s in common])
        y = np.array([float(p3[s]["KEGG_dehalogenase_KO_sum"]) for s in common])
        if len(common) > 2 and x.std() > 0 and y.std() > 0:
            rho, p = stats.spearmanr(x, y)
            report.append(f"\nSpearman correlation: candidate-gene TPM vs KEGG "
                          f"dehalogenase-relevant KO abundance, n={len(common)}: "
                          f"rho={rho:.3f}, p={p:.4f} (exploratory, n=18)")

    out_txt = OUT / "phase7_report.txt"
    out_txt.write_text("\n".join(report), encoding="utf-8")
    print("\n".join(report))
    print(f"\nWrote {out_txt}")

if __name__ == "__main__":
    main()
