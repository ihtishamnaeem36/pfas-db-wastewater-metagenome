"""
Generate Figure 11: 3D Structural Modeling and Foldseek TM-align Validation.
Panels:
  (A) Per-candidate mean and pocket pLDDT confidence (ESMFold).
  (B) Foldseek TM-scores against experimentally solved reference crystal structures.
  (C) TM-score vs. Sequence Identity (fold conservation despite divergence).
  (D) C-alpha RMSD distribution across verified structural homologs.
"""
import csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

ROOT = Path(__file__).resolve().parents[1]
FIG_DIR = ROOT / "figures"
EVID_TSV = ROOT / "06_structural_modeling" / "structural_evidence_table.tsv"
PLDDT_TSV = ROOT / "06_structural_modeling" / "plddt_summary.tsv"
FOLDSEEK_TSV = ROOT / "06_structural_modeling" / "foldseek_tmalign_results.tsv"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8.5,
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

def make_fig11():
    rows = []
    with open(EVID_TSV, encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["TM_score"] and float(r["TM_score"]) > 0:
                rows.append({
                    "gene_id": r["gene_id"],
                    "ref": r["reference_structure_compared"].split()[0],
                    "mean_plddt": float(r["mean_pLDDT"]),
                    "active_plddt": float(r["active_site_pLDDT"]),
                    "tm": float(r["TM_score"]),
                    "rmsd": float(r["RMSD_angstrom"]),
                    "notes": r["notes"],
                })

    # Parse foldseek hits to get sequence identity
    fident_map = {}
    with open(FOLDSEEK_TSV, encoding="utf-8") as f:
        for line in f:
            p = line.strip().split("\t")
            if len(p) >= 12:
                gid = p[0]
                fid = float(p[6]) * 100.0
                tm = float(p[2])
                if gid not in fident_map or tm > fident_map[gid][0]:
                    fident_map[gid] = (tm, fid)

    fig = plt.figure(figsize=(7.2, 5.8))
    gs = gridspec.GridSpec(2, 2, hspace=0.35, wspace=0.32)

    # Panel A: pLDDT Distribution (Overall vs. Pocket)
    axA = fig.add_subplot(gs[0, 0])
    y_pos = np.arange(len(rows))
    cands = [r["gene_id"].split("_contig_")[1] for r in rows]
    mean_plddt = [r["mean_plddt"] for r in rows]
    active_plddt = [r["active_plddt"] for r in rows]

    axA.barh(y_pos - 0.18, mean_plddt, height=0.35, color="#0072B2", label="Mean pLDDT", alpha=0.85)
    axA.barh(y_pos + 0.18, active_plddt, height=0.35, color="#D55E00", label="Active-site pLDDT", alpha=0.85)
    axA.axvline(70, color="gray", linestyle="--", linewidth=0.8, alpha=0.7)
    axA.axvline(90, color="#009E73", linestyle=":", linewidth=0.8, alpha=0.7)
    axA.set_yticks(y_pos)
    axA.set_yticklabels(cands, fontsize=7.5)
    axA.set_xlabel("pLDDT Score (ESMFold)", fontsize=8.5)
    axA.set_xlim(50, 100)
    axA.set_title("A  Structural Prediction Confidence", fontsize=9, fontweight="bold", loc="left")
    axA.legend(loc="lower left", fontsize=7.5, frameon=False)
    axA.text(71, -0.8, "Conf.", fontsize=7, color="gray")
    axA.text(91, -0.8, "V.High", fontsize=7, color="#009E73")

    # Panel B: TM-scores against Solved Reference PDBs
    axB = fig.add_subplot(gs[0, 1])
    tm_scores = [r["tm"] for r in rows]
    refs = [r["ref"] for r in rows]
    colors = ["#D55E00" if t >= 0.8 else ("#E69F00" if t >= 0.5 else "#999999") for t in tm_scores]

    bars = axB.barh(y_pos, tm_scores, color=colors, height=0.6, alpha=0.85)
    axB.axvline(0.5, color="red", linestyle="--", linewidth=0.9, label="Same fold (TM=0.5)")
    axB.axvline(0.8, color="#009E73", linestyle=":", linewidth=0.9, label="High fidelity (TM=0.8)")
    axB.set_yticks(y_pos)
    axB.set_yticklabels([f"{refs[i]}" for i in range(len(rows))], fontsize=7.5)
    axB.set_xlabel("Foldseek Alignment TM-score", fontsize=8.5)
    axB.set_xlim(0.3, 1.05)
    axB.set_title("B  Tertiary Fold Alignment (TM-align)", fontsize=9, fontweight="bold", loc="left")
    axB.legend(loc="lower left", fontsize=7.5, frameon=False)

    # Panel C: TM-score vs. Sequence Identity
    axC = fig.add_subplot(gs[1, 0])
    xs = [fident_map[r["gene_id"]][1] for r in rows]
    ys = [r["tm"] for r in rows]
    axC.scatter(xs, ys, color="#0072B2", s=45, alpha=0.85, edgecolors="k", linewidth=0.5)
    for i, r in enumerate(rows):
        if r["tm"] > 0.94 or r["gene_id"].startswith("MSLW1_contig_257141"):
            lbl = r["gene_id"].split("_")[0] + "_" + r["gene_id"].split("_")[-1]
            axC.annotate(lbl, (xs[i], ys[i]), textcoords="offset points", xytext=(5, -2), fontsize=6.8)
    axC.axhline(0.5, color="red", linestyle="--", linewidth=0.8, alpha=0.6)
    axC.set_xlabel("Sequence Identity to Reference (%)", fontsize=8.5)
    axC.set_ylabel("Structural TM-Score", fontsize=8.5)
    axC.set_ylim(0.4, 1.02)
    axC.set_xlim(15, 60)
    axC.set_title("C  Structure Conserved Despite Sequence Divergence", fontsize=9, fontweight="bold", loc="left")
    axC.text(17, 0.52, "Same global fold threshold", fontsize=7, color="red")

    # Panel D: RMSD (C-alpha) vs. Reference Length
    axD = fig.add_subplot(gs[1, 1])
    rmsds = [r["rmsd"] for r in rows]
    axD.hist(rmsds, bins=8, range=(0.5, 6.0), color="#009E73", edgecolor="white", alpha=0.85)
    axD.axvline(2.5, color="orange", linestyle="--", linewidth=0.8)
    axD.set_xlabel("Cα RMSD to Solved Reference (Å)", fontsize=8.5)
    axD.set_ylabel("Number of Candidates", fontsize=8.5)
    axD.set_title("D  Structural Superposition Coordinate Precision", fontsize=9, fontweight="bold", loc="left")
    axD.text(2.6, 3.5, "Sub-2.5 Å precision", fontsize=7.5, color="#D55E00")

    pdf_out = FIG_DIR / "Fig11_structural_validation.pdf"
    png_out = FIG_DIR / "Fig11_structural_validation.png"
    plt.tight_layout()
    plt.savefig(pdf_out, dpi=300)
    plt.savefig(png_out, dpi=300)
    plt.close()
    print(f"Saved Figure 11 to {pdf_out} and {png_out}")

if __name__ == "__main__":
    make_fig11()
