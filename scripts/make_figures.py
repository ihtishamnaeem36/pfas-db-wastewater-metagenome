"""
High-impact-journal-style figures for the PFAS candidate-gene discovery pipeline.
All output as vector PDF (infinite resolution). No default-matplotlib bar-chart
styling; each figure is purpose-built for the evidence type it shows.
"""
import csv
import re
import json
from pathlib import Path
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.patches import Circle
import matplotlib.gridspec as gridspec

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "figures"
FIG.mkdir(exist_ok=True)

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 9,
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "pdf.fonttype": 42,   # editable text in Illustrator, not outlines
    "ps.fonttype": 42,
    "svg.fonttype": "none",
})

# Okabe-Ito colorblind-safe qualitative palette (avoids the red/green pairing that
# is indistinguishable under deuteranopia/protanopia -- required for high-IF
# submission accessibility standards).
CITY_COLORS = {"Mardan": "#0072B2", "Peshawar": "#E69F00", "Swat": "#009E73"}
TIER_COLORS = {"High-confidence": "#D55E00", "Moderate": "#E69F00", "Exploratory": "#BBBBBB"}
SOURCE_MARKERS = {"Hospital": "o", "Community": "s", "Slaughterhouse": "^"}

def parse_sample_meta(s):
    m = re.match(r"^([MPS])(HW|CW|SLW)(\d)$", s)
    city = {"M": "Mardan", "P": "Peshawar", "S": "Swat"}[m.group(1)]
    source = {"HW": "Hospital", "CW": "Community", "SLW": "Slaughterhouse"}[m.group(2)]
    return city, source, m.group(3)


# ---------------------------------------------------------------------------
# Figure 1: Study design + Phase-0/1/2 screening funnel (schematic + Sankey-like)
# ---------------------------------------------------------------------------
def fig1_screening_funnel():
    fig = plt.figure(figsize=(7.2, 4.2))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1.3], wspace=0.35)

    # Panel A: 3x3x2 design grid
    axA = fig.add_subplot(gs[0])
    cities = ["Mardan", "Peshawar", "Swat"]
    sources = ["Hospital", "Community", "Slaughterhouse"]
    for i, city in enumerate(cities):
        for j, source in enumerate(sources):
            for t in range(2):
                x = j + (t - 0.5) * 0.28
                y = i
                axA.scatter(x, y, s=170, color=CITY_COLORS[city],
                            marker=SOURCE_MARKERS[source], edgecolor="white",
                            linewidth=0.8, zorder=3)
    axA.set_xticks(range(3)); axA.set_xticklabels(sources, fontsize=8)
    axA.set_yticks(range(3)); axA.set_yticklabels(cities, fontsize=8)
    axA.set_xlim(-0.6, 2.6); axA.set_ylim(-0.6, 2.6)
    axA.set_title("A. Study design (n = 18)\n3 cities x 3 sources x 2 timepoints",
                   fontsize=9, loc="left", fontweight="bold")
    axA.grid(alpha=0.15, zorder=0)
    for spine in axA.spines.values():
        spine.set_visible(False)
    legend_elems = [plt.Line2D([0], [0], marker=SOURCE_MARKERS[s], color="w",
                     markerfacecolor="grey", markersize=7, label=s) for s in sources]
    axA.legend(handles=legend_elems, loc="upper center", bbox_to_anchor=(0.5, -0.12),
               ncol=3, frameon=False, fontsize=7)

    # Panel B: funnel
    axB = fig.add_subplot(gs[1])
    stages = [
        ("Predicted proteins\n(corrected NR catalogue)", 777124),
        ("DIAMOND+HMMER hits\n(71 refs, e<1e-5)", 11907),
        ("High-specificity\nreference hits", 1819),
        ("High-confidence +\nModerate tier", 386),
        ("Priority 1\n(top validation targets)", 166),
    ]
    y = np.arange(len(stages))[::-1]
    widths = np.array([np.log10(v + 1) for _, v in stages])
    widths = widths / widths.max()
    colors = plt.cm.Blues(np.linspace(0.35, 0.9, len(stages)))
    for yi, (label, val), w, c in zip(y, stages, widths, colors):
        axB.barh(yi, w, height=0.6, color=c, edgecolor="black", linewidth=0.6)
        axB.text(w + 0.02, yi, f"{label}  (n={val:,})", va="center", ha="left", fontsize=7.5)
    axB.set_xlim(0, 2.6)
    axB.set_ylim(-0.7, len(stages) - 0.3)
    axB.axis("off")
    axB.set_title("B. Curated-database screening funnel", fontsize=9, loc="left", fontweight="bold")

    fig.savefig(FIG / "Fig1_study_design_and_screening_funnel.pdf", bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 2: Homology evidence -- identity vs coverage scatter, sized by TPM,
# colored by confidence tier, reference-tier shape
# ---------------------------------------------------------------------------
def fig2_homology_evidence():
    with open(ROOT / "02_homology_screening" / "phase2_ranked_candidates.tsv", encoding="utf-8") as f:
        all_rows = list(csv.DictReader(f, delimiter="\t"))
    # HMM-only hits (best_ref = "HMM:family_*.aln") have no DIAMOND alignment, so no
    # identity/coverage-of-reference-length to plot on these axes; excluded here and
    # noted in the caption/text (they are a small, real, separate class of evidence --
    # not omitted from the candidate table, just not plottable on this particular axis).
    rows = [r for r in all_rows if not r["best_ref"].startswith("HMM:")]
    n_hmm_only = len(all_rows) - len(rows)

    fig, ax = plt.subplots(figsize=(5.8, 4.8))
    tier_colors = TIER_COLORS
    tier_order = ["Exploratory", "Moderate", "High-confidence"]
    for tier in tier_order:
        sub = [r for r in rows if r["confidence_tier"] == tier]
        if not sub:
            continue
        x = [float(r["coverage_pct"]) for r in sub]
        y = [float(r["pident"]) for r in sub]
        s = [20 + 20 * np.log2(1 + float(r["total_tpm_across_18_samples"])) for r in sub]
        marker = ["v" if r["ref_tier"] == "B" else "o" for r in sub]
        alpha = 0.55 if tier == "Exploratory" else 0.85
        for xi, yi, si, mi in zip(x, y, s, marker):
            ax.scatter(xi, yi, s=si, marker=mi, color=tier_colors[tier],
                       edgecolor="black", linewidth=0.3, alpha=alpha, zorder=3,
                       label=tier if (xi, yi) == (x[0], y[0]) else None)

    ax.axvspan(70, 100, color="#D55E00", alpha=0.05, zorder=0)
    ax.axhspan(30, 100, color="#D55E00", alpha=0.05, zorder=0)
    ax.axvline(70, color="grey", lw=0.6, ls="--", zorder=1)
    ax.axhline(30, color="grey", lw=0.6, ls="--", zorder=1)
    ax.text(71, 3, "coverage >= 70%", fontsize=6.5, color="grey", rotation=90, va="bottom")
    ax.text(2, 31, "identity >= 30%", fontsize=6.5, color="grey")

    ax.set_xlabel("Alignment coverage of reference length (%)")
    ax.set_ylabel("Percent identity to best-matching curated reference (%)")
    ax.set_title("Sequence-homology evidence for PFAS candidate-gene screening hits",
                  fontsize=9.5, loc="left", fontweight="bold")

    from matplotlib.lines import Line2D
    handles = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=tier_colors[t],
               markeredgecolor="black", markersize=7, label=t) for t in tier_order
    ] + [
        Line2D([0], [0], marker="o", color="w", markerfacecolor="grey",
               markeredgecolor="black", markersize=7, label="Tier A reference"),
        Line2D([0], [0], marker="v", color="w", markerfacecolor="grey",
               markeredgecolor="black", markersize=7, label="Tier B reference"),
    ]
    ax.legend(handles=handles, loc="lower right", fontsize=6.5, frameon=False)
    ax.text(0.02, 0.98, f"Marker size ~ log2(summed TPM across 18 samples)\n"
            f"n={len(rows)} DIAMOND-aligned hits shown ({n_hmm_only} additional "
            f"HMM-only hits not plottable on these axes)",
            transform=ax.transAxes, fontsize=6, va="top", color="grey")

    fig.savefig(FIG / "Fig2_homology_evidence_scatter.pdf", bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 3: Candidate abundance heatmap across 18 samples (High-conf + Moderate)
# ---------------------------------------------------------------------------
def fig3_abundance_heatmap(top_n=20):
    with open(ROOT / "02_homology_screening" / "phase2_ranked_candidates.tsv", encoding="utf-8") as f:
        rows = [r for r in csv.DictReader(f, delimiter="\t")
                if r["confidence_tier"] in ("High-confidence", "Moderate")]
    n_total = len(rows)
    tpm_cols = [c for c in rows[0].keys() if c.startswith("TPM_")]
    samples = [c[4:] for c in tpm_cols]
    order = sorted(range(len(samples)), key=lambda i: (parse_sample_meta(samples[i])[0],
                                                          parse_sample_meta(samples[i])[1],
                                                          samples[i]))
    samples = [samples[i] for i in order]
    tpm_cols = [tpm_cols[i] for i in order]

    rows.sort(key=lambda r: (0 if r["confidence_tier"] == "High-confidence" else 1, -float(r["sw_score"])))
    rows = rows[:top_n]
    mat = np.array([[np.log1p(float(r[c])) for c in tpm_cols] for r in rows])
    labels = [f"{r['gene_id']}  ({r['confidence_tier'][:4]}, {r['ref_tier']})" for r in rows]

    fig, ax = plt.subplots(figsize=(8.2, 4.2))
    im = ax.imshow(mat, aspect="auto", cmap="YlOrRd", vmin=0)
    ax.set_xticks(range(len(samples)))
    ax.set_xticklabels(samples, rotation=90, fontsize=6.5)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=6.5)

    # city color strip above
    city_bounds = defaultdict(list)
    for i, s in enumerate(samples):
        city_bounds[parse_sample_meta(s)[0]].append(i)
    for city, idxs in city_bounds.items():
        ax.axvspan(min(idxs) - 0.5, max(idxs) + 0.5, ymin=1.0, ymax=1.03,
                   color=CITY_COLORS[city], clip_on=False, transform=ax.get_xaxis_transform())

    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label("log1p(TPM)", fontsize=7.5)
    cbar.ax.tick_params(labelsize=6.5)

    ax.set_title(f"Candidate-gene abundance across 18 wastewater metagenomes\n"
                 f"(top {len(rows)} of {n_total} High-confidence/Moderate candidates, by SW score)",
                  fontsize=9, loc="left", fontweight="bold")
    from matplotlib.patches import Patch
    city_handles = [Patch(facecolor=CITY_COLORS[c], edgecolor="none", label=c) for c in CITY_COLORS]
    ax.legend(handles=city_handles, loc="lower left", bbox_to_anchor=(0.0, 1.13),
              ncol=3, frameon=False, fontsize=7, handlelength=1.0, columnspacing=1.0)

    fig.savefig(FIG / "Fig3_candidate_abundance_heatmap.pdf", bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 4: PCoA ordination (community structure) with candidate-abundance overlay
# ---------------------------------------------------------------------------
def fig4_pcoa():
    with open(ROOT / "07_statistics" / "pcoa_results.json", encoding="utf-8") as f:
        pc = json.load(f)
    samples = pc["samples"]
    coords = np.array(pc["coords"])
    var_exp = pc["var_explained_pct"]

    with open(ROOT / "02_homology_screening" / "phase2_ranked_candidates.tsv", encoding="utf-8") as f:
        rows = [r for r in csv.DictReader(f, delimiter="\t")
                if r["confidence_tier"] in ("High-confidence", "Moderate")]
    tpm_cols = [c for c in rows[0].keys() if c.startswith("TPM_")]
    totals = defaultdict(float)
    for r in rows:
        for c in tpm_cols:
            totals[c[4:]] += float(r[c])

    fig, ax = plt.subplots(figsize=(5.4, 4.6))
    for s, (x, y) in zip(samples, coords[:, :2]):
        city, source, tp = parse_sample_meta(s)
        size = 60 + 400 * (totals[s] / (max(totals.values()) + 1e-9))
        ax.scatter(x, y, s=size, color=CITY_COLORS[city], marker=SOURCE_MARKERS[source],
                   edgecolor="black", linewidth=0.5, alpha=0.85, zorder=3)
        ax.annotate(s, (x, y), fontsize=5.5, xytext=(3, 3), textcoords="offset points",
                    color="#333333")

    ax.axhline(0, color="grey", lw=0.4, zorder=0)
    ax.axvline(0, color="grey", lw=0.4, zorder=0)
    ax.set_xlabel(f"PCoA axis 1 ({var_exp[0]:.1f}% var.)")
    ax.set_ylabel(f"PCoA axis 2 ({var_exp[1]:.1f}% var.)")
    ax.set_title("Community structure (Bray-Curtis PCoA) with candidate-gene\n"
                 "abundance overlay (marker size)", fontsize=9, loc="left", fontweight="bold")

    from matplotlib.lines import Line2D
    handles = [Line2D([0], [0], marker="o", color="w", markerfacecolor=c,
               markeredgecolor="black", markersize=7, label=city)
               for city, c in CITY_COLORS.items()]
    handles += [Line2D([0], [0], marker=m, color="w", markerfacecolor="grey",
                markeredgecolor="black", markersize=7, label=s)
                for s, m in SOURCE_MARKERS.items()]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.02, 1.0),
              fontsize=6.5, frameon=False)

    fig.savefig(FIG / "Fig4_PCoA_community_structure.pdf", bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 5: Multi-evidence consolidation matrix (dot-matrix / "oncoprint" style --
# the standard high-IF-journal format for multi-criteria evidence summaries;
# replaces an earlier radar-chart draft, which is discouraged in the data-
# visualization literature because area/angle encodings distort magnitude
# perception and radar panels are rare in Nature/Cell/Elsevier environmental
# journals).
# ---------------------------------------------------------------------------
def fig5_evidence_consolidation(top_n=20):
    with open(ROOT / "08_candidate_prioritization" / "final_candidate_ranking.tsv", encoding="utf-8") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    n_total = len(rows)
    rows.sort(key=lambda r: (-int(r["n_evidence_lines_of_5"]), -float(r["composite_score"])))
    rows = rows[:top_n]

    criteria = ["Sequence\nhomology", "Tier A\nreference", "OTU-supported\ntaxonomy",
                "Genomic\nneighborhood", "Recurrence\n(≥3/18)"]

    def evidence_vector(r):
        homology = 1
        tier_a = 1 if r["reference_tier"] == "A" else 0
        taxo = 1 if r["taxonomy_OTU_supported"] == "True" else 0
        neigh = 1 if (r["contig_relevant_neighbors"] not in ("NA", "0", "")) else 0
        recur = 1 if int(r["n_samples_detected_of_18"]) >= 3 else 0
        return [homology, tier_a, taxo, neigh, recur]

    n = len(rows)
    mat = np.array([evidence_vector(r) for r in rows])
    labels = [r["gene_id"] for r in rows]
    scores = [float(r["composite_score"]) for r in rows]
    tiers = [r["confidence_tier_phase2"] for r in rows]

    fig = plt.figure(figsize=(8.4, 0.34 * n + 1.6))
    gs = gridspec.GridSpec(1, 2, width_ratios=[len(criteria), 1.15], wspace=0.08,
                            left=0.30, right=0.97, top=0.80, bottom=0.22)
    axM = fig.add_subplot(gs[0])
    axS = fig.add_subplot(gs[1])

    for i in range(n):
        for j in range(len(criteria)):
            filled = mat[i, j] == 1
            axM.scatter(j, i, s=170, marker="o",
                        facecolor=TIER_COLORS[tiers[i]] if filled else "white",
                        edgecolor=TIER_COLORS[tiers[i]] if filled else "#BBBBBB",
                        linewidth=1.1, zorder=3)
    axM.set_xlim(-0.6, len(criteria) - 0.4)
    axM.set_ylim(-0.7, n - 0.3)
    axM.set_xticks(range(len(criteria)))
    axM.set_xticklabels(criteria, fontsize=7, rotation=0)
    axM.set_yticks(range(n))
    axM.set_yticklabels(labels, fontsize=7)
    axM.invert_yaxis()
    for spine in axM.spines.values():
        spine.set_visible(False)
    axM.tick_params(length=0)
    axM.set_axisbelow(True)
    axM.grid(axis="y", color="#EEEEEE", lw=0.6, zorder=0)
    axM.set_title("A. Evidence matrix", fontsize=9, loc="left", fontweight="bold")

    axS.barh(range(n), scores, height=0.55,
             color=[TIER_COLORS[t] for t in tiers], edgecolor="black", linewidth=0.4, zorder=3)
    axS.set_ylim(axM.get_ylim())
    axS.set_yticks([])
    axS.set_xlabel("Composite\nscore", fontsize=7)
    axS.tick_params(axis="x", labelsize=6.5)
    for spine in ["top", "right", "left"]:
        axS.spines[spine].set_visible(False)
    axS.set_title("B. Composite\nscore", fontsize=9, loc="left", fontweight="bold")

    from matplotlib.lines import Line2D
    handles = [Line2D([0], [0], marker="o", color="w", markerfacecolor=TIER_COLORS[t],
               markeredgecolor=TIER_COLORS[t], markersize=8, label=t)
               for t in ["High-confidence", "Moderate"]]
    handles.append(Line2D([0], [0], marker="o", color="w", markerfacecolor="white",
                    markeredgecolor="#BBBBBB", markersize=8, label="Criterion not met"))
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.55, 0.0),
               ncol=3, frameon=False, fontsize=7)

    fig.suptitle(f"Multi-evidence support for ranked PFAS candidate genes "
                 f"(top {len(rows)} of {n_total})",
                 fontsize=10.5, fontweight="bold", x=0.30, ha="left", y=0.97)
    fig.savefig(FIG / "Fig5_evidence_consolidation_matrix.pdf", bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 6: Dehalogenase-relevant KEGG/EggNOG/GO community-level context
# ---------------------------------------------------------------------------
def fig6_functional_context():
    path = ROOT / "03_functional_context" / "dehalogenase_relevant_term_abundance_by_sample.tsv"
    with open(path, encoding="utf-8") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))

    fig, axes = plt.subplots(1, 2, figsize=(8.5, 3.6))
    groups = defaultdict(list)
    for r in rows:
        groups[(r["city"], r["source"])].append(float(r["KEGG_dehalogenase_KO_sum"]))

    cities = ["Mardan", "Peshawar", "Swat"]
    sources = ["Hospital", "Community", "Slaughterhouse"]
    ax = axes[0]
    x = np.arange(len(sources))
    width = 0.25
    for i, city in enumerate(cities):
        means = [np.mean(groups.get((city, s), [0])) for s in sources]
        sems = [np.std(groups.get((city, s), [0])) / np.sqrt(max(1, len(groups.get((city, s), [1]))))
                for s in sources]
        ax.bar(x + (i - 1) * width, means, width, yerr=sems, capsize=2,
               color=CITY_COLORS[city], label=city, edgecolor="black", linewidth=0.4)
    ax.set_xticks(x); ax.set_xticklabels(sources, fontsize=8)
    ax.set_ylabel("Summed KEGG dehalogenase-relevant\nKO abundance (TPM-like units)", fontsize=7.5)
    ax.set_title("A. Community-level dehalogenase-KO\nbaseline (n=8 KOs; no PFAS-specific KO exists)",
                 fontsize=8, loc="left", fontweight="bold")
    ax.legend(fontsize=6.5, frameon=False)

    ax2 = axes[1]
    cand_path = ROOT / "02_homology_screening" / "phase2_ranked_candidates.tsv"
    with open(cand_path, encoding="utf-8") as f:
        crows = [r for r in csv.DictReader(f, delimiter="\t")
                 if r["confidence_tier"] in ("High-confidence", "Moderate")]
    tpm_cols = [c for c in crows[0].keys() if c.startswith("TPM_")]
    cand_totals = defaultdict(float)
    for r in crows:
        for c in tpm_cols:
            cand_totals[c[4:]] += float(r[c])
    kegg_by_sample = {r["sample"]: float(r["KEGG_dehalogenase_KO_sum"]) for r in rows}
    common = [s for s in cand_totals if s in kegg_by_sample]
    xs = [cand_totals[s] for s in common]
    ys = [kegg_by_sample[s] for s in common]
    colors = [CITY_COLORS[parse_sample_meta(s)[0]] for s in common]
    markers_ = [SOURCE_MARKERS[parse_sample_meta(s)[1]] for s in common]
    for xi, yi, ci, mi in zip(xs, ys, colors, markers_):
        ax2.scatter(xi, yi, color=ci, marker=mi, s=55, edgecolor="black", linewidth=0.4, zorder=3)
    ax2.set_xlabel("Curated-screen candidate-gene TPM (sum)", fontsize=7.5)
    ax2.set_ylabel("KEGG dehalogenase-relevant KO TPM (sum)", fontsize=7.5)
    ax2.set_title("B. Curated candidates vs. generic\nKEGG dehalogenase-KO abundance",
                 fontsize=8, loc="left", fontweight="bold")

    fig.tight_layout()
    fig.savefig(FIG / "Fig6_functional_context.pdf", bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 7: Permutation-null validation of the Phase 2 screening thresholds
# ---------------------------------------------------------------------------
def fig7_null_model():
    path = ROOT / "02_homology_screening" / "phase2_null_model_real_report.txt"
    text = path.read_text(encoding="utf-8")
    real = {}
    decoy = {}
    for line in text.splitlines():
        parts = line.split()
        if len(parts) >= 4 and parts[0] in ("High-confidence", "Moderate", "Exploratory"):
            real[parts[0]] = int(parts[1])
            decoy[parts[0]] = int(parts[2])

    tiers = ["High-confidence", "Moderate", "Exploratory"]
    fig, ax = plt.subplots(figsize=(4.6, 3.6))
    x = np.arange(len(tiers))
    w = 0.32
    real_vals = [max(real[t], 0.5) for t in tiers]  # floor for log scale display only
    decoy_vals = [max(decoy[t], 0.5) for t in tiers]
    ax.bar(x - w/2, real_vals, w, color=[TIER_COLORS[t] for t in tiers],
           edgecolor="black", linewidth=0.6, label="Real DIAMOND+HMMER screen")
    ax.bar(x + w/2, decoy_vals, w, color="white",
           edgecolor=[TIER_COLORS[t] for t in tiers], hatch="////", linewidth=0.9,
           label="Residue-shuffled\ndecoy catalogue (n=777,124)")
    ax.set_yscale("log")
    for xi, t in zip(x, tiers):
        fdr = decoy[t] / real[t] if real[t] else float("nan")
        ax.text(xi, max(real_vals[list(tiers).index(t)], decoy_vals[list(tiers).index(t)]) * 1.3,
                f"FDR={fdr:.3f}\n({decoy[t]}/{real[t]})", ha="center", fontsize=6.5)

    ax.set_xticks(x)
    ax.set_xticklabels(tiers, fontsize=7.5)
    ax.set_ylabel("Candidate genes reaching tier (log scale)")
    ax.set_ylim(0.3, 3000)
    ax.set_title("Whole-catalogue permutation-null validation\nof the real DIAMOND+HMMER screening thresholds",
                 fontsize=8.7, loc="left", fontweight="bold")
    ax.legend(fontsize=6.3, frameon=False, loc="upper left")
    fig.savefig(FIG / "Fig7_permutation_null_validation.pdf", bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 8: Sequence-based active-site (catalytic triad) conservation matrix
# ---------------------------------------------------------------------------
def fig8_active_site_conservation(top_n=15):
    path = ROOT / "06_structural_modeling" / "phase6a_active_site_seqcheck.tsv"
    if not path.exists():
        return
    with open(path, encoding="utf-8") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    n_total = len(rows)
    rows.sort(key=lambda r: -int(r["n_conservative_or_better_of_3"]))
    rows = rows[:top_n]

    residues = [("Asp104", "Asp104_candidate_residue", "Asp104_class", "D"),
                ("His271", "His271_candidate_residue", "His271_class", "H"),
                ("Asp128", "Asp128_candidate_residue", "Asp128_class", "D")]

    def color_for(cls):
        if cls == "IDENTICAL":
            return "#009E73"
        if cls.startswith("conservative"):
            return "#E69F00"
        if cls.startswith("NON-CONSERVATIVE"):
            return "#D55E00"
        return "#EEEEEE"

    n = len(rows)
    fig, ax = plt.subplots(figsize=(4.6, 0.5 * n + 1.2))
    for i, r in enumerate(rows):
        for j, (label, rescol, clscol, exp) in enumerate(residues):
            c = color_for(r[clscol])
            ax.add_patch(plt.Rectangle((j - 0.45, i - 0.4), 0.9, 0.8, facecolor=c,
                                        edgecolor="black", linewidth=0.5, zorder=2))
            obs = r[rescol] if r[rescol] != "" else "-"
            ax.text(j, i, obs, ha="center", va="center", fontsize=8, fontweight="bold",
                    color="white" if c != "#EEEEEE" else "black", zorder=3)

    ax.set_xlim(-0.7, len(residues) - 0.3)
    ax.set_ylim(-0.7, n - 0.3)
    ax.set_xticks(range(len(residues)))
    ax.set_xticklabels([f"{lbl}\n(ref: DehH1/Q1JU72)" for lbl, *_ in residues], fontsize=7)
    ax.set_yticks(range(n))
    ax.set_yticklabels([r["gene_id"] for r in rows], fontsize=7)
    ax.invert_yaxis()
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(length=0)
    ax.set_title(f"Sequence-based catalytic-triad conservation (top {len(rows)} of {n_total})\n"
                 f"(1-D alignment proxy; not a substitute for Phase 6 structural check)",
                 fontsize=8.3, loc="left", fontweight="bold")

    from matplotlib.patches import Patch
    handles = [Patch(facecolor="#009E73", edgecolor="black", label="Identical"),
               Patch(facecolor="#E69F00", edgecolor="black", label="Conservative substitution"),
               Patch(facecolor="#D55E00", edgecolor="black", label="Non-conservative mismatch"),
               Patch(facecolor="#EEEEEE", edgecolor="black", label="Not aligned (gap)")]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.18),
               ncol=2, frameon=False, fontsize=6.5)

    fig.savefig(FIG / "Fig8_active_site_conservation.pdf", bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 9: Computational biophysical property screening (candidates vs. curated
# reference calibration baseline)
# ---------------------------------------------------------------------------
def fig9_biophysical_properties():
    cand_path = ROOT / "09_biophysical_screening" / "candidate_biophysical_properties.tsv"
    ref_path = ROOT / "09_biophysical_screening" / "reference_biophysical_properties.tsv"
    if not cand_path.exists():
        return
    with open(cand_path, encoding="utf-8") as f:
        cand = list(csv.DictReader(f, delimiter="\t"))
    with open(ref_path, encoding="utf-8") as f:
        ref = list(csv.DictReader(f, delimiter="\t"))

    cand.sort(key=lambda r: -float(r["expressibility_proxy"]))
    n_total = len(cand)
    cand_top = cand[:12] + cand[-8:] if n_total > 20 else cand  # best + worst extremes

    fig, axes = plt.subplots(1, 3, figsize=(9.5, 4.6))

    ax = axes[0]
    y = np.arange(len(cand_top))
    colors = ["#009E73" if r["instability_class"] == "stable" else "#D55E00" for r in cand_top]
    ax.barh(y, [float(r["expressibility_proxy"]) for r in cand_top], color=colors,
            edgecolor="black", linewidth=0.4)
    ax.axvline(np.mean([float(r["expressibility_proxy"]) for r in ref]), color="grey",
               ls="--", lw=1, label="curated-reference mean")
    ax.set_yticks(y)
    ax.set_yticklabels([r["id"] for r in cand_top], fontsize=6)
    ax.invert_yaxis()
    ax.set_xlabel("Expressibility proxy score", fontsize=8)
    title_a = "A. Composite expressibility proxy" if n_total <= 20 else \
        f"A. Expressibility proxy (best 12 + worst 8 of {n_total})"
    ax.set_title(title_a, fontsize=8.7, loc="left", fontweight="bold")
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(facecolor="#009E73", edgecolor="black", label="stable (instab. <40)"),
                        Patch(facecolor="#D55E00", edgecolor="black", label="likely unstable (>=40)")],
               fontsize=6, frameon=False, loc="lower right")

    ax = axes[1]
    ax.scatter([float(r["gravy"]) for r in cand], [float(r["instability_index"]) for r in cand],
               s=45, color="#0072B2", edgecolor="black", linewidth=0.4, label="Candidates", zorder=3)
    ax.scatter([float(r["gravy"]) for r in ref], [float(r["instability_index"]) for r in ref],
               s=45, color="#BBBBBB", edgecolor="black", linewidth=0.4, marker="D",
               label="Curated references", zorder=3)
    ax.axhline(40, color="#D55E00", ls="--", lw=1)
    ax.text(ax.get_xlim()[0] if False else -0.4, 41, "instability threshold", fontsize=6, color="#D55E00")
    ax.set_xlabel("GRAVY (hydropathy)", fontsize=8)
    ax.set_ylabel("Instability index", fontsize=8)
    ax.set_title("B. Hydropathy vs. instability", fontsize=9, loc="left", fontweight="bold")
    ax.legend(fontsize=6.5, frameon=False, loc="upper left")

    ax = axes[2]
    ax.scatter([float(r["aliphatic_index"]) for r in cand], [float(r["mw_kDa"]) for r in cand],
               s=45, color="#0072B2", edgecolor="black", linewidth=0.4, label="Candidates", zorder=3)
    ax.scatter([float(r["aliphatic_index"]) for r in ref], [float(r["mw_kDa"]) for r in ref],
               s=45, color="#BBBBBB", edgecolor="black", linewidth=0.4, marker="D",
               label="Curated references", zorder=3)
    ax.set_xlabel("Aliphatic index (thermostability proxy)", fontsize=8)
    ax.set_ylabel("Molecular weight (kDa)", fontsize=8)
    ax.set_title("C. Aliphatic index vs. size", fontsize=9, loc="left", fontweight="bold")
    ax.legend(fontsize=6.5, frameon=False, loc="upper right")

    fig.suptitle("Computational biophysical property screening (sequence-only; Tier 1)",
                 fontsize=10.5, fontweight="bold", x=0.02, ha="left")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(FIG / "Fig9_biophysical_properties.pdf", bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 10: Family-level summary of the expanded 386-candidate screen
# ---------------------------------------------------------------------------
def fig10_family_summary():
    path = ROOT / "08_candidate_prioritization" / "family_level_summary.tsv"
    if not path.exists():
        return
    with open(path, encoding="utf-8") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    rows = [r for r in rows if int(r["n_candidate_genes"]) >= 2]
    rows.sort(key=lambda r: -int(r["n_candidate_genes"]))

    labels = [r["reference_family"].replace("OriginalCurated:", "Curated: ")
              .replace("Toolkit:", "Toolkit: ") for r in rows]
    n_genes = [int(r["n_candidate_genes"]) for r in rows]
    n_hc = [int(r["n_high_confidence"]) for r in rows]
    n_mod = [int(r["n_moderate"]) for r in rows]
    n_samples = [int(r["n_samples_with_detection"]) for r in rows]
    n_cities = [int(r["n_cities"]) for r in rows]

    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.4), gridspec_kw={"width_ratios": [1.4, 1]})

    ax = axes[0]
    y = np.arange(len(rows))
    ax.barh(y, n_hc, color=TIER_COLORS["High-confidence"], edgecolor="black", linewidth=0.4,
            label="High-confidence")
    ax.barh(y, n_mod, left=n_hc, color=TIER_COLORS["Moderate"], edgecolor="black", linewidth=0.4,
            label="Moderate")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel("Candidate genes")
    ax.set_title("A. Candidates per reference family", fontsize=9, loc="left", fontweight="bold")
    ax.legend(fontsize=7, frameon=False, loc="lower right")

    ax = axes[1]
    colors = ["#0072B2" if c == 3 else ("#56B4E9" if c == 2 else "#BBBBBB") for c in n_cities]
    ax.barh(y, n_samples, color=colors, edgecolor="black", linewidth=0.4)
    ax.set_yticks(y)
    ax.set_yticklabels([])
    ax.invert_yaxis()
    ax.axvline(18, color="grey", ls="--", lw=0.8)
    ax.text(18, -0.7, "all 18\nsamples", fontsize=6, ha="center", color="grey")
    ax.set_xlabel("Samples with detection (/18)")
    ax.set_title("B. Detection breadth\n(color = cities: 3 dark blue, 2 light blue, 1 grey)",
                 fontsize=8.3, loc="left", fontweight="bold")

    fig.suptitle("Family-level summary: 386 candidates collapse to a small number of "
                 "reference enzyme families", fontsize=10, fontweight="bold", x=0.02, ha="left")
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(FIG / "Fig10_family_level_summary.pdf", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    fig1_screening_funnel()
    print("Fig1 done")
    fig2_homology_evidence()
    print("Fig2 done")
    fig3_abundance_heatmap()
    print("Fig3 done")
    fig4_pcoa()
    print("Fig4 done")
    fig5_evidence_consolidation()
    print("Fig5 done")
    fig6_functional_context()
    print("Fig6 done")
    fig7_null_model()
    print("Fig7 done")
    fig8_active_site_conservation()
    print("Fig8 done")
    fig9_biophysical_properties()
    print("Fig9 done")
    fig10_family_summary()
    print("Fig10 done")
    print(f"\nAll figures written to {FIG}")
