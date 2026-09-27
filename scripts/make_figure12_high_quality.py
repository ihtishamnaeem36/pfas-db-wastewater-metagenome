# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "pymol-open-source-whl",
#     "matplotlib",
#     "numpy",
#     "pillow",
# ]
# ///
"""
Generate Publication-Quality Figure 12: 3D Structural Superposition of Priority Dehalogenase Candidates.
Uses PyMOL ray tracing (ray_trace_mode=1, smooth loops, fancy helices, ambient shadows)
to create crisp, publication-grade structural biology artwork.
Panels:
  (A) FAcD tertiary fold superposition: MSLW1_contig_257141_2 (cyan) vs. solved defluorinase PDB 8SDC (salmon).
  (B) Active-site catalytic triad: Asp104 - His271 - Asp128 with stick representation and labeled residues.
  (C) Chlorobenzoyl-CoA dehalogenase: PHW2_contig_78442_19 (marine) vs. 4-chlorobenzoyl-CoA dehalogenase PDB 1NZY (warm pink).
  (D) HAD dehalogenase: PHW1_contig_70865_19 (forest green) vs. Rhodococcus HAD PDB 3UMG (wheat).
"""
import os
import sys
from pathlib import Path
from PIL import Image

os.environ["PYOPENGL_PLATFORM"] = "osmesa"
import pymol
pymol.pymol_argv = ["pymol", "-cq"]
pymol.finish_launching()
from pymol import cmd

ROOT = Path(".").resolve()
FIG_DIR = ROOT / "figures"
PRED_DIR = ROOT / "06_structural_modeling" / "predictions"
REF_DIR = ROOT / "06_structural_modeling" / "reference_pdbs"
TMP_DIR = ROOT / "figures" / "_tmp_3d_hq"
TMP_DIR.mkdir(exist_ok=True)

def setup_scene():
    cmd.reinitialize()
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 1)
    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("cartoon_smooth_loops", 1)
    cmd.set("ray_shadows", 1)
    cmd.set("depth_cue", 1)
    cmd.set("antialias", 2)
    cmd.set("light_count", 2)
    cmd.set("spec_power", 250)
    cmd.set("ray_trace_mode", 1)  # Publication outline mode

# ---------------------------------------------------------------------------
# Panel A: FAcD Global Fold (MSLW1 vs 8SDC)
# ---------------------------------------------------------------------------
setup_scene()
cmd.load(str(PRED_DIR / "MSLW1_contig_257141_2.pdb"), "cand_facd")
cmd.load(str(REF_DIR / "8SDC.pdb"), "ref_facd")
cmd.align("cand_facd", "ref_facd and chain A")

cmd.hide("everything")
cmd.show("cartoon", "cand_facd or (ref_facd and chain A)")
cmd.color("cyan", "cand_facd")
cmd.color("salmon", "ref_facd and chain A")
cmd.set("cartoon_transparency", 0.25, "ref_facd and chain A")
cmd.orient("cand_facd")
cmd.zoom("cand_facd", buffer=1.5)

png_a = TMP_DIR / "panel_a.png"
cmd.ray(1600, 1300)
cmd.png(str(png_a))

# ---------------------------------------------------------------------------
# Panel B: Active-Site Catalytic Triad (Asp104, His271, Asp128)
# ---------------------------------------------------------------------------
cmd.hide("everything")
cmd.show("cartoon", "cand_facd or (ref_facd and chain A)")
cmd.set("cartoon_transparency", 0.75)

cmd.select("triad_c", "cand_facd and (resi 104 or resi 271 or resi 128)")
cmd.select("triad_r", "ref_facd and chain A and (resi 104 or resi 271 or resi 128)")

cmd.show("sticks", "triad_c or triad_r")
cmd.set("stick_radius", 0.32)

# Color carbons by protein, heteroatoms by element
cmd.color("cyan", "triad_c and elem C")
cmd.color("salmon", "triad_r and elem C")
cmd.color("red", "(triad_c or triad_r) and elem O")
cmd.color("blue", "(triad_c or triad_r) and elem N")

cmd.zoom("triad_c", buffer=3.0)
png_b = TMP_DIR / "panel_b.png"
cmd.ray(1600, 1300)
cmd.png(str(png_b))

# ---------------------------------------------------------------------------
# Panel C: Chlorobenzoyl-CoA dehalogenase (PHW2 vs 1NZY)
# ---------------------------------------------------------------------------
setup_scene()
cmd.load(str(PRED_DIR / "PHW2_contig_78442_19.pdb"), "cand_cba")
cmd.load(str(REF_DIR / "1NZY.pdb"), "ref_cba")
cmd.align("cand_cba", "ref_cba and chain A")

cmd.hide("everything")
cmd.show("cartoon", "cand_cba or (ref_cba and chain A)")
cmd.color("marine", "cand_cba")
cmd.color("warmpink", "ref_cba and chain A")
cmd.set("cartoon_transparency", 0.25, "ref_cba and chain A")
cmd.orient("cand_cba")
cmd.zoom("cand_cba", buffer=1.5)

png_c = TMP_DIR / "panel_c.png"
cmd.ray(1600, 1300)
cmd.png(str(png_c))

# ---------------------------------------------------------------------------
# Panel D: HAD dehalogenase (PHW1 vs 3UMG)
# ---------------------------------------------------------------------------
setup_scene()
cmd.load(str(PRED_DIR / "PHW1_contig_70865_19.pdb"), "cand_had")
cmd.load(str(REF_DIR / "3UMG.pdb"), "ref_had")
cmd.align("cand_had", "ref_had and chain A")

cmd.hide("everything")
cmd.show("cartoon", "cand_had or (ref_had and chain A)")
cmd.color("forest", "cand_had")
cmd.color("wheat", "ref_had and chain A")
cmd.set("cartoon_transparency", 0.25, "ref_had and chain A")
cmd.orient("cand_had")
cmd.zoom("cand_had", buffer=1.5)

png_d = TMP_DIR / "panel_d.png"
cmd.ray(1600, 1300)
cmd.png(str(png_d))

# Save PyMOL session
cmd.save(str(FIG_DIR / "Fig12_3d_superpositions.pse"))
cmd.quit()
print("High-quality ray tracing completed.")

# ---------------------------------------------------------------------------
# Combine Panels with Matplotlib at Publication Layout
# ---------------------------------------------------------------------------
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

fig = plt.figure(figsize=(8.0, 7.0))
gs = gridspec.GridSpec(2, 2, hspace=0.15, wspace=0.12)

panels = [
    (gs[0, 0], png_a, "A  FAcD Fold Superposition (MSLW1 vs. 8SDC)\n    TM-score = 0.985, RMSD = 0.99 Å (cyan / salmon)"),
    (gs[0, 1], png_b, "B  Active-Site Catalytic Triad Alignment\n    Asp104 - His271 - Asp128 (sticks; N blue, O red)"),
    (gs[1, 0], png_c, "C  Chlorobenzoyl-CoA Dehalogenase (PHW2 vs. 1NZY)\n    TM-score = 0.939, RMSD = 2.04 Å (marine / pink)"),
    (gs[1, 1], png_d, "D  HAD Dehalogenase Fold (PHW1 vs. 3UMG)\n    TM-score = 0.946, RMSD = 1.71 Å (green / wheat)"),
]

for g, path, title in panels:
    ax = fig.add_subplot(g)
    im = Image.open(path)
    ax.imshow(im)
    ax.axis("off")
    ax.set_title(title, fontsize=8.8, fontweight="bold", loc="left", pad=4)

fig_out_pdf = FIG_DIR / "Fig12_3d_structure_superposition.pdf"
fig_out_png = FIG_DIR / "Fig12_3d_structure_superposition.png"

plt.tight_layout()
plt.savefig(fig_out_pdf, dpi=600)
plt.savefig(fig_out_png, dpi=600)
plt.close()

# Cleanup temporary panel files
for p in [png_a, png_b, png_c, png_d]:
    if p.exists():
        p.unlink()
if TMP_DIR.exists():
    try:
        TMP_DIR.rmdir()
    except:
        pass

print(f"Generated High-Quality Figure 12 at 600 DPI: {fig_out_pdf} and {fig_out_png}")
