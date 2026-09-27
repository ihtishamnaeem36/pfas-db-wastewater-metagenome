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
Generate Figure 12: Atomic-Resolution 3D Structural Superposition of Priority Dehalogenase Candidates.
Panels:
  (A) Full tertiary fold superposition of FAcD candidate MSLW1_contig_257141_2 (cyan) on solved defluorinase PDB 8SDC (salmon).
  (B) Active-site zoom-in showing spatial alignment of the catalytic triad (Asp104, His271, Asp128).
  (C) Chlorobenzoyl-CoA dehalogenase candidate PHW2_contig_78442_19 (marine) on PDB 1NZY (yellow-orange).
  (D) HAD dehalogenase candidate PHW1_contig_70865_19 (green) on PDB 3UMG (wheat).
"""
import os
import sys
from pathlib import Path
from PIL import Image

# Headless rendering setup
os.environ["PYOPENGL_PLATFORM"] = "osmesa"

import pymol
pymol.pymol_argv = ["pymol", "-cq"]
pymol.finish_launching()

from pymol import cmd

ROOT = Path(".").resolve()
FIG_DIR = ROOT / "figures"
PRED_DIR = ROOT / "06_structural_modeling" / "predictions"
REF_DIR = ROOT / "06_structural_modeling" / "reference_pdbs"

TMP_DIR = ROOT / "figures" / "_tmp_3d"
TMP_DIR.mkdir(exist_ok=True)

# ---------------------------------------------------------------------------
# Panel A: FAcD global superposition (MSLW1 on 8SDC)
# ---------------------------------------------------------------------------
cmd.reinitialize()
cmd.bg_color("white")
cmd.set("ray_opaque_background", 1)
cmd.set("antialias", 2)

cand_a = PRED_DIR / "MSLW1_contig_257141_2.pdb"
ref_a = REF_DIR / "8SDC.pdb"

cmd.load(str(cand_a), "cand_facd")
cmd.load(str(ref_a), "ref_facd")
cmd.align("cand_facd", "ref_facd and chain A")

cmd.hide("everything")
cmd.show("cartoon", "cand_facd or (ref_facd and chain A)")
cmd.color("cyan", "cand_facd")
cmd.color("salmon", "ref_facd and chain A")
cmd.set("cartoon_transparency", 0.20, "ref_facd and chain A")
cmd.orient("cand_facd")
cmd.zoom("cand_facd", buffer=2.0)

png_a = TMP_DIR / "panel_a.png"
cmd.png(str(png_a), width=1200, height=1000, dpi=300, ray=0)

# ---------------------------------------------------------------------------
# Panel B: FAcD Active-Site Zoom-in (Catalytic Triad: Asp104, His271, Asp128)
# ---------------------------------------------------------------------------
cmd.hide("everything")
cmd.show("cartoon", "cand_facd or (ref_facd and chain A)")
cmd.set("cartoon_transparency", 0.70)

# Select triad residues in candidate and reference
cmd.select("triad_cand", "cand_facd and (resi 104 or resi 271 or resi 128)")
cmd.select("triad_ref", "ref_facd and chain A and (resi 104 or resi 271 or resi 128)")

cmd.show("sticks", "triad_cand or triad_ref")
cmd.color("cyan", "triad_cand and not elem N and not elem O")
cmd.color("salmon", "triad_ref and not elem N and not elem O")
cmd.set("stick_radius", 0.28)

cmd.zoom("triad_cand", buffer=3.5)
png_b = TMP_DIR / "panel_b.png"
cmd.png(str(png_b), width=1200, height=1000, dpi=300, ray=0)

# ---------------------------------------------------------------------------
# Panel C: Chlorobenzoyl-CoA dehalogenase (PHW2 on 1NZY)
# ---------------------------------------------------------------------------
cmd.reinitialize()
cmd.bg_color("white")
cand_c = PRED_DIR / "PHW2_contig_78442_19.pdb"
ref_c = REF_DIR / "1NZY.pdb"

cmd.load(str(cand_c), "cand_cba")
cmd.load(str(ref_c), "ref_cba")
cmd.align("cand_cba", "ref_cba and chain A")

cmd.hide("everything")
cmd.show("cartoon", "cand_cba or (ref_cba and chain A)")
cmd.color("marine", "cand_cba")
cmd.color("warmpink", "ref_cba and chain A")
cmd.set("cartoon_transparency", 0.25, "ref_cba and chain A")
cmd.orient("cand_cba")
cmd.zoom("cand_cba", buffer=2.0)

png_c = TMP_DIR / "panel_c.png"
cmd.png(str(png_c), width=1200, height=1000, dpi=300, ray=0)

# ---------------------------------------------------------------------------
# Panel D: HAD dehalogenase (PHW1 on 3UMG)
# ---------------------------------------------------------------------------
cmd.reinitialize()
cmd.bg_color("white")
cand_d = PRED_DIR / "PHW1_contig_70865_19.pdb"
ref_d = REF_DIR / "3UMG.pdb"

cmd.load(str(cand_d), "cand_had")
cmd.load(str(ref_d), "ref_had")
cmd.align("cand_had", "ref_had and chain A")

cmd.hide("everything")
cmd.show("cartoon", "cand_had or (ref_had and chain A)")
cmd.color("forest", "cand_had")
cmd.color("wheat", "ref_had and chain A")
cmd.set("cartoon_transparency", 0.25, "ref_had and chain A")
cmd.orient("cand_had")
cmd.zoom("cand_had", buffer=2.0)

png_d = TMP_DIR / "panel_d.png"
cmd.png(str(png_d), width=1200, height=1000, dpi=300, ray=0)

# Save master session
cmd.save(str(FIG_DIR / "Fig12_3d_superpositions.pse"))
cmd.quit()
print("PyMOL rendering complete.")

# ---------------------------------------------------------------------------
# Composite 4-Panel Figure via Matplotlib
# ---------------------------------------------------------------------------
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

fig = plt.figure(figsize=(7.5, 6.8))
gs = gridspec.GridSpec(2, 2, hspace=0.18, wspace=0.15)

panels = [
    (gs[0, 0], png_a, "A  FAcD Fold Superposition (MSLW1 vs. 8SDC)\n    TM-score = 0.985, RMSD = 0.99 Å"),
    (gs[0, 1], png_b, "B  Active-Site Catalytic Triad Alignment\n    Asp104 - His271 - Asp128 (cyan / salmon)"),
    (gs[1, 0], png_c, "C  Chlorobenzoyl-CoA Dehalogenase (PHW2 vs. 1NZY)\n    TM-score = 0.939, RMSD = 2.04 Å"),
    (gs[1, 1], png_d, "D  HAD Dehalogenase Fold (PHW1 vs. 3UMG)\n    TM-score = 0.946, RMSD = 1.71 Å"),
]

for g, path, title in panels:
    ax = fig.add_subplot(g)
    im = Image.open(path)
    ax.imshow(im)
    ax.axis("off")
    ax.set_title(title, fontsize=8.5, fontweight="bold", loc="left", pad=3)

fig_out_pdf = FIG_DIR / "Fig12_3d_structure_superposition.pdf"
fig_out_png = FIG_DIR / "Fig12_3d_structure_superposition.png"

plt.tight_layout()
plt.savefig(fig_out_pdf, dpi=600)
plt.savefig(fig_out_png, dpi=600)
plt.close()

# Cleanup temporary panel images
for p in [png_a, png_b, png_c, png_d]:
    if p.exists():
        p.unlink()
if TMP_DIR.exists():
    try:
        TMP_DIR.rmdir()
    except:
        pass

print(f"Generated Figure 12 at 600 DPI: {fig_out_pdf} and {fig_out_png}")
