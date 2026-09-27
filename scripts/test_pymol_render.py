# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "pymol-open-source-whl",
# ]
# ///
import os
import sys
from pathlib import Path

# Headless platform
os.environ["PYOPENGL_PLATFORM"] = "osmesa"

import pymol
pymol.pymol_argv = ["pymol", "-cq"]
pymol.finish_launching()

from pymol import cmd

cand_pdb = Path("06_structural_modeling/predictions/MSLW1_contig_257141_2.pdb").resolve()
ref_pdb = Path("06_structural_modeling/reference_pdbs/8SDC.pdb").resolve()

print(f"Loading candidate: {cand_pdb}")
cmd.load(str(cand_pdb), "cand")
cand_atoms = cmd.count_atoms("cand")
print(f"Candidate atoms: {cand_atoms}")

print(f"Loading reference: {ref_pdb}")
cmd.load(str(ref_pdb), "ref")
ref_atoms = cmd.count_atoms("ref")
print(f"Reference atoms: {ref_atoms}")

if cand_atoms == 0 or ref_atoms == 0:
    print("Error: empty structure")
    cmd.quit()
    sys.exit(1)

# Align candidate to reference
cmd.align("cand", "ref and chain A")

# Styling
cmd.bg_color("white")
cmd.hide("everything")
cmd.show("cartoon", "cand or (ref and chain A)")
cmd.color("cyan", "cand")
cmd.color("salmon", "ref and chain A")
cmd.set("cartoon_transparency", 0.15, "ref and chain A")

# Orient and zoom
cmd.orient("cand")
cmd.zoom("cand", buffer=2.0)

# Render test png
out_png = Path("figures/test_render_3d.png").resolve()
cmd.png(str(out_png), width=1200, height=900, dpi=300, ray=0)

# Save session
out_pse = Path("figures/test_render_3d.pse").resolve()
cmd.save(str(out_pse))

print(f"Saved test render to {out_png}")
cmd.quit()
