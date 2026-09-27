# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "pymol-open-source-whl",
# ]
# ///
import os
import sys
from pathlib import Path

os.environ["PYOPENGL_PLATFORM"] = "osmesa"
import pymol
pymol.pymol_argv = ["pymol", "-cq"]
pymol.finish_launching()
from pymol import cmd

cmd.reinitialize()
cmd.bg_color("white")

cand = Path("06_structural_modeling/predictions/MSLW1_contig_257141_2.pdb").resolve()
ref = Path("06_structural_modeling/reference_pdbs/8SDC.pdb").resolve()

cmd.load(str(cand), "cand")
cmd.load(str(ref), "ref")
cmd.align("cand", "ref and chain A")

cmd.hide("everything")
cmd.show("cartoon", "cand or (ref and chain A)")

# High quality styling
cmd.color("cyan", "cand")
cmd.color("salmon", "ref and chain A")
cmd.set("cartoon_transparency", 0.25, "ref and chain A")
cmd.set("cartoon_fancy_helices", 1)
cmd.set("cartoon_smooth_loops", 1)
cmd.set("ray_shadows", 1)
cmd.set("depth_cue", 1)
cmd.set("antialias", 2)
cmd.set("light_count", 2)
cmd.set("spec_power", 300)
cmd.set("ray_trace_mode", 1) # Good outline for publication

cmd.orient("cand")
cmd.zoom("cand", buffer=1.5)

cmd.ray(2400, 2000)
cmd.png("figures/test_beauty_ray.png")
print("Beauty ray complete!")
cmd.quit()
