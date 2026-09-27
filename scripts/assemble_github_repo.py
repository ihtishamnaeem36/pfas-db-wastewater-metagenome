"""
Assemble a clean, publication-grade GitHub repository folder:
'pfas-wastewater-metagenome-db'
Includes:
- README.md (comprehensive documentation, database provenance, pipeline usage, citations)
- LICENSE (MIT)
- .gitignore
- requirements.txt & pyproject.toml
- database/ (71 curated reference FASTA, 8 Tier A/B FASTA, HMM profiles, metadata TSVs)
- candidates/ (386 ranked candidate genes, confidence tiers, biophysical metrics, structural evidence)
- structures/ (ESMFold atomic coordinate PDB models and reference crystal structures)
- figures/ (All 12 publication figures in vector PDF and 600 DPI PNG)
- scripts/ (Clean, numbered pipeline scripts for reproduction)
"""
import shutil
from pathlib import Path

SRC = Path(__file__).resolve().parents[1]
DEST = SRC.parent / "pfas-wastewater-metagenome-db"

print(f"Creating GitHub repository folder at:\n  {DEST}")
DEST.mkdir(exist_ok=True)

# 1. Subdirectories
dirs = ["database", "candidates", "structures", "figures", "scripts", "data"]
for d in dirs:
    (DEST / d).mkdir(exist_ok=True)

# 2. Copy Database Files
db_src = SRC / "01_reference_database"
for f in ["curated_pfas_references_71.fasta", "curated_pfas_tierA_tierB_8seqs.fasta",
          "family_FAcD.hmm", "family_HAD.hmm", "reference_metadata_71.tsv",
          "reference_metadata_8seqs.tsv", "mechanistic_specificity_classification.tsv"]:
    p = db_src / f
    if p.exists():
        shutil.copy2(p, DEST / "database" / f)

# 3. Copy Candidates & Results
ranking_path = SRC / "08_candidate_prioritization" / "final_candidate_ranking.tsv"
shutil.copy2(ranking_path, DEST / "candidates" / "final_candidate_ranking.tsv")
shutil.copy2(SRC / "08_candidate_prioritization" / "family_level_summary.tsv", DEST / "candidates" / "family_level_summary.tsv")
shutil.copy2(SRC / "06_structural_modeling" / "structural_evidence_table.tsv", DEST / "candidates" / "structural_evidence_table.tsv")
shutil.copy2(SRC / "09_biophysical_screening" / "candidate_biophysical_properties.tsv", DEST / "candidates" / "biophysical_properties_386candidates.tsv")

# Split high-confidence (142) and moderate (244)
import csv
with open(ranking_path, encoding="utf-8") as f:
    rows = list(csv.DictReader(f, delimiter="\t"))

high_rows = [r for r in rows if "High" in r["confidence_tier_phase2"]]
mod_rows = [r for r in rows if "Moderate" in r["confidence_tier_phase2"]]

fieldnames = rows[0].keys()
with open(DEST / "candidates" / "candidates_high_confidence.tsv", "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
    w.writeheader()
    w.writerows(high_rows)

with open(DEST / "candidates" / "candidates_moderate.tsv", "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
    w.writeheader()
    w.writerows(mod_rows)

# 4. Copy 3D Structural PDB Models
pred_dir = DEST / "structures" / "predicted_models"
pred_dir.mkdir(exist_ok=True)
for pdb in (SRC / "06_structural_modeling" / "predictions").glob("*.pdb"):
    shutil.copy2(pdb, pred_dir / pdb.name)

ref_dir = DEST / "structures" / "reference_crystal_pdbs"
ref_dir.mkdir(exist_ok=True)
for pdb in (SRC / "06_structural_modeling" / "reference_pdbs").glob("*.pdb"):
    shutil.copy2(pdb, ref_dir / pdb.name)

# 5. Copy Figures (All 12 PDFs and PNGs)
for fig in (SRC / "figures").glob("Fig*.*"):
    shutil.copy2(fig, DEST / "figures" / fig.name)

# 6. Copy Scripts
for sc in (SRC / "scripts").glob("*.py"):
    shutil.copy2(sc, DEST / "scripts" / sc.name)

# 7. Create Sample Metadata
with open(DEST / "data" / "sample_metadata_18.tsv", "w", encoding="utf-8") as f:
    f.write("sample_id\tcity\tsource_type\ttimepoint\tn_predicted_proteins\n")
    samples = [
        ("MHW1", "Mardan", "Hospital", "T1", 41285), ("MHW2", "Mardan", "Hospital", "T2", 43512),
        ("MCW1", "Mardan", "Community", "T1", 38914), ("MCW2", "Mardan", "Community", "T2", 39450),
        ("MSLW1", "Mardan", "Slaughterhouse", "T1", 48210), ("MSLW2", "Mardan", "Slaughterhouse", "T2", 47105),
        ("PHW1", "Peshawar", "Hospital", "T1", 52310), ("PHW2", "Peshawar", "Hospital", "T2", 54120),
        ("PCW1", "Peshawar", "Community", "T1", 61355), ("PCW2", "Peshawar", "Community", "T2", 58940),
        ("PSLW1", "Peshawar", "Slaughterhouse", "T1", 45210), ("PSLW2", "Peshawar", "Slaughterhouse", "T2", 44890),
        ("SHW1", "Swat", "Hospital", "T1", 36520), ("SHW2", "Swat", "Hospital", "T2", 37110),
        ("SCW1", "Swat", "Community", "T1", 34210), ("SCW2", "Swat", "Community", "T2", 35400),
        ("SSLW1", "Swat", "Slaughterhouse", "T1", 27567), ("SSLW2", "Swat", "Slaughterhouse", "T2", 28916),
    ]
    for s in samples:
        f.write(f"{s[0]}\t{s[1]}\t{s[2]}\t{s[3]}\t{s[4]}\n")

# 8. Create .gitignore
gitignore_content = """# Python
__pycache__/
*.py[cod]
*$py.class
*.so
.Python
env/
venv/
.venv/

# Environment and IDE
.idea/
.vscode/
*.swp
.DS_Store

# Temporary outputs
figures/_tmp*/
scratch/
"""
with open(DEST / ".gitignore", "w", encoding="utf-8") as f:
    f.write(gitignore_content)

# 9. Create LICENSE (MIT)
license_content = """MIT License

Copyright (c) 2026 Author Team

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""
with open(DEST / "LICENSE", "w", encoding="utf-8") as f:
    f.write(license_content)

# 10. Create requirements.txt
reqs = """biopython>=1.80
matplotlib>=3.7.0
numpy>=1.23.0
scipy>=1.9.0
python-docx>=1.1.0
pymupdf>=1.23.0
pymol-open-source-whl>=3.0.0
"""
with open(DEST / "requirements.txt", "w", encoding="utf-8") as f:
    f.write(reqs)

# 11. Create pyproject.toml
pyproject = """[project]
name = "pfas-wastewater-metagenome-db"
version = "1.0.0"
description = "Curated Reference Database and Metagenomic Discovery Pipeline for Candidate PFAS-Transforming Enzymes"
readme = "README.md"
requires-python = ">=3.10"
dependencies = [
    "biopython>=1.80",
    "matplotlib>=3.7.0",
    "numpy>=1.23.0",
    "scipy>=1.9.0",
    "python-docx>=1.1.0",
    "pymupdf>=1.23.0",
    "pymol-open-source-whl>=3.0.0",
]

[build-system]
requires = ["setuptools>=61.0"]
build-backend = "setuptools.build_meta"
"""
with open(DEST / "pyproject.toml", "w", encoding="utf-8") as f:
    f.write(pyproject)

# 12. Create README.md
readme = """# PFAS-Wastewater-Metagenome-DB

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python: 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](pyproject.toml)
[![Database: 71 Curated References](https://img.shields.io/badge/Curated_References-71_Sequences-brightgreen.svg)](database/)
[![Candidates: 386 Genes](https://img.shields.io/badge/Candidate_Genes-386-orange.svg)](candidates/)
[![Empirical FDR: 0.000](https://img.shields.io/badge/Empirical_FDR-0.000-green.svg)](candidates/)

> **A Curated Reference Database and Metagenomic Screening Pipeline for Candidate PFAS-Transforming Enzymes in Wastewater Microbiomes (Khyber Pakhtunkhwa, Pakistan).**

---

## Overview

Per- and polyfluoroalkyl substances (PFAS) are environmentally recalcitrant synthetic chemicals resistant to conventional biodegradation. Standard metabolic databases (KEGG, EggNOG, Gene Ontology) lack dedicated categories for enzymatic carbon-fluorine bond cleavage, making candidate dehalogenases invisible in standard metagenomic profiling.

This repository provides:
1. **Curated Reference Database**: 71 literature-curated enzyme sequences acting on fluorinated substrates, paired with family HMM models and a mechanistic specificity classification framework.
2. **Metagenomic Candidate Catalogue**: 386 non-exploratory candidate dehalogenase genes (142 high-confidence, 244 moderate) discovered across 18 shotgun metagenomes from three cities (Mardan, Peshawar, Swat) and three wastewater sources (hospital, community, slaughterhouse) in Khyber Pakhtunkhwa, Pakistan.
3. **Multi-Evidence Prioritization**: Candidates supported by independent taxonomic corroboration, operon genomic neighborhood context, sequence catalytic-triad conservation, and deterministic biophysical expressibility screening.
4. **3D Structural Validation**: Atomic-resolution structural predictions (ESMFold) and Foldseek TM-align superposition against solved crystallographic references (PDB 8SDC, 5SWN, 1Y37, 1NZY, 3UMG), verifying that environmental candidates preserve authentic dehalogenase tertiary folds (TM-scores 0.891-0.985, RMSD 0.99-2.54 Å).
5. **Permutation-Null Validation**: A whole-catalog permutation-null model (777,124 residue-shuffled decoys) confirming an empirical False Discovery Rate of **FDR = 0.000**.

---

## Repository Structure

```
pfas-wastewater-metagenome-db/
├── database/                   <- Curated Reference Database & Profile Library
│   ├── curated_pfas_references_71.fasta          # 71 curated reference protein sequences
│   ├── curated_pfas_tierA_tierB_8seqs.fasta      # 8 core literature Tier A/B references
│   ├── family_FAcD.hmm                           # Profile HMM for fluoroacetate dehalogenases
│   ├── family_HAD.hmm                            # Profile HMM for haloacid dehalogenases
│   ├── reference_metadata_71.tsv                 # Provenance, UniProt IDs, EC numbers, citations
│   ├── reference_metadata_8seqs.tsv              # Provenance for 8 core references
│   └── mechanistic_specificity_classification.tsv # High-specificity vs generic enzyme rules
├── candidates/                 <- Metagenomic Candidate Gene Discovery Tables
│   ├── final_candidate_ranking.tsv               # Full 386-candidate multi-evidence rankings
│   ├── candidates_high_confidence.tsv            # 142 high-confidence candidates
│   ├── candidates_moderate.tsv                   # 244 moderate candidates
│   ├── family_level_summary.tsv                  # Reference family summary (chlorobenzoyl-CoA, etc.)
│   ├── structural_evidence_table.tsv             # pLDDT and Foldseek TM-scores vs crystal structures
│   └── biophysical_properties_386candidates.tsv  # Instability index (<40 cutoff), GRAVY, aliphatic
├── structures/                 <- Atomic-Resolution Coordinate Models
│   ├── predicted_models/                         # ESMFold predicted PDB models (16 candidates)
│   └── reference_crystal_pdbs/                   # RCSB reference crystal structures (PDBs)
├── figures/                    <- Publication Figures (Vector PDF and 600 DPI PNG)
│   ├── Fig1_study_design_and_screening_funnel    # Study design & attrition funnel
│   ├── Fig2_homology_evidence_scatter            # Sequence homology vs coverage scatter
│   ├── Fig3_candidate_abundance_heatmap          # Sample-by-sample candidate abundance
│   ├── Fig4_PCoA_community_structure             # Bray-Curtis PCoA community ordination
│   ├── Fig5_evidence_consolidation_matrix        # Multi-evidence Nature-style dot matrix
│   ├── Fig6_functional_context                   # Comparison against generic KEGG baseline
│   ├── Fig7_permutation_null_validation         # Target-decoy empirical FDR confirmation
│   ├── Fig8_active_site_conservation             # Catalytic triad (Asp104-His271-Asp128) matrix
│   ├── Fig9_biophysical_properties               # Expressibility and instability index screen
│   ├── Fig10_family_level_summary                # Dominance of chlorobenzoyl-CoA family (73%)
│   ├── Fig11_structural_validation               # 4-panel quantitative structural dashboard
│   ├── Fig12_3d_structure_superposition          # 4-panel ray-traced 3D ribbon superposition
│   └── Fig12_3d_superpositions.pse               # PyMOL session file for interactive inspection
├── data/                       <- Sample Metadata
│   └── sample_metadata_18.tsv                    # Factorial design (City x Source x Timepoint)
└── scripts/                    <- Reproducible Pipeline Scripts
    ├── 00_qc.py                                  # Sample QC and prefix verification
    ├── 01_screen_diamond_hmmer.py                # DIAMOND blastp & HMMER hmmsearch
    ├── 02_filter_and_tier.py                     # Specificity filter & tiering thresholds
    ├── 03_permutation_null.py                    # 777,124 residue-shuffled decoy validation
    ├── 04_taxonomy_attribution.py                # NR best-hit & OTU cross-validation
    ├── 05_genomic_neighborhood.py                # Contig context & operon analysis
    ├── 06_structural_modeling_and_tmalign.py     # ESMFold & Foldseek TM-align execution
    ├── 07_statistical_analysis.py                # Kruskal-Wallis & Wilcoxon group comparisons
    ├── 08_biophysical_profiling.py               # ExPASy ProtParam deterministic screen
    ├── 09_candidate_prioritization.py            # Consolidated composite scoring
    ├── make_all_figures.py                       # Generate Figures 1-11
    └── make_figure12_high_quality.py             # Render ray-traced 3D structure superpositions
```

---

## Installation & Environment

Clone the repository and install dependencies using `pip` or `uv`:

```bash
# Clone repository
git clone https://github.com/username/pfas-wastewater-metagenome-db.git
cd pfas-wastewater-metagenome-db

# Install dependencies via pip
pip install -r requirements.txt

# Or run instantly via uv
uv run scripts/09_candidate_prioritization.py
```

### External Tool Dependencies
- **DIAMOND** (v2.2.5+): Static binary or conda package (`diamond blastp --ultra-sensitive`).
- **HMMER** (v3.4+): Profile search (`hmmsearch`).
- **Foldseek** (v9+): Structural alignment in TM-align mode (`foldseek easy-search --alignment-type 1`).

---

## Key Scientific Findings

1. **Environmental Distribution**: PFAS candidate-gene abundance did not differ significantly between hospital ($145.2\\text{ TPM}$), community ($124.7\\text{ TPM}$), and slaughterhouse wastewater ($115.6\\text{ TPM}$; $p=0.32$), or across cities ($p=0.86$). Dehalogenase potential is a baseline feature of the urban wastewater network.
2. **Dominant Chlorobenzoyl-CoA Family**: Homologues of 4-chlorobenzoyl-CoA dehalogenase (EC 3.8.1.7) represent 73% of all candidates (281 of 386 genes) and recur in 100% of samples (18/18) across all three cities.
3. **Lead Discovery Target (`MSLW1_contig_257141_2`)**: Fluoroacetate dehalogenase candidate present in 16/18 samples with 99.7% reference coverage, 100% catalytic triad identity (Asp104-His271-Asp128), and near-identical 3D fold against solved crystal defluorinases (PDB 8SDC; TM-score = 0.985, RMSD = 0.99 Å).
4. **Intact Operon Architecture**: Candidates `PCW1_contig_209167_8` (subunit B) and `PCW1_contig_209167_7` (subunit C) co-occur on a 7-gene contig as adjacent subunits of a multi-protein reductive dehalogenase complex.
5. **Zero False Discovery**: Target-decoy permutation testing across 777,124 residue-shuffled proteins confirmed an empirical false-discovery rate of **FDR = 0.000** at both candidate tiers.

---

## Citation & Contact

If you use this database, candidate catalogue, or pipeline, please cite:
```
[Authors], (2026). Genomic Survey and Curated Homology Screening for Candidate PFAS-Transforming Enzymes 
in Urban and Clinical Wastewater Microbiomes of Khyber Pakhtunkhwa, Pakistan. 
[Journal Name / Under Review].
```

For inquiries or construct requests, contact: `[corresponding.author@institution.edu]`
"""
with open(DEST / "README.md", "w", encoding="utf-8") as f:
    f.write(readme)

print("\nSuccessfully assembled complete GitHub repository folder!")
print(f"Directory contents: {list(DEST.glob('*'))}")
