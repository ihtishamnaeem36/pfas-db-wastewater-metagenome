"""
Build the submission-ready manuscript docx from the actual pipeline outputs.
Structured per Science of the Total Environment / Journal of Hazardous Materials
conventions: Highlights, Abstract, Keywords, Introduction, Materials and Methods,
Results, Discussion, Limitations, Future Work, CRediT, Data Availability, References.

Strictly written in accordance with the Human Academic Writing Style Guide:
- Jagged sentence-length rhythm alternating short (<15 words) and long sentences.
- Direct copulas (is, are, was, were, has, had, shows, contains).
- Zero banned AI vocabulary (additionally, moreover, furthermore, notably, crucially, pivotal, etc.).
- Zero copula-avoidance verbs (represents, constitutes, occupies, serves as, reflects).
- Zero rhetorical em-dashes and parenthetical asides.
- Zero trailing -ing significance tails.
- Strict British spelling throughout (characterise, prioritisation, programme, catalogue, homologue, modelled).
- Preserves all scientific facts, numbers, taxa, and claim boundaries.
"""
import csv
from pathlib import Path
import docx
from docx import Document
from docx.shared import Inches, Pt, Cm
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT.parent / "PFAS_Candidate_Gene_Manuscript.docx"

d = Document()
for section in d.sections:
    section.left_margin = Cm(2.5)
    section.right_margin = Cm(2.5)

style = d.styles["Normal"]
style.font.name = "Times New Roman"
style.font.size = Pt(11)
style.paragraph_format.space_after = Pt(6)
style.paragraph_format.line_spacing = 1.15

def h_title(t):
    p = d.add_paragraph()
    r = p.add_run(t)
    r.bold = True
    r.font.size = Pt(16)
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    return p

def h1(t, num=None):
    p = d.add_heading((f"{num}. " if num else "") + t, level=1)
    for r in p.runs:
        r.font.name = "Times New Roman"
        r.font.size = Pt(13)
        r.font.color.rgb = docx.shared.RGBColor(0, 0, 0)
    return p

def h2(t, num=None):
    p = d.add_heading((f"{num} " if num else "") + t, level=2)
    for r in p.runs:
        r.font.name = "Times New Roman"
        r.font.size = Pt(11.5)
        r.font.color.rgb = docx.shared.RGBColor(0, 0, 0)
    return p

def para(t, italic=False, bold=False, justify=True, size=11):
    p = d.add_paragraph()
    r = p.add_run(t)
    r.italic = italic
    r.bold = bold
    r.font.size = Pt(size)
    if justify:
        p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
    return p

def bullet(t):
    p = d.add_paragraph(t, style="List Bullet")
    for r in p.runs:
        r.font.size = Pt(11)
    return p

def numbered(t):
    p = d.add_paragraph(t, style="List Number")
    return p

def caption(t):
    p = d.add_paragraph()
    r = p.add_run(t)
    r.font.size = Pt(9.5)
    r.italic = True
    p.alignment = WD_ALIGN_PARAGRAPH.LEFT
    return p

def reference(t):
    p = d.add_paragraph()
    p.paragraph_format.left_indent = Cm(0.75)
    p.paragraph_format.first_line_indent = Cm(-0.75)
    r = p.add_run(t)
    r.font.size = Pt(10)
    return p

def table_from_tsv(path, max_rows=30, col_widths=None, fontsize=7.5, cols=None,
                    header_labels=None):
    with open(path, encoding="utf-8") as f:
        rows = list(csv.reader(f, delimiter="\t"))
    header, data = rows[0], rows[1:max_rows+1]
    if cols is not None:
        idx = [header.index(c) for c in cols]
        header = [header_labels[i] if header_labels else header[idx[i]] for i in range(len(idx))]
        data = [[row[i] for i in idx] for row in data]
    t = d.add_table(rows=1, cols=len(header))
    t.style = "Light Grid Accent 1"
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    t.autofit = False
    t.allow_autofit = False
    tblPr = t._tbl.tblPr
    layout = docx.oxml.shared.OxmlElement("w:tblLayout")
    layout.set(docx.oxml.ns.qn("w:type"), "fixed")
    tblPr.append(layout)
    if col_widths:
        for i, w in enumerate(col_widths):
            t.columns[i].width = Inches(w)
        for row in t.rows:
            for cell, w in zip(row.cells, col_widths):
                cell.width = Inches(w)
    hdr_cells = t.rows[0].cells
    for i, hh in enumerate(header):
        hdr_cells[i].text = hh
        if col_widths:
            hdr_cells[i].width = Inches(col_widths[i])
        for p in hdr_cells[i].paragraphs:
            for r in p.runs:
                r.bold = True
                r.font.size = Pt(fontsize)
    for row in data:
        cells = t.add_row().cells
        for i, v in enumerate(row):
            cells[i].text = v
            if col_widths:
                cells[i].width = Inches(col_widths[i])
    for row in t.rows:
        for cell in row.cells:
            for p in cell.paragraphs:
                for r in p.runs:
                    r.font.size = Pt(fontsize)
    return t

def figure(fn, cap_text, width=6.3):
    png = ROOT / "figures" / fn.replace(".pdf", ".png")
    if png.exists():
        d.add_picture(str(png), width=Inches(width))
        d.paragraphs[-1].alignment = WD_ALIGN_PARAGRAPH.CENTER
    caption(cap_text)


# =============================================================================
# TITLE / AUTHORS / AFFILIATIONS
# =============================================================================
h_title("Genomic Survey and Curated Homology Screening for Candidate "
        "PFAS-Transforming Enzymes in Urban and Clinical Wastewater "
        "Microbiomes of Khyber Pakhtunkhwa, Pakistan")
p = d.add_paragraph()
p.alignment = WD_ALIGN_PARAGRAPH.CENTER
r = p.add_run("[Author names to be inserted]\n[Affiliations to be inserted]\n"
              "[Corresponding author email to be inserted]")
r.italic = True
r.font.size = Pt(10.5)

# =============================================================================
# HIGHLIGHTS
# =============================================================================
h1("Highlights")
for hl in [
    "Reference database expansion using published toolkits yielded 71 curated dehalogenase sequences.",
    "A mechanistic specificity filter excluded 9,921 hits to generic metabolic enzyme families.",
    "The screen identified 386 candidate dehalogenase genes (142 high-confidence, 244 moderate).",
    "Chlorobenzoyl-CoA dehalogenase homologues recurred in all 18 samples across three cities.",
    "Whole-catalogue permutation testing confirmed zero decoy hits at both candidate tiers.",
    "Atomic structural modelling and Foldseek TM-align verified dehalogenase tertiary folds (TM up to 0.985).",
]:
    bullet(hl)

# =============================================================================
# ABSTRACT
# =============================================================================
h1("Abstract")
para(
    "Per- and polyfluoroalkyl substances (PFAS) resist biodegradation because of the high "
    "stability of carbon-fluorine bonds. Few enzymes capable of cleaving these bonds have "
    "been described. Standard databases such as KEGG, EggNOG, and Gene Ontology lack "
    "specific categories for defluorination. Consequently, standard profiling cannot "
    "identify candidate genes. Here, we screened shotgun metagenomes from 18 wastewater "
    "samples from Khyber Pakhtunkhwa, Pakistan. The sample set spans three cities (Mardan, "
    "Peshawar, Swat), three source types (hospital, community, slaughterhouse), and two "
    "timepoints. We built a reference database of 71 sequences by combining eight "
    "literature-curated dehalogenases with 63 fluorinated-substrate-acting proteins from "
    "the published PFAS-Biodegradation-Toolkit. The non-redundant catalogue of 777,124 "
    "predicted proteins was searched using DIAMOND blastp and HMMER hmmsearch (E < 1e-5), "
    "yielding 11,907 hits. Most added references belong to core metabolic families tested "
    "against only one fluorinated substrate analog. We therefore applied a mechanistic "
    "specificity filter. We retained only enzymes whose sole known catalytic function is "
    "carbon-halogen bond cleavage (EC 3.8.1.x hydrolases and EC 1.3.7.8 reductive dehalogenases). "
    "Hits to non-specific references (9,921 genes) were excluded from candidate tiers. This "
    "filter left 386 candidate genes: 142 high-confidence and 244 moderate. Chlorobenzoyl-CoA "
    "dehalogenase homologues (EC 3.8.1.7) accounted for 281 candidates (73%). This family was "
    "detected in all 18 samples. Whole-catalogue permutation testing against 777,124 residue-shuffled "
    "decoys produced zero decoy hits reaching either candidate tier (empirical FDR = 0.000). "
    "Taxonomic assignments matched host lineages in the sample OTU table for 187 candidates. "
    "Three-dimensional atomic models generated with ESMFold and aligned with Foldseek confirmed "
    "dehalogenase folds across all priority families (TM-scores 0.891 to 0.985; RMSD 0.99 to 2.54 Å). "
    "The top fluoroacetate dehalogenase candidate (MSLW1_contig_257141_2) aligned to solved "
    "defluorinases (PDB 8SDC and 5SWN) with TM-score 0.985, RMSD 0.99 Å, and full catalytic triad "
    "conservation (Asp104-His271-Asp128). Deterministic biophysical profiling classified 291 of 386 "
    "candidates as sequence-stable. These findings deliver a curated reference database, a "
    "multi-evidence candidate catalogue, and verified structural targets from an understudied "
    "wastewater system."
)

h2("Keywords")
para("PFAS; per- and polyfluoroalkyl substances; metagenomics; wastewater; fluoroacetate "
     "dehalogenase; haloacid dehalogenase; curated homology screening; candidate-gene "
     "discovery; Khyber Pakhtunkhwa", justify=False)

d.add_page_break()

# =============================================================================
# 1. INTRODUCTION
# =============================================================================
h1("Introduction", 1)
para(
    "Per- and polyfluoroalkyl substances (PFAS) are synthetic chemicals containing alkyl "
    "chains with multiple carbon-fluorine bonds. These bonds are among the strongest single "
    "bonds in organic chemistry. Bond strength confers thermal and chemical stability, leading "
    "to widespread industrial use and environmental persistence. Wastewater systems are "
    "established point sources and environmental reservoirs for PFAS. However, the microbial "
    "genes capable of transforming these compounds in wastewater remain largely uncharacterised. "
    "Most published surveillance programmes have focused on high-income regions, leaving low- "
    "and middle-income settings unexamined."
)
para(
    "A major obstacle to metagenomic PFAS research is functional annotation. General pathway "
    "databases such as KEGG, EggNOG, and Gene Ontology lack dedicated categories for carbon-fluorine "
    "bond cleavage. Enzymes with confirmed defluorinating activity are catalogued unevenly across "
    "the primary literature. These include fluoroacetate dehalogenases (FAcD) and selected haloacid "
    "dehalogenase (HAD) superfamily members. Standard profiling tools can detect broad hydrolase "
    "and oxidoreductase families. They cannot identify the specific sequence variants relevant to "
    "organofluorine degradation."
)
para(
    "This study addresses that limitation through curated homology screening. We do not claim "
    "measured PFAS transformation. No paired chemical measurements, such as congener concentrations "
    "or fluoride release rates, were available for these samples. This is a discovery-stage study. "
    "We built a tiered reference database of PFAS-relevant enzymes from the biochemical literature. "
    "We then screened this database against 18 shotgun metagenomes from wastewater systems in "
    "Khyber Pakhtunkhwa, Pakistan. We report candidate genes supported by multiple independent lines "
    "of computational evidence: sequence homology, taxonomic placement, genomic context, "
    "deterministic biophysical stability, and tertiary structural alignment. These candidates serve "
    "as prioritised targets for future laboratory validation."
)
para(
    "This work provides four contributions. First, we assembled a tiered reference database of 71 "
    "enzymes with published activity on fluorinated substrates. Second, we established a screening "
    "protocol validated by a whole-catalogue permutation-null model. Third, we identified 386 candidate "
    "genes supported by taxonomic, genomic context, and atomic structural evidence. Fourth, we provide "
    "the first genomic survey of PFAS-candidate enzymes in wastewater from Khyber Pakhtunkhwa."
)

# =============================================================================
# 2. MATERIALS AND METHODS
# =============================================================================
h1("Materials and Methods", 2)

h2("Study design and sample metadata", "2.1")
para(
    "The study design comprised 18 shotgun metagenomes from Khyber Pakhtunkhwa, Pakistan. The sampling "
    "followed a 3 x 3 x 2 factorial design. Factors were city (Mardan, Peshawar, Swat), source type "
    "(hospital, community, slaughterhouse wastewater), and timepoint (n = 2 independent samples per cell). "
    "Upstream sequencing, assembly, non-redundant gene catalogue construction, per-gene TPM normalisation, "
    "and initial functional annotations were delivered prior to this analysis (Fig. 1A)."
)

h2("Data correction and quality control", "2.2")
para(
    "We corrected two sample-labelling errors before analysis. First, an unassigned sample (MHW3) with "
    "no design counterpart was removed from all data tables. Second, a mislabelled sample (SSLW3) was "
    "verified as the second timepoint for the Swat community source and renamed SCW2. These corrections "
    "were applied across abundance tables, the OTU table, the protein FASTA file, and the gene-level "
    "TPM table. The gene-level TPM table was corrected in this pipeline; other tables had been corrected "
    "upstream. Post-correction inspection verified 18 distinct samples conforming to the factorial design. "
    "No duplicate sample identifiers remained. Column sums for KEGG orthology abundance fell between "
    "617,302 and 810,659 TPM-equivalent units across samples. Contig numbering was consistent across all "
    "18 sample prefixes, giving 777,124 predicted proteins."
)

h2("Curated PFAS-relevant gene and enzyme reference database", "2.3")
para(
    "We compiled a reference set of enzymes with published experimental evidence for carbon-fluorine "
    "bond cleavage. Entries were stratified into two tiers based on evidence quality. Tier A (n = 11) "
    "comprises enzymes with structurally or biochemically confirmed defluorination activity. Tier B "
    "(n = 10) comprises enzymes with proposed, homology-based, or consortium-level activity. Seven "
    "Tier A sequences with public database accessions were retrieved (three FAcD, three HAD, and one DEF2). "
    "One Tier B sequence from Delftia acidovorans (PZP66635.1) had a resolvable accession and was "
    "included, yielding an initial eight-sequence database. We also built two multiple sequence "
    "alignments (FAcD, n = 3; HAD, n = 3; MAFFT --auto) to generate profile hidden Markov models."
)

h2("Reference-database expansion: the PFAS-Biodegradation-Toolkit", "2.3b")
para(
    "We expanded the reference set using the PFAS-Biodegradation-Toolkit (Wackett and Robinson, 2024). "
    "This resource catalogues 64 proteins with reported activity on fluorinated compounds, largely "
    "retrieved from BRENDA. One entry (Q1JU72, DehH1) matched our curated set and was retained once. "
    "The remaining 63 sequences were incorporated, yielding a 71-sequence search database."
)
para(
    "Most entries in this toolkit belong to general metabolic enzyme families. These include sugar "
    "epimerases, acyl-CoA dehydrogenases, and transaminases that were tested against isolated "
    "fluorinated analogues. In a preliminary trial, generic references matched hundreds of non-specific "
    "genes. For example, a single dTDP-sugar epimerase (P27830) matched 1,179 genes in our catalogue. "
    "These hits reflect background metabolic diversity rather than defluorinating capacity."
)
para(
    "We therefore applied a mechanistic specificity filter before screening. We classified references "
    "into four groups. High-specificity references were restricted to enzyme classes whose primary "
    "catalytic function is carbon-halogen bond cleavage: hydrolytic dehalogenases (EC 3.8.1.x) and "
    "reductive dehalogenases (EC 1.3.7.8). Perfluoroalkyl-acid active complexes formed a second group, "
    "containing the CarD/E/F-Etf complex from Acetobacterium woodii (EC 1.3.1.108). Fluoride transport "
    "proteins (CLC F-/H+ antiporter and Fluc channel) formed a third group. The remaining 57 references "
    "were classified as low-specificity metabolic enzymes. Only high-specificity references were used "
    "for candidate tiering."
)

h2("Sequence-similarity and profile-based homology screening", "2.4")
para(
    "We executed DIAMOND (Buchfink et al., 2015) and HMMER (Eddy, 2011) within a local Linux subsystem "
    "(WSL2, Ubuntu) without root privileges. Static binaries for DIAMOND v2.2.5 and HMMER v3.4 were "
    "obtained and executed directly. We built a DIAMOND database from the 71 reference sequences "
    "(diamond makedb). The 777,124 predicted proteins were searched against this database (diamond blastp, "
    "--ultra-sensitive, E < 1e-5, 12 threads). This search identified 18,191 pairwise alignments "
    "representing 11,907 unique genes. In parallel, the two curated family HMM profiles (family_FAcD "
    "and family_HAD) were searched against the catalogue using hmmsearch (E < 1e-5, 12 threads), "
    "recovering 434 hits."
)
para(
    "Hits were separated according to reference specificity. The DIAMOND search yielded 1,776 hits to "
    "high-specificity references, 1,471 hits to the CarD/E/F-Etf complex, 195 hits to fluoride transporters, "
    "and 9,921 hits to low-specificity references. Only high-specificity DIAMOND hits and HMMER hits "
    "were carried forward into candidate tiering. Hits to the remaining categories are reported in "
    "supplementary tables (02_homology_screening/perfluoro_complex_hits.tsv, fluoride_transport_hits.tsv, "
    "and broad_homology_lowspecificity.tsv)."
)
para(
    "Candidates were tiered using identity, coverage, and statistical significance thresholds. "
    "High-confidence candidates required a high-specificity DIAMOND match with >=70% query coverage and "
    ">=30% sequence identity, or an HMM match with E < 1e-30. Moderate candidates required 55-70% coverage "
    "and >=25% sequence identity, or an HMM match with E < 1e-15. These thresholds reflect the divergence "
    "expected for environmental homologues. Catalytic residue conservation and structural alignment were "
    "used as independent mechanistic checks."
)
para(
    "Before establishing the WSL environment, an initial screening was performed using an in-house k-mer "
    "seeded Smith-Waterman aligner in Python against the eight original references. All six high-confidence "
    "candidates identified by the Python implementation were confirmed by DIAMOND and HMMER. Full records "
    "of that preliminary screen are archived in 02_homology_screening/custom_substitute_archive/."
)

h2("Permutation-based null-model validation", "2.5")
para(
    "We estimated the empirical false-discovery rate using a whole-catalogue permutation test. Each of the "
    "777,124 predicted proteins was shuffled independently using the Fisher-Yates algorithm (random "
    "seed = 42). This procedure preserved sequence length and amino acid composition while disrupting "
    "sequence order. The shuffled decoy catalogue was searched against the 71 references using the "
    "identical DIAMOND and HMMER commands. Output alignments were evaluated with the same specificity "
    "filter and tiering criteria. The empirical false-discovery rate was calculated as the ratio of decoy "
    "candidates to real candidates at each tier."
)

h2("Taxonomic attribution", "2.6")
para(
    "Taxonomy was assigned using the pre-computed NR best-hit and lowest common ancestor (LCA) annotations "
    "provided with the catalogue. We cross-referenced these assignments with the sample-level OTU table. "
    "For each candidate, we checked whether taxonomic tokens from its lineage were present at non-zero "
    "abundance in the same sample. Candidates without corresponding entries in the taxonomy table were "
    "recorded as unannotated. We did not impute missing taxonomic data."
)

h2("Genomic neighbourhood analysis", "2.7")
para(
    "We examined contigs harbouring candidate genes to evaluate local genomic context. All predicted "
    "proteins located on the same contig as a candidate were extracted. These flanking genes were annotated "
    "using their NR best-hit descriptions. Genes annotated with functions related to halogen metabolism, "
    "transport, cofactor biosynthesis, or stress response were flagged as corroborating context. Contig "
    "length and gene count were recorded to assess assembly fragmentation."
)

h2("Sequence-based active-site conservation assessment", "2.8")
para(
    "The fluoroacetate dehalogenase DehH1 from Burkholderia sp. FA1 (UniProt Q1JU72) contains an established "
    "catalytic triad: Asp104, His271, and Asp128. For candidates matching the FAcD family, we mapped these "
    "reference positions onto candidate sequences using pairwise alignments. Each aligned residue was "
    "classified as identical, conservatively substituted, non-conservatively substituted, or unaligned. "
    "Conservative substitutions were defined as aspartate-to-glutamate or histidine-to-lysine/arginine. "
    "This sequence comparison is an initial proxy. It does not replace structural evaluation."
)

h2("Atomic-resolution structural modelling and Foldseek TM-align validation", "2.9")
para(
    "We generated three-dimensional models for priority candidates using ESMFold (Lin et al., 2023). "
    "Coordinates were saved in PDB format. Predicted local distance difference test (pLDDT) scores were "
    "recorded in the temperature factor field. Residues with pLDDT >= 70 were classified as confident. "
    "Residues with pLDDT >= 90 were classified as highly confident."
)
para(
    "Experimentally determined crystal structures of confirmed dehalogenases were retrieved from the "
    "RCSB Protein Data Bank. These comprised PDB 1Y37 (Burkholderia sp. FA1 DehH1, 1.80 Å), PDB 5SWN "
    "(Rhodopseudomonas palustris FAcD, 1.55 Å), PDB 8SDC (Dechloromonas aromatica DEF2, 1.85 Å), PDB 3UMG "
    "(Rhodococcus jostii RHA1 HAD, 1.85 Å), PDB 1NZY (Pseudomonas sp. CBS3 4-chlorobenzoyl-CoA dehalogenase, "
    "1.80 Å), PDB 2V4U (Pseudomonas sp. ADP AtzA, 2.60 Å), and PDB 4U3E (Thauera aromatica benzoyl-CoA "
    "reductase, 2.20 Å)."
)
para(
    "Structural superposition was conducted using Foldseek in TM-align mode (--alignment-type 1; "
    "van Kempen et al., 2024). We calculated alignment TM-scores, length-normalized TM-scores, and coordinate "
    "root-mean-square deviations (RMSD, Å). A TM-score >= 0.50 establishes a shared global fold (Zhang "
    "and Skolnick, 2004). A TM-score >= 0.80 establishes high structural fidelity. Active-site coordinates "
    "and pLDDT values were inspected directly against crystallographic residues."
)

h2("Statistical analysis", "2.10")
para(
    "Candidate gene abundance was calculated as summed TPM per sample. Values were compared across cities "
    "and source types using Kruskal-Wallis tests. Pairwise differences were evaluated using Wilcoxon "
    "rank-sum tests. Community differences were assessed by principal coordinates analysis (PCoA; Gower, 1966) "
    "using Bray-Curtis dissimilarity from the OTU table (Bray and Curtis, 1957). Correlations with generic "
    "dehalogenase KOs were tested using Spearman rank correlation. Because each cell contains two samples "
    "(n = 2), all statistical tests are exploratory."
)

h2("Computational biophysical property screening", "2.11")
para(
    "We calculated sequence-based biophysical properties to assess potential for heterologous expression. "
    "Hosted online tools (Protein-Sol and external structure servers) proved unavailable due to server-side "
    "connection limits. We therefore calculated deterministic properties using Biopython ProtParam "
    "(Gasteiger et al., 2005). Parameters included molecular weight, theoretical isoelectric point, the "
    "grand average of hydropathy (GRAVY; Kyte and Doolittle, 1982), the aliphatic index (Ikai, 1980), "
    "and the Guruprasad instability index (Guruprasad et al., 1990). Proteins with an instability index "
    "below 40 are classified as stable in vitro. These properties are closed-form compositional metrics "
    "rather than machine-learning predictions. To assist prioritisation, we computed an expressibility "
    "proxy by summing z-scores of GRAVY and instability index, normalised against the curated reference set."
)

d.add_page_break()

# =============================================================================
# 3. RESULTS
# =============================================================================
h1("Results", 3)

h2("Corrected dataset overview", "3.1")
para(
    "Quality control confirmed 18 complete samples conforming to the factorial design. The catalogue "
    "contained 777,124 predicted proteins across all samples. Per-sample gene counts ranged from 27,567 "
    "to 61,355 (Fig. 1A). KEGG orthology column sums ranged from 617,302 to 810,659 TPM-equivalent units. "
    "This variation matches normal differences in sequencing depth."
)

h2("Curated screening recovers candidates absent from generic annotations", "3.2")
para(
    "Generic functional databases provided poor resolution for dehalogenase discovery. In the delivered "
    "KEGG database, only eight orthology groups mention halogenated substrates. Zero EggNOG groups and "
    "three Gene Ontology terms matched halogen keywords (Fig. 6A). Standard annotation tools cannot "
    "distinguish PFAS-specific enzymes from broad metabolic background."
)
para(
    "Screening against the 71-sequence reference database recovered 11,907 DIAMOND hits from the 777,124 "
    "proteins. After applying the specificity filter, 1,819 genes remained eligible for tiering (Fig. 1B). "
    "Threshold filtering yielded 142 high-confidence and 244 moderate candidates, totalling 386 non-exploratory "
    "genes (Fig. 2). An additional 1,433 high-specificity hits fell below tiering thresholds and were recorded "
    "as exploratory."
)
para(
    "Chlorobenzoyl-CoA dehalogenase homologues (EC 3.8.1.7) dominated the candidate set (Fig. 10). This "
    "family accounted for 281 of 386 candidates (73%). Genes from this family were detected in all 18 samples "
    "across all three cities, with a combined abundance of 1,006 TPM units. The triazine chlorohydrolase "
    "family (EC 3.8.1.8) contributed 22 candidates, also detected in all 18 samples. Curated FAcD and HAD "
    "families accounted for 26 and 19 candidates, detected in 17 and 18 samples, respectively. Reductive "
    "dehalogenase subunits (EC 1.3.7.8) contributed 17 candidates."
)
para(
    "The top-ranked candidate overall was PCW1_contig_209167_8, a benzoyl-CoA reductase subunit B homologue. "
    "It met all six evaluation criteria (composite score 152.9). Its independent NR annotation matched its "
    "curated reference assignment. It was co-located on a seven-gene contig with candidate PCW1_contig_209167_7 "
    "(subunit C), consistent with an intact multi-subunit operon."
)
para(
    "High-confidence candidates showed 22.1% to 73.7% sequence identity and 55.7% to 99.7% coverage against "
    "references (Fig. 2). The candidate with the highest sequence match to FAcD was MSLW1_contig_257141_2. "
    "It occurred in 16 of 18 samples and covered 99.7% of the H-1 reference sequence (Q01398). All six "
    "candidates identified in the preliminary Python screen were recovered here, confirming consistency "
    "across alignment methods."
)

h2("Whole-catalogue permutation-null model confirms screening specificity", "3.3")
para(
    "Screening the 777,124 shuffled decoy sequences against all 71 references produced zero decoy hits in "
    "either candidate tier (empirical FDR = 0.000; Fig. 7). Across the entire decoy catalogue, only 10 "
    "alignments were detected. Only one decoy hit aligned to a high-specificity reference. That hit reached "
    "47.2% coverage, well below the 55% threshold. The remaining nine decoy alignments involved low-specificity "
    "or CarD/E/F references, which are excluded by the specificity filter. These results confirm that the "
    "sequence thresholds and the specificity filter eliminate composition-driven false discoveries."
)

h2("Taxonomic attribution and independent corroboration", "3.4")
para(
    "Taxonomy was resolved for 195 of the 386 candidate genes. The remaining 191 genes lacked entries in the "
    "delivered taxonomy table, indicating incomplete coverage in the upstream annotation file. Resolved "
    "candidates were predominantly Betaproteobacteria, particularly Burkholderiales and Rhodocyclales. This "
    "taxonomic profile matches known environmental dehalogenase hosts."
)
para(
    "The independent OTU table corroborated the assigned taxon in 187 of the 386 candidates. In several "
    "cases, independent annotations matched the screening results directly. Candidate MCW2_contig_261091_1 "
    "was annotated independently as fluoroacetate dehalogenase in the NR database. Candidate PCW1_contig_209167_8 "
    "was annotated independently as benzoyl-CoA reductase subunit B. These independent annotations confirm "
    "the curated assignments."
)

h2("Genomic neighbourhood analysis", "3.5")
para(
    "Several candidate genes occurred in genomic contexts supporting halogen metabolism. Candidate "
    "PHW1_contig_70865_19 (HAD family) was located on an 18-gene contig (5,594 amino acids). This contig "
    "carried two adjacent genes annotated with transport and halogen-metabolism functions. Candidates "
    "PCW1_contig_209167_8 and PCW1_contig_209167_7 occurred together on a seven-gene contig as adjacent "
    "subunits of a reductive dehalogenase complex. Most other candidates occurred on short contigs containing "
    "one to four genes. Contig fragmentation was incorporated as a penalty in the composite prioritisation score."
)

h2("Sequence-based active-site conservation", "3.6")
para(
    "We evaluated catalytic triad conservation (Asp104, His271, Asp128) across 37 FAcD-homologous candidates "
    "using reference DehH1 (Q1JU72). Candidate MSLW1_contig_257141_2 was the only gene retaining all three "
    "catalytic residues (Asp, His, Asp; Fig. 8). Several other candidates retained one or both catalytic "
    "aspartate residues. In those candidates, the histidine residue fell outside the aligned region due to "
    "shorter alignment lengths. This check applies only to the FAcD alpha/beta-hydrolase fold. It is not "
    "applicable to the chlorobenzoyl-CoA or reductive dehalogenase families, which utilise different "
    "catalytic mechanisms."
)

h2("Community-level functional context", "3.7")
para(
    "Summed abundance of the eight generic halogen-related KEGG KOs showed no systematic differences across "
    "cities or source types (Fig. 6A). Candidate gene abundance did not correlate with generic dehalogenase "
    "KO abundance (Spearman rho = -0.32, p = 0.20, n = 18; Fig. 6B). This lack of correlation indicates that "
    "curated screening captures a distinct signal not detected by standard KEGG profiling."
)

h2("Distribution patterns across source type, city, and timepoint", "3.8")
para(
    "Candidate gene abundance did not vary significantly by source type (Kruskal-Wallis H = 2.26, p = 0.32) "
    "or city (H = 0.29, p = 0.86). Bray-Curtis ordination of the OTU table separated communities along two "
    "main axes, explaining 26.8% and 20.1% of variance (Fig. 4). Candidate abundance did not correlate with "
    "community ordination scores. Because no chemical measurements exist for these samples, these distributions "
    "describe gene abundance rather than chemical contamination."
)

h2("Consolidated candidate prioritisation", "3.9")
para(
    "We integrated sequence homology, reference tier, taxonomic corroboration, genomic context, sample "
    "recurrence, and structural alignment into a composite score (Fig. 5). A total of 166 candidates met at "
    "least four of the six evaluation criteria. The top 20 candidates are listed in Table 1."
)
para(
    "The six highest-scoring individual candidates are priority targets for experimental follow-up. Candidate "
    "PCW1_contig_209167_8 ranked first (score 152.8). It is a benzoyl-CoA reductase subunit homologue supported "
    "by independent NR annotation and an adjacent subunit gene. Candidate PHW1_contig_70865_19 ranked second "
    "(score 130.9). It belongs to the HAD family and is situated on an 18-gene contig. Candidates "
    "PCW1_contig_150234_1, MHW2_contig_176972_3, PCW1_contig_8454_3, and PHW2_contig_78442_19 occupied the "
    "next four ranks (scores 117.3 to 120.9). All four belong to the chlorobenzoyl-CoA dehalogenase family."
)
para(
    "Beyond individual genes, the chlorobenzoyl-CoA dehalogenase family is the most widespread signal. "
    "Its 281 candidate genes occurred in all 18 samples across all three cities. This broad recurrence "
    "makes this family a priority for future biochemical assays."
)
_table1_path = ROOT / "08_candidate_prioritization" / "final_candidate_ranking.tsv"
with open(_table1_path, encoding="utf-8") as f:
    _t1rows = list(csv.DictReader(f, delimiter="\t"))
_tmp_path = ROOT / "figures" / "_table1_manuscript.tsv"
with open(_tmp_path, "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["gene_id", "confidence_tier_phase2", "pident", "coverage_pct",
                "n_samples_detected_of_18", "taxonomy_OTU_supported",
                "structural_evidence", "n_evidence_lines_of_5", "composite_score", "priority_tier"])
    for r in _t1rows:
        struct_short = r["structural_evidence"].split(";")[0].replace("Verified: ", "")
        if "Screened" in struct_short:
            struct_short = "Screened"
        elif len(struct_short) > 22:
            struct_short = struct_short[:22]
        w.writerow([
            r["gene_id"], r["confidence_tier_phase2"][:4], r["pident"], r["coverage_pct"],
            r["n_samples_detected_of_18"], "Y" if r["taxonomy_OTU_supported"] == "True" else "N",
            struct_short,
            r["n_evidence_lines_of_5"], r["composite_score"],
            r["priority_tier"].split(" (")[0],
        ])
table_from_tsv(
    _tmp_path, max_rows=20,
    header_labels=["Gene ID", "Tier", "% id.", "% cov.", "n det.\n(/18)", "OTU\nsupp.",
                   "3D Structure\n(PDB/TM)", "Evidence\n(/6)", "Score", "Priority"],
    cols=["gene_id", "confidence_tier_phase2", "pident", "coverage_pct",
          "n_samples_detected_of_18", "taxonomy_OTU_supported",
          "structural_evidence", "n_evidence_lines_of_5", "composite_score", "priority_tier"],
    col_widths=[1.45, 0.50, 0.50, 0.50, 0.50, 0.45, 1.25, 0.55, 0.55, 0.65],
    fontsize=7.5,
)
_tmp_path.unlink()
caption("Table 1. Top 20 ranked candidate genes, integrating sequence homology against high-specificity "
        "references, taxonomic corroboration, genomic neighbourhood, sample recurrence, and atomic-resolution "
        "structural alignment (ESMFold and Foldseek TM-align). Full 386-candidate rankings are provided in "
        "Supplementary Table S1 (08_candidate_prioritization/final_candidate_ranking.tsv).")

h2("Computational biophysical property screening", "3.10")
para(
    "Deterministic biophysical calculations classified 291 of the 386 candidates as sequence-stable "
    "(instability index < 40; Fig. 9). The remaining 95 candidates had instability indices above 40, "
    "indicating potential instability in heterologous expression systems. In comparison, seven of the "
    "eight Tier A/B reference enzymes were classified as stable (mean 32.1; DEF2 was 43.3). Hydropathy "
    "(GRAVY) and aliphatic indices were comparable between candidates and references (Fig. 9B, C)."
)
para(
    "Four top-ranked candidates classified as stable: PCW1_contig_209167_8 (39.8), PHW1_contig_70865_19 "
    "(33.7), PCW1_contig_150234_1 (34.2), and PHW2_contig_78442_19 (36.1). Candidate MSLW1_contig_257141_2 "
    "scored near the boundary (40.1). This borderline score suggests that chaperone co-expression or "
    "solubility tags may assist recombinant expression."
)

h2("Atomic-resolution 3D structural modelling and Foldseek tertiary fold verification", "3.11")
para(
    "ESMFold generated complete structural models for priority candidates (Fig. 11A). The mean global "
    "pLDDT was 84.5% (range: 69.7% to 94.0%). Residues in predicted active-site pockets regularly exceeded "
    "90.0% pLDDT. Two candidates reached mean pLDDT values above 93%: PHW1_contig_283624_2 (94.0%, with "
    "98% of residues >= 70) and PCW1_contig_226866_2 (93.0%)."
)
para(
    "Foldseek TM-align superposition against reference crystal structures confirmed that candidates adopt "
    "dehalogenase folds (Fig. 11B). All tested candidates exceeded the 0.50 same-fold threshold, with "
    "TM-scores from 0.891 to 0.985 (Fig. 11C). Coordinate RMSD values ranged from 0.99 to 2.54 Å (Fig. 11D)."
)
para(
    "Candidate MSLW1_contig_257141_2 aligned closely with solved defluorinases (Fig. 12A). It achieved a TM-score of "
    "0.985 (RMSD 0.99 Å) against Dechloromonas DEF2 (PDB 8SDC), 0.981 (RMSD 1.10 Å) against Rhodopseudomonas "
    "palustris FAcD (PDB 5SWN), and 0.973 (RMSD 1.36 Å) against Burkholderia DehH1 (PDB 1Y37). Its active-site "
    "pocket had a mean pLDDT of 88.8, with spatial alignment of the Asp104-His271-Asp128 catalytic triad (Fig. 12B)."
)
para(
    "Five chlorobenzoyl-CoA dehalogenase candidates aligned to 4-chlorobenzoyl-CoA dehalogenase (PDB 1NZY; Fig. 12C). "
    "These were PHW2_contig_78442_19, MSLW2_contig_7812_2, PCW1_contig_150234_1, PCW1_contig_8454_3, "
    "and MHW2_contig_176972_3. TM-scores ranged from 0.911 to 0.939, with RMSD values of 2.04 to 2.36 Å. "
    "Sequence identity to the reference was 26.1% to 32.3%. The HAD candidate PHW1_contig_70865_19 aligned to "
    "Rhodococcus jostii HAD (PDB 3UMG; Fig. 12D) with a TM-score of 0.946 (RMSD 1.71 Å; pocket pLDDT 89.9). These "
    "alignments demonstrate that candidate enzymes preserve the three-dimensional architecture of catalytic "
    "dehalogenases despite sequence divergence."
)

d.add_page_break()

# =============================================================================
# FIGURES
# =============================================================================
h1("Figures", None)
figure("Fig1_study_design_and_screening_funnel.pdf",
       "Fig. 1. Study design and screening funnel. (A) Factorial study design comprising three cities, "
       "three wastewater sources, and two timepoints (n = 18). (B) Screening workflow from 777,124 predicted "
       "proteins through DIAMOND, HMMER, and the mechanistic specificity filter to candidate tiers.")
figure("Fig2_homology_evidence_scatter.pdf",
       "Fig. 2. Sequence homology metrics for 1,776 high-specificity DIAMOND alignments. Percent identity "
       "is plotted against reference coverage. Symbols indicate reference tier (Tier A circle, Tier B triangle). "
       "Marker size scales with summed candidate abundance (TPM) across all 18 samples.")
figure("Fig3_candidate_abundance_heatmap.pdf",
       "Fig. 3. Abundance heatmap of the top 20 candidate genes across 18 wastewater metagenomes. Values are "
       "log1p-transformed TPM, grouped by city.")
figure("Fig4_PCoA_community_structure.pdf",
       "Fig. 4. Principal coordinates analysis (PCoA) of microbial community structure based on Bray-Curtis "
       "dissimilarity from the OTU table. Marker size indicates summed candidate gene abundance.")
figure("Fig5_evidence_consolidation_matrix.pdf",
       "Fig. 5. Evidence consolidation matrix for the top 20 candidate genes. (A) Matrix showing independent "
       "evidence criteria met per candidate. (B) Composite prioritisation score.")
figure("Fig6_functional_context.pdf",
       "Fig. 6. Community-level functional context. (A) Summed abundance of eight halogen-related KEGG "
       "orthology groups across cities and source types. (B) Candidate gene abundance plotted against generic "
       "dehalogenase KO abundance per sample.")
figure("Fig7_permutation_null_validation.pdf",
       "Fig. 7. Whole-catalogue permutation-null validation. Candidate counts in real versus residue-shuffled "
       "decoy catalogues (n = 777,124), with annotated empirical false-discovery rates (FDR).")
figure("Fig8_active_site_conservation.pdf",
       "Fig. 8. Active-site residue conservation for fluoroacetate dehalogenase candidates evaluated against "
       "the Burkholderia DehH1 catalytic triad (Asp104, His271, Asp128).")
figure("Fig9_biophysical_properties.pdf",
       "Fig. 9. Computational biophysical property screening. (A) Expressibility proxy scores for top and "
       "bottom candidates. (B) Hydropathy (GRAVY) plotted against the Guruprasad instability index for all "
       "386 candidates and reference sequences. (C) Aliphatic index plotted against molecular weight.")
figure("Fig10_family_level_summary.pdf",
       "Fig. 10. Family-level summary of the 386 candidate genes. (A) Candidate counts per reference family, "
       "stratified by confidence tier. (B) Detection breadth across 18 samples, coloured by number of cities "
       "represented.")
figure("Fig11_structural_validation.pdf",
       "Fig. 11. Atomic-resolution three-dimensional structural modelling and Foldseek TM-align validation. "
       "(A) Per-candidate mean and active-site pocket pLDDT scores from ESMFold. (B) Foldseek TM-scores against "
       "reference crystal structures; dashed line marks same-fold threshold (TM = 0.50), dotted line marks high "
       "fidelity (TM = 0.80). (C) TM-scores plotted against sequence identity. (D) Cα coordinate RMSD "
       "distribution against reference structures.")
figure("Fig12_3d_structure_superposition.pdf",
       "Fig. 12. Atomic-resolution three-dimensional structural superposition of representative dehalogenase candidates "
       "against solved crystallographic reference structures. (A) Ribbon diagram of fluoroacetate dehalogenase candidate "
       "MSLW1_contig_257141_2 (cyan) superimposed on Dechloromonas aromatica DEF2 defluorinase (salmon; PDB 8SDC; TM-score = 0.985, "
       "RMSD = 0.99 Å). (B) Zoom-in of the active-site catalytic triad showing spatial alignment of Asp104, His271, and Asp128 "
       "(sticks). (C) 4-Chlorobenzoyl-CoA dehalogenase candidate PHW2_contig_78442_19 (marine) superimposed on Pseudomonas sp. CBS3 "
       "crystal structure (warm pink; PDB 1NZY; TM-score = 0.939, RMSD = 2.04 Å). (D) Haloacid dehalogenase candidate "
       "PHW1_contig_70865_19 (forest green) superimposed on Rhodococcus jostii RHA1 HAD crystal structure (wheat; PDB 3UMG; "
       "TM-score = 0.946, RMSD = 1.71 Å).")

d.add_page_break()

# =============================================================================
# 4. DISCUSSION
# =============================================================================
h1("Discussion", 4)
para(
    "This study demonstrates that curated reference databases identify candidate dehalogenase genes that "
    "are undetectable with generic annotations. General databases such as KEGG and EggNOG lack terms for "
    "fluorinated substrate cleavage (Section 3.2). Consequently, broad functional profiling misses these "
    "enzymes. Purpose-built reference databases are necessary to discover narrow metabolic classes in "
    "metagenomic data. Permutation testing confirmed a low false-discovery rate at the high-confidence tier "
    "(Section 3.3). Decoy sequences produced zero candidates. This validation indicates that the candidate "
    "list does not reflect random sequence matches passing permissive thresholds."
)
para(
    "Expanding the database with the PFAS-Biodegradation-Toolkit increased the reference pool from eight to "
    "71 sequences (Wackett and Robinson, 2024). Naive screening against this expanded set produced 11,907 "
    "hits. Most references in the toolkit are core metabolic enzymes tested against a single fluorinated "
    "compound. In preliminary tests, a single sugar epimerase reference matched 1,179 catalogue genes. These "
    "alignments reflected generic homology rather than dehalogenation capacity. The specificity filter "
    "resolved this problem. Restricting candidates to enzyme classes with dedicated dehalogenase functions "
    "(EC 3.8.1.x and EC 1.3.7.8) reduced the candidate set to 386 genes. This filter kept empirical "
    "false-discovery rates at zero (Section 3.3). Future screens using fluorinated-substrate databases must "
    "account for the primary biochemical function of reference enzymes."
)
para(
    "The most frequent candidate family was chlorobenzoyl-CoA dehalogenase (EC 3.8.1.7). This family was "
    "absent from the initial eight-reference screen. It contributed 281 candidate genes detected across all "
    "18 samples in all three cities (Section 3.2). Whether these homologues act on polyfluorinated compounds "
    "remains untested. Reference enzymes in this class cleave carbon-chlorine bonds on aromatic rings rather "
    "than aliphatic carbon-fluorine bonds. However, they share the hydrolytic dehalogenation mechanism of "
    "fluoroacetate dehalogenases. Their ubiquity across diverse wastewater samples makes them clear targets "
    "for experimental testing."
)
para(
    "Multiple independent lines of evidence converged on a small set of priority candidates. Candidate "
    "PCW1_contig_209167_8 combined sequence homology, independent NR annotation, and operon co-localization. "
    "Candidate MSLW1_contig_257141_2 showed full catalytic triad conservation and a TM-score of 0.985 against "
    "solved defluorinase structures (Section 3.11). Structural alignments confirmed that candidate proteins "
    "adopt genuine dehalogenase folds despite sequence identity below 35%. TM-scores above 0.90 for both FAcD "
    "and chlorobenzoyl-CoA dehalogenase candidates indicate that tertiary structure is conserved."
)
para(
    "Taxonomically resolved candidates were assigned primarily to Betaproteobacteria, including "
    "Burkholderiales and Rhodocyclales. This distribution is consistent with known bacterial hosts of "
    "haloacid and fluoroacetate dehalogenases (Section 2.3). The curated screen thus identifies candidate "
    "genes from bacterial lineages known to carry dehalogenation systems."
)
para(
    "This study provides the first survey of PFAS-relevant candidate genes in wastewater metagenomes from "
    "Khyber Pakhtunkhwa. The sampled sites receive mixed hospital, community, and slaughterhouse waste. "
    "These inputs differ from the municipal treatment plants typically studied in high-income regions. "
    "Candidate abundance did not differ significantly by source type or city (Section 3.8). This even "
    "distribution indicates that candidate dehalogenase genes are widespread across wastewater environments "
    "in this region."
)

# =============================================================================
# 5. LIMITATIONS
# =============================================================================
h1("Limitations", 5)
bullet("No paired chemical measurements (PFAS congener concentrations or fluoride release) were available "
       "for these samples. This study makes no claim of measured chemical transformation.")
bullet("The two curated HMM profiles (family_FAcD and family_HAD) each contain only three sequences. These "
       "profiles detect close homologues but may miss divergent variants. Expanding these alignments with "
       "additional sequences is required.")
bullet("Delivered functional tables contained community-aggregated annotations rather than per-gene KEGG "
       "or EggNOG calls. Novelty relative to standard annotations was therefore assessed at the community level.")
bullet("Delivered taxonomy annotations covered only 195 of the 386 candidate genes. The remaining 191 genes "
       "lacked matching records, likely due to gene-identifier differences between pre- and post-clustering files.")
bullet("Tiering thresholds and the specificity filter were defined after observing the initial DIAMOND and "
       "HMMER output distributions. Both choices were subsequently validated using the permutation-null model, "
       "but they were not established a priori.")
bullet("The specificity filter relies on broad enzyme commission classes (EC 3.8.1.x and EC 1.3.7.8). It does "
       "not evaluate individual substrate profiles. Some excluded hits may possess defluorinating activity, "
       "while some retained hits may lack it.")
bullet("Atomic-resolution structural modelling was performed on a representative subset of priority candidates. "
       "Full structural prediction across all 386 candidates remains to be completed.")
bullet("Because each city-by-source cell contained two samples (n = 2), all statistical group comparisons "
       "are exploratory.")
bullet("Several literature-described defluorinating enzymes (such as DEF1 and Agaricus bisporus laccase) "
       "lack public database accessions and could not be included in the reference set.")
bullet("Biophysical property calculations were sequence-based. They do not substitute for experimental "
       "solubility or expression assays.")

# =============================================================================
# 6. FUTURE WORK
# =============================================================================
h1("Future Work", 6)
para(
    "We identify five priorities for future research. "
    "First, candidate binding pockets should be tested by molecular docking. Target ligands include "
    "trifluoroacetate, PFOA, PFOS, and fluorinated acyl-CoA analogues. "
    "Second, profile HMMs should be constructed for newly recovered families. These include "
    "chlorobenzoyl-CoA dehalogenases and reductive dehalogenase subunits. "
    "Third, reference chlorobenzoyl-CoA dehalogenases should be tested for defluorinating activity in vitro. "
    "Fourth, missing literature accessions should be incorporated as sequences become public. "
    "Fifth, top candidates should be expressed recombinantly. Fluoride release or 19F NMR assays should "
    "be performed to measure catalytic defluorination directly."
)

# =============================================================================
# 7. DATA AVAILABILITY
# =============================================================================
h1("Data Availability", 7)
para(
    "The corrected sample dataset, curated PFAS-relevant reference database (FASTA files, HMM profiles, "
    "and citation tables), analysis pipeline scripts, intermediate tables, and figures are available in the "
    "project repository (PFAS_Pipeline/). Permutation-null decoy sequences can be regenerated using the "
    "fixed random seed (42) specified in Section 2.5."
)

h1("CRediT authorship / Acknowledgments / Declarations", None)
para("[To be completed by the author team prior to submission: CRediT roles, funding "
     "sources, conflict-of-interest statement, and any ethics/permit declarations for "
     "wastewater sample collection.]", italic=True)

d.add_page_break()

# =============================================================================
# REFERENCES
# =============================================================================
h1("References", None)
para("Methodological references:", bold=True, justify=False)
for ref in [
    "Altschul, S.F., Gish, W., Miller, W., Myers, E.W., Lipman, D.J., 1990. Basic local "
    "alignment search tool. J. Mol. Biol. 215, 403-410.",
    "Bray, J.R., Curtis, J.T., 1957. An ordination of the upland forest communities of "
    "southern Wisconsin. Ecol. Monogr. 27, 325-349.",
    "Buchfink, B., Xie, C., Huson, D.H., 2015. Fast and sensitive protein alignment using "
    "DIAMOND. Nat. Methods 12, 59-60.",
    "Cock, P.J.A., Antao, T., Chang, J.T., et al., 2009. Biopython: freely available Python "
    "tools for computational molecular biology and bioinformatics. Bioinformatics 25, "
    "1422-1423.",
    "Eddy, S.R., 2011. Accelerated profile HMM searches. PLoS Comput. Biol. 7, e1002195.",
    "Gasteiger, E., Hoogland, C., Gattiker, A., et al., 2005. Protein identification and "
    "analysis tools on the ExPASy server, in: Walker, J.M. (Ed.), The Proteomics "
    "Protocols Handbook. Humana Press, pp. 571-607.",
    "Guruprasad, K., Reddy, B.V.B., Pandit, M.W., 1990. Correlation between stability of "
    "a protein and its dipeptide composition: a novel approach for predicting in vivo "
    "stability of a protein from its primary sequence. Protein Eng. 4, 155-161.",
    "Hebditch, M., Carballo-Amador, M.A., Charonis, S., Curtis, R., Warwicker, J., 2017. "
    "Protein-Sol: a web tool for predicting protein solubility from sequence. "
    "Bioinformatics 33, 3098-3100.",
    "Ikai, A., 1980. Thermostability and aliphatic index of globular proteins. J. Biochem. "
    "88, 1895-1898.",
    "Kyte, J., Doolittle, R.F., 1982. A simple method for displaying the hydropathic "
    "character of a protein. J. Mol. Biol. 157, 105-132.",
    "Gower, J.C., 1966. Some distance properties of latent root and vector methods used in "
    "multivariate analysis. Biometrika 53, 325-338.",
    "Henikoff, S., Henikoff, J.G., 1992. Amino acid substitution matrices from protein "
    "blocks. Proc. Natl. Acad. Sci. U.S.A. 89, 10915-10919.",
    "Jumper, J., Evans, R., Pritzel, A., et al., 2021. Highly accurate protein structure "
    "prediction with AlphaFold. Nature 596, 583-589.",
    "Kruskal, W.H., Wallis, W.A., 1952. Use of ranks in one-criterion variance analysis. J. "
    "Am. Stat. Assoc. 47, 583-621.",
    "Lin, Z., Akin, H., Rao, R., et al., 2023. Evolutionary-scale prediction of atomic-"
    "level protein structure with a language model. Science 379, 1123-1130.",
    "Mirdita, M., Schutze, K., Moriwaki, Y., Heo, L., Ovchinnikov, S., Steinegger, M., 2022. "
    "ColabFold: making protein folding accessible to all. Nat. Methods 19, 679-682.",
    "Smith, T.F., Waterman, M.S., 1981. Identification of common molecular subsequences. J. "
    "Mol. Biol. 147, 195-197.",
    "van Kempen, M., Kim, S.S., Tumescheit, C., et al., 2024. Fast and accurate protein "
    "structure search with Foldseek. Nat. Biotechnol. 42, 243-246.",
    "Wackett, L.P., Robinson, S.L., 2024. A prescription for engineering PFAS "
    "biodegradation. Biochem. J. 481, 1757-1770. doi:10.1042/BCJ20240283. Data: "
    "github.com/serina-robinson/PFAS-Biodegradation-Toolkit.",
    "Zhang, Y., Skolnick, J., 2004. Scoring function for automated assessment of protein "
    "structure template quality. Proteins 57, 702-710.",
]:
    reference(ref)

para("PFAS/dehalogenase-biochemistry references:", bold=True, justify=False)
for ref in [
    "Atashgahi, S., Lu, Y., Zheng, Y., Saccenti, E., Suarez-Diez, M., Ramiro-Garcia, J., Smidt, H., 2017. "
    "Geochemical and microbial diversity of organohalide-respiring environments. Environ. Microbiol. "
    "19, 4152-4167. doi:10.1111/1462-2920.13876.",
    "Benning, M.M., Wesenberg, G., Liu, R.Q., Dunaway-Mariano, D., Holden, H.M., 1996. The "
    "three-dimensional structure of 4-chlorobenzoyl-CoA dehalogenase: How to cleave an unactivated "
    "carbon-halogen bond. Biochemistry 35, 8103-8109. doi:10.1021/bi9605342.",
    "Chan, P.W.Y., Yakunin, A.F., Edwards, E.A., 2021. Microbial enzymes for defluorination of "
    "per- and polyfluoroalkyl substances (PFAS). Curr. Opin. Biotechnol. 73, 230-237. "
    "doi:10.1016/j.copbio.2021.08.016.",
    "Garcia-Ramos, J.C., et al., 2021. Engineered PLP-dependent enzymes for selective "
    "hydrodefluorination. ChemCatChem 13, 3105-3114. doi:10.1002/cctc.202100344.",
    "Kamachi, T., Nakayama, T., Esaki, N., Yoshizawa, K., 2009. Catalytic mechanism of fluoroacetate "
    "dehalogenase: A computational study. Chem. Eur. J. 15, 7394-7403. doi:10.1002/chem.200802611.",
    "Khusnutdinova, A.N., Flick, R., Popovic, A., Brown, G., Tchigvintsev, A., Nocek, B., "
    "Yakunin, A.F., 2023. Haloacid dehalogenase superfamily reveals hidden enzymatic potential for "
    "cleaving carbon-fluorine bonds. FEBS J. 290, 4812-4828. doi:10.1111/febs.16875.",
    "Kunze, C., Bommer, M., Hagen, W.R., Goris, T., Dobbek, H., 2017. Native mass spectrometric and "
    "EPR studies of reductive dehalogenase complexes. FEBS J. 284, 3920-3935. doi:10.1111/febs.14281.",
    "Kurihara, T., Esaki, N., 2008. Bacterial hydrolytic dehalogenases and related enzymes: occurrences, "
    "reaction mechanisms, and applications. Chem. Rec. 8, 67-74. doi:10.1002/tcr.20141.",
    "Liu, J., et al., 2017. Organohalide-respiring bacteria and reductive dehalogenases: Current "
    "status and perspectives. Environ. Microbiol. Rep. 9, 477-492. doi:10.1111/1758-2229.12574.",
    "Mekureyaw, M.F., Kallenbach, E.M., Fahlman, B.D., et al., 2025. Per- and polyfluoroalkyl substances "
    "(PFAS) degradation potential in freshwater cyanobacteria. J. Hazard. Mater. Adv. 17, 100523. "
    "doi:10.1016/j.hazadv.2024.100523.",
    "Nakayama, T., Kamachi, T., Jitsumori, K., Omi, R., Hirotsu, K., Esaki, N., Kurihara, T., 2012. "
    "Substrate recognition mechanism of fluoroacetate dehalogenase from Rhodopseudomonas palustris "
    "strain CGA009. Chem. Eur. J. 18, 8392-8401. doi:10.1002/chem.201103986.",
    "Wasmund, K., Cooper, M., Schreiber, L., et al., 2014. Single-cell genomics reveals the diversity "
    "and metabolic potential of organohalide-respiring Dehalococcoidia in marine sediments. "
    "Environ. Microbiol. 16, 2140-2156. doi:10.1111/1462-2920.12480.",
]:
    reference(ref)

d.save(str(OUT))
print(f"Wrote {OUT}")
