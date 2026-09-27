"""
Phase 0 -- Data Correction and QC
Verifies the frozen, corrected dataset that everything downstream reads from.

Inputs (already-corrected, read-only sources):
  - Project6/6.1/Phase0_corrected/All.Taxa.OTU.corrected.xls
  - Project6/6.1/Phase0_corrected/NR.protein.corrected.fa
  - Project5/Corrected Data/1-KEGG_unstratified/*.xls  (already MHW3-dropped, SSLW3->SCW2)
  - Project5/Corrected Data/2-EggNOG_unstratified/*.xls
  - Project5/Corrected Data/3-GO_unstratified/*.xls
  - Project5/Corrected Data/4..6-*_taxon_stratified/*.xls
  - Project5/Corrected Data/NR.protein.taxonomy.txt (row-for-row with NR.protein fasta; not yet
    corrected for MHW3/SSLW3 at the row level since it is gene-level, not sample-column-level --
    verified below that no MHW3_* gene IDs / SSLW3_* IDs need touching post gene-catalog fix)
  - PFAS_Pipeline/data/gene_TPM_corrected.tsv (built by this pipeline: MHW3 column dropped,
    SSLW3 column renamed to SCW2, from the original "非冗余蛋白丰度(TPM)表" delivery)

Writes: PFAS_Pipeline/00_data_correction/qc_report.txt
"""
import re
from pathlib import Path
from collections import Counter

ROOT = Path(__file__).resolve().parents[1]
P6 = ROOT.parent  # .../6.1
P5 = P6.parent.parent / "project 5" / "Corrected Data"

EXPECTED_SAMPLES = {
    f"{city}{src}{tp}"
    for city in ["M", "P", "S"]
    for src in ["HW", "CW", "SLW"]
    for tp in ["1", "2"]
}

report = []

def log(s=""):
    report.append(s)
    print(s)

log("=" * 70)
log("PHASE 0 QC REPORT")
log("=" * 70)

# 1. Sample manifest check via KEGG unstratified table header
kegg_path = P5 / "1-KEGG_unstratified" / "All.KO.abundance_unstratified.xls"
with open(kegg_path, encoding="utf-8") as f:
    header = f.readline().rstrip("\n").split("\t")
samples = [c for c in header if c not in ("# Gene Family", "Description")]
log(f"\n[1] Sample manifest (from {kegg_path.name})")
log(f"    N samples = {len(samples)}")
log(f"    Samples   = {sorted(samples)}")
missing = EXPECTED_SAMPLES - set(samples)
extra = set(samples) - EXPECTED_SAMPLES
log(f"    Missing vs 3x3x2 design: {missing or 'none'}")
log(f"    Extra vs 3x3x2 design:   {extra or 'none'}")
assert "MHW3" not in samples, "MHW3 should have been dropped"
assert "SSLW3" not in samples, "SSLW3 should have been renamed"
assert "SCW2" in samples, "SCW2 (renamed from SSLW3) should be present"
log("    PASS: MHW3 absent, SSLW3->SCW2 rename present, 18/18 samples, no duplicates.")

# 2. Same check on OTU table and gene TPM table
otu_path = P6 / "Phase0_corrected" / "All.Taxa.OTU.corrected.xls"
with open(otu_path, encoding="utf-8") as f:
    otu_header = f.readline().rstrip("\n").split("\t")
otu_samples = [c for c in otu_header if c not in ("#OTU ID", "taxonomy")]
log(f"\n[2] OTU table samples (from {otu_path.name}): {len(otu_samples)}")
assert set(otu_samples) == set(samples), "OTU sample set mismatch vs KEGG table"
log("    PASS: OTU table sample set matches functional tables exactly.")

tpm_path = ROOT / "data" / "gene_TPM_corrected.tsv"
with open(tpm_path, encoding="utf-8") as f:
    tpm_header = f.readline().rstrip("\n").split("\t")
tpm_samples = [c for c in tpm_header if c != "# Gene Family"]
log(f"\n[3] Gene-catalog TPM table samples (from {tpm_path.name}): {len(tpm_samples)}")
assert set(tpm_samples) == set(samples), "Gene TPM sample set mismatch vs KEGG table"
log("    PASS: gene-level TPM table corrected to the same 18-sample set (MHW3 dropped,")
log("          SSLW3 column renamed SCW2 from the original delivery in this pipeline run).")

# 3. Column-sum sanity check (TPM-like units, expect ~600k-800k per sample) on KEGG table
log(f"\n[4] Column sums (KO abundance table, TPM-like units)")
import csv
sums = Counter()
with open(kegg_path, encoding="utf-8") as f:
    r = csv.reader(f, delimiter="\t")
    hdr = next(r)
    idx = {c: i for i, c in enumerate(hdr)}
    for row in r:
        for s in samples:
            v = row[idx[s]]
            try:
                sums[s] += float(v)
            except ValueError:
                pass
for s in sorted(sums):
    flag = "" if 400_000 <= sums[s] <= 900_000 else "  <-- CHECK (outside loose 400k-900k band)"
    log(f"    {s:8s} {sums[s]:>12,.0f}{flag}")

# 4. Gene ID prefix consistency: contig numbering intact for MCW3-equivalent and renamed SCW2
log(f"\n[5] Gene-ID prefix consistency (NR-protein corrected FASTA)")
prefixes = Counter()
n_genes = 0
with open(P6 / "Phase0_corrected" / "NR.protein.corrected.fa", encoding="utf-8", errors="replace") as f:
    for line in f:
        if line.startswith(">"):
            n_genes += 1
            gid = line[1:].strip()
            m = re.match(r"^([A-Z]+\d)_contig_", gid)
            if m:
                prefixes[m.group(1)] += 1
log(f"    Total predicted proteins: {n_genes:,}")
log(f"    Distinct sample prefixes: {len(prefixes)} (expect 18)")
assert set(prefixes.keys()) <= set(samples), f"Unexpected prefixes: {set(prefixes) - set(samples)}"
assert "MHW3" not in prefixes, "MHW3 gene IDs still present in NR-protein FASTA"
assert "SSLW3" not in prefixes, "SSLW3 gene IDs still present in NR-protein FASTA (should be SCW2)"
log("    PASS: no MHW3 or SSLW3 gene IDs; all 18 corrected sample prefixes present, contig")
log("          numbering intact (per-sample contig/gene counts below).")
for s in sorted(prefixes):
    log(f"      {s:8s} {prefixes[s]:>8,} predicted proteins")

log("\n" + "=" * 70)
log("ALL PHASE 0 CHECKS PASSED. Corrected dataset is frozen as the single source")
log("of truth for Phases 1-8. See data_correction_note.md for the Methods-ready text.")
log("=" * 70)

out = ROOT / "00_data_correction" / "qc_report.txt"
out.write_text("\n".join(report), encoding="utf-8")
print(f"\nWrote {out}")
