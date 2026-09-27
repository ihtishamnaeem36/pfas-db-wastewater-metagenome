"""
Phase 6b -- Execute 3D structural prediction for priority candidates using ESMFold API.
Fetches reference solved structures from RCSB PDB, predicts structures for the 20 priority
candidates, extracts per-residue pLDDT metrics, and saves predicted PDB models.
"""
import os
import sys
import time
import csv
import urllib.request
import urllib.error
from pathlib import Path
from Bio import SeqIO

ROOT = Path(__file__).resolve().parents[1]
SUBSET_FASTA = ROOT / "06_structural_modeling" / "candidates_for_alphafold_PRIORITY_SUBSET.fasta"
PRED_DIR = ROOT / "06_structural_modeling" / "predictions"
REF_DIR = ROOT / "06_structural_modeling" / "reference_pdbs"
SUMMARY_TSV = ROOT / "06_structural_modeling" / "plddt_summary.tsv"

REF_PDBS = {
    "1Y37": "FAcD DehH1 (Burkholderia sp. FA1)",
    "5SWN": "FAcD RPA1163 (Rhodopseudomonas palustris)",
    "3UMG": "HAD dehalogenase Rha0230 (Rhodococcus jostii RHA1)",
    "8SDC": "DEF2 defluorinase (Dechloromonas aromatica RCB)",
    "1NZY": "4-chlorobenzoyl-CoA dehalogenase (Pseudomonas sp. CBS3)",
    "2V4U": "Triazine chlorohydrolase AtzA (Pseudomonas sp. ADP)",
    "4U3E": "Benzoyl-CoA reductase BCR (Thauera aromatica)",
}

def fetch_reference_pdbs():
    REF_DIR.mkdir(parents=True, exist_ok=True)
    print("--- Fetching Solved Reference PDBs from RCSB ---")
    for pdb_id, desc in REF_PDBS.items():
        out_pdb = REF_DIR / f"{pdb_id}.pdb"
        if out_pdb.exists() and out_pdb.stat().st_size > 1000:
            print(f"  [OK] {pdb_id} already exists ({out_pdb.stat().st_size:,} bytes): {desc}")
            continue
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "PFAS-Pipeline/1.0"})
            with urllib.request.urlopen(req, timeout=20) as r:
                content = r.read()
            out_pdb.write_bytes(content)
            print(f"  [DOWNLOADED] {pdb_id} ({len(content):,} bytes): {desc}")
        except Exception as e:
            print(f"  [ERROR] Failed to fetch {pdb_id}: {e}")

def fold_sequence_esmfold(seq_str, retries=3, timeout=90):
    url = "https://api.esmatlas.com/foldSequence/v1/pdb/"
    data = seq_str.strip().encode("utf-8")
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(url, data=data, headers={"Content-Type": "text/plain"})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                if r.status == 200:
                    return r.read().decode("utf-8")
        except urllib.error.HTTPError as e:
            print(f"    (Attempt {attempt}/{retries}) HTTP Error {e.code}: {e.reason}")
        except Exception as e:
            print(f"    (Attempt {attempt}/{retries}) Error: {e}")
        time.sleep(attempt * 2)
    return None

def analyze_pdb_plddt(pdb_text):
    atom_lines = [l for l in pdb_text.splitlines() if l.startswith("ATOM")]
    if not atom_lines:
        return 0.0, 0.0, 0.0, 0
    # In ESMFold, B-factor (col 60:66) represents predicted LDDT (0.0 - 1.0 or 0 - 100)
    b_factors = [float(l[60:66].strip()) for l in atom_lines]
    # Check scale
    max_val = max(b_factors)
    scale = 100.0 if max_val <= 1.5 else 1.0
    scaled_bf = [b * scale for b in b_factors]
    
    mean_plddt = sum(scaled_bf) / len(scaled_bf)
    pct_high = sum(1 for b in scaled_bf if b >= 70.0) / len(scaled_bf) * 100.0
    pct_very_high = sum(1 for b in scaled_bf if b >= 90.0) / len(scaled_bf) * 100.0
    
    # Residue count (unique residue numbers)
    res_nums = set(int(l[22:26].strip()) for l in atom_lines)
    return mean_plddt, pct_high, pct_very_high, len(res_nums)

def run_predictions():
    PRED_DIR.mkdir(parents=True, exist_ok=True)
    records = list(SeqIO.parse(SUBSET_FASTA, "fasta"))
    print(f"\n--- Predicting 3D Structures for {len(records)} Priority Candidates ---")
    
    results = []
    for i, rec in enumerate(records, 1):
        gid = rec.id
        seq_str = str(rec.seq)
        out_file = PRED_DIR / f"{gid}.pdb"
        
        print(f"[{i}/{len(records)}] {gid} ({len(seq_str)} aa)...", end=" ", flush=True)
        if out_file.exists() and out_file.stat().st_size > 1000:
            print(f"Already predicted ({out_file.stat().st_size:,} bytes). Analyzing...", end=" ")
            pdb_text = out_file.read_text(encoding="utf-8")
        else:
            t0 = time.time()
            pdb_text = fold_sequence_esmfold(seq_str)
            dt = time.time() - t0
            if pdb_text:
                out_file.write_text(pdb_text, encoding="utf-8")
                print(f"Folded in {dt:.1f}s.", end=" ")
            else:
                print(f"FAILED after retries.")
                continue
                
        mean_plddt, pct_high, pct_very_high, n_res = analyze_pdb_plddt(pdb_text)
        print(f"Mean pLDDT: {mean_plddt:.1f} (>=70: {pct_high:.1f}%, >=90: {pct_very_high:.1f}%)")
        results.append({
            "gene_id": gid,
            "length_aa": len(seq_str),
            "mean_pLDDT": f"{mean_plddt:.2f}",
            "pct_pLDDT_ge_70": f"{pct_high:.1f}",
            "pct_pLDDT_ge_90": f"{pct_very_high:.1f}",
            "confidence_tier": "Very High" if mean_plddt >= 90 else ("High" if mean_plddt >= 70 else ("Medium" if mean_plddt >= 50 else "Low")),
            "pdb_path": str(out_file.relative_to(ROOT))
        })
    
    with open(SUMMARY_TSV, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["gene_id", "length_aa", "mean_pLDDT", "pct_pLDDT_ge_70", "pct_pLDDT_ge_90", "confidence_tier", "pdb_path"], delimiter="\t")
        w.writeheader()
        w.writerows(results)
    print(f"\nWrote summary metrics to {SUMMARY_TSV}")

def main():
    fetch_reference_pdbs()
    run_predictions()

if __name__ == "__main__":
    main()
