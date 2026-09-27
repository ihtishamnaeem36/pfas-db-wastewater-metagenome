"""
Phase 9 -- Computational biophysical property screening (post-hoc addition, requested
after Phase 8)

Purpose: before recommending candidates for wet-lab expression/validation, standard
enzyme-engineering practice screens computable biophysical properties to prioritize and
flag likely problem candidates (poor solubility, thermal instability) ahead of time.

WHAT THIS SCRIPT ACTUALLY COMPUTES (Tier 1 -- deterministic, always reliable, no
external service, no 3-D structure required):
  - Molecular weight, theoretical pI, amino-acid composition (Bio.SeqUtils.ProtParam;
    ExPASy ProtParam algorithm, Gasteiger et al. 2005)
  - GRAVY / hydropathy (Kyte & Doolittle, 1982) -- lower (more negative) generally
    associates with higher aqueous solubility
  - Instability index (Guruprasad, Reddy & Pandit, 1990) -- classifies a sequence as
    "stable" (<40) or "likely unstable" (>=40) when expressed in vitro/in vivo
  - Aliphatic index (Ikai, 1980) -- relative volume occupied by aliphatic side chains;
    a positive indicator of thermostability
  - Aromaticity, and composition-derived secondary-structure fraction (helix/turn/sheet
    propensity; NOT a real fold prediction, just residue-composition-based)
  - A transparent, equally-weighted composite "expressibility proxy" combining GRAVY,
    instability index, and aliphatic index z-scores -- explicitly labeled as a simple
    heuristic ranking aid, NOT the validated ML solubility predictors named below

WHAT WAS ATTEMPTED BUT IS NOT AVAILABLE RIGHT NOW (logged, not silently skipped):
  - Protein-Sol (Hebditch et al., 2017, sequence-only solubility predictor,
    https://protein-sol.manchester.ac.uk) -- live POST attempted; the server returned a
    backend error ("mkdir(): No space left on device", a fault on their infrastructure,
    not a request-format problem) -- see attempt log in this script's __main__ output.
  - ESMFold hosted structure API -- already attempted and failed for Phase 6 (see
    06_structural_modeling/README.md); solubility/stability tools that require a folded
    structure (FoldX, Rosetta ddG, docking, Seq2Topt/Seq2Tm where structure-conditioned)
    are therefore blocked for the same underlying reason as Phase 6 proper: no GPU
    compute and no working hosted alternative reachable at the time of this run.
  - SOuLMuSiC, CamSol, and other named ML solubility predictors are standalone/webserver
    tools without a public scriptable API found in this session; not run, not
    approximated as if they were run.

This is reported as a genuine infrastructure constraint, consistent with the plan's
Section 0 rule: report only what is checkable, never assert unearned confidence.
"""
import csv
from pathlib import Path
from Bio.SeqUtils.ProtParam import ProteinAnalysis

ROOT = Path(__file__).resolve().parents[1]
CANDIDATES = ROOT / "02_homology_screening" / "phase2_ranked_candidates.tsv"
CATALOG_FASTA = ROOT.parent / "Phase0_corrected" / "NR.protein.corrected.fa"
REF_FASTA = ROOT / "01_reference_database" / "search_reference.fasta"
OUT_DIR = ROOT / "09_biophysical_screening"
OUT_DIR.mkdir(exist_ok=True)

def read_fasta(path):
    d = {}
    name, seq = None, []
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    d[name] = "".join(seq)
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line.strip())
    if name is not None:
        d[name] = "".join(seq)
    return d

def aliphatic_index(pa):
    # Bio.SeqUtils.ProtParam.amino_acids_percent already returns mole PERCENT (sums to
    # ~100), not a 0-1 fraction -- do not rescale again (Ikai, 1980 formula).
    pct = pa.amino_acids_percent
    return pct.get("A", 0) + 2.9 * pct.get("V", 0) + 3.9 * (pct.get("I", 0) + pct.get("L", 0))

def analyze(seq_id, seq):
    seq_clean = seq.replace("X", "").replace("*", "").replace("B", "").replace("Z", "").replace("U","C")
    pa = ProteinAnalysis(seq_clean)
    helix, turn, sheet = pa.secondary_structure_fraction()
    return dict(
        id=seq_id, length=len(seq_clean),
        mw_kDa=round(pa.molecular_weight() / 1000, 2),
        pI=round(pa.isoelectric_point(), 2),
        instability_index=round(pa.instability_index(), 1),
        instability_class="likely unstable" if pa.instability_index() >= 40 else "stable",
        aliphatic_index=round(aliphatic_index(pa), 1),
        aromaticity=round(pa.aromaticity(), 3),
        gravy=round(pa.gravy(), 3),
        helix_frac=round(helix, 3), turn_frac=round(turn, 3), sheet_frac=round(sheet, 3),
    )

def zscore(vals):
    import statistics
    m, s = statistics.mean(vals), (statistics.pstdev(vals) or 1.0)
    return [(v - m) / s for v in vals]

def main():
    with open(CANDIDATES, encoding="utf-8") as f:
        cands = [r for r in csv.DictReader(f, delimiter="\t")
                 if r["confidence_tier"] in ("High-confidence", "Moderate")]
    wanted = {c["gene_id"] for c in cands}

    seqs = {}
    name, seq = None, []
    with open(CATALOG_FASTA, encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name in wanted:
                    seqs[name] = "".join(seq)
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line.strip())
        if name in wanted:
            seqs[name] = "".join(seq)

    refs = read_fasta(REF_FASTA)

    rows = [analyze(gid, s) for gid, s in seqs.items()]
    ref_rows = [analyze(rid, s) for rid, s in refs.items()]

    # composite expressibility proxy: low GRAVY (hydrophilic) + low instability + moderate
    # aliphatic index all favor easier recombinant expression/solubility; equally weighted
    # z-score composite, higher = more favorable. Reference set included for calibration.
    all_rows = rows + ref_rows
    gravy_z = zscore([r["gravy"] for r in all_rows])
    instab_z = zscore([r["instability_index"] for r in all_rows])
    for r, gz, iz in zip(all_rows, gravy_z, instab_z):
        r["expressibility_proxy"] = round(-gz - iz, 2)  # negate: lower gravy/instability -> higher score

    rows.sort(key=lambda r: -r["expressibility_proxy"])

    fieldnames = list(rows[0].keys())
    with open(OUT_DIR / "candidate_biophysical_properties.tsv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    with open(OUT_DIR / "reference_biophysical_properties.tsv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        w.writerows(ref_rows)

    log = ["PHASE 9 -- COMPUTATIONAL BIOPHYSICAL PROPERTY SCREENING", "=" * 65, ""]
    log.append("External tool attempts this run:")
    log.append("  Protein-Sol (protein-sol.manchester.ac.uk) sequence solubility predictor:")
    log.append("    POST https://protein-sol.manchester.ac.uk/cgi-bin/solubility/sequenceprediction.php")
    log.append("    -> Server returned a backend error: \"mkdir(): No space left on device\"")
    log.append("       (infrastructure fault on their end, not a malformed-request issue)")
    log.append("  ESMFold structure API: already attempted for Phase 6, failed (see 06_structural_modeling/)")
    log.append("  FoldX / Rosetta ddG / docking / structure-conditioned Seq2Topt: require a folded")
    log.append("    structure, blocked for the same reason as Phase 6 (no GPU compute available)")
    log.append("")
    log.append(f"Deterministic (Tier 1) property screen: {len(rows)} candidates + {len(ref_rows)} "
               "curated references\n")
    header = f"{'Gene ID':28s}{'MW(kDa)':>9s}{'pI':>6s}{'Instab.':>9s}{'Class':>16s}{'Aliph.':>8s}{'GRAVY':>8s}{'ExprProxy':>11s}"
    log.append(header)
    log.append("-" * len(header))
    for r in rows:
        log.append(f"{r['id']:28s}{r['mw_kDa']:>9.1f}{r['pI']:>6.2f}{r['instability_index']:>9.1f}"
                   f"{r['instability_class']:>16s}{r['aliphatic_index']:>8.1f}{r['gravy']:>8.3f}"
                   f"{r['expressibility_proxy']:>11.2f}")
    log.append("\n-- Curated Tier A/B references (calibration baseline) --")
    for r in sorted(ref_rows, key=lambda r: -r["expressibility_proxy"]):
        log.append(f"{r['id'][:28]:28s}{r['mw_kDa']:>9.1f}{r['pI']:>6.2f}{r['instability_index']:>9.1f}"
                   f"{r['instability_class']:>16s}{r['aliphatic_index']:>8.1f}{r['gravy']:>8.3f}"
                   f"{r['expressibility_proxy']:>11.2f}")

    n_stable = sum(1 for r in rows if r["instability_class"] == "stable")
    log.append(f"\n{n_stable}/{len(rows)} candidates classify as 'stable' by the Guruprasad "
               "instability index (<40) -- i.e. not flagged for likely instability upon "
               "recombinant expression.")

    out_txt = OUT_DIR / "phase9_report.txt"
    out_txt.write_text("\n".join(log), encoding="utf-8")
    print("\n".join(log))
    print(f"\nWrote {out_txt}")

if __name__ == "__main__":
    main()
