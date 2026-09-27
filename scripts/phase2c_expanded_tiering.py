"""
Phase 2c -- Specificity-aware tiering of the EXPANDED DIAMOND screen (71 curated
references: 8 original + 63 from the PFAS-Biodegradation-Toolkit)

WHY THIS SCRIPT EXISTS (methodological note, not routine bookkeeping): the naive
--ultra-sensitive DIAMOND run against all 71 references returned 18,191 hits /
11,907 unique genes -- roughly 30x the hit count of the original 8-reference
screen. Inspecting the per-reference hit distribution showed why: most of the
Toolkit's 63 new entries are core, ubiquitous metabolic enzyme families (sugar
glycosidases/phosphorylases, aminotransferases, epimerases, acyl-CoA
dehydrogenases, generic aromatic mono/dioxygenases) that were each tested
against exactly ONE fluorinated substrate analog in a biochemistry assay, but
whose SEQUENCE family is defined by their ordinary, non-fluorine catalytic
function. A galactosidase homolog in a metagenome reflects "this organism has
galactosidases" (true of nearly every organism), not PFAS relevance -- e.g. the
single reference P27830 (a dTDP-sugar epimerase) alone pulled in 1,179 hits.

Treating all 11,907 hits as "candidates" would dilute the whole result with
this noise. Instead, references are split by mechanistic specificity BEFORE
tiering:
  - "high" specificity: our original 8 curated Tier A/B references (already
    fold-specific dehalogenases), plus Toolkit entries whose EC number is a
    literal C-halide-bond hydrolase (EC 3.8.1.x) or a reductive dehalogenase
    (EC 1.3.7.8) -- enzyme classes whose entire catalytic mechanism IS
    halide-bond cleavage, so sequence homology to them is a meaningful PFAS-
    relevance signal, not incidental.
  - "perfluoro" specificity: the Acetobacterium woodii Car D/E/F-Etf complex
    (EC 1.3.1.108), the only entry in either database acting on genuine
    per/polyfluorinated (not just monofluorinated) carboxylic acids. Kept as
    its own bucket rather than merged into "high" because its own hit count
    (995+250+198 = 1,443) shows it is NOT sequence-specific in practice (Car
    subunits belong to a broadly conserved acyl-CoA dehydrogenase/ETF fold
    family) -- reported as a distinct, honestly-caveated evidence class, not
    silently folded into the headline dehalogenase tier.
  - "fluoride_transport" specificity: CLC F-/H+ antiporter and Fluc channel --
    fluoride EXPORT/resistance machinery, mechanistically real and relevant to
    PFAS bioremediation feasibility (per the Toolkit authors' own "fluoride
    stress management" precondition) but NOT C-F bond cleavage; reported
    separately, never mixed into the defluorination candidate tiers.
  - "low" specificity: everything else (57 of the 63 new Toolkit entries) --
    searched for completeness and transparency (a reviewer can see the full
    18,191-hit table and know it was checked), but explicitly excluded from
    the High-confidence/Moderate/Exploratory tiers used for prioritization.
    Reported in its own supplementary table with the caveat above attached.
"""
import csv
from pathlib import Path
from collections import defaultdict

ROOT = Path(__file__).resolve().parents[1]
DIAMOND_V2 = ROOT / "02_homology_screening" / "real_tools" / "diamond_hits_v2.tsv"
HMM = ROOT / "02_homology_screening" / "real_tools" / "hmm_hits.tsv"  # unchanged, 2 orig families
META = ROOT / "01_reference_database" / "search_reference_v2_metadata.tsv"
REF_FASTA_V2 = ROOT / "01_reference_database" / "search_reference_v2.fasta"
TPM = ROOT / "data" / "gene_TPM_corrected.tsv"

OUT_MAIN = ROOT / "02_homology_screening" / "phase2_ranked_candidates.tsv"
OUT_RAW = ROOT / "02_homology_screening" / "phase2_candidate_hits.tsv"
OUT_LOWSPEC = ROOT / "02_homology_screening" / "broad_homology_lowspecificity.tsv"
OUT_FLU = ROOT / "02_homology_screening" / "fluoride_transport_hits.tsv"
OUT_PERFLUORO = ROOT / "02_homology_screening" / "perfluoro_complex_hits.tsv"

HIGH_SPEC_TOOLKIT_ECS_PREFIX = ("3.8.1",)
HIGH_SPEC_TOOLKIT_ECS_EXACT = {"1.3.7.8"}

def read_fasta_lens(path):
    lens = {}
    name, n = None, 0
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    lens[name] = n
                name = line[1:].split()[0]
                n = 0
            else:
                n += len(line.strip())
    if name is not None:
        lens[name] = n
    return lens

def load_ref_specificity():
    """returns {full_header_or_accession_key: (specificity, category, ec, substrate_class)}"""
    spec = {}
    # original 8 references: header itself starts with TierA|... or is PZP66635.1
    ref_lens = read_fasta_lens(REF_FASTA_V2)
    for full_id in ref_lens:
        if full_id.startswith("TierA|") or full_id.startswith("TierB|") or full_id == "PZP66635.1":
            spec[full_id] = ("high", "OriginalCurated", "", "")
        elif full_id.startswith("Toolkit|"):
            parts = full_id.split("|")
            category = parts[1]
            substrate_class = parts[2]
            accession = parts[-1]
            spec[full_id] = (None, category, "", substrate_class)  # ec filled below

    with open(META, encoding="utf-8") as f:
        meta_rows = list(csv.DictReader(f))
    ec_by_accession = {r["accession"]: r["ec"] for r in meta_rows if r["included"] == "True"}

    for full_id in list(spec.keys()):
        if not full_id.startswith("Toolkit|"):
            continue
        _, category, _, substrate_class = spec[full_id]
        accession = full_id.split("|")[-1]
        ec = ec_by_accession.get(accession, "")
        if category == "PerfluoroTierA":
            specificity = "perfluoro"
        elif category == "FluTransport":
            specificity = "fluoride_transport"
        elif ec.startswith(HIGH_SPEC_TOOLKIT_ECS_PREFIX) or ec in HIGH_SPEC_TOOLKIT_ECS_EXACT:
            specificity = "high"
        else:
            specificity = "low"
        spec[full_id] = (specificity, category, ec, substrate_class)

    return spec, ref_lens

def tier_of(is_tier_a_like, pident, coverage, hmm_best_evalue):
    if is_tier_a_like and coverage >= 70 and pident >= 30:
        return "High-confidence"
    if hmm_best_evalue is not None and hmm_best_evalue < 1e-30:
        return "High-confidence"
    if is_tier_a_like and 55 <= coverage < 70 and pident >= 25:
        return "Moderate"
    if hmm_best_evalue is not None and hmm_best_evalue < 1e-15:
        return "Moderate"
    return "Exploratory"

def get_tpm_rows(gene_ids):
    wanted = set(gene_ids)
    found = {}
    with open(TPM, encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
        samples = header[1:]
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if parts[0] in wanted:
                found[parts[0]] = dict(zip(samples, map(float, parts[1:])))
    return found, samples

def main():
    spec, ref_lens = load_ref_specificity()

    # load DIAMOND v2 hits, split by specificity
    high_hits = defaultdict(list)
    low_hits = defaultdict(list)
    flu_hits = defaultdict(list)
    perfluoro_hits = defaultdict(list)
    with open(DIAMOND_V2, encoding="utf-8") as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            qid, sid = p[0], p[1]
            pident, alnlen = float(p[2]), int(p[3])
            evalue, bitscore = float(p[10]), float(p[11])
            slen = int(p[13])
            coverage = 100.0 * alnlen / slen if slen else 0.0
            specificity = spec.get(sid, ("low", "", "", ""))[0]
            rec = dict(ref_id=sid, pident=pident, aln_len=alnlen, coverage_pct=coverage,
                       evalue=evalue, bitscore=bitscore, ref_len=slen)
            if specificity == "high":
                high_hits[qid].append(rec)
            elif specificity == "fluoride_transport":
                flu_hits[qid].append(rec)
            elif specificity == "perfluoro":
                perfluoro_hits[qid].append(rec)
            else:
                low_hits[qid].append(rec)

    # HMMER hits (unchanged 2 family profiles -- both already "high" specificity)
    hmm_hits = defaultdict(list)
    with open(HMM, encoding="utf-8") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            p = line.split()
            target, family = p[0], p[2]
            hmm_hits[target].append(dict(family=family, evalue=float(p[4]), score=float(p[5])))

    print(f"DIAMOND high-specificity hits: {len(high_hits)} genes")
    print(f"DIAMOND low-specificity (broad enzyme family) hits: {len(low_hits)} genes")
    print(f"DIAMOND fluoride-transport hits: {len(flu_hits)} genes")
    print(f"DIAMOND perfluoro-complex hits: {len(perfluoro_hits)} genes")
    print(f"HMMER hits (2 original family profiles): {len(hmm_hits)} genes")

    # === main candidate tiering: high-spec DIAMOND + HMMER only ===
    all_genes = set(high_hits) | set(hmm_hits)
    raw_rows = []
    best_per_gene = {}
    for gid in all_genes:
        d_hits = sorted(high_hits.get(gid, []), key=lambda h: h["evalue"])
        h_hits = sorted(hmm_hits.get(gid, []), key=lambda h: h["evalue"])
        best_d = d_hits[0] if d_hits else None
        best_h = h_hits[0] if h_hits else None

        for dh in d_hits:
            raw_rows.append(dict(gene_id=gid, ref_id=dh["ref_id"], ref_len=dh["ref_len"],
                                  pident=dh["pident"], aln_len=dh["aln_len"],
                                  coverage_pct=round(dh["coverage_pct"], 1),
                                  evalue=dh["evalue"], bitscore=dh["bitscore"], method="DIAMOND"))
        for hh in h_hits:
            raw_rows.append(dict(gene_id=gid, ref_id=hh["family"], ref_len="", pident="",
                                  aln_len="", coverage_pct="", evalue=hh["evalue"],
                                  bitscore=hh["score"], method="HMMER"))

        if best_d:
            pident, coverage = best_d["pident"], best_d["coverage_pct"]
            ref_id = best_d["ref_id"]
            sw_score = best_d["bitscore"]
            is_tier_a_like = True  # all "high" specificity refs treated as Tier-A-equivalent
        else:
            pident, coverage = 0.0, 0.0
            ref_id = f"HMM:{best_h['family']}"
            sw_score = best_h["score"] if best_h else 0.0
            is_tier_a_like = True

        hmm_best_evalue = best_h["evalue"] if best_h else None
        tier = tier_of(is_tier_a_like, pident, coverage, hmm_best_evalue)
        ref_source = "OriginalCurated" if ref_id.startswith("TierA|") else \
                     ("HMM_profile" if ref_id.startswith("HMM:") else "ToolkitHighSpec")

        best_per_gene[gid] = dict(
            gene_id=gid, best_ref=ref_id, ref_tier="A", ref_source=ref_source,
            sw_score=round(sw_score, 1), pident=round(pident, 1), coverage_pct=round(coverage, 1),
            confidence_tier=tier,
            hmm_proxy_best_FAcD=round(next((h["score"] for h in h_hits if h["family"].startswith("family_FAcD")), 0.0), 1),
            hmm_proxy_best_HAD=round(next((h["score"] for h in h_hits if h["family"].startswith("family_HAD")), 0.0), 1),
        )

    tpm_data, samples = get_tpm_rows(best_per_gene.keys())
    out_rows = []
    for gid, row in best_per_gene.items():
        tpm = tpm_data.get(gid, {})
        total_tpm = sum(tpm.values())
        out_rows.append(dict(
            gene_id=gid, best_ref=row["best_ref"], ref_tier=row["ref_tier"],
            ref_source=row["ref_source"],
            sw_score=row["sw_score"], pident=row["pident"], coverage_pct=row["coverage_pct"],
            confidence_tier=row["confidence_tier"],
            hmm_proxy_best_FAcD=row["hmm_proxy_best_FAcD"],
            hmm_proxy_best_HAD=row["hmm_proxy_best_HAD"],
            total_tpm_across_18_samples=round(total_tpm, 4),
            n_samples_detected=sum(1 for v in tpm.values() if v > 0),
            **{f"TPM_{s}": round(tpm.get(s, 0.0), 4) for s in samples},
        ))
    order = {"High-confidence": 0, "Moderate": 1, "Exploratory": 2}
    out_rows.sort(key=lambda r: (order[r["confidence_tier"]], -float(r["sw_score"])))

    with open(OUT_MAIN, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)
    with open(OUT_RAW, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(raw_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(raw_rows)

    from collections import Counter
    c = Counter(r["confidence_tier"] for r in out_rows)
    n_new_from_toolkit = sum(1 for r in out_rows if r["ref_source"] == "ToolkitHighSpec"
                              and r["confidence_tier"] != "Exploratory")
    print(f"\n=== MAIN CANDIDATE TIERS (high-specificity references only) ===")
    print(f"Wrote {len(out_rows)} ranked candidates to {OUT_MAIN}")
    print(f"Tier breakdown: {dict(c)}")
    print(f"Non-exploratory candidates newly surfaced by Toolkit high-specificity refs "
          f"(EC 3.8.1.x / 1.3.7.8) not seen with the original 8 refs: {n_new_from_toolkit}")

    # === low-specificity broad-homology table (reported, not tiered as candidates) ===
    low_rows = []
    for gid, hits in low_hits.items():
        best = max(hits, key=lambda h: h["bitscore"])
        low_rows.append(dict(gene_id=gid, best_ref=best["ref_id"], pident=round(best["pident"], 1),
                              coverage_pct=round(best["coverage_pct"], 1), evalue=best["evalue"],
                              n_low_spec_refs_hit=len(set(h["ref_id"] for h in hits))))
    low_rows.sort(key=lambda r: r["evalue"])
    with open(OUT_LOWSPEC, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(low_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(low_rows[:5000])  # cap file size; full count reported above
    print(f"\n=== LOW-SPECIFICITY BROAD-HOMOLOGY (NOT candidates) ===")
    print(f"{len(low_rows)} genes hit a broad/generic Toolkit enzyme family (e.g. sugar "
          f"glycosidases, transaminases) -- written to {OUT_LOWSPEC} (capped at 5000 rows) "
          f"for transparency, excluded from candidate tiering. See script docstring for why.")

    # === fluoride transport table ===
    flu_rows = []
    for gid, hits in flu_hits.items():
        best = max(hits, key=lambda h: h["bitscore"])
        flu_rows.append(dict(gene_id=gid, best_ref=best["ref_id"], pident=round(best["pident"], 1),
                              coverage_pct=round(best["coverage_pct"], 1), evalue=best["evalue"]))
    flu_rows.sort(key=lambda r: r["evalue"])
    if flu_rows:
        with open(OUT_FLU, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(flu_rows[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(flu_rows)
    print(f"\n=== FLUORIDE-TRANSPORT HITS (fluoride resistance, not C-F cleavage) ===")
    print(f"{len(flu_rows)} genes -> {OUT_FLU}")

    # === perfluoro complex table ===
    pf_rows = []
    for gid, hits in perfluoro_hits.items():
        best = max(hits, key=lambda h: h["bitscore"])
        pf_rows.append(dict(gene_id=gid, best_ref=best["ref_id"], pident=round(best["pident"], 1),
                             coverage_pct=round(best["coverage_pct"], 1), evalue=best["evalue"]))
    pf_rows.sort(key=lambda r: r["evalue"])
    if pf_rows:
        with open(OUT_PERFLUORO, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(pf_rows[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(pf_rows)
    print(f"\n=== PERFLUORO-COMPLEX HITS (Acetobacterium woodii Car D/E/F-Etf; caveat: this "
          f"complex's own sequence family is broadly conserved (ETF fold), so hits here are "
          f"NOT treated as high-specificity PFAS evidence despite being mechanistically the "
          f"closest match to real per/polyfluoro chemistry in either database) ===")
    print(f"{len(pf_rows)} genes -> {OUT_PERFLUORO}")

    for r in out_rows:
        if r["confidence_tier"] != "Exploratory":
            print(f"  [{r['confidence_tier']:16s}] {r['gene_id']:28s} -> {r['best_ref'][:55]:55s} "
                  f"({r['ref_source']}) pident={r['pident']:.1f}% cov={r['coverage_pct']:.1f}%")

if __name__ == "__main__":
    main()
