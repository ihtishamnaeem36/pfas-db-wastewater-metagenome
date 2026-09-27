"""
Phase 4 -- Taxonomic attribution of candidates

Method (b) from the plan: NR-protein best-hit / LCA taxonomy per gene (direct,
gene-ID-resolved -- used as the primary taxonomy source here).

Method (a) from the plan (stratified functional tables via g__/s__ regex) is not
applicable at single-gene resolution: the delivered stratified KEGG/EggNOG/GO tables
stratify community-level KO/COG/GO abundance by taxon, not individual predicted
proteins, so there is no per-candidate-gene row to regex out of them (see the caveat
already logged in phase2_2_filter_tier.py for why no gene-to-KO map exists in this
delivery). Two independent sources are used instead, matching the plan's intent of
avoiding a single-source taxonomy call:
  (b) NR best-hit + LCA taxonomy (direct, per gene ID)
  (c) OTU table cross-check: is the LCA genus also observed (any abundance > 0) in the
      OTU table for the same originating sample? Independent evidence the organism
      (or a close relative) is really present in that sample, not just an artifact of
      a single low-coverage contig's best BLAST-style hit.
Disagreement/no-OTU-support is flagged explicitly as lower confidence, per the plan.

KNOWN LIMITATION: NR.protein.taxonomy.txt (2,761,184 rows) is larger than the
non-redundant gene catalog it is paired with (777,124 predicted proteins), and a
subset of candidate gene IDs have no row in it at all. The most likely explanation is
that this taxonomy annotation was run pre-clustering (i.e. on all raw ORFs before
non-redundancy/representative-sequence collapse), so representative-gene IDs in the
post-clustering NR catalog do not all have a 1:1 row here. Genes with no taxonomy row
are reported as such (LCA/best-hit = "NA") rather than silently dropped or guessed --
this is itself flagged as lower-confidence evidence per the plan's instruction to
treat taxonomic disagreement/absence explicitly, not paper over it.
"""
import csv
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
P5 = ROOT.parents[2] / "project 5" / "Corrected Data"
TAXO = P5 / "NR.protein.taxonomy.txt"
OTU = ROOT.parent / "Phase0_corrected" / "All.Taxa.OTU.corrected.xls"
CANDIDATES = ROOT / "02_homology_screening" / "phase2_ranked_candidates.tsv"
OUT = ROOT / "04_taxonomic_attribution" / "candidate_taxonomy.tsv"

def load_candidates():
    rows = []
    with open(CANDIDATES, encoding="utf-8") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            if row["confidence_tier"] in ("High-confidence", "Moderate"):
                rows.append(row)
    return rows

def load_taxonomy(gene_ids):
    wanted = set(gene_ids)
    out = {}
    with open(TAXO, encoding="utf-8") as f:
        header = f.readline()
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if parts[0] in wanted:
                out[parts[0]] = dict(
                    lca=parts[1], best=parts[2], best_hit_id=parts[3],
                    pident=parts[4], evalue=parts[12], description=parts[14],
                )
                if len(out) == len(wanted):
                    break
    return out

def load_otu_genera_by_sample():
    """sample -> set of genus names with nonzero abundance"""
    with open(OTU, encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
    samples = header[1:-1]  # drop leading #OTU ID and trailing taxonomy
    genera = {s: set() for s in samples}
    with open(OTU, encoding="utf-8") as f:
        next(f)
        for line in f:
            parts = line.rstrip("\n").split("\t")
            tax = parts[-1]
            m = re.search(r"g__([^;]+)", tax)
            genus = m.group(1) if m else None
            if not genus or genus == "":
                continue
            for s, v in zip(samples, parts[1:-1]):
                try:
                    if float(v) > 0:
                        genera[s].add(genus)
                except ValueError:
                    pass
    return genera

def sample_of(gene_id):
    m = re.match(r"^([A-Z]+\d)_contig_", gene_id)
    return m.group(1) if m else None

def genus_from_lca(lca_string):
    parts = [p for p in lca_string.split(";") if p and p.lower() != "unknown"]
    return parts[-1] if parts else None

def main():
    cands = load_candidates()
    tax = load_taxonomy([c["gene_id"] for c in cands])
    otu_genera = load_otu_genera_by_sample()

    rows = []
    for c in cands:
        gid = c["gene_id"]
        t = tax.get(gid, {})
        sample = sample_of(gid)
        lca_genus_or_lower = genus_from_lca(t.get("lca", ""))
        best_taxonomy = t.get("best", "")
        m = re.search(r";([A-Za-z0-9_]+);([A-Za-z0-9_]+)$", best_taxonomy) or \
            re.search(r";([A-Za-z0-9_]+)$", best_taxonomy)
        supported = False
        checked_genus = None
        if sample and sample in otu_genera:
            # try matching any taxonomic token from LCA/best against OTU genera set
            tokens = set(re.split(r";|_", t.get("lca", "") + ";" + best_taxonomy))
            hits = tokens & otu_genera[sample]
            supported = len(hits) > 0
            checked_genus = ",".join(sorted(hits)) if hits else None

        rows.append(dict(
            gene_id=gid, sample=sample, confidence_tier=c["confidence_tier"],
            best_ref=c["best_ref"],
            LCA_taxonomy=t.get("lca", "NA"),
            best_hit_taxonomy=t.get("best", "NA"),
            best_hit_accession=t.get("best_hit_id", "NA"),
            best_hit_pident=t.get("pident", "NA"),
            best_hit_evalue=t.get("evalue", "NA"),
            best_hit_description=t.get("description", "NA"),
            OTU_table_supports_taxonomy=supported,
            OTU_matching_taxa=checked_genus or "none",
        ))

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    print(f"Wrote {len(rows)} taxonomy-resolved candidates to {OUT}\n")
    n_supported = sum(1 for r in rows if r["OTU_table_supports_taxonomy"])
    print(f"OTU-table cross-check: {n_supported}/{len(rows)} candidates have independent "
          f"OTU-table support for their assigned taxon in the originating sample.\n")
    for r in rows:
        flag = "OK " if r["OTU_table_supports_taxonomy"] else "!! "
        print(f"  {flag}{r['gene_id']:28s} [{r['confidence_tier'][:4]}] "
              f"best-hit={r['best_hit_description'][:45]:45s} "
              f"LCA={r['LCA_taxonomy'][:60]}")

if __name__ == "__main__":
    main()
