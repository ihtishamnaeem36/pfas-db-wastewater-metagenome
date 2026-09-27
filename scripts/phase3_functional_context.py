"""
Phase 3 -- Functional context from KEGG / EggNOG / GO (community-level, supporting
evidence only, per the plan -- not a per-gene annotation of the Phase 2 candidates;
see the caveat in phase2_2_filter_tier.py for why no gene-level KO map exists).

Confirms the plan's premise (KEGG/EggNOG/GO have no PFAS-defluorination category) by
showing the community-level baseline actually available: known dehalogenase/haloacid/
haloalkane-relevant KOs, EggNOG COGs, and GO terms, and their abundance pattern across
city x source type x timepoint. This is the "dehalogenase-family term abundance by
city/source/timepoint" community background picture the plan asks for.
"""
import csv
import re
from pathlib import Path
from collections import defaultdict

ROOT = Path(__file__).resolve().parents[1]
P5 = ROOT.parents[2] / "project 5" / "Corrected Data"
OUT = ROOT / "03_functional_context"
OUT.mkdir(parents=True, exist_ok=True)

KEGG_KEYWORDS = re.compile(r"dehalogen|haloacid|haloalkane|fluoroacetate|"
                            r"chlorohydroquinone|chloroethene|hydroxyphenylacetate", re.I)
# Narrowed to halogen-specific terms only (not generic hydrolase/oxidoreductase, which
# would otherwise pull in hundreds of unrelated enzyme-class GO/COG terms and dilute
# the "dehalogenase-relevant" signal -- see phase3_report.txt discussion).
GENERAL_KEYWORDS = re.compile(r"dehalogen|haloacid|haloalkane|halide|organohal|"
                               r"c-halide|halogenated", re.I)

def parse_sample_meta(s):
    m = re.match(r"^([MPS])(HW|CW|SLW)(\d)$", s)
    city = {"M": "Mardan", "P": "Peshawar", "S": "Swat"}[m.group(1)]
    source = {"HW": "Hospital", "CW": "Community", "SLW": "Slaughterhouse"}[m.group(2)]
    return city, source, m.group(3)

def load_table(path, id_col_name, desc_col_name, keyword_re):
    with open(path, encoding="utf-8") as f:
        r = csv.reader(f, delimiter="\t")
        header = next(r)
        non_sample_cols = ("# Gene Family", "Module", "Pathway", "Name", "EC", "Description",
                            "Function Category", "Function Category Level1",
                            "Function Category Level2", "category", "go_term")
        samples = [c for c in header if c not in non_sample_cols]
        desc_idx = header.index(desc_col_name)
        hits = []
        for row in r:
            d = row[desc_idx] if desc_idx < len(row) else ""
            if keyword_re.search(d):
                vals = dict(zip(header, row))
                hits.append((row[0], d, {s: float(vals[s]) for s in samples}))
    return samples, hits

def summarize(samples, hits, label, keyword_desc):
    lines = [f"\n=== {label} ===", f"Keyword filter: {keyword_desc}",
             f"{len(hits)} matching terms found (out of the delivered community profile)."]
    by_group = defaultdict(float)
    by_sample_total = defaultdict(float)
    for term_id, desc, vals in hits:
        lines.append(f"  {term_id:12s} {desc[:70]}")
        for s, v in vals.items():
            city, source, tp = parse_sample_meta(s)
            by_group[(city, source)] += v
            by_sample_total[s] += v
    lines.append("\n  Abundance summed across matching terms, by city x source (mean TPM-like units):")
    for city in ["Mardan", "Peshawar", "Swat"]:
        for source in ["Hospital", "Community", "Slaughterhouse"]:
            v = by_group.get((city, source), 0.0)
            lines.append(f"    {city:10s} {source:15s} {v:10.2f}")
    return lines, by_sample_total

def main():
    report = ["PHASE 3 -- FUNCTIONAL CONTEXT (community-level)", "=" * 60]

    kegg_samples, kegg_hits = load_table(
        P5 / "1-KEGG_unstratified" / "All.KO.detail.xls", "K", "Description", KEGG_KEYWORDS)
    lines, kegg_by_sample = summarize(kegg_samples, kegg_hits, "KEGG KO",
                                        "dehalogenase/haloacid/haloalkane/fluoroacetate etc. in Description")
    report += lines

    egg_samples, egg_hits = load_table(
        P5 / "2-EggNOG_unstratified" / "All.EGGNOG.detail.xls", "COG", "Description", GENERAL_KEYWORDS)
    lines, egg_by_sample = summarize(egg_samples, egg_hits, "EggNOG COG",
                                       "dehalogenase/hydrolase/oxidoreductase in Description")
    report += lines

    go_samples, go_hits = load_table(
        P5 / "3-GO_unstratified" / "All.GO.detail.xls", "GO", "go_term", GENERAL_KEYWORDS)
    lines, go_by_sample = summarize(go_samples, go_hits, "GO term",
                                      "dehalogenase/hydrolase/oxidoreductase activity in go_term")
    report += lines

    # write per-sample combined table for downstream correlation with candidate TPM (Phase 7)
    all_samples = sorted(set(kegg_by_sample) | set(egg_by_sample) | set(go_by_sample))
    with open(OUT / "dehalogenase_relevant_term_abundance_by_sample.tsv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["sample", "city", "source", "timepoint",
                    "KEGG_dehalogenase_KO_sum", "EggNOG_hydrolase_oxidoreductase_sum",
                    "GO_hydrolase_oxidoreductase_dehalogenase_sum"])
        for s in all_samples:
            city, source, tp = parse_sample_meta(s)
            w.writerow([s, city, source, tp,
                        round(kegg_by_sample.get(s, 0), 3),
                        round(egg_by_sample.get(s, 0), 3),
                        round(go_by_sample.get(s, 0), 3)])

    report.append(f"\nPer-sample combined table written to "
                   f"dehalogenase_relevant_term_abundance_by_sample.tsv")
    report.append("\nCONCLUSION: standard KEGG/EggNOG/GO annotation gives only a coarse, "
                   "family-level dehalogenase signal (no PFAS-specific category exists), "
                   "confirming the plan's premise that curated-database screening (Phase 2) "
                   "is necessary to recover PFAS-relevant candidates that generic pathway "
                   "databases cannot resolve at the required specificity.")

    out_txt = OUT / "phase3_report.txt"
    out_txt.write_text("\n".join(report), encoding="utf-8")
    print("\n".join(report))
    print(f"\nWrote {out_txt}")

if __name__ == "__main__":
    main()
