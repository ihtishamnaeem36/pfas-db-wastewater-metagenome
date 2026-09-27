"""
Phase 1b -- Merge the PFAS-Biodegradation-Toolkit into the curated reference database

Source: Wackett, L. and Robinson, S. "A prescription for engineering PFAS
biodegradation." Biochemical Journal (2024) 481(23):1757-1770.
doi:10.1042/BCJ20240283. Data: github.com/serina-robinson/PFAS-Biodegradation-Toolkit
(commit as of 2026-08-22; data file dated 2025-02-07).

64 unique, literature-curated proteins where fluorinated compounds/fluoride are
reaction participants (primarily sourced from BRENDA), spanning 11 substrate
classes (aromatic_F, sugar_F, acylCoA_F, aminoacid_F, othercarboxylate_F,
perfluoroalkenoate_F, organophosphate_F, acetate_F, nucleoside_F,
polysaccharide_F, anion_F). The toolkit's own authors state the same caveat our
Phase 1 citation table already carries: most of these act on monofluorinated
model substrates, not complex PFAS -- they are starting points for engineering,
not confirmed PFAS-degrading enzymes.

One entry (Q1JU72, DehH1/FA1) already exists in our Tier A set -- independent
convergence between two separately curated databases on the same enzyme, kept
once (not duplicated) and logged as a validation point.

Classification into our existing Tier A / Tier B / new categories:
  - All defluorination/dehalogenation enzyme entries -> Tier A (they meet our
    own Tier A definition: structurally/biochemically confirmed activity on a
    fluorinated substrate), same caveat as existing Tier A entries.
  - anion_F entries (CLC F-/H+ antiporter, Fluc channel) are NOT C-F bond-
    cleaving enzymes -- they are fluoride EXPORT/resistance machinery, a
    mechanistically distinct evidence class (relevant to the toolkit authors'
    "fluoride stress management" precondition for real bioremediation, not to
    C-F cleavage itself). Kept as a separate FluTransport category, searched
    but scored/reported separately from the defluorination candidate tiers.
  - perfluoroalkenoate_F entries (Acetobacterium woodii CarD/E/F-Etf reductive
    defluorination complex) are flagged specifically: this is the only entry
    in either database that acts on genuine per/polyfluorinated (not just
    monofluorinated) carboxylic acids -- the closest mechanistic match to
    actual PFAS chemistry available, called out separately in the search
    reference header for downstream reporting.
"""
import csv
from pathlib import Path
import openpyxl

ROOT = Path(__file__).resolve().parents[1]
TOOLKIT_DIR = ROOT / "01_reference_database" / "PFAS_Biodegradation_Toolkit"
TOOLKIT_XLSX = TOOLKIT_DIR / "fluorinated_deg_proteins_table.xlsx"
EXISTING_REF = ROOT / "01_reference_database" / "search_reference.fasta"
OUT_FASTA = ROOT / "01_reference_database" / "search_reference_v2.fasta"
OUT_META = ROOT / "01_reference_database" / "search_reference_v2_metadata.tsv"

EXISTING_ACCESSIONS = {"Q1JU72", "Q01398", "WP_178618037.1", "WP_139652913.1", "PZP66635.1"}
PERFLUORO_CLASS = "perfluoroalkenoate_F"  # Acetobacterium woodii Car complex (H6LGM6/7/8)
ANION_CLASS = "anion_F"

def read_existing():
    seqs = {}
    name, seq = None, []
    with open(EXISTING_REF, encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(seq)
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line.strip())
    if name is not None:
        seqs[name] = "".join(seq)
    return seqs

def main():
    existing = read_existing()
    wb = openpyxl.load_workbook(TOOLKIT_XLSX)
    ws = wb.active
    rows = list(ws.iter_rows(values_only=True))
    header = rows[0]
    data = [dict(zip(header, r)) for r in rows[1:] if r[0]]

    meta_rows = []
    new_records = []  # (header_line, sequence)
    n_skipped_dup = 0
    n_perfluoro = 0
    n_anion = 0
    n_added = 0

    for r in data:
        acc = r["uniprot_id"]
        if acc in EXISTING_ACCESSIONS:
            n_skipped_dup += 1
            meta_rows.append(dict(accession=acc, category="DUPLICATE_OF_EXISTING",
                                   substrate_class=r["major_substrate_class"],
                                   ec=r["ec_class_fourth"], other_id=r["other_id"],
                                   reference=r["reference"], included=False))
            continue

        substrate_class = r["major_substrate_class"] or "unclassified"
        if substrate_class == PERFLUORO_CLASS:
            category = "PerfluoroTierA"
            n_perfluoro += 1
        elif substrate_class == ANION_CLASS:
            category = "FluTransport"
            n_anion += 1
        else:
            category = "ToolkitTierA"

        seq = r["prot"].replace("*", "").strip()
        safe_id = str(r["other_id"]).replace(" ", "_").replace("|", "_")[:40]
        header_line = f"Toolkit|{category}|{substrate_class}|{safe_id}|{acc}"
        new_records.append((header_line, seq))
        n_added += 1
        meta_rows.append(dict(accession=acc, category=category,
                               substrate_class=substrate_class,
                               ec=r["ec_class_fourth"], other_id=r["other_id"],
                               reference=r["reference"], included=True))

    OUT_FASTA.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_FASTA, "w", encoding="utf-8") as f:
        for name, seq in existing.items():
            f.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + "\n")
        for header_line, seq in new_records:
            f.write(f">{header_line}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + "\n")

    with open(OUT_META, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["accession", "category", "substrate_class",
                                            "ec", "other_id", "reference", "included"])
        w.writeheader()
        w.writerows(meta_rows)

    print(f"Existing references kept: {len(existing)}")
    print(f"Toolkit entries: {len(data)} total, {n_skipped_dup} already in our database "
          f"(Q1JU72), {n_added} newly added")
    print(f"  ToolkitTierA (defluorination/dehalogenation enzymes): "
          f"{n_added - n_perfluoro - n_anion}")
    print(f"  PerfluoroTierA (genuine per/polyfluoro-acting complex): {n_perfluoro}")
    print(f"  FluTransport (fluoride export/resistance, not C-F cleavage): {n_anion}")
    print(f"\nTotal merged search database: {len(existing) + n_added} sequences")
    print(f"Wrote {OUT_FASTA}")
    print(f"Wrote {OUT_META}")

if __name__ == "__main__":
    main()
