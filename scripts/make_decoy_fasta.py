"""Generate a whole-catalog residue-shuffled decoy FASTA for permutation-null testing
of the REAL DIAMOND/HMMER screen (fixed seed=42, same approach as the earlier custom-
pipeline null model, now applied to validate the real tools instead)."""
import random
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CATALOG = ROOT.parent / "Phase0_corrected" / "NR.protein.corrected.fa"
OUT = ROOT / "data" / "decoy_shuffled_catalog.fa"

def main():
    rng = random.Random(42)
    n = 0
    with open(CATALOG, encoding="utf-8", errors="replace") as fin, \
         open(OUT, "w", encoding="utf-8") as fout:
        name, seq_lines = None, []
        def flush(nm, sq):
            nonlocal n
            if not nm:
                return
            chars = list(sq)
            rng.shuffle(chars)
            fout.write(f">{nm}\n")
            s = "".join(chars)
            for i in range(0, len(s), 60):
                fout.write(s[i:i+60] + "\n")
            n += 1
        for line in fin:
            line = line.rstrip("\n")
            if line.startswith(">"):
                flush(name, "".join(seq_lines))
                name = line[1:].split()[0]
                seq_lines = []
            else:
                seq_lines.append(line.strip())
        flush(name, "".join(seq_lines))
    print(f"Wrote {n:,} shuffled decoy sequences to {OUT}")

if __name__ == "__main__":
    main()
