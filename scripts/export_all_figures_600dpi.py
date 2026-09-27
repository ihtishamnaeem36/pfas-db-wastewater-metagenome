"""
Export all figures (Fig 1 through Fig 12) to true 600 DPI PNGs for ultra-high-resolution
publication artwork compliance (Nature / Science / Elsevier standards).
"""
from pathlib import Path
import fitz  # PyMuPDF

ROOT = Path(__file__).resolve().parents[1]
FIG_DIR = ROOT / "figures"

DPI = 600
ZOOM = DPI / 72.0  # standard PDF point size is 72 dpi
matrix = fitz.Matrix(ZOOM, ZOOM)

pdf_files = sorted(FIG_DIR.glob("Fig*.pdf"))
print(f"Found {len(pdf_files)} PDF figures to export at {DPI} DPI:")

for pdf_path in pdf_files:
    doc = fitz.open(str(pdf_path))
    page = doc[0]
    pix = page.get_pixmap(matrix=matrix, alpha=False)
    png_path = FIG_DIR / (pdf_path.stem + ".png")
    pix.save(str(png_path))
    print(f"  Exported {png_path.name}: {pix.width}x{pix.height} px ({pix.size / 1024:.1f} KB)")
    doc.close()

print(f"\nAll figures exported successfully at {DPI} DPI.")
