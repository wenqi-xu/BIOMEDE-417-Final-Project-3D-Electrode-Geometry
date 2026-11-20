"""Convert COMSOL exports in m to mm and save them to ../comsol_exports. Usage: python convert_to_mm.py file1 file2 ..."""
import os
import sys
import numpy as np

OUT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "comsol_exports")
os.makedirs(OUT_DIR, exist_ok=True)

for f in sys.argv[1:]:
    header = [line for line in open(f) if line.startswith("%")]
    if not any("Length unit" in line and line.split(":")[-1].strip() == "m" for line in header):
        print(f"skip {f}: length unit is not m")
        continue
    d = np.loadtxt(f, comments="%")
    d[:, :3] *= 1e3
    base = os.path.basename(f)
    out = os.path.join(OUT_DIR, base if base.endswith(".txt") else base + ".txt")
    np.savetxt(out, d, fmt="%.10g", delimiter="\t",
               header=f"converted from {base} (m -> mm)\nLength unit: mm\nx\ty\tz\tV (V)", comments="% ")
    print(f"{f} -> {out}")
