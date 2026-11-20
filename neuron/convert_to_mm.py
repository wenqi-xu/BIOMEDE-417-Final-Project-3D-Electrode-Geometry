"""Convert COMSOL exports whose length unit is m into mm, and save them as
comsol_exports/{name}.txt so run_threshold.py can read them.
V is not changed (still V for I_0 = 1 mA). The original files are not touched.

Usage (from the neuron folder):
    python convert_to_mm.py ../comsol_exports/raw_m_units/*
"""
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(os.path.dirname(HERE), "comsol_exports")
os.makedirs(OUT_DIR, exist_ok=True)

for f in sys.argv[1:]:
    header = [line for line in open(f) if line.startswith("%")]
    if not any("Length unit" in line and line.split(":")[-1].strip() == "m" for line in header):
        print(f"skip {f}: length unit is not m")
        continue
    d = np.loadtxt(f, comments="%")
    d[:, :3] *= 1e3                                            # m -> mm
    base = os.path.basename(f)
    name = base if base.endswith(".txt") else base + ".txt"    # names contain dots (d0.05mm), so no splitext
    out = os.path.join(OUT_DIR, name)
    np.savetxt(out, d, fmt="%.10g", delimiter="\t",
               header=f"converted from {os.path.basename(f)} (m -> mm)\nLength unit: mm\nx\ty\tz\tV (V)",
               comments="% ")
    print(f"{f} -> {out}")
